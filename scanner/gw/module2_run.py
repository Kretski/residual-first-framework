#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module2_run.py — Module 1 v2.0 driver (MODULE1_DESIGN_v2.0.md). The v1.0 files are imported
UNCHANGED; v1.0 remains stopped (tag v1.0-stopped-0a).

The single statistical change (design §3, amended before freeze): every belt and
expected-interval trial contains, besides Gaussian noise and calibration errors, the term
P[h(theta) - h_ref] of a posterior sample NOT used to build the subspace of that belt:
the belts of subspace A use the first 400 samples of set B, the belts of subspace B the first
400 samples of set A (cross-sets; the published LALInference labels have too few samples for
a separate set). This is done by replacing module1_catalog.event_trials (which the frozen
catalog stage and module1_offsource.event_c0 call) with a version that adds this term.

Stages:
  dset       per pair: the first 400 samples of set A and of set B (positions of
             module1_event.build: fixed permutation, seed 20261001), generated and projected as
             the off-source injections; stored in --cache2
  catalog    the FROZEN module1_run.stage_catalog with the cross-set term in every trial
             (belt trial j uses sample j mod 400, check trial j uses (j + 200) mod 400)
  offsource  as the frozen off-source stage, with NEW segments (ladder index k >= 8, none
             overlapping a v1.0 off-source segment of the pair) and NEW injections (set E:
             the 8 usable samples following set C of v1.0); c0 from the v2.0 event belt
  hash       v2.0 code hash (v1.0 code hash + this file)

The v1.0 cache is loaded with the v1.0 code hash, so it stays valid. Outputs carry both hashes
and the v2.0 status ('official-v2' only if FREEZE_v2.0.json matches, else 'trial-v2').
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time

import numpy as np

import module1_run as run                     # installs the on-source safeguard
import h5py                                   # noqa: E402
import module1_catalog as mc                  # noqa: E402
import module1_event as me                    # noqa: E402
import module1_offsource as mo                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
N_D, D_SHIFT, K_MIN = 400, 200, 8
_ORIG_LOAD = run.load
_ORIG_TRIALS = mc.event_trials


# ---------------------------------------------------------------- hashes / status ----
def v2_hash():
    v1, per = run.code_hash()
    own = run.sha256_file(os.path.join(HERE, "module2_run.py"), normalise_eol=True)
    return hashlib.sha256((v1 + own).encode()).hexdigest(), v1, own


def v2_status(h2, manifest):
    if not os.path.exists(manifest):
        return "trial-v2", "no v2.0 freeze manifest"
    m = json.load(open(manifest, encoding="utf-8"))
    if m.get("code_sha256_v2") != h2:
        return "trial-v2", f"code hash differs from {m.get('tag', '?')}"
    return "official-v2", f"code hash matches {m.get('tag', '?')}"


# ---------------------------------------------------------------- set D -------------
def d_paths(args, ev, label):
    stem = os.path.join(args.cache2, f"{ev}__{label.replace(':', '_')}")
    return stem + "_D.npz", stem + "_D.json"


def pair_setup(args, ev, label, r5, e):
    import f3_derivatives as f3
    meta, ifos, fb = e["meta"], e["meta"]["ifos"], e["fb"]
    with h5py.File(os.path.join(args.pe_dir, r5["file"].split(" ")[0]), "r") as fh:
        sample, cfg, _ = f3.read_label(fh, label)
        st = f3.settings(cfg)
        tab = me.posterior_table(fh, label)
    seg_on = sample["geocent_time"] + me.rc.POST_TRIGGER - st["duration"]
    f, _, _ = me.rc.polarizations(label, sample, st)
    band = (f >= st["f_an"]) & (f <= st["f_high"])
    if not np.allclose(f[band], fb):
        raise RuntimeError(f"{ev} {label}: frequency band differs from the cache")
    calf = {i: e["sb"][i].mean_factor(fb) for i in ifos}

    def network_h(s):
        _, sp, sc = me.rc.polarizations(label, s, st)
        return {i: me.rc.project(f, sp, sc, i, s, seg_on)[band] * calf[i] for i in ifos}
    idx = np.random.default_rng(me.SEED).permutation(len(tab["log_likelihood"]))
    return st, tab, idx, seg_on, band, network_h


def stage_dset(args, chash1, h2, status):
    """First N_D samples of set A and of set B, exactly at the positions used by
    module1_event.build (A = first n usable samples, B = next n usable)."""
    os.makedirs(args.cache2, exist_ok=True)
    for ev, label, r5 in run.select(args):
        npz, js = d_paths(args, ev, label)
        if os.path.exists(npz) and os.path.exists(js):
            m = json.load(open(js))
            if m.get("code_sha256_v1") == chash1 and m.get("n_D") == N_D and m.get("rule") == "cross-sets":
                print(f"=== {ev} {label}: cross-set samples cached, skipped", flush=True)
                continue
        t0 = time.time()
        e = _ORIG_LOAD(args, ev, label, r5, chash1)
        st, tab, idx, seg_on, band, network_h = pair_setup(args, ev, label, r5, e)
        n, fails0 = e["meta"]["n"], e["meta"]["fails"]
        href_s = e["net"].stack(e["href"])
        rows = len(href_s)
        LA = np.empty((N_D, rows), np.float32)    # from set A, used in the belts of subspace B
        LB = np.empty((N_D, rows), np.float32)    # from set B, used in the belts of subspace A
        pos, usable, fails, mA, mB = 0, 0, 0, 0, 0
        # walk the permutation as build() did; without build failures A = idx[:n], B = idx[n:2n]
        while (mA < N_D or mB < N_D) and pos < len(idx):
            j = int(idx[pos])
            pos += 1
            in_A = usable < n
            if fails0 == 0:
                if (in_A and mA >= N_D) or (not in_A and usable >= n + N_D):
                    usable += 1
                    continue
            s = {k: float(tab[k][j]) for k in me.NEED}
            try:
                vec = e["net"].stack(network_h(s)) - href_s
            except Exception:
                fails += 1
                continue
            if in_A and mA < N_D:
                LA[mA] = vec
                mA += 1
            elif not in_A and mB < N_D:
                LB[mB] = vec
                mB += 1
            usable += 1
            if (mA + mB) % 100 == 0:
                print(f"    cross-set samples {mA + mB}/{2 * N_D}  ({time.time() - t0:.0f} s)", flush=True)
        if mA < N_D or mB < N_D:
            raise RuntimeError(f"{ev} {label}: only {mA} set-A and {mB} set-B samples")
        rmsA = float(np.sqrt(np.mean(np.sum(LA.astype(np.float64) ** 2, axis=1))))
        rmsB = float(np.sqrt(np.mean(np.sum(LB.astype(np.float64) ** 2, axis=1))))
        np.savez_compressed(npz, LA=LA, LB=LB)
        json.dump({"event": ev, "label": label, "rule": "cross-sets", "n_D": N_D, "fails": fails,
                   "build_fails": fails0, "code_sha256_v1": chash1, "code_sha256_v2": h2, "status": status,
                   "rms_mismatch_A": rmsA, "rms_mismatch_B": rmsB}, open(js, "w"), indent=1)
        print(f"=== {ev} {label}: cross-set samples {N_D} + {N_D} (failures {fails}), rms mismatch "
              f"A {rmsA:.2f}  B {rmsB:.2f}  ({time.time() - t0:.0f} s)", flush=True)


# ---------------------------------------------------------------- patched pieces -----
def load_v2(args, ev, label, r5, chash):
    """v1.0 cache (loaded with the v1.0 hash) plus the set-D mismatch vectors."""
    e = _ORIG_LOAD(args, ev, label, r5, run.code_hash()[0])
    npz, js = d_paths(args, ev, label)
    if not (os.path.exists(npz) and os.path.exists(js)):
        raise RuntimeError(f"{ev} {label}: set D missing; run the dset stage first")
    d = np.load(npz)
    e["Dleak"] = {"A": d["LB"], "B": d["LA"]}      # cross-sets: never the subspace's own samples
    return e


def event_trials_v2(ev, grid, n_belt, n_exp, seed, which="A"):
    """module1_catalog.event_trials with the set-D mismatch term in every trial (design §3).
    The random-number sequence is exactly that of v1.0; the D term uses no random numbers."""
    L = ev["Dleak"][which]
    nD = len(L)
    A = ev[which]
    sc = me.ProfileScan(A, ev["K"], ev["href"], ev["net"], ev["fb"], grid)
    B = sc.b.astype(np.float64)
    G = sc.G
    del sc
    rng = np.random.default_rng(seed)
    cal = mc.cal_error_fn(ev, rng)
    nvec = A.Q.shape[0]
    out = []
    for ntr, shift in ((n_belt, 0), (n_exp, D_SHIFT)):
        V = np.empty((ntr, len(grid)))
        for j in range(ntr):
            n = rng.normal(size=nvec) + cal() + L[(j + shift) % nD].astype(np.float64)
            npj = n - A.Q @ (A.Q.T @ n)
            V[j] = B @ npj
        out.append(V)
    return G, out[0], out[1]


def install_patches():
    run.load = load_v2
    mc.event_trials = event_trials_v2


# ---------------------------------------------------------------- stages -------------
def stage_catalog(args, chash1, h2, status):
    install_patches()
    out = args.out or f"module2_catalog_{args.labels}.json"
    args.out = out
    run.stage_catalog(args, chash1, status)
    res = json.load(open(out))
    res.update({"design": "MODULE1_DESIGN_v2.0", "code_sha256_v2": h2, "status": status,
                "n_D": N_D, "D_shift_for_checks": D_SHIFT})
    json.dump(res, open(out, "w"), indent=1)
    print(f"v2.0: status {status}, code_sha256_v2 {h2}; written {out}")


def v1_segments(ev, label):
    out = []
    for fn in ("primary_offsource_search_eob.csv", "primary_offsource_search_xphm.csv",
               "primary_offsource_confirm_eob.csv", "primary_offsource_confirm_xphm.csv",
               "diag_0a_gate2.csv"):
        if os.path.exists(fn):
            for r in csv.DictReader(open(fn)):
                if r["event"] == ev and r["label"] == label:
                    out.append(float(r["seg_start"]))
    return sorted(set(out))


def stage_offsource(args, chash1, h2, status):
    from scipy.signal import resample_poly
    from scipy.signal.windows import tukey
    install_patches()
    out = args.out or f"module2_offsource_{args.labels}.csv"
    keys = ["event", "label", "status", "code_sha256_v1", "code_sha256_v2", "seg_start", "inj", "q0", "c0",
            "covered", "lambda_hat", "dt_ms", "abs_a", "e_index"]
    new = not os.path.exists(out)
    fout = open(out, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fout, fieldnames=keys)
    if new:
        w.writeheader()
    manifest = {(r["event"], r["ifo"]): r for r in
                csv.DictReader(open(os.path.join(args.strain_dir, "strain_manifest.csv"), encoding="utf-8"))}
    summary = {"C": [0, 0], "ref": [0, 0]}
    for k, (ev, label, r5) in enumerate(run.select(args)):
        t0 = time.time()
        e = load_v2(args, ev, label, r5, chash1)
        ifos, fb = e["meta"]["ifos"], e["fb"]
        st, tab, idx, seg_on, band, network_h = pair_setup(args, ev, label, r5, e)
        dur, srate = st["duration"], st["srate"]
        # set E: the 8 usable samples following set C of v1.0 (C = first 8 usable after A and B)
        pos, nCseen, hE = 2 * e["meta"]["n"] + e["meta"]["fails"], 0, []
        while len(hE) < mo.N_C and pos < len(idx):
            j = int(idx[pos])
            pos += 1
            s = {kk: float(tab[kk][j]) for kk in me.NEED}
            try:
                h = network_h(s)
            except Exception:
                continue
            if nCseen < mo.N_C:
                nCseen += 1
                continue
            hE.append((j, h))
        if len(hE) < mo.N_C:
            print(f"  {ev} {label}: only {len(hE)} set-E waveforms", flush=True)
        # NEW segments: ladder index k >= 8, no overlap with any v1.0 segment of the pair
        old = v1_segments(ev, label)
        windows = mo.forbidden_windows(args.events_list, dur)
        win = tukey(int(round(dur * srate)), alpha=2 * me.rc.ROLL_OFF / dur)
        strain = {}
        for i in ifos:
            x, dt, t0f = run._ORIG_READ_STRAIN(os.path.join(args.strain_dir, manifest[(ev, i)]["file"]))
            if abs(1 / dt - srate) > 1e-6:
                x = resample_poly(x, int(srate), int(round(1 / dt)))
                dt = 1.0 / srate
            strain[i] = (x, dt, t0f)
        segs = []
        for kk in range(K_MIN, K_MIN + 6 * mo.N_OFF):
            for side in (-1, +1):
                if len(segs) >= mo.N_OFF:
                    break
                s0 = seg_on + side * kk * (dur + mo.OFF_GAP)
                if any(abs(s0 - o) < dur for o in old):
                    continue
                ok, _ = mo.segment_allowed(s0, dur, windows)
                if not ok:
                    continue
                try:
                    dseg = {i: (np.fft.rfft(run._ORIG_SEGMENT(strain[i][0], strain[i][1], strain[i][2], s0, dur)
                                            * win) * strain[i][1])[band] for i in ifos}
                    segs.append((s0, dseg))
                except ValueError:
                    pass
        del strain
        own_grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * e["sigma_lin"]
        c0, scan, tz = mo.event_c0(e, own_grid, args.belt_trials, me.SEED + 11 + 1000 * k)
        Seff = {i: 4 * (fb[1] - fb[0]) / e["net"].sw[i] ** 2 for i in ifos}
        for m, (s0, dseg) in enumerate(segs):
            for inj in ("C", "ref"):
                if inj == "C":
                    if m >= len(hE):
                        continue
                    cj, hinj = hE[m]
                else:
                    cj, hinj = -1, e["href"]
                d = {i: dseg[i] + hinj[i] for i in ifos}
                r, tsh, a = mo.refine(d, e["href"], Seff, fb[1] - fb[0], fb=fb)
                q0, g = mo.q0_of_residual(e["net"].stack(r), e, scan, tz)
                cov = q0 <= c0
                summary[inj][0] += cov
                summary[inj][1] += 1
                w.writerow({"event": ev, "label": label, "status": status, "code_sha256_v1": chash1,
                            "code_sha256_v2": h2, "seg_start": s0, "inj": inj, "q0": q0, "c0": c0,
                            "covered": cov, "lambda_hat": own_grid[g], "dt_ms": tsh * 1e3, "abs_a": abs(a),
                            "e_index": cj})
        fout.flush()
        print(f"{ev:18s} {label:30s} new segments {len(segs)}  set E {len(hE)}  c0 {c0:.2f}  "
              f"({time.time() - t0:.0f} s)", flush=True)
    fout.close()
    for inj in ("C", "ref"):
        nc, n = summary[inj]
        if n:
            p = nc / n
            lo, hi = mo.binomial_interval(0.90, n)
            verdict = ("STOP: below the 99% binomial interval around 0.90" if p < lo else
                       "conservative (above the interval)" if p > hi else "within the interval")
            print(f"{'E' if inj == 'C' else inj:4s} coverage {nc}/{n} = {p:.3f}; 99% interval around 0.90: "
                  f"{lo:.3f}–{hi:.3f} -> {verdict if inj == 'C' else 'diagnostic'}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["dset", "catalog", "offsource", "hash"])
    ap.add_argument("--events", default="GW150914,GW230627_015337")
    ap.add_argument("--labels", default="both", choices=["both", "xphm", "eob"])
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    ap.add_argument("--residual-check", default="residual_check.csv")
    ap.add_argument("--cache", default=os.path.expanduser("~/gw/cache_module1"))
    ap.add_argument("--cache2", default=os.path.expanduser("~/gw/cache_module2"))
    ap.add_argument("--n", type=int, default=800)
    ap.add_argument("--belt-trials", type=int, default=5000)
    ap.add_argument("--exp-trials", type=int, default=2000)
    ap.add_argument("--half-width", type=float, default=8.0)
    ap.add_argument("--step", type=float, default=0.05)
    ap.add_argument("--strain-dir", default=os.path.expanduser("~/gwdata/strain"))
    ap.add_argument("--events-list", default="event_selection/event_list_v1.csv")
    ap.add_argument("--freeze-manifest-v2", default=os.path.join(HERE, "FREEZE_v2.0.json"))
    ap.add_argument("--write-manifest", default="", metavar="TAG")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    h2, chash1, own = v2_hash()
    status, why = v2_status(h2, args.freeze_manifest_v2)
    print(f"v1.0 code SHA-256 {chash1}\nmodule2_run.py SHA-256 {own}\nv2.0 code SHA-256 {h2}\n"
          f"status: {status} ({why})", flush=True)
    if args.stage == "hash":
        if args.write_manifest:
            if os.path.exists(args.freeze_manifest_v2):
                raise SystemExit(f"{args.freeze_manifest_v2} exists; refusing to overwrite")
            json.dump({"tag": args.write_manifest, "code_sha256_v2": h2, "code_sha256_v1": chash1,
                       "module2_run_sha256": own, "n_D": N_D, "D_rule": "cross-sets", "D_shift": D_SHIFT,
                       "k_min_segments": K_MIN, "design": "MODULE1_DESIGN_v2.0.md",
                       "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
                      open(args.freeze_manifest_v2, "w", encoding="utf-8"), indent=1)
            print(f"written {args.freeze_manifest_v2}")
        return
    if args.stage == "offsource" and args.belt_trials == 5000:
        args.belt_trials = 2000          # event-only belt at Lambda = 0, as in v1.0
    {"dset": stage_dset, "catalog": stage_catalog, "offsource": stage_offsource}[args.stage](
        args, chash1, h2, status)


if __name__ == "__main__":
    main()
