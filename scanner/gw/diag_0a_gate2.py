#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_0a_gate2.py — Gate 2 of DIAG_0a_SPEC.md (pre-specified diagnostic of the rule 0a
failure of Module 1 v1.0). It cannot change the v1.0 result (tag v1.0-stopped-0a).

Pairs (fixed by the specification): the four failing SEOBNRv4PHM search events, the
four controls (smallest |dt_ms| in residual_check.csv), and, as a cross-model check,
the IMRPhenomXPHM labels of the four failing events.

For every pair, on the SAME off-source segments, set-C injections, belts and seeds as
the v1.0 off-source stage (module1_run.stage_offsource):
  v1.0  baseline from the v1.0 cache (module1_run.load); its 'covered' flags must equal
        the stored ones in primary_offsource_search_{eob,xphm}.csv EXACTLY, otherwise
        the diagnostic stops;
  (a)   reference = posterior medoid: among the first 200 set-B samples, the one with the
        smallest mean raw whitened distance to the others; templates, calibration
        directions, free directions and the set-A/B subspaces are rebuilt around it;
  (b)   the set-A/B posterior samples are aligned in time and phase to the v1.0 reference
        (module1_offsource.refine, noise-free; amplitude untouched) before the
        principal-component subspaces are built; everything else as v1.0.
The set-A/B waveforms are generated once per pair, in the order and with the samples of
module1_event.build (fixed permutation, seed 20261001), and shared by (a) and (b).

Reading (fixed in DIAG_0a_SPEC.md): a variant supports its hypothesis if the pooled
set-C coverage of the four failing events reaches at least 24/32 while the pooled
coverage of the four controls does not fall by more than 3/32 from its recomputed
baseline. Both variants are reported whatever the outcome.

  python diag_0a_gate2.py --out diag_0a_gate2.csv
Resumable: pairs already complete in --out are skipped.
"""
import argparse
import csv
import os
import time

import h5py
import numpy as np
from scipy.signal import resample_poly
from scipy.signal.windows import tukey

import f3_derivatives as f3
import module1_event as me
import module1_offsource as mo
import module1_run as run

FAILING = ["GW150914", "GW170814", "GW190706_222641", "GW200129_065458"]
N_MEDOID = 200
THRESH_FAIL, MAX_CTRL_DROP = 24, 3
KEYS = ["event", "label", "role", "variant", "seg_start", "inj", "q0", "c0", "covered", "dt_ms", "sigma_lin", "k_A"]


def driver_args():
    return argparse.Namespace(
        pe_dir=os.path.expanduser("~/gwdata/pe"), f5="f5/f5_results.csv",
        residual_check="residual_check.csv", cache=os.path.expanduser("~/gw/cache_module1"), n=800,
        strain_dir=os.path.expanduser("~/gwdata/strain"), events_list="event_selection/event_list_v1.csv",
        half_width=8.0, step=0.05, belt_trials=2000)


def controls(args):
    prim = [e for e in open("primary_search_eob.txt").read().split(",") if e]
    rows = [r for r in csv.DictReader(open(args.residual_check, encoding="utf-8"))
            if "SEOBNRv4PHM" in r["label"] and r["event"] in prim and r["event"] not in FAILING]
    rows.sort(key=lambda r: abs(float(r["dt_ms"])))
    return [r["event"] for r in rows[:4]]


def stored_flags(ev, label):
    fn = "primary_offsource_search_eob.csv" if "SEOBNR" in label else "primary_offsource_search_xphm.csv"
    out = {}
    for r in csv.DictReader(open(fn)):
        if r["event"] == ev and r["label"] == label:
            out[(round(float(r["seg_start"]), 4), r["inj"])] = r["covered"] in ("True", "1")
    return out


def run_index(ev, label):
    """Position of the pair in the v1.0 off-source run (it fixes the belt seed)."""
    fn = "primary_search_eob.txt" if "SEOBNR" in label else "primary_search_xphm.txt"
    return [e for e in open(fn).read().split(",") if e].index(ev)


def geometry(href, sb, net, fb, K, ifos, D_A, D_B):
    """Estimators A and B around a reference, as module1_event.build does."""
    T = {p: net.stack({i: 1j * K * fb ** p * href[i] for i in ifos}) for p in (2, 3, 4)}
    tdir = net.stack({i: -2j * np.pi * fb * href[i] for i in ifos})
    adir = net.stack(href)
    pdir = net.stack({i: 1j * href[i] for i in ifos})
    J = np.column_stack([net.single(i, d) for i in ifos for d in sb[i].directions(fb, href[i])])
    est, ks = {}, {}
    for name, D in (("A", D_A), ("B", D_B)):
        U, k = me.pca_basis(D, me.FRAC)
        Q = me.orthonormal([U, tdir, adir, pdir])
        est[name], ks[name] = me.Estimator(Q, J, T), k
    return est, ks, T


def offsource(e, segs, hC, k, args, ifos, fb):
    """The v1.0 off-source chain for one geometry; returns rows (seg, inj, q0, c0, covered, dt)."""
    own_grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * e["sigma_lin"]
    c0, scan, tz = mo.event_c0(e, own_grid, args.belt_trials, me.SEED + 11 + 1000 * k)
    Seff = {i: 4 * (fb[1] - fb[0]) / e["net"].sw[i] ** 2 for i in ifos}
    out = []
    for m, (s0, dseg) in enumerate(segs):
        for inj in ("C", "ref"):
            if inj == "C":
                if m >= len(hC):
                    continue
                hinj = hC[m][1]
            else:
                hinj = e["href"]
            d = {i: dseg[i] + hinj[i] for i in ifos}
            r, tsh, _ = mo.refine(d, e["href"], Seff, fb[1] - fb[0], fb=fb)
            q0, _ = mo.q0_of_residual(e["net"].stack(r), e, scan, tz)
            out.append((s0, inj, q0, c0, bool(q0 <= c0), tsh * 1e3))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="diag_0a_gate2.csv")
    cli = ap.parse_args()
    args = driver_args()
    chash, _ = run.code_hash()
    rows5 = {r["commonName"]: r for r in csv.DictReader(open(args.f5, encoding="utf-8"))}
    ctrl = controls(args)
    pairs = [(e, rows5[e]["eob_label"], "failing") for e in FAILING] + \
            [(e, rows5[e]["eob_label"], "control") for e in ctrl] + \
            [(e, rows5[e]["xphm_label"], "failing-XPHM") for e in FAILING]
    print(f"code SHA-256 {chash}\ncontrols: {ctrl}", flush=True)
    done = set()
    if os.path.exists(cli.out):
        for r in csv.DictReader(open(cli.out)):
            done.add((r["event"], r["label"], r["variant"]))
    fo = open(cli.out, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fo, fieldnames=KEYS)
    if not done:
        w.writeheader()
    manifest = {(r["event"], r["ifo"]): r for r in
                csv.DictReader(open(os.path.join(args.strain_dir, "strain_manifest.csv"), encoding="utf-8"))}

    for ev, label, role in pairs:
        if all((ev, label, v) in done for v in ("v1.0", "a", "b")):
            print(f"\n=== {ev} {label}: done, skipped", flush=True)
            continue
        t0 = time.time()
        print(f"\n=== {ev} {label} ({role})", flush=True)
        r5 = rows5[ev]
        k = run_index(ev, label)
        e0 = run.load(args, ev, label, r5, chash)                       # v1.0 geometry (cache)
        meta, ifos, fb, net, sb, K = e0["meta"], e0["meta"]["ifos"], e0["fb"], e0["net"], e0["sb"], e0["K"]
        with h5py.File(os.path.join(args.pe_dir, r5["file"].split(" ")[0]), "r") as fh:
            sample, cfg, _ = f3.read_label(fh, label)
            st = f3.settings(cfg)
            tab = me.posterior_table(fh, label)
        dur, srate = st["duration"], st["srate"]
        seg_on = sample["geocent_time"] + me.rc.POST_TRIGGER - dur
        f, _, _ = me.rc.polarizations(label, sample, st)
        band = (f >= st["f_an"]) & (f <= st["f_high"])
        if not np.allclose(f[band], fb):
            raise RuntimeError(f"{ev} {label}: frequency band differs from the cache")
        calf = {i: sb[i].mean_factor(fb) for i in ifos}
        Seff = {i: 4 * (fb[1] - fb[0]) / net.sw[i] ** 2 for i in ifos}

        def network_h(s):
            fp, sp, sc = me.rc.polarizations(label, s, st)
            if len(fp) != len(f) or not np.allclose(fp, f):
                raise ValueError("frequency grid mismatch")
            return {i: me.rc.project(f, sp, sc, i, s, seg_on)[band] * calf[i] for i in ifos}

        # set-C injections and off-source segments, exactly as stage_offsource
        idx = np.random.default_rng(me.SEED).permutation(len(tab["log_likelihood"]))
        pos = 2 * meta["n"] + meta["fails"]
        hC = []
        while len(hC) < mo.N_C and pos < len(idx):
            j = int(idx[pos])
            pos += 1
            sC = {kk: float(tab[kk][j]) for kk in me.NEED}
            try:
                hC.append((j, network_h(sC)))
            except Exception:
                continue
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
        for kk in range(1, 4 * mo.N_OFF):
            for side in (-1, +1):
                if len(segs) >= mo.N_OFF:
                    break
                s0 = seg_on + side * kk * (dur + mo.OFF_GAP)
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

        def write(variant, rows, sigma_lin, kA):
            for s0, inj, q0, c0, cov, dtms in rows:
                w.writerow({"event": ev, "label": label, "role": role, "variant": variant, "seg_start": s0,
                            "inj": inj, "q0": q0, "c0": c0, "covered": cov, "dt_ms": dtms,
                            "sigma_lin": sigma_lin, "k_A": kA})
            fo.flush()

        # ---- v1.0 baseline and the reproduction gate
        base = offsource(e0, segs, hC, k, args, ifos, fb)
        stored = stored_flags(ev, label)
        mism = [(s0, inj) for s0, inj, _, _, cov, _ in base if stored.get((round(s0, 4), inj)) != cov]
        nC = sum(cov for _, inj, _, _, cov, _ in base if inj == "C")
        print(f"  v1.0 baseline: C coverage {nC}/{sum(inj == 'C' for _, inj, *_ in base)}; "
              f"reproduction mismatches {len(mism)}", flush=True)
        if mism or len(base) != len(stored):
            fo.close()
            raise SystemExit(f"STOP: v1.0 baseline does not reproduce the stored flags for {ev} {label} "
                             f"({len(mism)} mismatches, {len(base)} vs {len(stored)} rows)")
        write("v1.0", base, e0["sigma_lin"], meta.get("k", {}).get("A", "") if isinstance(meta.get("k"), dict) else "")

        # ---- set-A/B waveforms, generated once (order and samples of module1_event.build)
        href0 = e0["href"]
        href0_s = net.stack(href0)
        rows_n = len(href0_s)
        n = meta["n"]
        H = {"A": np.empty((rows_n, n), np.float32), "B": np.empty((rows_n, n), np.float32)}
        HL = {"A": np.empty((rows_n, n), np.float32), "B": np.empty((rows_n, n), np.float32)}
        Bsamples = []
        pos, fails, tg = 0, 0, time.time()
        for name in ("A", "B"):
            m = 0
            while m < n and pos < len(idx):
                j = int(idx[pos])
                pos += 1
                s = {kk: float(tab[kk][j]) for kk in me.NEED}
                try:
                    hs = network_h(s)
                except Exception:
                    fails += 1
                    continue
                H[name][:, m] = net.stack(hs)
                _, tsh, a = mo.refine(hs, href0, Seff, fb[1] - fb[0], fb=fb)
                rot = np.exp(2j * np.pi * fb * tsh - 1j * np.angle(a))      # time and phase only
                HL[name][:, m] = net.stack({i: hs[i] * rot for i in ifos})
                if name == "B" and m < N_MEDOID:
                    Bsamples.append(s)
                m += 1
                if (pos % 100) == 0:
                    el = time.time() - tg
                    print(f"    samples {pos}/{2 * n}  ({el:.0f} s)", flush=True)
            if m < n:
                raise RuntimeError(f"set {name}: only {m} samples")
        if fails != meta["fails"]:
            print(f"  WARNING: {fails} generation failures, cache recorded {meta['fails']}", flush=True)

        # ---- variant (b): aligned samples, v1.0 reference
        est_b, ks_b, T_b = geometry(href0, sb, net, fb, K, ifos,
                                    HL["A"] - href0_s[:, None].astype(np.float32),
                                    HL["B"] - href0_s[:, None].astype(np.float32))
        e_b = dict(e0, A=est_b["A"], B=est_b["B"], est=est_b, T=T_b, sigma_lin=est_b["A"].sigma(3))
        rows_b = offsource(e_b, segs, hC, k, args, ifos, fb)
        print(f"  (b) aligned: C coverage {sum(c for _, i, _, _, c, _ in rows_b if i == 'C')}/"
              f"{sum(i == 'C' for _, i, *_ in rows_b)}, k_A {ks_b['A']}, sigma_lin {e_b['sigma_lin']:.3e}", flush=True)
        write("b", rows_b, e_b["sigma_lin"], ks_b["A"])
        del HL

        # ---- variant (a): posterior medoid among the first 200 set-B samples
        X = H["B"][:, :N_MEDOID].astype(np.float64)
        G = X.T @ X
        d2 = np.diag(G)[:, None] + np.diag(G)[None, :] - 2 * G
        dist = np.sqrt(np.clip(d2, 0, None)).mean(axis=1)
        jm = int(np.argmin(dist))
        href_a = network_h(Bsamples[jm])
        href_a_s = net.stack(href_a).astype(np.float32)
        est_a, ks_a, T_a = geometry(href_a, sb, net, fb, K, ifos,
                                    H["A"] - href_a_s[:, None], H["B"] - href_a_s[:, None])
        e_a = dict(e0, A=est_a["A"], B=est_a["B"], est=est_a, T=T_a, href=href_a, sigma_lin=est_a["A"].sigma(3))
        rows_a = offsource(e_a, segs, hC, k, args, ifos, fb)
        print(f"  (a) medoid (set-B sample {jm}, distance of v1.0 reference to it "
              f"{np.linalg.norm(href0_s - href_a_s):.2f}): C coverage {sum(c for _, i, _, _, c, _ in rows_a if i == 'C')}/"
              f"{sum(i == 'C' for _, i, *_ in rows_a)}, k_A {ks_a['A']}, sigma_lin {e_a['sigma_lin']:.3e}", flush=True)
        write("a", rows_a, e_a["sigma_lin"], ks_a["A"])
        del H
        print(f"  ({time.time() - t0:.0f} s)", flush=True)
    fo.close()

    # ---- summary and the pre-specified reading
    res = list(csv.DictReader(open(cli.out)))
    def pooled(role, variant):
        r = [x for x in res if x["role"] == role and x["variant"] == variant and x["inj"] == "C"]
        return sum(x["covered"] == "True" for x in r), len(r)
    print("\nGate 2 summary (set-C coverage, pooled):")
    print(f"  {'variant':8s} {'failing':>10s} {'controls':>10s} {'failing-XPHM':>14s}")
    for v in ("v1.0", "a", "b"):
        f_, c_, x_ = pooled("failing", v), pooled("control", v), pooled("failing-XPHM", v)
        print(f"  {v:8s} {f_[0]:>4d}/{f_[1]:<5d} {c_[0]:>4d}/{c_[1]:<5d} {x_[0]:>6d}/{x_[1]:<7d}")
    cb = pooled("control", "v1.0")[0]
    for v, hyp in (("a", "(a) mis-centred reference"), ("b", "(b) parameter-dependent time offsets")):
        f_, c_ = pooled("failing", v), pooled("control", v)
        ok = f_[0] >= THRESH_FAIL and (cb - c_[0]) <= MAX_CTRL_DROP
        print(f"  {hyp}: failing {f_[0]}/{f_[1]} (needs >= {THRESH_FAIL}/32), controls {c_[0]} vs baseline {cb} "
              f"(drop <= {MAX_CTRL_DROP}) -> {'SUPPORTED' if ok else 'NOT SUPPORTED'}")
    print("Both variants are reported whatever the outcome; neither is adopted on this basis alone.")


if __name__ == "__main__":
    main()
