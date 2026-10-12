#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
supplement_validation.py — the two post-result validation checks of
MODULE1_SUPPLEMENT_VALIDATION_REGISTRATION.md. Off-source data only; cannot change the committed
on-source sets.

  V1  coverage at a non-zero true Lambda (+-2 sigma_lin of the event, alternating by segment):
      injection = set-F waveform x exp(i K f^3 Lambda_t), the signal model of the estimator;
      covered if q(Lambda_t) <= c(Lambda_t) of the event's own belt.
  V2  waveform-model leakage: injection = SEOBNRv5PHM waveform with the set-F parameters
      (Lambda = 0), analysed with the IMRPhenomXPHM reference and chain; covered if q(0) <= c(0);
      the scan estimate is recorded in units of the event's sigma_lin.

Common: IMRPhenomXPHM pairs of Addendum D; NEW segments (ladder index k >= 16, none overlapping a
v1.0 / v2.0 / diagnostic off-source segment of the pair); set F = the 8 usable samples after set
E (C = first 8 usable after A and B, E = next 8, F = next 8); chain and belts as the v2.0
off-source stage (own grid, 2000 belt trials, seed SEED + 11 + 1000 k), with the full Neyman
belt of the event so that c(Lambda_t) is available.

Stages:
  test   reproduces one stored v2.0 off-source row through this code path (q0 and c0), and
         probes SEOBNRv5PHM generation for every pair (one set-F sample); run before freezing
  v1     check V1 for both cohorts
  v2     check V2 for both cohorts
  hash   code hash; --write-manifest TAG writes FREEZE_supplement_validation.json
"""
import argparse
import copy
import csv
import hashlib
import json
import os
import time

import numpy as np

import supplement_onsource as so
import module2_run as m2                     # noqa: E402
import module1_run as run                    # noqa: E402
import module1_catalog as mc                 # noqa: E402
import module1_event as me                   # noqa: E402
import module1_offsource as mo               # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
K_MIN_VAL, N_SEG, N_F, T_SIGMA = 16, 8, 8, 2.0
COHORTS = ("search", "confirm")


def val_hash():
    hs, _, _ = so.supp_hash()
    own = run.sha256_file(os.path.join(HERE, "supplement_validation.py"), normalise_eol=True)
    return hashlib.sha256((hs + own).encode()).hexdigest(), hs, own


def val_status(hv, manifest):
    if not os.path.exists(manifest):
        return "trial-val", "no validation freeze manifest"
    m = json.load(open(manifest, encoding="utf-8"))
    if m.get("code_sha256_val") != hv:
        return "trial-val", f"code hash differs from {m.get('tag', '?')}"
    return "official-val", f"code hash matches {m.get('tag', '?')}"


def used_segments(ev, label):
    out = list(m2.v1_segments(ev, label))
    for fn in ("v2_offsource_search_xphm.csv", "v2_offsource_confirm_xphm.csv"):
        if os.path.exists(fn):
            for r in csv.DictReader(open(fn)):
                if r["event"] == ev and r["label"] == label:
                    out.append(float(r["seg_start"]))
    return sorted(set(out))


def set_F(idx, tab, e, network_h):
    """C = first 8 usable after A and B, E = next 8, F = next 8 (usable = generated in-process)."""
    pos, seen, F = 2 * e["meta"]["n"] + e["meta"]["fails"], 0, []
    while len(F) < N_F and pos < len(idx):
        j = int(idx[pos])
        pos += 1
        s = {kk: float(tab[kk][j]) for kk in me.NEED}
        try:
            h = network_h(s)
        except Exception:
            continue
        if seen < 16:
            seen += 1
            continue
        F.append((j, s, h))
    return F


def new_segments(args, ev, label, e, st, band, seg_on):
    old = used_segments(ev, label)
    dur = st["duration"]
    windows = mo.forbidden_windows(args.events_list, dur)
    segs = []
    for kk in range(K_MIN_VAL, K_MIN_VAL + 8 * N_SEG):
        for side in (-1, +1):
            if len(segs) >= N_SEG:
                break
            s0 = seg_on + side * kk * (dur + mo.OFF_GAP)
            if any(abs(s0 - o) < dur for o in old):
                continue
            ok, _ = mo.segment_allowed(s0, dur, windows)
            if not ok:
                continue
            try:
                segs.append((s0, so.strain_segment(args, ev, e["meta"]["ifos"], st, band, s0)))
            except ValueError:
                pass
    return segs


def event_belt(e, own_grid, k):
    """Full Neyman belt of the event (v2.0 trials; same seed as the v2.0 off-source stage)."""
    G, Vb, _ = mc.event_trials(e, own_grid, 2000, 1, me.SEED + 11 + 1000 * k, which="A")
    D = mc.delta_matrix(0.5 * (G + G.T))
    return mc.neyman(Vb, D)


def q_curve(e, R, own_grid):
    chi = so.chi2_curve(e, "A", R, own_grid)
    return chi - chi.min(), int(np.argmin(chi))


def eob_polarizations(st, s):
    """SEOBNRv5PHM with the label settings: every setting naming IMRPhenomXPHM is replaced and
    XPHM-specific waveform arguments are dropped."""
    st2 = copy.deepcopy(st)                       # keeps numeric types (no JSON round trip)
    for k, v in list(st2.items()):
        if isinstance(v, str) and "IMRPhenomXPHM" in v:
            st2[k] = "SEOBNRv5PHM"
    if isinstance(st2.get("wf_args"), dict):
        st2["wf_args"] = {k: v for k, v in st2["wf_args"].items()
                          if not any(t in k for t in ("PhenomX", "SpinTaylor", "Prec", "MB"))}
    return me.rc.polarizations("C00:SEOBNRv5PHM", s, st2)


def pair_context(args, k, ev, label, r5):
    e = run.load(args, ev, label, r5, run.code_hash()[0])
    st, tab, idx, seg_on, band, network_h = m2.pair_setup(args, ev, label, r5, e)
    own_grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * e["sigma_lin"]
    return e, st, tab, idx, seg_on, band, network_h, own_grid


def run_check(which, hv, status):
    m2.install_patches()
    out = f"supplement_validation_{which}.csv"
    if os.path.exists(out):
        raise SystemExit(f"refused: {out} exists")
    keys = ["cohort", "event", "label", "status", "code_sha256_val", "seg_start", "f_index", "lambda_true",
            "q_true", "c_true", "covered", "estimate_sigma", "note"]
    fo = open(out, "w", newline="", encoding="utf-8")
    w = csv.DictWriter(fo, fieldnames=keys)
    w.writeheader()
    tot = [0, 0]
    skipped = []
    for coh in COHORTS:
        args = so.driver_args(coh)
        for k, (ev, label, r5) in enumerate(run.select(args)):
            t0 = time.time()
            e, st, tab, idx, seg_on, band, network_h, own_grid = pair_context(args, k, ev, label, r5)
            ifos, fb = e["meta"]["ifos"], e["fb"]
            F = set_F(idx, tab, e, network_h)
            segs = new_segments(args, ev, label, e, st, band, seg_on)
            c = event_belt(e, own_grid, k)
            tz = int(np.argmin(np.abs(own_grid)))
            calf = {i: e["sb"][i].mean_factor(fb) for i in ifos}
            f_full, _, _ = me.rc.polarizations(label, F[0][1], st) if F else (None, None, None)
            for m, (s0, dseg) in enumerate(segs):
                if m >= len(F):
                    break
                j, s, hF = F[m]
                note = ""
                if which == "v1":
                    sign = +1 if m % 2 == 0 else -1
                    t = tz + sign * int(round(T_SIGMA / args.step))
                    lam = own_grid[t]
                    hinj = {i: hF[i] * np.exp(1j * e["K"] * fb ** 3 * lam) for i in ifos}
                else:
                    t, lam = tz, 0.0
                    try:
                        f2, sp, sc = eob_polarizations(st, s)
                        if len(f2) != len(f_full) or not np.allclose(f2, f_full):
                            raise ValueError("frequency grid differs")
                        hinj = {i: me.rc.project(f2, sp, sc, i, s, seg_on)[band] * calf[i] for i in ifos}
                    except Exception as exc:
                        skipped.append((coh, ev, f"{type(exc).__name__}: {exc}"[:120]))
                        note = "SEOBNRv5PHM generation failed; skipped"
                        w.writerow({"cohort": coh, "event": ev, "label": label, "status": status,
                                    "code_sha256_val": hv, "seg_start": s0, "f_index": j, "lambda_true": lam,
                                    "q_true": "", "c_true": "", "covered": "", "estimate_sigma": "", "note": note})
                        continue
                d = {i: dseg[i] + hinj[i] for i in ifos}
                R, _, _ = so.residual(e, d)
                q, g = q_curve(e, R, own_grid)
                cov = bool(q[t] <= c[t])
                tot[0] += cov
                tot[1] += 1
                w.writerow({"cohort": coh, "event": ev, "label": label, "status": status, "code_sha256_val": hv,
                            "seg_start": s0, "f_index": j, "lambda_true": lam, "q_true": float(q[t]),
                            "c_true": float(c[t]), "covered": cov,
                            "estimate_sigma": float((own_grid[g] - lam) / e["sigma_lin"]), "note": note})
            fo.flush()
            print(f"  {coh:7s} {ev:18s} segments {len(segs)}  set F {len(F)}  ({time.time() - t0:.0f} s)", flush=True)
            del e
    fo.close()
    nc, n = tot
    lo, hi = mo.binomial_interval(0.90, n)
    verdict = ("within the 99% interval" if lo <= nc / n <= hi else
               "above the interval (conservative)" if nc / n > hi else "BELOW the 99% interval")
    rows = [r for r in csv.DictReader(open(out)) if r["covered"] != ""]
    est = np.array([float(r["estimate_sigma"]) for r in rows])
    print(f"\n{which.upper()}: pooled coverage {nc}/{n} = {nc / n:.3f}; 99% interval around 0.90: {lo:.3f}-{hi:.3f} "
          f"-> {verdict}")
    print(f"     estimate minus truth, in sigma_lin: median {np.median(est):+.3f}, 16-84% "
          f"[{np.percentile(est, 16):+.3f}, {np.percentile(est, 84):+.3f}]")
    if skipped:
        print(f"     skipped (SEOBNRv5PHM not generated): {len(skipped)} segments, events "
              f"{sorted(set(x[1] for x in skipped))}")
    print(f"written {out} (status {status})")


def stage_test(hv, status):
    m2.install_patches()
    args = so.driver_args("search")
    pairs = run.select(args)
    ev, label, r5 = pairs[0]
    e, st, tab, idx, seg_on, band, network_h, own_grid = pair_context(args, 0, ev, label, r5)
    row = next(r for r in csv.DictReader(open("v2_offsource_search_xphm.csv"))
               if r["event"] == ev and r["label"] == label and r["inj"] == "C")
    s = {kk: float(tab[kk][int(row["e_index"])]) for kk in me.NEED}
    hE = network_h(s)
    dseg = so.strain_segment(args, ev, e["meta"]["ifos"], st, band, float(row["seg_start"]))
    R, _, _ = so.residual(e, {i: dseg[i] + hE[i] for i in e["meta"]["ifos"]})
    q, _ = q_curve(e, R, own_grid)
    c = event_belt(e, own_grid, 0)
    tz = int(np.argmin(np.abs(own_grid)))
    dq = abs(q[tz] - float(row["q0"])) / max(1.0, abs(float(row["q0"])))
    dc = abs(c[tz] - float(row["c0"])) / max(1.0, abs(float(row["c0"])))
    print(f"reproduction ({ev}): q0 {q[tz]:.6f} vs {float(row['q0']):.6f}; c0 {c[tz]:.6f} vs {float(row['c0']):.6f} "
          f"-> {'PASS' if max(dq, dc) < 1e-6 else 'FAIL: do not freeze'}", flush=True)
    print(f"settings keys of {label}: {sorted(st.keys())}", flush=True)
    print("SEOBNRv5PHM probe (one set-F sample per pair):", flush=True)
    for coh in COHORTS:
        a = so.driver_args(coh)
        for k, (ev, label, r5) in enumerate(run.select(a)):
            e, st, tab, idx, seg_on, band, network_h, own_grid = pair_context(a, k, ev, label, r5)
            F = set_F(idx, tab, e, network_h)
            try:
                f0, _, _ = me.rc.polarizations(label, F[0][1], st)
                f2, sp, sc = eob_polarizations(st, F[0][1])
                ok = len(f2) == len(f0) and np.allclose(f2, f0)
                msg = "ok" if ok else "frequency grid differs"
            except Exception as exc:
                msg = f"FAILS: {type(exc).__name__}: {exc}"[:100]
            print(f"  {coh:7s} {ev:18s} set F {len(F)}  SEOBNRv5PHM: {msg}", flush=True)
            del e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["test", "v1", "v2", "hash"])
    ap.add_argument("--freeze-manifest", default=os.path.join(HERE, "FREEZE_supplement_validation.json"))
    ap.add_argument("--write-manifest", default="", metavar="TAG")
    cli = ap.parse_args()
    hv, hs, own = val_hash()
    status, why = val_status(hv, cli.freeze_manifest)
    print(f"supplement code SHA-256 {hs}\nsupplement_validation.py SHA-256 {own}\nvalidation code SHA-256 {hv}\n"
          f"status: {status} ({why})", flush=True)
    if cli.stage == "hash":
        if cli.write_manifest:
            if os.path.exists(cli.freeze_manifest):
                raise SystemExit(f"{cli.freeze_manifest} exists; refusing to overwrite")
            json.dump({"tag": cli.write_manifest, "code_sha256_val": hv, "code_sha256_supp": hs,
                       "supplement_validation_sha256": own,
                       "registration": "MODULE1_SUPPLEMENT_VALIDATION_REGISTRATION.md",
                       "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
                      open(cli.freeze_manifest, "w", encoding="utf-8"), indent=1)
            print(f"written {cli.freeze_manifest}")
        return
    if cli.stage == "test":
        stage_test(hv, status)
        return
    if status != "official-val":
        raise SystemExit("refused: the validation code must be frozen first")
    run_check(cli.stage, hv, status)


if __name__ == "__main__":
    main()
