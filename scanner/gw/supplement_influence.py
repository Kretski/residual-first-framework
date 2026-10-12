#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
supplement_influence.py — DESCRIPTIVE drop-one influence table for the supplementary
IMRPhenomXPHM results (supplement_result_{search,confirm}.json). Not part of the registered
analysis: it was decided after the on-source results were committed, it changes no reported
set, and it is reported as a description of how much each event contributes.

For each cohort:
  - per event and subspace, the belt trials (module1_catalog.event_trials with the v2.0
    mismatch term, seeds SEED + 7 + 1000 k, 5000 trials) and the on-source chi2 curve on the
    cohort grid are computed ONCE, exactly as in supplement_onsource.py;
  - consistency gate: the set from all events must equal the stored on-source set, otherwise
    the script stops;
  - the set without event i uses the sums minus event i (same grid). For the event with the
    smallest sigma_lin, which fixes the grid, the stored leave-one-out set (own grid) is used.

  python supplement_influence.py --cohort search
  python supplement_influence.py --cohort confirm
"""
import argparse
import csv
import json
import os

import numpy as np

import supplement_onsource as so
import module2_run as m2                     # noqa: E402
import module1_run as run                    # noqa: E402
import module1_catalog as mc                 # noqa: E402
import module1_event as me                   # noqa: E402


def belt_c(Gs, Vbs):
    D = mc.delta_matrix(0.5 * (Gs + Gs.T))
    return mc.neyman(Vbs, D), D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", choices=["search", "confirm"], required=True)
    cli = ap.parse_args()
    res_file = f"supplement_result_{cli.cohort}.json"
    if not os.path.exists(res_file):
        raise SystemExit(f"refused: {res_file} does not exist (descriptive table only after the result)")
    stored = json.load(open(res_file))
    hs, h2, own = so.supp_hash()
    own_inf = run.sha256_file(os.path.abspath(__file__), normalise_eol=True)
    print(f"supplement code SHA-256 {hs}\nsupplement_influence.py SHA-256 {own_inf}\n"
          f"DESCRIPTIVE analysis (not registered; decided after the on-source results)", flush=True)
    if stored.get("code_sha256_supp") != hs:
        raise SystemExit("refused: the stored result was computed with different supplement code")

    m2.install_patches()
    args = so.driver_args(cli.cohort)
    pairs = run.select(args)
    metas = [json.load(open(run.cache_paths(args, ev, label)[1], encoding="utf-8")) for ev, label, _ in pairs]
    sig = np.array([m["sigma_lin"] for m in metas])
    smin = sig.min()
    grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * smin
    w = 1 / sig ** 2

    G = {s: [] for s in "AB"}
    V = {s: [] for s in "AB"}
    C = {s: [] for s in "AB"}
    for k, (ev, label, r5) in enumerate(pairs):
        e = run.load(args, ev, label, r5, run.code_hash()[0])
        for s in "AB":
            g, vb, _ = mc.event_trials(e, grid, args.belt_trials, 1, me.SEED + 7 + 1000 * k, which=s)
            G[s].append(g)
            V[s].append(vb)
        st, tab, idx, seg_on, band, _ = m2.pair_setup(args, ev, label, r5, e)
        d = so.strain_segment(args, ev, e["meta"]["ifos"], st, band, seg_on)
        R, _, _ = so.residual(e, d)
        for s in "AB":
            C[s].append(so.chi2_curve(e, s, R, grid))
        del e
        print(f"  {ev:18s} per-event belts and chi2 done", flush=True)

    Gt = {s: sum(G[s]) for s in "AB"}
    Vt = {s: sum(V[s]) for s in "AB"}
    Ct = {s: sum(C[s]) for s in "AB"}
    belts_all = {s: belt_c(Gt[s], Vt[s]) for s in "AB"}
    full = so.union_set(Ct, belts_all, grid)
    ref = stored["set"]
    if not (np.isclose(full["lower"], ref["lower"]) and np.isclose(full["upper"], ref["upper"])):
        raise SystemExit(f"STOP: recomputed full set {full['lower']:+.4e}..{full['upper']:+.4e} "
                         f"differs from the stored {ref['lower']:+.4e}..{ref['upper']:+.4e}")
    width0 = full["upper"] - full["lower"]
    print(f"\nconsistency: full set reproduced [{full['lower']:+.3e}, {full['upper']:+.3e}] m^2", flush=True)

    rows = []
    imin = int(np.argmin(sig))
    for i, (ev, label, _) in enumerate(pairs):
        if i == imin:
            s_set = stored["leave_one_out"]["set"]
            note = "own grid (stored leave-one-out)"
        else:
            belts = {s: belt_c(Gt[s] - G[s][i], Vt[s] - V[s][i]) for s in "AB"}
            s_set = so.union_set({s: Ct[s] - C[s][i] for s in "AB"}, belts, grid)
            note = ""
        if s_set is None:
            rows.append((ev, sig[i], w[i] / w.sum(), np.nan, np.nan, None, np.nan, "empty set"))
            continue
        wd = s_set["upper"] - s_set["lower"]
        rows.append((ev, sig[i], w[i] / w.sum(), s_set["lower"], s_set["upper"], s_set["contains_zero"],
                     100 * (wd / width0 - 1), note))

    out = f"supplement_influence_{cli.cohort}.csv"
    with open(out, "w", newline="", encoding="utf-8") as fo:
        wr = csv.writer(fo)
        wr.writerow(["event", "sigma_lin_m2", "weight_share", "lower_without_m2", "upper_without_m2",
                     "contains_zero", "width_change_percent", "note"])
        for r in rows:
            wr.writerow(r)
    print(f"\nDrop-one influence ({cli.cohort}, IMRPhenomXPHM, DESCRIPTIVE); full set width {width0:.3e} m^2")
    print(f"  {'event':18s} {'sigma_lin':>10s} {'share':>6s}   {'set without the event':>28s}  0 in  width")
    for r in sorted(rows, key=lambda r: -r[2]):
        if np.isnan(r[3]):
            print(f"  {r[0]:18s} {r[1]:10.2e} {r[2]:6.3f}   {'(empty)':>28s}")
        else:
            print(f"  {r[0]:18s} {r[1]:10.2e} {r[2]:6.3f}   [{r[3]:+.2e}, {r[4]:+.2e}]  {str(r[5]):5s} "
                  f"{r[6]:+6.1f}%  {r[7]}")
    print(f"written {out}")


if __name__ == "__main__":
    main()
