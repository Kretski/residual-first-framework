#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
supplement_onsource.py — supplementary IMRPhenomXPHM-only analysis of Module 1
(MODULE1_SUPPLEMENT_XPHM_REGISTRATION.md). Secondary analysis, selected after the EOB failures
of v1.0 and v2.0. Uses the v2.0 construction exactly as frozen (module2_run.py, tag
v2.0-module1); the frozen files are imported unchanged.

Stages:
  test      the on-source chain applied to an OFF-source segment with its set-E injection, for
            the first --n-test events of the cohort; q(0) and c(0) must reproduce the values
            stored by the v2.0 off-source stage (otherwise the code must not be frozen)
  gate      the stopping rules for IMRPhenomXPHM alone (registration §3):
              rule 0  catalog coverage of the v2.0 XPHM catalogs (lower end of the 99% interval
                      around 0.90 for the number of check trials) and grid edge <= 1%;
              rule 0a pooled off-source coverage of the v2.0 XPHM off-source stage;
              rule 3  calibration-only, noise-free injections (100 draws from the calibration
                      priors): fails if |median catalog estimate| > 1.0 sigma_cat or if more
                      than 10% of the draws give a confidence set that excludes Lambda = 0;
              rules 1, 1a, 2a are pair-level and already applied (Addendum D lists).
  onsource  the on-source confidence set of ONE cohort (--cohort search|confirm), only with
            --unblind, only with the frozen code (FREEZE_supplement_xphm.json), only after a
            passing gate file, and for the confirmation cohort only after the search result
            file exists. The v2.0 catalog critical values are recomputed with the frozen seeds
            and must reproduce the stored range, otherwise the run stops.
  hash      code hash; --write-manifest TAG writes FREEZE_supplement_xphm.json

  python supplement_onsource.py test --cohort search --n-test 3
  python supplement_onsource.py gate
  python supplement_onsource.py onsource --cohort search --unblind
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time

import numpy as np

import module2_run as m2                     # installs the v1.0 on-source safeguard
import module1_run as run                    # noqa: E402
import module1_catalog as mc                 # noqa: E402
import module1_event as me                   # noqa: E402
import module1_offsource as mo               # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
COHORT = {"search": "primary_search_xphm.txt", "confirm": "primary_confirm_xphm.txt"}
V2CAT = {"search": "v2_catalog_search_xphm.json", "confirm": "v2_catalog_confirm_xphm.json"}
V2OFF = {"search": "v2_offsource_search_xphm.csv", "confirm": "v2_offsource_confirm_xphm.csv"}
LVK = (-2.4e-11, 7.4e-12)
N_CAL, CAL_SIGMA_MAX, CAL_EXCL_MAX = 100, 1.0, 0.10


# ------------------------------------------------------------------ hashes / status ----
def supp_hash():
    h2, h1, own2 = m2.v2_hash()
    own = run.sha256_file(os.path.join(HERE, "supplement_onsource.py"), normalise_eol=True)
    return hashlib.sha256((h2 + own).encode()).hexdigest(), h2, own


def supp_status(hs, manifest):
    if not os.path.exists(manifest):
        return "trial-supp", "no supplement freeze manifest"
    m = json.load(open(manifest, encoding="utf-8"))
    if m.get("code_sha256_supp") != hs:
        return "trial-supp", f"code hash differs from {m.get('tag', '?')}"
    return "official-supp", f"code hash matches {m.get('tag', '?')}"


# ------------------------------------------------------------------ common pieces ------
def driver_args(cohort):
    return argparse.Namespace(
        events=open(COHORT[cohort]).read().strip(), labels="xphm",
        pe_dir=os.path.expanduser("~/gwdata/pe"), f5="f5/f5_results.csv",
        residual_check="residual_check.csv", cache=os.path.expanduser("~/gw/cache_module1"),
        cache2=os.path.expanduser("~/gw/cache_module2"), n=800,
        strain_dir=os.path.expanduser("~/gwdata/strain"), events_list="event_selection/event_list_v1.csv",
        half_width=8.0, step=0.05, belt_trials=5000, exp_trials=2000, out="")


def strain_segment(args, ev, ifos, st, band, s0):
    """Whitening-free rfft of the strain segment [s0, s0 + duration], as the off-source stage."""
    from scipy.signal import resample_poly
    from scipy.signal.windows import tukey
    manifest = {(r["event"], r["ifo"]): r for r in
                csv.DictReader(open(os.path.join(args.strain_dir, "strain_manifest.csv"), encoding="utf-8"))}
    dur, srate = st["duration"], st["srate"]
    win = tukey(int(round(dur * srate)), alpha=2 * me.rc.ROLL_OFF / dur)
    out = {}
    for i in ifos:
        x, dt, t0f = run._ORIG_READ_STRAIN(os.path.join(args.strain_dir, manifest[(ev, i)]["file"]))
        if abs(1 / dt - srate) > 1e-6:
            x = resample_poly(x, int(srate), int(round(1 / dt)))
            dt = 1.0 / srate
        out[i] = (np.fft.rfft(run._ORIG_SEGMENT(x, dt, t0f, s0, dur) * win) * dt)[band]
    return out


def residual(e, d):
    ifos, fb = e["meta"]["ifos"], e["fb"]
    Seff = {i: 4 * (fb[1] - fb[0]) / e["net"].sw[i] ** 2 for i in ifos}
    r, tsh, a = mo.refine(d, e["href"], Seff, fb[1] - fb[0], fb=fb)
    return e["net"].stack(r), tsh, a


def chi2_curve(e, wsp, R, grid):
    """chi2(g) up to a constant for a stacked residual R, subspace wsp, grid."""
    A = e[wsp]
    sc = me.ProfileScan(A, e["K"], e["href"], e["net"], e["fb"], grid)
    pr = R - A.Q @ (A.Q.T @ R)
    v = sc.b.astype(np.float64) @ pr
    return -2 * v + np.diag(sc.G)


def cohort_belts(args, pairs, grid, keep=None):
    """Critical values per subspace, exactly as module1_run.stage_catalog (seeds SEED+7+1000k)."""
    keep = range(len(pairs)) if keep is None else keep
    out = {}
    for wsp in ("A", "B"):
        Gs, Vbs = [], []
        for k in keep:
            ev, label, r5 = pairs[k]
            e = run.load(args, ev, label, r5, run.code_hash()[0])
            G, vb, _ = mc.event_trials(e, grid, args.belt_trials, 1, me.SEED + 7 + 1000 * k, which=wsp)
            Gs.append(G)
            Vbs.append(vb)
            del e
        D = mc.delta_matrix(0.5 * (sum(Gs) + sum(Gs).T))
        out[wsp] = (mc.neyman(sum(Vbs), D), D)
    return out


def union_set(chis, belts, grid):
    inside = np.zeros(len(grid), bool)
    qs = {}
    for wsp in ("A", "B"):
        q = chis[wsp] - chis[wsp].min()
        qs[wsp] = q
        inside |= q <= belts[wsp][0]
    idx = np.where(inside)[0]
    if len(idx) == 0:
        return None
    return {"lower": float(grid[idx[0]]), "upper": float(grid[idx[-1]]),
            "contiguous": bool(idx[-1] - idx[0] + 1 == len(idx)),
            "contains_zero": bool(inside[int(np.argmin(np.abs(grid)))]),
            "touches_edge": bool(idx[0] == 0 or idx[-1] == len(grid) - 1),
            "best_A": float(grid[int(np.argmin(chis["A"]))]), "best_B": float(grid[int(np.argmin(chis["B"]))])}


# ------------------------------------------------------------------ stages -------------
def stage_test(cli, hs, status):
    m2.install_patches()
    args = driver_args(cli.cohort)
    args.belt_trials = 2000
    pairs = run.select(args)
    stored = [r for r in csv.DictReader(open(V2OFF[cli.cohort])) if r["inj"] == "C"]
    worst = 0.0
    for k, (ev, label, r5) in enumerate(pairs[:cli.n_test]):
        e = run.load(args, ev, label, r5, run.code_hash()[0])
        st, tab, idx, seg_on, band, network_h = m2.pair_setup(args, ev, label, r5, e)
        row = next(r for r in stored if r["event"] == ev and r["label"] == label)
        s0, j = float(row["seg_start"]), int(row["e_index"])
        s = {kk: float(tab[kk][j]) for kk in me.NEED}
        hE = network_h(s)
        dseg = strain_segment(args, ev, e["meta"]["ifos"], st, band, s0)
        d = {i: dseg[i] + hE[i] for i in e["meta"]["ifos"]}
        R, _, _ = residual(e, d)
        own_grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * e["sigma_lin"]
        c0, scan, tz = mo.event_c0(e, own_grid, args.belt_trials, me.SEED + 11 + 1000 * k)
        q0, _ = mo.q0_of_residual(R, e, scan, tz)
        dq = abs(q0 - float(row["q0"])) / max(1.0, abs(float(row["q0"])))
        dc = abs(c0 - float(row["c0"])) / max(1.0, abs(float(row["c0"])))
        worst = max(worst, dq, dc)
        print(f"  {ev:18s} segment {s0:.1f}: q0 {q0:.6f} vs stored {float(row['q0']):.6f}; "
              f"c0 {c0:.6f} vs stored {float(row['c0']):.6f}", flush=True)
    ok = worst < 1e-6
    print(f"\nreproduction of the v2.0 off-source chain: worst relative difference {worst:.2e} -> "
          f"{'PASS' if ok else 'FAIL: do not freeze'}  (status {status})")


def stage_gate(cli, hs, status):
    res = {"status": status, "code_sha256_supp": hs, "rules": {}}
    fail = False
    # rule 0
    for coh in ("search", "confirm"):
        r = json.load(open(V2CAT[coh]))["all"]
        n = json.load(open(V2CAT[coh]))["exp_trials"]
        lo, _ = mo.binomial_interval(0.90, n)
        cmin = min(v["neyman"] for v in r["coverage"].values())
        ok = cmin >= lo and r["edge"] <= 0.01
        fail |= not ok
        res["rules"][f"rule0_{coh}"] = {"min_coverage": cmin, "lower_99": lo, "edge": r["edge"], "pass": ok}
        print(f"rule 0  {coh:8s}: min Neyman coverage {cmin:.3f} (>= {lo:.3f}), edge {r['edge']:.3f} -> "
              f"{'pass' if ok else 'FAIL'}")
    # rule 0a
    nc = n = 0
    for coh in ("search", "confirm"):
        for r in csv.DictReader(open(V2OFF[coh])):
            if r["inj"] == "C":
                nc += r["covered"] in ("True", "1")
                n += 1
    lo, _ = mo.binomial_interval(0.90, n)
    ok = nc / n >= lo
    fail |= not ok
    res["rules"]["rule0a"] = {"covered": nc, "n": n, "lower_99": lo, "pass": ok}
    print(f"rule 0a pooled   : {nc}/{n} = {nc / n:.3f} (>= {lo:.3f}) -> {'pass' if ok else 'FAIL'}")
    # rule 3: calibration-only, noise-free injections, per cohort
    m2.install_patches()
    for coh in ("search", "confirm"):
        args = driver_args(coh)
        pairs = run.select(args)
        metas = [json.load(open(run.cache_paths(args, ev, label)[1], encoding="utf-8")) for ev, label, _ in pairs]
        smin = min(m["sigma_lin"] for m in metas)
        sig_cat = 1.0 / np.sqrt(sum(1 / m["sigma_lin"] ** 2 for m in metas))
        grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * smin
        belts = cohort_belts(args, pairs, grid)
        est, excl = [], 0
        chis = [{"A": np.zeros(len(grid)), "B": np.zeros(len(grid))} for _ in range(N_CAL)]
        for k, (ev, label, r5) in enumerate(pairs):
            e = run.load(args, ev, label, r5, run.code_hash()[0])
            rng = np.random.default_rng(me.SEED + 31 + 1000 * k)
            cal = mc.cal_error_fn(e, rng)
            draws = [cal() for _ in range(N_CAL)]
            for wsp in ("A", "B"):
                A = e[wsp]
                sc = me.ProfileScan(A, e["K"], e["href"], e["net"], e["fb"], grid)
                B = sc.b.astype(np.float64)
                dG = np.diag(sc.G)
                for j, nvec in enumerate(draws):
                    pr = nvec - A.Q @ (A.Q.T @ nvec)
                    chis[j][wsp] += -2 * (B @ pr) + dG
            del e
        for j in range(N_CAL):
            u = union_set(chis[j], belts, grid)
            est.append(grid[int(np.argmin(chis[j]["A"]))])
            excl += (u is None) or (not u["contains_zero"])
        med = float(np.median(est))
        ok = abs(med) <= CAL_SIGMA_MAX * sig_cat and excl / N_CAL <= CAL_EXCL_MAX
        fail |= not ok
        res["rules"][f"rule3_{coh}"] = {"median_estimate": med, "sigma_cat": sig_cat,
                                       "fraction_excluding_zero": excl / N_CAL, "pass": ok}
        print(f"rule 3  {coh:8s}: median estimate {med:+.2e} m^2 ({med / sig_cat:+.2f} sigma_cat), "
              f"sets excluding 0: {excl}/{N_CAL} -> {'pass' if ok else 'FAIL'}")
    res["all_pass"] = not fail
    json.dump(res, open("supplement_gate.json", "w"), indent=1)
    print(f"\ngate: {'ALL RULES PASS' if not fail else 'A RULE FIRED: no on-source computation'}; "
          f"written supplement_gate.json (status {status})")


def stage_onsource(cli, hs, status):
    if not cli.unblind:
        raise SystemExit("refused: --unblind is required")
    if status != "official-supp":
        raise SystemExit(f"refused: status {status}; the code must be frozen first")
    if not os.path.exists("supplement_gate.json") or not json.load(open("supplement_gate.json"))["all_pass"]:
        raise SystemExit("refused: no passing supplement_gate.json")
    if json.load(open("supplement_gate.json")).get("code_sha256_supp") != hs:
        raise SystemExit("refused: the gate was evaluated with different code")
    if cli.cohort == "confirm" and not os.path.exists("supplement_result_search.json"):
        raise SystemExit("refused: the search result must exist (and be committed) first")
    out = f"supplement_result_{cli.cohort}.json"
    if os.path.exists(out):
        raise SystemExit(f"refused: {out} exists; an on-source result is computed once")
    m2.install_patches()
    args = driver_args(cli.cohort)
    pairs = run.select(args)
    metas = [json.load(open(run.cache_paths(args, ev, label)[1], encoding="utf-8")) for ev, label, _ in pairs]
    smin = min(m["sigma_lin"] for m in metas)
    grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * smin
    print(f"cohort {cli.cohort}: {len(pairs)} pairs; grid ±{args.half_width:g} x {smin:.3e} m^2", flush=True)
    belts = cohort_belts(args, pairs, grid)
    crange = [min(belts[w][0].min() for w in "AB"), max(belts[w][0].max() for w in "AB")]
    stored = json.load(open(V2CAT[cli.cohort]))["all"]["critical_values"]
    if not np.allclose(crange, stored, rtol=1e-9, atol=1e-12):
        raise SystemExit(f"STOP: critical values {crange} do not reproduce the v2.0 catalog {stored}")
    print(f"critical values reproduced: {crange[0]:.4f}–{crange[1]:.4f}", flush=True)
    chis_e = []
    for k, (ev, label, r5) in enumerate(pairs):
        e = run.load(args, ev, label, r5, run.code_hash()[0])
        st, tab, idx, seg_on, band, _ = m2.pair_setup(args, ev, label, r5, e)
        d = strain_segment(args, ev, e["meta"]["ifos"], st, band, seg_on)
        R, tsh, a = residual(e, d)
        chis_e.append({w: chi2_curve(e, w, R, grid) for w in "AB"})
        print(f"  {ev:18s} on-source residual: time shift {tsh * 1e3:+.2f} ms, |a| {abs(a):.3f}", flush=True)
        del e
    tot = {w: sum(c[w] for c in chis_e) for w in "AB"}
    res = {"cohort": cli.cohort, "status": status, "code_sha256_supp": hs, "pairs": [p[:2] for p in pairs],
           "grid": [float(grid[0]), float(grid[-1]), float(grid[1] - grid[0])],
           "critical_values": crange, "set": union_set(tot, belts, grid), "LVK_GWTC4": LVK}
    w = np.array([1 / m["sigma_lin"] ** 2 for m in metas])
    top = int(np.argmax(w))
    keep = [i for i in range(len(pairs)) if i != top]
    smin_loo = min(metas[i]["sigma_lin"] for i in keep)
    loo_grid = grid if smin_loo <= smin * (1 + 1e-9) else \
        np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * smin_loo
    if loo_grid is grid:
        belts_loo = cohort_belts(args, pairs, grid, keep)
        tot_loo = {w2: sum(chis_e[i][w2] for i in keep) for w2 in "AB"}
    else:
        belts_loo = cohort_belts(args, pairs, loo_grid, keep)
        tot_loo = {w2: np.zeros(len(loo_grid)) for w2 in "AB"}
        for i in keep:
            ev, label, r5 = pairs[i]
            e = run.load(args, ev, label, r5, run.code_hash()[0])
            st, tab, idx, seg_on, band, _ = m2.pair_setup(args, ev, label, r5, e)
            d = strain_segment(args, ev, e["meta"]["ifos"], st, band, seg_on)
            R, _, _ = residual(e, d)
            for w2 in "AB":
                tot_loo[w2] += chi2_curve(e, w2, R, loo_grid)
            del e
    res["leave_one_out"] = {"left_out": metas[top]["event"], "set": union_set(tot_loo, belts_loo, loo_grid)}
    json.dump(res, open(out, "w"), indent=1)
    s = res["set"]
    print(f"\n90% confidence set ({cli.cohort}, IMRPhenomXPHM, v2.0 construction, SECONDARY analysis):")
    if s is None:
        print("  empty set")
    else:
        print(f"  [{s['lower']:+.3e}, {s['upper']:+.3e}] m^2  contiguous {s['contiguous']}  "
              f"contains 0: {s['contains_zero']}  touches edge: {s['touches_edge']}")
    sl = res["leave_one_out"]["set"]
    if sl:
        print(f"  leave-one-out without {res['leave_one_out']['left_out']}: [{sl['lower']:+.3e}, {sl['upper']:+.3e}] m^2"
              f"  contains 0: {sl['contains_zero']}")
    print(f"  LVK GWTC-4.0 (90%): [{LVK[0]:+.1e}, {LVK[1]:+.1e}] m^2")
    print(f"written {out} (status {status}); commit it before any further on-source computation")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["test", "gate", "onsource", "hash"])
    ap.add_argument("--cohort", choices=["search", "confirm"], default="search")
    ap.add_argument("--n-test", type=int, default=3)
    ap.add_argument("--unblind", action="store_true")
    ap.add_argument("--freeze-manifest", default=os.path.join(HERE, "FREEZE_supplement_xphm.json"))
    ap.add_argument("--write-manifest", default="", metavar="TAG")
    cli = ap.parse_args()
    hs, h2, own = supp_hash()
    status, why = supp_status(hs, cli.freeze_manifest)
    print(f"v2.0 code SHA-256 {h2}\nsupplement_onsource.py SHA-256 {own}\nsupplement code SHA-256 {hs}\n"
          f"status: {status} ({why})", flush=True)
    if cli.stage == "hash":
        if cli.write_manifest:
            if os.path.exists(cli.freeze_manifest):
                raise SystemExit(f"{cli.freeze_manifest} exists; refusing to overwrite")
            json.dump({"tag": cli.write_manifest, "code_sha256_supp": hs, "code_sha256_v2": h2,
                       "supplement_onsource_sha256": own,
                       "registration": "MODULE1_SUPPLEMENT_XPHM_REGISTRATION.md",
                       "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
                      open(cli.freeze_manifest, "w", encoding="utf-8"), indent=1)
            print(f"written {cli.freeze_manifest}")
        return
    {"test": stage_test, "gate": stage_gate, "onsource": stage_onsource}[cli.stage](cli, hs, status)


if __name__ == "__main__":
    main()
