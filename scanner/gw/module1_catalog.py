#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_catalog.py — catalog combination and Neyman construction for Module 1
(MODULE1_DESIGN v0.9.10, §6 item 7 and §10). Synthetic noise only: no on-source
strain is read; the output is the EXPECTED interval and the coverage of the method.

Per event (geometry from module1_event.build, one model per run):
  ProfileScan on a COMMON absolute grid of Lambda [m^2]: b_g = C^-1 P s(Lambda_g),
  Gram matrix G_gh = (P s_g) . b_h, and for every noise trial v_g = b_g . (P n).
Because the noise enters linearly, the catalog needs only the SUMS over events:
  for data generated at Lambda_t,
    chi2_cat(g) - chi2_cat(t) = 2 (V_t - V_g) + Delta_t(g),
    V = sum_i v_i,  Delta_t(g) = sum_i (G_tt - 2 G_gt + G_gg)_i.
Test statistic (profile likelihood ratio, Feldman–Cousins ordering):
    q(t) = chi2_cat(t) - min_g chi2_cat(g) >= 0.
Neyman belt: c_t = 90% quantile of q(t) over belt trials generated at Lambda_t.
Interval for a data set: { t : q_obs(t) <= c_t }.
Reported:
  - expected interval under Lambda = 0 (median lower and upper ends over independent
    trials), compared with the LVK GWTC-4.0 bound;
  - coverage of the Neyman interval at several true values (independent trials), and
    of the Wilks interval (q <= 2.71) for comparison;
  - catalog joint-template matrix of f^2, f^3, f^4 (linear, sum over events): its
    condition number and the inflation of sigma(f^3) in a joint fit.

  python module1_catalog.py --events GW150914,GW200129_065458,GW230627_015337 \
      --labels xphm --n 800 --belt-trials 1000 --exp-trials 300
"""
import argparse
import csv
import json
import os
import sys
import time

import h5py
import numpy as np

import module1_event as me

LVK = (-2.4e-11, 7.4e-12)       # GWTC-4.0 combined, group velocity, 90%

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def cal_error_fn(ev, rng):
    net, sb, href, fb = ev["net"], ev["sb"], ev["href"], ev["fb"]

    def f():
        out = {}
        for i in net.ifos:
            pr = sb[i].pri
            out[i] = (sb[i].factor(fb, rng.normal(0, pr["sg_a"]), rng.normal(0, pr["sg_p"])) - 1.0) * href[i]
        return net.stack(out)
    return f


def event_trials(ev, grid, n_belt, n_exp, seed, which="A"):
    """Gram matrix and noise projections v (trials x grid) for one event and subspace.
    The SAME noise realisations are used for both subspaces (same seed), so that the A
    and B confidence sets of a trial come from the same data, as required by the union
    of v0.9.14."""
    A = ev[which]
    sc = me.ProfileScan(A, ev["K"], ev["href"], ev["net"], ev["fb"], grid)
    B = sc.b.astype(np.float64)
    G = sc.G
    del sc
    rng = np.random.default_rng(seed)
    cal = cal_error_fn(ev, rng)
    nvec = A.Q.shape[0]
    out = []
    for ntr in (n_belt, n_exp):
        V = np.empty((ntr, len(grid)))
        for j in range(ntr):
            n = rng.normal(size=nvec) + cal()
            npj = n - A.Q @ (A.Q.T @ n)
            V[j] = B @ npj
        out.append(V)
    return G, out[0], out[1]


def delta_matrix(Gsum):
    """Delta[t, g] = G_tt - 2 G_gt + G_gg (noiseless chi2 difference, truth t)."""
    d = np.diag(Gsum)
    return d[:, None] - 2 * Gsum + d[None, :]


def qstat(V, Delta, t):
    """q(t) for trials generated at truth t: -min_g [2 (V_t - V_g) + Delta_t(g)]."""
    D = 2 * (V[:, [t]] - V) + Delta[t][None, :]
    return -D.min(axis=1)


def q_obs_all(Vrow, Delta, t0):
    """q_obs(t) for every t, for one data set generated at truth t0."""
    E = 2 * (Vrow[t0] - Vrow) + Delta[t0]          # chi2(g) - chi2(t0)
    return E - E.min()


def neyman(V_belt, Delta, level=0.90):
    return np.array([np.quantile(qstat(V_belt, Delta, t), level) for t in range(Delta.shape[0])])


def union_interval(qs, cs, grid):
    """Union of the confidence sets of the subspaces (v0.9.14, §7)."""
    inside = np.zeros(len(grid), bool)
    for q, c in zip(qs, cs):
        inside |= q <= c
    idx = np.where(inside)[0]
    if len(idx) == 0:
        return np.nan, np.nan, True
    return grid[idx[0]], grid[idx[-1]], (idx[-1] - idx[0] + 1) == len(idx)


def interval(q, c, grid):
    inside = q <= c
    idx = np.where(inside)[0]
    if len(idx) == 0:
        return np.nan, np.nan, True
    contiguous = (idx[-1] - idx[0] + 1) == len(idx)
    return grid[idx[0]], grid[idx[-1]], contiguous


def catalog_fisher(events):
    F = np.zeros((3, 3))
    for ev in events:
        A = ev["A"]
        for a, p in enumerate((2, 3, 4)):
            for b, q in enumerate((2, 3, 4)):
                F[a, b] += float(A.Tp[p] @ A.CTp[q])
    F = 0.5 * (F + F.T)
    R = F / np.sqrt(np.outer(np.diag(F), np.diag(F)))
    w = np.linalg.eigvalsh(R)
    infl = np.sqrt(np.linalg.inv(R)[1, 1])
    return R, w[-1] / w[0], infl


def analyse(events, grid, n_belt, n_exp, seed, log=print):
    ng = len(grid)
    Gsum = np.zeros((ng, ng))
    Vb = np.zeros((n_belt, ng))
    Ve = np.zeros((n_exp, ng))
    for k, ev in enumerate(events):
        t = time.time()
        G, vb, ve = event_trials(ev, grid, n_belt, n_exp, seed + 1000 * k)
        Gsum += G
        Vb += vb
        Ve += ve
        log(f"  trials for {ev['name']} ({ev['label']}): {time.time() - t:.0f} s")
    Delta = delta_matrix(0.5 * (Gsum + Gsum.T))
    c = neyman(Vb, Delta)
    t0 = int(np.argmin(np.abs(grid)))
    res = {"grid_min": float(grid[0]), "grid_max": float(grid[-1]), "grid_step": float(grid[1] - grid[0])}
    lo, hi, cont = [], [], 0
    for j in range(n_exp):
        q = q_obs_all(Ve[j], Delta, t0)
        a, b, ok = interval(q, c, grid)
        lo.append(a)
        hi.append(b)
        cont += ok
    res["exp_lower_median"] = float(np.nanmedian(lo))
    res["exp_upper_median"] = float(np.nanmedian(hi))
    res["exp_contiguous_fraction"] = cont / n_exp
    res["exp_at_grid_edge_fraction"] = float(np.mean((np.array(lo) <= grid[0]) | (np.array(hi) >= grid[-1])))
    cov = {}
    for frac in (-0.5, -0.2, 0.0, 0.2, 0.5):
        t = int(np.argmin(np.abs(grid - frac * grid[-1])))
        qn = np.array([q_obs_all(Ve[j], Delta, t)[t] for j in range(n_exp)])
        cov[f"{grid[t]:+.3e}"] = {"neyman": float(np.mean(qn <= c[t])), "wilks": float(np.mean(qn <= 2.71))}
    res["coverage"] = cov
    res["critical_values_range"] = [float(c.min()), float(c.max())]
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="GW150914,GW200129_065458,GW230627_015337")
    ap.add_argument("--labels", default="xphm", choices=["xphm", "eob"])
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    ap.add_argument("--n", type=int, default=800)
    ap.add_argument("--belt-trials", type=int, default=1000)
    ap.add_argument("--exp-trials", type=int, default=300)
    ap.add_argument("--half-width", type=float, default=8.0, help="grid half-width in units of the smallest sigma_lin")
    ap.add_argument("--step", type=float, default=0.05, help="grid step in units of the smallest sigma_lin")
    ap.add_argument("--out", default="module1_catalog_check.json")
    args = ap.parse_args()
    rows5 = {r["commonName"]: r for r in csv.DictReader(open(args.f5, encoding="utf-8"))}
    events = []
    for ev in args.events.split(","):
        r5 = rows5[ev]
        ifos = [i for i in r5["ifos"].split(",") if i]
        label = r5["xphm_label"] if args.labels == "xphm" else r5["eob_label"]
        print(f"\n=== {ev}  {label}", flush=True)
        with h5py.File(os.path.join(args.pe_dir, r5["file"].split(" ")[0]), "r") as fh:
            est, ks, info, sb, href, net, fb = me.build(ev, label, r5["xphm_label"], fh, ifos, args.n,
                                                        me.FRAC, lambda m: print(m, flush=True))
        A = est["A"]
        print(f"  z median {info['z_med']:.4f}  k {ks['A']}  s3 {A.s[3]:.4f}  sigma_lin {A.sigma(3):.3e} m^2", flush=True)
        events.append({"name": ev, "label": label, "A": A, "K": info["K"], "href": href, "net": net,
                       "fb": fb, "sb": sb, "sigma_lin": A.sigma(3)})
        del est["B"]
    smin = min(e["sigma_lin"] for e in events)
    grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * smin
    R, cond, infl = catalog_fisher(events)
    sig_cat = 1.0 / np.sqrt(sum(1 / e["sigma_lin"] ** 2 for e in events))
    print(f"\ncatalog: {len(events)} events; linear sigma_cat {sig_cat:.3e} m^2; grid ±{args.half_width:g} x "
          f"{smin:.3e} m^2, {len(grid)} points")
    print(f"catalog joint f2/f3/f4: rho23 {R[0, 1]:+.4f} rho34 {R[1, 2]:+.4f} rho24 {R[0, 2]:+.4f}; "
          f"condition {cond:.0f}; sigma(f3) inflation in a joint fit x{infl:.1f}", flush=True)
    t = time.time()
    res = analyse(events, grid, args.belt_trials, args.exp_trials, me.SEED + 7)
    print(f"\nNeyman construction ({args.belt_trials} belt trials, {args.exp_trials} independent trials; "
          f"{time.time() - t:.0f} s); critical values {res['critical_values_range'][0]:.2f}–"
          f"{res['critical_values_range'][1]:.2f} (Wilks would be 2.71)")
    print(f"expected 90% interval under Lambda = 0 (median ends): [{res['exp_lower_median']:+.2e}, "
          f"{res['exp_upper_median']:+.2e}] m^2   (LVK GWTC-4.0: [{LVK[0]:+.1e}, {LVK[1]:+.1e}])")
    print(f"  contiguous intervals {res['exp_contiguous_fraction']:.2f}; touching the grid edge "
          f"{res['exp_at_grid_edge_fraction']:.2f}")
    print("coverage at true Lambda:   Neyman   Wilks")
    for k, v in res["coverage"].items():
        print(f"   {k} m^2      {v['neyman']:.3f}    {v['wilks']:.3f}")
    res.update({"events": [(e["name"], e["label"], e["sigma_lin"]) for e in events], "sigma_cat_linear": sig_cat,
                "joint_rho": R.tolist(), "joint_condition": cond, "joint_inflation_f3": infl,
                "belt_trials": args.belt_trials, "exp_trials": args.exp_trials, "n": args.n})
    json.dump(res, open(args.out, "w"), indent=1)
    print(f"\nwritten {args.out}")


if __name__ == "__main__":
    main()
