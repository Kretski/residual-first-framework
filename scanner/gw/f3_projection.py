#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f3_projection.py — Module 1 feasibility check F3, version 3.

The quantity that matters for Module 1 is not each waveform derivative but the
ORTHOGONALIZED template T3⊥ = T3 − P_tangent T3 (MODULE1_DESIGN §6). Noisy
derivatives in directions with little weight or nearly orthogonal to T3 do not
matter. This script measures, per event and model label:

  survival s = |T3⊥| / |T3|               (fraction of the f³ template that the
                                           GR parameters cannot absorb)
  stability d = |T3⊥(ε) − T3⊥(ε/2)| / |T3⊥(ε/2)|  for several step multipliers

and, for the EOB models, the cross-check of a tangent basis built from
IMRPhenomXPHM derivatives at the EOB reference point.

Tangent basis: derivatives with respect to the 11 PE parameters (as in
f3_derivatives.py) plus the analytic coalescence-time direction −2πif·h.
Templates: T_p = i f^p h (p = 2, 3, 4); weighted inner product with the H1 PSD.
Criterion (proposed for MODULE1_DESIGN): d < 0.05 for the chosen step.

  python f3_projection.py --pe-dir ~/gwdata/pe --out ~/gw/f3
"""
import argparse
import csv
import os
import sys
import time

import h5py
import numpy as np

import f3_derivatives as f3

MULTIPLIERS = {"IMRPhenomXPHM": [0.25, 1.0, 4.0], "EOB": [10.0, 30.0, 100.0]}
D_CRIT = 0.05

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def real_stack(v, sw):
    """Complex vector → real vector so that the real dot product equals
    Re Σ a* b w (the noise-weighted inner product)."""
    return np.concatenate([v.real * sw, v.imag * sw])


def project_out(T, B):
    """T − P_B T with an SVD-truncated least-squares fit (real-stacked)."""
    U, S, Vt = np.linalg.svd(B, full_matrices=False)
    keep = S > S[0] * 1e-8
    U = U[:, keep]
    return T - U @ (U.T @ T)


def basis(gen, p, h0, mult, half, f, band, sw):
    cols = []
    for name in f3.STEPS:
        d = f3.derivative(gen, p, name, mult * f3.STEPS[name][0] / 2 ** half, h0)
        cols.append(real_stack(d[band], sw))
    cols.append(real_stack((-2j * np.pi * f * h0)[band], sw))     # coalescence time
    return np.column_stack(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--out", default=os.path.expanduser("~/gw/f3"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    files = os.listdir(args.pe_dir)
    rows = []
    for ev, (sub, labels) in f3.EVENTS.items():
        fn = next((x for x in files if sub in x), None)
        if fn is None:
            continue
        with h5py.File(os.path.join(args.pe_dir, fn), "r") as fh:
            xphm_label = labels[0]
            xs, xcfg, _ = f3.read_label(fh, xphm_label)
            xst = f3.settings(xcfg)
            for label in labels:
                t0 = time.time()
                sample, cfg, psd = f3.read_label(fh, label)
                st = f3.settings(cfg)
                approx = ("IMRPhenomXPHM" if "XPHM" in label else
                          "SEOBNRv4PHM" if "v4PHM" in label else "SEOBNRv5PHM")
                gen = f3.Generator(approx, st, set())
                p = {k: sample[k] for k in f3.STEPS}
                f, h0 = gen.fd(p)
                band = (f >= st["f_an"]) & (f <= st["f_high"])
                Sn = np.interp(f[band], psd[:, 0], psd[:, 1], left=np.inf, right=np.inf)
                sw = np.sqrt(1.0 / Sn)
                T = {q: real_stack((1j * f ** q * h0)[band], sw) for q in (2, 3, 4)}
                print(f"\n=== {ev}  {label} ===", flush=True)
                mults = MULTIPLIERS["IMRPhenomXPHM" if approx == "IMRPhenomXPHM" else "EOB"]
                best = None
                for m in mults:
                    try:
                        Tp = []
                        for half in (0, 1):
                            B = basis(gen, p, h0, m, half, f, band, sw)
                            Tp.append(project_out(T[3], B))
                    except Exception as e:
                        print(f"  step ×{m:g}: failed ({type(e).__name__}: {str(e)[:120]})")
                        continue
                    s = np.linalg.norm(Tp[1]) / np.linalg.norm(T[3])
                    d = np.linalg.norm(Tp[0] - Tp[1]) / max(np.linalg.norm(Tp[1]), 1e-300)
                    print(f"  step ×{m:<5g} survival s = {s:.4f}   stability d = {d:.3e}   "
                          f"{'OK' if d < D_CRIT else 'NOT OK'}", flush=True)
                    rows.append({"event": ev, "label": label, "basis": "own", "multiplier": m,
                                 "survival": s, "stability": d, "ok": d < D_CRIT})
                    if best is None or d < best[1]:
                        best = (m, d, Tp[1], s)
                # cross-check for EOB: tangent basis from XPHM derivatives at the EOB point
                if approx != "IMRPhenomXPHM" and best is not None:
                    xgen = f3.Generator("IMRPhenomXPHM", dict(st, wf_args=xst["wf_args"]), set())
                    try:
                        xf, xh0 = xgen.fd(p)
                        if len(xf) != len(f) or not np.allclose(xf, f):
                            xd_list = []
                            for half in (0, 1):
                                cols = []
                                for name in f3.STEPS:
                                    d_ = f3.derivative(xgen, p, name, f3.STEPS[name][0] / 2 ** half, xh0)
                                    cols.append(real_stack(np.interp(f[band], xf, d_.real)
                                                           + 1j * np.interp(f[band], xf, d_.imag), sw))
                                cols.append(real_stack((-2j * np.pi * f * h0)[band], sw))
                                xd_list.append(project_out(T[3], np.column_stack(cols)))
                        else:
                            xd_list = [project_out(T[3], basis(xgen, p, xh0, 1.0, half, f, band, sw))
                                       for half in (0, 1)]
                        sx = np.linalg.norm(xd_list[1]) / np.linalg.norm(T[3])
                        dx = np.linalg.norm(xd_list[0] - xd_list[1]) / np.linalg.norm(xd_list[1])
                        cross = np.linalg.norm(xd_list[1] - best[2]) / np.linalg.norm(best[2])
                        print(f"  XPHM-derivative basis at this point: survival {sx:.4f}, stability {dx:.3e}, "
                              f"difference to own best basis (×{best[0]:g}) {cross:.3e}")
                        rows.append({"event": ev, "label": label, "basis": "XPHM", "multiplier": 1.0,
                                     "survival": sx, "stability": dx, "ok": dx < D_CRIT,
                                     "cross_difference": cross})
                    except Exception as e:
                        print(f"  XPHM cross-check failed: {type(e).__name__}: {str(e)[:120]}")
                # survival of the control templates with the best own basis
                if best is not None:
                    for q in (2, 4):
                        B = basis(gen, p, h0, best[0], 1, f, band, sw)
                        sq = np.linalg.norm(project_out(T[q], B)) / np.linalg.norm(T[q])
                        print(f"  control template f^{q}: survival {sq:.4f}")
                print(f"  ({gen.count} waveforms, {time.time() - t0:.0f} s)", flush=True)
    keys = ["event", "label", "basis", "multiplier", "survival", "stability", "ok", "cross_difference"]
    with open(os.path.join(args.out, "f3_projection.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows([{k: r.get(k, "") for k in keys} for r in rows])
    print(f"\nwritten {os.path.join(args.out, 'f3_projection.csv')}")


if __name__ == "__main__":
    main()
