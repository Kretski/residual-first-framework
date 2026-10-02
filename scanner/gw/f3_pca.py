#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f3_pca.py — Module 1 feasibility check F3, version 4.

Finite-difference derivatives gave a nearly degenerate, noisy tangent basis
(f3_projection.py: stable survival, unstable direction of T3⊥). Here the subspace
that the GR fit can absorb is built from the posterior itself:

  1. reference waveform h_ref at the maximum-likelihood sample;
  2. two DISJOINT random sets A and B of N posterior samples (fixed seed);
  3. Δ_i = h(θ_i) − h_ref, whitened with the released H1 PSD (SNR units:
     4 Δf Σ |·|² / S_n);
  4. PCA (SVD) of the Δ_i; the kept components plus the analytic
     coalescence-time direction −2πif·h_ref form the tangent subspace;
  5. T3⊥ = T3 − P T3 for T_p = i f^p h_ref.

Kept components are tested with several pre-declared rules:
  frac≥0.99 / 0.999 / 0.9999 of the explained variance, and
  rms≥0.1 (components whose RMS amplitude over the samples is at least 0.1
  noise units — directions the posterior actually explores).
Stability: d = |T3⊥(A) − T3⊥(B)| / |T3⊥(B)|, cosine overlap c between T3⊥(A)
and T3⊥(B), and relative survival difference |sA − sB| / sB.
Note: for a pure f³ signal r = Λ T3 the estimate ⟨r, T3⊥⟩/⟨T3⊥, T3⊥⟩ = Λ for ANY
subspace (⟨T3, T3⊥⟩ = |T3⊥|²); the subspace sets the estimator noise (via the
survival s) and the leakage of reference-point imperfections into T3⊥ (via the
direction, i.e. c). The rule to be frozen is chosen from these numbers.

  python f3_pca.py --pe-dir ~/gwdata/pe --out ~/gw/f3 --n 200
"""
import argparse
import csv
import os
import sys
import time

import h5py
import numpy as np

import f3_derivatives as f3

RULES = [("frac", 0.99), ("frac", 0.999), ("frac", 0.9999), ("rms", 0.1)]
D_CRIT = 0.05
SEED = 20261001

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def posterior_table(fh, label):
    ps = fh[label]["posterior_samples"]
    if isinstance(ps, h5py.Dataset) and ps.dtype.names:
        arr = ps[()]
        return {n: np.asarray(arr[n]) for n in arr.dtype.names}
    names = [f3.text(v) for v in ps["parameter_names"][()]]
    smp = ps["samples"][()]
    return {n: smp[:, i] for i, n in enumerate(names)}


def stack(v, sw):
    return np.concatenate([v.real * sw, v.imag * sw])


def subspace(D, rule, value):
    """Orthonormal basis (columns) of the kept principal components of D."""
    U, S, _ = np.linalg.svd(D, full_matrices=False)
    if rule == "frac":
        cum = np.cumsum(S ** 2) / np.sum(S ** 2)
        k = int(np.searchsorted(cum, value) + 1)
    else:                                   # RMS amplitude over samples
        rms = S / np.sqrt(D.shape[1])
        k = int(np.sum(rms >= value))
    return U[:, :max(k, 1)], k


def project_out(T, Q):
    return T - Q @ (Q.T @ T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--out", default=os.path.expanduser("~/gw/f3"))
    ap.add_argument("--n", type=int, default=200, help="samples per set (two disjoint sets)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    files = os.listdir(args.pe_dir)
    rows = []
    for ev, (sub, labels) in f3.EVENTS.items():
        fn = next((x for x in files if sub in x), None)
        if fn is None:
            continue
        with h5py.File(os.path.join(args.pe_dir, fn), "r") as fh:
            for label in labels:
                t0 = time.time()
                tab = posterior_table(fh, label)
                ref_sample, cfg, psd = f3.read_label(fh, label)
                st = f3.settings(cfg)
                approx = ("IMRPhenomXPHM" if "XPHM" in label else
                          "SEOBNRv4PHM" if "v4PHM" in label else "SEOBNRv5PHM")
                gen = f3.Generator(approx, st, set())
                p_ref = {k: ref_sample[k] for k in f3.STEPS}
                f, h_ref = gen.fd(p_ref)
                df = f[1] - f[0]
                band = (f >= st["f_an"]) & (f <= st["f_high"])
                Sn = np.interp(f[band], psd[:, 0], psd[:, 1], left=np.inf, right=np.inf)
                sw = np.sqrt(4 * df / Sn)                          # SNR units
                nsamp = len(tab["log_likelihood"])
                rng = np.random.default_rng(SEED)
                idx = rng.permutation(nsamp)
                need = 2 * args.n
                print(f"\n=== {ev}  {label}  ({nsamp} samples; using 2 × {args.n}) ===", flush=True)
                sets, fails = {"A": [], "B": []}, 0
                pos = 0
                for name in ("A", "B"):
                    while len(sets[name]) < args.n and pos < nsamp:
                        i = idx[pos]
                        pos += 1
                        p = {k: float(tab[k][i]) for k in f3.STEPS}
                        try:
                            _, h = gen.fd(p)
                        except Exception:
                            fails += 1
                            continue
                        sets[name].append(stack((h - h_ref)[band], sw))
                if min(len(sets["A"]), len(sets["B"])) < args.n:
                    print(f"  not enough usable samples (failures {fails}); skipped")
                    continue
                tdir = stack((-2j * np.pi * f * h_ref)[band], sw)
                T = {q: stack((1j * f ** q * h_ref)[band], sw) for q in (2, 3, 4)}
                snr_ref = float(np.linalg.norm(stack(h_ref[band], sw)))
                print(f"  reference optimal SNR (H1 PSD) {snr_ref:.1f}; generation failures {fails}")
                for rule, val in RULES:
                    res = {}
                    for name in ("A", "B"):
                        D = np.column_stack(sets[name])
                        Q, k = subspace(D, rule, val)
                        Qt, _ = np.linalg.qr(np.column_stack([Q, tdir]))
                        res[name] = (k, {q: project_out(T[q], Qt) for q in T})
                    a3, b3 = res["A"][1][3], res["B"][1][3]
                    d = np.linalg.norm(a3 - b3) / np.linalg.norm(b3)
                    s3 = np.linalg.norm(b3) / np.linalg.norm(T[3])
                    s3a = np.linalg.norm(a3) / np.linalg.norm(T[3])
                    cos = float(a3 @ b3 / (np.linalg.norm(a3) * np.linalg.norm(b3)))
                    ds = abs(s3a - s3) / s3
                    s2 = np.linalg.norm(res["B"][1][2]) / np.linalg.norm(T[2])
                    s4 = np.linalg.norm(res["B"][1][4]) / np.linalg.norm(T[4])
                    ok = d < D_CRIT
                    print(f"  {rule}≥{val:<7g} k = {res['A'][0]:3d}/{res['B'][0]:3d}   survival f³ {s3:.4f}  "
                          f"(f² {s2:.4f}, f⁴ {s4:.4f})   d = {d:.3e}  cos {cos:.5f}  Δs/s {ds:.3e}  "
                          f"{'OK' if ok else 'NOT OK'}",
                          flush=True)
                    rows.append({"event": ev, "label": label, "rule": f"{rule}>={val}",
                                 "k_A": res["A"][0], "k_B": res["B"][0], "survival_f3": s3,
                                 "survival_f2": s2, "survival_f4": s4, "stability": d, "cos": cos,
                                 "dsurv": ds, "ok": ok,
                                 "snr_ref": snr_ref, "fails": fails})
                print(f"  ({gen.count} waveforms, {time.time() - t0:.0f} s)", flush=True)
    keys = ["event", "label", "rule", "k_A", "k_B", "survival_f3", "survival_f2", "survival_f4",
            "stability", "cos", "dsurv", "ok", "snr_ref", "fails"]
    with open(os.path.join(args.out, "f3_pca.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print("\nrules passing (d < 0.05) for ALL event/model pairs:")
    for rule, val in RULES:
        sub = [r for r in rows if r["rule"] == f"{rule}>={val}"]
        print(f"  {rule}≥{val:<7g} {sum(r['ok'] for r in sub)}/{len(sub)}")


if __name__ == "__main__":
    main()
