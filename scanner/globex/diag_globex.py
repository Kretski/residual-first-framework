#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_globex.py — EXPLORATORY diagnostic of the v0.1.0 real result (not part of
the registered analysis). It asks whether the ~1.47× excess of the negative
deviation over linear theory is a depth error or a frequency-dependent effect.

For every measured point (ω, k) the effective depth at which linear theory
ω² = g k tanh(k h) reproduces the measurement is
    h_eff = atanh(ω² / (g k)) / k .
  • depth error only        → h_eff constant in frequency within each block;
  • nonlinearity / other    → h_eff depends on frequency.

  python diag_globex.py data\\A3 results\\v010_real_A3
"""
import os
import sys

import numpy as np
from scipy import stats

import globex_scanner as gs

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

folder, out = sys.argv[1], sys.argv[2]
case = os.path.basename(os.path.normpath(folder))
W = np.load(os.path.join(folder, f"WG_{case}.npy"), mmap_mode="r")
x = np.load(os.path.join(folder, f"x_{case}.npy"))
n = W.shape[0]
halves = {"search": slice(0, n // 2), "confirm": slice(n // 2, 2 * (n // 2))}
os.makedirs(out, exist_ok=True)

results = {}
for name, sl in halves.items():
    d = gs.collect(np.asarray(W[sl]), x)
    arg = d["omega"] ** 2 / (gs.G * d["k"])
    ok = arg < 1
    h_eff = np.full_like(d["k"], np.nan)
    h_eff[ok] = np.arctanh(arg[ok]) / d["k"][ok]
    f = d["omega"] / (2 * np.pi)
    print(f"\n=== {name} ({case}) ===")
    print(" block  h_nominal  median h_eff  ratio   slope dh_eff/df [m/Hz]   p(slope=0)   invalid")
    rows = []
    for b in np.unique(d["block"]):
        m = (d["block"] == b) & np.isfinite(h_eff)
        hn = d["h"][d["block"] == b][0]
        sl_ = stats.linregress(f[m], h_eff[m])
        rows.append((b, hn, np.median(h_eff[m]), sl_.slope, sl_.pvalue))
        print(f"   {b}     {hn:.3f}      {np.median(h_eff[m]):.3f}      {np.median(h_eff[m]) / hn:.3f}"
              f"      {sl_.slope:+.4f}                {sl_.pvalue:.1e}      {int((~ok & (d['block'] == b)).sum())}")
    results[name] = (d, h_eff, f)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit(0)

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
colors = plt.cm.viridis(np.linspace(0, 0.9, 4))
for name, mk in (("search", "o"), ("confirm", "x")):
    d, h_eff, f = results[name]
    for b in np.unique(d["block"]):
        m = d["block"] == b
        ax[0].plot(d["kh"][m], d["r"][m], mk, color=colors[b], ms=3,
                   label=f"block {b} ({name})" if name == "search" else None)
        ax[1].plot(f[m], h_eff[m], mk, color=colors[b], ms=3)
        ax[1].axhline(d["h"][m][0], color=colors[b], lw=0.8, ls="--")
kh = np.linspace(0.3, 0.72, 100)
ax[0].plot(kh, gs.r_linear(kh), "k-", lw=1.5, label="linear theory")
ax[0].set_xlabel("kh (nominal depth)")
ax[0].set_ylabel("r = ω / (k √(g h)) − 1")
ax[0].set_title(f"{case}: residual vs kh (o search, x confirm)")
ax[0].legend(fontsize=7)
ax[1].set_xlabel("f [Hz]")
ax[1].set_ylabel("h_eff [m]  (dashed: nominal)")
ax[1].set_title("effective depth from linear theory")
fig.tight_layout()
png = os.path.join(out, "diag_depth.png")
fig.savefig(png, dpi=110)
print(f"\nplot: {png}")
