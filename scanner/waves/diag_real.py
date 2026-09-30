#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_real.py — diagnostics of the real-data spectrum (development record only).

Reads the cached 3-D spectrum written by waves_scanner.py (real mode) and shows:
  1) the direction-averaged spectrum E(|k|, f) with the deep-water dispersion
     curve and the cells accepted by the ridge extraction;
  2) the ratio ω_meas / √(g k) of the accepted cells versus |k|;
  3) the directional energy distribution in a frequency band.

This is a DIAGNOSTIC of the measurement step, run after the v0.5.1 result was
recorded. It does not change the registered analysis.

  python diag_real.py results\\v051_real search
"""
import os
import sys

import numpy as np

import waves_scanner as ws

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

out = sys.argv[1]
half = sys.argv[2] if len(sys.argv) > 2 else "search"
c = np.load(os.path.join(out, f"spectrum_{half}.npz"))
P, f, ka, kb = c["P"].astype(np.float64), c["f"], c["ka"], c["kb"]
print(f"spectrum {P.shape}  (f, k_a, k_b)")

KA, KB = np.meshgrid(ka, kb, indexing="ij")
kx, ky = -KA, -KB                                  # propagation vector
kmag = np.hypot(kx, ky)
theta = np.degrees(np.arctan2(ky, kx))

# accepted ridge cells (same rule as the scanner)
d = ws.extract_ridge(P, f, ka, kb)
ratio = d["omega"] / np.sqrt(ws.G * d["k"])
print(f"accepted cells: {len(d['k'])} of {d['n_candidates']}")
for lo, hi in [(5.8, 10), (10, 15), (15, 20), (20, 28), (28, 36.5)]:
    m = (d["k"] >= lo) & (d["k"] < hi)
    if m.sum():
        print(f"  |k| {lo:4.1f}-{hi:4.1f}: n = {m.sum():4d}, median ω/√(gk) = {np.median(ratio[m]):.3f}, "
              f"IQR {np.percentile(ratio[m], 25):.3f}–{np.percentile(ratio[m], 75):.3f}")

# direction-averaged spectrum on a |k| grid
kbins = np.linspace(0, 50, 101)
kc = 0.5 * (kbins[1:] + kbins[:-1])
E = np.zeros((len(f), len(kc)))
flat = kmag.ravel()
idx = np.digitize(flat, kbins) - 1
Pf = P.reshape(P.shape[0], -1)
for j in range(len(kc)):
    sel = idx == j
    if sel.any():
        E[:, j] = Pf[:, sel].mean(axis=1)

# directional energy at 1.2–2 Hz, over the cells in the working band
band = (f >= 1.2) & (f <= 2.0)
kw = (kmag >= ws.k_deep(1.2)) & (kmag <= ws.k_deep(2.0))
Edir = P[band].sum(axis=0)[kw]
th = theta[kw]
hist, edges = np.histogram(th, bins=36, range=(-180, 180), weights=Edir)
peak_dir = 0.5 * (edges[np.argmax(hist)] + edges[np.argmax(hist) + 1])
print(f"dominant propagation direction (1.2–2 Hz): {peak_dir:.0f}°  "
      f"(share in the peak 10° bin: {hist.max() / hist.sum():.1%})")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("(matplotlib not available)")
    sys.exit(0)

fig, ax = plt.subplots(1, 3, figsize=(16, 5))
im = ax[0].pcolormesh(kc, f, np.log10(E + 1e-30), shading="auto", cmap="viridis")
kk = np.linspace(0.1, 50, 400)
ax[0].plot(kk, np.sqrt(ws.G * kk) / (2 * np.pi), "w-", lw=1, label="ω² = g k")
ax[0].plot(d["k"], d["omega"] / (2 * np.pi), "r.", ms=2, label="accepted cells")
ax[0].axvspan(ws.k_deep(1.2), ws.k_deep(3.0), color="w", alpha=0.08)
ax[0].set_xlabel("|k| [rad/m]")
ax[0].set_ylabel("f [Hz]")
ax[0].set_ylim(0, 6)
ax[0].set_title(f"log10 E(|k|, f), {half}")
ax[0].legend(loc="upper left", fontsize=8)
fig.colorbar(im, ax=ax[0])

ax[1].plot(d["k"], ratio, ".", ms=3)
ax[1].axhline(1, color="k", lw=0.8)
ax[1].set_xlabel("|k| [rad/m]")
ax[1].set_ylabel("ω_meas / √(g k)")
ax[1].set_title("accepted cells")

ax[2].bar(0.5 * (edges[1:] + edges[:-1]), hist / hist.sum(), width=10)
ax[2].set_xlabel("propagation direction [deg]")
ax[2].set_ylabel("energy share, 1.2–2 Hz")
ax[2].set_title("directional distribution")
fig.tight_layout()
png = os.path.join(out, f"diag_{half}.png")
fig.savefig(png, dpi=110)
print(f"plot: {png}")
