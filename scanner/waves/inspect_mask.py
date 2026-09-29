#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inspect_mask.py — geometry of the valid region (mask_Z) and the largest
hole-free rectangle. Reads only the mask, the coordinates and ONE frame.
Computes no spectra and does not look at the dispersion.

  python inspect_mask.py data\\Surfaces_20111004_113800_short.nc
"""
import sys

import numpy as np
import xarray as xr

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def largest_rectangle(mask):
    """Largest rectangle of ones (row-histogram algorithm)."""
    rows, cols = mask.shape
    heights = np.zeros(cols, dtype=int)
    best = (0, 0, 0, 0, 0)  # area, row0, row1, col0, col1
    for r in range(rows):
        heights = np.where(mask[r], heights + 1, 0)
        stack = []
        for c in range(cols + 1):
            h = heights[c] if c < cols else 0
            start = c
            while stack and stack[-1][1] >= h:
                s, sh = stack.pop()
                area = sh * (c - s)
                if area > best[0]:
                    best = (area, r - sh + 1, r, s, c - 1)
                start = s
            stack.append((start, h))
    return best


path = sys.argv[1]
ds = xr.open_dataset(path)
mask = ds["mask_Z"].values > 0.5
X, Y = ds["X"].values, ds["Y"].values
print(f"grid: {mask.shape}, valid points: {mask.sum()} ({mask.mean():.1%})")
print(f"X: {np.nanmin(X):.2f} … {np.nanmax(X):.2f} m,  Y: {np.nanmin(Y):.2f} … {np.nanmax(Y):.2f} m")
dx = np.nanmedian(np.abs(np.diff(X, axis=1)))
dy = np.nanmedian(np.abs(np.diff(Y, axis=0)))
print(f"step: dx = {dx:.3f} m, dy = {dy:.3f} m")

area, r0, r1, c0, c1 = largest_rectangle(mask)
ny, nx = r1 - r0 + 1, c1 - c0 + 1
print(f"\nlargest hole-free rectangle: rows {r0}…{r1}, columns {c0}…{c1}")
print(f"  size: {nx} × {ny} points = {nx * dx:.2f} × {ny * dy:.2f} m")
print(f"  k resolution: {2 * np.pi / (nx * dx):.3f} × {2 * np.pi / (ny * dy):.3f} rad/m")

t = ds["time"].values
dt = np.median(np.diff(t).astype("timedelta64[ns]").astype(float)) / 1e9
print(f"\ntime: {len(t)} frames, step {dt:.4f} s, duration {len(t) * dt / 60:.1f} min")
gaps = np.sum(np.abs(np.diff(t).astype(float) / 1e9 - dt) > 0.5 * dt)
print(f"irregular time steps: {gaps}")

z0 = ds["Z"].isel(time=0).values
inside = z0[r0:r1 + 1, c0:c1 + 1]
print(f"\nframe 0: NaN inside the rectangle: {np.isnan(inside).sum()}, "
      f"NaN outside the mask: {np.isnan(z0[~mask]).mean():.1%}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.imshow(mask, origin="lower", cmap="gray")
    ax.add_patch(plt.Rectangle((c0 - 0.5, r0 - 0.5), nx, ny, fill=False, ec="r", lw=2))
    ax.set_title("mask_Z (white = valid), red = rectangle")
    fig.savefig("mask.png", dpi=110)
    print("\nplot: mask.png")
except ImportError:
    pass
