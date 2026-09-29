#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inspect_mask.py — геометрия на валидната област (mask_Z) и най-големият
правоъгълник без дупки. Чете само маската, координатите и ЕДИН кадър.
Не прави спектри и не гледа дисперсията.

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
    """Най-голям правоъгълник от единици (алгоритъм с хистограма по редове)."""
    rows, cols = mask.shape
    heights = np.zeros(cols, dtype=int)
    best = (0, 0, 0, 0, 0)  # площ, ред0, ред1, кол0, кол1
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
print(f"решетка: {mask.shape}, валидни точки: {mask.sum()} ({mask.mean():.1%})")
print(f"X: {np.nanmin(X):.2f} … {np.nanmax(X):.2f} m,  Y: {np.nanmin(Y):.2f} … {np.nanmax(Y):.2f} m")
dx = np.nanmedian(np.abs(np.diff(X, axis=1)))
dy = np.nanmedian(np.abs(np.diff(Y, axis=0)))
print(f"стъпка: dx = {dx:.3f} m, dy = {dy:.3f} m")

area, r0, r1, c0, c1 = largest_rectangle(mask)
ny, nx = r1 - r0 + 1, c1 - c0 + 1
print(f"\nнай-голям правоъгълник без дупки: редове {r0}…{r1}, колони {c0}…{c1}")
print(f"  размер: {nx} × {ny} точки = {nx * dx:.2f} × {ny * dy:.2f} m")
print(f"  резолюция по k: {2 * np.pi / (nx * dx):.3f} × {2 * np.pi / (ny * dy):.3f} rad/m")

t = ds["time"].values
dt = np.median(np.diff(t).astype("timedelta64[ns]").astype(float)) / 1e9
print(f"\nвреме: {len(t)} кадъра, стъпка {dt:.4f} s, продължителност {len(t) * dt / 60:.1f} min")
gaps = np.sum(np.abs(np.diff(t).astype(float) / 1e9 - dt) > 0.5 * dt)
print(f"неравни стъпки във времето: {gaps}")

z0 = ds["Z"].isel(time=0).values
inside = z0[r0:r1 + 1, c0:c1 + 1]
print(f"\nкадър 0: NaN в правоъгълника: {np.isnan(inside).sum()}, "
      f"NaN извън маската: {np.isnan(z0[~mask]).mean():.1%}")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.imshow(mask, origin="lower", cmap="gray")
    ax.add_patch(plt.Rectangle((c0 - 0.5, r0 - 0.5), nx, ny, fill=False, ec="r", lw=2))
    ax.set_title("mask_Z (бяло = валидно), червено = правоъгълник")
    fig.savefig("mask.png", dpi=110)
    print("\nграфика: mask.png")
except ImportError:
    pass
