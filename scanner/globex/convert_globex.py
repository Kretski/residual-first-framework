#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_globex.py — converts WG_<case>.txt to binary .npy (float32) and
reports missing values per column. No wave analysis.

  python convert_globex.py data\\A1

Output in the same folder:
  WG_<case>.npy   — array (time, gauge), metres, 128 Hz
  x_<case>.npy    — gauge positions [m]
"""
import os
import sys
import time

import numpy as np
import pandas as pd

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

folder = sys.argv[1]
case = os.path.basename(os.path.normpath(folder))
fx = os.path.join(folder, f"xWG_{case}.txt")
fw = os.path.join(folder, f"WG_{case}.txt")

x = np.loadtxt(fx).ravel()
np.save(os.path.join(folder, f"x_{case}.npy"), x)
print(f"positions: {len(x)}")

t0 = time.time()
parts = []
for i, chunk in enumerate(pd.read_csv(fw, sep=r"\s+", header=None, dtype=np.float32,
                                       engine="c", chunksize=50_000)):
    parts.append(chunk.to_numpy(dtype=np.float32))
    print(f"  read {sum(p.shape[0] for p in parts)} lines ({time.time() - t0:.0f} s)", flush=True)
W = np.concatenate(parts, axis=0)
del parts
if W.shape[1] != len(x):
    print(f"WARNING: columns {W.shape[1]} ≠ positions {len(x)}")
out = os.path.join(folder, f"WG_{case}.npy")
np.save(out, W)
print(f"\nwritten {out}: shape {W.shape}, {os.path.getsize(out) / 1e6:.0f} MB")

fs = 128.0
bad = ~np.isfinite(W)
print(f"missing values in total: {bad.sum()} ({bad.mean():.2%})")
cols = np.flatnonzero(bad.any(axis=0))
if len(cols):
    print("columns with missing values (position → last valid row → minutes):")
    for c in cols:
        good = np.flatnonzero(~bad[:, c])
        last = good[-1] if len(good) else -1
        print(f"  column {c:3d}  x = {x[c]:7.3f} m  → row {last:6d}  ({(last + 1) / fs / 60:5.1f} min), "
              f"missing {bad[:, c].sum()}")
else:
    print("no columns with missing values")
