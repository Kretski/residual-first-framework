#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sessions_globex.py — column diagnostics (no dispersion analysis):
  1) padded segments: zeros or constant values at the end of the record;
  2) spacing between neighbouring positions (identifies the trolleys);
  3) MAGNITUDE of the coherence between neighbouring columns in the wave band
     and in a high-frequency band. Simultaneously recorded pairs are almost
     fully coherent; pairs from different sessions less so.
     The phase is NOT computed and not shown.

  python sessions_globex.py data\\A1
"""
import os
import sys

import numpy as np
from scipy import signal

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

folder = sys.argv[1]
case = os.path.basename(os.path.normpath(folder))
W = np.load(os.path.join(folder, f"WG_{case}.npy"), mmap_mode="r")
x = np.load(os.path.join(folder, f"x_{case}.npy"))
fs = 128.0
nt, nc = W.shape
print(f"shape {W.shape}, {nt / fs / 60:.1f} min")

# 1) padding: long tails of zeros or of unchanging values
print("\n1) columns with long constant segments (≥ 10 s):")
found = False
for c in range(nc):
    col = np.asarray(W[:, c])
    same = np.concatenate([[0], (np.diff(col) == 0).astype(np.int8), [0]])
    edges = np.flatnonzero(np.diff(same))
    starts, ends = edges[::2], edges[1::2]              # half-open segments
    if len(starts):
        k = int(np.argmax(ends - starts))
        best, best_end = int(ends[k] - starts[k]), int(ends[k] - 1)
    else:
        best, best_end = 0, 0
    if best >= 10 * fs:
        found = True
        start = best_end - best
        print(f"  column {c:3d} x = {x[c]:7.3f} m: constant {best / fs:7.1f} s "
              f"from {start / fs / 60:5.1f} to {best_end / fs / 60:5.1f} min (value {col[best_end]:.4g})")
if not found:
    print("  none")

# 2) + 3) neighbouring pairs
seg = slice(0, int(30 * 60 * fs))           # the first 30 min suffice for the diagnostic
nper = int(32 * fs)
print("\n2–3) neighbouring columns: spacing and coherence magnitude")
print("   i    x_i [m]   Δx [m]   |coh| 0,4–1,5 Hz   |coh| 5–30 Hz")
for c in range(nc - 1):
    a = np.asarray(W[seg, c], dtype=np.float64)
    b = np.asarray(W[seg, c + 1], dtype=np.float64)
    f, C = signal.coherence(a, b, fs=fs, nperseg=nper)
    band = (f >= 0.4) & (f <= 1.5)
    hi = (f >= 5) & (f <= 30)
    print(f"  {c:3d}  {x[c]:7.3f}  {x[c + 1] - x[c]:7.3f}      {np.sqrt(C[band]).mean():.3f}"
          f"             {np.sqrt(C[hi]).mean():.3f}", flush=True)
