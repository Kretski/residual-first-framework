#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
blocks_globex.py — checks the session boundaries of the inshore trolley
(blocks of 11 gauges at 0.37 m). Simultaneously recorded gauges share
electrical interference (50 Hz mains, data-acquisition noise); gauges from
different sessions do not. The coherence MAGNITUDE is computed in the 40–60 Hz
band, where there are no waves. The phase is not computed.

  python blocks_globex.py data\\A1
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
seg = slice(0, int(20 * 60 * fs))
nper = int(8 * fs)

# presumed blocks: 9 × 11 columns starting at 82
starts = [82 + 11 * j for j in range(9)]
block = -np.ones(len(x), int)
for j, s in enumerate(starts):
    block[s:s + 11] = j


def coh(a, b, lo, hi):
    f, C = signal.coherence(np.asarray(W[seg, a], float), np.asarray(W[seg, b], float),
                            fs=fs, nperseg=nper)
    m = (f >= lo) & (f <= hi)
    return np.sqrt(C[m]).mean(), np.sqrt(C[(f >= 49.5) & (f <= 50.5)]).max()


print("neighbouring pairs in the inshore-trolley zone")
print("  i   x_i [m]  block  boundary?  |coh| 40–60 Hz   max |coh| near 50 Hz")
within, across = [], []
for c in range(80, len(x) - 1):
    m, p50 = coh(c, c + 1, 40, 60)
    edge = block[c] != block[c + 1]
    (across if edge else within).append(m)
    print(f"{c:4d}  {x[c]:7.3f}   {block[c]:3d}   {'YES' if edge else '   '}        "
          f"{m:.3f}            {p50:.3f}", flush=True)
print(f"\nmean within blocks: {np.mean(within):.3f}   at boundaries: {np.mean(across):.3f}")

print("\ncontrol: pairs 11 columns apart (always in different blocks)")
for c in range(82, 82 + 11 * 8, 11):
    m, _ = coh(c, c + 11, 40, 60)
    print(f"  {c} ↔ {c + 11}: |coh| 40–60 Hz = {m:.3f}")
