#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inspect_globex.py — description of the GLOBEX files, no analysis.
Reads the positions (xWG) in full, and the large file (WG) only at the start,
counting its lines in chunks without loading it into memory.

  python inspect_globex.py data\\A1
"""
import os
import re
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

folder = sys.argv[1]
case = os.path.basename(os.path.normpath(folder))
fx = os.path.join(folder, f"xWG_{case}.txt")
fw = os.path.join(folder, f"WG_{case}.txt")


def split(line):
    return [t for t in re.split(r"[,\s;]+", line.strip()) if t]


print(f"=== {fx} ===")
with open(fx, encoding="utf-8", errors="replace") as fh:
    raw = fh.read()
print("first 300 characters:\n" + raw[:300])
vals = np.array([float(t) for t in split(raw.replace("\n", " "))
                 if re.fullmatch(r"[-+0-9.eEnaN]+", t)])
print(f"\nnumbers: {len(vals)}, from {np.nanmin(vals):.3f} to {np.nanmax(vals):.3f}")
if len(vals) > 1:
    d = np.diff(np.sort(vals))
    print(f"steps between neighbouring positions: min {d[d > 0].min():.3f}, median {np.median(d):.3f}, max {d.max():.3f}")

print(f"\n=== {fw} ({os.path.getsize(fw) / 1e9:.2f} GB) ===")
with open(fw, encoding="utf-8", errors="replace") as fh:
    head = [fh.readline() for _ in range(3)]
for i, line in enumerate(head):
    toks = split(line)
    print(f"line {i + 1}: {len(toks)} values; start: {' '.join(toks[:6])} …")

n_lines = 0
with open(fw, "rb") as fh:
    while True:
        chunk = fh.read(64 * 1024 * 1024)
        if not chunk:
            break
        n_lines += chunk.count(b"\n")
print(f"\nnumber of lines: {n_lines}")
print(f"at 128 Hz this is {n_lines / 128 / 60:.1f} min, if lines are time steps")

with open(fw, encoding="utf-8", errors="replace") as fh:
    sample = " ".join(fh.readline() for _ in range(2000))
toks = split(sample)
nan_like = sum(t.lower() in ("nan", "-999", "-9999", "9999") for t in toks)
print(f"in the first 2000 lines: {len(toks)} values, of which empty/NaN-like: {nan_like}")
