#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
null_effect_sizes.py — distribution of |a| (effect size at k_ref) in the
NULL synthetic realizations, per test and template. Used to calibrate the
minimum effect size (the systematic floor of the method).

  python null_effect_sizes.py results\\synth_null_A200\\synth_null.json results\\synth_null\\synth_null.json
"""
import json
import sys

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

vals = {}
stat_flags = {}
for path in sys.argv[1:]:
    with open(path, encoding="utf-8") as fh:
        runs = json.load(fh)["runs"]
    for run in runs:
        for test, r in run.items():
            for q, v in r["per"].items():
                key = (test, q)
                vals.setdefault(key, []).append(abs(v["a"]))
                stat_flags.setdefault(key, []).append(v["padj"] < 0.01)

print("test template |   N | |a| median   90 %     95 %     99 %     max   | statistically significant (without effect-size floor)")
for (test, q), a in sorted(vals.items()):
    a = np.array(a)
    qs = np.percentile(a, [50, 90, 95, 99])
    sf = np.mean(stat_flags[(test, q)])
    print(f"  {test}   k^{q}   | {len(a):3d} | {qs[0]:.1e}  {qs[1]:.1e}  {qs[2]:.1e}  {qs[3]:.1e}  {a.max():.1e} | {sf:.3f}")
