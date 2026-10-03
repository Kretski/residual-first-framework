#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check_calib_sigma.py — what do the envelope columns 3–6 mean (±1σ or a wider
interval)? Read-only. MODULE1_DESIGN §6, item 2.

priors/calibration/<IFO> has 7 columns: f, A_med, phi_med, A_lo, phi_lo, A_hi, phi_hi.
The analysis prior of each spline node is stored as text (priors/... |S150) and as
5000 samples (priors/samples/recalib_<IFO>_<amplitude|phase>_<k>). For each node:
  - node frequency from recalib_<IFO>_frequency_<k>;
  - std and mean of the prior samples;
  - (hi − lo)/2 and the median of the envelope interpolated at the node frequency.
If std ≈ (hi − lo)/2 the envelope columns are ±1σ; if std ≈ (hi − lo)/3.29 they are a
90% interval. The amplitude prior is for δA = A − 1.

  python check_calib_sigma.py --pe-dir ~/gwdata/pe --f5 f5/f5_results.csv
"""
import argparse
import csv
import os
import re
import sys

import h5py
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

EVENTS = ["GW150914", "GW200129_065458", "GW230627_015337"]


def txt(x):
    x = x[()] if isinstance(x, h5py.Dataset) else x
    if isinstance(x, np.ndarray):
        x = x.ravel()[0] if x.size else ""
    return x.decode(errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    args = ap.parse_args()
    rows = {r["commonName"]: r for r in csv.DictReader(open(args.f5, encoding="utf-8"))}
    for ev in EVENTS:
        r = rows[ev]
        with h5py.File(os.path.join(args.pe_dir, r["file"].split(" ")[0]), "r") as fh:
            for label in (r["xphm_label"], r["eob_label"]):
                if not label or label not in fh or "priors" not in fh[label]:
                    continue
                pr = fh[label]["priors"]
                texts, samples = {}, {}
                def visit(name, obj):
                    if not isinstance(obj, h5py.Dataset):
                        return
                    base = name.split("/")[-1]
                    if not base.startswith("recalib_"):
                        return
                    if obj.dtype.kind in "SOU":
                        texts[base] = txt(obj)
                    elif obj.ndim == 1 and obj.shape[0] > 10:
                        samples[base] = np.asarray(obj[()], float)
                pr.visititems(visit)
                print(f"\n=== {ev}  {label}")
                if not texts and not samples:
                    print("  no recalib priors found")
                    continue
                ex = [k for k in texts if "amplitude_0" in k or "frequency_0" in k][:2]
                for k in ex:
                    print(f"  prior text {k}: {texts[k]}")
                if "calibration" not in pr:
                    print("  no priors/calibration envelope")
                    continue
                for ifo in sorted(pr["calibration"].keys()):
                    env = np.asarray(pr["calibration"][ifo][()], float)
                    f = env[:, 0]
                    print(f"  {ifo}:  node   f[Hz]   | ampl: std(prior)  (hi-lo)/2  ratio   mean(prior)  A_med-1 "
                          f"| phase: std(prior)  (hi-lo)/2  ratio")
                    for k in range(30):
                        fk_key = f"recalib_{ifo}_frequency_{k}"
                        if fk_key not in texts and fk_key not in samples:
                            break
                        if fk_key in samples:
                            fk = float(np.median(samples[fk_key]))
                        else:
                            m = re.search(r"([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)", texts[fk_key].split("(")[-1])
                            fk = float(m.group(1)) if m else float("nan")
                        a = samples.get(f"recalib_{ifo}_amplitude_{k}")
                        p = samples.get(f"recalib_{ifo}_phase_{k}")
                        if a is None or p is None or not np.isfinite(fk):
                            print(f"       {k:2d}  {fk:8.1f}   (samples missing)")
                            continue
                        hw_a = 0.5 * (np.interp(fk, f, env[:, 5]) - np.interp(fk, f, env[:, 3]))
                        hw_p = 0.5 * (np.interp(fk, f, env[:, 6]) - np.interp(fk, f, env[:, 4]))
                        print(f"       {k:2d}  {fk:8.1f}   |  {a.std():.4f}     {hw_a:.4f}   {a.std() / hw_a:5.2f}"
                              f"   {a.mean():+.4f}     {np.interp(fk, f, env[:, 1]) - 1:+.4f} "
                              f"|  {p.std():.4f}     {hw_p:.4f}   {p.std() / hw_p:5.2f}")
    print("\nratio ≈ 1.00 → envelope columns are ±1σ;  ratio ≈ 0.61 → 90% interval.")


if __name__ == "__main__":
    main()
