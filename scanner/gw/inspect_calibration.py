#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inspect_calibration.py — read-only check of the calibration information in the PE
files, before the Module 1 wrapper uses it (MODULE1_DESIGN §6, item 2).

For one event per catalog (GWTC-2.1, GWTC-3, GWTC-4.0) and each label of the
primary test, prints:
  - every HDF5 object whose path contains 'calib' (shape, dtype);
  - for the envelope arrays: column count, first rows, last row, and the
    frequency range;
  - the posterior columns of the calibration parameters (recalib_*), which the
    alternative "maximum-likelihood calibration" rule would use.
No residual, no estimate. Output is only a description of the files.

  python inspect_calibration.py --pe-dir ~/gwdata/pe --f5 f5/f5_results.csv
"""
import argparse
import csv
import os
import sys

import h5py
import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

EVENTS = ["GW150914", "GW200129_065458", "GW230627_015337"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    args = ap.parse_args()
    rows = {r["commonName"]: r for r in csv.DictReader(open(args.f5, encoding="utf-8"))}
    np.set_printoptions(precision=5, suppress=False, linewidth=140)
    for ev in EVENTS:
        r = rows[ev]
        fn = r["file"].split(" ")[0]
        print(f"\n######## {ev}   file {fn}")
        with h5py.File(os.path.join(args.pe_dir, fn), "r") as fh:
            for label in (r["xphm_label"], r["eob_label"]):
                if not label or label not in fh:
                    continue
                g = fh[label]
                print(f"\n=== {label}   top-level keys: {sorted(g.keys())}")
                hits = []
                g.visititems(lambda name, obj: hits.append((name, obj)) if "calib" in name.lower() else None)
                for name, obj in hits:
                    if isinstance(obj, h5py.Dataset):
                        print(f"  dataset {name}: shape {obj.shape}, dtype {obj.dtype}")
                        a = obj[()]
                        if a.dtype.names:
                            print(f"    fields: {a.dtype.names}")
                            a = np.column_stack([np.asarray(a[n], float) for n in a.dtype.names])
                        a = np.asarray(a)
                        if a.ndim == 2 and a.shape[0] > 3 and np.issubdtype(a.dtype, np.number):
                            print(f"    columns: {a.shape[1]};  col0 range {a[:, 0].min():.3g} – {a[:, 0].max():.3g}")
                            print(f"    first rows:\n{a[:3]}")
                            print(f"    last row:\n{a[-1]}")
                    else:
                        print(f"  group   {name}: keys {sorted(obj.keys())}")
                ps = g["posterior_samples"]
                names = (ps.dtype.names if isinstance(ps, h5py.Dataset) and ps.dtype.names
                         else [v.decode() if isinstance(v, bytes) else str(v)
                               for v in ps["parameter_names"][()]])
                rc = [n for n in names if n.lower().startswith("recalib")]
                print(f"  posterior calibration columns: {len(rc)}"
                      + (f"  e.g. {rc[:4]} … {rc[-2:]}" if rc else ""))


if __name__ == "__main__":
    main()
