#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inspect_nc.py — shows the description of a NetCDF file (dimensions, variables,
units, attributes) WITHOUT loading the data into memory and without analysis.

  python inspect_nc.py data\\Surfaces_20111004_113800_short.nc
"""
import os
import sys

import xarray as xr

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

path = sys.argv[1]
print("file size:", round(os.path.getsize(path) / 1e9, 2), "GB\n")
ds = xr.open_dataset(path)          # lazy open: the data are not read
print(ds)
print("\n--- file attributes ---")
for k, v in ds.attrs.items():
    print(f"{k}: {v}")
print("\n--- variables ---")
for name, var in ds.variables.items():
    units = var.attrs.get("units", "")
    long = var.attrs.get("long_name", "")
    print(f"{name}: {var.dims} {var.shape} [{units}] {long}")
