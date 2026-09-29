#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inspect_nc.py — показва описанието на NetCDF файл (измерения, променливи,
единици, атрибути), БЕЗ да зарежда данните в паметта и без анализ.

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
print("размер на файла:", round(os.path.getsize(path) / 1e9, 2), "GB\n")
ds = xr.open_dataset(path)          # мързеливо отваряне: данните не се четат
print(ds)
print("\n--- атрибути на файла ---")
for k, v in ds.attrs.items():
    print(f"{k}: {v}")
print("\n--- променливи ---")
for name, var in ds.variables.items():
    units = var.attrs.get("units", "")
    long = var.attrs.get("long_name", "")
    print(f"{name}: {var.dims} {var.shape} [{units}] {long}")
