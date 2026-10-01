#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inspect_pe.py — Module 1 feasibility check F5 (structure only, no analysis).

Opens one LVK PESummary parameter-estimation release file and reports, per
analysis label: the approximant, the number of posterior samples, whether a
log-likelihood column exists (needed for the maximum-likelihood reference
point), the detectors with PSDs and calibration envelopes, and the waveform
settings stored in the configuration (reference/minimum frequency, waveform
arguments, mode array, duration, sampling frequency).

  python inspect_pe.py <file.h5> [--full-config]
"""
import sys

import h5py
import numpy as np

KEYS_OF_INTEREST = ("reference-frequency", "reference_frequency", "minimum-frequency",
                    "maximum-frequency", "sampling-frequency", "duration",
                    "waveform-approximant", "waveform-arguments-dict", "mode-array",
                    "frequency-domain-source-model", "calibration-model",
                    "spline-calibration-nodes", "distance-marginalization",
                    "phase-marginalization", "time-marginalization")


def text(x):
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode(errors="replace")
    if isinstance(x, np.ndarray):
        if x.dtype.kind in "SO" and x.size:
            return ", ".join(text(v) for v in x.ravel()[:10])
        return str(x.tolist() if x.size < 10 else f"array {x.shape}")
    return str(x)


def describe_label(f, label):
    g = f[label]
    print(f"\n=== {label} ===")
    print("  groups:", ", ".join(sorted(g.keys())))
    if "approximant" in g:
        print("  approximant:", text(g["approximant"][()]))
    ps = g["posterior_samples"] if "posterior_samples" in g else None
    if ps is not None:
        names, n = (), 0
        try:
            if isinstance(ps, h5py.Dataset) and ps.dtype.names:
                names = ps.dtype.names
                n = ps.shape[0]
            elif isinstance(ps, h5py.Group):
                keys = list(ps.keys())
                if "parameter_names" in keys and "samples" in keys:      # older PESummary layout
                    names = tuple(text(v) for v in ps["parameter_names"][()])
                    n = ps["samples"].shape[0]
                else:
                    names = tuple(keys)
                    n = len(ps[keys[0]]) if keys else 0
        except Exception as e:
            print(f"  posterior samples: unreadable layout ({type(e).__name__}: {e})")
        print(f"  posterior samples: {n}, parameters: {len(names)}")
        print("  has log_likelihood:", "log_likelihood" in names)
        wanted = ("mass_1", "mass_2", "chirp_mass", "mass_ratio", "luminosity_distance",
                  "redshift", "geocent_time", "a_1", "a_2", "tilt_1", "tilt_2",
                  "phi_12", "phi_jl", "theta_jn", "psi", "ra", "dec", "phase",
                  "network_matched_filter_snr")
        print("  key parameters present:", ", ".join(w for w in wanted if w in names))
        missing = [w for w in wanted if w not in names]
        if missing:
            print("  key parameters missing:", ", ".join(missing))
    for kind in ("psds", "calibration_envelope"):
        if kind in g:
            sub = g[kind]
            print(f"  {kind}: " + ", ".join(f"{k} {sub[k].shape}" for k in sub.keys()))
        else:
            print(f"  {kind}: NOT PRESENT")
    if "config_file" in g:
        cfg = g["config_file"]
        flat = {}

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                flat[name.split("/")[-1]] = obj[()]
        cfg.visititems(visit)
        print("  config entries:", len(flat))
        shown = [k for k in KEYS_OF_INTEREST if k in flat]
        for k in shown:
            print(f"    {k} = {text(flat[k])[:200]}")
        if not shown or FULL:
            print("  (non-standard or full config listing)")
            for k in sorted(flat):
                print(f"    {k} = {text(flat[k])[:160]}")
    else:
        print("  config_file: NOT PRESENT")


FULL = False


def main():
    global FULL
    path = sys.argv[1]
    FULL = "--full-config" in sys.argv
    with h5py.File(path, "r") as f:
        top = sorted(f.keys())
        print("top-level keys:", ", ".join(top))
        labels = [k for k in top if k not in ("version", "history")]
        for label in labels:
            try:
                describe_label(f, label)
            except Exception as e:
                print(f"\n=== {label} === could not be described: {type(e).__name__}: {e}")


if __name__ == "__main__":
    main()
