#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_diag_calibration.py — Addendum C diagnostic (pre-specified, separate from the
primary analysis).

Events: exactly those excluded from the search cohort solely by the section 3
calibration rule (commit 6f46913): GW170608, GW190707_093326, GW190728_064510,
GW190924_021846. Their labels have no recalibration priors.

What this wrapper changes, and nothing else:
  * module1_calibration.load is replaced by a version that first tries the frozen
    loader and, only if it fails, builds the same per-detector dictionary
    {freqs, mu_a, sg_a, mu_p, sg_p} from the calibration_envelope table of the label
    (or of the IMRPhenomXPHM label of the same event), the way bilby builds spline
    priors from envelope files:
        nodes  n log-spaced in [f_low, f_high] of the analysis (n from the PE config,
               default 10)
        mu_a = median amplitude − 1,  sg_a = (upper − lower) / 2
        mu_p = median phase,          sg_p = (upper − lower) / 2
    Envelope columns: frequency, amplitude median, phase median, amplitude −1σ,
    phase −1σ, amplitude +1σ, phase +1σ. If a detector has no envelope either, the
    calibration term is omitted for that detector (zero widths) and the result is
    marked as not including calibration uncertainty there.
  * the status of every output is 'diagnostic-C' (never 'official'), and the cache
    directory defaults to ~/gw/cache_diag_C, so nothing mixes with the primary cache.

The frozen files (module1_run.py, module1_event.py, ...) are imported unchanged.

Usage (same stages as module1_run.py):
  python module1_diag_calibration.py build    --labels both
  python module1_diag_calibration.py checks   --labels both --out diag_C_checks.csv
  python module1_diag_calibration.py catalog  --labels xphm --out diag_C_catalog_xphm.json
  python module1_diag_calibration.py catalog  --labels eob  --out diag_C_catalog_eob.json
  python module1_diag_calibration.py offsource --labels both --out diag_C_offsource.csv
"""

import os
import sys

import numpy as np

import module1_calibration as cal
import module1_run as run

DIAG_EVENTS = "GW170608,GW190707_093326,GW190728_064510,GW190924_021846"
STATUS = "diagnostic-C"

_frozen_load = cal.load


def envelope_priors(table, f_lo, f_hi, n_nodes):
    """Per-detector prior dictionary from a calibration-envelope table (n, 7)."""
    t = np.asarray(table, float)
    if t.ndim != 2 or t.shape[1] < 7:
        raise ValueError(f"unexpected envelope shape {t.shape}")
    fr = t[:, 0]
    nodes = np.geomspace(f_lo, f_hi, n_nodes)
    a_med, p_med = np.interp(nodes, fr, t[:, 1]), np.interp(nodes, fr, t[:, 2])
    a_lo, p_lo = np.interp(nodes, fr, t[:, 3]), np.interp(nodes, fr, t[:, 4])
    a_hi, p_hi = np.interp(nodes, fr, t[:, 5]), np.interp(nodes, fr, t[:, 6])
    # amplitude stored as the factor (1 + dA) in LVK envelope files; guard against
    # tables that already store dA
    offset = 1.0 if np.median(t[:, 1]) > 0.5 else 0.0
    return {"freqs": nodes, "mu_a": a_med - offset, "sg_a": (a_hi - a_lo) / 2.0,
            "mu_p": p_med, "sg_p": (p_hi - p_lo) / 2.0}


def _config_value(fh, label, keys, default=None):
    g = fh[label]
    if "config_file" not in g:
        return default
    found = {}

    def visit(name, obj):
        if hasattr(obj, "shape"):
            found[name.split("/")[-1].replace("_", "-")] = obj[()]
    g["config_file"].visititems(visit)
    for k in keys:
        if k in found:
            v = found[k]
            v = v.decode() if isinstance(v, bytes) else (v[0].decode() if hasattr(v, "__len__") and len(v) and isinstance(v[0], bytes) else v)
            try:
                return float(v)
            except (TypeError, ValueError):
                continue
    return default


def diag_load(fh, label, xphm_label, ifos):
    try:
        return _frozen_load(fh, label, xphm_label, ifos)
    except Exception as frozen_error:
        import f3_derivatives as f3
        _, cfg, _ = f3.read_label(fh, label)
        st = f3.settings(cfg)
        n_nodes = int(_config_value(fh, label, ["spline-calibration-nodes", "spcal-nodes"], 10))
        pri, sources = {}, []
        for i in ifos:
            table = None
            for lab in (label, xphm_label):
                g = fh[lab]
                if "calibration_envelope" in g and i in g["calibration_envelope"]:
                    cand = np.asarray(g["calibration_envelope"][i][()], float)
                    if cand.ndim == 2 and cand.shape[0] >= 2 and cand.shape[1] >= 7:   # empty tables occur
                        table = cand
                        sources.append(f"{i}: envelope of {lab}")
                        break
            if table is None:
                nodes = np.geomspace(st["f_an"], st["f_high"], n_nodes)
                z = np.zeros(n_nodes)
                pri[i] = {"freqs": nodes, "mu_a": z, "sg_a": z, "mu_p": z, "sg_p": z}
                sources.append(f"{i}: NO CALIBRATION (no priors, no envelope)")
                continue
            pri[i] = envelope_priors(table, st["f_an"], st["f_high"], n_nodes)
        src = f"Addendum C ({'; '.join(sources)}); frozen loader said: {type(frozen_error).__name__}"
        print(f"  calibration: {src}", flush=True)
        return pri, src


def main():
    cal.load = diag_load
    run.cal.load = diag_load
    run.me.cal.load = diag_load
    run.official_status = lambda *a, **k: (STATUS, "Addendum C diagnostic wrapper")
    argv = sys.argv[1:]
    if "--events" not in argv:
        argv += ["--events", DIAG_EVENTS]
    if "--cache" not in argv:
        argv += ["--cache", os.path.expanduser("~/gw/cache_diag_C")]
    sys.argv = [sys.argv[0]] + argv
    run.main()


if __name__ == "__main__":
    main()
