#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_offsource.py — check of the Neyman belts with real off-source detector noise
(MODULE1_DESIGN v0.9.12, §6 item 7 and §7). The ONLY part of Module 1 that reads strain,
and only off-source segments.

Per event and model (geometry from the module1_run cache):
  - off-source segments as in residual_check.py (duration of the PE analysis, spaced by
    64 s on both sides of the event, up to 8 valid ones); a segment is refused if it
    overlaps the window [t_ev - duration - 2 s, t_ev + 2 s + 2 s] of ANY event of the
    event list (selected or excluded), not only the event itself;
  - two injections on EVERY segment, on the same noise:
      C   a posterior sample of set C (the indices following those consumed by sets A and
          B in the fixed permutation of seed 20261001; disjoint from A and B), projected
          with its own sky position, polarization and time;
      ref the reference waveform itself (control);
    both with the calibration mean factor applied, as the reference in the cache;
  - the residual chain of the real analysis: common time shift (±100 ms, 0.05 ms steps)
    and complex amplitude refined against the reference, residual r = d - a h_ref;
  - profile scan of r on the event's own grid (±8 sigma_lin, step 0.05 sigma_lin) and
    q(0) = chi2(0) - min chi2, compared with the event's own Neyman critical value c_0
    at Lambda = 0 (Gaussian belt, same construction as module1_catalog).
Acceptance (catalog level, C injections only, per model): coverage = fraction with
q(0) <= c_0; below the lower end of the 99% binomial interval around 0.90 the analysis
stops before any on-source estimate. Per-event coverage and the ref injections are
diagnostics. Software injections do not pass through the detector calibration: this
check covers non-Gaussian noise and the signal–reference mismatch, not calibration.
"""
import csv
import os

import numpy as np
from scipy.signal.windows import tukey

import module1_catalog as mc
import module1_event as me

OFF_GAP = 64.0
N_OFF = 8
N_C = 8
GUARD_MARGIN = 2.0


def forbidden_windows(event_list_csv, duration):
    rows = list(csv.DictReader(open(event_list_csv, encoding="utf-8")))
    cols = rows[0].keys() if rows else []
    key = next((k for k in ("GPS", "gps", "geocent_time", "GPS_time") if k in cols), None)
    if key is None:
        raise RuntimeError(f"no GPS column in {event_list_csv}; columns: {list(cols)}")
    out = []
    for r in rows:
        try:
            t = float(r[key])
        except (TypeError, ValueError):
            continue
        out.append((t - duration - GUARD_MARGIN, t + me.rc.POST_TRIGGER + GUARD_MARGIN, r.get("commonName", "?")))
    return out


def segment_allowed(s0, duration, windows):
    for a, b, name in windows:
        if s0 < b and s0 + duration > a:
            return False, name
    return True, ""


def refine(d, h, Seff_band, df, tmax=0.1, step=5e-5, fb=None):
    """Common time shift and complex amplitude, as residual_check.check (band vectors)."""
    def ip(a, b, i):
        return 4 * df * np.sum(np.conj(a) * b / Seff_band[i])
    ifos = list(d)
    hh = sum(ip(h[i], h[i], i).real for i in ifos)
    best = (0.0, 0j)
    for t in np.arange(-tmax, tmax + 1e-9, step):
        ph = np.exp(-2j * np.pi * fb * t)
        z = sum(ip(h[i] * ph, d[i], i) for i in ifos)
        if abs(z) > abs(best[1]):
            best = (t, z)
    t, z = best
    a = z / hh
    ph = np.exp(-2j * np.pi * fb * t)
    return {i: d[i] - a * h[i] * ph for i in ifos}, t, a


def event_c0(e, grid, n_belt, seed):
    """Event-only Neyman critical value at Lambda = 0, and the scan objects for that grid."""
    G, Vb, _ = mc.event_trials(e, grid, n_belt, 1, seed)
    Delta = mc.delta_matrix(0.5 * (G + G.T))
    t0 = int(np.argmin(np.abs(grid)))
    c0 = float(np.quantile(mc.qstat(Vb, Delta, t0), 0.90))
    sc = me.ProfileScan(e["A"], e["K"], e["href"], e["net"], e["fb"], grid)
    return c0, sc, t0


def q0_of_residual(r_stacked, e, sc, t0):
    """q(0) and the grid estimate for a residual vector (not of the form s(L) + n)."""
    A = e["A"]
    pr = r_stacked - A.Q @ (A.Q.T @ r_stacked)
    v = sc.b.astype(np.float64) @ pr
    chi = -2 * v + np.diag(sc.G)          # chi2(g) up to a constant
    g = int(np.argmin(chi))
    return float(chi[t0] - chi[g]), g


def binomial_interval(p, n, z=2.576):
    s = np.sqrt(p * (1 - p) / n)
    return p - z * s, p + z * s
