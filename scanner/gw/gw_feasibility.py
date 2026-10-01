#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gw_feasibility.py — Module 1 feasibility checks F1–F2 (no GW data involved).

F1: prints the installed versions of lalsuite, pyseobnr and gwpy.
F2: generates one waveform per model for a typical binary black hole and
    measures the generation time (median of several repetitions):
      - IMRPhenomXPHM   (lalsimulation, frequency domain)
      - SEOBNRv4PHM     (lalsimulation, time domain)    — discovery cohort
      - SEOBNRv5PHM     (pyseobnr, time domain)         — confirmation cohort
F6 (rough): extrapolates the cost of the numerical derivatives needed for the
    orthogonalization (§6 of MODULE1_DESIGN): 2 waveforms per parameter.

  python gw_feasibility.py
"""
import statistics
import sys
import time

import numpy as np

REPEATS = 3
N_PARAMS = 15                       # parameters differentiated in §6
M1, M2 = 36.0, 29.0                 # GW150914-like, solar masses
CHI1 = (0.1, 0.0, 0.3)
CHI2 = (0.0, 0.1, -0.2)
DIST_MPC = 400.0
INCL = 0.4
F_LOW = 20.0
F_REF = 20.0
SRATE = 4096.0
DURATION = 8.0


def versions():
    print("=== F1: versions ===")
    for name in ("lal", "lalsimulation", "pyseobnr", "gwpy", "numpy", "scipy"):
        try:
            mod = __import__(name)
            v = getattr(mod, "__version__", None)
            if v is None and name == "lal":
                v = mod.VCSInfo.version
            print(f"  {name:14s} {v}")
        except Exception as e:
            print(f"  {name:14s} NOT AVAILABLE ({type(e).__name__}: {e})")
    print(f"  python         {sys.version.split()[0]}")


def timed(fn, label):
    times = []
    out = None
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        out = fn()
        times.append(time.perf_counter() - t0)
    med = statistics.median(times)
    print(f"  {label:16s} median {med:8.3f} s  (runs: {', '.join(f'{t:.3f}' for t in times)})")
    return med, out


def phenom_xphm():
    import lal
    import lalsimulation as ls
    deltaF = 1.0 / DURATION
    return ls.SimInspiralChooseFDWaveform(
        M1 * lal.MSUN_SI, M2 * lal.MSUN_SI, *CHI1, *CHI2,
        DIST_MPC * 1e6 * lal.PC_SI, INCL, 0.0, 0.0, 0.0, 0.0,
        deltaF, F_LOW, SRATE / 2, F_REF, None, ls.IMRPhenomXPHM)


def seob_v4phm():
    import lal
    import lalsimulation as ls
    return ls.SimInspiralChooseTDWaveform(
        M1 * lal.MSUN_SI, M2 * lal.MSUN_SI, *CHI1, *CHI2,
        DIST_MPC * 1e6 * lal.PC_SI, INCL, 0.0, 0.0, 0.0, 0.0,
        1.0 / SRATE, F_LOW, F_REF, None, ls.SEOBNRv4PHM)


def seob_v5phm():
    from pyseobnr.generate_waveform import GenerateWaveform
    params = {
        "mass1": M1, "mass2": M2,
        "spin1x": CHI1[0], "spin1y": CHI1[1], "spin1z": CHI1[2],
        "spin2x": CHI2[0], "spin2y": CHI2[1], "spin2z": CHI2[2],
        "deltaT": 1.0 / SRATE, "f22_start": F_LOW, "f_ref": F_REF,
        "phi_ref": 0.0, "distance": DIST_MPC, "inclination": INCL,
        "approximant": "SEOBNRv5PHM",
    }
    return GenerateWaveform(params).generate_td_polarizations()


def main():
    versions()
    print("\n=== F2: one waveform per model (GW150914-like BBH, f_low = 20 Hz) ===")
    cost = {}
    for label, fn in (("IMRPhenomXPHM", phenom_xphm), ("SEOBNRv4PHM", seob_v4phm),
                      ("SEOBNRv5PHM", seob_v5phm)):
        try:
            med, out = timed(fn, label)
            hp = out[0]
            n = hp.data.length if hasattr(hp, "data") else len(hp)
            print(f"  {'':16s} output length {n}, finite: {bool(np.all(np.isfinite(hp.data.data))) if hasattr(hp, 'data') else 'n/a'}")
            cost[label] = med
        except Exception as e:
            print(f"  {label:16s} FAILED: {type(e).__name__}: {e}")

    print("\n=== F6 (rough): cost of the tangent basis per event and model ===")
    print(f"  {N_PARAMS} parameters × 2 waveforms (central differences) + 1 reference")
    for label, t in cost.items():
        per_event = (2 * N_PARAMS + 1) * t
        print(f"  {label:16s} ≈ {per_event:7.1f} s per event; "
              f"× 1000 null/ladder realizations ≈ {per_event * 1000 / 3600:7.1f} h "
              f"(if the basis is recomputed each time)")
    print("\nNote: the tangent basis can be computed once per event at the reference point")
    print("and reused for all injections around it; then the dominant cost is the")
    print("injected waveforms themselves (1 per realization).")


if __name__ == "__main__":
    main()
