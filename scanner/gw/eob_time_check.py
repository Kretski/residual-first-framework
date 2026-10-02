#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eob_time_check.py — diagnose the time/phase convention of the time-domain (EOB)
path used in residual_check.py.

At the SAME parameters (the EOB label's maximum-likelihood sample) three
frequency-domain h_+ are compared:
  (a) IMRPhenomXPHM, generated directly in the frequency domain;
  (b) the EOB model through our path (TD waveform, rfft, epoch phase factor);
  (c) the EOB model through LALSimulation's SimInspiralFD (official TD→FD
      conditioning).
For each pair: noise-weighted match maximised over a time shift (±0.2 s) and
an overall phase, and the best time shift.
Expected if our path is right: match(b, c) ≈ 1 at Δt ≈ 0; match(b, a) high
(different models, same parameters) at small Δt.

  python eob_time_check.py --event GW150914
"""
import argparse
import csv
import os
import sys

import h5py
import numpy as np

import lal
import lalsimulation as ls

import f3_derivatives as f3
import residual_check as rc

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def via_simfd(label, sample, st, f):
    p = {k: sample[k] for k in f3.STEPS}
    m1, m2, s, iota = f3.physical(p, st)
    approx = ls.SEOBNRv4PHM if "v4PHM" in label else ls.SEOBNRv5PHM
    dur = st["duration"]
    hp, hc = ls.SimInspiralFD(m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s,
                              p["luminosity_distance"] * 1e6 * lal.PC_SI, iota, p["phase"], 0, 0, 0,
                              1.0 / dur, st["f_wf"], st["srate"] / 2, st["f_ref"],
                              f3.laldict(st["wf_args"], set()), approx)
    out = np.zeros(len(f), complex)
    n = min(hp.data.length, len(f))
    out[:n] = hp.data.data[:n]
    # SimInspiralFD places the waveform peak at t = epoch offset; bring it to t = 0
    out *= np.exp(-2j * np.pi * f * (-float(hp.epoch)))
    return out


def match(a, b, w, df, band, pad=16):
    x = np.zeros(len(a) * pad, complex)
    x[:len(a)] = np.where(band, np.conj(a) * b * w, 0)
    z = np.fft.ifft(x) * len(x)
    k = int(np.argmax(np.abs(z)))
    na = np.sqrt(np.sum((np.abs(a) ** 2 * w)[band]))
    nb = np.sqrt(np.sum((np.abs(b) ** 2 * w)[band]))
    dt = 1.0 / (len(x) * df)
    shift = k * dt if k < len(x) // 2 else (k - len(x)) * dt
    return float(np.abs(z[k]) / (na * nb)), shift


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--event", default="GW150914")
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    args = ap.parse_args()
    r = next(x for x in csv.DictReader(open(args.f5, encoding="utf-8")) if x["commonName"] == args.event)
    pe_file = r["file"].split(" ")[0]
    with h5py.File(os.path.join(args.pe_dir, pe_file), "r") as fh:
        label = r["eob_label"]
        sample, cfg, psd = f3.read_label(fh, label)
        st = f3.settings(cfg)
        f, hp_b, _ = rc.polarizations(label, sample, st)
        st_x = dict(st, wf_args={})
        _, hp_a, _ = rc.polarizations("IMRPhenomXPHM", sample, st_x)
        try:
            hp_c = via_simfd(label, sample, st, f)
        except Exception as e:
            print(f"SimInspiralFD failed: {type(e).__name__}: {e}")
            hp_c = None
    band = (f >= st["f_an"]) & (f <= st["f_high"])
    w = 1.0 / np.interp(f, psd[:, 0], psd[:, 1], left=np.inf, right=np.inf)
    df = f[1] - f[0]
    print(f"{args.event} {label}: f_wf {st['f_wf']} Hz, f_ref {st['f_ref']} Hz, srate {st['srate']}, duration {st['duration']} s")
    pairs = [("ours (b) vs XPHM (a)", hp_b, hp_a)]
    if hp_c is not None:
        pairs += [("ours (b) vs SimInspiralFD (c)", hp_b, hp_c), ("SimInspiralFD (c) vs XPHM (a)", hp_c, hp_a)]
    for name, x, y in pairs:
        m, sh = match(x, y, w, df, band)
        print(f"  {name:34s} match {m:.4f}   best time shift {sh * 1e3:+8.2f} ms")


if __name__ == "__main__":
    main()
