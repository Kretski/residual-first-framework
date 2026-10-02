#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f3_derivatives.py — Module 1 feasibility check F3 (no strain data, no residual).

Question: are numerical derivatives of the waveform with respect to the PE
parameters stable at the real reference points? They form the tangent basis
used to orthogonalize the f^p templates (MODULE1_DESIGN §6).

For each test event and model label:
  1. reference point = maximum-likelihood posterior sample of that label;
  2. waveform settings from that label's configuration (f_ref, waveform start
     frequency, band, sampling rate, duration, waveform arguments);
  3. h_+(f) generated with the same model (IMRPhenomXPHM: FD; SEOBNRv4PHM:
     LALSimulation TD; SEOBNRv5PHM: pyseobnr TD; TD waveforms are FFT'd on a
     fixed grid with the epoch phase-corrected);
  4. central differences D(ε), D(ε/2), D(ε/4) for each parameter; convergence
     ratios r1 = |D(ε) − D(ε/2)| / |D(ε/2)| and r2 = |D(ε/2) − D(ε/4)| / |D(ε/4)|
     in the noise-weighted norm with the released H1 PSD.
     Stable: r2 < 1e-2. Second-order truncation gives r2 ≈ r1 / 4.

This checks the numerics of the waveform derivatives only (the detector
projection is linear and does not change their stability).

  python f3_derivatives.py --pe-dir ~/gwdata/pe --out ~/gw/f3
"""
import argparse
import ast
import csv
import os
import re
import sys
import time

import h5py
import numpy as np

import lal
import lalsimulation as ls

EVENTS = {   # event → (file substring, [labels])
    "GW150914": ("GW150914_095045", ["C01:IMRPhenomXPHM", "C01:SEOBNRv4PHM"]),
    "GW190412": ("GW190412_053044", ["C01:IMRPhenomXPHM", "C01:SEOBNRv4PHM"]),
    "GW230814_230901": ("GW230814_230901", ["C00:IMRPhenomXPHM-SpinTaylor", "C00:SEOBNRv5PHM"]),
    "GW231123_135430": ("GW231123_135430", ["C00:IMRPhenomXPHM-SpinTaylor", "C00:SEOBNRv5PHM"]),
}
# parameter → (base step, kind: 'rel' relative / 'abs' absolute, lower, upper)
STEPS = {   # v2: larger angular and spin steps (v1 used 2e-3: noise-dominated, r1 ≈ r2)
    "chirp_mass": (1e-4, "rel", 0.0, np.inf),
    "mass_ratio": (4e-3, "abs", 0.02, 1.0),
    "a_1": (1e-2, "abs", 0.0, 0.99),
    "a_2": (1e-2, "abs", 0.0, 0.99),
    "tilt_1": (2e-2, "abs", 0.0, np.pi),
    "tilt_2": (2e-2, "abs", 0.0, np.pi),
    "phi_12": (2e-2, "abs", -np.inf, np.inf),
    "phi_jl": (2e-2, "abs", -np.inf, np.inf),
    "theta_jn": (2e-2, "abs", 0.0, np.pi),
    "phase": (2e-2, "abs", -np.inf, np.inf),
    "luminosity_distance": (1e-4, "rel", 0.0, np.inf),
}
# EOB models integrate ODEs with adaptive steps and are not smooth on very small
# parameter changes: their base steps are multiplied by this factor (capped below).
EOB_STEP_FACTOR = 10.0
STABLE = 1e-2
NEGLIGIBLE = 1e-2      # |D| / max|D| of the event below this: direction negligible

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def text(x):
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode(errors="replace")
    if isinstance(x, np.ndarray):
        return ", ".join(text(v) for v in x.ravel()[:20])
    return str(x)


def read_label(f, label):
    g = f[label]
    ps = g["posterior_samples"]
    if isinstance(ps, h5py.Dataset) and ps.dtype.names:
        arr = ps[()]
        cols = {n: np.asarray(arr[n]) for n in arr.dtype.names}
    else:
        names = [text(v) for v in ps["parameter_names"][()]]
        smp = ps["samples"][()]
        cols = {n: smp[:, i] for i, n in enumerate(names)}
    i = int(np.argmax(cols["log_likelihood"]))
    sample = {}
    for k, v in cols.items():                    # numeric columns only (some files carry text columns)
        if np.ndim(v) == 1 and np.issubdtype(np.asarray(v).dtype, np.number):
            sample[k] = float(v[i])
    cfg = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            cfg[name.split("/")[-1].replace("_", "-")] = text(obj[()])
    if "config_file" in g:
        g["config_file"].visititems(visit)
    psd = np.asarray(g["psds"]["H1"][()]) if "H1" in g["psds"] else np.asarray(
        g["psds"][list(g["psds"].keys())[0]][()])
    return sample, cfg, psd


def nums(s):
    return [float(v) for v in re.findall(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", s)]


def settings(cfg):
    """Waveform settings from bilby or LALInference configuration keys."""
    def num(*keys, default=None):
        for k in keys:
            v = cfg.get(k)
            if v is not None and str(v).strip() != "":
                try:
                    return float(v)
                except ValueError:
                    raise ValueError(f"config key '{k}' not numeric: {v!r}")
        if default is not None:
            return default
        raise ValueError(f"config keys {keys} missing or empty")
    flow_s = cfg.get("minimum-frequency") or cfg.get("flow") or "20"
    m = re.search(r"waveform\s*['\"]?\s*:\s*([0-9.]+)", flow_s)
    det_vals = [float(v) for k, v in re.findall(r"([HLVK]1)['\"]?\s*:\s*([0-9.]+)", flow_s)]
    f_an = min(det_vals) if det_vals else min(nums(flow_s))
    f_wf = float(m.group(1)) if m else num("fmin-template", default=f_an)
    fh_s = cfg.get("maximum-frequency") or cfg.get("fhigh") or ""
    fh_vals = [float(v) for k, v in re.findall(r"([HLVK]1)['\"]?\s*:\s*([0-9.]+)", fh_s)]
    srate = num("sampling-frequency", "srate")
    return {
        "f_ref": num("reference-frequency", "fref"),
        "f_wf": f_wf, "f_an": f_an,
        "f_high": min(fh_vals) if fh_vals else srate / 2,
        "srate": srate,
        "duration": num("duration", "seglen"),
        "wf_args": ast.literal_eval(cfg["waveform-arguments-dict"])
        if cfg.get("waveform-arguments-dict", "").strip().startswith("{") else {},
    }


def laldict(wf_args, unknown):
    d = lal.CreateDict()
    for k, v in wf_args.items():
        fn = getattr(ls, f"SimInspiralWaveformParamsInsert{k}", None)
        if fn is None:
            unknown.add(k)
        else:
            fn(d, int(v))
    return d


def physical(p, st):
    """Sampling parameters → masses (detector frame), Cartesian spins, inclination."""
    mc, q = p["chirp_mass"], p["mass_ratio"]
    m1 = mc * (1 + q) ** 0.2 / q ** 0.6
    m2 = m1 * q
    iota, s1x, s1y, s1z, s2x, s2y, s2z = ls.SimInspiralTransformPrecessingNewInitialConditions(
        p["theta_jn"], p["phi_jl"], p["tilt_1"], p["tilt_2"], p["phi_12"], p["a_1"], p["a_2"],
        m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, st["f_ref"], p["phase"])
    return m1, m2, (s1x, s1y, s1z, s2x, s2y, s2z), iota


class Generator:
    def __init__(self, approx, st, unknown):
        self.approx, self.st = approx, st
        self.dict = laldict(st["wf_args"], unknown)
        self.N = None
        self.count = 0

    def fd(self, p):
        self.count += 1
        st = self.st
        m1, m2, s, iota = physical(p, st)
        dist = p["luminosity_distance"] * 1e6 * lal.PC_SI
        if self.approx == "IMRPhenomXPHM":
            df = 1.0 / st["duration"]
            hp, hc = ls.SimInspiralChooseFDWaveform(
                m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s, dist, iota, p["phase"], 0, 0, 0,
                df, st["f_wf"], st["srate"] / 2, st["f_ref"], self.dict, ls.IMRPhenomXPHM)
            f = hp.f0 + df * np.arange(hp.data.length)
            return f, hp.data.data.copy()
        dt = 1.0 / st["srate"]
        if self.approx == "SEOBNRv4PHM":
            hp, hc = ls.SimInspiralChooseTDWaveform(
                m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s, dist, iota, p["phase"], 0, 0, 0,
                dt, st["f_wf"], st["f_ref"], self.dict, ls.SEOBNRv4PHM)
        else:
            from pyseobnr.generate_waveform import GenerateWaveform
            gp = {"mass1": m1, "mass2": m2, "spin1x": s[0], "spin1y": s[1], "spin1z": s[2],
                  "spin2x": s[3], "spin2y": s[4], "spin2z": s[5], "deltaT": dt,
                  "f22_start": st["f_wf"], "f_ref": st["f_ref"], "phi_ref": p["phase"],
                  "distance": p["luminosity_distance"], "inclination": iota,
                  "approximant": "SEOBNRv5PHM"}
            gp.update(st["wf_args"])            # the PE waveform arguments (e.g. lmax_nyquist)
            hp, hc = GenerateWaveform(gp).generate_td_polarizations()
        x = hp.data.data.copy()
        if self.N is None:
            self.N = int(2 ** np.ceil(np.log2(max(len(x) * 2, st["srate"] * st["duration"]))))
        if len(x) > self.N:          # never truncate silently
            raise ValueError(f"time-domain waveform longer than the FFT window ({len(x)} > {self.N})")
        X = np.fft.rfft(x, n=self.N) * dt
        f = np.fft.rfftfreq(self.N, dt)
        X *= np.exp(-2j * np.pi * f * float(hp.epoch))
        return f, X


def step_pair(p, name, eps):
    base, kind, lo, hi = STEPS[name]
    h = eps * (abs(p[name]) if kind == "rel" else 1.0)
    if kind == "abs" and np.isfinite(hi - lo):
        h = min(h, 0.25 * (hi - lo))
    up, dn = dict(p), dict(p)
    up[name] = p[name] + h
    dn[name] = p[name] - h
    if up[name] > hi or dn[name] < lo:                 # one-sided near a bound
        sgn = -1 if up[name] > hi else 1
        p1, p2 = dict(p), dict(p)
        p1[name] = p[name] + sgn * h
        p2[name] = p[name] + 2 * sgn * h
        return ("one-sided", sgn * h, p1, p2)
    return ("central", h, up, dn)


def derivative(gen, p, name, eps, h0):
    kind, h, a, b = step_pair(p, name, eps)
    if kind == "central":
        return (gen.fd(a)[1] - gen.fd(b)[1]) / (2 * h)
    return (-3 * h0 + 4 * gen.fd(a)[1] - gen.fd(b)[1]) / (2 * h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--out", default=os.path.expanduser("~/gw/f3"))
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    files = os.listdir(args.pe_dir)
    rows = []
    for ev, (sub, labels) in EVENTS.items():
        fn = next((x for x in files if sub in x), None)
        if fn is None:
            print(f"{ev}: PE file not found")
            continue
        with h5py.File(os.path.join(args.pe_dir, fn), "r") as f:
            for label in labels:
                t0 = time.time()
                sample, cfg, psd = read_label(f, label)
                st = settings(cfg)
                approx = ("IMRPhenomXPHM" if "XPHM" in label else
                          "SEOBNRv4PHM" if "v4PHM" in label else "SEOBNRv5PHM")
                fac = 1.0 if approx == "IMRPhenomXPHM" else EOB_STEP_FACTOR
                unknown = set()
                print(f"\n=== {ev}  {label}  ({approx}) ===")
                print(f"  f_ref {st['f_ref']}, waveform start {st['f_wf']} Hz, band {st['f_an']}–{st['f_high']} Hz, "
                      f"srate {st['srate']}, duration {st['duration']} s, waveform args {st['wf_args']}, "
                      f"step factor {fac:g}", flush=True)
                try:
                    gen = Generator(approx, st, unknown)
                    p = {k: sample[k] for k in STEPS}
                    fr, h0 = gen.fd(p)
                except Exception as e:
                    print(f"  GENERATION FAILED at the reference point: {type(e).__name__}: {str(e)[:300]}")
                    rows.append({"event": ev, "label": label, "parameter": "(reference)", "step": "",
                                 "r1": "", "r2": "", "norm_D": "", "weight": "", "stable": False,
                                 "class": "generation failed", "unknown_args": ""})
                    continue
                if unknown:
                    print(f"  UNKNOWN WAVEFORM ARGS (not passed to LAL): {sorted(unknown)}")
                band = (fr >= st["f_an"]) & (fr <= st["f_high"])
                Sn = np.interp(fr[band], psd[:, 0], psd[:, 1], left=np.inf, right=np.inf)
                w = 1.0 / Sn

                def norm(v):
                    return float(np.sqrt(np.sum(np.abs(v[band]) ** 2 * w)))
                print(f"  reference: Mc {p['chirp_mass']:.3f}, q {p['mass_ratio']:.3f}, a1 {p['a_1']:.2f}, "
                      f"a2 {p['a_2']:.2f}, D {p['luminosity_distance']:.0f} Mpc;  |h| (weighted) {norm(h0):.3e}")
                res = []
                for name, (base, kind, lo, hi) in STEPS.items():
                    try:
                        D = [derivative(gen, p, name, fac * base / 2 ** j, h0) for j in range(3)]
                    except Exception as e:
                        res.append((name, None, None, None, f"failed: {type(e).__name__}"))
                        continue
                    nD = [norm(d) for d in D]
                    r1 = norm(D[0] - D[1]) / max(nD[1], 1e-300)
                    r2 = norm(D[1] - D[2]) / max(nD[2], 1e-300)
                    res.append((name, nD[2], r1, r2, ""))
                dmax = max(r[1] for r in res if r[1] is not None and r[0] != "luminosity_distance")
                for name, nd, r1, r2, err in res:
                    if err:
                        cls, wgt = err, ""
                        print(f"  {name:20s} {err}")
                    else:
                        wgt = nd / dmax if name != "luminosity_distance" else float("nan")
                        if r2 < STABLE:
                            cls = "stable"
                        elif name != "luminosity_distance" and wgt < NEGLIGIBLE:
                            cls = "unstable but negligible"
                        else:
                            cls = "NOT STABLE"
                        print(f"  {name:20s} |D| {nd:.3e}  weight {wgt:8.2e}   r1 {r1:.2e}   r2 {r2:.2e}   {cls}")
                    rows.append({"event": ev, "label": label, "parameter": name,
                                 "step": fac * STEPS[name][0], "r1": r1, "r2": r2, "norm_D": nd,
                                 "weight": wgt, "stable": cls == "stable", "class": cls,
                                 "unknown_args": ";".join(sorted(unknown))})
                print(f"  waveforms generated: {gen.count}  ({time.time() - t0:.0f} s)", flush=True)
    with open(os.path.join(args.out, "f3_results.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print("\nsummary by model and class:")
    c = Counter(("XPHM" if "XPHM" in r["label"] else r["label"].split(":")[1], r["class"]) for r in rows)
    for k, v in sorted(c.items()):
        print(f"  {k[0]:14s} {k[1]:26s} {v}")
    bad = [r for r in rows if r["class"] not in ("stable", "unstable but negligible")]
    print(f"\n{len(rows)} entries, {len(bad)} problematic:")
    for r in bad:
        print(f"  {r['event']} {r['label']} {r['parameter']}: {r['class']}"
              + (f" (r2 = {r['r2']:.2e}, weight {r['weight']:.2e})" if isinstance(r['r2'], float) else ""))


if __name__ == "__main__":
    main()
