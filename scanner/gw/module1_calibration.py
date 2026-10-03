#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_calibration.py — detector calibration for the Module 1 estimator
(MODULE1_DESIGN §5 and §6, item 2; rule fixed in v0.9.8).

Source of the calibration model: the bilby priors of the PE label
(priors/.../recalib_<IFO>_amplitude_k, _phase_k: Gaussian(mu, sigma);
recalib_<IFO>_frequency_k: DeltaFunction(peak)). Labels without recalib priors
(the LALInference SEOBNRv4PHM labels) take them from the IMRPhenomXPHM label of the
same event (same detectors and strain). The plotting envelope priors/calibration/<IFO>
is NOT used: its median has the opposite sign convention in GWTC-4.0 (C00) labels
compared with GWTC-2.1/3 (C01) labels (check_calib_sigma.py, v0.9.8 record).

Model (bilby 2.8.2, bilby.gw.detector.calibration.CubicSpline, LIGO-T2300140):
  factor(f) = (1 + dA(f)) (2 + i dphi(f)) / (2 - i dphi(f)),
  dA, dphi = cubic splines in log10 f through the node values (linear in them).
Rule:
  - the reference waveform of each detector is multiplied by factor(f) evaluated at
    the prior means mu;
  - the spread enters J: amplitude node k -> sigma_A,k * B_k(f) * h,
    phase node k -> i * sigma_phi,k * B_k(f) * h  (first order of the factor),
    B_k = bilby spline response to a unit value at node k.
"""
import re

import h5py
import numpy as np
from bilby.gw.detector.calibration import CubicSpline

_NUM = r"([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)"


def _txt(x):
    x = x[()] if isinstance(x, h5py.Dataset) else x
    if isinstance(x, np.ndarray):
        x = x.ravel()[0] if x.size else ""
    return x.decode(errors="replace") if isinstance(x, (bytes, np.bytes_)) else str(x)


def _prior_texts(group):
    out = {}

    def visit(name, obj):
        base = name.split("/")[-1]
        if isinstance(obj, h5py.Dataset) and base.startswith("recalib_") and obj.dtype.kind in "SOU":
            out[base] = _txt(obj)
    if "priors" in group:
        group["priors"].visititems(visit)
    return out


def parse_priors(texts, ifo):
    """Node frequencies and Gaussian (mu, sigma) of amplitude and phase for one detector."""
    freqs, mu_a, sg_a, mu_p, sg_p = [], [], [], [], []
    k = 0
    while f"recalib_{ifo}_frequency_{k}" in texts:
        m = re.search(r"peak\s*=\s*" + _NUM, texts[f"recalib_{ifo}_frequency_{k}"])
        freqs.append(float(m.group(1)))
        for kind, mus, sgs in (("amplitude", mu_a, sg_a), ("phase", mu_p, sg_p)):
            t = texts[f"recalib_{ifo}_{kind}_{k}"]
            if not t.startswith("Gaussian"):
                raise ValueError(f"recalib_{ifo}_{kind}_{k}: not Gaussian: {t[:60]}")
            mus.append(float(re.search(r"mu\s*=\s*" + _NUM, t).group(1)))
            sgs.append(float(re.search(r"sigma\s*=\s*" + _NUM, t).group(1)))
        k += 1
    if k < 4:
        raise ValueError(f"{ifo}: {k} calibration nodes found (need >= 4)")
    freqs = np.array(freqs)
    expect = np.logspace(np.log10(freqs[0]), np.log10(freqs[-1]), k)
    if not np.allclose(freqs, expect, rtol=1e-6):
        raise ValueError(f"{ifo}: nodes are not log-spaced as bilby's CubicSpline assumes")
    return {"freqs": freqs, "mu_a": np.array(mu_a), "sg_a": np.array(sg_a),
            "mu_p": np.array(mu_p), "sg_p": np.array(sg_p)}


def load(fh, label, fallback_label, ifos):
    """Calibration priors per detector: from the label, else from fallback_label (XPHM)."""
    texts = _prior_texts(fh[label])
    source = label
    if not any(t.startswith("recalib_") for t in texts):
        texts = _prior_texts(fh[fallback_label])
        source = fallback_label
    cal = {ifo: parse_priors(texts, ifo) for ifo in ifos}
    return cal, source


class SplineBasis:
    """bilby CubicSpline for one detector, with its node -> f linear map."""

    def __init__(self, pri, ifo="X"):
        self.pri = pri
        self.n = len(pri["freqs"])
        self.spline = CubicSpline(f"recalib_{ifo}_", pri["freqs"][0], pri["freqs"][-1], self.n)

    def factor(self, f, amp, ph):
        params = {f"recalib_X_amplitude_{k}": amp[k] for k in range(self.n)}
        params.update({f"recalib_X_phase_{k}": ph[k] for k in range(self.n)})
        self.spline.prefix = "recalib_X_"
        return self.spline.get_calibration_factor(f, **params)

    def basis(self, f):
        """B[k, :] = dA(f) for a unit value at node k (all others 0)."""
        B = np.empty((self.n, len(f)))
        z = np.zeros(self.n)
        for k in range(self.n):
            e = z.copy()
            e[k] = 1.0
            B[k] = (self.factor(f, e, z) - 1.0).real
        return B

    def mean_factor(self, f):
        return self.factor(f, self.pri["mu_a"], self.pri["mu_p"])

    def directions(self, f, h):
        """Scaled first-order calibration directions (2n, len(f)) for waveform h."""
        B = self.basis(f)
        amp = self.pri["sg_a"][:, None] * B * h[None, :]
        ph = 1j * self.pri["sg_p"][:, None] * B * h[None, :]
        return np.concatenate([amp, ph], axis=0)


def _self_test():
    from scipy.interpolate import CubicSpline as SciCS
    rng = np.random.default_rng(1)
    pri = {"freqs": np.logspace(np.log10(20.0), np.log10(896.0), 10),
           "mu_a": rng.normal(0, 0.01, 10), "sg_a": np.full(10, 0.03),
           "mu_p": rng.normal(0, 0.01, 10), "sg_p": np.full(10, 0.03)}
    sb = SplineBasis(pri)
    f = np.linspace(20.0, 896.0, 3000)
    # 1. bilby spline == not-a-knot cubic spline in log10 f
    vals = rng.normal(0, 0.05, 10)
    ours = (sb.factor(f, vals, np.zeros(10)) - 1).real
    ref = SciCS(np.log10(pri["freqs"]), vals, bc_type="not-a-knot")(np.log10(f))
    print(f"1. bilby spline vs not-a-knot cubic spline in log10 f: max |diff| = {np.max(np.abs(ours - ref)):.2e}")
    # 2. linearity: basis reconstruction of dA
    B = sb.basis(f)
    print(f"2. basis reconstruction of dA: max |diff| = {np.max(np.abs(vals @ B - ours)):.2e}")
    # 3. first-order directions vs exact factor for a draw from the priors
    h = np.exp(1j * 0.01 * f)
    da = rng.normal(0, 1, 10)
    dp = rng.normal(0, 1, 10)
    J = sb.directions(f, h)
    lin = np.concatenate([da, dp]) @ J
    exact = sb.factor(f, pri["sg_a"] * da, pri["sg_p"] * dp) * h - h
    rel = np.linalg.norm(lin - exact) / np.linalg.norm(exact)
    print(f"3. first-order directions vs exact factor (1-sigma draw): relative difference {rel:.2e}"
          f"  (second order in 3%: expect ~1e-2)")
    # 4. prior text parsing
    texts = {}
    for k in range(10):
        texts[f"recalib_H1_frequency_{k}"] = f"DeltaFunction(peak={float(pri['freqs'][k])!r}, name='x', unit=None)"
        texts[f"recalib_H1_amplitude_{k}"] = f"Gaussian(mu={float(pri['mu_a'][k])!r}, sigma={float(pri['sg_a'][k])!r}, name='x')"
        texts[f"recalib_H1_phase_{k}"] = f"Gaussian(mu={float(pri['mu_p'][k])!r}, sigma={float(pri['sg_p'][k])!r}, name='x')"
    got = parse_priors(texts, "H1")
    ok = all(np.allclose(got[key], pri[key]) for key in pri)
    print(f"4. prior-text parsing round trip: {'ok' if ok else 'FAILED'}")


if __name__ == "__main__":
    _self_test()
