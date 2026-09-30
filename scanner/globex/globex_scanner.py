#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
globex_scanner.py — Module 0 (shallow-water laboratory waves) of the
residual-first scanner: the NEGATIVE-sign test.

Data: GLOBEX data base (Michallet et al., Zenodo doi:10.5281/zenodo.4009405,
CC BY 4.0). Scheldt flume, Deltares: 0.85 m still-water depth over a flat
section, 1:80 concrete beach with its toe at x = 16.57 m. Wave gauges at
128 Hz. The inshore trolley carried 11 gauges at 0.37 m and was placed at
9 positions (sessions); only gauges of the SAME session (a "block") are
recorded simultaneously and only within-block phase differences are used.

Measurement: for each block and frequency f, the cross-spectral phase of the
11 gauges relative to the central gauge is fitted with a straight line in x;
the slope gives the local wavenumber k at the block centre, with depth h
from the nominal bed profile.

Baseline (shallow water, NON-dispersive):
    ω = k √(g h) (1 + s_b)          s_b: free offset per block
                                    (absorbs depth/set-down error, mean current)
Residual:
    r = ω / (k √(g h)) − 1
Templates (one global amplitude a):
    r = s_b + a (kh)^q ,  q ∈ {1, 2, 3}
Known effect (linear theory, ω² = g k tanh kh):
    r = √(tanh(kh)/kh) − 1 ≈ −(kh)²/6 + 19 (kh)⁴/360
    → q = 2 must win with a < 0.

Modes:
  synth-null     linear random waves WITHOUT dispersion (calibration barrier)
  synth-inject   dispersion scale ×m, wrong shapes, systematics (reflection,
                 bound second harmonics)
  real           GLOBEX case (primary A3), search/confirm split in time;
                 includes a reflection diagnostic (pre-registered flag at R > 0.05)

Dependencies: numpy, scipy.
"""

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy import signal, stats

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

G = 9.81

# ----------------------------------------------------------------------------
# FIXED CONFIGURATION — its hash goes into the registration
# ----------------------------------------------------------------------------
CONFIG = {
    "version": "0.1.0",
    "geometry": {
        "h0": 0.85,                 # still-water depth over the flat section [m]
        "toe_x": 16.57,             # toe of the 1:80 slope [m]
        "slope": 1.0 / 80.0,
    },
    "blocks": {
        "starts": [82, 93, 104, 115],   # inshore-trolley sessions 0–3 (x ≈ 51–67 m,
                                        # before breaking); boundaries confirmed by
                                        # the 50 Hz coherence diagnostic
        "size": 11,
        "ref_index": 5,             # central gauge of the block
    },
    "spectrum": {
        "fs": 128.0,
        "nperseg": 32768,           # 256 s → Δf = 0.0039 Hz; the spectral-slope
                                    # (leakage) bias of the phase scales as 1/T²:
                                    # with 64 s it was +0.5 % in k at 0.45 Hz
        "overlap": 0.5,
        "freq_stride": 2,           # every other bin (Hann: neighbouring bins correlated)
    },
    "case": "A3",                   # primary case: JONSWAP γ = 20, Tp = 2.25 s, Hs = 0.10 m.
                                    # A1 (Tp = 1.58 s) reaches kh ≈ 0.5–1, where the
                                    # (kh)^4 term makes linear theory locally almost
                                    # linear in kh and the (kh)^2 shape cannot be
                                    # identified (synthetic calibration); A1 is used
                                    # only for the magnitude, not for the shape.
    "selection": {
        "band_hz": [0.33, 0.62],    # around the A3 peak (0.44 Hz), below 2 f_p
        "coh_min": 0.9,             # mean coherence of the block gauges with the reference
        "kh_max": 0.7,              # (kh)^2 dominates: quartic term ≤ 13 % of it
    },
    "synthetic": {"hs": 0.10, "fp": 1 / 2.25, "gamma": 20.0,
                  "reflector_x": 80.0},     # reflected waves originate at this point
    "reflection": {
        "flag_R": 0.05,             # median reflection coefficient in the band above
                                    # this → result flagged, no shape claim
    },
    "templates": [1, 2, 3],
    "kh_ref": 0.7,                  # where the effect size is measured
    "stats": {
        "alpha": 0.01,              # after Holm over templates
        "min_effect_rel": 5e-3,     # |a (kh_ref)^q| ≥ 0.5 %
        "barrier_test_alpha": 0.05,
    },
    "split": "first_half_search_second_half_confirm",
}


def config_hash(cfg=CONFIG):
    s = json.dumps(cfg, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(s).hexdigest()[:16]


def depth(x, cfg=CONFIG):
    g = cfg["geometry"]
    return np.where(x <= g["toe_x"], g["h0"], g["h0"] - (x - g["toe_x"]) * g["slope"])


# ----------------------------------------------------------------------------
# Dispersion models (used by the synthetic generator)
# ----------------------------------------------------------------------------
def r_linear(kh):
    """Relative deviation of full linear theory from the shallow-water speed."""
    kh = np.maximum(kh, 1e-9)
    return np.sqrt(np.tanh(kh) / kh) - 1.0


def solve_k(omega, h, r_func, iters=40):
    """Solve ω = k √(g h) (1 + r(kh)) for k (vectorised Newton, numeric derivative)."""
    omega, h = np.broadcast_arrays(np.asarray(omega, float), np.asarray(h, float))
    c0 = np.sqrt(G * h)
    k = omega / c0
    for _ in range(iters):
        F = k * c0 * (1 + r_func(k * h)) - omega
        dk = 1e-6 * np.maximum(k, 1e-6)
        dF = ((k + dk) * c0 * (1 + r_func((k + dk) * h)) - omega - F) / dk
        step = F / dF
        k = np.maximum(k - step, 1e-6)
        if np.max(np.abs(step / k)) < 1e-12:
            break
    return k


# ----------------------------------------------------------------------------
# Measurement
# ----------------------------------------------------------------------------
def measure_block(W, xb, cfg=CONFIG):
    """
    W: (time, 11) elevations of one block, xb: positions.
    Returns per-frequency k (from the phase slope), its standard error,
    the mean coherence and the block-centre depth.
    scipy.signal.csd(x, y) returns conj(X)·Y, so for an onshore-propagating wave
    cos(ωt − k x) the phase of csd(η_i, η_ref) is +k (x_i − x_ref): k = +slope.
    """
    sp, sel = cfg["spectrum"], cfg["selection"]
    ref = cfg["blocks"]["ref_index"]
    fs, nper = sp["fs"], sp["nperseg"]
    nover = int(nper * sp["overlap"])
    W = np.asarray(W, float) - np.asarray(W, float).mean(axis=0)
    f, Prr = signal.welch(W[:, ref], fs=fs, nperseg=nper, noverlap=nover)
    band = np.flatnonzero((f >= sel["band_hz"][0]) & (f <= sel["band_hz"][1]))
    band = band[::sp["freq_stride"]]
    n = W.shape[1]
    G_ = np.empty((n, len(band)), complex)
    coh = np.empty((n, len(band)))
    for i in range(n):
        _, Gir = signal.csd(W[:, i], W[:, ref], fs=fs, nperseg=nper, noverlap=nover)
        _, Pii = signal.welch(W[:, i], fs=fs, nperseg=nper, noverlap=nover)
        G_[i] = Gir[band]
        coh[i] = np.abs(Gir[band]) ** 2 / (Pii[band] * Prr[band] + 1e-300)
    order = np.argsort(xb)
    xs = xb[order] - xb[ref]
    phase = np.unwrap(np.angle(G_[order]), axis=0)
    phase -= phase[np.flatnonzero(order == ref)[0]]
    nseg = max(1, int((W.shape[0] - nper) / (nper - nover)) + 1)
    out = []
    for j, fi in enumerate(band):
        c = np.clip(coh[order, j], 1e-6, 1 - 1e-9)
        others = np.arange(n) != np.flatnonzero(order == ref)[0]
        mc = c[others].mean()
        var = (1 - c) / (2 * nseg * c)
        var[~others] = var[others].min()
        w = 1.0 / var
        A = np.column_stack([xs, np.ones_like(xs)])
        Aw = A * np.sqrt(w)[:, None]
        yw = phase[:, j] * np.sqrt(w)
        coef, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
        cov = np.linalg.inv(Aw.T @ Aw)
        k = coef[0]
        out.append({"f": f[fi], "k": k, "se_k": float(np.sqrt(cov[0, 0])), "coh": mc})
    return out


def reflection_block(W, xb, k_of_f, cfg=CONFIG):
    """
    Reflection coefficient from the standing-wave modulation of the spectral
    density along x: P_i(f) = c0 + c_t (x_i − x̄) + c1 cos(2 k x_i) + c2 sin(2 k x_i),
    R ≈ √(c1² + c2²) / (2 c0). k_of_f: dict f → measured k.
    """
    sp = cfg["spectrum"]
    fs, nper = sp["fs"], sp["nperseg"]
    nover = int(nper * sp["overlap"])
    W = np.asarray(W, float) - np.asarray(W, float).mean(axis=0)
    P = []
    for i in range(W.shape[1]):
        f, Pii = signal.welch(W[:, i], fs=fs, nperseg=nper, noverlap=nover)
        P.append(Pii)
    P = np.array(P)                                           # (gauge, f)
    xm = xb - xb.mean()
    R = {}
    for fi, k in k_of_f.items():
        j = int(np.argmin(np.abs(f - fi)))
        A = np.column_stack([np.ones_like(xm), xm, np.cos(2 * k * xb), np.sin(2 * k * xb)])
        c, *_ = np.linalg.lstsq(A, P[:, j], rcond=None)
        R[fi] = float(np.hypot(c[2], c[3]) / (2 * max(c[0], 1e-30)))
    return R


def collect(W, x, cfg=CONFIG):
    """All blocks → arrays of (block, ω, k, se, h, kh, r, se_r) after selection."""
    b = cfg["blocks"]
    sel = cfg["selection"]
    rows = []
    refl = []
    for bi, s in enumerate(b["starts"]):
        cols = np.arange(s, s + b["size"])
        xb = x[cols]
        hc = float(depth(xb).mean())
        ms = measure_block(W[:, cols], xb, cfg)
        Rb = reflection_block(W[:, cols], xb, {m["f"]: m["k"] for m in ms if m["k"] > 0}, cfg)
        refl.append(float(np.median(list(Rb.values()))) if Rb else float("nan"))
        for m in ms:
            if m["k"] <= 0 or m["coh"] < sel["coh_min"]:
                continue
            kh = m["k"] * hc
            if kh > sel["kh_max"]:
                continue
            om = 2 * np.pi * m["f"]
            r = om / (m["k"] * np.sqrt(G * hc)) - 1
            se_r = (1 + r) * m["se_k"] / m["k"]
            rows.append((bi, om, m["k"], m["se_k"], hc, kh, r, se_r, m["coh"]))
    a = np.array(rows) if rows else np.empty((0, 9))
    keys = ["block", "omega", "k", "se_k", "h", "kh", "r", "se_r", "coh"]
    d = {k_: a[:, i] for i, k_ in enumerate(keys)}
    d["block"] = d["block"].astype(int)
    d["reflection_median"] = refl
    return d


# ----------------------------------------------------------------------------
# Test (linear weighted least squares)
# ----------------------------------------------------------------------------
def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (m - rank) * p[i])
        adj[i] = min(1.0, run)
    return adj


def _design(d, qs):
    blocks = np.unique(d["block"])
    X = [(d["block"] == b).astype(float) for b in blocks]
    X += [d["kh"] ** q for q in qs]
    return np.column_stack(X)


def _wls(d, qs):
    X = _design(d, qs)
    w = 1.0 / np.maximum(d["se_r"], 1e-12)
    c, *_ = np.linalg.lstsq(X * w[:, None], d["r"] * w, rcond=None)
    rss = float(np.sum(((d["r"] - X @ c) * w) ** 2))
    return c, rss, X.shape[1]


def _ftest(rss_small, rss_big, dof):
    if rss_big <= 0 or dof <= 0:
        return np.inf, 0.0
    F = max(0.0, (rss_small - rss_big) / (rss_big / dof))
    return float(F), float(stats.f.sf(F, 1, dof))


def run_test(d, cfg=CONFIG):
    templates = cfg["templates"]
    N = len(d["r"])
    if N < 10:
        return {"N": N, "error": "too few points", "significant": False, "any_flag": False,
                "winner": None, "shape_unique": False, "per": {}}
    _, rss0, p0 = _wls(d, [])
    per = {}
    for q in templates:
        c, rss1, p1 = _wls(d, [q])
        F, p = _ftest(rss0, rss1, N - p1)
        per[q] = {"a": float(c[-1]), "F": F, "p": p, "rss": rss1,
                  "effect": float(c[-1] * cfg["kh_ref"] ** q)}
    padj = holm([per[q]["p"] for q in templates])
    alpha, mine = cfg["stats"]["alpha"], cfg["stats"]["min_effect_rel"]
    for q, pa in zip(templates, padj):
        per[q]["padj"] = float(pa)
        per[q]["flag"] = bool(pa < alpha and abs(per[q]["effect"]) >= mine)
    winner = min(templates, key=lambda q: per[q]["rss"])
    others = [q for q in templates if q != winner]
    part = {}
    for q in others:
        _, rss_both, pb = _wls(d, [winner, q])
        F, p = _ftest(per[q]["rss"], rss_both, N - pb)
        part[q] = {"F_winner_given_q": F, "p": p}
    padj2 = holm([part[q]["p"] for q in others])
    for q, pa in zip(others, padj2):
        part[q]["padj"] = float(pa)
    unique = all(part[q]["padj"] < alpha for q in others)
    offsets, _, _ = _wls(d, [])
    return {"N": N, "per": per, "winner": winner, "significant": per[winner]["flag"],
            "shape_unique": unique, "partial": part,
            "any_flag": any(per[q]["flag"] for q in templates),
            "block_offsets": [float(v) for v in offsets],
            "kh_range": [float(d["kh"].min()), float(d["kh"].max())]}


def expected_a(d, m=1.0):
    """Amplitude of the q = 2 template fitted to NOISELESS linear theory
    (dispersion scale m) at the same kh points and weights — the reference
    value for the recovered amplitude."""
    dd = dict(d)
    dd["r"] = m * r_linear(d["kh"])
    c, _, _ = _wls(dd, [2])
    return float(c[-1])


def declared(r):
    if not r["significant"]:
        return "none"
    sign = "+" if r["per"][r["winner"]]["a"] > 0 else "-"
    return f"k^{r['winner']}{sign}" + (" distinguishable" if r["shape_unique"] else " undetermined")


def print_test(r, label):
    if "error" in r:
        print(f"  [{label}] {r['error']} (N = {r['N']})")
        return
    print(f"\n  [{label}] points {r['N']}, kh {r['kh_range'][0]:.2f}–{r['kh_range'][1]:.2f}, "
          f"block offsets " + ", ".join(f"{v:+.4f}" for v in r["block_offsets"]))
    for q, v in r["per"].items():
        mark = "  ← significant" if v["flag"] else ""
        print(f"     template (kh)^{q}: a = {v['a']:+.4f}  effect at kh_ref {v['effect']:+.4f}  "
              f"p(Holm) = {v['padj']:.2e}{mark}")
    print(f"     declared: {declared(r)}")


# ----------------------------------------------------------------------------
# Synthetic records (linear random waves on the GLOBEX beach, WKB shoaling)
# ----------------------------------------------------------------------------
def jonswap(f, hs=0.10, fp=1 / 1.58, gamma=3.3):
    sig = np.where(f <= fp, 0.07, 0.09)
    S = f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * gamma ** np.exp(-(f - fp) ** 2 / (2 * sig ** 2 * fp ** 2))
    df = f[1] - f[0]
    S *= (hs / 4) ** 2 / np.sum(S * df)
    return S


def synth_record(rng, xg, duration, fs=128.0, r_func=r_linear, reflection=0.0,
                 bound=0.0, noise=5e-4, hs=None, fp=None, gamma=None, f_lo=0.2, f_hi=2.0,
                 cfg=CONFIG):
    """
    Linear random waves generated at the slope toe, propagated with the
    dispersion ω = k √(g h)(1 + r(kh)) and WKB shoaling (energy flux).
    Options: a reflected wave (coefficient R, reflected at x_r and travelling
    offshore, phase-locked to the incoming wave) and
    phase-locked second harmonics (Stokes self-interaction, scale 'bound').
    Returns (time, gauges).
    """
    n = int(duration * fs)
    f = np.fft.rfftfreq(n, 1 / fs)
    band = (f >= f_lo) & (f <= f_hi)
    fb = f[band]
    om = 2 * np.pi * fb
    df = f[1] - f[0]
    sc = cfg["synthetic"]
    hs = sc["hs"] if hs is None else hs
    fp = sc["fp"] if fp is None else fp
    gamma = sc["gamma"] if gamma is None else gamma
    amp0 = np.sqrt(2 * jonswap(fb, hs, fp, gamma) * df)
    phi = rng.uniform(0, 2 * np.pi, len(fb))
    x0 = cfg["geometry"]["toe_x"]
    x_r = cfg["synthetic"]["reflector_x"]
    xgrid = np.linspace(x0, max(max(xg) + 0.5, x_r), 3000)
    hgrid = depth(xgrid, cfg)
    K = solve_k(om[None, :], hgrid[:, None], r_func)                # (x, f)
    dx = xgrid[1] - xgrid[0]
    phase_int = np.concatenate([np.zeros((1, len(fb))),
                                np.cumsum(0.5 * (K[1:] + K[:-1]) * dx, axis=0)])
    # group velocity by finite difference in ω
    dom = 1e-4
    K2 = solve_k((om + dom)[None, :], hgrid[:, None], r_func)
    cg = dom / np.maximum(K2 - K, 1e-12)
    shoal = np.sqrt(cg[0][None, :] / cg)
    theta_r = np.array([np.interp(x_r, xgrid, phase_int[:, j]) for j in range(len(fb))])
    out = np.empty((n, len(xg)), np.float32)
    for gi, xgi in enumerate(xg):
        # exact position: interpolate phase, shoaling and k between grid points
        th = np.array([np.interp(xgi, xgrid, phase_int[:, j]) for j in range(len(fb))])
        a = amp0 * np.array([np.interp(xgi, xgrid, shoal[:, j]) for j in range(len(fb))])
        kk_x = np.array([np.interp(xgi, xgrid, K[:, j]) for j in range(len(fb))])
        hh_x = float(depth(np.array([xgi]), cfg)[0])
        spec = np.zeros(len(f), complex)
        spec[band] = a * np.exp(1j * (phi - th))                     # incoming
        if reflection > 0:
            # reflected at x_r: travelled to x_r and back → phase φ − 2θ(x_r) + θ(x)
            spec[band] += reflection * a * np.exp(1j * (phi - 2 * theta_r + th))
        if bound > 0:
            kk, hh = kk_x, hh_x
            coef = kk * np.cosh(kk * hh) * (2 + np.cosh(2 * kk * hh)) / (4 * np.sinh(kk * hh) ** 3)
            a2 = bound * coef * a ** 2
            idx2 = np.searchsorted(f, 2 * fb)
            ok = idx2 < len(f)
            np.add.at(spec, idx2[ok], (a2 * np.exp(2j * (phi - th)))[ok])
        out[:, gi] = np.fft.irfft(spec * n / 2, n)
    out += rng.normal(0, noise, out.shape).astype(np.float32)
    return out


def synth_gauges():
    """Gauge geometry of the analysed blocks: 11 gauges at 0.37 m per block,
    blocks starting at the GLOBEX positions of columns 82, 93, 104, 115."""
    starts_x = [50.986, 55.066, 59.137, 63.207]
    x = np.full(181, np.nan)
    for s, x0 in zip(CONFIG["blocks"]["starts"], starts_x):
        x[s:s + 11] = x0 + 0.37 * np.arange(11)
    return x


def synth_measure(rng, args, **kw):
    x = synth_gauges()
    cols = np.concatenate([np.arange(s, s + 11) for s in CONFIG["blocks"]["starts"]])
    W = np.zeros((int(args.duration * 128), 181), np.float32)
    W[:, cols] = synth_record(rng, x[cols], args.duration, **kw)
    return collect(W, x)


def barrier(flags):
    M, k = len(flags), int(sum(flags))
    pb = float(stats.binom.sf(k - 1, M, CONFIG["stats"]["alpha"])) if k > 0 else 1.0
    return k / max(M, 1), pb, pb >= CONFIG["stats"]["barrier_test_alpha"]


def r_power(q, a):
    return lambda kh: a * np.maximum(kh, 1e-9) ** q


def r_scaled(m):
    return lambda kh: m * r_linear(kh)


# ----------------------------------------------------------------------------
# Modes
# ----------------------------------------------------------------------------
def run_synth_null(args):
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    print(f"config hash: {config_hash()}  (seed {args.seed})")
    flags, rows = [], []
    for i in range(args.m):
        t0 = time.time()
        d = synth_measure(rng, args, r_func=lambda kh: 0.0 * kh)       # non-dispersive
        r = run_test(d)
        flags.append(r["any_flag"])
        rows.append(r)
        print(f"  {i + 1:3d}/{args.m}: points {r['N']:3d}  flag {r['any_flag']!s:5}  "
              f"{declared(r)}  ({time.time() - t0:.0f} s)", flush=True)
    fpr, pb, ok = barrier(flags)
    print(f"null: {sum(flags)}/{len(flags)}, FPR {fpr:.3f} (target ≤ {CONFIG['stats']['alpha']}), "
          f"binomial p = {pb:.3f} → {'PASSED' if ok else 'NOT PASSED → STOP'}")
    with open(os.path.join(args.out, "synth_null.json"), "w", encoding="utf-8") as fh:
        json.dump({"config_hash": config_hash(), "args": vars(args), "runs": rows}, fh,
                  indent=1, default=float)
    return 0 if ok else 2


def run_synth_inject(args):
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    print(f"config hash: {config_hash()}  (seed {args.seed})")
    a_ref = -1.0 / 6.0 * CONFIG["kh_ref"] ** 2        # leading-order effect at kh_ref
    cases = []
    if "scale" in args.parts:
        cases += [(f"dispersion ×{m:g}", "k^2-", dict(r_func=r_scaled(m)), m) for m in args.mults]
    if "wrong" in args.parts:
        cases += [(f"wrong shape (kh)^{q}", f"k^{q}-",
                   dict(r_func=r_power(q, a_ref / CONFIG['kh_ref'] ** q)), None) for q in (1, 3)]
    if "syst" in args.parts:
        cases += [("null + reflection R=0.03", None, dict(r_func=lambda kh: 0 * kh, reflection=0.03), 0.0),
                  ("null + reflection R=0.1", None, dict(r_func=lambda kh: 0 * kh, reflection=0.1), 0.0),
                  ("null + bound harmonics ×1", None, dict(r_func=lambda kh: 0 * kh, bound=1.0), 0.0),
                  ("null + bound harmonics ×3", None, dict(r_func=lambda kh: 0 * kh, bound=3.0), 0.0),
                  ("dispersion ×1 + reflection R=0.03", "k^2-", dict(r_func=r_linear, reflection=0.03), 1.0),
                  ("dispersion ×1 + reflection R=0.1", "k^2-", dict(r_func=r_linear, reflection=0.1), 1.0),
                  ("dispersion ×1 + bound ×3", "k^2-", dict(r_func=r_linear, bound=3.0), 1.0)]
    table = []
    from collections import Counter
    for label, truth, kw, mult in cases:
        decl, ratios, refls = [], [], []
        for _ in range(args.m):
            d = synth_measure(rng, args, **kw)
            r = run_test(d)
            decl.append(declared(r))
            refls.append(np.nanmedian(d["reflection_median"]))
            if mult and "per" in r and 2 in r["per"]:
                ratios.append(r["per"][2]["a"] / expected_a(d, mult))
        n = len(decl)
        c = Counter(decl)
        det = sum(v for k_, v in c.items() if k_ != "none") / n
        ident = c.get(f"{truth} distinguishable", 0) / n if truth else float("nan")
        false = (sum(v for k_, v in c.items() if k_.endswith("distinguishable") and not k_.startswith(truth))
                 / n) if truth else sum(v for k_, v in c.items() if k_ != "none") / n
        print(f"  {label}: " + ", ".join(f"{k_}: {v / n:.2f}" for k_, v in sorted(c.items())))
        print(f"     detection {det:.2f} | form_identification {ident:.2f} | false_form {false:.2f}"
              + (f" | median a/a_expected {np.median(ratios):.3f}" if ratios else "")
              + f" | estimated R {np.median(refls):.3f}", flush=True)
        table.append({"case": label, "counts": dict(c), "detection_rate": det,
                      "form_identification_rate": ident, "false_form_rate": false,
                      "a_ratio_median": float(np.median(ratios)) if ratios else None})
        with open(os.path.join(args.out, "synth_inject.json"), "w", encoding="utf-8") as fh:
            json.dump({"config_hash": config_hash(), "args": vars(args), "table": table}, fh,
                      indent=1)


def run_real(args):
    os.makedirs(args.out, exist_ok=True)
    print(f"config hash: {config_hash()}")
    folder = args.data
    case = os.path.basename(os.path.normpath(folder))
    W = np.load(os.path.join(folder, f"WG_{case}.npy"), mmap_mode="r")
    x = np.load(os.path.join(folder, f"x_{case}.npy"))
    n = W.shape[0]
    halves = {"search": slice(0, n // 2), "confirm": slice(n // 2, 2 * (n // 2))}
    res = {}
    for name, sl in halves.items():
        d = collect(np.asarray(W[sl]), x)
        r = run_test(d)
        r["reflection_median_per_block"] = d["reflection_median"]
        r["reflection_flag"] = bool(np.nanmedian(d["reflection_median"]) > CONFIG["reflection"]["flag_R"])
        print(f"\n  [{name}] reflection (median per block): "
              + ", ".join(f"{v:.3f}" for v in d["reflection_median"])
              + ("  → FLAG (R > %.2f)" % CONFIG["reflection"]["flag_R"] if r["reflection_flag"] else ""))
        r["a_expected_q2"] = expected_a(d) if r["N"] >= 10 else None
        print_test(r, name)
        if r.get("a_expected_q2") is not None and 2 in r["per"]:
            print(f"     [interpretation] a(q=2) = {r['per'][2]['a']:+.4f}, linear theory at the same "
                  f"points {r['a_expected_q2']:+.4f}, ratio {r['per'][2]['a'] / r['a_expected_q2']:.3f}")
        res[name] = r

    def ok(r):
        return r["significant"] and r["winner"] == 2 and r["per"][2]["a"] < 0
    found = ok(res["search"])
    conf = found and ok(res["confirm"])
    rflag = res["search"]["reflection_flag"] or res["confirm"]["reflection_flag"]
    shape = (conf and res["search"]["shape_unique"] and res["confirm"]["shape_unique"]
             and not rflag)
    print("\n=== RESULT ===")
    print("finite-depth term (kh)^2, a < 0: " + ("found and confirmed" if conf
                                                   else "found but not confirmed" if found else "not found"))
    if rflag:
        print("   REFLECTION FLAG: median R above the pre-registered limit → no shape claim")
    if conf:
        print("   the (kh)^2 shape is " + ("DISTINGUISHABLE in both halves" if shape
                                           else "NOT DISTINGUISHABLE (or reflection-flagged) → reported as "
                                                "'negative deviation, undetermined shape'"))
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"config_hash": config_hash(), "case": case, "results": res,
                   "found": found, "confirmed": conf, "shape_unique": shape,
                   "reflection_flag": rflag},
                  fh, indent=1, default=float)


def main():
    ap = argparse.ArgumentParser(description="Module 0 (GLOBEX shallow-water waves) of the residual-first scanner")
    sub = ap.add_subparsers(dest="mode", required=True)
    for name in ("synth-null", "synth-inject"):
        p = sub.add_parser(name)
        p.add_argument("--m", type=int, default=(200 if name == "synth-null" else 30),
                       help="number of realizations")
        p.add_argument("--duration", type=float, default=2070.0,
                       help="synthetic record length [s] (= one half of a GLOBEX A-series record)")
        p.add_argument("--seed", type=int, default=1)
        p.add_argument("--out", default=f"results/{name.replace('-', '_')}")
        if name == "synth-inject":
            p.add_argument("--mults", type=float, nargs="+", default=[0.0, 0.25, 0.5, 1.0, 2.0])
            p.add_argument("--parts", nargs="+", default=["scale", "wrong", "syst"],
                           choices=["scale", "wrong", "syst"])
    r = sub.add_parser("real")
    r.add_argument("--data", required=True, help="folder with WG_<case>.npy and x_<case>.npy")
    r.add_argument("--out", default="results/real")
    args = ap.parse_args()
    if args.mode == "synth-null":
        return run_synth_null(args)
    if args.mode == "synth-inject":
        return run_synth_inject(args)
    return run_real(args)


if __name__ == "__main__":
    sys.exit(main())
