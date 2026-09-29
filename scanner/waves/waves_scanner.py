#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
waves_scanner.py — Module 0 (sea waves, stereo field η(x, y, t)) of the residual-first scanner
=================================================================================

Data: IFREMER stereo data set (Guimarães et al. 2020, CC-BY), e.g.
  BS_2011/2011-10-04_11-38-00_12Hz/nc/Surfaces_20111004_113800_short.nc
  grid 0.05 m, 12 Hz, 20 min; working band 0.6–3 Hz (deep water, kh > 40).

Measurement: for every cell k = (k_a, k_b) of the 3-D spectrum S(f, k_a, k_b)
the ridge frequency ω_meas(k) is found (maximum over f + parabolic refinement).

Baseline model (parameters fitted jointly, nonlinear least squares):
  ω = √(g_eff · k · (1 + C k²)) + k·U
    g_eff  — free (absorbs any scale error of the stereo calibration)
    U      — free vector (current / Doppler)
    C      — capillary coefficient ℓ_c² = γ/(ρg) (≈ 7.2·10⁻⁶ m²)

Residual templates (relative effect a at k_ref):
  δω = a · √(g k) · (k / k_ref)^q
  q = 0 is absorbed by g_eff and is not tested.
  q = 2 is the shape of the capillary term and of the Λ model (ω² = gk(1 + Λk²)).

Active test (v0.5.x):
  TEST B (known effect): baseline without C; templates q ∈ {1, 2, 3}.
                         Expectation: q = 2 with a > 0 (surface tension).
  The null calibration is synthetic; the real data are NOT a null test.
  Test A (baseline with free C) was retired in v0.5.0 — reason in CONFIG["retired"].
Split: first half of the record — search, second half — confirmation.

Modes:
  synth-null     synthetic fields (calibration barrier for test B)
  synth-inject   synthetic fields: scale of capillary amplitudes + wrong shapes
  real           real NetCDF file

Dependencies: numpy, scipy; for real: xarray, netCDF4.
"""

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy import fft as sfft
from scipy import optimize, stats

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

G = 9.81
RHO = 1025.0
GAMMA_NOMINAL = 0.0727                       # N/m, clean sea water ~20 °C
C_NOMINAL = GAMMA_NOMINAL / (RHO * G)        # ℓ_c² ≈ 7,23e-6 m²

# ----------------------------------------------------------------------------
# FIXED CONFIGURATION — its hash goes into the registration
# ----------------------------------------------------------------------------
CONFIG = {
    "version": "0.5.1",
    "note": "English translation of v0.5.0 (config hash 069ba514c71704e5); no functional change.",
    "spectrum": {
        "seg_frames": 512,          # 42.7 s at 12 Hz → Δf = 0.0234 Hz
        "time_window": "hann",
        "space_window": "hann",
        "prewhiten": "laplacian",   # power × ~k⁴: flattens the k⁻⁴ slope and
                                    # removes the frequency bias from leakage in k
    },
    "ridge": {
        "band_hz": [1.2, 3.0],      # cells with k from 5.8 to 36 rad/m; below 1.2 Hz
                                    # the k step exceeds 8 % of k and leakage
                                    # biases the frequency (see calibration)
        "search_hz": [0.2, 5.8],    # peak search range in f (Nyquist 6 Hz)
        "snr_db": 10.0,             # peak power / median over f ≥ 10 dB
        "cell_stride": 2,           # only cells with even indices (reduces the
                                    # correlation introduced by the spatial window)
        "weights": "inv_cg",        # weight ∝ 1/c_g ∝ √k: the ω error comes from
                                    # the cell width, σ ≈ c_g·Δk (calibrated)
    },
    "templates_B": [1, 2, 3],
    "retired": {
        "test_A": ("Retired in v0.5.0 as a model-design decision made BEFORE the "
                   "new B barrier series. Reason: with a free capillary coefficient C, "
                   "g_eff and the k^1/k^3 templates are nearly collinear for "
                   "k = 5.8–36 rad/m; the systematic ridge bias (~1e-4) is amplified "
                   "to |a(k^1)| ~ 1e-3–4e-3. Observed v0.4.0 failures: "
                   "2/40 (seed 1) and 8/200 (seed 1001, p = 0.001)."),
    },
    "k_ref": 20.0,                  # rad/m, where the effect size is measured
    "stats": {
        "alpha": 0.01,              # after Holm over templates
        "min_effect_rel": 5e-4,     # |δω/ω| at k_ref ≥ 0.05 %: systematic floor
                                    # (ridge bias up to ~8e-4 per ring;
                                    # 99th percentile of |a(k^1)| in B nulls: 4.7e-4)
        "barrier_test_alpha": 0.05,
    },
    "split": "first_half_search_second_half_confirm",
}


def config_hash(cfg=CONFIG):
    s = json.dumps(cfg, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(s).hexdigest()[:16]


def k_deep(f_hz):
    return (2 * np.pi * f_hz) ** 2 / G


# ----------------------------------------------------------------------------
# Spectrum
# ----------------------------------------------------------------------------
class SpectrumAccumulator:
    """Averages |FFT|² over segments (Welch without overlap).
    Index [f, a, b]; f ≥ 0. A wave cos(k·x − ωt) with ω > 0 appears at
    spatial index −k, so the propagation vector is −(k_a, k_b)."""

    def __init__(self, nt, n0, n1, dt, dx0, dx1, prewhiten=None):
        self.prewhiten = CONFIG["spectrum"]["prewhiten"] if prewhiten is None else prewhiten
        if self.prewhiten == "laplacian":
            n0, n1 = n0 - 2, n1 - 2
        self.nt, self.n0, self.n1 = nt, n0, n1
        self.wt = np.hanning(nt).astype(np.float32)[:, None, None]
        self.ws = np.outer(np.hanning(n0), np.hanning(n1)).astype(np.float32)[None]
        self.f = sfft.rfftfreq(nt, dt)
        self.ka = 2 * np.pi * sfft.fftfreq(n0, dx0)
        self.kb = 2 * np.pi * sfft.fftfreq(n1, dx1)
        self.P = None
        self.count = 0

    def add(self, z):
        z = np.asarray(z, dtype=np.float32)
        z = z - z.mean(axis=0, keepdims=True)
        if self.prewhiten == "laplacian":
            z = (z[:, 2:, 1:-1] + z[:, :-2, 1:-1] + z[:, 1:-1, 2:] + z[:, 1:-1, :-2]
                 - 4 * z[:, 1:-1, 1:-1])
        F = sfft.rfft(z * self.wt * self.ws, axis=0, workers=-1)
        F = sfft.fft2(F, axes=(1, 2), workers=-1)
        p = (F.real ** 2 + F.imag ** 2).astype(np.float64)
        self.P = p if self.P is None else self.P + p
        self.count += 1

    def mean(self):
        return self.P / self.count


def extract_ridge(P, f, ka, kb, cfg=CONFIG):
    """Returns arrays: kx, ky (propagation vector), kmag, ω_meas, snr."""
    r = cfg["ridge"]
    KA, KB = np.meshgrid(ka, kb, indexing="ij")
    kx, ky = -KA, -KB
    kmag = np.hypot(kx, ky)
    k_lo, k_hi = k_deep(r["band_hz"][0]), k_deep(r["band_hz"][1])
    sel = (kmag >= k_lo) & (kmag <= k_hi)
    s = r["cell_stride"]
    if s > 1:
        ia = np.arange(len(ka))[:, None] % s == 0
        ib = np.arange(len(kb))[None, :] % s == 0
        sel &= ia & ib
    i0 = int(np.searchsorted(f, r["search_hz"][0]))
    i1 = int(np.searchsorted(f, r["search_hz"][1]))
    Pc = P[i0:i1][:, sel]                                 # (nf, ncells)
    j = np.argmax(Pc, axis=0)
    med = np.median(Pc, axis=0) + 1e-300
    pk = Pc[j, np.arange(Pc.shape[1])]
    snr = pk / med
    ok = (j > 0) & (j < Pc.shape[0] - 1) & (10 * np.log10(snr) >= r["snr_db"])
    cols = np.arange(Pc.shape[1])[ok]
    jj = j[ok]
    la = np.log(Pc[jj - 1, cols] + 1e-300)
    lb = np.log(Pc[jj, cols] + 1e-300)
    lc = np.log(Pc[jj + 1, cols] + 1e-300)
    den = la - 2 * lb + lc
    delta = np.where(den < 0, 0.5 * (la - lc) / np.where(den < 0, den, -1), 0.0)
    df = f[1] - f[0]
    f_pk = f[i0 + jj] + delta * df
    return {
        "kx": kx[sel][ok], "ky": ky[sel][ok], "k": kmag[sel][ok],
        "omega": 2 * np.pi * f_pk, "snr": snr[ok],
        "n_candidates": int(sel.sum()), "df_hz": float(df),
    }


# ----------------------------------------------------------------------------
# Models and tests
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


def _model(theta, d, with_C, qs, k_ref):
    g, ux, uy = theta[0], theta[1], theta[2]
    i = 3
    C = 0.0
    if with_C:
        C = theta[i]
        i += 1
    om = np.sqrt(np.maximum(g * d["k"] * (1 + C * d["k"] ** 2), 1e-12)) + d["kx"] * ux + d["ky"] * uy
    for q in qs:
        om = om + theta[i] * np.sqrt(G * d["k"]) * (d["k"] / k_ref) ** q
        i += 1
    return om


def _fit(d, w, with_C, qs, k_ref):
    th0 = [G, 0.0, 0.0] + ([C_NOMINAL] if with_C else []) + [0.0] * len(qs)
    res = optimize.least_squares(
        lambda th: (_model(th, d, with_C, qs, k_ref) - d["omega"]) * w, th0,
        x_scale="jac", method="lm", max_nfev=4000)
    return res.x, float(np.sum(res.fun ** 2))


def _ftest(rss_small, rss_big, dof):
    if rss_big <= 0:
        return np.inf, 0.0
    F = max(0.0, (rss_small - rss_big) / (rss_big / dof))
    return F, float(stats.f.sf(F, 1, dof))


def run_test(d, test, cfg=CONFIG):
    """test = 'B'. Returns a dictionary with the results."""
    if test != "B":
        raise ValueError("only test B is active since v0.5.0 (test A retired, see CONFIG['retired'])")
    k_ref = cfg["k_ref"]
    with_C = False
    templates = cfg["templates_B"]
    wt = cfg["ridge"]["weights"]
    w = (np.sqrt(d["k"]) if wt == "inv_cg" else
         np.sqrt(d["snr"]) if wt == "sqrt_snr" else np.ones_like(d["k"]))
    w = w / np.median(w)
    N = len(d["k"])
    th0, rss0 = _fit(d, w, with_C, [], k_ref)
    n_base = len(th0)
    per = {}
    for q in templates:
        th, rss1 = _fit(d, w, with_C, [q], k_ref)
        F, p = _ftest(rss0, rss1, N - n_base - 1)
        per[q] = {"a": float(th[-1]), "F": float(F), "p": p, "rss": rss1,
                  "theta": [float(x) for x in th]}
    padj = holm([per[q]["p"] for q in templates])
    alpha, mine = cfg["stats"]["alpha"], cfg["stats"]["min_effect_rel"]
    for q, pa in zip(templates, padj):
        per[q]["padj"] = float(pa)
        per[q]["flag"] = bool(pa < alpha and abs(per[q]["a"]) >= mine)
    winner = min(templates, key=lambda q: per[q]["rss"])
    # Is the shape distinguishable? The winner must be needed even in the presence
    # of every other template (partial F-test, Holm over the alternatives).
    others = [q for q in templates if q != winner]
    part = {}
    for q in others:
        _, rss_both = _fit(d, w, with_C, [winner, q], k_ref)
        F, p = _ftest(per[q]["rss"], rss_both, N - n_base - 2)
        part[q] = {"F_winner_given_q": float(F), "p": p}
    padj2 = holm([part[q]["p"] for q in others]) if others else []
    for q, pa in zip(others, padj2):
        part[q]["padj"] = float(pa)
    shape_unique = all(part[q]["padj"] < alpha for q in others)
    base = {"g_eff": float(th0[0]), "Ux": float(th0[1]), "Uy": float(th0[2])}
    if with_C:
        base["C"] = float(th0[3])
    significant = per[winner]["flag"]
    shape = (winner if shape_unique else
             "undetermined: " + ", ".join(f"k^{q}" for q in [winner] + [q for q in others
                                                                       if part[q]["padj"] >= alpha]))
    return {"test": test, "N": N, "base": base, "per": per, "winner": winner,
            "significant": significant, "shape_unique": shape_unique, "shape": shape,
            "partial": part, "any_flag": any(per[q]["flag"] for q in templates)}


def print_test(r, label):
    print(f"\n  [{label}] test {r['test']}: cells {r['N']}, "
          + ", ".join(f"{k} = {v:.4g}" for k, v in r["base"].items()))
    for q, v in r["per"].items():
        mark = "  ← significant" if v["flag"] else ""
        print(f"     template k^{q}: a = {v['a']:+.2e}   p(Holm) = {v['padj']:.2e}{mark}")
    print(f"     best: k^{r['winner']}, significant: {r['significant']}, shape: {r['shape']}")


# ----------------------------------------------------------------------------
# Synthetic field
# ----------------------------------------------------------------------------
def synth_segments(rng, n_seg, nt, n0, n1, dt, dx, hs=0.4, fp=0.35, theta0=None,
                   spread_s=4, U=(0.0, 0.0), C=0.0, q_inj=None, a_inj=0.0,
                   noise_m=0.01, k_ref=20.0):
    """
    Generator of segments η(t, x, y) of size (nt, n0, n1).
    Modes lie on a grid with twice the k density (field 2× larger, centre
    cropped), so the wavenumbers do not coincide with the analysis grid —
    as in a real sea.
    Dispersion: ω = √(g k (1 + C k²)) + k·U + a_inj √(g k)(k/k_ref)^q_inj.
    """
    N0, N1 = 2 * n0, 2 * n1
    ka = 2 * np.pi * np.fft.fftfreq(N0, dx)
    kb = 2 * np.pi * np.fft.fftfreq(N1, dx)
    KA, KB = np.meshgrid(ka, kb, indexing="ij")
    k = np.hypot(KA, KB)
    k[0, 0] = 1e-9
    th = np.arctan2(KB, KA)
    if theta0 is None:
        theta0 = rng.uniform(-np.pi, np.pi)
    kp = k_deep(fp)
    Fk = k ** -4 * np.exp(-1.25 * (kp / k) ** 2)
    Fk *= np.cos((th - theta0) / 2) ** (2 * spread_s)
    Fk[0, 0] = 0.0
    dk2 = (ka[1] - ka[0]) * (kb[1] - kb[0])
    amp = np.sqrt(2 * Fk * dk2)
    amp *= (hs / 4) / np.sqrt(np.sum(amp ** 2) / 2)
    om = np.sqrt(G * k * (1 + C * k ** 2)) + KA * U[0] + KB * U[1]
    if q_inj is not None:
        om += a_inj * np.sqrt(G * k) * (k / k_ref) ** q_inj
    phase0 = rng.uniform(0, 2 * np.pi, k.shape)
    keep = amp > amp.max() * 1e-6
    A = np.where(keep, amp, 0).astype(np.complex64) * (N0 * N1)
    t_glob = 0.0
    for _ in range(n_seg):
        z = np.empty((nt, n0, n1), dtype=np.float32)
        Gt = (A * np.exp(1j * (phase0 - om * t_glob))).astype(np.complex64)
        R = np.exp(-1j * om * dt).astype(np.complex64)
        for it in range(nt):
            full = sfft.ifft2(Gt, workers=-1)
            z[it] = full.real[n0 // 2:n0 // 2 + n0, n1 // 2:n1 // 2 + n1]
            Gt *= R
        z += rng.normal(0, noise_m, z.shape).astype(np.float32)
        t_glob += nt * dt
        yield z


def synth_run(rng, args, **kw):
    nt = CONFIG["spectrum"]["seg_frames"]
    acc = SpectrumAccumulator(nt, args.n0, args.n1, 1 / 12.0, 0.05, 0.05)
    for z in synth_segments(rng, args.segments, nt, args.n0, args.n1, 1 / 12.0, 0.05,
                            noise_m=args.noise, **kw):
        acc.add(z)
    return extract_ridge(acc.mean(), acc.f, acc.ka, acc.kb)


def barrier(flags):
    M, k = len(flags), int(sum(flags))
    pb = float(stats.binom.sf(k - 1, M, CONFIG["stats"]["alpha"])) if k > 0 else 1.0
    return k / max(M, 1), pb, pb >= CONFIG["stats"]["barrier_test_alpha"]


def random_U(rng):
    return tuple(rng.normal(0, 0.1, 2))


def run_synth_null(args):
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    print(f"config hash: {config_hash()}  (seed {args.seed}, tests {args.tests})")
    flags = {t: [] for t in args.tests}
    rows = []
    for i in range(args.m):
        t0 = time.time()
        U = random_U(rng)
        row, line = {}, []
        if "B" in args.tests:
            dB = synth_run(rng, args, U=U, C=0.0)          # no capillarity → null for B
            rB = run_test(dB, "B")
            row["B"] = rB
            flags["B"].append(rB["any_flag"])
            fl = [q for q, v in rB["per"].items() if v["flag"]]
            line.append(f"B flag {rB['any_flag']!s:5}" + (f" {['k^%d' % q for q in fl]}" if fl else ""))
        rows.append(row)
        print(f"  {i + 1:3d}/{args.m}: " + "  ".join(line) + f"  ({time.time() - t0:.0f} s)", flush=True)
        if (i + 1) % 10 == 0:                              # intermediate save
            with open(os.path.join(args.out, "synth_null.json"), "w", encoding="utf-8") as fh:
                json.dump({"config_hash": config_hash(), "args": vars(args), "runs": rows}, fh,
                          indent=1, default=float)
    for name, fl in flags.items():
        fpr, pb, ok = barrier(fl)
        print(f"test {name}: {sum(fl)}/{len(fl)}, FPR {fpr:.3f} (target ≤ {CONFIG['stats']['alpha']}), "
              f"binomial p = {pb:.3f} → {'PASSED' if ok else 'NOT PASSED → STOP'}")
    with open(os.path.join(args.out, "synth_null.json"), "w", encoding="utf-8") as fh:
        json.dump({"config_hash": config_hash(), "args": vars(args), "runs": rows}, fh,
                  indent=1, default=float)


def _declared(r):
    """What the scanner declares: 'none', 'k^q distinguishable' or 'k^q undetermined'."""
    if not r["significant"]:
        return "none"
    return f"k^{r['winner']}" + (" distinguishable" if r["shape_unique"] else " undetermined")


def _summary_line(label, decl, n):
    from collections import Counter
    c = Counter(decl)
    parts = ", ".join(f"{k}: {v / n:.2f}" for k, v in sorted(c.items()))
    print(f"  {label}: {parts}", flush=True)
    return dict(c)


def run_synth_inject(args):
    """
    Three separate metrics for every injection:
      detection_rate            — significant deviation (any shape)
      form_identification_rate  — the true shape wins AND is distinguishable
      false_form_rate           — a DISTINGUISHABLE but wrong shape is declared
    """
    os.makedirs(args.out, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    print(f"config hash: {config_hash()}  (seed {args.seed})")
    a_cap = C_NOMINAL * CONFIG["k_ref"] ** 2 / 2
    print(f"expected capillary effect at k_ref = {CONFIG['k_ref']} rad/m: δω/ω ≈ {a_cap:.2e}\n")
    cases = []
    if "cap" in args.parts:
        cases += [(f"capillarity ×{m:g}", 2, dict(C=C_NOMINAL * m), m) for m in args.mults]
    if "wrong" in args.parts:
        cases += [(f"wrong shape k^{q}", q, dict(C=0.0, q_inj=q, a_inj=a_cap), 1.0) for q in (1, 3)]
    table = []
    print("test B; distribution of declared outcomes (fraction of realizations)")
    for label, q_true, kw, mult in cases:
        decl, ratios = [], []
        for _ in range(args.m):
            r = run_test(synth_run(rng, args, U=random_U(rng), **kw), "B")
            decl.append(_declared(r))
            if mult > 0 and q_true == 2:
                ratios.append(r["per"][2]["a"] / (a_cap * mult))
        n = len(decl)
        counts = _summary_line(label, decl, n)
        det = sum(d != "none" for d in decl) / n
        ident = counts.get(f"k^{q_true} distinguishable", 0) / n
        false = sum(v for k, v in counts.items()
                    if k.endswith("distinguishable") and not k.startswith(f"k^{q_true} ")) / n
        if q_true == 2 and mult == 0:
            ident = false = float("nan")
        print(f"     detection {det:.2f} | form_identification {ident:.2f} | false_form {false:.2f}"
              + (f" | median a/a_true {np.median(ratios):.3f}" if ratios else ""), flush=True)
        table.append({"case": label, "q_true": q_true, "counts": counts, "detection_rate": det,
                      "form_identification_rate": ident, "false_form_rate": false,
                      "a_ratio_median": float(np.median(ratios)) if ratios else None})
        with open(os.path.join(args.out, "synth_inject.json"), "w", encoding="utf-8") as fh:
            json.dump({"config_hash": config_hash(), "args": vars(args), "table": table}, fh,
                      indent=1, ensure_ascii=False)


# ----------------------------------------------------------------------------
# Real data
# ----------------------------------------------------------------------------
def spectra_real(path, out):
    import xarray as xr
    ds = xr.open_dataset(path)
    Z = ds["Z"]
    X, Y = ds["X"].values, ds["Y"].values
    d0 = float(np.nanmedian(np.abs(np.gradient(X, axis=0)) + np.abs(np.gradient(Y, axis=0))))
    d1 = float(np.nanmedian(np.abs(np.gradient(X, axis=1)) + np.abs(np.gradient(Y, axis=1))))
    t = ds["time"].values
    dt = float(np.median(np.diff(t).astype("timedelta64[ns]").astype(float)) / 1e9)
    nt = CONFIG["spectrum"]["seg_frames"]
    n_total = Z.shape[0]
    n_seg = n_total // nt
    half = n_seg // 2
    print(f"grid {Z.shape[1]}×{Z.shape[2]}, steps {d0:.3f}/{d1:.3f} m, dt {dt:.4f} s; "
          f"{n_seg} segments of {nt} frames (Δf = {1 / (nt * dt):.4f} Hz)")
    spectra = {}
    for name, rng_ in (("search", range(0, half)), ("confirm", range(half, 2 * half))):
        cache = os.path.join(out, f"spectrum_{name}.npz")
        if os.path.exists(cache):
            c = np.load(cache)
            spectra[name] = (c["P"], c["f"], c["ka"], c["kb"])
            print(f"  {name}: from cache {cache}")
            continue
        acc = SpectrumAccumulator(nt, Z.shape[1], Z.shape[2], dt, d0, d1)
        for s in rng_:
            t0 = time.time()
            z = Z.isel(time=slice(s * nt, (s + 1) * nt)).values
            if np.isnan(z).any():
                z = np.nan_to_num(z, nan=0.0)
            acc.add(z)
            print(f"  {name}: segment {s + 1}/{n_seg} ({time.time() - t0:.0f} s)", flush=True)
        P = acc.mean()
        np.savez_compressed(cache, P=P.astype(np.float32), f=acc.f, ka=acc.ka, kb=acc.kb)
        spectra[name] = (P, acc.f, acc.ka, acc.kb)
    return spectra


def run_real(args):
    os.makedirs(args.out, exist_ok=True)
    print(f"config hash: {config_hash()}")
    spectra = spectra_real(args.file, args.out)
    summary = {"config_hash": config_hash(), "file": os.path.basename(args.file)}
    res = {}
    for name in ("search", "confirm"):
        d = extract_ridge(*spectra[name])
        print(f"\n=== {name.upper()}: candidate cells {d['n_candidates']}, above threshold {len(d['k'])}, "
              f"Δf = {d['df_hz']:.4f} Hz ===")
        rB = run_test(d, "B")
        print_test(rB, name)
        res[name] = {"B": rB, "n_cells": len(d["k"])}
    sB, cB = res["search"]["B"], res["confirm"]["B"]

    def b_ok(r):
        return r["significant"] and r["winner"] == 2 and r["per"][2]["a"] > 0
    b_found = b_ok(sB)
    b_conf = b_found and b_ok(cB)
    b_shape = b_conf and sB["shape_unique"] and cB["shape_unique"]
    print("\n=== RESULT (v0.5.x, test B only; the null test is synthetic) ===")
    print("capillary term k^2, a > 0: " + ("found and confirmed" if b_conf
                                          else "found but not confirmed" if b_found else "not found"))
    if b_conf:
        print("   the k^2 shape is " + ("DISTINGUISHABLE from k^1 and k^3 in both halves" if b_shape
                                     else "NOT DISTINGUISHABLE from a neighbouring shape → reported as "
                                          "'deviation with undetermined shape'"))
    for name in ("search", "confirm"):
        a = res[name]["B"]["per"][2]["a"]
        C_est = 2 * a / CONFIG["k_ref"] ** 2
        print(f"[interpretation, {name}] ℓ_c² ≈ {C_est:.2e} m² → γ ≈ {C_est * RHO * G:.4f} N/m "
              f"(clean water {GAMMA_NOMINAL})")
    summary.update({"results": res, "B_found": b_found, "B_confirmed": b_conf,
                    "B_shape_unique": b_shape})
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1, default=float, ensure_ascii=False)


def main():
    ap = argparse.ArgumentParser(description="Module 0 (sea waves) of the residual-first scanner")
    sub = ap.add_subparsers(dest="mode", required=True)
    for name in ("synth-null", "synth-inject"):
        p = sub.add_parser(name)
        p.add_argument("--m", type=int, default=(40 if name == "synth-null" else 10), help="number of realizations")
        p.add_argument("--segments", type=int, default=14, help="segments per realization")
        p.add_argument("--n0", type=int, default=271)
        p.add_argument("--n1", type=int, default=281)
        p.add_argument("--noise", type=float, default=0.01, help="stereo noise [m]")
        p.add_argument("--seed", type=int, default=1)
        if name == "synth-null":
            p.add_argument("--tests", default="B", choices=["B"])
        p.add_argument("--out", default=f"results/{name.replace('-', '_')}")
        if name == "synth-inject":
            p.add_argument("--mults", type=float, nargs="+", default=[0.0, 0.25, 0.5, 1.0, 2.0])
            p.add_argument("--parts", nargs="+", default=["cap", "wrong"], choices=["cap", "wrong"])
    r = sub.add_parser("real")
    r.add_argument("--file", required=True)
    r.add_argument("--out", default="results/real")
    args = ap.parse_args()
    if args.mode == "synth-null":
        return run_synth_null(args)
    if args.mode == "synth-inject":
        return run_synth_inject(args)
    return run_real(args)


if __name__ == "__main__":
    sys.exit(main())
