#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piano_scanner.py — Module 0 (piano) of the residual-first anomaly scanner
=========================================================================

Idea:
  Theory subtracted:        ideal string, f_n = n * f0 (f0 is fitted).
  Residuals:                Δf_n = f_n - n * f0_fit   [Hz]
  Templates (fixed):        Δf_n ∝ n^p,  p ∈ {0, 2, 3, 4}
                            (p = 1 is fully absorbed by f0 and is not tested)
  Expectation for a piano:  stiff string → f_n = n f0 √(1 + B n²)
                            → Δf_n ≈ (f0 B / 2) n³  → p = 3 must win.
  Null test:                sound with exactly harmonic partials (synthetic,
                            or bowed violin/cello)
                            → the scanner must report "nothing".

Modes:
  synth-null    calibration barrier on synthetic harmonic tones
  synth-inject  sensitivity and shape identification (correct and wrong shapes)
  real          real recordings (e.g. University of Iowa MIS) with a
                search / confirmation split and an optional real null test

Examples:
  python piano_scanner.py synth-null   --out out_null
  python piano_scanner.py synth-inject --out out_inject
  python piano_scanner.py real --data Piano_mf --null-data Violin_nonvib --out out_real

Dependencies: numpy, scipy; for AIFF files also soundfile (pip install soundfile).
The physical interpretation (estimate of B) is kept separate and printed ONLY
as secondary information — the scanner is not a "B detector".
"""

import argparse
import csv
import hashlib
import json
import os
import re
import sys

import numpy as np
from scipy import stats

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# ----------------------------------------------------------------------------
# FIXED CONFIGURATION — not changed after pre-registration.
# Its hash is written to every output so that it can be cited in the registration.
# ----------------------------------------------------------------------------
CONFIG = {
    "version": "0.4.1",
    "analysis": {
        "skip_after_peak_s": 0.10,     # skip after tone onset (transient)
        "window_s": 1.0,               # length of the analysed segment
        "min_window_s": 0.3,           # shorter segment → note excluded
        "zero_pad": 8,                 # zero padding factor for the FFT
        "max_partials": 60,            # highest partial number searched
        "fit_max_n": 25,               # only n ≤ fit_max_n enter the fit
        "min_partials": 8,             # fewer usable partials → note excluded
        "seed_tol_rel": 0.06,          # search for the first partials: ±6 %
        "seed_max_n": 3,               # the starting partial may be n = 1..3
        "track_halfwidth_rel": 0.25,   # tracking window: ±0.25 partial spacing
        "snr_db": 20.0,                # peak threshold above the median floor
        "max_consecutive_misses": 3,   # stop after this many misses in a row
        "f_max_hz": 12000.0,           # upper frequency limit
        "weights": "snr",              # fit weights: σ_n ∝ 1/SNR
        "onset_rel": 0.1,              # tone onset: 10 % of the maximum
        "freq_method": "phase_demod",  # frequencies by phase demodulation (vibrato-robust)
        "ref_max_n": 6,                # common phase φ(t) from the 3 strongest partials with n ≤ 6
        "phase_halfband_rel": 0.4,     # band around each partial: ±0.4·f0
        "phase_max_rms": 1.0,          # partial with noisier phase [rad] is dropped
    },
    "templates": [0, 2, 3, 4],
    "stats": {
        "alpha_note": 0.01,            # per-note threshold (after Holm over templates)
        "min_effect_cents": 0.1,       # minimum effect: |deviation at n = 10| ≥ 0.1 cent
        "barrier_test_alpha": 0.05,    # barrier: STOP if the false-alarm count is
                                       # significantly above α (one-sided binomial test)
        "group_alpha": 0.001,          # threshold for structure in the search group
        "confirm_alpha": 0.01,         # threshold in the confirmation group
        "shape_majority": 0.6,         # the shape must win in ≥ 60 % of the notes
    },
    "groups": {"rule": "midi_parity", "search": "even", "confirm": "odd"},
    "primary_dynamic": "mf",
}


def config_hash(cfg=CONFIG):
    s = json.dumps(cfg, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(s).hexdigest()[:16]


# ----------------------------------------------------------------------------
# File names and notes
# ----------------------------------------------------------------------------
NOTE_INDEX = {"C": 0, "D": 2, "E": 4, "F": 5, "G": 7, "A": 9, "B": 11}
NOTE_RE = re.compile(r"(?<![A-Za-z])([A-G])(b|#)?(-?\d)(?!\d)")
DYNAMICS = ("pp", "mf", "ff")


def note_to_midi(letter, accidental, octave):
    idx = NOTE_INDEX[letter] + (-1 if accidental == "b" else 1 if accidental == "#" else 0)
    return 12 * (int(octave) + 1) + idx


def midi_to_freq(midi):
    return 440.0 * 2.0 ** ((midi - 69) / 12.0)


def parse_filename(path):
    """Returns (midi, note name, dynamic) or None."""
    base = os.path.basename(path)
    tokens = re.split(r"[.\s_\-]+", base)
    dyn = next((t for t in tokens if t in DYNAMICS), None)
    matches = NOTE_RE.findall(base)
    if not matches:
        return None
    letter, acc, octv = matches[-1]
    midi = note_to_midi(letter, acc, octv)
    return midi, f"{letter}{acc}{octv}", dyn


def load_audio(path):
    """Mono float signal and sampling rate."""
    try:
        import soundfile as sf
        x, sr = sf.read(path, always_2d=True)
        x = x.mean(axis=1)
    except ImportError:
        if path.lower().endswith(".wav"):
            from scipy.io import wavfile
            sr, x = wavfile.read(path)
            x = x.astype(float)
            if x.ndim == 2:
                x = x.mean(axis=1)
        else:
            raise RuntimeError("AIFF files need soundfile: pip install soundfile")
    x = np.asarray(x, dtype=float)
    peak = np.max(np.abs(x))
    if peak > 0:
        x = x / peak
    return x, int(sr)


# ----------------------------------------------------------------------------
# Spectrum and partial tracking
# ----------------------------------------------------------------------------
def analysis_segment(x, sr, a):
    """The segment starts skip_after_peak_s after the tone ONSET (first
    crossing of onset_rel of the envelope maximum). For a struck string (piano)
    the onset practically coincides with the maximum; for a bowed tone it does not."""
    k = max(1, int(0.01 * sr))
    env = np.convolve(np.abs(x), np.ones(k) / k, mode="same")
    i_on = int(np.argmax(env >= a["onset_rel"] * env.max()))
    start = i_on + int(a["skip_after_peak_s"] * sr)
    seg = x[start:start + int(a["window_s"] * sr)]
    if len(seg) < int(a["min_window_s"] * sr):
        return None
    return seg


def spectrum(seg, sr, pad):
    w = np.hanning(len(seg))
    nfft = 1 << int(np.ceil(np.log2(len(seg) * pad)))
    mag = np.abs(np.fft.rfft(seg * w, nfft))
    return mag, sr / nfft


def find_peak(mag, df, lo, hi, floor_lo, floor_hi, snr_db):
    i0, i1 = int(np.ceil(lo / df)), int(np.floor(hi / df))
    if i0 < 1 or i1 >= len(mag) - 1 or i1 - i0 < 3:
        return None
    k = i0 + int(np.argmax(mag[i0:i1 + 1]))
    if k == i0 or k == i1:                         # at the edge → not a true peak
        return None
    j0, j1 = max(1, int(floor_lo / df)), min(len(mag) - 1, int(floor_hi / df))
    floor = np.median(mag[j0:j1 + 1]) + 1e-300
    snr = mag[k] / floor
    if 20 * np.log10(snr) < snr_db:
        return None
    a_, b_, c_ = np.log(mag[k - 1:k + 2] + 1e-300)
    denom = a_ - 2 * b_ + c_
    delta = 0.5 * (a_ - c_) / denom if denom < 0 else 0.0
    return (k + delta) * df, snr


def _band(X, f, f_c, halfband):
    """Isolates a band ±halfband around f_c (smooth mask) → complex signal."""
    u = (f - f_c) / halfband
    mask = np.where(np.abs(u) < 1, 0.5 * (1 + np.cos(np.pi * u)), 0.0)
    mask[f < 0] = 0.0                                   # analytic signal
    return np.fft.ifft(X * mask) * 2


def refine_phase_demod(seg, sr, found, a):
    """
    Partial frequencies by phase demodulation.
    1) The common phase φ(t) = mean of ψ_k(t)/k is measured from the strongest
       low partials (n ≤ ref_max_n); it contains the vibrato and pitch drift.
    2) For each partial the residual phase θ_n(t) = arg[y_n · e^{-i n φ(t)}]
       is fitted with a straight line (identical time weights for all n).
       f_n = n·f_ref + slope/2π.
    For an exactly harmonic sound θ_n is constant for every n, even with vibrato.
    A partial whose θ_n rms exceeds phase_max_rms (noise / lost phase) is dropped.
    """
    N = len(seg)
    X = np.fft.fft(seg)
    f = np.fft.fftfreq(N, 1.0 / sr)
    t = np.arange(N) / sr
    f1 = np.median([fn / n for n, fn, _ in found[:5]])
    hb = a["phase_halfband_rel"] * f1
    i0, i1 = int(0.1 * N), int(0.9 * N)
    tt = t[i0:i1] - t[i0:i1].mean()
    A = np.column_stack([tt, np.ones_like(tt)])

    low = [p for p in found if p[0] <= a["ref_max_n"]]
    ref = sorted(low, key=lambda p: -p[2])[:3]
    phis = []
    for n, fn, _ in ref:
        psi = np.unwrap(np.angle(_band(X, f, fn, hb)))
        phis.append((psi - psi[i0]) / n)
    phi = np.mean(phis, axis=0)
    c_ref, *_ = np.linalg.lstsq(A, phi[i0:i1], rcond=None)
    f_ref = c_ref[0] / (2 * np.pi)

    out = []
    for n, fn, snr in found:
        z = _band(X, f, fn, hb) * np.exp(-1j * n * phi)
        th = np.unwrap(np.angle(z))[i0:i1]
        c, *_ = np.linalg.lstsq(A, th, rcond=None)
        rms = float(np.sqrt(np.mean((th - A @ c) ** 2)))
        if rms <= a["phase_max_rms"]:
            out.append((n, n * f_ref + c[0] / (2 * np.pi), snr))
    return out


def track_partials(seg, sr, f_nom, a):
    """
    Tracks partials without assuming a model of the deviation: the next partial
    is searched around a linear extrapolation from the last two found.
    """
    mag, df = spectrum(seg, sr, a["zero_pad"])
    f_max = min(a["f_max_hz"], 0.45 * sr)
    found = []                      # (n, f, snr)
    misses = 0
    for n in range(1, a["max_partials"] + 1):
        if len(found) >= 2:
            (n1, f1, _), (n2, f2, _) = found[-2], found[-1]
            s = (f2 - f1) / (n2 - n1)
            f_pred, half = f2 + s * (n - n2), a["track_halfwidth_rel"] * s
        elif len(found) == 1:
            n1, f1, _ = found[0]
            s = f1 / n1
            f_pred, half = s * n, a["track_halfwidth_rel"] * s
        else:
            if n > a["seed_max_n"]:
                break
            s = f_nom
            f_pred, half = n * f_nom, a["seed_tol_rel"] * n * f_nom
        if f_pred + half > f_max:
            break
        pk = find_peak(mag, df, f_pred - half, f_pred + half,
                       f_pred - 0.5 * s, f_pred + 0.5 * s, a["snr_db"])
        if pk is None:
            misses += 1
            if found and misses >= a["max_consecutive_misses"]:
                break
            continue
        misses = 0
        found.append((n, pk[0], pk[1]))
    if a["freq_method"] == "phase_demod" and len(found) >= 2:
        found = refine_phase_demod(seg, sr, found, a)
    return found


# ----------------------------------------------------------------------------
# Residuals, templates and statistics
# ----------------------------------------------------------------------------
def holm(pvals):
    p = np.asarray(pvals, dtype=float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


def fit_templates(n, f, templates, w=None):
    """
    Baseline f = f0·n; for each template: f = f0·n + a·n^p and an F-test.
    w = weights 1/σ_n. The peak-frequency error is ∝ 1/SNR (verified on
    synthetic data), therefore w = peak SNR.
    """
    n = np.asarray(n, float)
    f = np.asarray(f, float)
    N = len(n)
    w = np.ones(N) if w is None else np.asarray(w, float)
    w = w / np.median(w)
    nmax = n.max()
    x = n / nmax                                       # scaling for numerical stability
    fw = f * w
    c0, *_ = np.linalg.lstsq((x * w)[:, None], fw, rcond=None)
    rss0 = float(np.sum((fw - x * w * c0[0]) ** 2))
    f0_fit = c0[0] / nmax
    out = {"f0_fit": f0_fit, "rss0": rss0, "resid": f - n * f0_fit, "per": {}}
    for p in templates:
        X = np.column_stack([x, x ** p]) * w[:, None]
        c, *_ = np.linalg.lstsq(X, fw, rcond=None)
        rss1 = float(np.sum((fw - X @ c) ** 2))
        dof = N - 2
        if rss1 <= 1e-18 * max(rss0, 1e-300):
            F, pval = np.inf, 0.0
        else:
            F = max(0.0, (rss0 - rss1) / (rss1 / dof))
            pval = float(stats.f.sf(F, 1, dof))
        out["per"][p] = {"F": F, "p": pval, "rss": rss1,
                         "a": c[1] / nmax ** p, "f0": c[0] / nmax}
    padj = holm([out["per"][p]["p"] for p in templates])
    for p, pa in zip(templates, padj):
        out["per"][p]["padj"] = float(pa)
    out["winner"] = min(templates, key=lambda p: out["per"][p]["rss"])
    return out


def analyze_signal(x, sr, f_nom, cfg=CONFIG):
    a = cfg["analysis"]
    seg = analysis_segment(x, sr, a)
    if seg is None:
        return None, "signal too short"
    partials = track_partials(seg, sr, f_nom, a)
    n_all = len(partials)
    partials = [p for p in partials if p[0] <= a["fit_max_n"]]
    if len(partials) < a["min_partials"]:
        return None, f"too few partials ({len(partials)})"
    n = [p[0] for p in partials]
    f = [p[1] for p in partials]
    snr = [p[2] for p in partials]
    fit = fit_templates(n, f, cfg["templates"], w=snr if cfg["analysis"]["weights"] == "snr" else None)
    alpha = cfg["stats"]["alpha_note"]
    w = fit["winner"]
    for p, v in fit["per"].items():
        # effect size: deviation of the 10th partial relative to 10·f0, in cents
        rel = v["a"] * 10.0 ** p / (10.0 * v["f0"])
        v["cents10"] = float(1200 * np.log2(max(1e-12, 1 + rel)))
    min_c = cfg["stats"]["min_effect_cents"]
    flag = {p: bool(fit["per"][p]["padj"] < alpha and abs(fit["per"][p]["cents10"]) >= min_c)
            for p in cfg["templates"]}
    res = {
        "n_partials": n_all, "n_fit": len(n), "f0_fit": fit["f0_fit"],
        "winner": w,
        "significant": flag[w],
        "flags": flag,
        "per": fit["per"], "n": n, "f": f, "resid": fit["resid"],
        # secondary interpretation (separate module): stiff string, Δf ≈ (f0 B/2) n³
        "B_est": 2.0 * fit["per"][3]["a"] / fit["per"][3]["f0"] if 3 in fit["per"] else np.nan,
    }
    return res, "ok"


# ----------------------------------------------------------------------------
# Synthetic tones
# ----------------------------------------------------------------------------
def synth_note(f0, rng, sr=44100, dur=1.3, model="harmonic", B=0.0,
               p=3, cents_at_10=0.0, noise_db=-60.0, jitter_cents=0.0):
    """
    model = 'harmonic' : f_n = n f0
            'stiff'    : f_n = n f0 √(1 + B n²)         (exact stiff string)
            'power'    : f_n = n f0 + a n^p, where the deviation at n = 10
                         equals cents_at_10 cents
    """
    nmax = min(80, int(0.45 * sr / f0))
    n = np.arange(1, nmax + 1, dtype=float)
    if model == "stiff":
        fn = n * f0 * np.sqrt(1.0 + B * n ** 2)
    elif model == "power":
        d10 = 10 * f0 * (2 ** (cents_at_10 / 1200.0) - 1)
        fn = n * f0 + d10 * (n / 10.0) ** p
    else:
        fn = n * f0
    if jitter_cents > 0:
        fn = fn * 2 ** (rng.normal(0, jitter_cents, len(fn)) / 1200.0)
    keep = (fn > 0) & (fn < 0.45 * sr)
    n, fn = n[keep], fn[keep]
    amps = n ** -0.8 * rng.lognormal(0, 0.3, len(n))
    tau0 = float(np.clip(3.0 * np.sqrt(261.6 / f0), 0.3, 8.0))
    t = np.arange(int(dur * sr)) / sr - 0.05
    on = t >= 0
    tt = t[on]
    x = np.zeros_like(t)
    for k in range(len(n)):
        tau = tau0 / (1 + 0.08 * n[k])
        env = np.exp(-tt / tau) * (1 - np.exp(-tt / 0.003))
        x[on] += amps[k] * env * np.sin(2 * np.pi * fn[k] * tt + rng.uniform(0, 2 * np.pi))
    x /= np.max(np.abs(x))
    x += rng.normal(0, 10 ** (noise_db / 20.0), len(x))
    return x, sr


def random_f0s(rng, m, lo=27.5, hi=1100.0):
    return np.exp(rng.uniform(np.log(lo), np.log(hi), m))


def detune(rng, f0, cents=10.0):
    """The nominal frequency differs from the true one (as for a real piano)."""
    return f0 * 2 ** (rng.normal(0, cents) / 1200.0)


# ----------------------------------------------------------------------------
# Group analysis
# ----------------------------------------------------------------------------
def barrier(results, cfg=CONFIG):
    """Calibration barrier. Returns (FPR, binomial p, passed)."""
    M = len(results)
    k = int(sum(r["significant"] for r in results))
    if M == 0:
        return np.nan, np.nan, False
    pb = float(stats.binom.sf(k - 1, M, cfg["stats"]["alpha_note"])) if k > 0 else 1.0
    return k / M, pb, pb >= cfg["stats"]["barrier_test_alpha"]


def group_summary(results, cfg=CONFIG):
    st = cfg["stats"]
    M = len(results)
    out = {"M": M, "per_template": {}, "shape": None, "shape_share": 0.0}
    if M == 0:
        return out
    for p in cfg["templates"]:
        k = sum(r["flags"][p] for r in results)
        pb = float(stats.binom.sf(k - 1, M, st["alpha_note"])) if k > 0 else 1.0
        out["per_template"][p] = {"k": k, "binom_p": pb}
    min_p = min(v["binom_p"] for v in out["per_template"].values())
    out["group_p"] = min(1.0, min_p * len(cfg["templates"]))
    sig = [r for r in results if r["significant"]]
    out["n_significant"] = len(sig)
    if sig:
        winners = [r["winner"] for r in sig]
        vals, counts = np.unique(winners, return_counts=True)
        best = int(vals[np.argmax(counts)])
        share = counts.max() / len(sig)
        out["shape_share"] = float(share)
        out["shape"] = best if share >= st["shape_majority"] else None
    return out


def print_group(title, g, cfg=CONFIG):
    print(f"\n=== {title} ===")
    print(f"notes analysed: {g['M']}")
    if g["M"] == 0:
        return
    for p, v in g["per_template"].items():
        print(f"  template n^{p}: significant notes {v['k']:3d}/{g['M']}   binomial p = {v['binom_p']:.2e}")
    print(f"  significant notes in total: {g['n_significant']}")
    print(f"  group p (Bonferroni over templates): {g['group_p']:.2e}")
    shape = "undetermined" if g["shape"] is None else f"n^{g['shape']}"
    print(f"  shape: {shape} (share among significant notes: {g['shape_share']:.0%})")


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------
def ensure_out(out):
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "config.json"), "w", encoding="utf-8") as fh:
        json.dump({"config": CONFIG, "config_hash": config_hash()}, fh,
                  indent=2, ensure_ascii=False)


def write_notes_csv(path, rows, cfg=CONFIG):
    cols = ["label", "midi", "note", "dynamic", "group", "f_nom", "n_partials",
            "n_fit", "f0_fit", "winner", "significant", "B_est"]
    for p in cfg["templates"]:
        cols += [f"F_{p}", f"p_{p}", f"padj_{p}", f"a_{p}", f"cents10_{p}"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for meta, r in rows:
            line = [meta.get(c, "") for c in cols[:6]]
            line += [r["n_partials"], r["n_fit"], f"{r['f0_fit']:.6f}", r["winner"],
                     int(r["significant"]), f"{r['B_est']:.3e}"]
            for p in cfg["templates"]:
                v = r["per"][p]
                line += [f"{v['F']:.4g}", f"{v['p']:.3e}", f"{v['padj']:.3e}", f"{v['a']:.4e}",
                         f"{v['cents10']:.4f}"]
            w.writerow(line)


def save_plots(out, rows, k):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib not available — plots skipped)")
        return
    for meta, r in rows[:k]:
        n = np.array(r["n"], float)
        w = r["winner"]
        v = r["per"][w]
        model = v["f0"] * n + v["a"] * n ** w - r["f0_fit"] * n
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(n, r["resid"], "o", label="residuals Δf_n")
        ax.plot(n, model, "-", label=f"best template n^{w}")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel("partial number n")
        ax.set_ylabel("Δf_n = f_n − n·f0 [Hz]")
        ax.set_title(f"{meta.get('label', '')}  (significant: {r['significant']})")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(out, f"resid_{meta.get('label', 'note')}.png"), dpi=120)
        plt.close(fig)


# ----------------------------------------------------------------------------
# Modes
# ----------------------------------------------------------------------------
def run_synth_null(args):
    ensure_out(args.out)
    rng = np.random.default_rng(args.seed)
    rows, skipped = [], 0
    for i, f0 in enumerate(random_f0s(rng, args.m)):
        x, sr = synth_note(f0, rng, model="harmonic", noise_db=args.noise_db,
                           jitter_cents=args.jitter)
        r, why = analyze_signal(x, sr, detune(rng, f0))
        if r is None:
            skipped += 1
            continue
        rows.append(({"label": f"null{i:04d}", "f_nom": f"{f0:.3f}"}, r))
    results = [r for _, r in rows]
    fpr, pb, passed = barrier(results)
    print(f"config hash: {config_hash()}")
    print(f"synthetic null test: {len(results)} notes (excluded {skipped})")
    for p in CONFIG["templates"]:
        print(f"  FPR template n^{p}: {np.mean([r['flags'][p] for r in results]):.3f}")
    print(f"per-note FPR (Holm): {fpr:.3f}  (expected ≤ {CONFIG['stats']['alpha_note']}), "
          f"binomial p = {pb:.3f}")
    print("BARRIER: PASSED" if passed else "BARRIER: NOT PASSED → STOP")
    write_notes_csv(os.path.join(args.out, "notes_synth_null.csv"), rows)
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"config_hash": config_hash(), "fpr": fpr, "passed": passed,
                   "M": len(results), "skipped": skipped}, fh, indent=2)
    return 0 if passed else 2


def run_synth_inject(args):
    ensure_out(args.out)
    rng = np.random.default_rng(args.seed)
    templates = CONFIG["templates"]
    table = []
    print(f"config hash: {config_hash()}")
    print("\nInjections with shape n^p (amplitude = deviation at n = 10, in cents)")
    print("true shape | cents | detected | correct shape among detected")
    for p_true in templates:
        for c in args.cents:
            det = correct = total = 0
            for f0 in random_f0s(rng, args.m):
                x, sr = synth_note(f0, rng, model="power", p=p_true, cents_at_10=c,
                                   noise_db=args.noise_db, jitter_cents=args.jitter)
                r, _ = analyze_signal(x, sr, detune(rng, f0))
                if r is None:
                    continue
                total += 1
                if r["significant"]:
                    det += 1
                    correct += int(r["winner"] == p_true)
            dr = det / total if total else np.nan
            cr = correct / det if det else np.nan
            table.append({"model": "power", "p": p_true, "cents": c, "N": total,
                          "detect": dr, "correct_shape": cr})
            print(f"   n^{p_true}      | {c:5.1f} | {dr:6.2f}  | {cr:6.2f}")
    print("\nExact stiff string f_n = n f0 √(1+B n²) — n^3 must win")
    print("      B     | detected | correct shape (n^3) | median B_est / B")
    for B in args.B:
        det = correct = total = 0
        ratios = []
        for f0 in random_f0s(rng, args.m):
            x, sr = synth_note(f0, rng, model="stiff", B=B,
                               noise_db=args.noise_db, jitter_cents=args.jitter)
            r, _ = analyze_signal(x, sr, detune(rng, f0))
            if r is None:
                continue
            total += 1
            if r["significant"]:
                det += 1
                correct += int(r["winner"] == 3)
                ratios.append(r["B_est"] / B)
        dr = det / total if total else np.nan
        cr = correct / det if det else np.nan
        med = float(np.median(ratios)) if ratios else np.nan
        table.append({"model": "stiff", "B": B, "N": total, "detect": dr,
                      "correct_shape": cr, "B_ratio_median": med})
        print(f"  {B:9.1e} | {dr:6.2f}  | {cr:6.2f}               | {med:6.3f}")
    with open(os.path.join(args.out, "injections.json"), "w", encoding="utf-8") as fh:
        json.dump({"config_hash": config_hash(), "table": table}, fh, indent=2)
    return 0


def collect_files(folder, dynamic):
    rows = []
    for root, _, files in os.walk(folder):
        for fn in sorted(files):
            if not fn.lower().endswith((".aif", ".aiff", ".wav", ".flac")):
                continue
            info = parse_filename(fn)
            if info is None:
                continue
            midi, name, dyn = info
            if dynamic != "all" and dyn is not None and dyn != dynamic:
                continue
            rows.append((os.path.join(root, fn), midi, name, dyn))
    return rows


def analyze_files(files, label_prefix=""):
    rows, skipped = [], []
    for path, midi, name, dyn in files:
        x, sr = load_audio(path)
        r, why = analyze_signal(x, sr, midi_to_freq(midi))
        meta = {"label": f"{label_prefix}{name}_{dyn}", "midi": midi, "note": name,
                "dynamic": dyn, "f_nom": f"{midi_to_freq(midi):.3f}"}
        if r is None:
            skipped.append((os.path.basename(path), why))
            continue
        rows.append((meta, r))
    return rows, skipped


def run_real(args):
    ensure_out(args.out)
    print(f"config hash: {config_hash()}")
    summary = {"config_hash": config_hash()}

    # 1) Real null test (calibration barrier)
    if args.null_data:
        nfiles = collect_files(args.null_data, "all")
        nrows, nskip = analyze_files(nfiles, "null_")
        nres = [r for _, r in nrows]
        fpr, pb, passed = barrier(nres)
        from collections import Counter
        reasons = Counter(why.split(" (")[0] for _, why in nskip)
        if reasons:
            print("  exclusion reasons:", dict(reasons))
        print(f"\nreal null test: {len(nres)} notes (excluded {len(nskip)})")
        print(f"per-note FPR: {fpr:.3f}  (expected ≤ {CONFIG['stats']['alpha_note']}), "
              f"binomial p = {pb:.3f}")
        write_notes_csv(os.path.join(args.out, "notes_null.csv"), nrows)
        summary["null"] = {"M": len(nres), "fpr": fpr, "binom_p": pb}
        if not passed:
            print("BARRIER: NOT PASSED → STOP. The piano analysis is not run.")
            summary["stopped"] = True
            with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as fh:
                json.dump(summary, fh, indent=2, default=float)
            return 2
        print("BARRIER: PASSED")
        if args.null_only:
            print("(--null-only: the piano is not analysed)")
            with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as fh:
                json.dump(summary, fh, indent=2, default=float)
            return 0
    else:
        print("\n(no real null test — the synthetic barrier must have passed)")

    if args.null_only:
        print("--null-only requires --null-data")
        return 1
    if not args.data:
        print("--data is missing")
        return 1

    # 2) Piano: split by a fixed rule
    dyn = args.dynamic or CONFIG["primary_dynamic"]
    files = collect_files(args.data, dyn)
    rows, skipped = analyze_files(files)
    for meta, _ in rows:
        meta["group"] = "search" if meta["midi"] % 2 == 0 else "confirm"
    for name, why in skipped:
        print(f"  excluded: {name} ({why})")
    search = [r for m, r in rows if m["group"] == "search"]
    confirm = [r for m, r in rows if m["group"] == "confirm"]
    gs, gc = group_summary(search), group_summary(confirm)
    print_group("SEARCH GROUP (even MIDI)", gs)
    print_group("CONFIRMATION GROUP (odd MIDI)", gc)

    st = CONFIG["stats"]
    found = gs["M"] > 0 and gs["group_p"] < st["group_alpha"] and gs["shape"] is not None
    confirmed = False
    if found:
        s = gs["shape"]
        confirmed = (gc["M"] > 0 and gc["shape"] == s
                     and gc["per_template"][s]["binom_p"] < st["confirm_alpha"])
    print("\n=== RESULT ===")
    if not found:
        print("no significant structure with a clear shape in the search group")
    else:
        print(f"search: structure with shape n^{gs['shape']}")
        print("confirmation: " + ("YES — the same shape is significant" if confirmed else "NO"))

    # 3) Secondary interpretation (separate module): inharmonicity coefficient
    B = [r["B_est"] for r in search + confirm if r["significant"] and r["winner"] == 3]
    if B:
        print(f"\n[interpretation] stiff-string B: median {np.median(B):.2e}, "
              f"range {np.min(B):.1e} … {np.max(B):.1e} ({len(B)} notes)")

    write_notes_csv(os.path.join(args.out, "notes_piano.csv"), rows)
    if args.plots:
        save_plots(args.out, rows, args.plots)
    summary.update({"dynamic": dyn, "search": gs, "confirm": gc,
                    "found": found, "confirmed": confirmed,
                    "skipped": skipped})
    with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, default=float, ensure_ascii=False)
    return 0


def main():
    ap = argparse.ArgumentParser(description="Module 0 (piano) of the residual-first scanner")
    sub = ap.add_subparsers(dest="mode", required=True)

    a = sub.add_parser("synth-null", help="calibration barrier on synthetic tones")
    a.add_argument("--m", type=int, default=300)
    a.add_argument("--noise-db", type=float, default=-60.0)
    a.add_argument("--jitter", type=float, default=0.0, help="random detuning [cents]")
    a.add_argument("--seed", type=int, default=1)
    a.add_argument("--out", default="out_synth_null")

    b = sub.add_parser("synth-inject", help="sensitivity and shape identification")
    b.add_argument("--m", type=int, default=60)
    b.add_argument("--cents", type=float, nargs="+", default=[0.5, 1, 2, 5, 20])
    b.add_argument("--B", type=float, nargs="+", default=[1e-5, 1e-4, 4e-4, 1e-3])
    b.add_argument("--noise-db", type=float, default=-60.0)
    b.add_argument("--jitter", type=float, default=0.0)
    b.add_argument("--seed", type=int, default=2)
    b.add_argument("--out", default="out_synth_inject")

    c = sub.add_parser("real", help="real recordings")
    c.add_argument("--data", help="folder with the piano notes")
    c.add_argument("--null-only", action="store_true",
                   help="null test only; the piano is not analysed")
    c.add_argument("--null-data", help="folder with a harmonic instrument (null test)")
    c.add_argument("--dynamic", choices=list(DYNAMICS) + ["all"])
    c.add_argument("--plots", type=int, default=0, help="number of residual plots")
    c.add_argument("--out", default="out_real")

    args = ap.parse_args()
    if args.mode == "synth-null":
        return run_synth_null(args)
    if args.mode == "synth-inject":
        return run_synth_inject(args)
    return run_real(args)


if __name__ == "__main__":
    sys.exit(main())
