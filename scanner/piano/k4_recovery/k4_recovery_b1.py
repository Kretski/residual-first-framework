#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_recovery_b1.py — B1 wrapper for the blind k⁴-equivalent recovery test.

Specification: RFF_K4_RECOVERY_TEST_SPEC_B0_v0.4_CORRECTED.md
               SHA-256 662b2303aff00bb1000d8a0ed93be64fc5f7fdc72aae811567fb2f5810a4d042

This file is an EXTERNAL layer over scanner/piano/piano_scanner.py (e7c7520).
The scanner is imported and used unchanged; its file hash and config hash are
verified at start-up (B0 §2). Additions (B0 §3):
  W1  sign record                    -> analyse_case(): sign_winner, sign_B
  W2  discrimination statistic       -> discrimination()
  W3  coefficient B̂ and interval     -> b_hat_interval()
  W4  declared outcome               -> analyse_case(): declared
  W5  sealed generation              -> make_signal(), synth_mixture(), run()
  W6  manifest, seeds, counts        -> build_manifest(), assert_manifest()
  parity guard                       -> analyse_case()

The analysis path receives only (x, sr, detuned f_nom). The truth is written to
truth.jsonl and is read only by k4_score_b1.py.

Modes:
  b2     the single pre-registered run: seed 20261122, the §4 manifest (5900 cases)
  smoke  technical check only (B0 §6.4): seed 1, 3 cases per cell (108 cases),
         plus the W3 covariance cross-check against scipy.optimize.curve_fit

Usage (only with the corresponding approval):
  python k4_recovery_b1.py smoke --out <dir>
  python k4_recovery_b1.py b2    --out <dir>
"""

import argparse
import hashlib
import json
import os
import sys
import time

import numpy as np
from scipy import optimize, stats

HERE = os.path.dirname(os.path.abspath(__file__))
SCANNER_DIR = os.path.dirname(HERE)                       # scanner/piano
sys.path.insert(0, SCANNER_DIR)
import piano_scanner as ps                                # noqa: E402  (unchanged scanner)

# ----------------------------------------------------------------------------
# Frozen inputs (B0 §2)
# ----------------------------------------------------------------------------
SPEC_NAME = "RFF_K4_RECOVERY_TEST_SPEC_B0_v0.4_CORRECTED.md"
SPEC_SHA256 = "662b2303aff00bb1000d8a0ed93be64fc5f7fdc72aae811567fb2f5810a4d042"
SCANNER_SHA256 = "335ba8b8abaef4b424b36a516f55df568b9345f83cc8b6bc1c0dc1cc33c9e5d4"
SCANNER_CONFIG_HASH = "5bdf34b190a4b9ef"

SEED_B2 = 20261122
SEED_SMOKE = 1
N_SPAWN_B2 = 5901                                         # index 5900 = shuffle stream
M_SMOKE = 3

PARITY_RTOL = 1e-12
CROSSCHECK_RTOL = 1e-6

# Generator settings (B0 §4): synth_note defaults, except the R stratum
NOISE_DEFAULT = -60.0
JITTER_DEFAULT = 0.0
NOISE_R = -40.0
JITTER_R = 0.3


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_frozen_inputs():
    """Abort if the scanner or its configuration differ from e7c7520 / v0.4.1."""
    scanner_path = os.path.join(SCANNER_DIR, "piano_scanner.py")
    got = sha256_file(scanner_path)
    if got != SCANNER_SHA256:
        sys.exit(f"ABORT: piano_scanner.py sha256 {got} != {SCANNER_SHA256}")
    ch = ps.config_hash()
    if ch != SCANNER_CONFIG_HASH:
        sys.exit(f"ABORT: config_hash {ch} != {SCANNER_CONFIG_HASH}")
    if ps.CONFIG["templates"] != [0, 2, 3, 4] or ps.CONFIG["stats"]["alpha_note"] != 0.01:
        sys.exit("ABORT: unexpected templates or alpha in the scanner CONFIG")
    return {"piano_scanner_sha256": got, "config_hash": ch}


# ----------------------------------------------------------------------------
# W6: manifest (built before any signal is generated)
# ----------------------------------------------------------------------------
def cents_to_B(c):
    """B of a pure n³ injection with deviation c cents at n = 10:
    f_n = n f0 + d10 (n/10)³, d10 = 10 f0 (2^(c/1200) − 1)
    → a3 = d10/1000, B = 2 a3 / f0 = 0.02 (2^(c/1200) − 1)  (independent of f0)."""
    return 0.02 * (2.0 ** (c / 1200.0) - 1.0)


def cell_definitions():
    """The 36 cells of B0 §4, in fixed order. Each cell: dict with
    cell id, block, model parameters, injected sign, truth fields, m (B2)."""
    cells = []

    def add(cell, block, m, sign, gen, p_true=None, B_true=None, level=None,
            noise=NOISE_DEFAULT, jitter=JITTER_DEFAULT):
        cells.append({"cell": cell, "block": block, "m": m, "sign": sign, "gen": gen,
                      "p_true": p_true, "B_true": B_true, "level": level,
                      "noise_db": noise, "jitter_cents": jitter})

    # N: harmonic null
    add("N", "N", 2000, 0, {"model": "harmonic"})
    # P: pure power n^p, p ∈ {2, 3, 4}, d10 ∈ {0.5, 2, 5} cents, ±
    for p in (2, 3, 4):
        for c in (0.5, 2.0, 5.0):
            for s in (+1, -1):
                cs = s * c
                add(f"P_p{p}_c{c:g}_{'+' if s > 0 else '-'}", "P", 100, s,
                    {"model": "power", "p": p, "cents_at_10": cs},
                    p_true=p, B_true=(cents_to_B(cs) if p == 3 else None), level=c)
    # S: exact stiff string, |B| ∈ {1e-5, 1e-4, 4e-4, 1e-3}, ±
    for b in (1e-5, 1e-4, 4e-4, 1e-3):
        for s in (+1, -1):
            add(f"S_B{b:g}_{'+' if s > 0 else '-'}", "S", 100, s,
                {"model": "stiff", "B": s * b}, p_true=3, B_true=s * b, level=b)
    # O: out of family, pure n⁵, d10 = 2 cents, ±
    for s in (+1, -1):
        add(f"O_p5_c2_{'+' if s > 0 else '-'}", "O", 100, s,
            {"model": "power", "p": 5, "cents_at_10": s * 2.0}, p_true=5, level=2.0)
    # C: contamination: n³ (5 cents, sign s) + n² with opposite sign, half the
    #    n = 10 deviation (in Hz)
    for s in (+1, -1):
        add(f"C_mix_{'+' if s > 0 else '-'}", "C", 100, s,
            {"model": "mixture", "cents3": s * 5.0, "ratio2": -0.5},
            p_true=3, B_true=cents_to_B(s * 5.0), level=5.0)
    # R: robustness stratum (jitter 0.3 cent, noise −40 dB)
    add("R_N", "R", 500, 0, {"model": "harmonic"}, noise=NOISE_R, jitter=JITTER_R)
    for s in (+1, -1):
        add(f"R_P_p3_c5_{'+' if s > 0 else '-'}", "R", 100, s,
            {"model": "power", "p": 3, "cents_at_10": s * 5.0},
            p_true=3, B_true=cents_to_B(s * 5.0), level=5.0, noise=NOISE_R, jitter=JITTER_R)
    for s in (+1, -1):
        add(f"R_S_B0.0004_{'+' if s > 0 else '-'}", "R", 100, s,
            {"model": "stiff", "B": s * 4e-4}, p_true=3, B_true=s * 4e-4, level=4e-4,
            noise=NOISE_R, jitter=JITTER_R)
    return cells


EXPECTED_B2 = {"N": 2000, "P": 1800, "S": 800, "O": 200, "C": 200, "R": 900}
EXPECTED_CELLS = 36
EXPECTED_TOTAL_B2 = 5900


def build_manifest(mode):
    cells = cell_definitions()
    rows = []
    for c in cells:
        m = c["m"] if mode == "b2" else M_SMOKE
        for _ in range(m):
            rows.append(c)
    return cells, rows


def assert_manifest(mode, cells, rows):
    """W6: counts asserted against the B0 §4 table before any signal exists."""
    if len(cells) != EXPECTED_CELLS:
        sys.exit(f"ABORT: {len(cells)} cells != {EXPECTED_CELLS}")
    if len({c["cell"] for c in cells}) != EXPECTED_CELLS:
        sys.exit("ABORT: duplicate cell ids")
    per_block = {}
    for r in rows:
        per_block[r["block"]] = per_block.get(r["block"], 0) + 1
    if mode == "b2":
        if per_block != EXPECTED_B2 or len(rows) != EXPECTED_TOTAL_B2:
            sys.exit(f"ABORT: manifest counts {per_block} (total {len(rows)}) "
                     f"!= {EXPECTED_B2} (total {EXPECTED_TOTAL_B2})")
    else:
        if len(rows) != EXPECTED_CELLS * M_SMOKE:
            sys.exit(f"ABORT: smoke manifest has {len(rows)} cases")
    return per_block


# ----------------------------------------------------------------------------
# W5: generators (the scanner's synth_note unchanged; one added mixture generator)
# ----------------------------------------------------------------------------
def synth_mixture(f0, rng, cents3, ratio2, sr=44100, dur=1.3, noise_db=-60.0, jitter_cents=0.0):
    """n³ term (deviation cents3 at n = 10) plus an n² term whose n = 10 deviation
    in Hz is ratio2 times that of the n³ term. Everything after the frequency
    law is a line-by-line copy of piano_scanner.synth_note (same RNG order)."""
    nmax = min(80, int(0.45 * sr / f0))
    n = np.arange(1, nmax + 1, dtype=float)
    d10_3 = 10 * f0 * (2 ** (cents3 / 1200.0) - 1)
    d10_2 = ratio2 * d10_3
    fn = n * f0 + d10_3 * (n / 10.0) ** 3 + d10_2 * (n / 10.0) ** 2
    # --- copied from synth_note (e7c7520, l. 387–403) ---
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


def make_signal(cell, rng):
    """Per case: f0 draw, signal, then detuned nominal f0 — all from the case's own stream."""
    f0 = float(ps.random_f0s(rng, 1)[0])
    g = cell["gen"]
    kw = dict(noise_db=cell["noise_db"], jitter_cents=cell["jitter_cents"])
    with np.errstate(invalid="ignore"):          # B < 0: √(1 + B n²) is NaN above the cut-off;
        if g["model"] == "harmonic":             # those partials are dropped by synth_note's keep filter
            x, sr = ps.synth_note(f0, rng, model="harmonic", **kw)
        elif g["model"] == "power":
            x, sr = ps.synth_note(f0, rng, model="power", p=g["p"], cents_at_10=g["cents_at_10"], **kw)
        elif g["model"] == "stiff":
            x, sr = ps.synth_note(f0, rng, model="stiff", B=g["B"], **kw)
        elif g["model"] == "mixture":
            x, sr = synth_mixture(f0, rng, g["cents3"], g["ratio2"], **kw)
        else:
            raise ValueError(g["model"])
    f_nom = float(ps.detune(rng, f0))
    return x, sr, f0, f_nom


# ----------------------------------------------------------------------------
# Fits on exactly the scanner's design (fit_templates, l. 295–329)
# ----------------------------------------------------------------------------
def design(n, w, powers):
    n = np.asarray(n, float)
    x = n / n.max()
    return np.column_stack([x] + [x ** p for p in powers]) * w[:, None]


def wls(n, fw, w, powers):
    X = design(n, w, powers)
    c, *_ = np.linalg.lstsq(X, fw, rcond=None)
    rss = float(np.sum((fw - X @ c) ** 2))
    return c, rss, X


def ftest_globex(rss_small, rss_big, dof):
    """Exact copy of globex_scanner._ftest (e7c7520, l. 281–285)."""
    if rss_big <= 0 or dof <= 0:
        return np.inf, 0.0
    F = max(0.0, (rss_small - rss_big) / (rss_big / dof))
    return float(F), float(stats.f.sf(F, 1, dof))


def discrimination(n, fw, w, winner, rss_single, alpha):
    """W2 (B0 §3): T(p*|q) for every competitor q, Holm over competitors;
    distinguishable iff all adjusted p < alpha. Reverse tests T(q|p*) are
    recorded; mixture_like iff for some q both directions are significant
    (each family Holm-adjusted; the flag does not change the declared outcome)."""
    N = len(n)
    others = [q for q in ps.CONFIG["templates"] if q != winner]
    fwd, rev = {}, {}
    for q in others:
        _, rss_both, X = wls(n, fw, w, [winner, q])
        dof = N - X.shape[1]                                  # N − 3
        fwd[q] = ftest_globex(rss_single[q], rss_both, dof)
        rev[q] = ftest_globex(rss_single[winner], rss_both, dof)
    padj_f = ps.holm([fwd[q][1] for q in others])
    padj_r = ps.holm([rev[q][1] for q in others])
    out = {}
    for i, q in enumerate(others):
        out[str(q)] = {"F_winner_given_q": fwd[q][0], "p": fwd[q][1], "padj": float(padj_f[i]),
                       "F_q_given_winner": rev[q][0], "p_rev": rev[q][1], "padj_rev": float(padj_r[i])}
    distinguishable = all(out[str(q)]["padj"] < alpha for q in others)
    mixture_like = any(out[str(q)]["padj"] < alpha and out[str(q)]["padj_rev"] < alpha for q in others)
    return distinguishable, mixture_like, out


def b_hat_interval(n, fw, w):
    """W3: B̂ = 2 c1 / (c0 n_max²) from the p = 3 model, delta-method se,
    t interval with N − 2 dof. Returns the fit pieces used for parity."""
    n = np.asarray(n, float)
    nmax = n.max()
    c, rss, X = wls(n, fw, w, [3])
    N = len(n)
    s2 = rss / (N - 2)
    cov = s2 * np.linalg.inv(X.T @ X)
    c0, c1 = float(c[0]), float(c[1])
    B = 2.0 * c1 / (c0 * nmax ** 2)
    g = np.array([-2.0 * c1 / (c0 ** 2 * nmax ** 2), 2.0 / (c0 * nmax ** 2)])
    se = float(np.sqrt(g @ cov @ g))
    t = float(stats.t.ppf(0.975, N - 2))
    return {"B_hat": B, "se": se, "lo": B - t * se, "hi": B + t * se, "t": t,
            "a3": c1 / nmax ** 3, "f0_3": c0 / nmax}


def crosscheck_se(n, f, w, B_hat, se):
    """Smoke only (B0 W3 checks): independent se via scipy.optimize.curve_fit on
    the unscaled model f = A n + C n³, sigma = 1/w, absolute_sigma=False."""
    n = np.asarray(n, float)
    f = np.asarray(f, float)

    def model(nn, A, C):
        return A * nn + C * nn ** 3

    def jac(nn, A, C):
        # analytic Jacobian of the (linear) model; B1 v0.3: the smoke run showed that the
        # default finite-difference Jacobian limits the reference se to ~1e-6 … 1e-1 relative
        return np.column_stack([nn, nn ** 3])

    p0 = [float(np.median(f / n)), 0.0]          # generic start, not the WLS solution
    popt, pcov = optimize.curve_fit(model, n, f, p0=p0, sigma=1.0 / w, absolute_sigma=False,
                                    jac=jac, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15,
                                    maxfev=10000)
    A, C = popt
    Bc = 2.0 * C / A
    g = np.array([-2.0 * C / A ** 2, 2.0 / A])
    sec = float(np.sqrt(g @ pcov @ g))
    ok = (abs(Bc - B_hat) <= CROSSCHECK_RTOL * max(abs(B_hat), 1e-300)
          and abs(sec - se) <= CROSSCHECK_RTOL * max(abs(se), 1e-300))
    return {"B_curvefit": float(Bc), "se_curvefit": sec, "ok": bool(ok)}


# ----------------------------------------------------------------------------
# Analysis of one case (the analysis path sees only x, sr, f_nom)
# ----------------------------------------------------------------------------
class ParityError(Exception):
    pass


def _close(a, b, rtol=PARITY_RTOL):
    if a == b:
        return True
    if not (np.isfinite(a) and np.isfinite(b)):
        return False
    return abs(a - b) <= rtol * max(abs(a), abs(b))


def analyse_case(x, sr, f_nom, crosscheck=False):
    cfg = ps.CONFIG
    a = cfg["analysis"]
    st = cfg["stats"]
    templates = cfg["templates"]
    alpha = st["alpha_note"]

    # Original scanner (parity reference)
    ref, why = ps.analyze_signal(x, sr, f_nom)

    # Replicated path (identical steps of analyze_signal, l. 332–345)
    seg = ps.analysis_segment(x, sr, a)
    if seg is None:
        mine, my_why = None, "signal too short"
        partials = []
    else:
        partials = ps.track_partials(seg, sr, f_nom, a)
        n_all = len(partials)
        partials = [p for p in partials if p[0] <= a["fit_max_n"]]
        if len(partials) < a["min_partials"]:
            mine, my_why = None, f"too few partials ({len(partials)})"
        else:
            mine, my_why = True, "ok"

    if (ref is None) != (mine is None) or (ref is None and why != my_why):
        raise ParityError(f"exclusion mismatch: scanner '{why}' vs wrapper '{my_why}'")
    if ref is None:
        return {"excluded": True, "reason": why}

    n = [p[0] for p in partials]
    f = [p[1] for p in partials]
    snr = [p[2] for p in partials]
    fit = ps.fit_templates(n, f, templates, w=snr if a["weights"] == "snr" else None)

    # significance exactly as analyze_signal (l. 346–358)
    flags = {}
    for p in templates:
        v = fit["per"][p]
        rel = v["a"] * 10.0 ** p / (10.0 * v["f0"])
        cents10 = float(1200 * np.log2(max(1e-12, 1 + rel)))
        flags[p] = bool(v["padj"] < alpha and abs(cents10) >= st["min_effect_cents"])
    winner = fit["winner"]
    significant = flags[winner]

    # parity guard
    if winner != ref["winner"] or significant != ref["significant"]:
        raise ParityError(f"winner/significant mismatch: {winner},{significant} vs "
                          f"{ref['winner']},{ref['significant']}")
    for p in templates:
        for key in ("p", "padj"):
            if not _close(fit["per"][p][key], ref["per"][p][key]):
                raise ParityError(f"template {p} {key}: {fit['per'][p][key]} vs {ref['per'][p][key]}")

    # weights and weighted data exactly as fit_templates (l. 302–310)
    w = np.asarray(snr, float)
    w = w / np.median(w)
    fw = np.asarray(f, float) * w
    rss_single = {p: fit["per"][p]["rss"] for p in templates}

    # W2
    disting, mixture_like, disc = discrimination(n, fw, w, winner, rss_single, alpha)

    # W3 (+ parity of B̂ with the scanner's B_est and per[3])
    bh = b_hat_interval(n, fw, w)
    if not (_close(bh["a3"], ref["per"][3]["a"]) and _close(bh["f0_3"], ref["per"][3]["f0"])
            and _close(bh["B_hat"], ref["B_est"])):
        raise ParityError("W3 coefficients differ from the scanner's per[3] / B_est")

    # W1 / W4
    a_win = float(fit["per"][winner]["a"])
    # three-valued signs (B1 review A): an exact zero is neither + nor −,
    # so it can never count as a correct sign for either injection
    sign_winner = int(np.sign(a_win))
    sign_B = int(np.sign(bh["B_hat"]))
    if not significant:
        declared = "none"
    else:
        declared = (f"n^{winner}{'+' if sign_winner > 0 else '-' if sign_winner < 0 else '0'} "
                    + ("distinguishable" if disting else "undetermined"))

    out = {"excluded": False, "reason": "ok",
           "n_partials": n_all, "n_fit": len(n),
           "winner": int(winner), "significant": bool(significant),
           "distinguishable": bool(disting), "mixture_like": bool(mixture_like),
           "declared": declared,
           "a_winner": a_win, "sign_winner": sign_winner,
           "B_hat": bh["B_hat"], "B_se": bh["se"], "B_lo": bh["lo"], "B_hi": bh["hi"],
           "sign_B": sign_B,
           "per": {str(p): {"a": float(fit["per"][p]["a"]), "p": float(fit["per"][p]["p"]),
                            "padj": float(fit["per"][p]["padj"]), "flag": flags[p]}
                   for p in templates},
           "discrimination": disc}
    if crosscheck:
        try:
            out["crosscheck"] = crosscheck_se(n, f, w, bh["B_hat"], bh["se"])
        except Exception as e:                       # smoke only: record, do not abort
            out["crosscheck"] = {"ok": False, "error": repr(e)}
    return out


# ----------------------------------------------------------------------------
# The run (B0 §6.2)
# ----------------------------------------------------------------------------
def run(mode, out_dir):
    if os.path.exists(out_dir) and os.listdir(out_dir):
        sys.exit(f"ABORT: output directory {out_dir} is not empty")
    os.makedirs(out_dir, exist_ok=True)
    frozen = verify_frozen_inputs()

    seed = SEED_B2 if mode == "b2" else SEED_SMOKE
    cells, rows = build_manifest(mode)
    per_block = assert_manifest(mode, cells, rows)           # before any signal exists
    n_cases = len(rows)
    n_spawn = N_SPAWN_B2 if mode == "b2" else n_cases + 1
    children = np.random.SeedSequence(seed).spawn(n_spawn)
    order = np.random.default_rng(children[n_cases]).permutation(n_cases)

    manifest = {
        "spec": SPEC_NAME, "spec_sha256": SPEC_SHA256, "mode": mode, "seed": seed,
        "n_cases": n_cases, "n_spawn": n_spawn, "shuffle_spawn_index": n_cases,
        "wrapper_sha256": sha256_file(os.path.abspath(__file__)),
        "per_block": per_block,
        "cells": [dict({k: c[k] for k in ("cell", "block", "sign", "p_true", "B_true", "level",
                                          "noise_db", "jitter_cents")},
                       gen=c["gen"], m=(c["m"] if mode == "b2" else M_SMOKE)) for c in cells],
        **frozen,
        "numpy": np.__version__, "python": sys.version.split()[0],
    }
    with open(os.path.join(out_dir, "manifest.json"), "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=1)

    log = open(os.path.join(out_dir, "run.log"), "w", encoding="utf-8")

    def logp(msg):
        print(msg, flush=True)
        log.write(msg + "\n")
        log.flush()

    logp(f"mode {mode}  seed {seed}  cases {n_cases}  spec {SPEC_SHA256[:16]}  "
         f"scanner {frozen['piano_scanner_sha256'][:16]}  config {frozen['config_hash']}")
    ftruth = open(os.path.join(out_dir, "truth.jsonl"), "w", encoding="utf-8")
    fana = open(os.path.join(out_dir, "analysis.jsonl"), "w", encoding="utf-8")
    t_start = time.time()
    n_excl = 0
    n_xc_fail = 0
    try:
        for pos, idx in enumerate(order):
            cell = rows[idx]
            case_id = f"k{pos:05d}"
            rng = np.random.default_rng(children[idx])
            x, sr, f0, f_nom = make_signal(cell, rng)
            ftruth.write(json.dumps({
                "case_id": case_id, "spawn_index": int(idx), "cell": cell["cell"],
                "block": cell["block"], "sign": cell["sign"], "p_true": cell["p_true"],
                "B_true": cell["B_true"], "level": cell["level"], "f0": f0,
                "noise_db": cell["noise_db"], "jitter_cents": cell["jitter_cents"]}) + "\n")
            t0 = time.perf_counter()
            try:
                res = analyse_case(x, sr, f_nom, crosscheck=(mode == "smoke"))
            except ParityError as e:
                logp(f"ABORT: parity guard failed at {case_id}: {e}")
                raise
            res["case_id"] = case_id
            res["f_nom"] = f_nom
            res["seconds"] = time.perf_counter() - t0
            n_excl += int(res["excluded"])
            if mode == "smoke" and not res["excluded"] and not res["crosscheck"]["ok"]:
                n_xc_fail += 1
            fana.write(json.dumps(res, default=float) + "\n")
            if (pos + 1) % 100 == 0 or pos + 1 == n_cases:
                el = time.time() - t_start
                logp(f"  {pos + 1}/{n_cases}  excluded {n_excl}  elapsed {el:.0f} s  "
                     f"({el / (pos + 1):.2f} s/case)")
    finally:
        ftruth.close()
        fana.close()
        if sys.exc_info()[0] is not None:
            log.close()
    logp(f"done: {n_cases} cases, excluded {n_excl}, wall {time.time() - t_start:.0f} s")
    if mode == "smoke":
        logp(f"W3 covariance cross-check failures: {n_xc_fail}")
    log.close()

    hashes = {name: sha256_file(os.path.join(out_dir, name))
              for name in ("manifest.json", "truth.jsonl", "analysis.jsonl", "run.log")}
    with open(os.path.join(out_dir, "hashes.json"), "w", encoding="utf-8") as fh:
        json.dump(hashes, fh, indent=1)
    for k, v in hashes.items():
        print(f"sha256 {v}  {k}")
    return 0


def main():
    ap = argparse.ArgumentParser(description="B1 wrapper: blind k⁴-equivalent recovery test")
    ap.add_argument("mode", choices=["smoke", "b2"])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    return run(args.mode, args.out)


if __name__ == "__main__":
    sys.exit(main())
