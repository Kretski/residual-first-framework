#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
piano_scanner.py — Модул 0 (пиано) на търсачката на аномалии в остатъци
=========================================================================

Идея:
  Теория, която изваждаме: идеална струна, f_n = n * f0 (f0 се фитва).
  Остатъци:                 Δf_n = f_n - n * f0_fit   [Hz]
  Шаблони (фиксирани):      Δf_n ∝ n^p,  p ∈ {0, 2, 3, 4}
                            (p = 1 се поглъща изцяло от f0 и не се тества)
  Очакване при пиано:       твърда струна → f_n = n f0 √(1 + B n²)
                            → Δf_n ≈ (f0 B / 2) n³  → трябва да спечели p = 3.
  Нулев тест:               звук с точно кратни обертонове (синтетичен, или
                            цигулка/виолончело със смичок, без вибрато)
                            → търсачката трябва да каже „нищо няма".

Режими:
  synth-null    калибровъчна бариера върху синтетични хармонични тонове
  synth-inject  чувствителност и разпознаване на формата (правилна и грешна)
  real          реални записи (напр. University of Iowa MIS) с разделяне
                търсене / потвърждение и незадължителен реален нулев тест

Примери:
  python piano_scanner.py synth-null   --out out_null
  python piano_scanner.py synth-inject --out out_inject
  python piano_scanner.py real --data Piano_mf --null-data Violin_nonvib --out out_real

Зависимости: numpy, scipy; за AIFF файлове и soundfile (pip install soundfile).
Физическото тълкуване (оценка на B) е отделено и се отпечатва САМО като
вторична информация — търсачката не е „детектор на B".
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
# ФИКСИРАНА КОНФИГУРАЦИЯ — не се променя след предварителната регистрация.
# Хешът ѝ се записва във всеки изход, за да може да се цитира в регистрацията.
# ----------------------------------------------------------------------------
CONFIG = {
    "version": "0.3.0",
    "analysis": {
        "skip_after_peak_s": 0.10,     # пропуск след удара (преходен процес)
        "window_s": 1.0,               # дължина на анализирания участък
        "min_window_s": 0.3,           # по-къс участък → нотата се изключва
        "zero_pad": 8,                 # допълване с нули при FFT
        "max_partials": 60,            # максимален номер на обертон
        "fit_max_n": 25,               # във фита влизат само n ≤ fit_max_n
        "min_partials": 8,             # по-малко обертона → нотата се изключва
        "seed_tol_rel": 0.06,          # търсене на първите обертони: ±6 %
        "seed_max_n": 3,               # стартов обертон може да е n = 1..3
        "track_halfwidth_rel": 0.25,   # прозорец при проследяване: ±0.25 интервал
        "snr_db": 20.0,                # праг за връх над медианния фон
        "max_consecutive_misses": 3,   # спиране след толкова пропуснати поред
        "f_max_hz": 12000.0,           # горна граница на честотата
        "weights": "snr",              # тегла във фита: σ_n ∝ 1/SNR
    },
    "templates": [0, 2, 3, 4],
    "stats": {
        "alpha_note": 0.01,            # праг на ниво нота (след Holm по шаблоните)
        "barrier_test_alpha": 0.05,    # бариера: СТОП, ако броят фалшиви аларми е
                                       # значимо над α (едностранен биномен тест)
        "group_alpha": 0.001,          # праг за структура в групата за търсене
        "confirm_alpha": 0.01,         # праг в групата за потвърждение
        "shape_majority": 0.6,         # формата трябва да печели в ≥ 60 % от нотите
    },
    "groups": {"rule": "midi_parity", "search": "even", "confirm": "odd"},
    "primary_dynamic": "mf",
}


def config_hash(cfg=CONFIG):
    s = json.dumps(cfg, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(s).hexdigest()[:16]


# ----------------------------------------------------------------------------
# Имена на файлове и ноти
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
    """Връща (midi, име на нотата, динамика) или None."""
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
    """Моно float сигнал и честота на дискретизация."""
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
            raise RuntimeError("За AIFF файлове е нужен soundfile: pip install soundfile")
    x = np.asarray(x, dtype=float)
    peak = np.max(np.abs(x))
    if peak > 0:
        x = x / peak
    return x, int(sr)


# ----------------------------------------------------------------------------
# Спектър и проследяване на обертоновете
# ----------------------------------------------------------------------------
def analysis_segment(x, sr, a):
    i_peak = int(np.argmax(np.abs(x)))
    start = i_peak + int(a["skip_after_peak_s"] * sr)
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
    if k == i0 or k == i1:                         # на ръба → не е истински връх
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


def track_partials(seg, sr, f_nom, a):
    """
    Проследява обертоновете без да предполага модел на отклонението:
    следващият се търси около линейна екстраполация на последните два намерени.
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
    return found


# ----------------------------------------------------------------------------
# Остатъци, шаблони и статистика
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
    Базов модел f = f0·n; за всеки шаблон: f = f0·n + a·n^p и F-тест.
    w = тегла 1/σ_n. Грешката по честота на връх е ∝ 1/SNR (проверено на
    синтетични данни), затова w = SNR на върха.
    """
    n = np.asarray(n, float)
    f = np.asarray(f, float)
    N = len(n)
    w = np.ones(N) if w is None else np.asarray(w, float)
    w = w / np.median(w)
    nmax = n.max()
    x = n / nmax                                       # мащабиране за стабилност
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
        return None, "кратък сигнал"
    partials = track_partials(seg, sr, f_nom, a)
    n_all = len(partials)
    partials = [p for p in partials if p[0] <= a["fit_max_n"]]
    if len(partials) < a["min_partials"]:
        return None, f"малко обертонове ({len(partials)})"
    n = [p[0] for p in partials]
    f = [p[1] for p in partials]
    snr = [p[2] for p in partials]
    fit = fit_templates(n, f, cfg["templates"], w=snr if cfg["analysis"]["weights"] == "snr" else None)
    alpha = cfg["stats"]["alpha_note"]
    w = fit["winner"]
    res = {
        "n_partials": n_all, "n_fit": len(n), "f0_fit": fit["f0_fit"],
        "winner": w,
        "significant": bool(fit["per"][w]["padj"] < alpha),
        "flags": {p: bool(fit["per"][p]["padj"] < alpha) for p in cfg["templates"]},
        "per": fit["per"], "n": n, "f": f, "resid": fit["resid"],
        # вторично тълкуване (отделен модул): твърда струна, Δf ≈ (f0 B/2) n³
        "B_est": 2.0 * fit["per"][3]["a"] / fit["per"][3]["f0"] if 3 in fit["per"] else np.nan,
    }
    return res, "ok"


# ----------------------------------------------------------------------------
# Синтетични тонове
# ----------------------------------------------------------------------------
def synth_note(f0, rng, sr=44100, dur=1.3, model="harmonic", B=0.0,
               p=3, cents_at_10=0.0, noise_db=-60.0, jitter_cents=0.0):
    """
    model = 'harmonic' : f_n = n f0
            'stiff'    : f_n = n f0 √(1 + B n²)         (точна твърда струна)
            'power'    : f_n = n f0 + a n^p, където отклонението при n = 10
                         е cents_at_10 цента
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
    """Номиналната честота се различава от истинската (както при реално пиано)."""
    return f0 * 2 ** (rng.normal(0, cents) / 1200.0)


# ----------------------------------------------------------------------------
# Групов анализ
# ----------------------------------------------------------------------------
def barrier(results, cfg=CONFIG):
    """Калибровъчна бариера. Връща (FPR, биномно p, преминала ли е)."""
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
    print(f"ноти в анализа: {g['M']}")
    if g["M"] == 0:
        return
    for p, v in g["per_template"].items():
        print(f"  шаблон n^{p}: значими ноти {v['k']:3d}/{g['M']}   биномно p = {v['binom_p']:.2e}")
    print(f"  значими ноти общо: {g['n_significant']}")
    print(f"  групово p (Бонферони по шаблоните): {g['group_p']:.2e}")
    shape = "неопределена" if g["shape"] is None else f"n^{g['shape']}"
    print(f"  форма: {shape} (дял сред значимите: {g['shape_share']:.0%})")


# ----------------------------------------------------------------------------
# Изход
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
        cols += [f"F_{p}", f"p_{p}", f"padj_{p}", f"a_{p}"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for meta, r in rows:
            line = [meta.get(c, "") for c in cols[:6]]
            line += [r["n_partials"], r["n_fit"], f"{r['f0_fit']:.6f}", r["winner"],
                     int(r["significant"]), f"{r['B_est']:.3e}"]
            for p in cfg["templates"]:
                v = r["per"][p]
                line += [f"{v['F']:.4g}", f"{v['p']:.3e}", f"{v['padj']:.3e}", f"{v['a']:.4e}"]
            w.writerow(line)


def save_plots(out, rows, k):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("(matplotlib липсва — графиките се пропускат)")
        return
    for meta, r in rows[:k]:
        n = np.array(r["n"], float)
        w = r["winner"]
        v = r["per"][w]
        model = v["f0"] * n + v["a"] * n ** w - r["f0_fit"] * n
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(n, r["resid"], "o", label="остатъци Δf_n")
        ax.plot(n, model, "-", label=f"най-добър шаблон n^{w}")
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel("номер на обертона n")
        ax.set_ylabel("Δf_n = f_n − n·f0 [Hz]")
        ax.set_title(f"{meta.get('label', '')}  (значимо: {r['significant']})")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(out, f"resid_{meta.get('label', 'note')}.png"), dpi=120)
        plt.close(fig)


# ----------------------------------------------------------------------------
# Режими
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
    print(f"хеш на конфигурацията: {config_hash()}")
    print(f"синтетичен нулев тест: {len(results)} ноти (изключени {skipped})")
    for p in CONFIG["templates"]:
        print(f"  FPR шаблон n^{p}: {np.mean([r['flags'][p] for r in results]):.3f}")
    print(f"FPR на ниво нота (Holm): {fpr:.3f}  (очаквано ≤ {CONFIG['stats']['alpha_note']}), "
          f"биномно p = {pb:.3f}")
    print("БАРИЕРА: ПРЕМИНАТА" if passed else "БАРИЕРА: НЕ Е ПРЕМИНАТА → СТОП")
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
    print(f"хеш на конфигурацията: {config_hash()}")
    print("\nИнжекции с форма n^p (амплитуда = отклонение при n = 10, в цента)")
    print("вярна форма | цента | открити | правилна форма сред откритите")
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
    print("\nТочна твърда струна f_n = n f0 √(1+B n²) — трябва да печели n^3")
    print("      B     | открити | правилна форма (n^3) | медиана B_est / B")
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
    print(f"хеш на конфигурацията: {config_hash()}")
    summary = {"config_hash": config_hash()}

    # 1) Реален нулев тест (калибровъчна бариера)
    if args.null_data:
        nfiles = collect_files(args.null_data, "all")
        nrows, nskip = analyze_files(nfiles, "null_")
        nres = [r for _, r in nrows]
        fpr, pb, passed = barrier(nres)
        print(f"\nреален нулев тест: {len(nres)} ноти (изключени {len(nskip)})")
        print(f"FPR на ниво нота: {fpr:.3f}  (очаквано ≤ {CONFIG['stats']['alpha_note']}), "
              f"биномно p = {pb:.3f}")
        write_notes_csv(os.path.join(args.out, "notes_null.csv"), nrows)
        summary["null"] = {"M": len(nres), "fpr": fpr, "binom_p": pb}
        if not passed:
            print("БАРИЕРА: НЕ Е ПРЕМИНАТА → СТОП. Анализът на пианото не се пуска.")
            summary["stopped"] = True
            with open(os.path.join(args.out, "summary.json"), "w", encoding="utf-8") as fh:
                json.dump(summary, fh, indent=2, default=float)
            return 2
        print("БАРИЕРА: ПРЕМИНАТА")
    else:
        print("\n(няма реален нулев тест — синтетичната бариера трябва да е минала)")

    # 2) Пиано: разделяне по фиксирано правило
    dyn = args.dynamic or CONFIG["primary_dynamic"]
    files = collect_files(args.data, dyn)
    rows, skipped = analyze_files(files)
    for meta, _ in rows:
        meta["group"] = "search" if meta["midi"] % 2 == 0 else "confirm"
    for name, why in skipped:
        print(f"  изключена: {name} ({why})")
    search = [r for m, r in rows if m["group"] == "search"]
    confirm = [r for m, r in rows if m["group"] == "confirm"]
    gs, gc = group_summary(search), group_summary(confirm)
    print_group("ГРУПА ЗА ТЪРСЕНЕ (четни MIDI)", gs)
    print_group("ГРУПА ЗА ПОТВЪРЖДЕНИЕ (нечетни MIDI)", gc)

    st = CONFIG["stats"]
    found = gs["M"] > 0 and gs["group_p"] < st["group_alpha"] and gs["shape"] is not None
    confirmed = False
    if found:
        s = gs["shape"]
        confirmed = (gc["M"] > 0 and gc["shape"] == s
                     and gc["per_template"][s]["binom_p"] < st["confirm_alpha"])
    print("\n=== РЕЗУЛТАТ ===")
    if not found:
        print("в групата за търсене няма значима структура с ясна форма")
    else:
        print(f"търсене: структура с форма n^{gs['shape']}")
        print("потвърждение: " + ("ДА — същата форма е значима" if confirmed else "НЕ"))

    # 3) Вторично тълкуване (отделен модул): коефициент на нехармоничност
    B = [r["B_est"] for r in search + confirm if r["significant"] and r["winner"] == 3]
    if B:
        print(f"\n[тълкуване] B при твърда струна: медиана {np.median(B):.2e}, "
              f"диапазон {np.min(B):.1e} … {np.max(B):.1e} ({len(B)} ноти)")

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
    ap = argparse.ArgumentParser(description="Модул 0 (пиано) на търсачката на остатъци")
    sub = ap.add_subparsers(dest="mode", required=True)

    a = sub.add_parser("synth-null", help="калибровъчна бариера върху синтетични тонове")
    a.add_argument("--m", type=int, default=300)
    a.add_argument("--noise-db", type=float, default=-60.0)
    a.add_argument("--jitter", type=float, default=0.0, help="случайно разместване [цента]")
    a.add_argument("--seed", type=int, default=1)
    a.add_argument("--out", default="out_synth_null")

    b = sub.add_parser("synth-inject", help="чувствителност и разпознаване на формата")
    b.add_argument("--m", type=int, default=60)
    b.add_argument("--cents", type=float, nargs="+", default=[0.5, 1, 2, 5, 20])
    b.add_argument("--B", type=float, nargs="+", default=[1e-5, 1e-4, 4e-4, 1e-3])
    b.add_argument("--noise-db", type=float, default=-60.0)
    b.add_argument("--jitter", type=float, default=0.0)
    b.add_argument("--seed", type=int, default=2)
    b.add_argument("--out", default="out_synth_inject")

    c = sub.add_parser("real", help="реални записи")
    c.add_argument("--data", required=True, help="папка с нотите на пианото")
    c.add_argument("--null-data", help="папка с хармоничен инструмент (нулев тест)")
    c.add_argument("--dynamic", choices=list(DYNAMICS) + ["all"])
    c.add_argument("--plots", type=int, default=0, help="брой графики на остатъци")
    c.add_argument("--out", default="out_real")

    args = ap.parse_args()
    if args.mode == "synth-null":
        return run_synth_null(args)
    if args.mode == "synth-inject":
        return run_synth_inject(args)
    return run_real(args)


if __name__ == "__main__":
    sys.exit(main())
