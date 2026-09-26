"""
EPTA DR2 replication of the frozen V2.1 HD test — follows PREREGISTRATION_EPTA.md.

Run INSIDE WSL/Linux (80-bit long double), from the project folder:

  python run_epta.py check J1909-3744     # 1. load ONE pulsar's TOAs with PINT (no residuals)
  python run_epta.py gate                 # 2. calibration gate: sampling/errors/positions only
  python run_epta.py real                 # 3. only if the gate PASSED: residuals + frozen V2.1
  python run_epta.py real --variant DR2new    (secondary dataset, after the primary)

The frozen analysis module (common_residual_search_v2.py, tag v2.1-epta-prereg)
is imported, never modified. Nothing here changes its statistics or defaults.

Author: Dimitar Kretski
License: MIT
"""
import os
import sys
import json
import hashlib
import datetime
import contextlib
from pathlib import Path

import numpy as np

ROOT = Path("EPTA_DR2")
CACHE = Path("epta_cache")
GATE_FILE = "epta_gate_{variant}.json"
GATE_THRESHOLD = 0.09            # pre-registered (Section 4)
PRIMARY = "hd_perp_mono_dipole:scramble"
FROZEN_SHA = "089a94e073c962eb332779e0cc50807ebbb30a6e3f02bb49bced0226497ec2b2"


# ---------------------------------------------------------------- environment
def check_precision():
    eps = float(np.finfo(np.longdouble).eps)
    if eps > 1e-18:
        sys.exit(f"long double eps = {eps:.3g}: NOT extended precision. "
                 f"Run this inside WSL/Linux (pre-registration, Amendment 1.2).")
    return eps


def frozen_hash_ok():
    """Compare with line endings normalised (git on Windows may add CRLF)."""
    raw = Path("common_residual_search_v2.py").read_bytes().replace(b"\r\n", b"\n")
    h = hashlib.sha256(raw).hexdigest()
    return h == FROZEN_SHA, h


def find_variant_dir(variant):
    hits = [p for p in ROOT.rglob(variant) if p.is_dir() and p.name == variant]
    if not hits:
        sys.exit(f"Folder '{variant}' not found under {ROOT} — run fetch_epta_dr2.py first.")
    return hits[0]


def setup_pint():
    clk = [p for p in ROOT.rglob("*") if p.is_dir() and p.name.lower().startswith("clock")]
    if clk:
        os.environ["PINT_CLOCK_OVERRIDE"] = str(clk[0].resolve())   # corrected NRT clock (README)
    import pint.logging
    pint.logging.setup(level="WARNING")
    return str(clk[0]) if clk else None


@contextlib.contextmanager
def cwd(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


# ---------------------------------------------------------------- PINT loading
def psr_files(psr_dir):
    name = psr_dir.name
    par = psr_dir / f"{name}.par"
    tims = sorted(psr_dir.glob("*_all.tim")) or sorted(psr_dir.glob("*.tim"))
    if not par.exists() or not tims:
        raise FileNotFoundError(f"{name}: par={par.exists()}, tim files={len(tims)}")
    return par, tims[0]


TEMPO2_COMMANDS = {"FORMAT", "MODE", "JUMP", "TIME", "PHASE", "EFAC", "EQUAD", "EMIN", "EMAX",
                   "FMIN", "FMAX", "INFO", "SKIP", "NOSKIP", "TRACK"}


def _is_comment(tok):
    """'C', 'c', 'C??', 'c!!', '#…' — C/c followed only by punctuation, or '#'."""
    t = tok[0]
    return t.startswith("#") or (t[0] in "Cc" and all(not ch.isalnum() for ch in t[1:]))


def _is_toa_line(tok):
    """A FORMAT 1 TOA line: name freq MJD error site ... with numeric freq/MJD/error."""
    if _is_comment(tok) or len(tok) < 5:
        return False
    try:
        float(tok[1]); float(tok[2]); float(tok[3])
        return True
    except ValueError:
        return False


def split_tim_tempo2(tim, out_dir):
    """
    tempo2 semantics, one file at a time (Amendments 4 / 4c):
      * INCLUDE is followed recursively; END ends only the file it is in.
      * FORMAT is INHERITED from the including file (EPTA sub-files often have
        no FORMAT line of their own) — every copy gets 'FORMAT 1' if any
        ancestor or itself declared it.
      * A line is kept if it is a valid FORMAT 1 TOA line or a known tempo2
        command. Everything else (comments 'C', 'c', 'C??', '#', stray text)
        is dropped — tempo2 cannot parse such lines as TOAs either — and COUNTED.
    Returns [(copy path, original name, TOAs kept, lines dropped, lines after END)].
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    leaves = []
    state = {"fmt1": False}     # tempo2/PINT: FORMAT is GLOBAL state in reading order (Amendment 4d)

    def read(path, depth=0):
        fmt1_at_start = state["fmt1"]
        own, n_toa, dropped, skipped, ended = [], 0, 0, 0, False
        for raw in open(path, errors="replace"):
            line = raw.rstrip("\n")
            tok = line.split()
            if not tok:
                continue
            key = tok[0].upper()
            if ended:
                skipped += 1
                continue
            if key == "END":
                ended = True
                continue
            if key == "FORMAT":
                state["fmt1"] = len(tok) > 1 and tok[1] == "1"
                own.append(line)                      # kept in place: applies to the lines after it
                continue
            if key == "INCLUDE" and len(tok) > 1 and depth < 5:
                p = Path(tok[1])
                read(p if p.is_absolute() else (path.parent / p).resolve(), depth + 1)
                continue
            if _is_toa_line(tok):
                # FORMAT 1 is whitespace-delimited, so leading blanks carry no
                # meaning — but PINT classifies ' …' lines with '.' at column 41
                # as fixed-column Parkes format before honouring FORMAT 1.
                # Strip them (Amendment 4e). Non-FORMAT-1 lines are kept raw.
                own.append(line.strip() if state["fmt1"] else line)
                n_toa += 1
            elif key in TEMPO2_COMMANDS - {"FORMAT"}:
                own.append(line)
            else:
                dropped += 1
        if n_toa:
            dst = out_dir / f"{len(leaves):03d}_{path.name}"
            dst.write_text(("FORMAT 1\n" if fmt1_at_start else "") + "\n".join(own) + "\n")
            leaves.append((dst, path.name, n_toa, dropped, skipped))

    read(Path(tim).resolve())
    return leaves


def read_leaf(pint_toa, path, model):
    try:
        return pint_toa.get_TOAs(str(path), model=model)
    except Exception as e1:
        txt = Path(path).read_text()
        bad = str(e1).split(": ")[-1].strip("'\" ")
        where = next((f"line {i + 1}: {l.strip()[:160]}" for i, l in enumerate(txt.splitlines())
                      if bad and bad in l), "offending line not located")
        raise RuntimeError(f"{Path(path).name}: {type(e1).__name__}: {e1}\n      → {where}")


def load_pint(psr_dir, with_residuals, verbose=False):
    import pint.models
    import pint.toa
    par, tim = psr_files(psr_dir)
    leaves = split_tim_tempo2(psr_dir / tim.name,
                              CACHE / f"{psr_dir.parent.name}_{psr_dir.name}_tims")
    if verbose:
        for _, name, k, dr, sk in leaves:
            print(f"    {name:<24s} TOAs {k:5d}" + (f", dropped (comment/non-TOA) {dr}" if dr else "")
                  + (f", after END {sk}" if sk else ""))
    with cwd(psr_dir):
        try:
            model = pint.models.get_model(par.name, allow_tcb=True, allow_T2=True)
        except TypeError:                     # older PINT without these keywords
            model = pint.models.get_model(par.name)
    parts = [read_leaf(pint.toa, p, model) for p, *_ in leaves]
    toas = parts[0] if len(parts) == 1 else pint.toa.merge_TOAs(parts)
    out = {
        "mjd": np.asarray(toas.get_mjds().value, dtype=float),
        "err_s": np.asarray(toas.get_errors().to_value("s"), dtype=float),
        "freq": np.asarray(toas.get_freqs().to_value("MHz"), dtype=float),
        # backend = '-group' flag (NUPPI '-sys' differs per frequency channel)
        "backend": np.array([g if g is not None else (sy if sy is not None else o)
                             for g, sy, o in zip(toas.get_flag_value("group")[0],
                                                 toas.get_flag_value("sys")[0],
                                                 toas.get_obss())]),
    }
    if with_residuals:
        import pint.residuals
        r = pint.residuals.Residuals(toas, model, subtract_mean=True)
        out["res_s"] = np.asarray(r.time_resids.to_value("s"), dtype=float)
    return out


def epoch_average(d):
    """Weighted average per (observing day, backend) — Amendment 1.5."""
    key = np.array([f"{int(np.floor(m))}|{b}" for m, b in zip(d["mjd"], d["backend"])])
    uniq, inv = np.unique(key, return_inverse=True)
    w = 1.0 / d["err_s"] ** 2
    W = np.bincount(inv, weights=w)
    out = {"mjd": np.bincount(inv, weights=w * d["mjd"]) / W,
           "err_s": 1.0 / np.sqrt(W),
           "freq": np.bincount(inv, weights=w * d["freq"]) / W,
           "n_raw": np.bincount(inv)}
    if "res_s" in d:
        out["res_s"] = np.bincount(inv, weights=w * d["res_s"]) / W
    order = np.argsort(out["mjd"])
    return {k: v[order] for k, v in out.items()}


class Psr:                                    # interface expected by the frozen V2.1 module
    def __init__(self, name, toas, residuals, errors, ra, dec):
        self.name, self.toas, self.residuals, self.errors = name, toas, residuals, errors
        self.ra, self.dec = ra, dec
        self.time_span = (toas.max() - toas.min()) / 365.25


def load_variant(variant, with_residuals):
    from nanograv_loader import read_par_file
    vdir = find_variant_dir(variant)
    CACHE.mkdir(exist_ok=True)
    tag = "res" if with_residuals else "samp"
    pulsars, failed = {}, {}
    for psr_dir in sorted(p for p in vdir.iterdir() if p.is_dir()):
        name = psr_dir.name
        cache = CACHE / f"{variant}_{name}_{tag}_a4e.npz"     # a4e: Amendment-4e reader
        try:
            if cache.exists():
                d = dict(np.load(cache))
            else:
                print(f"  PINT: {name} …", flush=True)
                d = epoch_average(load_pint(psr_dir, with_residuals))
                np.savez(cache, **d)
            par = read_par_file(psr_files(psr_dir)[0])
            res = d["res_s"] if with_residuals else np.zeros_like(d["mjd"])
            pulsars[name] = Psr(name, d["mjd"], res, d["err_s"], par["ra"], par["dec"])
        except Exception as e:                # pre-registered: list and exclude
            failed[name] = f"{type(e).__name__}: {e}"
    return pulsars, failed


# ---------------------------------------------------------------- gate simulator
def synthetic_h0_array(rng, template):
    """
    H_0 simulator for the gate (Amendment 2). Same model as the frozen module's
    synthetic_array(): white noise at the real TOA errors, intrinsic red noise
    (1 µs RMS, β = 3) on every third pulsar, weighted quadratic removal — but
    with the time grid spanning the ACTUAL data span (the frozen helper assumes
    16 yr; EPTA spans up to ~25 yr, beyond which its red noise would be flat).
    """
    import common_residual_search_v2 as v
    names = sorted(template)
    t0 = min(p.toas.min() for p in template.values())
    T = max(p.toas.max() for p in template.values()) - t0
    grid = np.linspace(0, T, 8192)
    out = {}
    for i, n in enumerate(names):
        p = template[n]
        t, e = p.toas - t0, p.errors
        r = e * rng.standard_normal(len(t))
        if i % 3 == 0:
            r += v._red(rng, t, 1.0e-6, 3.0, grid)
        V = np.vander((t - t.mean()) / T, 3)
        r -= V @ np.linalg.lstsq(V / e[:, None], r / e, rcond=None)[0]
        out[n] = Psr(n, t, r, e, p.ra, p.dec)
    return out


# ---------------------------------------------------------------- stages
def stage_check(psr):
    check_precision()
    setup_pint()
    psr_dir = find_variant_dir("DR2full") / psr
    d = load_pint(psr_dir, with_residuals=False, verbose=True)
    a = epoch_average(d)
    print(f"{psr}: {len(d['mjd'])} TOAs → {len(a['mjd'])} epoch-averaged "
          f"(day × backend); MJD {d['mjd'].min():.0f}–{d['mjd'].max():.0f} "
          f"({(d['mjd'].max() - d['mjd'].min()) / 365.25:.1f} yr)")
    print(f"  backends: {sorted(set(d['backend']))}")
    print(f"  median TOA error: {np.median(d['err_s']) * 1e6:.3f} µs; "
          f"frequencies {d['freq'].min():.0f}–{d['freq'].max():.0f} MHz")
    print("  (no residuals computed — sampling only)")


def stage_gate(variant, trials):
    import common_residual_search_v2 as v
    eps = check_precision()
    ok, h = frozen_hash_ok()
    clk = setup_pint()
    print(f"Calibration gate — EPTA {variant} (sampling, errors, positions only)")
    print(f"  long double eps {eps:.3g}; frozen V2.1 hash {'OK' if ok else 'MISMATCH ' + h}")
    pulsars, failed = load_variant(variant, with_residuals=False)
    pulsars = {n: p for n, p in pulsars.items() if np.isfinite(p.ra) and np.isfinite(p.dec)}
    print(f"  pulsars loaded: {len(pulsars)}; excluded: {len(failed)}")
    for n, e in failed.items():
        print(f"    - {n}: {e}")

    rng = np.random.default_rng(0)
    rej = []
    for k in range(trials):
        s = v.CommonResidualSearchV2(synthetic_h0_array(rng, pulsars),
                                     seed=k, verbose=False, n_null=200, n_scramble=200).run()
        rej.append(np.array(list(s.p_values) + list(s.exploratory_p)) <= 0.05)
        if (k + 1) % 20 == 0:
            print(f"    {k + 1}/{trials} trials", flush=True)
    names = list(s.test_names) + list(s.exploratory_names)
    rates = dict(zip(names, np.mean(rej, axis=0).round(3).tolist()))
    passed = rates[PRIMARY] <= GATE_THRESHOLD
    print("\n  H_0 rates p≤0.05:")
    for n, r in rates.items():
        print(f"    {n:<32s} {r:.3f}{'   ← PRIMARY' if n == PRIMARY else ''}")
    print(f"\n  GATE {'PASSED' if passed else 'FAILED'}: primary H_0 rate {rates[PRIMARY]:.3f} "
          f"(threshold {GATE_THRESHOLD}, pre-registered)")
    Path(GATE_FILE.format(variant=variant)).write_text(json.dumps({
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "variant": variant, "trials": trials, "rates": rates, "primary": PRIMARY,
        "threshold": GATE_THRESHOLD, "passed": bool(passed), "pulsars": sorted(pulsars),
        "excluded": failed, "longdouble_eps": eps, "frozen_hash_ok": ok,
        "clock_override": clk}, indent=2))


def stage_real(variant):
    import common_residual_search_v2 as v
    gate = Path(GATE_FILE.format(variant=variant))
    if not gate.exists() or not json.loads(gate.read_text())["passed"]:
        sys.exit(f"Gate for {variant} missing or FAILED — real-data run not allowed "
                 f"(pre-registration, Section 4).")
    eps = check_precision()
    ok, h = frozen_hash_ok()
    if not ok:
        sys.exit(f"Frozen V2.1 file changed (hash {h}) — aborting.")
    setup_pint()
    print(f"EPTA {variant} — REAL residuals, frozen V2.1\n")
    pulsars, failed = load_variant(variant, with_residuals=True)
    pulsars = {n: p for n, p in pulsars.items() if np.isfinite(p.ra) and np.isfinite(p.dec)}
    print(f"\nPer-pulsar QC (residual RMS after epoch averaging):")
    for n, p in sorted(pulsars.items()):
        print(f"  {n}: {len(p.toas):5d} epochs, {p.time_span:5.1f} yr, "
              f"RMS {np.std(p.residuals) * 1e6:8.3f} µs")
    for n, e in failed.items():
        print(f"  EXCLUDED {n}: {e}")
    print()
    s = v.CommonResidualSearchV2(pulsars).run()        # frozen defaults
    i = s.exploratory_names.index(PRIMARY)
    p = s.exploratory_p[i]
    print(f"\n=== PRE-REGISTERED PRIMARY RESULT ({variant}) ===")
    print(f"  {PRIMARY}: p = {p:.4g}  →  "
          f"{'p ≤ 0.05' if p <= 0.05 else 'NOT replicated at α = 0.05'}")
    s.save_json(f"epta_{variant}_results.json",
                extra={"dataset": f"EPTA DR2 {variant}, doi:10.5281/zenodo.8300645",
                       "primary_test": PRIMARY, "primary_p": p,
                       "excluded_pulsars": failed, "longdouble_eps": eps,
                       "gate_file": str(gate)})


if __name__ == "__main__":
    args = sys.argv[1:]
    variant = args[args.index("--variant") + 1] if "--variant" in args else "DR2full"
    if not args or args[0] not in ("check", "gate", "real"):
        sys.exit(__doc__)
    if args[0] == "check":
        stage_check(args[1] if len(args) > 1 and not args[1].startswith("--") else "J1909-3744")
    elif args[0] == "gate":
        trials = int(args[args.index("--trials") + 1]) if "--trials" in args else 200
        stage_gate(variant, trials)
    else:
        stage_real(variant)
