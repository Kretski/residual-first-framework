"""
Common Residual Search V2 — no interpolation, no zero-filling.

Pipeline
  1. Each pulsar's epoch-averaged residuals are BINNED (weighted mean per
     bin, default 30 d). Empty bins stay NaN — nothing is interpolated.
  1b. A quadratic is removed from each binned series over its own span.
  2. Pairwise correlation ρ_ij is computed ONLY over bins where both pulsars
     have data (pairs with < min_overlap common bins are dropped).
     Fisher z_ij = atanh(ρ_ij)·sqrt(n_ij − 3) puts pairs with different
     overlaps on a common scale.
  3. Null for "is there any cross-pulsar correlation?": each pulsar's
     series is MIRROR-SHIFTED within its own observed span (random circular
     shift of the even extension [x, reversed x], first L samples kept),
     then the same quadratic as in step 1b is re-fitted and removed.
     This keeps (a) every pulsar's own spectrum incl. intrinsic red noise,
     (b) the start/end of every data set (the dominant missing-data pattern)
     and (c) the effect of the spin-down fit.
     Calibration on IDEALISED synthetic arrays (independent random TOAs per
     pulsar, constant per-pulsar errors, red noise, quadratic removal;
     200 trials, H_0), FPR at α = 0.05:
        plain circular shift ........ 0.14   (anti-conservative: wrap-around jump)
        phase-randomised surrogate .. 0.00   (over-conservative)
        mirror shift (default) ...... 0.065  (consistent with 0.05)
     WARNING — these numbers do NOT carry over to real PTA data. With the real
     NANOGrav 15-yr sampling (--validate --real-gaps) the omnibus z² test
     (P1 cross_correlation) has FPR ≈ 0.67–0.70. Diagnosis
     (diagnostics_shared_structure.py): TOA uncertainties that change at the
     same dates for all pulsars (backend/receiver upgrades) make the binned
     series share a time-varying variance; independent per-pulsar shifts
     destroy that alignment, so the null z² is too small (FPR → 1.00 in a
     synthetic test with a common 5× error step; the same happens with TOA
     density that changes at the same dates — FPR 1.00 → 0.10 after thinning
     to one TOA per bin). Shared data GAPS alone are harmless (FPR 0.07).
     A second, non-synchronous mechanism: when intrinsic red noise dominates
     the white noise, the mirror-shift null is itself inaccurate (FPR ≈ 0.28
     with no shared structure at all). The idealised calibration above was in
     a white-dominated regime. Real-gaps ablation chain: 0.69 → 0.43 (constant
     errors) → 0.31 (+ one TOA per bin) ≈ red-dominance level.
     Signed projections and all sky-scramble tests remain calibrated on real
     sampling. P1 must therefore NOT be interpreted on real data.
     A plain time permutation would destroy red noise and is not a valid null.
  4. Angular structure: overlap-weighted least-squares projections of z_ij
     on monopole (clock), dipole (ephemeris) and Hellings–Downs (GWB).
     Two nulls: in-span shifts (as above), and SKY SCRAMBLES — the data
     are kept exactly and only pulsar positions are permuted. The sky
     scramble asks "does the correlation depend on angle as HD predicts",
     given whatever correlation is present.
  5. Benjamini–Hochberg over the registered tests.

  python common_residual_search_v2.py --real
  python common_residual_search_v2.py --real --min-span 10 --bin-days 30
  python common_residual_search_v2.py --validate          # synthetic, realistic gaps

This is NOT the NANOGrav optimal statistic or Bayesian analysis: it uses
post-fit residuals without the full noise covariance, so it is expected to
be less sensitive, and quadratic spin-down / astrometric fits are already
imprinted in the residuals.

Author: Dimitar Kretski
License: MIT
"""

import sys
import numpy as np
from typing import Dict, List, Optional, Tuple

from common_residual_search import benjamini_hochberg, hellings_downs, radec_to_unit


def permutation_p_value(t_obs, t_null) -> float:
    t_null = np.asarray(t_null)
    return (1.0 + np.sum(t_null >= t_obs)) / (1.0 + t_null.size)


# ============================================================================
# Binning (no interpolation)
# ============================================================================

def bin_residuals(toas, residuals, errors, edges, return_sigma: bool = False):
    """Inverse-variance weighted mean of residuals in each bin; NaN if empty.
    With return_sigma=True also returns the formal bin uncertainty 1/sqrt(Σw)."""
    idx = np.digitize(toas, edges) - 1
    ok = (idx >= 0) & (idx < len(edges) - 1)
    w = 1.0 / errors[ok] ** 2
    num = np.bincount(idx[ok], weights=w * residuals[ok], minlength=len(edges) - 1)
    den = np.bincount(idx[ok], weights=w, minlength=len(edges) - 1)
    out = np.full(len(edges) - 1, np.nan)
    has = den > 0
    out[has] = num[has] / den[has]
    if return_sigma:
        sig = np.full(len(edges) - 1, np.nan)
        sig[has] = 1.0 / np.sqrt(den[has])
        return out, sig
    return out


def build_binned_matrix(pulsars: Dict, bin_days: float = 30.0, return_sigma: bool = False):
    names = sorted(pulsars)
    t0 = min(p.toas.min() for p in pulsars.values())
    t1 = max(p.toas.max() for p in pulsars.values())
    edges = np.arange(t0, t1 + bin_days, bin_days)
    rows = [bin_residuals(pulsars[n].toas, pulsars[n].residuals, pulsars[n].errors, edges, True)
            for n in names]
    X = np.array([r[0] for r in rows])
    S = np.array([r[1] for r in rows])
    return (X, names, edges, S) if return_sigma else (X, names, edges)


# ============================================================================
# Pairwise statistics with missing data
# ============================================================================

def pairwise_correlations(X: np.ndarray, min_overlap: int = 12):
    """
    ρ_ij and n_ij over common non-NaN bins, fully vectorised
    (per-pair means computed on the overlap only).
    """
    M = np.isfinite(X).astype(float)
    Z = np.where(M > 0, X, 0.0)
    n = M @ M.T
    sx = Z @ M.T                     # Σ x_i over overlap with j
    sy = M @ Z.T                     # Σ x_j over overlap with i
    sxx = (Z ** 2) @ M.T
    syy = M @ (Z ** 2).T
    sxy = Z @ Z.T
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = sxy - sx * sy / n
        vx = sxx - sx ** 2 / n
        vy = syy - sy ** 2 / n
        rho = cov / np.sqrt(vx * vy)
    iu = np.triu_indices(X.shape[0], k=1)
    r, nn = rho[iu], n[iu]
    keep = (nn >= min_overlap) & np.isfinite(r)
    return r, nn, iu, keep


def fisher_z(r, n):
    r = np.clip(r, -0.999999, 0.999999)
    return np.arctanh(r) * np.sqrt(np.maximum(n - 3, 1))


def spans(X: np.ndarray) -> List[Tuple[int, int]]:
    out = []
    for row in X:
        ok = np.nonzero(np.isfinite(row))[0]
        out.append((ok[0], ok[-1] + 1) if len(ok) else (0, 0))
    return out


def detrend_rows(X: np.ndarray, deg: int = 2) -> np.ndarray:
    """
    Remove a polynomial (default quadratic) from each row over its occupied bins.
    Applied identically to the data and to every null realisation, so both have
    been through the same 'spin-down fit' (the timing model removes one per pulsar).
    """
    Y = X.copy()
    for i in range(X.shape[0]):
        ok = np.isfinite(X[i])
        if ok.sum() > deg + 2:
            t = np.nonzero(ok)[0].astype(float)
            t = (t - t.mean()) / (np.ptp(t) + 1.0)
            V = np.vander(t, deg + 1)
            Y[i, ok] = X[i, ok] - V @ np.linalg.lstsq(V, X[i, ok], rcond=None)[0]
    return Y


def phase_surrogate_within_span(X: np.ndarray, sp, rng) -> np.ndarray:
    """
    Fourier phase randomisation inside each pulsar's span (keeps the
    segment's power spectrum, no wrap-around jump); internal gaps are
    filled only to generate the surrogate and then re-imposed.
    """
    Y = np.full_like(X, np.nan)
    for i, (a, b) in enumerate(sp):
        seg = X[i, a:b]
        ok = np.isfinite(seg)
        if ok.sum() < 4:
            continue
        idx = np.arange(b - a)
        filled = np.interp(idx, idx[ok], seg[ok])
        F = np.fft.rfft(filled - filled.mean())
        ph = np.exp(1j * rng.uniform(0, 2 * np.pi, len(F)))
        ph[0] = 1.0
        if (b - a) % 2 == 0:
            ph[-1] = 1.0
        s = np.fft.irfft(F * ph, n=b - a)
        s[~ok] = np.nan
        Y[i, a:b] = s
    return Y


def mirror_shift_within_span(X: np.ndarray, sp, rng) -> np.ndarray:
    """
    Circular shift of the even (mirror) extension [seg, seg[::-1]] inside each
    span, keeping the first L samples. The extension is continuous at the wrap
    point, so — unlike a plain circular shift — no artificial jump (extra
    high-frequency power) is introduced into red-noise series. Random reversal
    of direction is included implicitly.
    """
    Y = np.full_like(X, np.nan)
    for i, (a, b) in enumerate(sp):
        L = b - a
        if L > 0:
            ext = np.concatenate([X[i, a:b], X[i, a:b][::-1]])
            Y[i, a:b] = np.roll(ext, rng.integers(0, 2 * L))[:L]
    return Y


def shift_within_span(X: np.ndarray, sp, rng) -> np.ndarray:
    """Circular shift of each row inside [first, last] observed bin; outside stays NaN."""
    Y = np.full_like(X, np.nan)
    for i, (a, b) in enumerate(sp):
        if b - a > 0:
            Y[i, a:b] = np.roll(X[i, a:b], rng.integers(0, b - a))
    return Y


# ============================================================================
# The search
# ============================================================================

def null_summary(obs: float, null, kind: str) -> Dict:
    """Everything needed to reproduce / audit a permutation p-value."""
    null = np.asarray(null, dtype=float)
    exceed = int(np.sum(null >= obs))
    return {'observed': float(obs), 'null_type': kind, 'n_null': int(null.size),
            'null_mean': float(null.mean()), 'null_std': float(null.std()),
            'null_quantiles': {q: float(np.quantile(null, float(q))) for q in
                               ('0.05', '0.5', '0.95', '0.99')},
            'n_exceed': exceed, 'p': (exceed + 1) / (null.size + 1),
            'p_formula': '(n_exceed + 1) / (n_null + 1), one-sided, null >= observed'}


class CommonResidualSearchV2:
    def __init__(self, pulsars: Dict, bin_days: float = 30.0, min_overlap: int = 12,
                 n_null: int = 2000, n_scramble: int = 2000, seed: Optional[int] = 0,
                 verbose: bool = True, null: str = "mirror", detrend_deg: int = 2,
                 normalize_errors: bool = False):
        """
        null: 'mirror' (default, calibrated — see module docstring),
              'shift'  (plain circular shift within span; anti-conservative, FPR ≈ 0.14),
              'phase'  (phase-randomised surrogate; over-conservative, FPR ≈ 0.00).
        """
        self.pulsars = pulsars
        X, self.names, self.edges, S = build_binned_matrix(pulsars, bin_days, return_sigma=True)
        self.detrend_deg = detrend_deg
        self.X = detrend_rows(X, detrend_deg) if detrend_deg >= 0 else X
        # V3 candidate (NOT part of the frozen V2.1 analysis): divide each bin by its
        # formal uncertainty so that error changes shared by all pulsars (backend
        # upgrades) no longer give the series a common time-varying variance.
        self.normalize_errors = normalize_errors
        if normalize_errors:
            self.X = self.X / S
        self.null = null
        self.bin_days, self.min_overlap = bin_days, min_overlap
        self.n_null, self.n_scramble = n_null, n_scramble
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.verbose = verbose
        self.provenance: Dict[str, Dict] = {}
        self.spans = spans(self.X)

        coords = np.array([(pulsars[n].ra, pulsars[n].dec) for n in self.names])
        self.nhat = radec_to_unit(coords[:, 0], coords[:, 1])

        r, nn, self.iu, self.keep = pairwise_correlations(self.X, min_overlap)
        self.rho, self.n_overlap = r, nn
        self.z = fisher_z(r, nn)
        self.test_names, self.p_values, self.results = [], [], {}
        self._null_z = None

    # ---------------------------------------------------------------- helpers
    def _templates(self, nhat):
        cosz = (nhat @ nhat.T)[self.iu]
        return {'monopole': np.ones_like(cosz), 'dipole': cosz,
                'hellings_downs': hellings_downs(cosz)}

    def _projections(self, z, keep, nhat):
        """Weighted LS amplitude of z on each template (weights = n_ij − 3 via z-scaling)."""
        out = {}
        for k, g in self._templates(nhat).items():
            gk = g[keep]
            w = np.sqrt(np.maximum(self.n_overlap[keep] - 3, 1))   # z already ∝ sqrt(n-3)
            out[k] = float(np.sum(w * gk * z[keep]) / np.sqrt(np.sum((w * gk) ** 2)))
        return out

    def _nulls(self):
        """In-span-shift null realisations of (z, keep) — generated once, reused."""
        if self._null_z is None:
            self._null_z = []
            for _ in range(self.n_null):
                gen = {"phase": phase_surrogate_within_span, "shift": shift_within_span,
                       "mirror": mirror_shift_within_span}[self.null]
                Y = gen(self.X, self.spans, self.rng)
                if self.detrend_deg >= 0:
                    Y = detrend_rows(Y, self.detrend_deg)
                r, nn, _, keep = pairwise_correlations(Y, self.min_overlap)
                self._null_z.append((fisher_z(r, nn), keep))
        return self._null_z

    def _register(self, name, p, res):
        self.test_names.append(name)
        self.p_values.append(float(p))
        self.results[name] = res

    # ---------------------------------------------------------------- tests
    def describe(self):
        occ = np.isfinite(self.X)
        span_frac = np.mean([(b - a) / self.X.shape[1] for a, b in self.spans])
        within = np.mean([np.mean(occ[i, a:b]) for i, (a, b) in enumerate(self.spans) if b > a])
        if self.verbose:
            print(f"Pulsars: {len(self.names)}, bins: {self.X.shape[1]} × {self.bin_days:g} d")
            print(f"  Mean observed span: {span_frac:.1%} of the full grid "
                  f"(outside-span bins are NOT filled)")
            print(f"  Occupied bins inside spans: {within:.1%}")
            print(f"  Pairs with ≥ {self.min_overlap} common bins: {self.keep.sum()} / {len(self.keep)} "
                  f"(median overlap {np.median(self.n_overlap[self.keep]):.0f} bins)")

    def test_cross_correlation(self):
        """Omnibus: mean z² over usable pairs vs in-span-shift null."""
        t_obs = np.mean(self.z[self.keep] ** 2)
        null = np.array([np.mean(z[k] ** 2) for z, k in self._nulls()])
        p = permutation_p_value(t_obs, null)
        self.provenance['P1_cross_correlation'] = null_summary(t_obs, null, 'mirror-shift within span')
        res = {'mean_z2': t_obs, 'null_mean': null.mean(), 'null_99': np.quantile(null, 0.99),
               'mean_rho': float(np.mean(self.rho[self.keep])), 'p_value': p}
        self._register('cross_correlation', p, res)
        if self.verbose:
            print("\nCross-pulsar correlation (no interpolation):")
            print(f"  mean z² = {t_obs:.3f}  (null mean {null.mean():.3f}, 99%: {res['null_99']:.3f}); "
                  f"mean ρ = {res['mean_rho']:+.4f}")
            print(f"  p = {p:.4g}   (floor 1/{self.n_null + 1} = {1/(self.n_null + 1):.2g})")
        return res

    def test_angular_structure(self):
        obs = self._projections(self.z, self.keep, self.nhat)

        shift_null = [self._projections(z, k, self.nhat) for z, k in self._nulls()]
        p_shift = {k: permutation_p_value(obs[k], [d[k] for d in shift_null]) for k in obs}

        scr_null = []
        for _ in range(self.n_scramble):
            scr_null.append(self._projections(self.z, self.keep,
                                              self.nhat[self.rng.permutation(len(self.names))]))
        p_scr = {k: permutation_p_value(obs[k], [d[k] for d in scr_null])
                 for k in ('dipole', 'hellings_downs')}

        g = self._templates(self.nhat)
        overlap = float(np.dot(g['hellings_downs'][self.keep], g['monopole'][self.keep]) /
                        (np.linalg.norm(g['hellings_downs'][self.keep]) *
                         np.linalg.norm(g['monopole'][self.keep])))

        self.provenance['P2_monopole'] = null_summary(
            obs['monopole'], [d['monopole'] for d in shift_null], f'{self.null}-shift within span')
        self.provenance['P3_hd_sky_scramble'] = null_summary(
            obs['hellings_downs'], [d['hellings_downs'] for d in scr_null], 'sky-position scramble')
        self.provenance['info_hd_shift'] = null_summary(
            obs['hellings_downs'], [d['hellings_downs'] for d in shift_null], f'{self.null}-shift within span')
        self.provenance['info_dipole_trace_scramble'] = null_summary(
            obs['dipole'], [d['dipole'] for d in scr_null], 'sky-position scramble')
        res = {'projections': obs, 'p_shift': p_shift, 'p_sky_scramble': p_scr,
               'hd_monopole_overlap': overlap}
        self._register('monopole', p_shift['monopole'], res)
        self._register('hd_sky_scramble', p_scr['hellings_downs'], res)
        if self.verbose:
            print("\nAngular structure (one-sided, positive amplitude):")
            print(f"  {'template':<15s} {'proj':>8s} {'p(in-span shift)':>18s} {'p(sky scramble)':>17s}")
            for k in obs:
                ps = f"{p_scr[k]:.4g}" if k in p_scr else "—"
                print(f"  {k:<15s} {obs[k]:+8.3f} {p_shift[k]:18.4g} {ps:>17s}")
            print(f"  HD–monopole template overlap: {overlap:.2f}")
        return res

    def _pair_weights(self, keep):
        return np.sqrt(np.maximum(self.n_overlap[keep] - 3, 1))

    def _nuisance_basis(self, nhat, with_dipole: bool):
        """
        Pair-space nuisance patterns. Monopole: 1. Dipole along an arbitrary
        direction d gives ρ_ab ∝ (n_a·d)(n_b·d) = Σ_ij d_i d_j · ½(n_ai n_bj + n_aj n_bi),
        so the full dipole subspace is spanned by the 6 symmetric products
        (their trace is cos ζ).
        """
        cols = [np.ones(len(self.iu[0]))]
        if with_dipole:
            a, b = nhat[self.iu[0]], nhat[self.iu[1]]
            for i in range(3):
                for j in range(i, 3):
                    cols.append(0.5 * (a[:, i] * b[:, j] + a[:, j] * b[:, i]))
        return np.column_stack(cols)

    def _hd_perp_projection(self, z, keep, nhat, with_dipole):
        w = self._pair_weights(keep)
        g = w * hellings_downs((nhat @ nhat.T)[self.iu])[keep]
        B = w[:, None] * self._nuisance_basis(nhat, with_dipole)[keep]
        g_perp = g - B @ np.linalg.lstsq(B, g, rcond=None)[0]
        norm = np.linalg.norm(g_perp)
        frac = norm / np.linalg.norm(g)
        return float(np.dot(g_perp, z[keep]) / norm), frac

    def test_orthogonalized_hd(self):
        """
        EXPLORATORY (added after seeing V2 results): HD with the monopole
        (and monopole + full dipole) subspace projected out.
        Reported as a separate family with its own BH correction.
        """
        out = {}
        for label, with_d in (("hd_perp_mono", False), ("hd_perp_mono_dipole", True)):
            obs, frac = self._hd_perp_projection(self.z, self.keep, self.nhat, with_d)
            n_shift = [self._hd_perp_projection(z, k, self.nhat, with_d)[0] for z, k in self._nulls()]
            n_scr = [self._hd_perp_projection(self.z, self.keep,
                                              self.nhat[self.rng.permutation(len(self.names))],
                                              with_d)[0] for _ in range(self.n_scramble)]
            self.provenance[f'E_{label}_shift'] = null_summary(obs, n_shift, f'{self.null}-shift within span')
            self.provenance[f'E_{label}_scramble'] = null_summary(obs, n_scr, 'sky-position scramble')
            out[label] = {'proj': obs, 'hd_fraction_kept': frac,
                          'p_shift': permutation_p_value(obs, n_shift),
                          'p_scramble': permutation_p_value(obs, n_scr)}
        self.exploratory = out
        if self.verbose:
            print("\nEXPLORATORY — HD after removing nuisance subspaces (separate BH family):")
            print(f"  {'test':<22s} {'proj':>7s} {'HD kept':>8s} {'p(shift)':>10s} {'p(scramble)':>12s}")
            for k, r in out.items():
                print(f"  {k:<22s} {r['proj']:+7.3f} {r['hd_fraction_kept']:8.0%} "
                      f"{r['p_shift']:10.4g} {r['p_scramble']:12.4g}")
        return out

    def save_json(self, path: str, extra: Optional[Dict] = None) -> str:
        """Write a self-contained, reproducible record of this run."""
        import json, hashlib, platform, datetime
        from pathlib import Path
        script = Path(__file__).read_bytes()
        rec = {
            'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'),
            'script': Path(__file__).name,
            'script_sha256': hashlib.sha256(script).hexdigest(),
            'python': platform.python_version(), 'numpy': np.__version__,
            'parameters': {'bin_days': self.bin_days, 'min_overlap': self.min_overlap,
                           'n_null': self.n_null, 'n_scramble': self.n_scramble,
                           'temporal_null': self.null, 'detrend_deg': self.detrend_deg,
                           'seed': self.seed},
            'pulsars': self.names,
            'binned_data_sha256': hashlib.sha256(np.nan_to_num(self.X, nan=-9e99).tobytes()).hexdigest(),
            'n_pairs_used': int(self.keep.sum()),
            'primary_family': {
                'tests': self.test_names, 'p_raw': self.p_values,
                'p_bh': [float(x) for x in self.p_adj], 'reject': [bool(x) for x in self.reject],
                'status': 'pre-specified before the V2 real-data run (not formally preregistered)'},
            'exploratory_family': {
                'tests': getattr(self, 'exploratory_names', []),
                'p_raw': getattr(self, 'exploratory_p', []),
                'reject': [bool(x) for x in getattr(self, 'exploratory_reject', [])],
                'status': 'added after seeing V2 results — exploratory only'},
            'null_distributions': self.provenance,
            'hd_monopole_overlap': self.results.get('monopole', {}).get('hd_monopole_overlap'),
        }
        if extra:
            rec.update(extra)
        Path(path).write_text(json.dumps(rec, indent=2, ensure_ascii=False), encoding='utf-8')
        if self.verbose:
            print(f"\nProvenance record written to {path}")
        return path

    def binned_hd_curve(self, n_bins: int = 8):
        """Mean ρ in angular-separation bins (for plotting / eyeballing HD)."""
        ang = np.degrees(np.arccos(np.clip((self.nhat @ self.nhat.T)[self.iu], -1, 1)))
        a, r = ang[self.keep], self.rho[self.keep]
        e = np.linspace(0, 180, n_bins + 1)
        rows = []
        for lo, hi in zip(e[:-1], e[1:]):
            m = (a >= lo) & (a < hi)
            if m.sum() > 2:
                rows.append(((lo + hi) / 2, r[m].mean(), r[m].std() / np.sqrt(m.sum()), int(m.sum())))
        if self.verbose:
            print("\nMean ρ vs angular separation (compare shape with HD):")
            print(f"  {'ζ [deg]':>8s} {'mean ρ':>9s} {'± se':>8s} {'pairs':>6s} {'HD(ζ)':>8s}")
            for c, m, s, k in rows:
                print(f"  {c:8.1f} {m:+9.4f} {s:8.4f} {k:6d} {hellings_downs(np.cos(np.radians(c))):+8.3f}")
        return rows

    def run(self, exploratory: bool = True):
        self.describe()
        self.test_cross_correlation()
        self.test_angular_structure()
        self.binned_hd_curve()
        p_adj, rej = benjamini_hochberg(self.p_values, 0.05)
        self.p_adj, self.reject = p_adj, rej
        if self.verbose:
            print("\nPRIMARY family (pre-specified in V2), Benjamini–Hochberg α = 0.05:")
            for n, p, q, r in zip(self.test_names, self.p_values, p_adj, rej):
                print(f"  {n:<18s} p_raw={p:.4g}  p_BH={q:.4g}  {'← significant' if r else ''}")
            print("  (monopole projection = overlap-weighted SIGNED mean of z_ij)")

        if exploratory:
            ex = self.test_orthogonalized_hd()
            names = [f"{k}:{t}" for k in ex for t in ("shift", "scramble")]
            ps = [ex[k][f"p_{t}"] for k in ex for t in ("shift", "scramble")]
            q, r = benjamini_hochberg(ps, 0.05)
            self.exploratory_names, self.exploratory_p, self.exploratory_reject = names, ps, r
            if self.verbose:
                print("  Exploratory BH:")
                for n, p, qq, rr in zip(names, ps, q, r):
                    print(f"    {n:<30s} p_raw={p:.4g}  p_BH={qq:.4g}  {'← significant' if rr else ''}")
        return self


# ============================================================================
# Synthetic validation with realistic spans/gaps
# ============================================================================

class _P:  # minimal stand-in for PulsarTimingData
    def __init__(self, toas, residuals, errors, ra, dec):
        self.toas, self.residuals, self.errors, self.ra, self.dec = toas, residuals, errors, ra, dec


def _red(rng, t_days, amp, beta, grid):
    f = np.fft.rfftfreq(len(grid))
    a = np.zeros_like(f); a[1:] = f[1:] ** (-beta / 2)
    x = np.fft.irfft(a * (rng.standard_normal(len(f)) + 1j * rng.standard_normal(len(f))), n=len(grid))
    return np.interp(t_days, grid, amp * x / x.std())


def synthetic_array(rng, template: Optional[Dict] = None, n_psr=68, hd_amp=0.0, mono_amp=0.0,
                    dipole_amp=0.0, flat_errors: bool = False, one_per_bin: float = 0.0):
    """
    Independent white + intrinsic red noise (1/3 of pulsars), optional HD-correlated
    or monopole common red process; quadratic removed per pulsar (mimics spin-down fit).
    If `template` (real pulsars) is given, their TOAs, errors and positions are reused.
    """
    T = 16 * 365.25
    grid = np.linspace(0, T, 4096)
    if template:
        names = sorted(template)
        toas = [template[n].toas - min(p.toas.min() for p in template.values()) for n in names]
        errs = [template[n].errors for n in names]
        if one_per_bin > 0:  # ablation: at most one TOA per analysis bin (removes TOA-density structure)
            thinned_t, thinned_e = [], []
            for t, e in zip(toas, errs):
                _, first = np.unique(np.floor(t / one_per_bin), return_index=True)
                thinned_t.append(t[first]); thinned_e.append(e[first])
            toas, errs = thinned_t, thinned_e
        if flat_errors:     # ablation: keep sampling & gaps, remove any time variation of errors
            errs = [np.full(len(e), np.median(e)) for e in errs]
        ra = np.array([template[n].ra for n in names]); dec = np.array([template[n].dec for n in names])
    else:
        names = [f"P{i:02d}" for i in range(n_psr)]
        toas, errs = [], []
        for _ in names:
            start = rng.uniform(0, 13) * 365.25
            n = int((T - start) / 25)
            toas.append(np.sort(rng.uniform(start, T, n)))
            errs.append(np.full(n, rng.uniform(0.3, 2.0) * 1e-6))
        ra = rng.uniform(0, 360, n_psr); dec = np.degrees(np.arcsin(rng.uniform(-0.6, 1, n_psr)))

    common = np.zeros((len(names), len(grid)))
    if hd_amp > 0:
        nh = radec_to_unit(ra, dec)
        G = hellings_downs(nh @ nh.T); np.fill_diagonal(G, 1.0)
        L = np.linalg.cholesky(G + 1e-9 * np.eye(len(names)))
        base = np.array([_red(rng, grid, 1.0, 13 / 3, grid) for _ in names])
        common += hd_amp * (L @ base)
    if mono_amp > 0:
        common += mono_amp * _red(rng, grid, 1.0, 13 / 3, grid)
    if dipole_amp > 0:                                     # ephemeris-like: s(t) · (n̂·d)
        d = rng.standard_normal(3); d /= np.linalg.norm(d)
        common += dipole_amp * np.outer(radec_to_unit(ra, dec) @ d, _red(rng, grid, 1.0, 13 / 3, grid))

    out = {}
    for i, n in enumerate(names):
        t, e = toas[i], errs[i]
        r = e * rng.standard_normal(len(t))
        if i % 3 == 0:
            r += _red(rng, t, 1.0e-6, 3.0, grid)
        r += np.interp(t, grid, common[i])
        V = np.vander((t - t.mean()) / T, 3)                    # remove quadratic (spin-down fit)
        r -= V @ np.linalg.lstsq(V / e[:, None], r / e, rcond=None)[0]
        out[n] = _P(t, r, e, ra[i], dec[i])
    return out


def validate(template=None, n_trials=100, seed=0, null="mirror", n_null=200,
             flat_errors=False, normalize_errors=False, h0_only=False, one_per_bin=0.0):
    rng = np.random.default_rng(seed)
    for label, kw, nt in [("H_0 (no common signal)", {}, n_trials),
                          ("HD common process 1.0 µs", {'hd_amp': 1.0e-6}, 20),
                          ("monopole 0.5 µs", {'mono_amp': 0.5e-6}, 20),
                          ("dipole 1.0 µs", {'dipole_amp': 1.0e-6}, 20)][:1 if h0_only else 4]:
        rej = []
        for k in range(nt):
            ps = synthetic_array(rng, template, flat_errors=flat_errors, one_per_bin=one_per_bin, **kw)
            s = CommonResidualSearchV2(ps, n_null=n_null, n_scramble=n_null, seed=k,
                                       verbose=False, null=null,
                                       normalize_errors=normalize_errors).run()
            rej.append(np.array(list(s.p_values) + list(s.exploratory_p)) <= 0.05)
        rej = np.array(rej)
        allnames = list(s.test_names) + list(s.exploratory_names)
        print(f"{label}:  rate p≤0.05 →\n    " +
              "\n    ".join(f"{n:<32s} {r:.2f}" for n, r in zip(allnames, rej.mean(axis=0))))


# ============================================================================

if __name__ == "__main__":
    args = sys.argv[1:]

    def opt(name, default, cast=float):
        return cast(args[args.index(name) + 1]) if name in args else default

    min_span = opt("--min-span", 0.0)
    bin_days = opt("--bin-days", 30.0)
    n_null = opt("--n-null", 2000, int)

    if "--real" in args or ("--validate" in args and "--real-gaps" in args):
        import io, contextlib
        from nanograv_loader import NANOGravDR15Loader
        with contextlib.redirect_stdout(io.StringIO()):          # silence per-file column messages
            pulsars = NANOGravDR15Loader(use_synthetic=False).load_all_pulsars()
        pulsars = {n: p for n, p in pulsars.items()
                   if np.isfinite(p.ra) and np.isfinite(p.dec) and p.time_span >= min_span}
    else:
        pulsars = None

    if "--validate" in args:
        print("Synthetic validation" + (" (REAL sampling, errors, positions)" if pulsars else "") + "\n")
        flat, norm = "--flat-errors" in args, "--normalize-errors" in args
        if flat:
            print("ABLATION: per-pulsar constant errors (real sampling, gaps, positions kept)")
        if norm:
            print("V3 CANDIDATE: residuals divided by bin uncertainty")
        opb = bin_days if "--one-per-bin" in args else 0.0
        if opb:
            print(f"ABLATION: at most one TOA per {opb:g}-d bin (TOA-density structure removed)")
        validate(template=pulsars, n_trials=int(opt("--trials", 100, int)),
                 flat_errors=flat, normalize_errors=norm, h0_only="--h0-only" in args,
                 one_per_bin=opb)
        print("\nExpected: H_0 rates ≈ 0.05 (±0.04 with 100 trials); "
              "HD signal → HD tests high; monopole/dipole → hd_perp_mono_dipole ≈ 0.05.")
    elif pulsars:
        print(f"Common Residual Search V2 — REAL NANOGrav 15-yr, "
              f"{len(pulsars)} pulsars with span ≥ {min_span:g} yr\n")
        seed = opt("--seed", 0, int)
        search = CommonResidualSearchV2(pulsars, bin_days=bin_days, n_null=n_null, n_scramble=n_null,
                                        null=opt("--null", "mirror", str), seed=seed).run()
        tag = f"span{min_span:g}_bin{bin_days:g}_seed{seed}"
        family = "PRIMARY" if (min_span == 0 and bin_days == 30) else "EXPLORATORY (sensitivity cut)"
        print(f"\nThis configuration belongs to the {family} analysis.")
        search.save_json(opt("--out", f"v2_results_{tag}.json", str),
                         extra={'configuration_family': family, 'min_span_yr': min_span,
                                'data': 'NANOGrav 15-yr v2.1.0 narrowband epoch-averaged residuals, '
                                        'doi:10.5281/zenodo.16051178 (archive MD5 '
                                        '557d42dd8486a5f8272d90dec9b228a8)'})
    else:
        print("Usage: python common_residual_search_v2.py --real [--min-span 10] [--bin-days 30]\n"
              "       python common_residual_search_v2.py --validate [--real-gaps] [--trials 100]")
