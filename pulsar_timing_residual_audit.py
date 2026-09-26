"""
Pulsar Timing Residual Audit — Residual-First Framework

Stages 0–3 adapted for pulsar timing:
  Stage 0: Load post-fit residuals
  Stage 1: Characterize residual (red/white noise, spectral lines, DM, χ²)
  Stage 2: Benjamini–Hochberg multiple-testing correction
  Stage 3: Replication / common signal across pulsars

Changes vs. previous version (why the numbers changed):
  * PSD: pulsar TOAs are UNEVENLY sampled, so Welch (which assumes uniform
    sampling) is replaced by a Lomb–Scargle periodogram on the grid f_k = k/T.
  * Red/white decomposition: the old log-log line fit had the sign of β
    inverted (slope = -β for S ∝ f^-β), had no white-noise term and evaluated
    the model at f = 0 → divide-by-zero / NaN → "Flagged 0/57, p = 1".
    Now: Whittle-likelihood fit of S(f) = A (f/f_yr)^-β + C.
  * Stage 1.1: max-ratio test P_k / S_model,k, calibrated by parametric
    bootstrap from the fitted red+white model (not a hard-coded σ = 1/√2).
  * Stage 1.2: weighted regression on ν^-2 (DM delay law), t-test on slope;
    skipped automatically for single-band data.
  * Stage 1.3: redness = likelihood-ratio (red+white vs white), calibrated by
    permuting the normalized residuals r_i/σ_i (no hard-coded se_β = 0.4).
  * Stage 1.4: χ² with dof = N − n_fit_params, two-sided.
  * BH-adjusted p-values fixed (and kept in original test order);
    summary() no longer crashes on a numpy array.
  * Stage 3: PC1 variance tested against random circular shifts of each
    pulsar (keeps each pulsar's own autocorrelation, destroys cross-pulsar
    correlation) instead of an ad-hoc z-score.

Reference: NANOGrav Collaboration (2023), "The NANOGrav 15 yr Data Set", ApJL.

Author: Dimitar Kretski
License: MIT
"""

import numpy as np
from scipy import stats, optimize
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field

SEC_PER_DAY = 86400.0
F_YR = 1.0 / (365.25 * SEC_PER_DAY)          # 1/yr in Hz


# ============================================================================
# Helpers
# ============================================================================

def benjamini_hochberg(p_values, alpha: float = 0.05):
    """BH step-up. Returns (adjusted p, reject flags), both in input order."""
    p = np.asarray(p_values, dtype=float)
    m = p.size
    order = np.argsort(p)
    ranks = np.arange(1, m + 1)
    q = np.minimum.accumulate((p[order] * m / ranks)[::-1])[::-1]
    p_adj = np.empty(m)
    p_adj[order] = np.clip(q, 0, 1)
    below = p[order] <= ranks / m * alpha
    reject = np.zeros(m, dtype=bool)
    if below.any():
        reject[order[:np.max(np.nonzero(below)[0]) + 1]] = True
    return p_adj, reject


def permutation_p_value(t_obs, t_null) -> float:
    t_null = np.asarray(t_null)
    return (1.0 + np.sum(t_null >= t_obs)) / (1.0 + t_null.size)


class LombScargleBasis:
    """
    Classic Lomb–Scargle periodogram as a one-sided PSD (s^2/Hz), with the
    trigonometric matrices cached for a fixed (t, f) — re-evaluating for new y
    (permutations, simulations) is then just two matrix-vector products.
    For white noise of variance σ²: E[S] = 2σ²T/N and S/E[S] ~ Exp(1).
    Written out explicitly so it works on any SciPy version.
    """

    def __init__(self, t_sec: np.ndarray, freqs_hz: np.ndarray):
        w = 2 * np.pi * freqs_hz[:, None]
        wt = w * t_sec[None, :]
        tau = np.arctan2(np.sum(np.sin(2 * wt), axis=1),
                         np.sum(np.cos(2 * wt), axis=1)) / (2 * 2 * np.pi * freqs_hz)
        arg = w * (t_sec[None, :] - tau[:, None])
        self.c, self.s = np.cos(arg), np.sin(arg)
        self.cc = np.sum(self.c ** 2, axis=1)
        self.ss = np.sum(self.s ** 2, axis=1)
        self.norm = 2.0 * (t_sec[-1] - t_sec[0]) / len(t_sec)

    def psd(self, y: np.ndarray) -> np.ndarray:
        y = y - np.mean(y)
        return self.norm * 0.5 * ((self.c @ y) ** 2 / self.cc + (self.s @ y) ** 2 / self.ss)


def lomb_scargle_psd(t_sec: np.ndarray, y: np.ndarray, freqs_hz: np.ndarray) -> np.ndarray:
    """One-off convenience wrapper around LombScargleBasis."""
    return LombScargleBasis(t_sec, freqs_hz).psd(y)


def _model_psd(theta, f):
    log_a, beta, log_c = theta
    return 10 ** log_a * (f / F_YR) ** (-beta) + 10 ** log_c


def fit_red_white_whittle(f: np.ndarray, S: np.ndarray) -> Dict:
    """
    Whittle fit of S(f) = A (f/f_yr)^-β + C.
    Returns best-fit params, negative log-likelihoods, and the LR statistic
    2·(nll_white − nll_red+white) ≥ 0.
    """
    c_white = np.mean(S)
    nll_white = np.sum(np.log(c_white) + S / c_white)

    lnf = np.log(f / F_YR)
    LN10 = np.log(10.0)

    def nll_and_grad(theta):
        log_a, beta, log_c = theta
        g = np.exp(-beta * lnf)
        A, C = 10 ** log_a, 10 ** log_c
        m = A * g + C
        val = np.sum(np.log(m) + S / m)
        d = (m - S) / m ** 2
        grad = np.array([np.sum(d * LN10 * A * g),
                         np.sum(d * (-A * g * lnf)),
                         np.sum(d * LN10 * C)])
        return val, grad

    lc0 = np.log10(c_white)
    la0 = np.log10(np.mean(S[:3]))
    bounds = [(lc0 - 12, lc0 + 12), (0.0, 8.0), (lc0 - 6, lc0 + 2)]
    best = None
    for beta0 in (1.0, 3.0, 5.0):
        res = optimize.minimize(nll_and_grad, [la0 - 1.0, beta0, lc0 - 0.3], jac=True,
                                method="L-BFGS-B", bounds=bounds)
        if best is None or res.fun < best.fun:
            best = res
    # White-only is nested (A → 0); never report LR < 0 because of optimizer noise
    nll_rw = min(best.fun, nll_white)
    log_a, beta, log_c = best.x
    beta_lo, beta_hi = bounds[1]
    at_bound = bool(min(beta - beta_lo, beta_hi - beta) < 0.02)
    return {
        'beta_at_bound': at_bound,
        'red_amplitude': 10 ** log_a, 'beta': beta, 'white_level': 10 ** log_c,
        'theta': best.x, 'nll_red_white': nll_rw, 'nll_white': nll_white,
        'white_only_level': c_white, 'lr_stat': 2 * (nll_white - nll_rw),
    }


# ============================================================================
# Result container (unchanged interface)
# ============================================================================

@dataclass
class TimingResidualAuditResult:
    pulsar_name: str
    n_observations: int
    time_span: float  # years
    p_values_raw: List[float] = field(default_factory=list)
    p_values_corrected: List[float] = field(default_factory=list)
    test_names: List[str] = field(default_factory=list)
    residual_rms: float = 0.0
    red_noise_amplitude: float = 0.0
    red_noise_spectral_index: float = 0.0
    white_noise_level: float = 0.0
    anomalies: List[Dict] = field(default_factory=list)
    h0_rejected: bool = False


# ============================================================================
# Single-pulsar audit
# ============================================================================

class PulsarTimingResidualAudit:
    """
    Null hypothesis: post-fit residuals = white measurement noise
    + intrinsic power-law red noise, no deterministic/common signal.
    """

    def __init__(self, pulsar_name: str, verbose: bool = False,
                 seed: Optional[int] = None, rng: Optional[np.random.Generator] = None):
        self.name = pulsar_name
        self.verbose = verbose
        self.rng = rng if rng is not None else np.random.default_rng(seed)

        self.toas = self.residuals = self.errors = self.frequencies = None
        self.dm = self.ra = self.dec = None
        self._t_sec = None

        self.psd_freq = self.psd_data = None
        self.psd_model_white = self.psd_model_red = None
        self.red_noise_amplitude = self.red_noise_spectral_index = None
        self.white_noise_level = None
        self.fit = None

        self.test_results: Dict[str, Dict] = {}
        self.test_names: List[str] = []
        self.p_values_raw: List[float] = []
        self.p_values_corrected = None
        self.reject_flags = None
        self.alpha_global = None
        self.anomalies: List[Dict] = []

    def _register(self, name, result):
        self.test_results[name] = result
        self.test_names.append(name)
        self.p_values_raw.append(float(result['p_value']))

    # ------------------------------------------------------------------ Stage 0
    def stage0_load_pulsar_data(self, toas, residuals, errors, frequencies, dm,
                                ra: float = None, dec: float = None):
        """toas [MJD], residuals/errors [s], frequencies [MHz], dm [pc/cm³]."""
        toas, residuals, errors = map(np.asarray, (toas, residuals, errors))
        if not (len(toas) == len(residuals) == len(errors)):
            raise ValueError("Mismatched array lengths")
        order = np.argsort(toas)
        self.toas, self.residuals, self.errors = toas[order], residuals[order], errors[order]
        self.frequencies = None if frequencies is None else np.asarray(frequencies)[order]
        self.dm, self.ra, self.dec = dm, ra, dec
        self._t_sec = (self.toas - self.toas[0]) * SEC_PER_DAY

        if self.verbose:
            print(f"[{self.name}] Stage 0: Data loaded")
            print(f"  Observations: {len(self.toas)}")
            print(f"  Time span: {self.toas[-1] - self.toas[0]:.1f} days")
            print(f"  Residual RMS: {np.std(self.residuals)*1e6:.3f} µs")

    # ------------------------------------------------------------------ Stage 1
    def _frequency_grid(self, n_freq: Optional[int] = None):
        """k/T, k = 1..n. Default n = (number of distinct observing days)/2 —
        many TOAs on the same day (multi-channel data) add no higher frequencies."""
        T = self._t_sec[-1]
        if n_freq is None:
            n_epochs = len(np.unique(np.round(self.toas)))
            n_freq = max(8, min(len(self.residuals), n_epochs) // 2)
        return np.arange(1, n_freq + 1) / T

    def stage1_psd_red_white_decomposition(self, n_freq: Optional[int] = None,
                                           oversample: int = 4
                                           ) -> Tuple[np.ndarray, np.ndarray, Dict]:
        """Lomb–Scargle PSD + Whittle fit of S = A (f/f_yr)^-β + C."""
        if self.residuals is None:
            raise RuntimeError("Call stage0_load_pulsar_data first")
        self._n_freq = n_freq
        self.psd_freq = self._frequency_grid(n_freq)
        self._ls = LombScargleBasis(self._t_sec, self.psd_freq)
        self._sim_phase = 2 * np.pi * np.outer(self._t_sec, self.psd_freq)
        self._sim_cos, self._sim_sin = np.cos(self._sim_phase), np.sin(self._sim_phase)
        self.psd_data = self._ls.psd(self.residuals)
        # Oversampled grid, used only by the spectral-line search (lines between
        # the k/T bins lose up to ~60% of their power on the coarse grid).
        # The bootstrap in stage 1.1 uses the same grid, so calibration is exact.
        max_elems = 3e7                       # ~250 MB per trig matrix
        while oversample > 1 and len(self._t_sec) * len(self.psd_freq) * oversample > max_elems:
            oversample -= 1
        if len(self._t_sec) * len(self.psd_freq) > max_elems:
            raise MemoryError(f"{len(self._t_sec)} TOAs × {len(self.psd_freq)} frequencies is too "
                              f"large — use epoch-averaged residuals or pass n_freq=")
        df = self.psd_freq[0] / oversample
        self._f_fine = np.arange(self.psd_freq[0], self.psd_freq[-1] + df / 2, df)
        self._ls_fine = LombScargleBasis(self._t_sec, self._f_fine)

        self.fit = fit_red_white_whittle(self.psd_freq, self.psd_data)
        self.red_noise_amplitude = self.fit['red_amplitude']
        self.red_noise_spectral_index = self.fit['beta']
        self.white_noise_level = self.fit['white_level']
        self.psd_model_white = np.full_like(self.psd_freq, self.white_noise_level)
        self.psd_model_red = _model_psd(self.fit['theta'], self.psd_freq)   # total model

        if self.verbose:
            print(f"[{self.name}] Stage 1.0: Red/white decomposition (Lomb–Scargle + Whittle)")
            print(f"  Frequencies: {len(self.psd_freq)} bins, "
                  f"{self.psd_freq[0]:.3e}–{self.psd_freq[-1]:.3e} Hz")
            print(f"  White level C: {self.white_noise_level:.3e} s²/Hz")
            print(f"  Red amplitude A (at 1/yr): {self.red_noise_amplitude:.3e} s²/Hz, "
                  f"β = {self.red_noise_spectral_index:.2f}")
            if self.fit['beta_at_bound']:
                print("  ⚠ β at fit boundary — spectral index UNCONSTRAINED "
                      "(red power confined to the lowest 1–2 bins, or no red noise at all); "
                      "do not quote this β.")

        return self.psd_freq, self.psd_data, {
            'beta': self.red_noise_spectral_index,
            'white_noise_amplitude': self.white_noise_level,
            'red_noise_amplitude': self.red_noise_amplitude,
            'red_noise_model': self.psd_model_red,
            'beta_at_bound': self.fit['beta_at_bound'],
        }

    def _simulate_from_model(self, theta) -> np.ndarray:
        """Gaussian residuals with PSD = A f^-β (Fourier basis) + white C."""
        N, T = len(self._t_sec), self._t_sec[-1]
        f = self.psd_freq
        log_a, beta, log_c = theta
        s_red = 10 ** log_a * (f / F_YR) ** (-beta)
        sd = np.sqrt(s_red / T)                       # var per cos/sin coeff = S Δf
        a = self.rng.standard_normal(len(f)) * sd
        b = self.rng.standard_normal(len(f)) * sd
        red = self._sim_cos @ a + self._sim_sin @ b
        sigma_w = np.sqrt(10 ** log_c * N / (2 * T))  # E[S_white] = 2σ²T/N
        return red + sigma_w * self.rng.standard_normal(N)

    def stage1_spectral_anomaly_test(self, n_sigma: float = 3.0, n_sim: int = 200) -> Dict:
        """
        Stage 1.1: excess power at any frequency relative to the fitted
        red+white model. Statistic: max_k P_k / S_model,k.
        Calibrated by parametric bootstrap (simulate → re-fit → recompute).
        n_sigma only sets the "flagged bins" diagnostic count.
        """
        if self.fit is None:
            raise RuntimeError("Call stage1_psd_red_white_decomposition first")

        ratio = self._ls_fine.psd(self.residuals) / _model_psd(self.fit['theta'], self._f_fine)
        t_obs = np.max(ratio)

        t_null = np.empty(n_sim)
        for i in range(n_sim):
            y = self._simulate_from_model(self.fit['theta'])
            fit_i = fit_red_white_whittle(self.psd_freq, self._ls.psd(y))
            t_null[i] = np.max(self._ls_fine.psd(y) / _model_psd(fit_i['theta'], self._f_fine))
        p_value = permutation_p_value(t_obs, t_null)

        tail = 2 * stats.norm.sf(n_sigma)
        n_flagged = int(np.sum(ratio > -np.log(tail)))
        result = {
            'test': 'spectral_anomaly',
            'max_ratio': t_obs,
            'peak_frequency_hz': self._f_fine[np.argmax(ratio)],
            'peak_period_days': 1 / self._f_fine[np.argmax(ratio)] / SEC_PER_DAY,
            'n_flagged': n_flagged, 'n_freq': len(self._f_fine),
            'n_sim': n_sim, 'p_value': p_value,
        }
        self._register('spectral_anomaly', result)
        if self.verbose:
            print(f"[{self.name}] Stage 1.1: Spectral anomaly test")
            print(f"  max P/S_model = {t_obs:.2f} at {result['peak_period_days']:.1f} d, "
                  f"bins above {n_sigma}σ-equiv.: {n_flagged}/{len(self._f_fine)}, p = {p_value:.4g}")
        return result

    def stage1_frequency_dependence_test(self) -> Dict:
        """
        Stage 1.2: residual ∝ ν^-2 (DM-like chromatic delay).
        Weighted least squares r = a + b·ν^-2, t-test on b.
        Caveat: red noise correlates residuals in time; if different bands are
        observed at systematically different epochs this p-value is optimistic.
        """
        if self.residuals is None:
            raise RuntimeError("Call stage0_load_pulsar_data first")
        nu = self.frequencies
        if nu is None or np.ptp(nu) < 1e-6 * np.max(np.abs(nu)):
            result = {'test': 'frequency_dependence', 'p_value': None,
                      'note': 'single observing frequency — test skipped'}
            self.test_results['frequency_dependence'] = result
            if self.verbose:
                print(f"[{self.name}] Stage 1.2: skipped (single band)")
            return result

        x = (nu / 1400.0) ** -2
        w = 1 / self.errors ** 2
        X = np.column_stack([np.ones_like(x), x])
        XtW = X.T * w
        cov = np.linalg.inv(XtW @ X)
        coef = cov @ (XtW @ self.residuals)
        resid = self.residuals - X @ coef
        dof = len(x) - 2
        scale = np.sum(w * resid ** 2) / dof          # absorbs EFAC-like mis-scaling
        se_b = np.sqrt(cov[1, 1] * scale)
        t_stat = coef[1] / se_b
        p_value = 2 * stats.t.sf(abs(t_stat), dof)

        result = {'test': 'frequency_dependence', 'slope_s_at_1400MHz': coef[1],
                  'slope_se': se_b, 't_stat': t_stat, 'p_value': p_value,
                  'correlation': np.corrcoef(nu, self.residuals)[0, 1]}
        self._register('frequency_dependence', result)
        if self.verbose:
            print(f"[{self.name}] Stage 1.2: Frequency (ν^-2) dependence")
            print(f"  slope = {coef[1]*1e6:.4f} ± {se_b*1e6:.4f} µs, p = {p_value:.4g}")
        return result

    def stage1_autocorrelation_redness_test(self, n_permutations: int = 200) -> Dict:
        """
        Stage 1.3: is there red noise at all?  LR statistic
        2(nll_white − nll_red+white), calibrated by permuting r_i/σ_i across
        TOAs (exact under 'white noise with the quoted errors').
        NOTE: the stated H_0 of this class INCLUDES intrinsic red noise, so a
        rejection here is characterization, not necessarily an anomaly —
        consider excluding this test from the BH family (count_in_bh=False).
        """
        if self.fit is None:
            raise RuntimeError("Call stage1_psd_red_white_decomposition first")
        t_obs = self.fit['lr_stat']
        z = self.residuals / self.errors
        t_null = np.empty(n_permutations)
        for i in range(n_permutations):
            y = self.rng.permutation(z) * self.errors
            S = self._ls.psd(y)
            t_null[i] = fit_red_white_whittle(self.psd_freq, S)['lr_stat']
        p_value = permutation_p_value(t_obs, t_null)

        result = {'test': 'red_noise_spectral_index', 'spectral_index': self.fit['beta'],
                  'beta_at_bound': self.fit['beta_at_bound'],
                  'lr_stat': t_obs, 'n_permutations': n_permutations, 'p_value': p_value}
        self._register('red_noise', result)
        if self.verbose:
            print(f"[{self.name}] Stage 1.3: Red noise presence (LR, permutation)")
            flag = " (β at bound — unconstrained)" if self.fit['beta_at_bound'] else ""
            print(f"  β = {self.fit['beta']:.2f}{flag}, LR = {t_obs:.2f}, p = {p_value:.4g}")
        return result

    def stage1_measurement_error_consistency_test(self, n_fit_params: int = 0) -> Dict:
        """
        Stage 1.4: Σ(r/σ)² vs χ²(N − n_fit_params), two-sided
        (errors under- OR over-estimated). Set n_fit_params to the number of
        timing-model parameters that were fitted.
        """
        if self.residuals is None:
            raise RuntimeError("Call stage0_load_pulsar_data first")
        chi2_value = np.sum((self.residuals / self.errors) ** 2)
        dof = len(self.residuals) - n_fit_params
        p_value = min(1.0, 2 * min(stats.chi2.sf(chi2_value, dof), stats.chi2.cdf(chi2_value, dof)))
        result = {'test': 'measurement_error_consistency', 'chi2': chi2_value, 'dof': dof,
                  'reduced_chi2': chi2_value / dof, 'p_value': p_value}
        self._register('chi2', result)
        if self.verbose:
            print(f"[{self.name}] Stage 1.4: Measurement error consistency")
            print(f"  χ²/dof = {chi2_value/dof:.3f} (dof={dof}), p = {p_value:.4g}")
        return result

    # ------------------------------------------------------------------ Stage 2
    def stage2_benjamini_hochberg_correction(self, alpha_global: float = 0.05) -> Dict:
        if not self.p_values_raw:
            raise RuntimeError("No tests run yet")
        p_adj, reject = benjamini_hochberg(self.p_values_raw, alpha_global)
        self.p_values_corrected, self.reject_flags = p_adj, reject
        self.alpha_global = alpha_global
        result = {'test': 'benjamini_hochberg', 'n_tests': len(p_adj),
                  'alpha_global': alpha_global, 'reject_h0': bool(reject.any()),
                  'test_names': list(self.test_names),
                  'p_values_raw': np.array(self.p_values_raw), 'p_values_corrected': p_adj}
        if self.verbose:
            print(f"[{self.name}] Stage 2: Benjamini–Hochberg correction")
            print(f"  Tests: {len(p_adj)}, Reject H_0: {result['reject_h0']}")
        return result

    def identify_anomalies(self) -> List[Dict]:
        if self.p_values_corrected is None:
            raise RuntimeError("Call stage2_benjamini_hochberg_correction first")
        self.anomalies = [
            {'test': n, 'p_raw': p, 'p_corrected': float(q), 'details': self.test_results[n]}
            for n, p, q, r in zip(self.test_names, self.p_values_raw,
                                  self.p_values_corrected, self.reject_flags) if r
        ]
        return self.anomalies

    def to_result(self) -> TimingResidualAuditResult:
        return TimingResidualAuditResult(
            pulsar_name=self.name, n_observations=len(self.toas),
            time_span=(self.toas[-1] - self.toas[0]) / 365.25,
            p_values_raw=list(self.p_values_raw),
            p_values_corrected=[] if self.p_values_corrected is None else list(self.p_values_corrected),
            test_names=list(self.test_names), residual_rms=float(np.std(self.residuals)),
            red_noise_amplitude=self.red_noise_amplitude or 0.0,
            red_noise_spectral_index=self.red_noise_spectral_index or 0.0,
            white_noise_level=self.white_noise_level or 0.0,
            anomalies=self.anomalies, h0_rejected=bool(self.anomalies))

    def summary(self) -> str:
        lines = ["=" * 60, f"Pulsar Timing Residual Audit: {self.name}", "=" * 60,
                 f"Observations: {len(self.toas) if self.toas is not None else 'N/A'}",
                 f"Tests: {len(self.p_values_raw)}", ""]
        if self.p_values_corrected is not None:
            lines.append(f"After Benjamini–Hochberg (α={self.alpha_global}):")
            for n, p, q in zip(self.test_names, self.p_values_raw, self.p_values_corrected):
                lines.append(f"  {n:<22s} p_raw={p:.4g}  p_BH={q:.4g}")
            lines.append(f"Anomalies: {len(self.anomalies)}")
            if self.anomalies:
                for a in self.anomalies:
                    lines.append(f"  - {a['test']}: p_corr={a['p_corrected']:.4g}")
            else:
                lines.append("✓ No anomalies (H_0 not rejected)")
        return "\n".join(lines)


# ============================================================================
# Stage 3: common residual across pulsars
# ============================================================================

def test_common_residual_detection(residual_matrix: np.ndarray, psr_names: List[str],
                                   n_components: int = 1, verbose: bool = False,
                                   n_shifts: int = 500, seed: Optional[int] = None) -> Dict:
    """
    PC1 variance fraction vs. a null built from independent random circular
    shifts of each pulsar's series (preserves each pulsar's own spectrum,
    destroys cross-pulsar coherence). residual_matrix: (n_pulsars, n_epochs)
    on a COMMON epoch grid; NaN = missing.

    This detects *any* common/correlated component (clock, ephemeris, GWB…).
    It is NOT a Hellings–Downs test — it does not use sky positions.
    """
    rng = np.random.default_rng(seed)
    X = np.array(residual_matrix, dtype=float)
    X = (X - np.nanmean(X, axis=1, keepdims=True)) / (np.nanstd(X, axis=1, keepdims=True) + 1e-30)
    X = np.nan_to_num(X, nan=0.0)
    n_psr, n_ep = X.shape

    def pc_fractions(M):
        s = np.linalg.svd(M, compute_uv=False)
        return s ** 2 / np.sum(s ** 2)

    U, S, _ = np.linalg.svd(X, full_matrices=False)
    var_exp = S ** 2 / np.sum(S ** 2)

    t_null = np.empty(n_shifts)
    for i in range(n_shifts):
        shifts = rng.integers(0, n_ep, size=n_psr)
        Xs = np.stack([np.roll(X[k], shifts[k]) for k in range(n_psr)])
        t_null[i] = pc_fractions(Xs)[0]
    p_value = permutation_p_value(var_exp[0], t_null)

    result = {'test': 'common_residual', 'variance_pc1': var_exp[0],
              'variance_pc2': var_exp[1] if len(var_exp) > 1 else 0.0,
              'variance_explained': var_exp[:max(n_components, 2)],
              'null_pc1_mean': float(np.mean(t_null)), 'null_pc1_95': float(np.quantile(t_null, 0.95)),
              'n_pulsars': n_psr, 'pc1_pattern': U[:, 0], 'n_shifts': n_shifts,
              'p_value': p_value}
    if verbose:
        print("Common residual detection (circular-shift null):")
        print(f"  Variance in PC1: {var_exp[0]:.2%}  (null mean {np.mean(t_null):.2%}, "
              f"95% {np.quantile(t_null, 0.95):.2%})")
        if len(var_exp) > 1:
            print(f"  Variance in PC2: {var_exp[1]:.2%}")
        print(f"  p-value: {p_value:.4g}")
    return result


# ============================================================================
# Synthetic data + validation
# ============================================================================

def make_synthetic_pulsar(n_obs=451, span_days=5478.8, sigma_us=0.3, red_amp=0.0, beta=4.33,
                          sine_amp_us=0.0, sine_period_days=200.0, dm_slope_us=0.0,
                          rng=None) -> Dict:
    """Irregular TOAs, two bands, heteroscedastic errors; optional red noise / sine / DM."""
    rng = rng if rng is not None else np.random.default_rng()
    toas = 53000 + np.sort(rng.uniform(0, span_days, n_obs))
    errors = sigma_us * 1e-6 * rng.uniform(0.7, 1.3, n_obs)
    freqs = rng.choice([820.0, 1400.0, 2300.0], n_obs)
    res = errors * rng.standard_normal(n_obs)
    t = (toas - toas[0]) * SEC_PER_DAY
    if red_amp > 0:
        T = t[-1]; f = np.arange(1, n_obs // 2 + 1) / T
        sd = np.sqrt(red_amp * (f / F_YR) ** (-beta) / T)
        ph = 2 * np.pi * np.outer(t, f)
        res += np.cos(ph) @ (rng.standard_normal(len(f)) * sd) + np.sin(ph) @ (rng.standard_normal(len(f)) * sd)
    if sine_amp_us > 0:
        res += sine_amp_us * 1e-6 * np.sin(2 * np.pi * t / (sine_period_days * SEC_PER_DAY))
    if dm_slope_us:
        res += dm_slope_us * 1e-6 * (freqs / 1400.0) ** -2
    return dict(toas=toas, residuals=res, errors=errors, frequencies=freqs, dm=2.64, ra=None, dec=None)


def _audit(d, rng, n_sim=100, verbose=False):
    a = PulsarTimingResidualAudit("sim", verbose=verbose, rng=rng)
    a.stage0_load_pulsar_data(**d)
    a.stage1_psd_red_white_decomposition()
    a.stage1_spectral_anomaly_test(n_sim=n_sim)
    a.stage1_frequency_dependence_test()
    a.stage1_autocorrelation_redness_test(n_permutations=n_sim)
    a.stage1_measurement_error_consistency_test()
    a.stage2_benjamini_hochberg_correction()
    a.identify_anomalies()
    return a


def synthetic_validation(n_trials: int = 100, seed: int = 0):
    """
    FPR under white H_0 + sensitivity to red noise, a sinusoid and DM.
    Note: the red-noise LR statistic is exactly 0 in ~30% of white-noise
    realisations (fit collapses to white), giving a point mass at p = 1.
    That is valid (conservative) but makes a KS-uniformity check meaningless.
    """
    rng = np.random.default_rng(seed)
    P = []; n_rej = 0
    for _ in range(n_trials):
        a = _audit(make_synthetic_pulsar(rng=rng), rng)
        P.append(a.p_values_raw); n_rej += bool(a.anomalies)
    P = np.array(P)
    ci = stats.binomtest(n_rej, n_trials).proportion_ci()
    print(f"H_0 (white): global BH FPR = {n_rej/n_trials:.3f} "
          f"(95% CI {ci.low:.3f}–{ci.high:.3f}; target ≤ 0.05)")
    for j, n in enumerate(a.test_names):
        k = int(np.sum(P[:, j] <= 0.05))
        ci = stats.binomtest(k, n_trials).proportion_ci()
        print(f"  {n:<22s} raw FPR={k/n_trials:.3f}  (95% CI {ci.low:.3f}–{ci.high:.3f})")
    cases = {'red noise (β=4.33)': dict(red_amp=1e-8),
             'sine 0.15 µs, 200 d': dict(sine_amp_us=0.15),
             'DM ν^-2 0.1 µs': dict(dm_slope_us=0.1)}
    for label, kw in cases.items():
        a = _audit(make_synthetic_pulsar(rng=rng, **kw), rng)
        print(f"  {label:<20s} → anomalies: {[x['test'] for x in a.anomalies] or 'none'}")


if __name__ == "__main__":
    import sys
    from nanograv_loader import NANOGravDR15Loader

    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--real" in sys.argv:
        psr = args[0] if args else "J1909-3744"
        print(f"Pulsar Timing Residual Audit — REAL NANOGrav 15-yr data: {psr}\n")
        data = NANOGravDR15Loader(use_synthetic=False).load_pulsar(psr)
    else:
        psr = "J0437−4715"
        print("Pulsar Timing Residual Audit — Synthetic Test (seed=1)\n")
        data = NANOGravDR15Loader(use_synthetic=True, seed=1).load_pulsar(psr)

    d = dict(toas=data.toas, residuals=data.residuals, errors=data.errors,
             frequencies=data.frequencies, dm=data.dm, ra=data.ra, dec=data.dec)

    audit = PulsarTimingResidualAudit(psr, verbose=True, seed=1)
    audit.stage0_load_pulsar_data(**d)
    print()
    audit.stage1_psd_red_white_decomposition()
    audit.stage1_spectral_anomaly_test()
    audit.stage1_frequency_dependence_test()
    audit.stage1_autocorrelation_redness_test()
    audit.stage1_measurement_error_consistency_test()
    print()
    audit.stage2_benjamini_hochberg_correction()
    audit.identify_anomalies()
    print()
    print(audit.summary())

    if "--validate" in sys.argv:
        print("\nSynthetic validation\n" + "-" * 60)
        synthetic_validation()
