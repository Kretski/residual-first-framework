"""
Residual-First Discovery Framework — Common Statistical Methods

Permutation-based anomaly tests + Benjamini–Hochberg correction.

Null hypothesis of every Stage-1 test here: the residual is an exchangeable
(i.e. white, i.i.d.) sequence. Permuting the samples destroys temporal
structure while preserving the marginal distribution, so the permutation
p-values are exact under that null. NOTE: coloured (red) noise violates the
null and WILL be flagged — whiten the residual first if red noise is part of
the expected background.

Changes vs. previous version:
  * IndexError fix: the permutation PSD now uses its own frequency array and
    exactly the same Welch/periodogram settings as stage1_psd_fitting.
  * Permutation p-values use (1 + b) / (1 + n)  (Phipson & Smyth 2010),
    never exactly 0 — the old mean(T_perm >= T_obs) is anti-conservative.
  * BH-adjusted p-values were computed incorrectly (p reversed, ranks not);
    now correct and returned in the ORIGINAL test order.
  * summary(): no more "truth value of an array is ambiguous".
  * FPR self-test: under the global null BH gives P(any rejection) <= alpha,
    not ≈ alpha, so the pass criterion is now a one-sided binomial check.
  * Optional seed / Generator for reproducibility.

Author: Dimitar Kretski
License: MIT
"""

import numpy as np
from scipy import signal, stats
from typing import Dict, List, Optional


# ============================================================================
# Shared helpers
# ============================================================================

def permutation_p_value(t_obs: float, t_perm: np.ndarray) -> float:
    """Unbiased permutation p-value (Phipson & Smyth 2010): (1 + b) / (1 + n)."""
    t_perm = np.asarray(t_perm)
    return (1.0 + np.sum(t_perm >= t_obs)) / (1.0 + t_perm.size)


def benjamini_hochberg(p_values, alpha: float = 0.05):
    """
    Benjamini–Hochberg step-up procedure.

    Returns
    -------
    p_adj  : BH-adjusted p-values, in the SAME order as the input
    reject : boolean array, same order as the input
    """
    p = np.asarray(p_values, dtype=float)
    m = p.size
    order = np.argsort(p)
    ranks = np.arange(1, m + 1)

    q_sorted = p[order] * m / ranks
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]   # monotone from the top
    q_sorted = np.clip(q_sorted, 0.0, 1.0)

    p_adj = np.empty(m)
    p_adj[order] = q_sorted

    below = p[order] <= ranks / m * alpha
    reject = np.zeros(m, dtype=bool)
    if below.any():
        k = np.max(np.nonzero(below)[0])          # largest rank passing
        reject[order[:k + 1]] = True
    return p_adj, reject


# ============================================================================
# Audit class
# ============================================================================

class ResidualAudit:
    """Blind residual analysis with permutation-based anomaly detection."""

    def __init__(self, name: str, dt: float, verbose: bool = False,
                 seed: Optional[int] = None,
                 rng: Optional[np.random.Generator] = None):
        self.name = name
        self.dt = dt
        self.verbose = verbose
        self.rng = rng if rng is not None else np.random.default_rng(seed)

        self.residual = None
        self.psd_freq = None
        self.psd_data = None
        self._psd_method = None
        self._psd_nperseg = None

        self.test_results: Dict[str, Dict] = {}
        self.test_names: List[str] = []
        self.p_values_raw: List[float] = []
        self.p_values_corrected = None      # np.ndarray after stage 2, same order as p_values_raw
        self.reject_flags = None
        self.alpha_global = None
        self.anomalies: List[Dict] = []

    def _register(self, name: str, result: Dict):
        self.test_results[name] = result
        self.test_names.append(name)
        self.p_values_raw.append(float(result['p_value']))

    # ------------------------------------------------------------------ Stage 0
    def stage0_baseline_subtraction(self, data: np.ndarray, model: np.ndarray) -> np.ndarray:
        """Stage 0: Residual = data - model."""
        data, model = np.asarray(data, float), np.asarray(model, float)
        if data.shape != model.shape:
            raise ValueError("Data and model length mismatch")
        self.residual = data - model
        if self.verbose:
            print(f"[{self.name}] Stage 0: Residual computed, RMS = {np.std(self.residual):.6e}")
        return self.residual

    # ------------------------------------------------------------------ Stage 1
    def _psd(self, x: np.ndarray):
        """PSD with the settings fixed in stage1_psd_fitting; f = 0 removed."""
        if self._psd_method == "welch":
            f, p = signal.welch(x, fs=1 / self.dt, nperseg=self._psd_nperseg,
                                scaling='density')
        else:
            f, p = signal.periodogram(x, fs=1 / self.dt, scaling='density')
        keep = f > 0
        return f[keep], p[keep]

    def stage1_psd_fitting(self, method: str = "welch", nperseg: Optional[int] = None):
        """Compute the PSD and freeze its settings for the permutation test."""
        if self.residual is None:
            raise RuntimeError("Call stage0_baseline_subtraction first")
        if method not in ("welch", "periodogram"):
            raise ValueError("method must be 'welch' or 'periodogram'")
        self._psd_method = method
        self._psd_nperseg = nperseg if nperseg is not None else max(len(self.residual) // 8, 8)
        self.psd_freq, self.psd_data = self._psd(self.residual)
        if self.verbose:
            print(f"[{self.name}] Stage 1.0: PSD computed ({len(self.psd_freq)} bins)")
        return self.psd_freq, self.psd_data

    def stage1_spectral_anomaly_test_permutation(self, n_permutations: int = 200) -> Dict:
        """Max |PSD - mean(PSD)| vs. its permutation distribution."""
        if self.psd_data is None:
            raise RuntimeError("Call stage1_psd_fitting first")

        def stat(p):
            return np.max(np.abs(p - np.mean(p)))

        t_obs = stat(self.psd_data)
        t_perm = np.empty(n_permutations)
        for i in range(n_permutations):
            _, p_perm = self._psd(self.rng.permutation(self.residual))
            t_perm[i] = stat(p_perm)

        result = {
            'test': 'spectral_anomaly_permutation',
            'test_statistic': t_obs,
            'p_value': permutation_p_value(t_obs, t_perm),
            'n_permutations': n_permutations,
            'permutation_dist': t_perm,
            'peak_frequency': self.psd_freq[np.argmax(np.abs(self.psd_data - np.mean(self.psd_data)))],
        }
        self._register('spectral_anomaly', result)
        if self.verbose:
            print(f"[{self.name}] Stage 1.1: Spectral anomaly (permutation), p = {result['p_value']:.4g}")
        return result

    def stage1_timefreq_anomaly_test_permutation(self, n_permutations: int = 200) -> Dict:
        """Max spectrogram power vs. its permutation distribution."""
        if self.residual is None:
            raise RuntimeError("Call stage0_baseline_subtraction first")
        nperseg = min(len(self.residual) // 4, 256)

        def stat(x):
            _, _, sxx = signal.spectrogram(x, fs=1 / self.dt, nperseg=nperseg)
            return np.max(sxx)

        t_obs = stat(self.residual)
        t_perm = np.array([stat(self.rng.permutation(self.residual))
                           for _ in range(n_permutations)])
        result = {
            'test': 'timefreq_anomaly_permutation',
            'test_statistic': t_obs,
            'p_value': permutation_p_value(t_obs, t_perm),
            'n_permutations': n_permutations,
        }
        self._register('timefreq_anomaly', result)
        if self.verbose:
            print(f"[{self.name}] Stage 1.2: Time-frequency (permutation), p = {result['p_value']:.4g}")
        return result

    def stage1_autocorr_anomaly_test_permutation(self, max_lag: Optional[int] = None,
                                                  n_permutations: int = 200) -> Dict:
        """Sum of squared ACF (lags 1..max_lag-1) vs. its permutation distribution."""
        if self.residual is None:
            raise RuntimeError("Call stage0_baseline_subtraction first")
        N = len(self.residual)
        if max_lag is None:
            max_lag = N // 4

        def stat(x):
            x = x - np.mean(x)
            acf = np.fft.irfft(np.abs(np.fft.rfft(x, n=2 * N)) ** 2)[:N]
            acf = acf / acf[0]
            return np.sum(acf[1:max_lag] ** 2)

        t_obs = stat(self.residual)
        t_perm = np.array([stat(self.rng.permutation(self.residual))
                           for _ in range(n_permutations)])
        result = {
            'test': 'autocorr_anomaly_permutation',
            'test_statistic': t_obs,
            'p_value': permutation_p_value(t_obs, t_perm),
            'n_permutations': n_permutations,
        }
        self._register('autocorr_anomaly', result)
        if self.verbose:
            print(f"[{self.name}] Stage 1.3: Autocorrelation (permutation), p = {result['p_value']:.4g}")
        return result

    # ------------------------------------------------------------------ Stage 2
    def stage2_benjamini_hochberg_correction(self, alpha_global: float = 0.05) -> Dict:
        """Apply Benjamini–Hochberg FDR control across all Stage-1 tests."""
        if not self.p_values_raw:
            raise RuntimeError("No tests run yet")
        p_adj, reject = benjamini_hochberg(self.p_values_raw, alpha_global)
        self.p_values_corrected = p_adj
        self.reject_flags = reject
        self.alpha_global = alpha_global

        result = {
            'test': 'benjamini_hochberg',
            'n_tests': len(p_adj),
            'alpha_global': alpha_global,
            'reject_h0': bool(reject.any()),
            'test_names': list(self.test_names),
            'p_values_raw': np.array(self.p_values_raw),
            'p_values_corrected': p_adj,
            'n_significant': int(reject.sum()),
        }
        if self.verbose:
            print(f"[{self.name}] Stage 2: BH correction — tests: {result['n_tests']}, "
                  f"reject H_0: {result['reject_h0']}, significant: {result['n_significant']}")
        return result

    def identify_anomalies(self) -> List[Dict]:
        """Tests rejected by BH (uses the stage-2 alpha)."""
        if self.p_values_corrected is None:
            raise RuntimeError("Call stage2_benjamini_hochberg_correction first")
        self.anomalies = [
            {'test': name, 'p_raw': p_raw, 'p_corrected': float(p_corr)}
            for name, p_raw, p_corr, rej in zip(self.test_names, self.p_values_raw,
                                                self.p_values_corrected, self.reject_flags)
            if rej
        ]
        return self.anomalies

    def summary(self) -> str:
        rms = f"{np.std(self.residual):.6e}" if self.residual is not None else "N/A"
        lines = [
            "=" * 60,
            f"Residual Audit: {self.name}",
            "=" * 60,
            f"Residual RMS: {rms}",
            f"Tests: {len(self.p_values_raw)}",
            "",
        ]
        if self.p_values_corrected is not None:
            lines.append(f"After BH (α={self.alpha_global}):")
            for name, p, q in zip(self.test_names, self.p_values_raw, self.p_values_corrected):
                lines.append(f"  {name:<22s} p_raw={p:.4g}  p_BH={q:.4g}")
            lines.append(f"Anomalies: {len(self.anomalies)}")
            if self.anomalies:
                for a in self.anomalies:
                    lines.append(f"  - {a['test']}: p_BH={a['p_corrected']:.4g}")
            else:
                lines.append("✓ No anomalies (H_0 not rejected)")
        return "\n".join(lines)


# ============================================================================
# Synthetic test suite
# ============================================================================

def _run_full_audit(x, dt, rng, n_permutations):
    audit = ResidualAudit("synthetic", dt=dt, rng=rng)
    audit.stage0_baseline_subtraction(x, np.zeros_like(x))
    audit.stage1_psd_fitting()
    audit.stage1_spectral_anomaly_test_permutation(n_permutations)
    audit.stage1_timefreq_anomaly_test_permutation(n_permutations)
    audit.stage1_autocorr_anomaly_test_permutation(n_permutations=n_permutations)
    bh = audit.stage2_benjamini_hochberg_correction(alpha_global=0.05)
    audit.identify_anomalies()
    return audit, bh


def test_gaussian_noise_false_positive_rate(n_trials: int = 200, n_permutations: int = 200,
                                            alpha: float = 0.05, seed: int = 0,
                                            verbose: bool = True):
    """
    Under the global null, BH controls P(any rejection) <= alpha.
    PASS = observed FPR is not significantly ABOVE alpha
    (one-sided binomial test at the 1% level). A value below alpha is fine.
    """
    rng = np.random.default_rng(seed)
    dt, n_samples = 1e-3, 1000
    n_rejected = 0
    per_test_p = []

    for _ in range(n_trials):
        x = rng.standard_normal(n_samples)
        audit, bh = _run_full_audit(x, dt, rng, n_permutations)
        per_test_p.append(audit.p_values_raw)
        n_rejected += bh['reject_h0']

    fp_rate = n_rejected / n_trials
    p_excess = stats.binomtest(n_rejected, n_trials, alpha, alternative='greater').pvalue
    per_test_p = np.array(per_test_p)
    per_test_fpr = np.mean(per_test_p <= alpha, axis=0)
    ks_p = [stats.kstest(per_test_p[:, j], 'uniform').pvalue for j in range(per_test_p.shape[1])]

    if verbose:
        print(f"Trials: {n_trials}, permutations/test: {n_permutations}")
        print(f"Global FPR (BH, any rejection): {fp_rate:.3f}  "
              f"(95% CI {stats.binomtest(n_rejected, n_trials).proportion_ci().low:.3f}–"
              f"{stats.binomtest(n_rejected, n_trials).proportion_ci().high:.3f}; must be ≤ {alpha})")
        print(f"Per-test raw FPR at α={alpha}: " +
              ", ".join(f"{v:.3f}" for v in per_test_fpr))
        print(f"KS uniformity of raw p-values:  " + ", ".join(f"{v:.3f}" for v in ks_p))
        print(f"  → {'✓ PASS' if p_excess > 0.01 else '✗ FAIL (FPR significantly above α)'}")
    return fp_rate


def test_injected_anomaly_detection(anomaly_type: str = "spectral", snr: float = 5.0,
                                    n_permutations: int = 200, seed: int = 1) -> bool:
    """
    Sensitivity to an injected anomaly.
    snr is the per-sample amplitude ratio (signal amplitude / noise RMS).
    """
    rng = np.random.default_rng(seed)
    dt, n_samples = 1e-3, 1000
    t = np.arange(n_samples) * dt
    noise = rng.standard_normal(n_samples)
    amp = snr * np.std(noise)

    if anomaly_type == "spectral":
        anomaly = amp * np.sin(2 * np.pi * 100 * t)
    elif anomaly_type == "coherent":
        anomaly = amp * signal.chirp(t, 50, t[-1], 150)
    else:
        raise ValueError(f"Unknown anomaly type: {anomaly_type}")

    audit, _ = _run_full_audit(noise + anomaly, dt, rng, n_permutations)
    return len(audit.anomalies) > 0


def power_curve(kind: str = "spectral", snrs=(0.1, 0.2, 0.3, 0.5, 1.0),
                n_trials: int = 40, n_permutations: int = 200):
    """Detection probability vs. SNR (optional, not run by default)."""
    for snr in snrs:
        rate = np.mean([test_injected_anomaly_detection(kind, snr, n_permutations, seed=s)
                        for s in range(n_trials)])
        print(f"  {kind:<9s} SNR={snr:<4}  detection rate = {rate:.2f}")


if __name__ == "__main__":
    print("Residual-First Common Module — Synthetic Tests\n")

    print("Test 1: False positive rate under H_0")
    print("-" * 60)
    test_gaussian_noise_false_positive_rate(n_trials=200, n_permutations=200)
    print()

    for label, kind, snr in [("Test 2: Injected spectral anomaly (SNR=0.3)", "spectral", 0.3),
                             ("Test 3: Injected coherent chirp (SNR=1.0)", "coherent", 1.0),
                             ("Test 4: Strong spectral anomaly (SNR=5)", "spectral", 5.0)]:
        print(label)
        print("-" * 60)
        print(f"  Result: {'✓ DETECTED' if test_injected_anomaly_detection(kind, snr) else '✗ MISSED'}")
        print()

