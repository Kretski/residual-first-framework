"""
Common Residual Search Across Pulsar Ensembles (Stage 3)

Tests (all calibrated with the same null):
  1. PCA        — is PC1 larger than expected for independent pulsars?
  2. Cross-corr — is there ANY inter-pulsar correlation?
  3. Sky / ORF  — do pairwise correlations follow Hellings–Downs
                  (also reports monopole = clock-like, dipole = ephemeris-like)?
  4. Spectrum   — is the common (PC1) mode redder than for independent pulsars?

NULL: every pulsar is independent. Realised by random CIRCULAR SHIFTS of each
pulsar's series, which keeps each pulsar's own spectrum (incl. intrinsic red
noise) and destroys only cross-pulsar coherence. Plain permutation of epochs
(previous version) destroys the red noise itself, so any red noise looked
like a "common signal".

Changes vs. previous version:
  * PCA null: circular shifts instead of per-pulsar shuffling.
  * Cross-correlation null: previous code applied the SAME permutation to all
    pulsars, which preserves the cross-correlations → the test could never fire.
  * Sky test: previous monopole/dipole regression used PC1 loadings whose sign
    is arbitrary, and treated RA linearly. Now: pairwise correlations projected
    on the monopole, dipole and Hellings–Downs overlap-reduction functions.
  * Spectrum test: previous code took the PSD of U[:,0] (one number PER PULSAR,
    length 5), not of the PC1 time series (Vt[0]) — hence the RankWarning and
    β = −1.97. Also the sign of β was inverted.
  * BH-adjusted p-values fixed ("Significant: 1, Reject H_0: False" was a
    symptom of that bug).
  * Permutation p-values use (1 + b)/(1 + n); NaN-aware standardisation.

Reference: NANOGrav Collaboration (2023), ApJL 951, L8.

Author: Dimitar Kretski
License: MIT
"""

import numpy as np
from scipy import stats
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass


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


def hellings_downs(cos_zeta: np.ndarray) -> np.ndarray:
    """Hellings–Downs ORF for distinct pulsars (Γ = 0.5 at ζ → 0, −0.125 minimum)."""
    x = np.clip((1 - cos_zeta) / 2, 1e-12, 1.0)
    return 0.5 - x / 4 + 1.5 * x * np.log(x)


def radec_to_unit(ra_deg, dec_deg) -> np.ndarray:
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


@dataclass
class CommonResidualResult:
    n_pulsars: int
    n_epochs: int
    variance_explained: np.ndarray
    pc1_pattern: np.ndarray
    pc1_amplitude: float
    p_value_pc1: float
    p_value_sky_correlation: float
    p_value_cross_coherence: float
    h0_rejected: bool
    significant_components: int


# ============================================================================
# Ensemble search
# ============================================================================

class PulsarEnsembleCommonResidualSearch:
    """
    residual_matrix: (n_pulsars, n_epochs) on a common, uniform epoch grid.
    psr_names MUST be in the same order as the rows of residual_matrix.
    """

    def __init__(self, residual_matrix: np.ndarray, psr_names: List[str],
                 psr_coords: Dict[str, Tuple[float, float]], verbose: bool = False,
                 n_null: int = 500, seed: Optional[int] = None):
        X = np.array(residual_matrix, dtype=float)
        if X.shape[0] != len(psr_names):
            raise ValueError("psr_names length must equal number of matrix rows")
        mu = np.nanmean(X, axis=1, keepdims=True)
        sd = np.nanstd(X, axis=1, keepdims=True)
        self.missing_fraction = float(np.mean(np.isnan(X)))
        self.X = np.nan_to_num((X - mu) / np.where(sd > 0, sd, 1.0), nan=0.0)
        self.residual_matrix = np.nan_to_num(X, nan=0.0)   # kept for compatibility

        self.psr_names = list(psr_names)
        self.psr_coords = psr_coords
        self.verbose = verbose
        self.n_null = n_null
        self.rng = np.random.default_rng(seed)
        self.n_pulsars, self.n_epochs = X.shape

        self.test_names: List[str] = []
        self.p_values: List[float] = []
        self.pca_results = None
        self.sky_correlation_matrix = None
        self.bh_result = None
        self._null_cache = None

    # ------------------------------------------------------------------ utils
    def _register(self, name, p):
        self.test_names.append(name)
        self.p_values.append(float(p))

    def _null_matrices(self):
        """Independent random circular shift per pulsar; generated once, reused by all tests."""
        if self._null_cache is None:
            n = self.n_epochs
            idx = np.arange(n)
            self._null_cache = []
            for _ in range(self.n_null):
                shifts = self.rng.integers(0, n, size=self.n_pulsars)
                self._null_cache.append(self.X[np.arange(self.n_pulsars)[:, None],
                                               (idx[None, :] - shifts[:, None]) % n])
        return self._null_cache

    @staticmethod
    def _corr_matrix(M):
        Z = M - M.mean(axis=1, keepdims=True)
        Z /= np.linalg.norm(Z, axis=1, keepdims=True) + 1e-30
        return Z @ Z.T

    def _pair_indices(self):
        return np.triu_indices(self.n_pulsars, k=1)

    # ------------------------------------------------------------------ Test 1
    def test_principal_components(self, n_components: int = 3) -> Dict:
        """PC1 variance fraction vs. circular-shift null."""
        U, S, Vt = np.linalg.svd(self.X, full_matrices=False)
        var = S ** 2 / np.sum(S ** 2)

        def pc1_frac(M):
            s = np.linalg.svd(M, compute_uv=False)
            return s[0] ** 2 / np.sum(s ** 2)

        null = np.array([pc1_frac(M) for M in self._null_matrices()])
        p = permutation_p_value(var[0], null)

        self.pca_results = {
            'test': 'pca_principal_components',
            'n_components_analyzed': n_components,
            'variance_explained': var[:n_components],
            'singular_values': S[:n_components],
            'pc1_pattern': U[:, 0],               # per-pulsar loadings
            'pc1_timeseries': Vt[0] * S[0],      # common time series
            'pc1_amplitude': S[0],
            'pc1_bootstrap_dist': null,
            'null_mean': float(null.mean()), 'null_95': float(np.quantile(null, 0.95)),
            'p_value': p, 'pca_U': U, 'pca_S': S, 'pca_Vt': Vt,
        }
        self._register('pca_pc1', p)
        if self.verbose:
            print("PCA Test:")
            print(f"  Variance PC1: {var[0]:.2%}  (null mean {null.mean():.2%}, "
                  f"95%: {np.quantile(null, 0.95):.2%})")
            if len(var) > 1:
                print(f"  Variance PC2: {var[1]:.2%}")
            print(f"  p-value (circular-shift null, n={self.n_null}): {p:.4g}")
        return self.pca_results

    # ------------------------------------------------------------------ Test 2
    def test_cross_pulsar_correlations(self) -> Dict:
        """Omnibus: mean squared pairwise correlation vs. circular-shift null."""
        iu = self._pair_indices()
        C = self._corr_matrix(self.X)
        rho = C[iu]
        t_obs = np.mean(rho ** 2)
        null = np.array([np.mean(self._corr_matrix(M)[iu] ** 2) for M in self._null_matrices()])
        p = permutation_p_value(t_obs, null)
        self.sky_correlation_matrix = C

        result = {'test': 'cross_pulsar_correlations', 'n_pairs': len(rho),
                  'mean_correlation': float(np.mean(rho)), 'std_correlation': float(np.std(rho)),
                  'mean_sq_correlation': float(t_obs), 'p_value': p,
                  'p_value_permutation': p, 'all_correlations': rho}
        self._register('cross_correlation', p)
        if self.verbose:
            print("Cross-Pulsar Correlation Test:")
            print(f"  Mean ρ: {np.mean(rho):.4f}, mean ρ²: {t_obs:.4f} "
                  f"(null mean {null.mean():.4f})")
            print(f"  p-value: {p:.4g}")
        return result

    # ------------------------------------------------------------------ Test 3
    def test_sky_correlation_pattern(self) -> Dict:
        """
        Project pairwise correlations onto monopole (clock), dipole (ephemeris)
        and Hellings–Downs (GWB) templates. The registered p-value is the
        one-sided HD test; the other two are reported as diagnostics.
        With few pulsars the three templates are strongly degenerate.
        """
        have = [n in self.psr_coords for n in self.psr_names]
        if sum(have) < 3:
            if self.verbose:
                print("Sky Correlation Test: skipped (<3 pulsars with coordinates)")
            return {'test': 'sky_correlation', 'p_value': None, 'note': 'insufficient coordinates'}

        rows = np.nonzero(have)[0]
        coords = np.array([self.psr_coords[self.psr_names[i]] for i in rows])
        n_hat = radec_to_unit(coords[:, 0], coords[:, 1])
        cosz_full = n_hat @ n_hat.T
        iu = np.triu_indices(len(rows), k=1)
        cosz = cosz_full[iu]
        templates = {'monopole': np.ones_like(cosz), 'dipole': cosz,
                     'hellings_downs': hellings_downs(cosz)}

        def projections(M):
            rho = self._corr_matrix(M[rows])[iu]
            return {k: np.dot(rho, g) / np.sqrt(np.dot(g, g)) for k, g in templates.items()}

        obs = projections(self.X)
        null = [projections(M) for M in self._null_matrices()]
        pvals = {k: permutation_p_value(obs[k], [d[k] for d in null]) for k in templates}
        hd_mono_overlap = np.dot(templates['hellings_downs'], templates['monopole']) / (
            np.linalg.norm(templates['hellings_downs']) * np.linalg.norm(templates['monopole']))

        result = {'test': 'sky_correlation_orf', 'n_pairs': len(cosz),
                  'angular_separations_deg': np.degrees(np.arccos(np.clip(cosz, -1, 1))),
                  'projections': obs, 'p_monopole': pvals['monopole'],
                  'p_dipole': pvals['dipole'], 'p_hellings_downs': pvals['hellings_downs'],
                  'hd_monopole_template_overlap': hd_mono_overlap,
                  'p_value': pvals['hellings_downs']}
        self._register('sky_hellings_downs', pvals['hellings_downs'])
        if self.verbose:
            print("Sky Correlation (ORF) Test:")
            for k in templates:
                print(f"  {k:<15s} projection = {obs[k]:+.3f}, p = {pvals[k]:.4g}")
            print(f"  (HD–monopole template overlap: {hd_mono_overlap:.2f}; "
                  f"close to ±1 ⇒ the two cannot be distinguished with these pulsars)")
        return result

    # ------------------------------------------------------------------ Test 4
    def test_common_residual_frequency_structure(self, low_fraction: float = 0.1) -> Dict:
        """
        Is the common mode (PC1 time series) redder than the PC1 of
        INDEPENDENT pulsars with the same individual spectra?
        Statistic: fraction of PC1 periodogram power in the lowest
        `low_fraction` of frequencies. Null: same statistic on PC1 of the
        circular-shifted matrices (intrinsic red noise is therefore part of
        the null — permuting PC1 in time would flag it as a false positive).
        """
        if self.pca_results is None:
            raise RuntimeError("Call test_principal_components first")
        ts = self.pca_results['pc1_timeseries']
        n_low = max(1, int(round(low_fraction * (len(ts) // 2))))

        def low_power_frac(x):
            P = np.abs(np.fft.rfft(x - x.mean()))[1:] ** 2
            return P[:n_low].sum() / P.sum()

        def pc1_ts(M):
            _, S, Vt = np.linalg.svd(M, full_matrices=False)
            return Vt[0] * S[0]

        t_obs = low_power_frac(ts)
        null = np.array([low_power_frac(pc1_ts(M)) for M in self._null_matrices()])
        p = permutation_p_value(t_obs, null)

        # descriptive log-log slope (β = −slope, S ∝ f^−β)
        f = np.fft.rfftfreq(len(ts))[1:]
        P = np.abs(np.fft.rfft(ts - ts.mean()))[1:] ** 2
        beta = -np.polyfit(np.log10(f), np.log10(P + 1e-300), 1)[0]

        result = {'test': 'common_residual_frequency_structure', 'spectral_index_beta': beta,
                  'low_freq_power_fraction': t_obs, 'null_mean': float(null.mean()),
                  'p_value': p, 'psd_frequencies': f, 'psd_data': P}
        self._register('pc1_redness', p)
        if self.verbose:
            print("Frequency Structure Test (PC1 time series):")
            print(f"  Power in lowest {n_low} bins: {t_obs:.2%} "
                  f"(independent-pulsar null: {null.mean():.2%}), descriptive β ≈ {beta:.2f}")
            print(f"  p-value: {p:.4g}")
        return result

    # ------------------------------------------------------------------ BH
    def stage3_benjamini_hochberg_correction(self, alpha_global: float = 0.05) -> Dict:
        if not self.p_values:
            raise RuntimeError("No tests run yet")
        p_adj, reject = benjamini_hochberg(self.p_values, alpha_global)
        self.bh_result = {'test': 'ensemble_benjamini_hochberg', 'n_tests': len(p_adj),
                          'alpha_global': alpha_global, 'reject_h0': bool(reject.any()),
                          'test_names': list(self.test_names),
                          'p_values_raw': np.array(self.p_values), 'p_values_corrected': p_adj,
                          'reject': reject, 'n_significant': int(reject.sum())}
        if self.verbose:
            print("\nEnsemble BH Correction:")
            for n, p, q, r in zip(self.test_names, self.p_values, p_adj, reject):
                print(f"  {n:<20s} p_raw={p:.4g}  p_BH={q:.4g}  {'← significant' if r else ''}")
            print(f"  Reject H_0: {self.bh_result['reject_h0']}")
        return self.bh_result

    def summary(self) -> str:
        lines = ["=" * 70, "Pulsar Ensemble Common Residual Search", "=" * 70,
                 f"Pulsars: {self.n_pulsars}", f"Epochs: {self.n_epochs}",
                 f"Missing (interpolation gaps): {self.missing_fraction:.1%}",
                 f"Tests: {len(self.p_values)}", ""]
        if self.pca_results is not None:
            lines.append("PCA Results:")
            lines.append(f"  PC1 variance: {self.pca_results['variance_explained'][0]:.2%}")
            lines.append(f"  PC1 p-value: {self.pca_results['p_value']:.4g}")
        if self.bh_result is not None:
            lines.append(f"After BH: {self.bh_result['n_significant']} significant, "
                         f"reject H_0 = {self.bh_result['reject_h0']}")
        return "\n".join(lines)


# ============================================================================
# Replication across arrays
# ============================================================================

def compare_nanograv_ipta(nanograv_matrix: np.ndarray, ipta_matrix: np.ndarray,
                          verbose: bool = True) -> Dict:
    """Independent PCA tests on two arrays (PC1 patterns are not compared:
    different pulsar sets)."""
    ng = PulsarEnsembleCommonResidualSearch(
        nanograv_matrix, [f"NANOGrav_{i}" for i in range(nanograv_matrix.shape[0])], {}, verbose)
    ip = PulsarEnsembleCommonResidualSearch(
        ipta_matrix, [f"IPTA_{i}" for i in range(ipta_matrix.shape[0])], {}, verbose)
    ng.test_principal_components()
    ip.test_principal_components()
    result = {'test': 'replication_across_arrays',
              'nanograv_pc1_pval': ng.pca_results['p_value'],
              'ipta_pc1_pval': ip.pca_results['p_value'],
              'both_significant': ng.pca_results['p_value'] < 0.05 and ip.pca_results['p_value'] < 0.05}
    if verbose:
        print("\n" + "=" * 70 + "\nReplication Test (NANOGrav vs IPTA)\n" + "=" * 70)
        print(f"NANOGrav PC1 p-value: {result['nanograv_pc1_pval']:.4g}")
        print(f"IPTA PC1 p-value: {result['ipta_pc1_pval']:.4g}")
        print(f"Both significant: {result['both_significant']}")
    return result


def _run_all(matrix, names, coords, verbose=False, seed=None, n_null=500):
    s = PulsarEnsembleCommonResidualSearch(matrix, names, coords, verbose=verbose,
                                           seed=seed, n_null=n_null)
    s.test_principal_components()
    s.test_cross_pulsar_correlations()
    s.test_sky_correlation_pattern()
    s.test_common_residual_frequency_structure()
    s.stage3_benjamini_hochberg_correction()
    return s


if __name__ == "__main__":
    import sys
    from nanograv_loader import NANOGravDR15Loader, PulsarTimingEnsemble

    if "--real" in sys.argv:
        print("Common Residual Search — REAL NANOGrav 15-yr data\n")
        loader = NANOGravDR15Loader(use_synthetic=False)
        pulsars = loader.load_all_pulsars()
        pulsars = {n: d for n, d in pulsars.items() if np.isfinite(d.ra) and np.isfinite(d.dec)}
        print(f"Loaded {len(pulsars)} pulsars with sky positions")
    else:
        print("Common Residual Search — Synthetic Ensemble Test (seed=1)\n")
        psr_list = ['J0023+0923', 'J0030+0451', 'J0340+4130', 'J0437−4715', 'J0613−0200']
        loader = NANOGravDR15Loader(use_synthetic=True, seed=1)
        pulsars = {n: loader.load_pulsar(n) for n in psr_list}

    ensemble = PulsarTimingEnsemble(pulsars)
    matrix = ensemble.common_residual_matrix(n_bins=100)
    names = ensemble.pulsar_names()                     # same order as matrix rows
    coords = {n: (d.ra, d.dec) for n, d in pulsars.items()}

    search = _run_all(matrix, names, coords, verbose=True, seed=0)
    print()
    print(search.summary())

    if "--validate" in sys.argv:
        print("\n" + "=" * 70 + "\nValidation\n" + "=" * 70)
        for label, amp, n_trials in [("H_0: no common signal", 0.0, 200),
                                     ("common signal 0.5 µs", 0.5e-6, 20)]:
            rej = np.zeros(4)
            glob = 0
            psr_list = ['J0023+0923', 'J0030+0451', 'J0340+4130', 'J0437−4715', 'J0613−0200']
            for k in range(n_trials):
                ld = NANOGravDR15Loader(use_synthetic=True, seed=100 + k, common_amplitude=amp)
                ps = {n: ld.load_pulsar(n) for n in psr_list}
                ens = PulsarTimingEnsemble(ps)
                s = _run_all(ens.common_residual_matrix(n_bins=100), ens.pulsar_names(),
                             {n: (d.ra, d.dec) for n, d in ps.items()}, seed=k, n_null=200)
                rej += np.array(s.p_values) <= 0.05
                glob += s.bh_result['reject_h0']
            print(f"{label}: BH global rejection rate = {glob/n_trials:.2f}")
            print("  raw rate p≤0.05 per test: " +
                  ", ".join(f"{n}={r/n_trials:.2f}" for n, r in zip(s.test_names, rej)))
