# Pulsar Timing Residual-First Audit
## Complete Workflow for NANOGrav DR15 Independent Validation

**Objective:** Independent verification of NANOGrav common residual detection using residual-first framework.

**Status:** Synthetic test suite ready; replace with real NANOGrav DR15 data when available.

**Timeline:** 3–4 weeks for Stages 0–3.

---

## Architecture Overview

```
┌──────────────────────────────────────────┐
│ STAGE 0: Load NANOGrav DR15 Pulsar Data  │
│ nanograv_loader.py                       │
│ Input: PSR names                         │
│ Output: PulsarTimingData objects         │
└─────────────────┬────────────────────────┘
                  │
┌─────────────────▼────────────────────────┐
│ STAGE 1: Independent Pulsar Audits       │
│ pulsar_timing_residual_audit.py          │
│ • PSD fit (red + white noise)            │
│ • Spectral anomaly test                  │
│ • Frequency dependence (DM variation)    │
│ • Red noise spectral index test          │
│ • χ² measurement error consistency       │
│ Output: p-values per test                │
└─────────────────┬────────────────────────┘
                  │
┌─────────────────▼────────────────────────┐
│ STAGE 2: Benjamini–Hochberg Correction   │
│ residual_first_common.py                 │
│ • Multiple testing control (FDR)         │
│ • Identify anomalies (p_corr < 0.05)     │
│ • Document null results                  │
└─────────────────┬────────────────────────┘
                  │
┌─────────────────▼────────────────────────┐
│ STAGE 3: Ensemble Replication            │
│ common_residual_search.py                │
│ • PCA (common signal detection)          │
│ • Cross-pulsar correlations              │
│ • Sky isotropy test                      │
│ • Red noise in common component          │
│ • Replication: NANOGrav vs IPTA          │
│ Output: H_0 rejection decision           │
└──────────────────────────────────────────┘
```

---

## Stage 0: Data Loading

### Single Pulsar

```python
from nanograv_loader import NANOGravDR15Loader

loader = NANOGravDR15Loader(use_synthetic=False)
# Requires manual download from https://data.nanoGrav.org/

psr_data = loader.load_pulsar('J0437−4715')

print(f"Pulsar: {psr_data.name}")
print(f"Observations: {psr_data.n_observations}")
print(f"Time span: {psr_data.time_span} years")
print(f"Residual RMS: {psr_data.residual_rms()*1e6:.3f} µs")
```

### Multiple Pulsars (Ensemble)

```python
from nanograv_loader import NANOGravDR15Loader, PulsarTimingEnsemble

loader = NANOGravDR15Loader(use_synthetic=True)

psr_names = ['J0023+0923', 'J0030+0451', 'J0340+4130', ...]
pulsars = {name: loader.load_pulsar(name) for name in psr_names}

ensemble = PulsarTimingEnsemble(pulsars)
residual_matrix = ensemble.common_residual_matrix(n_bins=100)
# Shape: (n_pulsars, 100 common epochs)
```

---

## Stage 1: Individual Pulsar Audits

### Null Hypothesis (Written Before Opening Data)

```yaml
# pulsar_audit_hypothesis.yaml
null_hypothesis: |
  Timing residuals are consistent with:
  - Standard timing model (TTM)
  - White measurement noise
  - Individual red noise (spin-down related)
  - NO common signal correlated across pulsars

baseline_model: "TEMPO2/TEMPO standard timing model"

tests:
  - name: "Red noise spectral index"
    method: "PSD fitting in log-log space"
    null: "β = 0 (white noise only)"
    
  - name: "Frequency dependence"
    method: "Correlation(residuals, observation frequency)"
    null: "r = 0 (no DM variation)"
    
  - name: "Spectral anomaly"
    method: "Deviation from red+white fit"
    null: "No excess power in any frequency band"
    
  - name: "χ² consistency"
    method: "(residuals / errors)² vs chi-square"
    null: "χ²/dof ≈ 1"
```

### Run Single Audit

```python
from pulsar_timing_residual_audit import PulsarTimingResidualAudit

audit = PulsarTimingResidualAudit('J0437−4715', verbose=True)

# Load data
audit.stage0_load_pulsar_data(
    toas=psr_data.toas,
    residuals=psr_data.residuals,
    errors=psr_data.errors,
    frequencies=psr_data.frequencies,
    dm=psr_data.dm,
    ra=psr_data.ra,
    dec=psr_data.dec
)

# Run predefined tests (blind)
audit.stage1_psd_red_white_decomposition()
audit.stage1_spectral_anomaly_test()
audit.stage1_frequency_dependence_test()
audit.stage1_autocorrelation_redness_test()
audit.stage1_measurement_error_consistency_test()

# Apply multiple testing correction
audit.stage2_benjamini_hochberg_correction(alpha_global=0.05)

# Identify anomalies
anomalies = audit.identify_anomalies()

# Print summary
print(audit.summary())
```

**Expected Output (under H_0):**

```
============================================================
Pulsar Timing Residual Audit: J0437−4715
============================================================
Observations: 451
Tests: 5

After Benjamini–Hochberg (α=0.05):
Anomalies: 0
✓ No anomalies (H_0 not rejected)
```

---

## Stage 2: Catalog Scan

**For each of 70 pulsars:**

```python
results = {}

for psr_name in psr_list:
    try:
        data = loader.load_pulsar(psr_name)
        audit = PulsarTimingResidualAudit(psr_name, verbose=False)
        audit.stage0_load_pulsar_data(
            data.toas, data.residuals, data.errors,
            data.frequencies, data.dm, data.ra, data.dec
        )
        audit.stage1_psd_red_white_decomposition()
        audit.stage1_spectral_anomaly_test()
        audit.stage1_frequency_dependence_test()
        audit.stage1_autocorrelation_redness_test()
        audit.stage1_measurement_error_consistency_test()
        audit.stage2_benjamini_hochberg_correction()
        
        results[psr_name] = {
            'h0_rejected': len(audit.identify_anomalies()) > 0,
            'n_anomalies': len(audit.identify_anomalies()),
            'anomalies': audit.identify_anomalies()
        }
    except Exception as e:
        print(f"Failed: {psr_name}: {e}")

# Aggregate
n_with_anomalies = sum(1 for r in results.values() if r['h0_rejected'])
print(f"Pulsars with anomalies: {n_with_anomalies}/{len(psr_list)}")
```

**Expected Result (under H_0):**
- Most pulsars: 0 anomalies
- 5% of pulsars: 1 anomaly (false positive by design)
- This is consistent with α=0.05

---

## Stage 3: Ensemble Replication

### Test 1: Principal Component Analysis

```python
from common_residual_search import PulsarEnsembleCommonResidualSearch

# Ensemble residual matrix (already constructed in Stage 0)
# Shape: (n_pulsars, n_epochs)

psr_coords = {name: (data.ra, data.dec) for name, data in pulsars.items()}

search = PulsarEnsembleCommonResidualSearch(
    residual_matrix, list(pulsars.keys()), psr_coords,
    verbose=True
)

# Null Hypothesis (written before running):
# First principal component amplitude is consistent with
# expected noise level given number of pulsars.

pca_result = search.test_principal_components()

print(f"Variance in PC1: {pca_result['variance_explained'][0]:.2%}")
print(f"PC1 p-value: {pca_result['p_value']:.6e}")
```

**Expected under H_0:** p > 0.05 (no common signal)

**Expected if true common signal:** p < 0.01 (strong detection)

### Test 2: Cross-Pulsar Correlations

```python
corr_result = search.test_cross_pulsar_correlations()

print(f"Mean pairwise correlation: {corr_result['mean_correlation']:.4f}")
print(f"p-value: {corr_result['p_value_permutation']:.6e}")
```

### Test 3: Sky Isotropy

```python
sky_result = search.test_sky_correlation_pattern()

print(f"Monopole SNR: {sky_result['monopole_snr']:.3f}")
print(f"Dipole F-stat: {sky_result['dipole_f_stat']:.3f}")
print(f"p-value: {sky_result['p_value']:.6e}")
```

**Interpretation:**
- p > 0.05 → Common signal is isotropic (consistent with stochastic GW background)
- p < 0.05 → Common signal is directional (possibly pulsar-noise artifact)

### Test 4: Frequency Structure of Common Signal

```python
freq_result = search.test_common_residual_frequency_structure()

print(f"PC1 spectral index β: {freq_result['spectral_index_beta']:.2f}")
print(f"p-value: {freq_result['p_value']:.6e}")
```

**Interpretation:**
- β ≈ 0 → White noise (unlikely for real GW background)
- β ≈ 2–3 → Red noise (consistent with GW background)
- β > 3 → Very red (possible AGN or other low-frequency source)

### Ensemble-Level Multiple Testing Correction

```python
ensemble_result = search.stage3_benjamini_hochberg_correction(alpha_global=0.05)

print(f"Significant tests: {ensemble_result['n_significant']}")
print(f"Reject H_0: {ensemble_result['reject_h0']}")
```

---

## Stage 3 Replication: Cross-Array Validation

**Key historical precedent:** NANOGrav found common residual in their data; IPTA independently confirmed it.

```python
from common_residual_search import compare_nanograv_ipta

# Load both arrays (currently synthetic)
nanograv_data = loader.load_all_pulsars()  # NANOGrav DR15
nanograv_matrix = ensemble_nanograv.common_residual_matrix()

ipta_data = ipta_loader.load_all_pulsars()  # IPTA DR1
ipta_matrix = ensemble_ipta.common_residual_matrix()

# Compare
replication = compare_nanograv_ipta(nanograv_matrix, ipta_matrix)

print(f"NANOGrav PC1 p-value: {replication['nanograv_pc1_pval']:.6e}")
print(f"IPTA PC1 p-value: {replication['ipta_pc1_pval']:.6e}")
print(f"Both significant: {replication['both_significant']}")
```

**Success Criterion (for "finding something"):**
- NANOGrav: PC1 p-value < 0.05 ✓
- IPTA: PC1 p-value < 0.05 ✓
- Same morphology (visual inspection of PC1 patterns)

**If true:**
- Publish: "Independent replication of NANOGrav common residual"
- Characterize: Spectral index, sky correlation, polarization (if data allows)
- Interpret: GW background? Systematics? Other?

---

## Stopping Rules and Null Results

### Single Pulsar (Stage 1–2)

**Stop if:**
- All tests pass BH correction with p_corr > 0.05
- Report: "[PSR name]: No timing anomalies detected"
- Null result is publishable as methodology validation

### Ensemble (Stage 3)

**Stop if:**
- PC1 p-value > 0.05 (NANOGrav)
- Report: "Independent audit of GWTC-4 timing residuals finds no common signal"

**Continue if:**
- PC1 p-value < 0.05 (NANOGrav AND IPTA)
- Proceed to replication and characterization

### Lambda Interpretation (Only if Stage 3 succeeds)

**Stage 4 (conditional):**
- Do not test Lambda unless ensemble shows robust, replicated common residual
- Null hypothesis: "Can standard GW dispersion relation explain the morphology?"
- Independent parameter: Use spectral index (β) from common signal
- Test: Does Lambda dispersion equation predict observed β?

**Only if independent prediction matches observation → Lambda constraint**

---

## Expected Timeline

| Week | Task | Deliverable |
|------|------|-------------|
| 1 | Synthetic validation | Test suite passes; false positive rate ≈ 5% |
| 2 | NANOGrav catalog scan (70 PSRs) | Table: anomalies per pulsar |
| 3 | Ensemble analysis | PC1 p-value, cross-correlations, sky test |
| 4 | Cross-array replication | NANOGrav vs IPTA comparison |
| 5+ | (If needed) Characterization & Lambda test | Spectral index, Lambda constraint |

---

## Key Files

```
/home/claude/
├── RESIDUAL_FIRST_DESIGN.md              # Framework specification
├── residual_first_common.py              # Common statistical module
│
├── PULSAR_TIMING_AUDIT_GUIDE.md          # This file
├── nanograv_loader.py                    # NANOGrav data access
├── pulsar_timing_residual_audit.py       # Single-pulsar audit
└── common_residual_search.py             # Ensemble & replication
```

---

## Historical Context

This audit strategy follows the scientific path that led to the GW background discovery:

```
2021: NANOGrav (15 yr dataset) finds common red noise signal
       ↓
2023: IPTA independently confirms signal
       ↓
2023: NANOGRAV + IPTA + CPTA coordinate detection announcement
       ↓
2024: Detailed characterization (spectral index, sky correlation, etc.)
       ↓
202X: Planned: Cross-frequency constraints (Lambda, modified gravity, etc.)
```

Your audit is **independent verification** at the 2021–2023 stage, before interpretations.

**Why this matters scientifically:**
- Reproduces discovery with blind methodology
- Validates NANOGrav's finding from first principles
- Identifies any systematic errors in their approach
- Provides independent basis for Lambda testing (if appropriate)

---

## Next Concrete Step

Run synthetic test suite:

```bash
python /home/claude/pulsar_timing_residual_audit.py
python /home/claude/common_residual_search.py
```

Both should complete without errors and show:
- `✓ No anomalies (H_0 not rejected)` for individual pulsars
- PC1 p-value > 0.05 for ensemble (no common signal in noise)

Once synthetic tests pass, acquire real NANOGrav DR15 data and repeat.

---

**End of Guide.**

Questions? Check RESIDUAL_FIRST_DESIGN.md for methodology details.

