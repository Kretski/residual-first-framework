# Pre-registration: independent replication of the V2.1 HD test on EPTA DR2

Author: Dimitar Kretski
Status: **to be committed and tagged BEFORE any EPTA residual is computed or inspected.**
Frozen analysis code: `common_residual_search_v2.py` (V2.1), SHA-256 recorded at tag time:
`<fill in: certutil -hashfile common_residual_search_v2.py SHA256>`

## 1. Background (what was already seen)

On NANOGrav 15-yr (v2.1.0, epoch-averaged, 68 pulsars, 30-d bins) the V2.1 pipeline gave:

| test | p (NANOGrav) | status |
|---|---|---|
| cross_correlation (z²) | 0.045 | **invalid** — FPR 0.69 on real sampling (see diagnostics_shared_structure.py) |
| monopole | 0.035 | calibrated |
| hd_sky_scramble | 0.034 | calibrated; leaks under dipole (0.60) |
| hd_perp_mono_dipole:scramble | 0.046 | calibrated (0.04); robust to monopole (0.05); partial dipole leakage (0.10–0.20) |

The last test was added after seeing V2 results (exploratory on NANOGrav).
Its choice as the primary EPTA test is therefore informed by NANOGrav, but EPTA
data have NOT been seen.

## 2. Data

EPTA DR2, doi:10.5281/zenodo.8300645, file EPTA-DR2.zip, MD5 307b2b99cd352876b07409899448af8b.
Primary dataset: the full 25-pulsar release as distributed. If the release contains
separate variants (e.g. "DR2full" and "DR2new"), DR2full is primary and DR2new is
secondary; both are reported.

Residuals: computed with PINT from the released .par/.tim files, using the released
timing model unchanged (no re-fit beyond what PINT needs to evaluate residuals).
Pulsars that PINT cannot load are listed and excluded; no other exclusions.

## 3. Analysis (frozen)

`common_residual_search_v2.py` V2.1 with defaults: 30-d bins, min_overlap = 12,
mirror-shift temporal null, quadratic detrend, n_null = n_scramble = 2000, seed = 0.
No `--min-span`, no other options.

**Primary test (single, no multiplicity correction):**
`hd_perp_mono_dipole:scramble`, one-sided, α = 0.05.

Secondary (reported, not used for the decision): all other V2.1 tests.
`cross_correlation` is reported but declared uninterpretable in advance.

## 4. Calibration gate (run BEFORE the real-data run)

`--validate --real-gaps` on EPTA sampling/errors/positions, H₀ only, 200 trials.
If the H₀ rate of the primary test is > 0.09 (≈ 3σ above 0.05 for 200 trials),
the real-data run is NOT performed and the calibration failure is reported instead.

## 5. Decision rule and reporting

- p ≤ 0.05: "EPTA DR2 shows an HD-like correlation beyond monopole + linear dipole,
  consistent with the NANOGrav result" — with the dipole-leakage caveat.
- p > 0.05: "not replicated at α = 0.05". Reported with equal prominence.
- EPTA and NANOGrav share many pulsars and the same GWB realisation; agreement tests
  independence of *systematics*, not of the signal.
- All outputs (JSON provenance, validation logs) are published regardless of outcome.
- No further variants are run on EPTA before this result is published.
