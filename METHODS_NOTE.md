# Residual-level Hellings–Downs diagnostics on NANOGrav 15-yr and EPTA DR2: calibration failures and a pre-registered replication

Dimitar Kretski · 26 September 2026

## Abstract

A simple residual-level test for Hellings–Downs (HD) correlations, calibrated by sky-position scrambling and orthogonalised against monopole and dipole patterns, gives p = 0.0135 on EPTA DR2full and p = 0.0065 on EPTA DR2new in a pre-registered replication. Both calibration gates passed (H₀ rates 0.045 and 0.040).

The more general finding is negative. On real pulsar-timing-array (PTA) sampling, an omnibus cross-correlation statistic has a false-positive rate of 0.69 (NANOGrav) to 1.00 (EPTA), although it is calibrated on idealised simulations (0.065). Temporal-shift nulls also fail for monopole and angular tests (0.10–0.15). Only sky-scramble tests remain calibrated across all three data sets.

We trace the omnibus failure to two mechanisms: time-varying variance shared across pulsars (TOA uncertainties and TOA density) and inaccuracy of the shift null for red-noise-dominated series. The HD result recovers the known PTA signal with a simple method. It is not a new detection.

## 1 Introduction and scope

The aim is methodological: to find which residual-level cross-pulsar statistics stay calibrated on real PTA data, and whether a calibrated one recovers the published HD signal. The major PTAs already report HD evidence with full Bayesian and optimal-statistic analyses that model each pulsar's noise.

Residual-level tests are attractive because they are simple and fast. They need only post-fit timing residuals, sky positions and permutation nulls. The risk is that real sampling, shared instrumental history and red noise violate the assumptions behind those nulls. This note documents where that happens.

We make no claim about new physics. In particular, the results do not bear on modified-gravity or dispersion models, which would need a quantitative prediction of their effect at nanohertz frequencies first.

## 2 Data

Three data sets were analysed with the same frozen analysis code.

| Data set | Pulsars used | Span | Residuals |
| --- | --- | --- | --- |
| NANOGrav 15-yr v2.1.0 | 68 | up to 16 yr | released epoch-averaged, un-whitened, post-fit |
| EPTA DR2full | 24 of 25 | 14–25 yr | PINT, released model re-fitted once (WLS) |
| EPTA DR2new | 25 | ~10 yr | PINT, released model re-fitted once (WLS) |

NANOGrav residuals were read directly from the release tables (`*_nb.avg.res`). The 23 pulsars with an extra whitened-residual column match the 23 with detected red noise in the NANOGrav paper, which confirms the column mapping.

EPTA DR2 provides only .par/.tim files. Residuals were computed with PINT 1.1.7 under WSL with 80-bit long double; native Windows PINT runs at reduced precision. The released models use `UNITS TCB` and `BINARY T2`; PINT's conversions are approximate, so each model was re-fitted once. J0030+0451 was excluded from DR2full because PINT rejected a JUMP with no TOAs.

Reading the EPTA .tim files with tempo2 semantics required four fixes: END ends only its own file; FORMAT is global state in reading order; comment markers such as `C??` must be recognised; and leading blanks must be stripped. Without the last one PINT classifies a FORMAT 1 line as fixed-column Parkes format whenever column 41 holds a decimal point.

## 3 Method (V2.1, frozen)

Residuals are binned in 30-d bins by inverse-variance weighted mean, with no interpolation; empty bins stay empty. A quadratic is removed from each binned series over its own span. For each pulsar pair the correlation ρ is computed only over bins both pulsars occupy (pairs with fewer than 12 common bins are dropped) and converted to a Fisher score:

```latex
z_{ab} = \operatorname{artanh}(\rho_{ab})\,\sqrt{n_{ab}-3}
```

Test statistics are weighted projections of z on pair-space templates: monopole (1), a dipole trace (cos ζ) and the HD curve. The orthogonalised tests first remove from the HD template its projection on the monopole alone, or on the monopole plus the full 6-dimensional dipole subspace spanned by the symmetric products of the two pulsars' unit vectors. That subspace covers correlations (n̂_a·d)(n̂_b·d) for any direction d.

Two nulls are used:

- **Temporal ("shift")**: each pulsar's series is mirror-shifted within its own span (a random circular shift of the even extension), then detrended like the data.
- **Sky scramble**: the data are kept and pulsar positions are permuted. This is exact when the correlations do not depend on position.

The omnibus test is the mean of z² under the temporal null. Tests form a primary family (omnibus, monopole, HD sky-scramble) and an exploratory family (orthogonalised HD, both nulls), each with its own Benjamini–Hochberg correction. Permutation p-values are (b + 1)/(n + 1), one-sided.

## 4 Calibration results

Only the sky-scramble tests are calibrated on real sampling in all three data sets. Rates are the fraction of H₀ simulations with p ≤ 0.05, using real TOA times, uncertainties and positions with synthetic noise (white noise plus 1 µs, β = 3 red noise on one pulsar in three).

| Test | Null | NANOGrav (200) | EPTA DR2full (200) | EPTA DR2new (200) |
| --- | --- | --- | --- | --- |
| omnibus z² | shift | 0.69 | 1.000 | 1.000 |
| monopole | shift | 0.07 | 0.135 | 0.150 |
| HD ⟂ monopole | shift | 0.07 | 0.120 | 0.105 |
| HD ⟂ monopole + dipole | shift | 0.04 | 0.120 | 0.115 |
| HD | sky scramble | 0.04 | 0.035 | 0.060 |
| HD ⟂ monopole | sky scramble | 0.04 | 0.045 | 0.065 |
| HD ⟂ monopole + dipole | sky scramble | 0.04 | 0.045 | 0.040 |

The omnibus test looked calibrated on idealised arrays (independent random TOAs, constant errors): mirror shift 0.065, versus 0.14 for a plain circular shift and 0.00 for phase-randomised surrogates. Its failure on real data has two sources, isolated with synthetic arrays (40 trials each) and real-sampling ablations (NANOGrav, 100 trials):

| Scenario | Omnibus H₀ rate |
| --- | --- |
| synthetic: independent gaps, constant errors | 0.03 |
| synthetic: shared observing gaps | 0.07 |
| synthetic: shared 5× error step | 1.00 |
| synthetic: shared high-cadence campaigns, constant errors | 1.00 |
| synthetic: same, thinned to one TOA per bin | 0.10 |
| synthetic: red-noise dominated (errors 0.05–0.2 µs), no shared structure | 0.25–0.28 |
| real sampling, real errors | 0.69 |
| real sampling, constant per-pulsar errors | 0.43 |
| real sampling, constant errors, one TOA per bin | 0.31 |

First, variance that changes at the same dates for many pulsars, through backend upgrades or campaign cadence, is destroyed by independent per-pulsar shifts, so the null z² is too small. Second, for short red-noise-dominated series the shift null is itself inaccurate. Dividing residuals by bin uncertainty did not help (0.98).

Nuisance leakage was measured on NANOGrav geometry (20 trials per injection). The shift-null HD tests fired at 0.25–0.30 under a 0.5 µs monopole and 0.35–0.75 under a 1 µs dipole. The sky-scramble test orthogonalised to monopole + dipole fired at 0.05 and 0.10–0.20.

## 5 Pre-registered EPTA replication

The single primary test, HD ⟂ monopole + dipole with the sky-scramble null, gave p ≤ 0.05 on both EPTA variants. It was chosen after seeing NANOGrav but before any EPTA residual.

Protocol. The analysis code was frozen and tagged (`v2.1-epta-prereg`) before EPTA data were downloaded. The pre-registration fixed the data variant (DR2full primary, DR2new secondary), the primary test at α = 0.05 without multiplicity correction, and a calibration gate: the real run was allowed only if the primary test's H₀ rate on EPTA sampling was ≤ 0.09. Technical amendments (TOA reading, precision, re-fit decision) were committed before the runs they affect. One gate run, corrupted by filesystem errors, is recorded as invalid.

| Variant | Gate H₀ rate | Pulsars | Primary p | HD sky scramble p | HD ⟂ mono (scramble) p |
| --- | --- | --- | --- | --- | --- |
| DR2full (primary) | 0.045 | 24 | **0.0135** | 0.0095 | 0.0080 |
| DR2new (secondary) | 0.040 | 25 | **0.0065** | 0.011 | 0.0095 |
| NANOGrav (exploratory, for comparison) | 0.04 | 68 | 0.046 | 0.034 | 0.026 |

Shift-null and omnibus results are not reported as evidence: they failed calibration on EPTA sampling. The DR2full monopole projection is negative (−1.32), which argues against a clock-type origin. The signal is stronger in the cleaner DR2new subset (reduced χ² 0.7–5.8 versus up to 7208 in DR2full), which argues against the DR2full processing problems as its source. The three results share pulsars and the same background realisation, and DR2new is a subset of DR2full, so no combined p-value is given.

Method check M1 measured leakage in EPTA geometry after the real runs, with a criterion fixed beforehand (rate ≤ 0.10 under both dipole amplitudes and the monopole). Rates of p ≤ 0.05 (100 trials; 50 for monopole and HD):

| Injection | DR2full: primary / HD sky / HD ⟂ mono | DR2new: primary / HD sky / HD ⟂ mono |
| --- | --- | --- |
| none (H₀) | 0.040 / 0.040 / 0.050 | 0.050 / 0.040 / 0.050 |
| dipole 0.5 µs | 0.060 / 0.090 / 0.070 | 0.020 / 0.110 / 0.080 |
| dipole 1.0 µs | 0.010 / 0.100 / 0.110 | 0.000 / 0.160 / 0.160 |
| monopole 0.5 µs | 0.000 / 0.020 / 0.020 | 0.020 / 0.020 / 0.040 |
| HD 1.0 µs (power) | 0.28 / 0.28 / 0.32 | 0.52 / 0.48 / 0.50 |

The primary test met the criterion on both variants. Removing the dipole subspace is what matters: without it, leakage under a 1 µs dipole reaches 0.10–0.16. On NANOGrav geometry the same test still leaked at 0.10–0.20, so robustness depends on the array.

## 6 Limitations

- **No noise modelling.** DM, chromatic and solar-wind noise are not modelled, and EFAC/EQUAD from the EPTA noise files are not applied. Several DR2full pulsars have reduced χ² of 124–7208; DR2new pulsars reach 5.8.
- **Approximate model conversion.** PINT's TCB→TDB conversion is approximate; one WLS re-fit per pulsar was used instead of tempo2.
- **Clock coverage.** Some Effelsberg and LEAP TOAs fall outside the shipped clock files.
- **Modest power.** The primary test detects a 1 µs HD process in 28–52% of EPTA simulations; it is far less sensitive than the collaborations' analyses.
- **Stress-test amplitudes.** Injected dipole and monopole amplitudes were chosen for testing, not estimated from real ephemeris or clock errors. Leakage estimates rest on 20–100 trials.
- **Garden of forking paths on NANOGrav.** The orthogonalised sky-scramble test was introduced after seeing NANOGrav results; only the EPTA runs are pre-registered.
- **Not independent.** NANOGrav and EPTA share pulsars and the same gravitational-wave background realisation; DR2new is a subset of DR2full.

## 7 Conclusions and recommendations

A calibrated, nuisance-orthogonalised sky-scramble statistic recovers HD-like correlations in NANOGrav 15-yr (p = 0.046, exploratory) and in a pre-registered EPTA DR2 replication (p = 0.0135 primary, 0.0065 secondary). This is consistent with the published PTA evidence and adds no new physics.

For anyone building residual-level PTA diagnostics:

1. Validate on the real TOA times, uncertainties and positions of the array, not on idealised simulations. The omnibus test passed idealised calibration and failed on real data at rates of 0.69–1.00.
2. Do not use per-pulsar temporal shifts as the null for statistics that sum many pairs. They destroy shared time-varying variance and misrepresent red-noise-dominated series.
3. Prefer sky-scramble nulls for angular tests. They keep every temporal property of the data.
4. Orthogonalise the HD template against the monopole and the full 6-dimensional dipole subspace. In EPTA geometry this cut dipole leakage from 0.10–0.16 to 0.00–0.06.
5. Measure nuisance leakage for each array's geometry; it differed between NANOGrav and EPTA.
6. Freeze code and gate each real-data run on a pre-specified calibration check.

## Reproducibility and data

All code, the pre-registration with its amendments and results, gate and provenance JSON files, and run logs are in [github.com/Kretski/residual-first-framework](https://github.com/Kretski/residual-first-framework).

| Item | Value |
| --- | --- |
| Frozen analysis | `common_residual_search_v2.py`, tag `v2.1-epta-prereg` (commit ff05b3e) |
| SHA-256 (LF line endings) | 089a94e073c962eb332779e0cc50807ebbb30a6e3f02bb49bced0226497ec2b2 |
| Pre-registration | `PREREGISTRATION_EPTA.md`, first commit e17dca0; Amendment 6 in ca75402 |
| Environment | PINT 1.1.7, Python 3.13, WSL Ubuntu, long double eps 1.08e-19 |
| Seeds | analysis seed 0; gates seed 0; method check M1 seed 1 |

Commands (run from a Linux filesystem):

```bash
python fetch_nanograv15.py
python common_residual_search_v2.py --real
python common_residual_search_v2.py --validate --real-gaps --trials 200
python diagnostics_shared_structure.py
python fetch_epta_dr2.py
python run_epta.py gate  [--variant DR2new]
python run_epta.py real  [--variant DR2new]
python run_epta.py dipole [--variant DR2new]
```

Data and software:

- NANOGrav Collaboration (2025), The NANOGrav 15-Year Data Set v2.1.0, [doi:10.5281/zenodo.16051178](https://doi.org/10.5281/zenodo.16051178); Agazie et al. 2023, ApJL 951, L8 and L9.
- EPTA Collaboration (2023), EPTA DR2, [doi:10.5281/zenodo.8300645](https://doi.org/10.5281/zenodo.8300645); A&A 678, A48 and A50.
- PINT pulsar timing software (Luo et al. 2021, ApJ 911, 45).
