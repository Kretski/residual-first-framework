# Shammass R0 — post-gate design calculation v2.1

**Status:** post-gate design analysis. Does not change R0. R0 verdict stays `NOT_TESTABLE` for both series.

**Files:** `r0_design_calc_v2_1.py` (SHA-256 `e8d82811dc770fe63d8ef9ca207e620b42fecf40051ff0ab0880f20dad234315`), `r0_design_results_v2_1.json` (SHA-256 `ff2f7f7086699ea87e8102b72f1c5b0c861278befd3789b1ad3f4126418a120d`). Seed 20261125.

**Inputs read:** the two GPE curve CSVs (k, f) and the `k_um_inv` column of the two experimental CSVs from `realdata/shammass2012/unified_v1/exports/`. No experimental frequency is read and no residual is formed.

## What it computes

1. **R0 rule:** full partial F with candidate-dependent denominator, Holm at 0.01 over q = 1, 2, 3, q* = smallest raw p among rejected. Full confusion matrix including "no detection".
2. **R1 (separate question):** one pre-declared test of the F·k² template at α = 0.01. Analytic m for 80% power and simulated power for true q = 1, 2, 3.
3. **R1 sensitivity to misspecification:** analysis assumes independent, equal-SEM errors; data generated with AR(1) correlation (ρ = 0.3, 0.6), SEM rising 1→3 across k, and log-normal SEM (σ = 0.5) both as a single fixed profile and as an ensemble with a new profile per replicate.

Every simulated rate carries its count, n and exact 95% Clopper–Pearson interval.

## Conditions that must be stated with any use of these numbers

- Scenarios `raw_*` assume the exact GPE curve and per-point frequency SEM only. `raw_N32_idealized` also assumes 32 equally spaced, independent, equally precise points over the same k range. They are mathematical benchmarks, not forecasts of an experiment.
- The AR(1) and SEM scenarios are sensitivity cases for chosen models, not estimates of the real correlation or SEM structure.
- The single fixed log-normal SEM profile is **conditional on one realization**. The ensemble version characterizes the distribution; the two differ markedly for the small series (FPR 0.04% vs 0.89%).
- m is a dimensionless amplitude in units of the per-point noise. It is not converted to a physical Λ or ΔΛ.
- A positive R1 result would mean "an excess with a significant projection on F·k²". It does not identify the generating form: R1 has higher power for a true q = 3 excess than for q = 2.

## Main findings

- Shape identification under the R0 rule needs m ≈ 64 (16 points, exact data) or 32–64 (32 points) for about 80–90% correct q* = 2; at the detection threshold it is near chance.
- R1 m for 80% power (large / small): gate-like 17.0 / 26.3; exact data, observed points 9.6 / 15.9; 32 idealized points 6.5 / 7.2.
- R1 false-positive rate under independent errors ≈ 1%. Under AR(1) correlation it rises to about 3–5% (ρ = 0.3) and 6–17% (ρ = 0.6); under SEM rising across k to about 2–4%. A future R1 protocol therefore needs per-point SEM and a covariance model or bound fixed in advance.
