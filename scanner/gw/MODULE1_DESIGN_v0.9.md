# Module 1 — Residual-first test of Λ-type dispersion in gravitational-wave data

**Status: DRAFT v0.9.1 — not frozen.** (v0.9.1: §2 coefficient and sign fixed by F4.) Nothing in this document has been run on
real gravitational-wave events. Version 1.0 will be frozen (commit, tag and hash)
only after the feasibility checks in §12 are complete. The first run after the
freeze is the calibration suite (§7–§9), never the real events.

## 0. Purpose and scope

This is **not** an attempt to refute the LIGO–Virgo–KAGRA (LVK) results. LVK have
tested the same class of modified dispersion with a full Bayesian framework and
found no evidence for it (GWTC-3 tests of General Relativity; GWTC-4.0 parameterized
tests). Module 1 is an **independent, residual/template-based cross-check** of the
same physical class, with a pre-calibrated estimator of our own. The realistic
expected outcome is an independent bound on Λ, weaker than the LVK bound.

## 1. Relation to Module 0 and the pending basin validation

Module 0 (`scanner/README.md`, tag `v0-module0-closed`) validated the scanner
logic on real data: positive sign recovered quantitatively (piano), negative sign
recovered qualitatively but not quantitatively (GLOBEX A3, amplitude distorted by
nonlinear shoaling), and two real-data confounds identified (Black Sea).

A clean, linear-regime negative-sign test (flume or basin data with low Ursell
number and several simultaneously recorded probes at known spacing) is still
pending; data have been requested (LEGI Grenoble; BSHC Varna). **Decision:** Module 1
preparation proceeds in parallel; when suitable basin data arrive they are analysed
under their own pre-registration, and their result is recorded here. The freeze of
v1.0 records the status of the basin validation (completed, or pending as a stated
limitation). The GW-specific parts of the pipeline (noise weighting, orthogonalization
against waveform parameters, detector calibration) cannot be validated by basin data
and are validated in §7–§9 in any case.

## 2. Hypothesis and signal template

Model: ω² = c²k²(1 + Λk²), the α = 4 member of the modified-dispersion class
E² = p²c² + A_α p^α c^α used by LVK.

For a compact-binary signal the dispersion adds a frequency-domain phase
(verified by F4, see below):

    δΨ(f) = sign(Λ) · (4π³/3) · |Λ| · I₄(z) · f³ / c³ ,
    I₄(z) = (c / H₀) ∫₀^z (1 + z')² / E(z') dz' ,

in the LALSimulation Fourier convention (positive Λ gives a positive f³
coefficient). The factor 1/3 is the 1/(α − 1) of the Mirshekari–Yunes–Will form.
Conversion to the LVK parameterization:

    Λ = A₄ (ħc)² = (λ_A / 2π)² ,    Λ [m²] = A₄ [eV⁻²] × 3.894·10⁻¹⁴ .

**F4 result (lal 7.7.1, lalsimulation 6.2.1, IMRPhenomXPHM, D_L = 400 Mpc,
λ_eff = 10⁻³·⁵ m):** the LALSimulation dispersion term is a pure f³ phase
(fitted f² coefficient ~10⁻²⁰, residual 7·10⁻¹⁶ rad), changes sign with sign(A),
and equals the expression above to a ratio of 0.9983 (Planck15 and Planck18 give
the same z to 10⁻³); the expression without the 1/3 factor is off by a factor of 3.
LALSimulation's λ_eff equals λ_A. The residual 0.17 % is attributed to the
cosmology used internally by LALSimulation and is negligible here; it is recorded.
In Module 1, Λ is always the coefficient defined through ω² = c²k²(1 + Λk²); any
comparison with other normalizations of the Λ model goes through this relation.

To first order in δΨ the residual after subtracting the GR waveform is

    r(f) ≈ i δΨ(f) h_GR(f) = Λ · T₃(f),   T_p(f) = i · (4π³/3) · I₄(z) / c³ · f^p · h_GR(f).

**Primary template:** p = 3. **False-shape controls:** p = 2 and p = 4.
p = 0 and p = 1 are fully degenerate with coalescence phase and time and are not tested.
The linear approximation requires |δΨ| ≪ 1 rad over the band; this is checked for
every injection amplitude.

## 3. Cohorts

- **Discovery cohort:** GWTC-1, GWTC-2.1 (O3a) and GWTC-3 (O3b).
- **Confirmation cohort:** GWTC-4.0 (O4a), fixed in advance. Later data (e.g. O4b)
  are not added after results are seen; they may form a separate, later blind extension.
- **Selection rule (applied from catalog tables only, before any residual is computed):**
  binary black holes (both component masses above 3 M☉, catalog median values) with
  median network matched-filter SNR ≥ 12. The resulting event list, with GWOSC
  identifiers and file checksums, is committed before the first run.

## 4. Data and preprocessing (fixed)

- Strain: GWOSC event data, the detectors, segment duration, sampling rate and
  frequency range used in the corresponding LVK parameter-estimation (PE) run.
- Noise power spectral densities: those released with the PE results, not re-estimated.
- Events with glitch subtraction: the LVK-released cleaned frames are used, as in PE.
- Inner product: ⟨a, b⟩ = Σ_detectors 4 Re ∫ a*(f) b(f) / S_n(f) df over the PE band.
- Any deviation from these rules for a specific event is decided and recorded before
  its residual is computed.

## 5. Waveform models and reference parameters

- Two models per event, each evaluated **with its own PE samples** (parameters
  from one model must not be used to generate the other; otherwise the residual
  contains model mismatch).
- Discovery cohort: IMRPhenomXPHM and SEOBNRv4PHM (the effective-one-body model
  released with GWTC-2.1/3 PE). Confirmation cohort: IMRPhenomXPHM and SEOBNRv5PHM.
  **Availability of these PE sample sets per event is verified in F5;** the
  pairing is fixed in v1.0.
- Reference point: the maximum-likelihood sample of each model's PE release.
- The same event is always analysed with both models by the same scanner; the model
  is never chosen according to the result.

## 6. Orthogonalization and estimator

For each event, detector network and model:

1. Tangent basis: numerical derivatives ∂h/∂θ_j at the reference point for all
   parameters sampled in PE (coalescence time and phase, masses, spins, distance,
   inclination, polarization, sky position), plus the detector-calibration
   directions (amplitude and phase spline nodes) as nuisance directions.
2. T_p^⊥ = T_p − P_tangent T_p, the part of the template that cannot be
   reproduced by a change of standard GR parameters or of the calibration.
3. Per-event estimate: Λ̂_i = ⟨r, T₃^⊥⟩ / ⟨T₃^⊥, T₃^⊥⟩, with Fisher uncertainty
   σ_i = ⟨T₃^⊥, T₃^⊥⟩^(−1/2); the empirical σ_i from §7 replaces it if they differ.
4. The same for p = 2 and p = 4 (control amplitudes).
5. The fraction of T₃ that survives orthogonalization, ‖T₃^⊥‖/‖T₃‖, is reported per event.

## 7. Null calibration (real detector noise)

For every event: off-source segments of real noise from the same detectors and the
same observing period (excluding times near known candidates and data-quality
vetoed times), with the event's own GR reference waveform injected. The full
pipeline is run and the distribution of Λ̂ under Λ = 0 is recorded: bias, scatter,
empirical σ_i, and the per-event and catalog-level false-positive rates. The number
of off-source injections per event is fixed in v1.0 from the compute budget (F6).

## 8. Injection ladder

| Injection | Checks |
|---|---|
| GR, Λ = 0 | null calibration (§7) |
| f³ at pre-fixed amplitudes of both signs (in units of the catalog σ) | recovery of sign and amplitude |
| f² at matched amplitude | wrong-shape rejection |
| f⁴ at matched amplitude | wrong-shape rejection |

Amplitudes are fixed in v1.0 and not changed after the first results.

## 9. Detector-calibration systematics

GR injections multiplied by calibration perturbations drawn from the published
LVK calibration-uncertainty models (O1–O3; O4), without dispersion. Criterion:
calibration-only injections must not produce a catalog-level f³ signal with a
consistent sign and the K(z) dependence of the Λ hypothesis.

## 10. Catalog-level tests

- Combined estimate: Λ̂ = Σ w_i Λ̂_i / Σ w_i, w_i = 1/σ_i² (K(z) is inside the template,
  so Λ is common to all events).
- **Sign consistency:** a real effect gives Λ̂_i of one sign across events.
- **Heterogeneity:** Q = Σ (Λ̂_i − Λ̂)²/σ_i², calibrated by the null suite; a real
  common Λ gives Q consistent with the null, systematics tend to give excess Q.
- **Permutation control:** K(z) is shuffled among events while the raw per-event
  f³ amplitudes are kept; this destroys the propagation relation but keeps waveform
  systematics, noise, SNR distribution and uncertainties. The observed catalog
  statistic is compared with the permutation distribution.
- **Waveform robustness:** the result must be consistent between the two models.

## 11. Decision rules (fixed in v1.0)

1. **Null failure:** GR injections give a significant systematic Λ̂ → calibration
   failure; real events are not interpreted.
2. **Wrong-shape failure:** f² or f⁴ injections are systematically recovered as f³
   → shape discrimination failure; no physical interpretation.
3. **Calibration failure:** calibration-only injections produce a catalog-level f³
   signal with the Λ sign and scaling → calibration confound unresolved.
4. **All controls pass:**
   - no significant f³ → bound on Λ (reported with the conversion to A₄ from §2);
   - significant f³ → **candidate only**, which must pass the confirmation cohort,
     the waveform and calibration controls and the permutation test; a p-value
     alone is not called a discovery.

Every result, including failures, is reported.

## 12. Feasibility checks before v1.0

- **F1** DONE: lalsuite (pip) with lal 7.7.1 / lalsimulation 6.2.1, pyseobnr 0.3.7,
  gwpy 4.0.2, Python 3.11 in WSL (`requirements_gw2.txt`, `environment_gw2.yml`).
  The conda-forge lal/lalsimulation builds showed a SWIG type mismatch.
- **F2** DONE: one waveform: IMRPhenomXPHM 0.014 s, SEOBNRv5PHM 0.1 s (first call
  ~9 s), SEOBNRv4PHM 2.3 s.
- **F3** Numerical derivatives: step sizes and stability for all parameters.
- **F4** DONE (see §2): coefficient (4π³/3), sign convention, λ_eff = λ_A,
  Λ = A₄(ħc)²; numerical agreement 0.9983 with LALSimulation.
- **F5** Verify, per event, the availability of PE samples for both models and of
  the released PSDs and calibration envelopes.
- **F6** Rough budget DONE: tangent basis per event 0.4 / 3.3 / 72 s
  (XPHM / v5PHM / v4PHM); the full calibration suite is feasible on a laptop.
  Final numbers of injections are fixed in v1.0. Compute budget: estimate the total cost of §7–§9 and fix the numbers of
  off-source injections and ladder realizations; if needed, run the large
  calibration suite with IMRPhenomXPHM only and a smaller pre-fixed control suite
  with the effective-one-body model.

## 13. Provenance

Software versions, data-file checksums, the event list, the configuration hash,
and git tags are recorded for every step. Code, output and documentation are in
English.
