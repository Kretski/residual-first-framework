# Module 1 — Residual-first test of Λ-type dispersion in gravitational-wave data

**Status: DRAFT v0.9.6 — not frozen.** (v0.9.6: residual chain verified on all
events; reference-fit quality gate; extrinsic refinement window; EOB generation rules.)
Previous: **DRAFT v0.9.5.** (v0.9.5: strain data, missing-data rule,
per-detector consistency check, primary counts after the detector rule.)
Previous: **DRAFT v0.9.4.** (v0.9.4: tangent subspace from posterior
samples (F3), subspace systematic in the calibration, exclusions following LVK.) (v0.9.1: §2 coefficient checked by F4;
v0.9.2: event list frozen, F5 results and rules in §3 and §5; v0.9.3: §2 CORRECTED to the
group-velocity phase used by LVK from GWTC-4.0 on — the v0.9.1 conclusion was wrong.) Nothing in this document has been run on
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

For a compact-binary signal the dispersion adds a frequency-domain phase. Two
prescriptions exist (GWTC-4.0 Tests of GR II, arXiv:2603.19020, §3.1):

- particle velocity (LVK up to GWTC-3; Mirshekari–Yunes–Will; implemented in
  LALSimulation's LIV option):  δΨ = +(4π³/3) Λ I₄(z) f³ / c³
- **group velocity (LVK from GWTC-4.0; consistent with the WKB treatment, Ezquiaga
  et al. 2022): δΨ = −4π³ Λ I₄(z) f³ / c³**

with I₄(z) = (c/H₀) ∫₀^z (1 + z')² / E(z') dz' and the same Fourier convention
(h̃ = h̃_GR e^{iδΨ}). The ratio is (1 − α) = −3: the particle velocity pc²/E =
c(1 − A₄p²c²/2) and the group velocity dE/dp = c(1 + 3A₄p²c²/2) differ in sign and
by a factor 3. **Module 1 uses the group-velocity form.** Injections must add this
phase explicitly; LALSimulation's built-in LIV option must NOT be used (it implements
the particle-velocity form).

Conversion: Λ = A₄ (ħc)² = (λ_A / 2π)², Λ [m²] = A₄ [eV⁻²] × 3.894·10⁻¹⁴.
LVK GWTC-4.0 combined bound (83 events, group velocity, Table 5): A₄ ∈ [−620, +190]
eV⁻² (90%), i.e. Λ ∈ [−2.4·10⁻¹¹, +7.4·10⁻¹²] m².

**F4 record (corrected):** F4 verified numerically that LALSimulation's LIV term is a
pure f³ phase equal to +(4π³/3) Λ I₄ f³/c³ (ratio 0.9983). v0.9.1 concluded from this
that the expression without the 1/3 factor was "wrong by a factor 3"; that conclusion
was incorrect — the 1/3 belongs to the particle-velocity prescription, and the
group-velocity form without it (and with the opposite sign) is the physically
appropriate one. In Module 1, Λ is always defined through ω² = c²k²(1 + Λk²).

To first order in δΨ the residual after subtracting the GR waveform is

    r(f) ≈ i δΨ(f) h_GR(f) = Λ · T₃(f),   T_p(f) = −i · 4π³ · I₄(z) / c³ · f^p · h_GR(f).

**Primary template:** p = 3. **False-shape controls:** p = 2 and p = 4.
p = 0 and p = 1 are fully degenerate with coalescence phase and time and are not tested.
The linear approximation requires |δΨ| ≪ 1 rad over the band; this is checked for
every injection amplitude.

## 3. Cohorts

- **Discovery cohort:** GWTC-1, GWTC-2.1 (O3a) and GWTC-3 (O3b).
- **Confirmation cohort:** GWTC-4.0 (O4a), fixed in advance. Later data (e.g. O4b)
  are not added after results are seen; they may form a separate, later blind extension.
- **Selection rule (applied from catalog tables only, before any residual is computed):**
  GWOSC catalogs GWTC-2.1-confident (includes the O1/O2 events reanalysed with
  IMRPhenomXPHM and SEOBNRv4PHM) and GWTC-3-confident for discovery, GWTC-4.0 for
  confirmation; mass_2_source (catalog median, source frame) > 3 M☉;
  network_matched_filter_snr (catalog PE median) ≥ 12.
- **Frozen list (tag `v1-module1-eventlist`, `scanner/gw/event_selection/`):**
  discovery 37 events, confirmation 24 events; raw catalog snapshots with SHA-256 and
  access time. 43 GWTC-4.0 candidates without public PE are excluded (no reference
  point); see `NOTES.md`.
- **F5 eligibility (`scanner/gw/f5/`):**
  discovery 32 BOTH_OK, 4 SINGLE_MODEL (GW151226, GW190521, GW190602_175927,
  GW190828_063405: no SEOBNRv4PHM in the public release), 1 excluded
  (GW191204_171526: no calibration envelope in the public file, for either model);
  confirmation 24 BOTH_OK.
- **Further exclusions from the primary test (decided in v0.9.4, before any residual):**
  - GW231123_135430: strong waveform-model systematics; LVK exclude it from their
    combined dispersion bound (arXiv:2603.19020, Appendix D: IMRPhenomXPHM and
    NRSur7dq4 disagree strongly). It also showed the least stable T₃⊥ in F3.
  - Events recorded by fewer than two detectors (GW230814_230901), as in the LVK
    selection, so that the primary sample is comparable with the LVK one. The
    detector list of every event is taken from its PE file (F5 output) and the rule
    is applied mechanically to both cohorts.
- **Primary catalog test: 31 discovery + 20 confirmation events** (v0.9.6: GW200225_060421
  removed by the reference-fit gate of §5).
- Earlier count in v0.9.5: **32 discovery + 20 confirmation events.** The detector rule
  removes three confirmation events whose PE used a single detector
  (GW230814_230901: L1; GW231231_154016: H1; GW240104_164932: H1); none of the
  three is in the LVK tests-of-GR selection either. SINGLE_MODEL events and the two excluded events are analysed
  separately with the same pipeline and reported, but are not part of the primary test.

## 4. Data and preprocessing (fixed)

- Strain: GWOSC event data, the detectors, segment duration, sampling rate and
  frequency range used in the corresponding LVK parameter-estimation (PE) run.
- Noise power spectral densities: those released with the PE results, not re-estimated.
- Events with glitch subtraction: the LVK-released cleaned frames are used, as in PE.
- Inner product: ⟨a, b⟩ = Σ_detectors 4 Re ∫ a*(f) b(f) / S_n(f) df over the PE band.
- Strain files (`scanner/gw/strain_manifest.csv`, SHA-256 per file): GWOSC 4096 s
  files at 4096 Hz for every detector used in the event's PE (134 files, 16.7 GB).
  One file provides both the on-source segment and the off-source segments of the
  null suite.
- Missing-data rule: if the standard GWOSC file is unavailable for a detector, the
  event-specific file of the catalog release is used; if none exists, the event is
  analysed with the remaining detectors and the detector rule of §3 applies.
  Applied once: GW170608 H1 (not in the bulk archive because H1 was not in nominal
  observing mode) is taken from the GWTC-1-confident v3 event release (4096 s,
  centred on the event). Its null suite uses off-source segments of the same file,
  and its null distribution is reported separately.
- Any deviation from these rules for a specific event is decided and recorded before
  its residual is computed.

## 5. Waveform models and reference parameters

- Two models per event, each evaluated **with its own PE samples** (parameters
  from one model must not be used to generate the other; otherwise the residual
  contains model mismatch).
- Discovery cohort: IMRPhenomXPHM and SEOBNRv4PHM (the effective-one-body model
  released with GWTC-2.1/3 PE). Confirmation cohort: IMRPhenomXPHM and SEOBNRv5PHM.
  Availability verified in F5 (see §3).
- Labels: discovery `C01:IMRPhenomXPHM`, `C01:SEOBNRv4PHM`; confirmation
  `C00:IMRPhenomXPHM-SpinTaylor`, `C00:SEOBNRv5PHM` (all 24 confirmation events use
  exactly these). `Mixed` labels are never used.
- Waveform settings are taken **per event and per label** from the PE configuration:
  reference frequency (20 Hz for most events, 10 Hz e.g. for GW231123_135430),
  analysis band, waveform starting frequency (e.g. `waveform: 10.0`, or LALInference
  `fmin-template`), sampling rate, segment duration, and the waveform-argument
  dictionary (e.g. `PhenomXPrecVersion: 320` for XPHM-SpinTaylor). A waveform
  generated without these settings is not the waveform of the PE analysis and its
  residual would contain the difference.
- Reference point: the maximum-likelihood sample of each model's PE release, followed
  by a local re-optimisation of the extrinsic parameters only (coalescence time,
  phase, distance) with the released PSDs. Reason: bilby analyses marginalise over
  time and distance, so the stored log-likelihood is marginalised and the time and
  distance of a sample are draws, not maxima; LALInference analyses have few samples.
  Intrinsic parameters are not changed.
- Calibration envelope: from the label itself; if absent, from the XPHM label of the
  same event (same detectors and strain); if absent in both, the event is excluded.
- **Extrinsic refinement (verified, `residual_check.py`):** a common time shift within
  ±100 ms (step 0.05 ms) and a common complex amplitude for the network. Bilby-based
  labels need |Δt| < 1 ms; LALInference-based labels (SEOBNRv4PHM) need event-dependent
  shifts up to 39 ms and arbitrary phases, i.e. their stored reference times are not at
  the likelihood maximum; after refinement both models fit equally well.
- **EOB generation:** time-domain waveforms are cropped to the analysis segment
  (merger 2 s before its end) with the same 0.4 s taper as the data; if the model
  refuses the PE sampling rate (ringdown above Nyquist), it is generated at 2× or 4×
  the rate and only the PE frequency grid up to f_high is used. Our TD→FD path agrees
  with IMRPhenomXPHM at the same parameters (match 0.993, Δt −0.4 ms, GW150914).
- **Reference-fit quality gate (fixed before any Λ estimate):** an event–model pair
  whose refined reference gives SNR_mf / catalog SNR outside 0.8–1.2 is moved to the
  separately analysed set for that model; since the primary test needs both models, the
  event leaves the primary test. Result over all 116 pairs (`residual_check.csv`):
  0 errors, 2 pairs outside: GW200225_060421 SEOBNRv4PHM (0.66; our EOB waveform matches
  XPHM only at 0.884 at that point) and GW231123_135430 SEOBNRv5PHM (0.71; event already
  excluded, the same waveform systematics LVK report in their Appendix D).
- The same event is always analysed with both models by the same scanner; the model
  is never chosen according to the result.

## 6. Orthogonalization and estimator

For each event, detector network and model:

1. Tangent subspace from the posterior (F3): finite-difference derivatives were
   found unstable (near-degenerate parameter directions; numerical noise of the EOB
   and XPHM-SpinTaylor models). Instead, two DISJOINT random sets A and B of n = 800
   posterior samples of the label (fixed seed 20261001; all samples if fewer) give
   waveforms h(θ_i); the differences Δ_i = h(θ_i) − h_ref, whitened with the released
   PSD, are decomposed by SVD, and the principal components explaining 99% of the
   variance (rule frac ≥ 0.99) form the subspace. The analytic coalescence-time
   direction −2πif·h_ref and the detector-calibration directions (amplitude and
   phase spline nodes) are always added. **Set A defines the primary subspace;**
   set B is used only for the subspace systematic (§7).
2. T_p^⊥ = T_p − P_A T_p, the part of the template that cannot be
   reproduced by a change of standard GR parameters or of the calibration.
3. Per-event estimate: Λ̂_i = ⟨r, T₃^⊥⟩ / ⟨T₃^⊥, T₃^⊥⟩, with Fisher uncertainty
   σ_i = ⟨T₃^⊥, T₃^⊥⟩^(−1/2); the empirical σ_i from §7 replaces it if they differ.
4. The same for p = 2 and p = 4 (control amplitudes).
5. Reported per event as diagnostics (no thresholds): the survival s = ‖T₃^⊥‖/‖T₃‖
   (it sets the estimator noise, σ ∝ 1/s), the cosine between T₃^⊥ from sets A and B
   (it controls the leakage of reference-point imperfections), and |s_A − s_B|/s_B.
   For a pure f³ signal the estimate is unbiased for any subspace, since
   ⟨T₃, T₃^⊥⟩ = ‖T₃^⊥‖².

## 7. Null calibration (real detector noise)

For every event: off-source segments of real noise from the same detectors and the
same observing period (excluding times near known candidates and data-quality
vetoed times), with the event's own GR reference waveform injected. The full
pipeline is run and the distribution of Λ̂ under Λ = 0 is recorded: bias, scatter,
empirical σ_i, and the per-event and catalog-level false-positive rates. The number
of off-source injections per event is fixed in v1.0 from the compute budget (F6).
Off-source segments are taken on both sides of the event, skipping invalid data (gaps,
file edges); at least 4 valid segments are required per detector.

**Subspace systematic.** Every null and ladder injection is analysed twice, with the
set-A and the set-B subspace. The difference Λ̂_A − Λ̂_B over the null suite measures
the cost of the subspace choice in the units of the estimate. Rule (fixed now): if the
median |Λ̂_A − Λ̂_B| over the null suite of an event exceeds 0.5 σ_i, the event is
moved from the primary test to the separately reported set; otherwise the RMS of
Λ̂_A − Λ̂_B is added in quadrature to σ_i.

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
- **Per-detector consistency:** Λ̂ is also estimated separately for each detector
  (same reference waveform, the detector's own projection and PSD). A propagation
  effect is the same in all detectors (same Λ, same source distance); detector
  calibration errors, glitches and noise are not. Statistic per event: χ² of the
  per-detector estimates about their weighted mean; its null distribution comes from
  the null suite.

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
     the waveform and calibration controls, the permutation test and the
     per-detector consistency test (catalog-level χ² not in the upper 1% tail of its
     null distribution); a p-value alone is not called a discovery.

Every result, including failures, is reported.

## 12. Feasibility checks before v1.0

- **F1** DONE: lalsuite (pip) with lal 7.7.1 / lalsimulation 6.2.1, pyseobnr 0.3.7,
  gwpy 4.0.2, Python 3.11 in WSL (`requirements_gw2.txt`, `environment_gw2.yml`).
  The conda-forge lal/lalsimulation builds showed a SWIG type mismatch.
- **F2** DONE: one waveform: IMRPhenomXPHM 0.014 s, SEOBNRv5PHM 0.1 s (first call
  ~9 s), SEOBNRv4PHM 2.3 s.
- **F3** DONE (`scanner/gw/f3*`): finite-difference derivatives unstable even with
  per-parameter step choice (XPHM-SpinTaylor and EOB numerical noise; nearly
  degenerate directions). Posterior-sample PCA (n = 200/400/800 per set, 4 events ×
  2 models): survival of f³ 0.07–0.53; the A/B direction agreement is best for
  frac ≥ 0.99 (cos up to 0.999) and improves with n for XPHM; higher thresholds pick
  up model noise (up to ~380 components for SEOBNRv5PHM, GW231123). Rule frac ≥ 0.99,
  n = 800, fixed seed adopted; residual instability handled as a measured systematic.
- **F4** DONE, corrected in v0.9.3 (see §2): LALSimulation implements the
  particle-velocity phase (+4π³/3, agreement 0.9983); Module 1 uses the group-velocity
  phase (−4π³) as LVK do from GWTC-4.0; λ_eff = λ_A; Λ = A₄(ħc)².
- **F5** DONE (see §3 and §5): all 61 PE files downloaded and MD5-verified (8.7 GB);
  content checked per label (samples with log-likelihood, parameters, PSDs,
  calibration envelopes, configuration). Two checker defects found and fixed before
  any residual: remote reading hit Zenodo rate limits (replaced by verified local
  files), and configuration keys occur with '_' or '-' (both accepted).
- **Residual chain** DONE (v0.9.6): see §5. Diagnostics note: with a least-squares
  complex amplitude the χ² drop equals SNR_mf² identically, so SNR_mf is compared with
  the catalog SNR instead; the reduced χ² of the residual is a weak test (about 7000
  degrees of freedom in the band), reported only as a diagnostic.
- **Version policy for extensions:** later data (O4b, GWTC-5.0) are analysed only with
  the code at the v1.0 tag, unchanged; any change makes a new version with its own
  registration before the new data are looked at.
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
