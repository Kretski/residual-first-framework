# Module 1 — Residual-first test of Λ-type dispersion in gravitational-wave data

**Status: v1.0 FROZEN (tag `v1.0-module1`, code SHA-256 22a5127d…, manifest
`FREEZE_v1.0.json`).** Changes after the tag are limited to recording what the official
run finds; the code is not changed. v1.0a: four more search-cohort events excluded under
the existing §3 calibration rule, found at the start of the official build (§3).
Previous: **DRAFT v0.9.15 — ready to freeze.** (v0.9.15: §11 completed — thresholds made
numerical, construction rule one-sided for the union, order of the stages and the opening
of the confirmation cohort fixed.)
Previous: **DRAFT v0.9.14.**
(v0.9.14: the A/B systematic enters as the union of the two confidence sets, per event and
for the catalog, instead of a decision-level threshold that was measured and rejected.)
Previous: **DRAFT v0.9.13.**
(v0.9.13: Gaussian-belt coverage becomes a stopping rule; the A/B systematic moves to the
scan, with a decision-level criterion; numbers of §11 collected.)
Previous: **DRAFT v0.9.12.** (v0.9.12: driver with verified cache and strain
guard (`module1_run.py`); off-source check of the belts specified (set C, paired reference
control, catalog-level binomial acceptance, guard around all listed events).)
Previous: **DRAFT v0.9.11.** (v0.9.11: catalog/Neyman machinery checked on
three events with synthetic noise; reporting rules for concentration (leave-one-out,
N_eff) and for non-contiguous intervals; belt sizes; belts from Gaussian noise with real
off-source segments as a check.)
Previous: **DRAFT v0.9.10.** (v0.9.10: GR-leakage rule and cohort roles fixed
by decision, after the O3 XPHM real-geometry runs; joint-template conditioning recorded.)
Previous: **DRAFT v0.9.9.** (v0.9.9: primary statistic changed to a profile
scan over Λ with the exact dispersive phase; the linear estimator fails at the scale of its
own σ (real-geometry check, §12). Neyman construction required; GR-leakage test added;
z from the posterior median distance; shape identifiability recorded as a limitation.)
Previous: **DRAFT v0.9.8.** (v0.9.8: calibration source fixed — bilby
recalib priors, not the plotting envelope, whose median has the opposite sign in
GWTC-4.0 files; envelope columns verified to be ±1σ; bilby 2.8.2 added to the
environment.)
Previous: **DRAFT v0.9.7.** (v0.9.7: §6 CORRECTED — detector-calibration
directions are no longer projected out freely; they enter as a constrained noise
covariance (generalized least squares). The v0.9.6 estimator could not see any f³
effect; found by a synthetic test of the estimator core before any real residual.)
Previous: **DRAFT v0.9.6.** (v0.9.6: residual chain verified on all
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
**Linear approximation: not valid at the scale that matters (v0.9.9).** The linear
estimator (§6, items 3–4) measures only the part of T₃ that survives the projection
(20–28% for GW230627_015337). Second-order terms of e^{iδΨ} project into the same small
subspace and dominate early: with the exact phase injected (no noise), the linear
estimate returns 0.93–0.97 of Λ at 0.1σ, 0.69–0.82 at 0.3σ, 0.02 to −0.03 at 1σ and
−0.52 at 3σ, at a signal-weighted phase ‖δΨ h‖/‖h‖ of only 0.13–0.17 rad at 1σ
(`module1_event.py`, GW230627_015337, both models). A bound from the linear σ would
exclude values the estimator cannot see. Therefore the **primary statistic is a profile
scan over Λ with the exact phase** (§6, item 7); the linear estimate and T_p are kept as
diagnostics and for the subspace geometry.

**Shape identifiability (limitation, v0.9.9).** After the projection the templates
T₂, T₃, T₄ are correlated at 0.90–0.98 in the C⁻¹ metric (GW230627_015337:
ρ(f²,f³) = 0.975, ρ(f³,f⁴) = 0.970–0.981), and f² or f⁴ injections at 5σ give
z(f³) ≈ +4.8 to +5.0. Within one event the method responds to any smooth phase
deviation growing at high frequency, not specifically to f³. The bound on Λ remains
valid because it assumes the α = 4 model. A significant result would be a candidate
for a dispersion of the class α ≈ 2–4, not specifically for the Λ model; shape
discrimination is attempted only at catalog level (§10) and its power is reported.
Joint-template conditioning (v0.9.10, from the measured correlations): the 3×3
correlation matrix of T₂, T₃, T₄ has condition number 2100–4800 and smallest eigenvalue
6·10⁻⁴–1.4·10⁻³ (GW150914, GW200129_065458, GW230627_015337); a joint fit of f², f³, f⁴
would inflate σ(f³) by 22–33 (by 4.4–5.4 for f², f³ only). Per-event power
identification is therefore not attempted; the catalog-level joint matrix and its
conditioning are computed and reported (§10).

## 3. Cohorts

- **Discovery cohort:** GWTC-1, GWTC-2.1 (O3a) and GWTC-3 (O3b).
- **Confirmation cohort:** GWTC-4.0 (O4a), fixed in advance.
- **Cohort roles (fixed in v0.9.10; limitation recorded).** The cohort is determined by
  the observing run only: O1–O3 = search, O4 = confirmation. No event is moved between
  cohorts according to SNR, waveform model, reference-fit quality or any estimate.
  Recorded limitation: the per-event spread of the scan is about 1·10⁻⁹ m² for the O3
  events checked and about 3.5·10⁻¹¹ m² for GW230627_015337 (O4), so the search cohort is
  expected to be roughly 5–10 times less sensitive than the confirmation cohort. The
  two-stage candidate logic is therefore unbalanced: a candidate first appearing at the
  O4 sensitivity cannot be produced by the search stage. The combined bound uses both
  cohorts; the candidate procedure is kept as registered and its imbalance is reported. Later data (e.g. O4b)
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
- **Primary catalog test: 27 discovery + 20 confirmation events** (v1.0a).
- **Four further discovery events excluded (v1.0a, found during the official build, before
  any estimate): GW170608, GW190707_093326, GW190728_064510, GW190924_021846.** No label
  of these files contains calibration priors (`recalib_<IFO>_*`) for any detector —
  checked over every label present, not only the two of the primary test. This is the rule
  already applied to GW191204_171526 above ("no calibration envelope in the public file,
  for either model"); it was not caught by F5 because F5 tested the plotting table
  `priors/calibration`, whereas v0.9.8 moved the calibration source to the bilby priors
  (§6, item 2). The exclusion is mechanical and was applied before any residual estimate.
  These events cannot be analysed with the frozen calibration model at all, so unlike the
  other excluded events they are not analysed separately either.
- Earlier count before v1.0a: 31 discovery + 20 confirmation.
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
- **Redshift for the template prefactor (v0.9.9):** z is computed from the posterior
  median luminosity distance of the label (Planck15, `module1_cosmo.py`), not from the
  maximum-likelihood sample. Reason: the maximum-likelihood distance differs between the
  two models of the same event by up to 55% (GW200129_065458: z = 0.160 vs 0.103), and
  Λ scales as 1/I₄(z); the median is far more stable (GW230627_015337: 0.066 vs 0.065).
- Calibration model: the bilby recalib priors of the label (Gaussian mu and sigma of
  the amplitude and phase at each spline node; node frequencies from
  recalib_<IFO>_frequency_k); if absent (LALInference SEOBNRv4PHM labels), those of the
  XPHM label of the same event (same detectors and strain); if absent in both, the
  event is excluded. Details and the reason for not using the envelope table: §6, item 2.
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

For each event, detector network and model. Implementation of the linear algebra:
`scanner/gw/module1_estimator_core.py` (no lalsuite dependency; synthetic self-test
in its `__main__`).

1. **Free nuisance subspace (projected out exactly).** Directions that the GR fit can
   absorb without any prior bound:
   - tangent subspace from the posterior (F3): finite-difference derivatives were
     found unstable (near-degenerate parameter directions; numerical noise of the EOB
     and XPHM-SpinTaylor models). Instead, two DISJOINT random sets A and B of n = 800
     posterior samples of the label (fixed seed 20261001; all samples if fewer) give
     waveforms h(θ_i); the differences Δ_i = h(θ_i) − h_ref, whitened with the
     released PSD, are decomposed by SVD, and the principal components explaining 99%
     of the variance (rule frac ≥ 0.99) are kept;
   - the analytic coalescence-time direction −2πif·h_ref;
   - the overall amplitude and phase directions h_ref and i·h_ref.

   Q_A is an orthonormal basis of these directions for set A. **Set A defines the
   primary subspace;** set B is used only for the subspace systematic (§7).
2. **Detector calibration (constrained, NOT projected out).** Calibration errors are
   unknown but bounded. Model and numbers are those of the LVK PE run itself: bilby's
   cubic-spline calibration (bilby 2.8.2, `CubicSpline`, LIGO-T2300140): nodes
   log-spaced in frequency (20–896 Hz in the O3 releases checked, 20–1792 Hz in O4),
   δA(f) and δφ(f) cubic splines in log₁₀ f through the node values (not-a-knot), and
   factor (1 + δA)(2 + iδφ)/(2 − iδφ). Per detector:
   - the reference waveform is multiplied by this factor at the prior means μ of the
     nodes;
   - the spread enters J: amplitude node k gives σ_A,k · B_k(f) · h_ref, phase node k
     gives i · σ_φ,k · B_k(f) · h_ref, with B_k the spline response to a unit value at
     node k (first order of the factor);
   - C = I + J₁J₁ᵀ, J₁ = (I − Q_A Q_Aᵀ) J, applied through the Woodbury identity
     C⁻¹v = v − J₁(I + J₁ᵀJ₁)⁻¹J₁ᵀv.
   Implementation: `module1_calibration.py` (self-test: bilby spline equals a
   not-a-knot spline in log₁₀ f to 3·10⁻¹⁶; first-order directions agree with the exact
   factor to 1.7% for a 1σ draw of 3% calibration, as expected at second order).
   **Calibration protocol (fixed in v0.9.8, before any Λ estimate).** The calibration
   priors describe the detectors and the strain, not the waveform model. For a label
   without recalib priors, "the XPHM label of the same event" means the
   IMRPhenomXPHM PE run of the same event on the same detectors and strain; it is used
   only as the place where these priors are stored. Only the priors (μ, σ, node
   frequencies) are used, never posterior calibration samples of any label, so the
   calibration of the test does not depend on the posterior of the model being tested.
   The rule was chosen from the file format alone (`inspect_calibration.py`,
   `check_calib_sigma.py`), before any residual projection or Λ estimate, and is not
   changed according to the F3/PCA results or any later result.
   **Why not the envelope table priors/calibration/<IFO>** (`check_calib_sigma.py`,
   GW150914, GW200129_065458, GW230627_015337): its columns are f, A_med, φ_med, A_lo,
   φ_lo, A_hi, φ_hi with lo/hi = ±1σ (prior std / half-width = 0.97–1.03; a few O4 H1
   nodes 0.91–1.08), but in the GWTC-4.0 (C00) labels the prior mean equals −(A_med − 1)
   at every node, while in GWTC-2.1/3 (C01) labels it equals +(A_med − 1). Applying the
   envelope median would have doubled the calibration offset in O4 instead of removing
   it. The recalib priors are the input of the PE run and have one convention in all
   catalogs. The maximum-likelihood sample's calibration was not chosen because the
   SEOBNRv4PHM labels have no calibration columns, so the rule would differ between
   models.
3. T_p^⊥ = (I − Q_A Q_Aᵀ) T_p, r^⊥ = (I − Q_A Q_Aᵀ) r.
4. Per-event estimate (generalized least squares):
   Λ̂_i = ⟨T₃^⊥, C⁻¹ r^⊥⟩ / ⟨T₃^⊥, C⁻¹ T₃^⊥⟩,  σ_i = ⟨T₃^⊥, C⁻¹ T₃^⊥⟩^(−1/2);
   the empirical σ_i from §7 replaces σ_i if they differ. Per-detector estimates
   (§10) use the same formula with the detector's own vectors.
5. The same for p = 2 and p = 4 (control amplitudes).
6. Reported per event as diagnostics (no thresholds): the survival
   s = ⟨T₃^⊥, C⁻¹ T₃^⊥⟩^(1/2) / ‖T₃‖ (it sets the estimator noise, σ ∝ 1/s), the cosine
   between T₃^⊥ from sets A and B (it controls the leakage of reference-point
   imperfections), and |s_A − s_B|/s_B. For r = Λ T₃ + noise the estimate is unbiased
   for any free subspace, since (I − Q_A Q_Aᵀ) T₃ = T₃^⊥. Unbiased in expectation does
   not mean usable: when s → 0, σ grows without bound and a single estimate carries no
   information.

7. **Primary statistic (v0.9.9): profile scan over Λ with the exact phase.** For Λ on a
   grid, the residual model is s(Λ) = h_ref (e^{iδΨ(Λ)} − 1) with the group-velocity
   phase of §2, and
   χ²(Λ) = ⟨P(r − s(Λ)), C⁻¹ P(r − s(Λ))⟩,  P = I − Q_A Q_Aᵀ.
   The dispersive phase is exact at every Λ; the GR parameters stay linearised (free
   subspace) and calibration constrained (C), as above. The event estimate is the grid
   minimum. For the catalog the grid is common and in absolute units (m²), and
   χ²_cat(Λ) = Σ_i χ²_i(Λ). **Intervals come from a Neyman construction** with
   exact-phase injections in noise (synthetic first, then real off-source noise, §7);
   the Δχ² ≤ 2.71 (Wilks) interval is not used: its coverage on the real geometry is
   0.74–0.87 instead of 0.90 (§12). Implementation: `ProfileScan` in `module1_event.py`
   (one projection per grid point; a noise trial costs one matrix-vector product) and
   `module1_catalog.py` (the noise enters linearly, so the catalog needs only the sums
   over events of the Gram matrices and of the noise projections; test statistic
   q(Λ) = χ²_cat(Λ) − min χ²_cat, Feldman–Cousins ordering).
   **Belt sizes (fixed in v0.9.11):** 5000 trials per belt; coverage checked with 2000
   independent trials at five true values; grid ±8 × the smallest σ_lin of the catalog,
   step 0.05 × that σ; the fraction of intervals touching the grid edge is reported and
   must be below 1%, otherwise the grid is widened before any real estimate.
   **Belt noise (fixed in v0.9.11):** a belt needs thousands of noise realisations, while
   each event has only 4–8 valid off-source segments. Belts are therefore built from
   Gaussian noise with the released PSDs plus calibration errors drawn from the priors.
   The real off-source segments (§7) are scanned as a check of the belts
   (`module1_offsource.py`, stage `offsource` of `module1_run.py`; v0.9.12):
   - **Injections:** on every segment two injections on the same noise: (C) a posterior
     sample of **set C** — the indices following those consumed by sets A and B in the
     fixed permutation (seed 20261001; 2n attempts plus failures, recorded in the
     cache), hence disjoint from both and deterministic — projected with its own sky
     position, polarization and time; (ref) the reference waveform itself, as a paired
     control. Set B is not used because it defines the second subspace of the A/B
     systematic. The difference between C and ref on the same noise separates the
     signal–reference mismatch from non-Gaussian noise.
   - **Chain:** the residual chain of the real analysis (time shift ±100 ms, complex
     amplitude, against the reference), then the profile scan on the event's own grid;
     q(0) is compared with the event's own Gaussian-belt critical value at Λ = 0.
   - **Acceptance (catalog level, per model, C injections only):** coverage over all
     segments of all event–model pairs (about 400 per model). Below the lower end of
     the 99% binomial interval around 0.90 (≈ 0.861 for 400) the analysis stops before
     any on-source estimate; above the upper end (≈ 0.939) the belts are reported as
     conservative. Per-event coverage (8 segments: 0.90 × 8 = 7.2 ± 0.85) and the ref
     injections are diagnostics only.
   - **Guard:** a segment is refused if it overlaps the window
     [t − duration − 2 s, t + 4 s] of any event of `event_list_v1.csv`, selected or
     excluded, not only the event analysed.
   - **Limitation:** software injections do not pass through the detector calibration;
     this check covers non-Gaussian noise and the signal–reference mismatch, not
     calibration (§9).
   - **Observed in the trial (v0.9.13, GW230627_015337 XPHM-SpinTaylor, 8 segments):**
     q(0) of the C and of the reference injection differ by less than 0.2 on every
     segment, and the two fail on the same segment. The posterior subspace absorbs the
     signal–reference mismatch, so in this case the statistic is set by the noise. The
     refinement behaves as expected: the time shift is one grid step (0.05 ms) for C and
     zero for the reference, |a| 0.89–0.99 for C and 0.96–1.02 for the reference.
   **Non-contiguous intervals (fixed in v0.9.11):** because of the sign ambiguity of the
   nonlinear phase, a confidence set may consist of separate pieces (12% of the
   expected-interval trials in the three-event check). The reported interval is the
   envelope (smallest to largest included Λ), which is conservative; the full set is
   reported alongside.

**Why calibration is not projected out (v0.9.7 correction).** v0.9.6 added the
calibration spline directions to the free subspace. A calibration phase error and
the f³ dispersion phase are both smooth functions of frequency over the band, so ten
unconstrained spline nodes per detector can reproduce almost all of T₃. The synthetic
test of the estimator core showed it: with free calibration the surviving fraction
of f³ was 0.000 and an injected signal was recovered at about 15% in a single
realisation, i.e. the test could not have seen any effect. In the toy self-test of
`module1_estimator_core.py` (flat envelope 5% / 3°, 200 realisations): free
calibration s = 0.003, σ/Λ_inj = 83; constrained calibration s = 0.036, σ/Λ_inj = 7,
with the empirical scatter equal to the analytic σ in both cases; null (Λ = 0,
calibration drawn from the envelope, 300 realisations) mean z = −0.001, std z = 0.987.
These toy numbers only demonstrate the mechanism; the survival with the real,
frequency-dependent envelopes is measured on the real geometry before v1.0 (§12).

## 7. Null calibration (real detector noise)

For every event: off-source segments of real noise from the same detectors and the
same observing period (excluding times near known candidates and data-quality
vetoed times), with the event's own GR reference waveform injected. The full
pipeline is run and the distribution of Λ̂ under Λ = 0 is recorded: bias, scatter,
empirical σ_i, and the per-event and catalog-level false-positive rates. The number
of off-source injections per event is fixed in v1.0 from the compute budget (F6).
**Gaussian-belt coverage (stopping rule, v0.9.13).** Independently of the off-source
check (§6, item 7), the coverage of the Neyman intervals is measured on the same
Gaussian noise model the belts are built from (2000 independent trials at five true
values, as in `module1_catalog.py`). By construction it must be ≈ 0.90: it cannot fail
for a physical reason, only through an error in the code, the grid or the number of
trials. Rule: each of the five coverages must lie inside the 99% binomial interval
around 0.90 for the number of trials used (≈ 0.883–0.917 for 2000); outside it, the
analysis stops and the cause is found before any real estimate. For the union of the A
and B sets (§7) the rule is one-sided: stop only below the lower end; above the upper end
the construction is reported as conservative, which is the expected behaviour of a union. The two checks have
different meanings: the Gaussian one tests the internal consistency of the construction,
the off-source one tests whether the method is valid in real noise.
Off-source segments are taken on both sides of the event, skipping invalid data (gaps,
file edges); at least 4 valid segments are required per detector.

**Subspace systematic (moved to the scan in v0.9.13).** Every null and ladder injection
is analysed twice, with the set-A and the set-B subspace, on the same noise. Since the
scan is the primary statistic, the systematic is measured on it and not on the linear
estimate (which is kept as a diagnostic and for comparison with the pre-v0.9.9 numbers).
**How it enters the result (fixed in v0.9.14): the reported confidence set is the union
of the 90% sets obtained with subspace A and with subspace B**, both per event and at
catalog level, each from its own Neyman belt with the same number of trials. The union
costs a second catalog belt (the belts are not interchangeable: cos(T₃⊥ A, B) is
0.88–0.95 on the real geometry), roughly 6–10 h of the ≈ 45 h budget; the geometry is
not rebuilt, since both subspaces are built anyway. No threshold is needed for the
disagreement: where A and B differ, the interval simply widens.
Consequences, measured on the real geometry (GW150914 XPHM with GW230627_015337
XPHM-SpinTaylor, 2000 trials, ±0.007): (i) the coverage of the union is 0.913–0.949 at the
five true values, against 0.895–0.909 with subspace A alone — the largest value, 0.949 at
Λ = +5.9·10⁻¹¹ m², is about 7 standard errors above 0.90 and is the measure of how
conservative the union is, not an anomaly; the expected interval is essentially unchanged
([−4.99·10⁻¹¹, +4.80·10⁻¹¹] m² against ±4.8·10⁻¹¹), so the widening affects only the few
trials in which A and B differ; (ii) the fraction of contiguous sets falls from 0.85 to
0.86 at catalog level and from 0.79 to 0.83 in the leave-one-out, i.e. the effect is small
here — the construction rule of §7 is therefore one-sided for the union
(stop only below the lower end of the interval; above the upper end is reported as
conservative). The envelope rule of §6 item 7 applies to the non-contiguous sets.
**Stopping rule (quantitative, kept):** the median |Λ̂_A − Λ̂_B| of the scan estimates
over the null suite must be ≤ 0.5 σ_lin; an event–model pair exceeding it moves to the
separately reported set.
**Rejected threshold, recorded:** a decision-level criterion (fraction of null trials in
which A and B disagree about whether Λ = 0 is inside the 90% interval ≤ 5%) was written
in v0.9.13 and then measured on the real geometry: 0.09 (GW150914 XPHM) and 0.07
(GW230627_015337 XPHM-SpinTaylor), 100 trials each (±0.03). Both would have failed from
the geometry alone, while the quantitative criterion passed with a large margin (median
difference 0.05 and 0.125 σ_lin). The threshold measures the sensitivity of a binary
decision near its critical value, not the instability of the subspace: at Λ = 0 the belt
rejects in 10% of trials by definition, so a marginal difference flips the decision, and
non-contiguous sets add jumps between minima; independent subspaces would disagree in
≈ 18%, identical ones in 0%. The threshold was therefore dropped rather than raised
until it passed, and the disagreement is reported as a diagnostic.

**GR-leakage test (v0.9.9).** Real GR waveforms h(θ_B) of set B (not used for the subspace
of set A) are used as data, h(θ_B) − h_ref plus noise, and scanned. The quantity compared
(`leak_std` in the output) is the standard deviation of the **scan estimate**, in units of
σ_lin, against the same quantity for pure noise (`scan_+0sig_std`) — not a z-score and not
a comparison with 1. On the real geometry it is 0.95 vs 1.13 (GW150914 XPHM) and 1.48 vs
1.66 (GW230627_015337 XPHM-SpinTaylor): the leakage distribution is 11–16% *narrower*
than pure noise, consistent with the per-segment check of §6 item 7 (q(0) of the set-C and
of the reference injection differ by less than 0.2). Building the belts from set-C
injections instead of Gaussian noise was considered and rejected on these numbers: the
mismatch does not widen the distribution, and it would make the belt depend on the
particular samples and mix the check into the construction. A GR-only signal must
give the same distribution of the scan estimate as pure noise; a shift or a wider spread
means that GR nonlinearity outside the linearised subspace leaks into Λ. Rule (fixed
in v0.9.10): the **primary test** is a two-sided two-sample Kolmogorov–Smirnov test
between the leakage and pure-noise scan estimates, at least 60 waveforms each, at
p < 0.01. **Action follows the primary test only:** a significant KS result in either
direction moves the event–model pair from the primary test to the separately reported
set. Reported with it, as interpretation (no effect on the action): KS not significant →
no detected difference; significant with a wider leakage distribution or a larger tail
fraction P(|est| ≥ 1σ_lin) → potential false-positive inflation; significant with a
narrower leakage distribution → a statistical difference without evidence of
false-positive inflation. The tail fraction is a secondary diagnostic, not a
replacement of the KS test. Note: this rule was fixed after the O3 XPHM runs, in which
GW150914 XPHM showed a narrower leakage distribution (std 0.83 vs 1.13, 60 trials); the
two-sided form was kept so that the rule is not adjusted to that result, accepting that
such a pair may leave the primary test. This tests posterior-scale GR variations only; waveform-model
systematics are tested by the two-model comparison (§10).

## 8. Injection ladder

| Injection | Checks |
|---|---|
| GR, Λ = 0 | null calibration (§7) |
| f³ at pre-fixed amplitudes of both signs, absolute Λ on the common grid | recovery of sign and amplitude; Neyman belts |
| f² at matched amplitude | how a wrong shape maps into Λ (reported; see §2, shape identifiability) |
| f⁴ at matched amplitude | how a wrong shape maps into Λ (reported) |

All injections use the exact phase factor e^{iδΨ} (for f² and f⁴ the same form with f^p),
never the linear template. Amplitudes include the LVK GWTC-4.0 bound values
(−2.4·10⁻¹¹ and +7.4·10⁻¹² m²). They are fixed in v1.0 and not changed after the first
results.

## 9. Detector-calibration systematics

GR injections multiplied by calibration perturbations drawn from the published
LVK calibration-uncertainty models (O1–O3; O4), without dispersion. Criterion:
calibration-only injections must not produce a catalog-level f³ signal with a
consistent sign and the K(z) dependence of the Λ hypothesis.
Since v0.9.7 the estimator models the calibration spread through C (§6); these
injections test that model with real envelopes and real noise, including the
part of the calibration error that the spline-node model does not describe.

## 10. Catalog-level tests

- Combined result (v0.9.9): χ²_cat(Λ) = Σ_i χ²_i(Λ) on the common absolute grid (K(z) is
  inside the model, so Λ is common to all events); estimate at its minimum, interval
  from the Neyman construction (§6, item 7). The weighted mean of linear estimates is
  reported only as a diagnostic.
- **Sign consistency:** a real effect gives Λ̂_i of one sign across events.
- **Heterogeneity:** Q = Σ (Λ̂_i − Λ̂)²/σ_i², calibrated by the null suite; a real
  common Λ gives Q consistent with the null, systematics tend to give excess Q.
- **Permutation control:** K(z) is shuffled among events while the raw per-event
  f³ amplitudes are kept; this destroys the propagation relation but keeps waveform
  systematics, noise, SNR distribution and uncertainties. The observed catalog
  statistic is compared with the permutation distribution.
- **Concentration (fixed in v0.9.11):** the effective number of events
  N_eff = (Σ w_i)² / Σ w_i², w_i = 1/σ_i² (linear σ), is reported, and the result is
  recomputed leaving out the event with the largest weight (leave-one-out) and reported
  next to the main result. In the three-event check one event carried more than 99% of
  the weight (N_eff ≈ 1.01).
- **Waveform robustness:** the result must be consistent between the two models.
- **Per-detector consistency:** Λ̂ is also estimated separately for each detector
  (same reference waveform, the detector's own projection and PSD). A propagation
  effect is the same in all detectors (same Λ, same source distance); detector
  calibration errors, glitches and noise are not. Statistic per event: χ² of the
  per-detector estimates about their weighted mean; its null distribution comes from
  the null suite.

## 11. Decision rules (fixed in v1.0)

**Order of the stages (fixed in v0.9.15).** build → checks → catalog (rules 0, 1, 1a, 2a)
→ offsource (rule 0a) → **only then** the real estimates, and only for the search cohort;
the result of the search cohort is committed before the confirmation cohort is run. No
on-source estimate is made before every stopping rule below has been evaluated and
recorded.

0. **Construction failure (v0.9.13; one-sided for the union, v0.9.15):** any of the five
   Gaussian-belt coverages **below** the lower end of the 99% binomial interval around
   0.90 (≈ 0.883 for 2000 trials) → stop. Values above the upper end are expected for the
   union of the A and B sets (§7) and are reported as conservative. More than 1% of the
   expected intervals touching the grid edge → stop and widen the grid (§6, item 7). In
   both cases the cause is found before any real estimate.
0a. **Off-source failure (v0.9.12):** catalog-level coverage of the set-C injections
   below the lower end of the 99% binomial interval around 0.90 (≈ 0.861 for 400
   segments) → stop before any on-source estimate (§6, item 7).
1. **Null failure (made numerical in v0.9.15):** over the null suite of an event–model
   pair, |mean Λ̂| > 0.3 σ_lin (GR injections, scan estimate) → the pair leaves the
   primary test; if this happens for more than one third of the pairs, the catalog result
   is not interpreted. The threshold is set from the trial geometry, where the null mean
   of the linear estimator is within ±0.16 σ (§12), i.e. it flags a bias about twice the
   largest observed, not a typical fluctuation.
1a. **Subspace systematic (v0.9.14):** the reported confidence set is the union of the A
   and B sets (§7); an event–model pair with median |Λ̂_A − Λ̂_B| > 0.5 σ_lin leaves the
   primary test. The decision-level disagreement is a diagnostic (the 5% threshold of
   v0.9.13 was measured and rejected; see §7).
2. **Shape (changed in v0.9.9):** per-event shape discrimination is not possible (§2);
   the former rule "f² or f⁴ recovered as f³ → failure" would always trigger and is
   replaced by reporting how f² and f⁴ injections map into Λ. A significant result is
   interpreted only as a candidate for the class α ≈ 2–4.
2a. **GR-leakage failure (v0.9.9; two-sided rule fixed in v0.9.10):** an event–model pair
   with a significant two-sided KS test (§7) leaves the primary test, in either
   direction; if more than one third of the pairs fail, the method is
   not interpreted.
3. **Calibration failure (made numerical in v0.9.15):** calibration-only injections
   (§9) give a catalog-level |Λ̂| > 1.0 σ_cat, or a catalog-level set that excludes
   Λ = 0 → the calibration confound is unresolved and no bound is reported.
4. **All controls pass:**
   - the search-cohort set contains Λ = 0 → **bound** on Λ (the union set of §7, reported
     as the envelope with the full set alongside, with N_eff, the leave-one-out result and
     the conversion to A₄ of §2), after which the confirmation cohort is run and reported
     in the same way;
   - the search-cohort set excludes Λ = 0 → **candidate only**. It must then pass, in
     this order: the confirmation cohort (its set must also exclude Λ = 0, with a
     consistent sign and a value within the search-cohort set), the two waveform models,
     the calibration controls (rule 3), the K(z) permutation test and the per-detector
     consistency test (catalog-level χ² not in the upper 1% tail of its null
     distribution). A p-value alone is not called a discovery, and a significant result
     remains a candidate for the class α ≈ 2–4 (rule 2), not specifically for the Λ
     model.

Every result, including failures, is reported.

## 12. Feasibility checks before v1.0

- **F1** DONE: lalsuite (pip) with lal 7.7.1 / lalsimulation 6.2.1, pyseobnr 0.3.7,
  gwpy 4.0.2, bilby 2.8.2 (added in v0.9.8), Python 3.11 in WSL (`requirements_gw2.txt`, `environment_gw2.yml`).
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
- **Estimator core** DONE (v0.9.7, `module1_estimator_core.py`): synthetic test
  found that freely projected calibration directions absorb the f³ template (see §6);
  corrected to generalized least squares with the envelope as covariance. TO DO before
  v1.0: the same test on the real geometry (reference waveforms, PSDs, PE samples and
  calibration envelopes of 2–3 events, Gaussian noise from the released PSD), which
  measures the real survival s and therefore the reachable sensitivity.
- **Real-geometry check** DONE (v0.9.9, `module1_event.py` v3, synthetic noise only, no
  on-source strain read; `module1_geometry_check.csv`, `module1_geometry_check_v3.csv`).
  Linear estimator, n = 800, frac ≥ 0.99: σ_lin = 8.5·10⁻¹⁰ / 8.3·10⁻¹⁰ m² (GW150914
  XPHM / SEOBNRv4PHM), 6.5·10⁻¹⁰ / 3.5·10⁻¹⁰ (GW200129_065458), 3.6·10⁻¹¹ / 2.1·10⁻¹¹
  (GW230627_015337, XPHM-SpinTaylor / SEOBNRv5PHM); nulls mean z within ±0.16, std
  0.96–1.08; A/B median 0.21–0.38σ. Exact-phase response: see §2 (linear estimator fails
  at ~1σ). Profile scan (GW230627_015337): noiseless Δχ² from Λ = 0 at ±1σ_lin 1.50 /
  0.55 and ≥ 10.8 / 4.1 for |Λ| ≥ 2σ_lin (no blind region within ±8σ_lin); spread of
  the estimate at Λ = 0: 0.97 / 1.66 σ_lin, i.e. ≈ 3.5·10⁻¹¹ m² for both models (the
  model dependence of σ_lin disappears); recovery at +3σ_lin 2.79 / 2.64; Wilks 90%
  coverage 0.78–0.91 / 0.74–0.84; GR leakage (60 set-B waveforms) mean +0.19 / −0.11,
  std 0.95 / 1.45 vs pure noise 0.97 / 1.66. Open: the XPHM null mean +0.25 at about
  2.6 standard errors (100 trials) — to be re-checked with more trials; the
  SEOBNRv4PHM → XPHM calibration branch is exercised in the GW150914 and GW200129 runs.
  O3 with the scan (v0.9.10 record, IMRPhenomXPHM; `module1_geometry_check_v3.csv`):
  GW150914 / GW200129_065458 — exact-phase linear response at 10⁻¹¹ m² 0.997 / 0.986
  (linear regime at the LVK scale for O3); scan Δχ² at ±1σ_lin 0.67 / 1.53, smallest at
  |Λ| ≥ 2σ_lin 6.09 / 5.30 (no blind region); scan spread at Λ = 0 1.13 / 1.21 σ_lin;
  Wilks coverage 0.75–0.83 / 0.69–0.83; GR leakage std 0.83 / 1.08 vs pure noise
  1.13 / 1.21, mean +0.01 / −0.12. z median vs maximum likelihood: 0.096 vs 0.107 and
  0.189 vs 0.160. The current code still uses the maximum-likelihood z in the scan; it
  is switched to the median (§5) before v1.0.
- **A/B union** DONE (v0.9.14, `trial_catalog_union.json`, code 7969f087, status trial):
  numbers above in §7. Wilks under the union is still below nominal (down to 0.835 in the
  leave-one-out), confirming that the Neyman construction is required. Joint f², f³, f⁴
  conditioning with two events: 1506, σ(f³) inflation ×18.5.
- **Driver and cache** DONE (v0.9.12, `module1_run.py`): stages build / checks / catalog /
  offsource / hash; cache entries carry the SHA-256 of the code (all Module 1 modules,
  line endings normalised), seeds, n, PCA rule, MD5 and SHA-256 of the PE file and
  SHA-256 of the input tables, and are used only if all match; strain readers are
  disabled in every stage except offsource; outputs are "official" only if the code hash
  equals the freeze manifest `FREEZE_v1.0.json`, otherwise "trial". Leave-one-out uses its
  own grid when the left-out event set the common grid. Trial (code 68b00bc6, status
  trial): GW150914 XPHM with GW230627_015337 XPHM-SpinTaylor — all numbers reproduced from
  the cache; catalog Neyman coverage 0.895–0.909 (2000 trials), Wilks 0.83–0.90;
  expected interval ±4.8·10⁻¹¹ m²; leave-one-out (GW150914 alone, own grid)
  [−1.2·10⁻⁹, +1.7·10⁻⁹] m², coverage 0.898–0.918, Wilks down to 0.756.
- **Catalog / Neyman check** DONE (v0.9.11, `module1_catalog.py`, synthetic noise only;
  GW150914 and GW200129_065458 (C01 XPHM) with GW230627_015337 (C00 XPHM-SpinTaylor), mixed
  cohorts for the machinery check only; 1000 belt and 300 independent trials;
  `module1_catalog_check.json`). σ_lin with z median: 9.5·10⁻¹⁰, 5.4·10⁻¹⁰,
  3.7·10⁻¹¹ m² (N_eff ≈ 1.01). Critical values 1.49–4.25 (Wilks 2.71). Expected 90%
  interval under Λ = 0 (median ends) [−5.2·10⁻¹¹, +4.6·10⁻¹¹] m², about 3 times wider
  than the LVK GWTC-4.0 interval [−2.4·10⁻¹¹, +7.4·10⁻¹²] (83 events); 88% of the
  intervals contiguous; none at the grid edge. Coverage at five true values
  (−1.48·10⁻¹⁰ … +1.48·10⁻¹⁰ m²): Neyman 0.873–0.933 (target 0.90; ±0.017 standard
  error, the five values share the same 300 realisations), Wilks 0.840–0.893. Catalog
  joint f², f³, f⁴: condition number 949, σ(f³) inflation ×14.7 (per event 22–33).
  Formula check on a toy problem: sum-over-events χ² equals the direct computation to
  6·10⁻⁶ of a 1.3·10³ range.
- **Calibration source** DONE (v0.9.8): see §6, item 2. bilby 2.8.2 added to `gw2`
  (pip dry run: no change to numpy, scipy, lalsuite, h5py, gwpy; `requirements_gw2.txt`
  updated).
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
