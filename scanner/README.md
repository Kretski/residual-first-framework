# Module 0 — validation of the residual-first scanner

Module 0 tests the scanner on systems where a dispersion-type deviation of known
form exists, before it is applied to gravitational-wave data (Module 1). The
scanner looks for reproducible structure in residuals after subtracting a baseline
model; physical interpretation is a separate step.

Every analysis followed the same protocol: fixed configuration with a hash,
registration by git commit/tag **before** the real data, synthetic null barrier,
injections (amplitude scale, wrong shapes, systematics), split into search and
confirmation, and reporting of every result, including failures.

## Summary

| Environment | Known effect | Result |
|---|---|---|
| Piano (University of Iowa MIS) | stiff-string inharmonicity, positive sign, Δf ∝ n³ | **Detected and confirmed**, shape identified, amplitude consistent with published B values |
| Black Sea stereo field (IFREMER) | capillary term, positive sign | **Not usable quantitatively**: two confounds identified (see below) |
| GLOBEX flume, case A3 (Zenodo) | finite-depth dispersion, **negative** sign | **Negative residual reproducible** in both halves; shape undetermined; amplitude ~1.47× linear theory, not interpretable quantitatively |

**Module 0 does not establish a new physical dispersion term. It establishes that
the scanner responds reproducibly to real dispersion deviations, while real-world
systematics and omitted known physics can dominate the inferred amplitude.**

## 1. Piano — `scanner/piano`, config `5bdf34b190a4b9ef` (v0.4.1)

- Final null barrier (cello, untouched before the final run): 2/89, binomial p = 0.224.
- Search group: structure in 21/21 notes, shape n³; confirmation: 18/19 (95 %).
- Stiff-string coefficient B ≈ 1.1·10⁻⁴ … 4.9·10⁻⁴ (secondary interpretation).
- Limitation: 40 of 86 notes (bass and middle register) passed the exclusion rule.
- v0.3 failed its violin null (FPR 0.368, vibrato and onset detection); reported.

Details: `scanner/piano/README.md`.

## 2. Black Sea — `scanner/waves`, config `069ba514c71704e5` (v0.5.0) / `c66f827f9b3b8c56` (v0.5.1, English translation, identical output)

- Test A (free capillary coefficient) **retired** before the new barrier series:
  ill-conditioned (C, g_eff and k¹/k³ nearly collinear); failures 2/40 and 8/200.
- Test B synthetic null: 0/200. Injections: capillary term ×1 detected 100 %,
  shape distinguishable 60 %; wrong shape k¹ never declared as distinguishable k².
- **Real record** (BS_2011, 2011-10-04 11:38): 'not found', with a non-physical
  baseline fit (g_eff < 0). Diagnostic (`diag_real.py`, development record):
  - |k| ≲ 17 rad/m: ridge visible but ω/√(gk) ≈ 1.06, i.e. a ~6 % shift from
    currents/drift, 40–80× the capillary term and not separable from it without
    an independent current measurement — **the main negative result**;
  - |k| ≳ 20 rad/m: accepted cells at 0.3–0.4 Hz, **probably** leakage of the
    dominant waves amplified by the Laplacian prewhitening (diagnosed, not confirmed).

## 3. GLOBEX A3 — `scanner/globex`, config `4b2d81580df7981b` (v0.1.0)

Data: GLOBEX data base (Michallet et al., doi:10.5281/zenodo.4009405, CC BY 4.0).
Only simultaneously recorded gauges (inshore-trolley blocks, boundaries confirmed
by 50 Hz coherence) are used.

- Synthetic null: 0/200. Recovery ×0.25–×1: detection 100 %, (kh)² identified
  93–100 %, amplitude 0.93–0.98 of expected. ×2: shape read as k³ (higher-order terms).
- Systematics: bound harmonics ×1, ×3 without effect; reflection R = 0.03 keeps the
  amplitude (0.99) but shape distinguishable only 17 %; R = 0.10 null gives 60 %
  false alarms, **all with positive sign**.
- A1 not used for the shape: at kh ≈ 0.5–1 linear theory is locally almost linear in kh.
- **Real A3**: reflection ≈ 1 % (no flag). Negative deviation significant in both
  halves (p ≈ 10⁻⁴⁸ … 10⁻⁵²).
  - Code verdict: 'found but not confirmed' (search winner k²−, confirmation
    winner k³−, both undetermined).
  - Rule in the calibration commit: negative sign detected and confirmed, shape undetermined.
  - The two formulations differ in labelling, not in substance; both are reported.
  - a(q=2) = −0.202 / −0.205 vs linear theory −0.138 / −0.139: ratio **1.47**,
    stable between halves, outside all calibrated systematics (0.90–0.99).
    **Not interpreted as a measured coefficient.**
- Exploratory depth inversion (`diag_globex.py`, not registered): h_eff/h rises
  shoreward (1.02–1.04 → 1.09–1.10); a negative frequency slope of h_eff is
  reproducible in the shallowest block (p = 1.1·10⁻³ and 3.9·10⁻⁴), other blocks
  mixed. **Consistent with nonlinear shoaling; a bottom/depth mismatch alone is
  insufficient to explain the frequency dependence.**

## Lessons carried to Module 1

1. Known physics missing from the baseline (currents, reflection, nonlinearity,
   bed geometry) appears in the residuals with a size comparable to or larger than
   the target effect, and can share its shape.
2. A synthetic null barrier is necessary but not sufficient: every real data set
   revealed a mechanism the synthetic model did not contain.
3. Sign and reproducibility are more robust than amplitude and polynomial order.
4. For gravitational waves the corresponding requirement is a complete baseline of
   known effects: waveform-model systematics, detector calibration errors, and the
   partial absorption of a phase deviation by the parameter fit.

## Open

- Linear-regime (low Ursell number) flume data requested from LEGI Grenoble and
  BSHC Varna would allow a quantitative negative-sign test.
