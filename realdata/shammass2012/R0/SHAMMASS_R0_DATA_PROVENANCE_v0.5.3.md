# Shammass et al. (2012) — R0 Data Provenance and Pre-analysis Protocol

**Status:** v0.5.3 — provenance record.  
**Analysis status:** closed. Synthetic gate verdict `NOT_TESTABLE` for both series (§6.1); no experiment-minus-GPE residual has been calculated, and none will be under this R0.  
**Protocol approval status:** NOT APPROVED. Spec v1.0.3 and this document carry no formal approval record; this is a separate open item (§12) and does not change the analysis status.  
**Repository location:** `realdata/shammass2012/`  
**Framework:** `Kretski/residual-first-framework`  
**Primary source:** I. Shammass, S. Rinott, A. Berkovitz, R. Schley, and J. Steinhauer, "The phonon dispersion relation of a Bose-Einstein condensate," arXiv:1207.3440v2 (2012).  
**Figure used:** Figure 2, PDF page 7.

## 0. Changes from v0.5.2 (documentation-only revision)

- Status line and new §6.1: the synthetic gate of `SHAMMASS_R0_PREANALYSIS_SPEC_v1.0.3` was run; verdict `NOT_TESTABLE` for the large and the small series. Run history, deviations and code audit are referenced, not repeated.
- §4.1: WebPlotDigitizer provenance recorded as `RECORDED (version as reported by the application; exact source not pinned)`; the pixel-to-value map is validated empirically (reproduction to ≤ 1.8e-15 on all 91 points).
- §5: blindness record updated: the gate read experimental `k` only; no residual was opened.
- §10: journal citation recorded; Figure 2 comparison still open.
- §11.2 (new): files excluded from the public repository for copyright reasons (source PDF and every file that embeds the figure image), with their hashes kept here.
- §12: LDA line closed with its qualifier; gate and interpolation/error-model lines closed by reference to spec v1.0.3; approval line updated.
- Status line, §11.0 and §12 distinguish the closed analysis from the open protocol approval, and the local archive from the public repository (after review of the first v0.5.3 draft).
- Section numbering in §0: the two earlier change lists that were both numbered 0.1 are renumbered 0.1 and 0.2. Their text is unchanged.
- No data file, coordinate, calibration, protocol rule or result was changed.

## 0.0 Changes from v0.5.1 (in v0.5.2)

- Corrected §3.5 to record the completed visual count and point-by-point check: 9 hollow circles counted in the original PDF main panel; points 1–9 confirmed; 0 excluded. The legacy count of 7 remains a post-check comparison only.
- Corrected §11.0 to state that the source PDF and four legacy CSVs are present under `source/` and `legacy/`; their recorded byte counts and SHA-256 values were checked against the files in this package.
- Corrected §12 to mark the visual check, source/legacy placement, and manifest verification as complete. Scientific data, WPD project, calibration, and exported coordinate/value fields were not changed.
- Corrected the S35 README to describe the completed checklist and the actual package contents; no generator script is claimed to be included.
- Regenerated `MANIFEST.sha256` after these documentation changes. No GPE residuals, interpolation, fit, or quantitative estimate was performed. R0 remains DRAFT and is not approved.

## 0.1 Changes from v0.4 (in v0.5)

- §3.5: completed the independent visual count and point-by-point marker identity check for the nine `exp_small` points. The original PDF main panel was counted first: 9 hollow circles. Dimitar Kretski then checked points 1–9 against Figure 2; all nine were confirmed, none excluded.
- Added `S35_VISUAL_CHECK/` with the numbered overview, zoom tiles, completed checklist, and procedure/record. The digitized coordinate/value fields were preserved; only `decision` and `reason` were populated.
- Updated the top-level, source, and legacy README files and regenerated `MANIFEST.sha256`.
- No residuals, interpolation, fit, or quantitative estimate was performed. R0 remains DRAFT and is not approved.

## 0.2 Changes from v0.3 (in v0.4)

- §4: operator roles corrected. Digitization, axis calibration and the calibration correction were performed by Dimitar Kretski in WebPlotDigitizer 5.2. The click-order dataset reassignment and its programmatic checks were performed by the assistant (Claude). v0.3 had these roles reversed.
- §4: calibration points are now recorded as read from the archived `wpd.json`, not as user-reported values.
- §3, §6, §8: the unified single-calibration digitization (`unified_v1`) is now the primary dataset; it is complete. Point counts and ranges updated to 16 / 9 / 38 / 28.
- §3.5: the small-condensate series has 9 points in the unified set and 7 in the legacy CSV. This is resolved by the visual check recorded in `S35_VISUAL_CHECK/README_S35_VISUAL_CHECK.md`: 9 hollow circles counted in the original PDF main panel; points 1–9 confirmed; 0 excluded. The legacy count remains 7 and is recorded only as a post-check comparison. No digitized coordinates or values were changed.
- §4: the "if re-digitization is performed" conditionals are removed; re-digitization has been done.
- §7: pixel-resolution lower bound added.
- §11.1: full SHA-256 values and byte counts for the unified archive, its members, the source project and the four exports. All four CSVs were exported by the assistant from `wpd.json`; v0.3's statement that `gpe_large.csv` was "user supplied" is corrected.
- §12: checklist updated.

No other control is relaxed. The power gate, error model, interpolation rules, blindness record, quantitative restriction and journal-version check are carried over unchanged in substance.

## 1. Purpose and scope

This document freezes data identity, provenance limitations, comparison hierarchy, and pre-analysis controls for applying the existing residual-first framework to the Shammass et al. (2012) Figure 2 data.

It is not a results report. No experiment-minus-GPE residual, residual-shape search, coefficient recovery, or \(k^4\) fit is performed or authorized by this draft.

The data directory is `realdata/shammass2012/`, not the repository root. Its layout is given in `realdata/shammass2012/README.md`. Do not modify the frozen B0/B1 or published B2 records. Nothing is to be committed until review and explicit approval.

## 2. Source-paper facts and physical interpretation

Figure 2 contains:
- measured dispersion points for the large condensate, \(\mu/h=2340\) Hz (filled circles);
- measured dispersion points for the small condensate, \(\mu/h=1200\) Hz (hollow circles);
- black curves from a cylindrically symmetric 2D Gross–Pitaevskii equation (GPE) simulation (solid: large; dashed: small);
- green curves for the local-density approximation (LDA);
- a dotted free-particle spectrum.

The paper states that the 2D GPE simulation agrees excellently with the measured values. It also states that the measured and simulated dispersions are depressed relative to LDA at low \(k\), interpreting this through the 3D-to-1D transition and transverse/radial mode discretization. Consequently, experiment-minus-LDA is a contextual comparison to a known approximation, not the primary unexplained residual.

### 2.1 Meaning of a \(k^4\) contribution

A \(k^4\) contribution is not, by itself, evidence of new physics. In the standard Bogoliubov expansion, the finite-healing-length correction appears in the expansion of \(\omega^2\), schematically:
\[
\omega^2 = c^2 k^2\left(1+\frac{\xi^2 k^2}{4}+\cdots\right).
\]
The standard GPE description already contains this physics. Therefore, a residual relative to the 2D GPE could at most test for an additional contribution beyond the standard model represented by that GPE (for example, an excess parameter such as \(\Delta\Lambda\), if its mapping and normalization are separately defined). It cannot be interpreted as a measurement of the full standard \(k^4\) contribution or of \(\Lambda\) without subtracting the standard contribution.

The LDA comparison is not suitable for claiming a new \(k^4\) effect because the paper documents a physical low-\(k\) departure from LDA.

## 3. Data inventory and verified metadata

### 3.0 Comparison hierarchy

1. **Primary:** `unified_v1` — four datasets digitized in one WebPlotDigitizer 5.2 session with one corrected axis calibration (§3.1–§3.4, §4).
2. **Comparison only:** legacy CSVs (§3.6). They are preserved unchanged and are not used as the primary analysis input.
3. **Context only:** LDA comparison table (§3.7).

### 3.1 Unified experimental points — large condensate, 16 points

**File:** `unified_v1/exports/exp_large.csv`  
**Columns:** `k_um_inv`, `omega_over_2pi_kHz`, `px`, `py` (decimal point `.`, separator `,`)  
**Range:** \(k=0.2361\)–\(3.8900\ \mu\mathrm{m}^{-1}\); \(f=0.0818\)–\(1.7827\) kHz.  
**Label:** `unified_v1_digitized_experimental_large`.

### 3.2 Unified experimental points — small condensate, 9 points

**File:** `unified_v1/exports/exp_small.csv`  
**Columns:** as §3.1.  
**Range:** \(k=0.3636\)–\(3.8964\ \mu\mathrm{m}^{-1}\); \(f=0.0913\)–\(1.4301\) kHz.  
**Label:** `unified_v1_digitized_experimental_small`. See §3.5 for the point-count finding.

### 3.3 Unified 2D GPE curve — large condensate, 38 points

**File:** `unified_v1/exports/gpe_large.csv`  
**Columns:** as §3.1.  
**Range:** \(k=0.1280\)–\(4.4115\ \mu\mathrm{m}^{-1}\); \(f=0.0578\)–\(1.9753\) kHz.  
**Label:** `unified_v1_digitized_published_2D_GPE_large_curve`. This is a digitized published simulation curve, not raw output from the authors' GPE solver and not experimental data.

### 3.4 Unified 2D GPE curve — small condensate, 28 points

**File:** `unified_v1/exports/gpe_small.csv`  
**Columns:** as §3.1.  
**Range:** \(k=0.1279\)–\(4.4177\ \mu\mathrm{m}^{-1}\); \(f=0.0433\)–\(1.6033\) kHz.  
**Label:** `unified_v1_digitized_published_2D_GPE_small_curve`. Same caveat as §3.3.

All four exports total 91 points. Every exported \(k\) and \(f\) value is computed by the archived calibration from the pixel coordinates stored in `wpd.json`; the `px`, `py` columns are those stored coordinates.

### 3.5 Finding: small-condensate point count (9 vs 7)

The unified session contains 9 hollow-circle (small-condensate) points. The legacy CSV contains 7 small-condensate points. The primary-source visual identity/count check is complete and recorded in `S35_VISUAL_CHECK/README_S35_VISUAL_CHECK.md`: the original PDF main panel was counted before opening the numbered review sheet and contained 9 hollow circles; points 1–9 were then checked sequentially and all 9 were marked `CONFIRMED`, with 0 `EXCLUDED`.

The difference between 9 primary unified points and 7 legacy points is therefore resolved for the identity/count of the primary series. The legacy count remains 7 and is retained only as a post-check comparison; the legacy file was not consulted when making the point decisions and remains unchanged. No digitized coordinates or values were changed. This visual check does not validate scientific interpretation, estimate digitization uncertainty, or authorize GPE residual analysis.

### 3.6 Legacy datasets — comparison only

These files are preserved unchanged under `legacy/`. Their hashes are recorded in §11.

**`shammass2012_fig2_digitized_v2(3).csv`** — 23 × 12; 16 large and 7 small experimental points. Range \(k=0.2385\)–\(3.8739\ \mu\mathrm{m}^{-1}\); \(f=0.0784\)–\(1.7922\) kHz. Recorded extraction method `manual marker-center digitization`, assistant-assisted, from the published plot. Tool, version, calibration coordinates and operator audit trail are not recovered. Do not claim WebPlotDigitizer for this set. Label: `legacy_digitized_experimental_points`. It stores constant digitization estimates \(\sigma_k=0.0158\ \mu\mathrm{m}^{-1}\) and \(\sigma_f=0.0078\) kHz of undocumented derivation; `experimental_sem_khz` is empty.

**`GPE_large_2D_clean(1).csv`** — 28 × 2; \(k=0.1344\)–\(4.2832\ \mu\mathrm{m}^{-1}\). Label: `legacy_digitized_published_2D_GPE_large_curve`.

**`GPE_small_2D_clean(1).csv`** — 29 × 2; \(k=0.1596\)–\(4.3152\ \mu\mathrm{m}^{-1}\). Label: `legacy_digitized_published_2D_GPE_small_curve`.

The legacy GPE curves share one calibration with each other (user statement), but that calibration is not recorded and is not shown to match the legacy experimental points. None of the legacy calibrations is established to match `unified_v1`.

### 3.7 LDA comparison table — context only

**File:** `legacy/shammass2012_lda_comparison(3).csv` — 23 × 17. Built on the legacy 23 experimental rows. It defines the LDA quantity as the 3D Thomas–Fermi LDA first moment \(m_1/m_0\) of the local Bogoliubov dynamic structure factor and contains a derived `residual_exp_minus_lda_khz` column.

**Provenance gap:** the generating script, its hash, physical parameters and numerical settings are not documented. Label: `derived_LDA_comparison_provenance_incomplete`. This table is not an independent experimental data set and is not rebuilt on `unified_v1` under this R0.

## 4. Unified digitization: calibration and dataset assignment

### 4.1 Tool and operators

- **Tool:** WebPlotDigitizer 5.2 (automeris.io). Project file format version 4.2 (`wpd.json`); archive `info.json` version 4.0. Provenance status: `RECORDED (version as reported by the application; exact source not pinned)`. Whether the web or desktop build was used is not recorded. The pixel-to-value map used by the R0 code reproduces all 91 exported values from their `px`, `py` to ≤ 1.8e-15 (audit stage), which validates the conversion empirically for this project.
- **Source image:** `Screenshot_1.png`, 795 × 694 px, a screenshot of Figure 2 (PDF page 7). SHA-256 in §11.1.
- **Digitization, axis calibration and calibration correction:** Dimitar Kretski, 2026-11-22.
- **Dataset reassignment by click order, programmatic identity check, CSV export and archive assembly:** the assistant (Claude), 2026-11-22, by script, from the user's saved project.

### 4.2 Calibration correction

The initial X2 calibration point (\(k=4\)) was placed at the outer edge of the plot frame, compressing \(k\) by about 8–9%. The initial Y1 point was about 4.5 px above the horizontal axis. Dimitar Kretski corrected both with WebPlotDigitizer's Tweak Calibration function. Values were recomputed from the same stored pixel coordinates. The earlier calibration must not be used.

### 4.3 Calibration points (read from archived `wpd.json`)

Axis type: linear XY (`XYAxes`). Values below are rounded to 0.1 px; full precision is in `wpd.json`.

| Point | Pixel (x, y) | Axis value |
|---|---:|---:|
| X1 | (122.1, 661.2) | \(k=0\) |
| X2 | (729.7, 659.7) | \(k=4\ \mu\mathrm{m}^{-1}\) |
| Y1 | (120.6, 661.2) | \(f=0\) |
| Y2 | (117.7, 44.7) | \(f=2\ \mathrm{kHz}\) |

The same calibration applies to all four unified datasets because they are in one project with one axis set.

### 4.4 Dataset reassignment

During digitization the points went into WebPlotDigitizer's `gpe_small` dataset because the intended datasets had not been selected. The saved source project (`source_project/exp_large_3.tar`, §11.1) therefore contains `gpe_large` 38, `gpe_small` 53, `exp_large` 0, `exp_small` 0, with the corrected calibration.

The assistant split the 53 points by click order, without changing any pixel coordinate:
- points 1–28 → `gpe_small` (checked by script: their pixel coordinates are identical, point by point, to the 28-point `gpe_small` dataset in the earlier saved project `source_project/wshammass_fig2_unified.tarpd_project.tar`);
- points 29–44 → `exp_large` (16 points);
- points 45–53 → `exp_small` (9 points).

Dataset colours were set to `gpe_large` (200,0,0), `gpe_small` (0,150,0), `exp_large` (255,165,0), `exp_small` (160,32,240). Only dataset membership, names and colours changed. The calibration block and the image are unchanged. The result is `shammass_fig2_unified_v1.tar`; the four CSVs were exported from its `wpd.json` by the assistant.

The split relies on the click order being: small-GPE curve, then large-condensate markers, then small-condensate markers. That order is as reported by the operator and is consistent with the point positions. The §3.5 visual check also serves as a check of this split.

## 5. Blindness and analysis history

The LDA comparison and experiment-minus-LDA residual column have already been generated and viewed. The overall workflow is therefore not blind with respect to LDA residuals. Any statement or decision informed by them is exploratory and must be disclosed.

The experiment-minus-2D-GPE residual has not been calculated or inspected, for either the legacy or the unified data. It is the only remaining target that can be kept unopened until the analysis specification is frozen.

Both GPE curves and the experimental markers were digitized with the figure in view. This is unavoidable for plot extraction and is recorded as prior exposure to the plotted shapes. The assistant's handling of the unified data was limited to reassignment, identity checks, export, counts, and min/max ranges; no interpolation or residual was computed.

Accurate description: `LDA comparison previously inspected; experiment-minus-GPE residual held unopened pending R0 freeze`. No claim of a fully blind study is permitted.

Update (v0.5.3): the synthetic gate read only the experimental `k` column; experimental frequencies were not parsed. Because both series are `NOT_TESTABLE`, spec v1.0.3 §9.1 does not authorise a detection analysis, and the experiment-minus-GPE residual remains unopened. Post-gate design calculations used only the GPE curves and the experimental `k`.

## 6. Statistical power and identifiability gate

The primary unified set has 16 large-condensate points and 9 small-condensate points; the §3.5 visual identity/count check is complete. The small series is limited; nearby candidate powers (for example \(k^2\), \(k^3\), \(k^4\)) may be poorly identifiable there, depending on the \(k\)-grid, noise and baseline model.

Before opening experiment-minus-GPE residuals, run a synthetic power/identifiability study using:
- exactly the observed `unified_v1` \(k\) coordinates for each series (after the §3.5 decision);
- the predeclared baseline and candidate model family;
- a predeclared noise model and a sensitivity range for its scale;
- the same interpolation and fitting procedure intended for the real data;
- fixed random seeds and a documented number of simulations.

Report selection/confusion rates, coefficient coverage where meaningful, false-positive rates under the null, and power under prespecified effect sizes. The effect-size grid and noise assumptions must be fixed before the real GPE residual is opened.

### 6.1 Gate result (added in v0.5.3)

The gate specified in `SHAMMASS_R0_PREANALYSIS_SPEC_v1.0.3` (SHA-256 `24d487ab6799cd90868cc5fabe378101b0bd6d705f1b60e6de1022be95fd6914`) was run once on 2026-11-25 with code `007d56ad92cca87959cab349231c751c2532f0c921c64f4974a1ae40036c06d9`. Result file SHA-256 `fd94f4b92bd49e446e1b1c8df513c36bfb2eb327de3d0900b0f4e4b8480ff4e3`.

| Series | G1 | G2 | Verdict |
|---|---|---|---|
| large (N = 16) | PASS | FAIL (largest CP lower bound of q = 2 detection 0.122; required 0.80) | `NOT_TESTABLE` |
| small (N = 9) | PASS | FAIL (0.034) | `NOT_TESTABLE` |

Statement: **not testable with available data.** This is the predeclared inconclusive outcome below, not a failure of the method. The main reason is that the target template \(F(k)k^2\) is strongly collinear with the base term \(F\) and with the neighbouring shapes over the observed \(k\) range.

Recorded deviations (gate run before formal approval of spec v1.0.3; run code not committed before the run; three minor implementation deviations) and the code audit are in `audit/R0_GATE_RUN_RECORD.md` and `audit/R0_GATE_CODE_AUDIT_v1.md`. The verdict does not depend on any of them. A separate post-gate design calculation (`design/`) is an accepted artifact; it does not change this verdict.

**Predeclared inconclusive outcome:** if the study shows inadequate power or poor identifiability for a series, report `not testable with available data` for that series. Do not describe this as a failure of the residual-first method, and do not pool the two series post hoc to manufacture power.

## 7. Error model and uncertainty limits

Experimental SEM values are not available numerically. The paper describes plotted error bars as SEM, mostly too small to see.

The `unified_v1` exports carry no per-point uncertainty. The legacy constants \(\sigma_k=0.0158\ \mu\mathrm{m}^{-1}\), \(\sigma_f=0.0078\) kHz belong to the legacy extraction, have undocumented derivation, and do not transfer to `unified_v1`.

**Pixel-resolution lower bound.** From the §4.3 calibration on the 795 × 694 px screenshot, one pixel corresponds to about \(4/607.6 \approx 0.0066\ \mu\mathrm{m}^{-1}\) in \(k\) and \(2/616.6 \approx 0.0032\) kHz in \(f\). This is a floor on digitization resolution, not an uncertainty estimate. Real marker-centre and curve-tracing error is expected to be several pixels, and is additional to calibration-point placement error and to experimental SEM.

Before residual inspection:
1. define a digitization-uncertainty model for `unified_v1` that is no smaller than the pixel floor, justified independently of the target residual (for example, by repeat placement or marker-size arguments);
2. define a prespecified sensitivity range over plausible total frequency-error scales;
3. keep experimental SEM, experimental-marker digitization error, GPE-curve digitization error and calibration error as separate components;
4. do not report calibrated p-values, confidence intervals, or quantitative physical coefficients unless the uncertainty model supports them.

At most, exploratory sign/shape descriptions may be made under explicitly stated assumptions, accounting for digitization and calibration limits.

## 8. GPE interpolation and domain coverage

`unified_v1` ranges (min/max only):
- `exp_large` \(0.2361\)–\(3.8900\ \mu\mathrm{m}^{-1}\) within `gpe_large` \(0.1280\)–\(4.4115\);
- `exp_small` \(0.3636\)–\(3.8964\ \mu\mathrm{m}^{-1}\) within `gpe_small` \(0.1279\)–\(4.4177\).

Each GPE curve brackets its own experimental series by range. This is a range check only; it does not verify monotonicity, local gaps, point ordering, or interpolation error.

Before residual computation, freeze:
- the pairing: `exp_large` with `gpe_large` and `exp_small` with `gpe_small` only;
- interpolation method (linear or a specified shape-preserving alternative);
- no extrapolation outside the digitized curve domain;
- sorting and handling of duplicate or non-monotonic \(k\) entries in the curve exports;
- interpolation-uncertainty estimation, for example by comparing prespecified schemes or resampling digitized curve points;
- propagation of experimental \(k\)-uncertainty into the matched GPE frequency, where applicable.

No interpolation or residual values are calculated in this draft.

## 9. Quantitative-claim restriction following B2

The synthetic B2 recovery result was classified V4: detection, sign and shape were recovered, but quantitative coefficient recovery was not calibrated (C4 failed). R0 must not use the Shammass application to claim a quantitative \(\Delta\Lambda\) estimate or a calibrated confidence interval unless a separately justified and independently validated recovery procedure is developed first.

For this R0, quantitative \(\Delta\Lambda\) estimation is **out of scope**. The protocol may assess whether a residual is distinguishable from zero and whether a prespecified shape is identifiable, and only if the power and uncertainty gates pass. No physical claim follows from a shape match alone.

## 10. Journal-version check

The source is arXiv:1207.3440v2. Before R0 approval, verify whether the journal-published version has the same Figure 2, labels and plotted data. Record the journal citation and compare the figure. Journal citation (recorded in v0.5.3): Phys. Rev. Lett. 109, 195301 (2012), DOI 10.1103/PhysRevLett.109.195301, title "Phonon dispersion relation of an atomic Bose-Einstein condensate". The figure comparison is still open (subscription access). If it differs, designate the chosen source figure and say why. Do not silently substitute one version for another.

## 11. File integrity record

### 11.0 Legacy and source files

In the **local archive** the source PDF is present under `source/` and the four legacy CSVs under `legacy/`; the byte counts and SHA-256 values below were checked against those files. In the **public repository** the four legacy CSVs are published, but the PDF is not: `source/` contains only a README with the citation and this hash (§11.2). These files are preserved unchanged and remain source/comparison material; their presence does not make the legacy datasets primary analysis inputs.

| File | Bytes | SHA-256 |
|---|---:|---|
| `1207.3440v2(1).pdf` | 456741 | `ef5416e243b0863dada51b68483b72edf83e3b3870bf417b9215b236ea61207c` |
| `shammass2012_fig2_digitized_v2(3).csv` | 4301 | `6243ea439753d532e2fcf319b9478bb84eedf37e1bbf73f96b26a43b1c542bb0` |
| `GPE_large_2D_clean(1).csv` | 857 | `7cf4c44fb5d37cf17a7232dfd3f4b8184f5f8e6bfffe9c1c61c94acebb812bec` |
| `GPE_small_2D_clean(1).csv` | 888 | `5fd5a8ad5bf4c10863d9a34ed688a20014dadc6f0506827c6600f4063ab5cff3` |
| `shammass2012_lda_comparison(3).csv` | 7264 | `b15c8a7f8d07c1242d25d4a4c2e4c0f93762adfb26d36210b38d2f94445ad9cc` |

### 11.1 Unified session

| File (under `unified_v1/`) | Bytes | SHA-256 |
|---|---:|---|
| `shammass_fig2_unified_v1.tar` | 102400 | `c9013218dea7eb40b7e14a44c5f0af476245a40712438b488fe8b6a2fa56b192` |
| `project/wpd.json` | 9495 | `291f83d0a7699f8c867eab6ff7448b2e867342a25009e10cae18356decc98beb` |
| `project/Screenshot_1.png` | 88873 | `3c868d835c788db677b8ec8d03eaf5785a3ce6fc624fd49a658c668d4789ed2e` |
| `project/info.json` | 65 | `577cc2929fdc0b691437162e3d6debe1024f9eff74b61fa55658ad656d44c159` |
| `exports/exp_large.csv` | 1243 | `203637ca80adee49d01bc01a2c9456db5fb5bcbf9e8efac768e8e2537adbb39b` |
| `exports/exp_small.csv` | 716 | `a39239136a7011a794a8e72e2f50f0176fef4e04eb5362076d6ed512e190aa02` |
| `exports/gpe_large.csv` | 2875 | `ddb5e53199d8048aba8ccd08c6905b44151c63e573d11ff7514657f703d2e4d9` |
| `exports/gpe_small.csv` | 2136 | `c19213b54d9c5433265140f8cad1517102733bc5031a581b3bc7e7c26b500066` |
| `source_project/exp_large_3.tar` | 102400 | `94a9e2cc253e4fb50fa5cedbd0cf3fa8aa012cfcebef856ed35acc18131125ac` |
| `source_project/wshammass_fig2_unified.tarpd_project.tar` | 102400 | `d2ecfb8ef9734b5a77696f4aeac3e0d1b35e7e8cafffecf5640200a362baa7ad` |

`project/` holds the three members of `shammass_fig2_unified_v1.tar`, extracted unchanged.

`source_project/` holds the user's two saved WebPlotDigitizer projects, unchanged, as audit evidence only:
- `exp_large_3.tar` — the project after calibration correction and before reassignment (§4.4); the direct input to the reassignment. Its `wpd.json` has SHA-256 `3f3db546b33a1217f7a12fe3600e63c367c722eb9f63e5bac7e93ceb50d2208f`.
- `wshammass_fig2_unified.tarpd_project.tar` — an earlier save (`gpe_large` 38, `gpe_small` 28) that still carries the **superseded, uncorrected** calibration (X2 at about x = 789 px, Y1 at about y = 649 px). Its `wpd.json` has SHA-256 `7e1366173496684eb0229139dc8924711158537dc6baaf27683221d54283f318`. It must not be used for values; it is kept because it documents the calibration correction and is the reference for the `gpe_small` identity check.

Both source images are byte-identical to `project/Screenshot_1.png`.

All four CSVs were exported by the assistant from `project/wpd.json`. The `gpe_large.csv` hash matches the value the user reported for the separately supplied copy. A complete manifest is in `MANIFEST.sha256` at the archive root; verify it on the laptop copy before analysis.

### 11.2 Files excluded from the public repository (added in v0.5.3)

The arXiv abstract page shows no open licence for 1207.3440, so the arXiv default non-exclusive licence applies and does not permit redistribution. The following are therefore kept only in the local archive and are not committed to a public repository. Their hashes above (§11.0, §11.1) and in the local `MANIFEST.sha256` remain the reference.

- `source/1207.3440v2(1).pdf` (the paper);
- `unified_v1/project/Screenshot_1.png` (image of Figure 2);
- `unified_v1/shammass_fig2_unified_v1.tar` and both `unified_v1/source_project/*.tar` (each embeds the same image);
- `S35_VISUAL_CHECK/s35_overview_numbered.png` and `S35_VISUAL_CHECK/s35_zoom_tiles.png` (derived from the image).

The digitized numbers, `wpd.json` (calibration and pixel coordinates, no image) and `info.json` are published. The public copy carries its own manifest listing only the published files.

## 12. Approval checklist

Analysis status: closed (`NOT_TESTABLE`, §6.1). Protocol approval status: NOT APPROVED until the last item below is recorded. The items are:

- [x] `realdata/shammass2012/` created in the intended working checkout and `MANIFEST.sha256` verified on Lenovo (23/23 MATCH; 0 errors).
- [x] Source PDF and the four legacy CSVs present unchanged under `source/` and `legacy/` in the local archive; hashes match §11.0. Public repository: legacy CSVs published, PDF excluded (§11.2).
- [x] Unified re-digitization performed; tool recorded as WebPlotDigitizer 5.2.
- [x] Corrected calibration recorded from the archived `wpd.json`.
- [x] Click-order reassignment, operators, and point counts (16 / 9 / 38 / 28) recorded.
- [x] §3.5 visual count and point-identity check completed: 9 hollow circles counted in the original PDF main panel; 9 CONFIRMED, 0 EXCLUDED; decisions recorded in `S35_VISUAL_CHECK/s35_checklist.csv`.
- [x] Legacy experimental extraction provenance retained as incomplete; no WebPlotDigitizer claim for it.
- [x] LDA provenance explicitly remains incomplete: generating script not recovered; table not used as an input (context only).
- [ ] Journal-version Figure 2 checked against arXiv v2.
- [x] Error model, digitization-uncertainty model and sensitivity scales fixed in spec v1.0.3 §5; no digitization sigma treated as SEM.
- [x] Synthetic power/identifiability gate completed before opening GPE residuals: `NOT_TESTABLE` for both series (§6.1).
- [x] Interpolation, pairing, domain, uncertainty propagation and exclusion rules fixed in spec v1.0.3 §1, §3, §4.
- [x] Blindness limitation recorded.
- [x] Quantitative \(\Delta\Lambda\) estimation excluded from this R0 scope.
- [ ] Explicit R0 approval recorded. (Spec v1.0.3 was not formally approved before the gate run; recorded as a deviation.)

**Current decision:** R0 is closed with the gate verdict `NOT_TESTABLE` for both series. No GPE residual calculation, residual-shape search, coefficient recovery, or \(k^4\) fitting has been performed, and none is authorised under this R0. A successor test (R1) would need a new pre-registered protocol and numerical source data. Open items: journal Figure 2 comparison (§10) and formal approval record.
