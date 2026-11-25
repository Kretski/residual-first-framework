# Shammass et al. (2012) — R0 Pre-analysis Specification

**Document:** `SHAMMASS_R0_PREANALYSIS_SPEC_v1.0.3` — **DRAFT / NOT APPROVED**  
**Supersedes:** `SHAMMASS_R0_PREANALYSIS_SPEC_v1.0.2_DRAFT.md` (SHA-256 `755689a674ec3bf7b915868fcdb7929675e1fa6b897d43889393d72be05d3235`), which superseded `SHAMMASS_R0_PREANALYSIS_SPEC_v1.0.1_DRAFT.md` (SHA-256 `c55680209d3de61b68e94f85cc1a3c6f787c6b1f73a040bdad9297582fe09da8`), which superseded `SHAMMASS_R0_PREANALYSIS_SPEC_v1.0_DRAFT.md` (SHA-256 `ac6c18acc6ed2a9cb42d189c7fc43c4d3003006690881e79e47723a99570b154`)  
**Companion to:** `R0/SHAMMASS_R0_DATA_PROVENANCE_DRAFT_v0.5.2.md` (SHA-256 `5658eb3d07b823afa5e477b56e8eb0f4f80f47ffd4364a8302a018421618e5f4`, 25620 bytes, as received)  
**S35 record:** `S35_VISUAL_CHECK/s35_checklist.csv` (SHA-256 `4ad5234ed0554ff21a76cdf49b696007124e0cd9ecd6225817906cf5ec2a414e`, 2367 bytes, as received; must equal the copy listed in `MANIFEST.sha256`)  
**Location when approved:** `realdata/shammass2012/R0/`

## 0. Changes from v1.0.2

Made after review of v1.0.2 (reviewer HOLD, two items). All other text is unchanged from v1.0.2.

- **§7 and §10.** The Holm family is defined per series. The nominal family-wise error rate of 0.01 applies to each series separately, not jointly; no combined detection claim is made from the two series.
- **§8.5.** G3 has a pre-fixed rule for zero detections: if there are none for a non-target shape, G3 is `NOT_ESTIMABLE`, which does not count as a pass. Such a series can be at most `SIGN_ONLY`.
- **§8.2.** Added a recorded fact on export precision relevant to the reproduction check (no tolerance changed).

## 0.00 Changes from v1.0.1 (in v1.0.2)

Made after review of v1.0.1 (reviewer HOLD, three precision items). Corrections A–E of v1.0.1 are unchanged. All other text is unchanged from v1.0.1.

- **§10.** The first allowed statement now names the Holm procedure at *nominal* family-wise \(\alpha=0.01\) and says the empirical null rejection rate is reported separately. It no longer reads as if the simulation had confirmed an empirical rate of at most 1%.
- **§5.2.** Added a note that the weights are a diagonal variance model only. Interpolated values at different \(k_i\) share curve nodes, so their errors can be correlated. The weights do not describe this; the gate, which perturbs the shared nodes, tests its effect empirically through G1–G3.
- **§8.2.** The calibration-reproduction check uses a fixed mixed tolerance \(|x_{\rm rec}-x_{\rm exp}|\le\epsilon_{\rm abs}+\epsilon_{\rm rel}|x_{\rm exp}|\), separately for \(k\) and \(f\), instead of a purely relative one, which is ill-conditioned near zero.

## 0.0 Changes from v1.0 (in v1.0.1)

Made after review of v1.0 (reviewer HOLD, five required corrections) and the decision of Dimitar Kretski on the false-positive criterion. All other text is unchanged from v1.0.

- **§2.4 (review item C).** The claim that calibration error leaves sign and shape unchanged is withdrawn. A shared calibration preserves the residual up to a common frequency scale, but an error in the origin of the \(k\) axis turns the target \(F(k)k^2\) into \(F(k)k^2\) plus a small \(F(k)k\) admixture. The effect is now quantified and simulated (§8.2) instead of assumed away.
- **§4.3 (item B).** The leave-one-out RMS is described as a curve-only diagnostic, not an upper bound. The interpolation position \(t_i\) is defined.
- **§5.2 (items A and E).** The constant 0.75 is removed; the mean of \((1-t)^2+t^2\) over uniform \(t\) is 2/3, and a constant is not needed. Each point uses its own \(t_i\). Curve-point click error is now two-axis, which adds the slope factor to the curve term.
- **§6 (review remark).** Weighted pairwise correlations between design columns are reported as diagnostics.
- **§8.2 (items C and E).** Noise is injected in pixel space on both axes for experimental and curve points, and the four calibration points are perturbed in every replicate.
- **§8.4, §8.5, §10 (item D; decision of Dimitar Kretski).** Null replicates increased from 2000 to 10,000. G1 now requires the exact two-sided 95% Clopper–Pearson interval to have lower bound ≤ 0.010 and upper bound ≤ 0.015. §10 states the nominal and empirical rates and does not claim that G1 proves the family-wise error rate equals 0.01. Total: 176,000 replicates (each fits the base model and the three extended models).
- **§8.6.** Outputs extended with the design correlations and the calibration-jitter record.
- Former §0 heading renumbered §0.1 to make room for this list; its text is unchanged.

No code has been written or run. No synthetic gate has been executed. No experiment-minus-GPE residual has been computed.

## 0.1 What this document is

It fixes, before any experiment-minus-GPE residual is computed, every analytical choice the R0 provenance document left open in §6, §7 and §8: data preparation, interpolation, uncertainty model, residual model family, test statistic, decision rules, the synthetic power/identifiability gate, and the allowed wording of outcomes.

It contains no residual values, no interpolated GPE values at experimental \(k\), and no simulation output. No code has been run for it.

Order of work after approval of this spec:
1. implement the gate and analysis code (separate review, as B1 was for B0);
2. compute the curve-only quantities of §4.3 (they use no experimental frequencies);
3. run the synthetic gate (§8) — one run, after explicit permission;
4. record the per-series gate verdict;
5. R0 freeze (commit + tag) with this spec, the gate output and the code hashes;
6. only then open the real residual (§9), one run, after explicit permission.

No step may be repeated with changed settings after seeing its output without a recorded deviation (§11).

## 1. Inputs and pairing

| Role | File (under `realdata/shammass2012/unified_v1/exports/`) | Points | SHA-256 |
|---|---|---:|---|
| exp, large | `exp_large.csv` | 16 | `203637ca80adee49d01bc01a2c9456db5fb5bcbf9e8efac768e8e2537adbb39b` |
| exp, small | `exp_small.csv` | 9 | `a39239136a7011a794a8e72e2f50f0176fef4e04eb5362076d6ed512e190aa02` |
| GPE, large | `gpe_large.csv` | 38 | `ddb5e53199d8048aba8ccd08c6905b44151c63e573d11ff7514657f703d2e4d9` |
| GPE, small | `gpe_small.csv` | 28 | `c19213b54d9c5433265140f8cad1517102733bc5031a581b3bc7e7c26b500066` |

- Pairing is fixed: `exp_large` ↔ `gpe_large`, `exp_small` ↔ `gpe_small`. No cross-pairing, no pooling of series.
- Variables used: `k_um_inv` (\(k\)) and `omega_over_2pi_kHz` (\(f\)). The `px`, `py` columns are not used except in §5 for the pixel-to-value scale.
- Legacy CSVs and the LDA table are not inputs.
- Code must verify all four SHA-256 values and the point counts before reading values, and abort on any mismatch.

## 2. Physical target and model family

### 2.1 Target form

With the excess dispersion \(\omega^2 = \omega_{\rm GPE}^2\,(1+\Delta\Lambda\,k^2)\), first order in \(\Delta\Lambda\) gives
\[
f_{\rm exp}(k) - F(k) \approx \tfrac12\,\Delta\Lambda\,k^2\,F(k),
\]
where \(F\) is the GPE curve. The target template is therefore \(T_2(k)=F(k)\,k^2\). At low \(k\), where \(F\approx ck/2\pi\), this is \(\propto k^3\), matching the \(n^3\leftrightarrow k^4\) correspondence of B0.

Sign convention: the sign of the fitted \(T_2\) coefficient is the sign of \(\Delta\Lambda\); \(\Delta\Lambda>0\) is the superluminal (upward) excess, as in the project's \(\Lambda\) convention. No magnitude of \(\Delta\Lambda\) is reported (§10).

### 2.2 Nuisance (base) terms

Always present in every fitted model:
- \(B_0(k)=1\): constant offset. Covers a systematic vertical offset between marker centres and the traced curve line.
- \(B_1(k)=F(k)\): proportional mismatch. Covers a small difference in sound speed or \(\mu\) between experiment and simulation, which is standard physics and must not be read as an excess.

### 2.3 Candidate excess shapes

\(T_q(k)=F(k)\,k^q\), \(q\in\{1,2,3\}\). \(q=2\) is the target; \(q=1\) and \(q=3\) are the neighbouring shapes used to measure identifiability. No other shapes are tested.

### 2.4 Calibration error: what the shared calibration does and does not remove

All four datasets share one calibration (R0 §4.3). An error in the calibration points changes the pixel-to-value map by an affine transform that is the same for the experimental points and for the curve.

**What it removes.** Linear interpolation commutes with affine maps. The residual measured at a given pixel position is therefore changed only by a common frequency-scale factor of about \(1\pm0.005\) (a 1–3 px placement error on a 616 px axis), plus a negligible rotation coupling. A frequency offset cancels in the difference. No independent calibration nuisance term between experiment and curve is needed. This holds only for `unified_v1`; it would not hold for the legacy files.

**What it does not remove.** The templates are evaluated at the measured \(k\). If the true wavenumber is \(k=a\,k'+b\), a true residual \(\gamma F k^2\) appears in measured coordinates as
\[
\gamma F\,(a^2k'^2 + 2ab\,k' + b^2).
\]
The \(b^2\) part is absorbed by the base term \(F\). The \(2ab\,k'\) part is an admixture of the competing shape \(q=1\), with relative size about \(2b/k'\). For an origin error of 2 px (\(b\approx0.013\ \mu{\rm m}^{-1}\)) this is about 7% at \(k'=0.36\), 2.6% at \(k'=1\), and 0.7% at \(k'=3.9\). A scale error \(a\ne1\) only rescales \(\gamma\).

This effect is small but not zero, and it acts on exactly the \(q=1\) versus \(q=2\) distinction. It is therefore not assumed away. It is included in every synthetic replicate through calibration-point jitter (§8.2), so its influence on detection, shape selection and sign is measured by G1–G3.

## 3. Data preparation (fixed)

1. GPE curve: sort by \(k\) ascending. If two curve points have \(|\Delta k|<10^{-9}\ \mu{\rm m}^{-1}\), replace them by one point at their mean \(f\); record the count.
2. After sorting, record (do not edit) any curve segment where \(f\) does not increase. Report the count. The curve is used as digitized.
3. Experimental points: kept in file order; no sorting is needed for the fit.
4. **Domain rule.** An experimental point is used only if \(k_{\min}^{\rm curve}\le k\le k_{\max}^{\rm curve}\) and the bracketing curve segment has length \(\le 0.5\ \mu{\rm m}^{-1}\). No extrapolation. By R0 §8 all 25 points are inside the ranges; the segment-length condition is checked by code.
5. **No other exclusion.** No outlier removal, no point dropping after any residual is seen. Influence diagnostics may be reported but never change the decision.

## 4. Interpolation (fixed)

### 4.1 Primary
Piecewise-linear interpolation of the prepared curve, \(F(k)\), evaluated at each experimental \(k\) (`numpy.interp` semantics, inside the domain only).

### 4.2 Sensitivity
PCHIP (`scipy.interpolate.PchipInterpolator`, no extrapolation) on the same prepared curve. Used only in §9.3.

### 4.3 Curve-only quantities
Computed from the curve alone; no experimental \(f\) is read:
- local slope \(s_i\) (kHz per \(\mu{\rm m}^{-1}\)) of the linear segment containing \(k_i\);
- interpolation position \(t_i=(k_i-k_L)/(k_R-k_L)\in[0,1]\), where \(k_L,k_R\) are the curve nodes bracketing \(k_i\);
- leave-one-out interpolation error for each interior curve point \(j\): \(e_j=f_j-\hat F_{-j}(k_j)\), linear interpolation without point \(j\). Report \(\mathrm{RMS}(e)\) per curve, in kHz and in pixels. This is a curve-only diagnostic estimate of interpolation and tracing sensitivity. It is not a guaranteed upper bound: it mixes local curvature, digitization error and the effect of the wider gap left by removing a point, and with uneven spacing that gap is not simply twice the original. It is reported, and used in §8 only as a check that the noise grid covers it (§5.3).

### 4.4 Propagation of \(k\)-error
An error \(\delta k\) in an experimental \(k\) maps to an apparent residual \(-s_i\,\delta k\). It enters through the variance model of §5.2.

## 5. Uncertainty model (fixed)

### 5.1 Scale
Pixel scales from R0 §4.3: \(a_k=0.0065830\ \mu{\rm m}^{-1}\)/px, \(a_f=0.0032438\) kHz/px, ratio \(a_k/a_f=2.0294\).

Absolute experimental SEM is not available and is **not invented**. The overall noise scale \(\sigma\) is treated as unknown and estimated from the fit residual. All tests (§6) are scale-free.

### 5.2 Relative variance (weights)
For point \(i\), assume isotropic click error \(\delta\) px on both axes, independently for the experimental marker and for each curve node. An error displaced along \(k\) by \(\delta_x\) px on a line of pixel slope \(\sigma_i=2.0294\,s_i\) is equivalent to a vertical error \(\sigma_i\delta_x\). Hence:
- experimental marker: vertical variance \(\delta^2(1+\sigma_i^2)\);
- linear interpolant between two independent nodes at position \(t_i\): \(\delta^2(1+\sigma_i^2)\,[(1-t_i)^2+t_i^2]\), using the slope of the same segment for both nodes.

So
\[
\mathrm{Var}(r_i)\ \propto\ v_i = \bigl(1+(2.0294\,s_i)^2\bigr)\,\bigl(1+(1-t_i)^2+t_i^2\bigr),
\]
and the weights are \(w_i=1/v_i\). They depend only on the curve and on the experimental \(k_i\), not on any experimental frequency.

No constant replaces \((1-t_i)^2+t_i^2\). For reference only, its mean over uniform \(t\) is \(2/3\).

**Scope of this model.** It is a diagonal variance model. Neighbouring experimental points that fall in the same or adjacent curve segments share curve nodes, so their interpolation errors are correlated; the weights do not describe these covariances, and the model is not presented as a full covariance matrix. The synthetic gate perturbs the shared curve nodes (§8.2), so G1–G3 test the effect of these correlations empirically, together with everything else the diagonal weights omit.

Experimental SEM is assumed to scale in the same way; it is absorbed into \(\sigma\). This is an assumption and is stated with every result.

### 5.3 Noise grid for the synthetic gate
\(\delta\in\{1,2,4,8\}\) px.
- 1 px is the resolution floor (R0 §7).
- 8 px exceeds the marker radius (about 6 px at screenshot scale). Error bars larger than a marker would be visible in Figure 2, and the paper states most are too small to see.

If the §4.3 leave-one-out RMS of a curve exceeds 8 px, the grid is extended to include the next power of two. This is a rule fixed now, not a choice made later.

## 6. Fit and test statistic (fixed)

For a series with \(N\) points, residual \(r_i=f_i-F(k_i)\), weighted least squares with \(w_i\):
- base model \(M_B\): \(r=\beta_0 B_0+\beta_1 B_1\);
- extended model \(M_q\): \(M_B + \gamma_q T_q\), for \(q=1,2,3\).

Statistic: partial F for adding \(T_q\) to \(M_B\), with \((1,\,N-3)\) degrees of freedom; p-value \(p_q\). Large: \(N-3=13\); small: \(N-3=6\).

Design matrices are column-scaled to unit weighted norm before solving; solution by QR (`numpy.linalg.lstsq` on \(\sqrt{w}\)-scaled rows). A design whose condition number after scaling exceeds \(10^{8}\) is reported as `ill_conditioned`, and that series is `NOT_TESTABLE`.

For each series, the weighted pairwise correlations between the columns \(1, F, F k, F k^2, F k^3\) at the observed \(k_i\) are reported, together with the condition number of each \(M_q\). They are diagnostics only: identifiability is decided empirically by G2 and G3, not by a correlation threshold. They are computed from the curve and the experimental \(k\) only.

## 7. Decision rule (fixed; identical in the gate and in the real analysis)

Per series:
1. Holm correction over \(\{p_1,p_2,p_3\}\) at family-wise \(\alpha=0.01\).
2. No rejection → `NO_DETECTION`.
3. Otherwise `DETECTED`, with selected shape \(q^*\) = the rejected \(q\) with the smallest raw \(p_q\) (ties → smaller \(q\)), and sign \(s^*=\mathrm{sign}(\hat\gamma_{q^*})\) (three-valued; exact zero counts as no sign).
4. `DETECTED` with \(q^*=2\) is the only outcome that refers to the target form. \(q^*\in\{1,3\}\) is reported as "excess residual of a non-target shape".

**Scope of error control.** The Holm family is defined separately within each series. The nominal family-wise error rate of 0.01 applies per series, not jointly across both series. No combined project-wide detection claim is made from the two series.

## 8. Synthetic power/identifiability gate

### 8.1 Synthetic truth
- Truth curve: the prepared digitized GPE curve of the same series, linearly interpolated (§3–§4). Curve-only.
- Experimental \(k\): the observed `unified_v1` \(k_i\) of the series. Observed experimental \(f_i\) are **not read**.
- True experimental frequency: \(f_i^{\rm true}=F(k_i)+A\,T_q(k_i)\).

### 8.2 Noise injection, per replicate

Noise is injected in pixel space, the space in which the digitization errors occur.

1. **To pixels.** Map the true experimental points \((k_i,f_i^{\rm true})\) and the curve nodes to pixel coordinates with the archived calibration (R0 §4.3, full precision from `wpd.json`).
2. **Click error.** Add independent \(\mathcal N(0,\delta^2)\) px to both pixel coordinates of every experimental point and every curve node.
3. **Calibration jitter.** Add independent \(\mathcal N(0,1^2)\) px to both pixel coordinates of each of the four calibration points X1, X2, Y1, Y2. The jitter is fixed at 1 px for all \(\delta\): calibration points are placed on axis ticks and were refined with Tweak Calibration.
4. **Back to values.** Convert all perturbed pixel coordinates to \((k,f)\) with the jittered calibration. The same jittered map is used for the experimental points and for the curve, as in the real data.
5. Re-prepare and re-interpolate the curve exactly as in §3–§4, recompute \(s_i\), \(t_i\) and \(w_i\), and apply §6–§7 unchanged.

The pixel↔value map must reproduce WebPlotDigitizer's `XYAxes` conversion. Before any replicate, the code checks that it reproduces the exported `gpe_large.csv` and `gpe_small.csv` values from their `px`, `py` columns within the mixed tolerance
\[
|x_{\rm rec}-x_{\rm exp}| \le \epsilon_{\rm abs}+\epsilon_{\rm rel}\,|x_{\rm exp}|,
\]
applied separately to \(k\) and \(f\), with \(\epsilon_{\rm rel}=10^{-9}\), \(\epsilon_{\rm abs}=10^{-9}\ \mu{\rm m}^{-1}\) for \(k\) and \(\epsilon_{\rm abs}=10^{-9}\) kHz for \(f\). These are about \(10^{-7}\) px, far below any physical scale. Recorded fact on the inputs: the `k_um_inv` and `omega_over_2pi_kHz` values in `gpe_large.csv` and `gpe_small.csv` are written with 15–17 significant digits (counted on the text only, no values compared), so export rounding is about \(10^{-15}\) relative and cannot by itself cause the check to fail. This check uses curve values only. If it fails, the run stops and the failure is recorded; the tolerances are not loosened without a new spec version (§11).

Random draws within one replicate are taken in this fixed order: calibration jitter (8 values, X1x, X1y, X2x, X2y, Y1x, Y1y, Y2x, Y2y); experimental points (\(x\) then \(y\), in file order); curve nodes (\(x\) then \(y\), in file order).

Base-term content is irrelevant: the partial F statistic is invariant to any \(\beta_0B_0+\beta_1B_1\) added to the truth, so none is injected.

### 8.3 Effect-size grid
Amplitude \(A\) is set so that \(\max_i |A\,T_q(k_i)| = m\cdot\delta\cdot a_f\), with \(m\in\{2,3,5,8\}\), both signs, for \(q\in\{1,2,3\}\). \(m\) is in units of the per-point frequency noise.

### 8.4 Replicates and seeds
- Null (\(A=0\)): 10,000 replicates per (series, \(\delta\)).
- Signal: 500 replicates per (series, \(\delta\), \(q\), \(m\), sign).
- Total: \(2\times4\times(10{,}000+3\times4\times2\times500)=176{,}000\) replicates; each replicate fits \(M_B\) and \(M_1,M_2,M_3\).
- Seeds: `numpy.random.SeedSequence(20261125)`, `spawn` in this fixed nested order: series (large, small) → \(\delta\) (1, 2, 4, 8) → case (null; then \(q\)=1,2,3 × \(m\)=2,3,5,8 × sign −,+). One child per case; each replicate draws from that child's `default_rng` sequentially.

### 8.5 Gate criteria, per series

CP = exact two-sided 95% Clopper–Pearson interval, from beta quantiles (`scipy.stats.beta.ppf(0.025, x, n-x+1)` and `beta.ppf(0.975, x+1, n-x)`, with 0 and 1 at the edges). No normal approximation.

- **G1 false positives.** Under the null, at every \(\delta\in\{1,2,4\}\), the CP interval of the `DETECTED` rate must satisfy both:
  - lower bound ≤ 0.010 — no statistical basis to conclude the rate exceeds the nominal 1%;
  - upper bound ≤ 0.015 — the data exclude a rate above 1.5%.

  With 10,000 replicates, a rate of exactly 1% gives about [0.0081, 0.0121]; the upper condition allows at most 126 rejections. G1 shows compatibility with the nominal level and bounds an excess. It does not show that the true rate equals 0.01.
- **G2 target power.** For \(q=2\), there is a smallest \(m^*\in\{2,3,5,8\}\) such that, for both signs and every \(\delta\in\{1,2,4\}\):
  - CP lower bound of `DETECTED` rate ≥ 0.80;
  - among detections, CP lower bound of \(q^*=2\) ≥ 0.80;
  - among detections, CP lower bound of correct sign ≥ 0.98.
- **G3 confusion (reported, gating).** For true \(q=1\) and \(q=3\) at \(m^*\), the rate of \(q^*=2\) among detections must have CP upper bound ≤ 0.20. Otherwise target-shape claims are not allowed for that series; only sign-of-excess claims are.
  - G3 is evaluated for each true \(q\in\{1,3\}\), each sign and each \(\delta\in\{1,2,4\}\); it passes only if every such cell passes.
  - If a cell has at least one `DETECTED`, the CP interval is computed on the detections in that cell.
  - If a cell has zero `DETECTED`, its G3 status is `NOT_ESTIMABLE`. `NOT_ESTIMABLE` is not a pass: absence of detections is not evidence of good identifiability. A series with any `NOT_ESTIMABLE` G3 cell can be at most `SIGN_ONLY`.
- \(\delta=8\) results are reported but do not gate.

**Gate verdict per series:**
- `TESTABLE(m*)` — G1, G2 and G3 pass;
- `SIGN_ONLY(m*)` — G1 and G2 pass, G3 fails or is `NOT_ESTIMABLE`;
- `NOT_TESTABLE` — G1 or G2 fails, or `ill_conditioned`.

A `NOT_TESTABLE` small series is the expected possible outcome with \(N=9\). It is reported as "not testable with available data", not as a failure of the method. The two series are never pooled.

### 8.6 Gate outputs
Per (series, \(\delta\), case): counts of `DETECTED`, \(q^*\) histogram, sign counts, CP intervals; the §3 preparation counts; §4.3 leave-one-out RMS; condition numbers and the §6 design correlations (at the observed \(k\)); the result of the §8.2 conversion check; summary of the calibration jitter actually drawn (mean and SD per coordinate); code and input hashes. Written to a results file before anything is printed.

## 9. Real-data analysis (after freeze only)

### 9.1 Primary
Apply §3–§7 once to each series that is `TESTABLE` or `SIGN_ONLY`. A `NOT_TESTABLE` series is not analysed for detection; its residual may be plotted only as a descriptive figure labelled "not testable".

### 9.2 Reported quantities
Outcome (`NO_DETECTION` / `DETECTED`), \(q^*\), sign, the three raw and Holm-adjusted p-values, \(\hat\sigma\) in kHz and in pixels, and \(N\). No \(\hat\gamma\) magnitude is converted to \(\Delta\Lambda\).

### 9.3 Sensitivity checks (do not change the primary outcome)
- PCHIP instead of linear;
- unweighted (OLS) instead of §5.2 weights.
If either changes the outcome, \(q^*\) or sign, the result is labelled `interpolation-sensitive` or `weight-sensitive` and must not be stated as a detection without that label.

### 9.4 Consistency with the gate
If \(\hat\sigma/a_f > 8\) px, the result is outside the validated noise range and is labelled `outside_gate_range`.

## 10. Allowed statements

- `DETECTED`, \(q^*=2\), series `TESTABLE`: "An excess residual of the \(F(k)k^2\) form with sign \(s^*\) is distinguishable from zero using a Holm procedure at nominal family-wise \(\alpha=0.01\), relative to the digitized published 2D GPE curve, under the stated uncertainty assumptions. The empirical null rejection rate and its two-sided 95% Clopper–Pearson interval are reported separately."
- `DETECTED` in a `SIGN_ONLY` series: sign of excess only, no shape statement.
- `NO_DETECTION`: "no excess residual distinguishable from zero; the gate sensitivity was \(m^*\) noise units". Not an upper limit on \(\Delta\Lambda\).
- Never: a value or interval for \(\Delta\Lambda\) or \(\Lambda\) (B2 V4, R0 §9); "evidence for new physics"; any statement based on experiment-minus-LDA.

Error control is per series (§7). Results for the large and small series are reported separately; no statement combines them into a single detection claim.

Error-rate wording, used verbatim with the gate result of the series: "Nominal family-wise error rate: 0.01. Empirical null rejection rate: \(\hat r\), with its two-sided 95% Clopper–Pearson confidence interval. G1 requires the interval lower bound to be ≤ 0.010 and the upper bound to be ≤ 0.015." It is never stated that G1 proves the family-wise error rate equals 0.01.

Every statement carries: digitized data from a published figure; SEM not available, absorbed into \(\sigma\); GPE curve is the published simulation, not the authors' solver output.

## 11. Deviations

Any change to §1–§10 after approval requires a new version of this spec with a recorded reason, made before the affected step is run. A change made after seeing gate output must be labelled `post-gate deviation`. A change after the real residual is opened makes the affected result exploratory.

## 12. Items closed or carried in R0 by this spec

- R0 §6 (gate design): specified here; execution pending.
- R0 §7 (error model): specified here (§5).
- R0 §8 (interpolation, pairing, domain, propagation, exclusion): specified here (§1, §3, §4).
- R0 §12 LDA line: proposed to mark complete as "generating script not recovered; provenance remains incomplete; table not used as an input".
- R0 §10 journal version: still open. The journal article is Phys. Rev. Lett. 109, 195301 (2012), DOI 10.1103/PhysRevLett.109.195301, title "Phonon dispersion relation of an atomic Bose-Einstein condensate". Its Figure 2 must be compared with arXiv v2 Figure 2. If it differs, the arXiv v2 figure remains the designated source unless a recorded decision says otherwise, and the difference is stated with every result.

## 13. Approval

- [ ] Reviewer approval of this spec.
- [ ] Dimitar Kretski approval of this spec.
- [ ] Code implementing §3–§8 reviewed (separate step).
- [ ] Permission for the single synthetic gate run.
