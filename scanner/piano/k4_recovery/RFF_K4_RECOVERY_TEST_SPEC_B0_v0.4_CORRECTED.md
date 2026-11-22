# B0: Blind recovery of a k⁴-equivalent coefficient — test specification v0.4 (consolidated)

**Status:** consolidated final draft, for the internal-consistency review before freezing by SHA-256. It incorporates the v0.4 amendment (C2a/C2b) and the decision on its strictness. Specification only: no code written, nothing run.
**Date:** 2026-11-22.
**Basis:** audit A, `RFF_K4_RECOVERY_AUDIT_v0.1.md` (sha256 eb7cb961…), on residual-first-framework **e7c7520**.
**Location (decided):** `residual-first-framework`, next to the piano scanner. Code goes in `scanner/piano/k4_recovery/`. Only the verdict and its scope are later entered in the Λ map.
**Seed (decided, final):** `20261122`. It is not changed after the smoke run or after any review of outputs.

**Changes from v0.3 (v0.4 amendment and its review, 2026-11-22):**
- **C2 split into C2a and C2b.**
  - **C2a:** sign of the winning deviation a_{p*}, over all scored P and S cells. **Diagnostic only.**
  - **C2b:** sign of **B̂**, over the **target set** (P n³ and S, scored), with 500 cases per sign. This is the **primary sign gate**.
  - Reason: in v0.3 the sign gate counted n² and n⁴ cells, whose winning-coefficient sign is not the sign of B.
- §7 ladder: C1 → **C2b** → (C3a, C3c) → C3b → C4 → C5. V1, V2 and V5 are reworded for the sign of B.
- §10: the target sample of 500 per sign is accepted knowingly (variant 1). It is not raised to 900, and the stricter effective demand is recorded.
- Appendix B:
  - C2b rules and operating characteristics for 500 and 400 cases per sign, plus thresholds for every smaller target size;
  - C2a's v0.3 figures kept and labelled diagnostic.
- A pre-declared consequence: below 183 detected cases a sign cannot pass C2b; N = 100 is explicitly marked unattainable.
- V3 has a separately defined four-cell C4-PURE rule; standard six-cell C4 is unchanged.
- **No numerical threshold changed.** Seed, levels, generator and manifest are unchanged.

**Changes from v0.2 (review of 2026-11-22, second round):**
- C2: the denominator for sign accuracy is the number of **detected** cases for that sign, not the 900 cases;
- Appendix B: C2 recomputed for that rule, together with the detection criterion, per sign and for both signs jointly;
- §5.1 added: what PASS and FAIL mean;
- §10: all decisions resolved, including keeping C2 and C3a strict. **No threshold changed.**

**Changes from v0.1 (review of 2026-11-22, first round):**
- §7 rebuilt as an ordered ladder; a C3b outcome added;
- W2 defined exactly (null hypothesis, both-significant case);
- W3 defined exactly (parametrisation, delta method), with mis-specification stated;
- §6 separates technical from statistical independence; per-case parameters recorded; counts asserted before analysis;
- smoke run defined as technical only;
- C3b and C3c interpretations narrowed;
- the binomial justification is now computed for the exact acceptance rules (Appendix B). **No threshold was changed.**

**Purpose:** give a final, checkable answer for **one** pipeline to the question:

> Can the method detect a deviation, decide that it has the k⁴-equivalent form, and recover its coefficient B — including **the sign of B when that sign is unknown** — with a calibrated interval?

**Phases:**
- **B0** (this document) fixes the questions, the criteria and the stopping rule.
- **B1** writes the minimal code (a wrapper; the original scanner is not changed).
- **B2** is one run.

Nothing in B1, the smoke run or B2 starts without separate approval.

---

## §1 What this test can and cannot establish

- **Within a closed family.** The test checks recognition of n³ (the k⁴-equivalent) against the fixed competitors {n⁰, n², n⁴}, plus one out-of-family signal and one contaminated signal. It does **not** test non-parametric discovery of an arbitrary dispersion law. Any positive verdict is worded as "within the family {0, 2, 3, 4}".
- **Synthetic only.** No real recording is used, so the real piano, GLOBEX and Black Sea results are not touched, re-used or re-interpreted.
- **One pipeline.** A verdict applies to the piano pipeline at config v0.4.1 plus the B1 wrapper, not to "the residual-first method" in general.
- **One generator.** All cases come from one synthetic generator. They are Monte Carlo replicates, not independent physical experiments (§6.3).
- **Injection recovery is not physical extraction.** A pass shows the estimator behaves as stated under the simulated conditions. It does not show that a real residual is free of systematics: GLOBEX's opposite-sign reflection false alarms are the reason.

---

## §2 Choice of pipeline: piano v0.4.1

**Why the piano pipeline:**
- it is the only one whose real-data analysis selected the k⁴-equivalent shape;
- its generator (`synth_note`) already produces both the pure power law and the exact stiff-string law;
- it needs no external data;
- the configuration that touched real data (5bdf34b190a4b9ef) is exactly the one the audit found **unvalidated** (audit F2, A2).

GLOBEX and Black Sea are not included.

**Frozen inputs (verified by the wrapper at start-up; abort on mismatch):**
- `scanner/piano/piano_scanner.py` sha256 `335ba8b8abaef4b424b36a516f55df568b9345f83cc8b6bc1c0dc1cc33c9e5d4` (e7c7520);
- `piano_scanner.config_hash()` = `5bdf34b190a4b9ef`.

`piano_scanner.py`, its CONFIG, templates and thresholds are **not modified**. B1 is an external layer only.

---

## §3 Minimal code changes (B1): a wrapper, the scanner untouched

**Files:**
- `scanner/piano/k4_recovery/k4_recovery_b1.py` (generator and analysis);
- `scanner/piano/k4_recovery/k4_score_b1.py` (scorer, the only code that reads the truth).

The wrapper imports `piano_scanner` and reuses `analysis_segment`, `track_partials`, `fit_templates`, `analyze_signal`, `synth_note`, `random_f0s` and `detune` unchanged.

**Notation (as in `fit_templates`, l. 295–329):**
- x = n / n_max;
- w = SNR / median(SNR);
- the two-column model for template p is [x, x^p], fitted by weighted least squares on f;
- RSS_p is its weighted residual sum of squares;
- N is the number of fitted partials.

### W1 Sign record

Per case, record sign(a_{p*}) of the winner (used by C2a) and sign(B̂) (used by C2b and C5).

### W2 Discrimination statistic (exact copy of the GLOBEX rule, `globex_scanner.run_test` l. 306–316, adapted to the piano design)

1. **Winner:** p* = argmin_p RSS_p over p ∈ {0, 2, 3, 4}. This is the original winner (l. 328), taken from `fit_templates`.
2. **For each competitor q ≠ p*:** fit the three-column model M(p*, q) = [x, x^{p*}, x^q] with the same x and w, giving RSS_{p*q}.
3. **Test T(p* | q):**
   - H₀: the coefficient of x^{p*} in M(p*, q) is zero, i.e. the model [x, x^q] is adequate.
   - F = (RSS_q − RSS_{p*q}) / (RSS_{p*q} / (N − 3)), with df (1, N − 3). If RSS_{p*q} ≤ 0 or N − 3 ≤ 0, F = ∞ and p = 0, as in GLOBEX `_ftest`.
4. **Holm** over the three competitors. p* is **distinguishable** if every Holm-adjusted p < 0.01.
5. **Meaning:** "distinguishable" means **p\* cannot be replaced by any single competitor**. It does not mean that q is absent, and it is not a comparison of two alternative models.
6. **Both-significant case:** the reverse tests T(q | p*) are also computed and recorded. A case where T(p* | q) and T(q | p*) are both significant for some q is flagged `mixture_like`.
   - The flag does **not** change the declared outcome; this keeps the GLOBEX definition exact.
   - It is reported per cell.

### W3 Coefficient and interval

- **Fit:** the p = 3 model with weighted design X_w = [x·w, x³·w], solved by least squares (exactly as `fit_templates`), giving coefficients (c₀, c₁).
- **Scale:** s² = RSS₃ / (N − 2).
- **Covariance:** Cov(c) = s² (X_wᵀ X_w)⁻¹.
- **Estimate:** a₃ = c₁ / n_max³ and f0 = c₀ / n_max, so **B̂ = 2a₃/f0 = 2c₁ / (c₀ n_max²)**. This is identical to `B_est` (l. 362).
- **Delta method:** gradient g = (−2c₁ / (c₀² n_max²), 2 / (c₀ n_max²)), se = √(g Cov gᵀ).
- **Interval:** B̂ ± t₀.₉₇₅,N−2 · se.
- B̂ and its interval are computed for **every analysed case**, whatever the winner. Sign and amplitude recovery are therefore not conditional on the shape decision.

**Assumptions:**
- the partial errors are independent with σ_n ∝ 1/SNR up to one common factor (the scanner's own assumption, l. 297–299);
- the linear p = 3 model is correct.

**Mis-specification, stated in advance:**
- **Pure n³ cases (block P):** the model is correct, so the empirical coverage there tests calibration.
- **Stiff-string cases (S):** the law has higher-order terms (≈ −B²n⁵/8 and beyond), so the p = 3 model is **mis-specified**. Coverage there is the **empirical coverage of the p = 3 interval under the stiff law**, measured against the injected B. It is never presented as calibrated by construction.
- The same applies to blocks O, C and R.

**Checks in B1 (static review and smoke run, never on B2 data):**
- wrapper c₀, c₁ reproduce `per[3]["f0"]`, `per[3]["a"]` and `B_est` (relative 1e-12);
- se agrees with an independent computation (`scipy.optimize.curve_fit` on the unscaled model with σ = 1/w and `absolute_sigma=False`), relative tolerance 1e-6.

### W4 Declared outcome per case

`none` / `n^p± distinguishable` / `n^p± undetermined`:
- "significant" is exactly v0.4.1's `significant` (Holm over 4 templates, α = 0.01, |cents10| ≥ 0.1);
- "distinguishable" is W2;
- the sign is that of a_{p*}.

### W5 Sealed generation

- **New generators:**
  - `mixture`: an n³ term plus an n² term, built as in `synth_note`'s `power` branch;
  - pure n⁵ via `synth_note(model="power", p=5)`, which is already supported.
- **Case order** is shuffled.
- **The truth** goes to a separate file that only the scorer reads (§6).

### W6 Manifest, seeds and counts

- **Manifest first:** before any signal is generated, the full case manifest (block, cell, sign, level, m) is built and its counts are **asserted against the §4 table**. On mismatch, abort before analysis.
- **Per-case seeds:**
  - each case gets its own child seed from `numpy.random.SeedSequence(20261122).spawn(5901)`; index 5900 is reserved for the shuffle;
  - each case's f0 draw, detuning, phases, partial amplitudes, jitter and noise come from its own stream.
- **Recorded per case:** case id, spawn index, block, cell, sign, level, f0, detuned f_nom, n_partials, n_fit, exclusion reason.
- **Scorer checks:** it asserts unique spawn indices and the per-cell counts before scoring.

### Parity guard

- For every case, the original `analyze_signal` is also called on the same signal.
- Equal `winner`, `significant` and per-template p-values are asserted (relative tolerance 1e-12).
- A mismatch is a wrapper defect and **aborts the run**.

---

## §4 Cases (all synthetic; f0 log-uniform 27.5–1100 Hz as `random_f0s`; nominal f0 detuned as `detune`)

Amplitudes are given as the deviation of partial 10 in cents, d₁₀.

| block | truth | levels | sign | m per cell | cells | role |
|---|---|---|---|---|---|---|
| N | harmonic | — | — | 2000 | 1 | null FPR |
| P | pure power n^p, p ∈ {2, 3, 4} | d₁₀ ∈ {0.5, 2, 5} cents | ± | 100 | 18 | sign, shape, amplitude (n³) |
| S | exact stiff string, f_n = n f0 √(1+Bn²) | \|B\| ∈ {1e-5, 1e-4, 4e-4, 1e-3} (d₁₀ ≈ 0.9, 8.6, 34, 165 cents for B > 0) | ± | 100 | 8 | sign, shape, amplitude for the physical law |
| O | out of family: pure n⁵ | d₁₀ = 2 cents | ± | 100 | 2 | must not be claimed as n³ |
| C | contamination: n³ (d₁₀ = 5 cents) plus an n² term of **opposite sign** and half the n = 10 deviation | — | ± (sign of the n³ part) | 100 | 2 | a second mechanism of opposite sign (GLOBEX lesson) |
| R | robustness stratum (jitter 0.3 cent per partial, noise −40 dB): null; n³ d₁₀ = 5 cents; stiff \|B\| = 4e-4 | — | ± | 500 null, 100 others | 5 | degradation under worse measurement |

**Total:** 2000 + 1800 + 800 + 200 + 200 + 900 = 5900 analyses in 36 cells. The other generator parameters are the `synth_note` defaults.

**Notes:**
- The S grid is the scanner's own pre-real-data default (`--B`, set at v0.3, before any piano note was analysed). The cents grid lies inside the default `--cents` grid. **No level is chosen from the real piano results.**
- For B < 0, partials with 1 + Bn² ≤ 0 do not exist. For |B| = 1e-3 that is n > 31, above `fit_max_n` = 25.
- **Scored cells per sign (C2a, diagnostic):** P n², n³, n⁴ at 2 and 5 cents (6 cells) plus S at 1e-4, 4e-4 and 1e-3 (3 cells), which gives 900 cases per sign.
- **Target set per sign (C2b, primary):** P n³ at 2 and 5 cents (2 cells, 200 cases) plus S at 1e-4, 4e-4 and 1e-3 (3 cells, 300 cases), which gives **500 cases per sign**. These are exactly the scored cells whose injected law has a k⁴-equivalent n³ coefficient.
- **Sensitivity cells:** S at 1e-5 (0.86 cent) and P at 0.5 cent.

---

## §5 Criteria (fixed; justified from design and statistics, not from real data; operating characteristics in Appendix B)

**Conventions:**
- CP is the Clopper–Pearson 95 % interval.
- A cell is **scored** if |d₁₀| ≥ 2 cents (≥ 20× the scanner's minimum effect of 0.1 cent). Lower levels are sensitivity cells: reported, not scored.
- A cell is **evaluable** if ≥ 80 % of its cases pass the scanner's exclusion rule. Non-evaluable cells are reported and not scored.

| ID | criterion | pass condition | basis |
|---|---|---|---|
| **C1 null** | FPR = fraction of block-N cases declared anything other than `none` | CP upper ≤ 0.02 (= 2α) | nominal α = 0.01 after Holm; 2α allows the F approximation under SNR weights. This turns the old "not rejected" barrier into a measured FPR |
| **C2a sign of the winning deviation (diagnostic only)** | pooled over the scored, evaluable P and S cells, **separately for + and −**. Detection: n_analysed = cases that pass the exclusion rule (≤ 900 per sign); detected = declared ≠ `none`. Sign accuracy: n = n_detected for that sign; correct = sign of a_{p*} equals the injected sign of the deviation | reported against the same bars: detection CP lower ≥ 0.95, sign-correct CP lower ≥ 0.98, \|acc₊ − acc₋\| ≤ 0.02. **It does not set the verdict** | it measures whether the winning template's coefficient follows the injected deviation across the family. It is **not** evidence that the sign of B was recovered |
| **C2b sign of B (primary gate)** | the **target set** (§4): the scored, evaluable cells P n³ (2, 5 cents) and S (\|B\| = 1e-4, 4e-4, 1e-3), **separately for + and −** (nominal 500 each). (1) n_analysed = target cases that pass the exclusion rule; (2) detected = the original v0.4.1 outcome is significant (declared ≠ `none`), **whatever the winning template**; (3) n_B = n_detected; undetected cases carry no inferred sign; (4) correct = sign(B̂) equals the injected sign of B, with B̂ = 2c₁/(c₀ n_max²) from W3 | for **each** sign: detection CP lower ≥ 0.95 on n_analysed **and** sign-correct CP lower ≥ 0.98 on n_B. C2b passes only if **both** signs pass. \|acc₊ − acc₋\| ≤ 0.02 is reported (it is implied by the two bounds; Appendix B) | S3 for the target coefficient. Each sign is evaluated with **its own actual denominators**. If a target cell is non-evaluable for one sign, the two signs have different sizes; they are **not pooled**, and the missing cell is **not replaced** by other data. Thresholds for every possible size are fixed in Appendix B. If n_B < 183 for a sign, CP lower ≥ 0.98 is unattainable and C2b fails by construction |
| **C3a shape, pure** | per p ∈ {2, 3, 4}, pooled over sign and scored levels (400): `n^p distinguishable` with the correct p | correct CP lower ≥ 0.80; wrong p distinguishable CP upper ≤ 0.02 | F2, F3, within the family |
| **C3b shape, stiff law** | scored S cells pooled (600) | wrong p distinguishable CP upper ≤ 0.02; n³ distinguishable CP lower ≥ 0.50 | **a minimum practical criterion, not a statistically derived bound.** A pass means recognition in at least half of the **pooled** evaluable stiff-law cases, not at every value of B. Per-\|B\| rates are reported, but no per-value claim is made from them |
| **C3c out of family** | block O (200): declared `n³± distinguishable` | CP upper ≤ 0.05 | **a limited test:** it checks only that an n⁵ residual is not claimed as k⁴. Whether n⁵ is declared n² or n⁴ is reported, not scored. This is the closed-family limitation of §1 |
| **C4 amplitude** | per scored evaluable cell of P (n³ only) and S: bias = median(B̂/B) − 1; coverage = fraction of 95 % intervals containing B | **C4-core (standard):** P n³ ±2, ±5 cents and S ±1e-4 (six cells): \|bias\| ≤ 0.10 in each cell; coverage ≥ 0.85 in each cell; pooled coverage ∈ [0.92, 0.99] (552–594/600). **C4-PURE (V3 branch only):** the four P n³ cells (±2, ±5 cents): same per-cell bias and coverage limits; pooled coverage ∈ [0.92, 0.99] (368–396/400). **C4-range:** under standard C4, the largest \|B\| for which the per-cell conditions hold for both signs is reported; under C4-PURE, no stiff-law range is claimed | For S the coverage is empirical under a mis-specified model (W3). The exact joint binomial operating characteristics for both rules are in Appendix B. The bias threshold has no binomial derivation: its false-fail probability depends on the estimator's spread, which is unknown before the run, and this is stated, not estimated |
| **C5 robustness** | block C: correct sign of B̂; block R: C1–C4 rates recomputed | **C5-sign:** block C sign-correct CP lower ≥ 0.90 for each sign. The rest of C and R is reported | only a sign flip by an opposite-sign second mechanism is disqualifying. Bias under a mis-specified model is expected and reported |

All thresholds are written before any B1 code exists. None refers to the real piano values (B 1.1e-4 … 4.9e-4, n³ in 39/40).

### §5.1 What PASS and FAIL mean

- **PASS:** the pipeline met the pre-registered, deliberately strict criteria.
- **FAIL:** it did not meet them. **A FAIL does not by itself show that the method is useless.** It may mean that its accuracy or power is insufficient for the task as stated (for example, with detection 0.99 and true B-sign accuracy 0.99 in both signs, C2b fails with probability ≈ 0.93; see Appendix B).
- Each verdict is reported with the measured rates and their CP intervals, so the size of any shortfall is visible.
- Neither outcome is extended beyond the piano pipeline v0.4.1 + B1 (§1, §8).

---

## §6 Hidden truth, confirmation discipline and independence

### §6.1 Order

1. B0 approved.
2. B1 code written.
3. Static review of B1.
4. Smoke run (separate approval).
5. B0 and B1 frozen by hash.
6. **One** B2 run (separate approval).
7. Scoring.

### §6.2 Technical separation (hidden truth)

- The generator creates the signal and writes the truth row to `truth.jsonl`.
- The analysis receives only `(x, sr, detuned f_nom)` and writes `analysis.jsonl` keyed by case id.
- Case order is shuffled.
- Both files and `run.log` are hashed **before** the scorer is executed.
- The scorer is frozen in B1 and is the only code that joins them.

### §6.3 Statistical independence (separate from §6.2)

- All 5900 cases come from **one generator and one pipeline**. They are Monte Carlo replicates, not independent physical experiments.
- The CP intervals in §5 treat the cases of a cell as independent, identically distributed draws. This is justified only because each case has its own child seed (W6), and therefore its own f0, detuning, phases, amplitudes, jitter and noise.
- The intervals quantify Monte Carlo uncertainty about rates **under this generator's distribution**: f0 range, noise model, partial amplitudes and decay of `synth_note`. Any generalisation beyond that distribution is outside the test.
- Shared factors (one generator code, one analysis code) are common to all cases by design. They are not a source of within-cell dependence, but they limit generality.

### §6.4 Smoke run (requires its own approval)

- **Purpose:** a **technical check only. It is not statistical calibration.**
- **Settings:** seed 1 (not 20261122); 3 cases per cell (108 cases).
- **Checks:**
  - completion;
  - the parity guard;
  - the W3 covariance cross-check;
  - output schema;
  - the W6 count assertions;
  - that the scorer runs.
- **Scorer mode:** `--smoke`. It prints only structural checks (counts, schema, hashes) and **no rates**.
- **Use of the outputs:** they are not used to change thresholds, levels, generators or parameters. A code defect found there is fixed as a documented B1 revision before freezing.

### §6.5 No second run

There is no second B2 under the same B0. Thresholds, levels and the seed cannot change after the run.

---

## §7 Verdict (ordered ladder; unambiguous)

**Evaluation order:** C1 → **C2b** → (C3a and C3c) → C3b → C4 → C5. C2a is always reported as a diagnostic and never sets the verdict.
- The verdict is set by the **first level that fails**.
- All later levels are still computed and appended as **qualifiers**; they never change the verdict.
- If no level fails, the verdict is V5.

| first failing level | verdict |
|---|---|
| C1 | **V0 — NOT VALIDATED (null control):** the false-alarm rate is not controlled |
| C2b | **V1 — NOT VALIDATED (target coefficient sign):** the pipeline has not demonstrated reliable blind recovery of the sign of B in the target n³ / stiff-string cases |
| C3a or C3c | **V2 — SIGN OF B ONLY:** deviations are detected and the sign of B is recovered blind in the target cases, but the n³ shape is not validated; no shape claim is allowed |
| C3b | **V3 — PURE-POWER SHAPE ONLY:** n³ recognition is validated within the family for the pure power law, **not** for the exact stiff-string law. In this branch only **C4-PURE** is evaluated: the four P n³ cells (±2, ±5 cents), with per-cell \|bias\| ≤ 0.10 and coverage ≥ 0.85, plus pooled coverage 368–396/400. The six-cell standard C4-core and C4-range are not claimed in this branch. The suffix "+ quantitative (pure n³)" is added if C4-PURE passes, or "not quantitative" if it fails |
| C4 | **V4 — DETECTION, SIGN OF B AND SHAPE; NOT QUANTITATIVE:** B may be quoted as sign and shape only, not as a value with an interval |
| none | **V5 — EXTRACTION VALIDATED (synthetic, within the family {0, 2, 3, 4}):** detects, recovers the sign of B blind, recognises n³ without confidently mislabelling the tested out-of-family residual as n³, with the empirically measured interval coverage over the reported \|B\| range |

**Qualifiers, always appended:**
- **C2a:** the rates and CP intervals of the winning-deviation sign (diagnostic).
- **C5-sign failed:** "an opposite-sign second mechanism can flip the inferred sign". This becomes a standing caveat for any real-data use.
- **C3b passed:** "n³ recognised in ≥ 50 % of pooled stiff-law cases; per-|B| rates reported; no per-value claim".
- **C3c:** the declared shapes for n⁵ (n², n⁴, undetermined).
- The `mixture_like` rate per cell (W2).
- The C4-range.
- Non-evaluable cells.
- The R-stratum rates.

**Stopping rule:**
- At most **one** revision B0 → v0.5 after B2, and only if the failure is traced to a **code defect** in the wrapper, shown by a test independent of the B2 results. Threshold changes are never a valid reason.
- After that, the verdict is final for this pipeline and is entered in the Λ map as the status of "independent extraction method, Module 0, piano pipeline".

---

## §8 What a pass would and would not allow next

- **Would allow (separately proposed, not part of B):** a blind application of the frozen wrapper to the **untouched `pp` and `ff` piano recordings** reserved by the piano README. That is a real-data test with a known sign, where agreement between dynamics can be checked.
- **Would not allow:**
  - any statement about Λ in a physical sector;
  - any extension to GW;
  - citing the piano B as a measurement;
  - any claim that "the residual-first method" in general is validated.

---

## §9 Effort and runtime

- **B1:** about 300–400 lines (wrapper with W1–W6 and the parity guard, plus the scorer).
- **B2:** 5900 analyses, plus 5900 parity calls, on Lenovo (CPU). The per-note time of v0.4.1 is unknown and is logged.
- **Inputs hashed:** `piano_scanner.py`, the wrapper, the scorer, B0. **Outputs hashed:** `manifest.json`, `truth.jsonl`, `analysis.jsonl`, `run.log`, `score.json`.

---

## §10 Decisions

**Resolved by the review:**
- location: `residual-first-framework`;
- seed: 20261122, final;
- smoke run: yes, technical only, separate approval;
- C3b = 0.50: kept as a practical criterion.

**Resolved by the second review (v0.3; "C2" there is the rule now called C2a, which is diagnostic. For the primary gate, see the next block):**
- **C2** (sign CP lower ≥ 0.98) and **C3a** (wrong-shape CP upper ≤ 0.02) are **kept strict**, knowingly. The effective demands are sign accuracy ≈ 0.995 per sign and a wrong-shape rate ≈ 0.002 (Appendix B).
- **Reason:** the risk of not validating an adequate method is preferred to the risk of validating one that sometimes gets the sign or the shape wrong.
- The sample sizes are not increased, because that would change the pre-registered experiment.
- **The only correction then:** the C2 (now C2a) sign denominator is n_detected per sign. The same denominator rule is carried into C2b.

**Resolved by the v0.4 amendment and its review (variant 1):**
- The primary sign gate is **C2b** (sign of B̂ on the target set). C2a is diagnostic only.
- **A target sample of 500 cases per sign is accepted.** It is not increased to 900, and C2b keeps detection CP lower ≥ 0.95 and sign CP lower ≥ 0.98.
- **The consequence, recorded as a consciously accepted risk and not an omission:**
  - the effective demand is B-sign accuracy of about 0.998–0.999 per sign;
  - at a true detection probability of 0.97, the detection part alone passes with probability ≈ 0.57 at 500 cases, even with perfect sign recovery (≈ 0.46 at 400).
- **Reason:** the experiment is not changed in order to raise the probability of a PASS. The required accuracy is known in advance instead.
- **Unequal sign groups:** if a target cell is non-evaluable, each sign is evaluated with its own denominators. The rules for 400, 300 and 200 cases are pre-written in Appendix B, so no rule is chosen after seeing which cell dropped. Losing three of five target cells leaves N = 200 for that sign, so C2b remains mathematically possible but requires at least 197 detections and zero sign errors. Losing four cells leaves N = 100; C2b is then impossible even with 100 detections and 100 correct signs. More generally, C2b fails by construction whenever fewer than 183 cases are detected.
- **The `pp`/`ff` recordings (§8):** they carry a known positive sign. They can test transfer, shape, amplitude and consistency between dynamics on real recordings, but **not** blind sign recovery. That belongs to a later physical channel.

**Nothing open.** After the internal-consistency review of v0.4, B0 is frozen by SHA-256. No criterion is changed because of B2 results.

---

## Appendix B: operating characteristics of the exact acceptance rules

**Method:**
- exact binomial computation, for each rule, of the probability of passing as a function of the true rate;
- CP bounds via the beta quantiles;
- cases assumed i.i.d. within a cell (§6.3) and all cells evaluable with the nominal n.

This is arithmetic on the rules, not a simulation of the scanner.

| rule | n | pass region | P(pass) at the true rate |
|---|---|---|---|
| C1: FPR CP upper ≤ 0.02 | 2000 | ≤ 27 false alarms | 0.005: 1.000 · 0.010: 0.948 · 0.0125: 0.70 · 0.015: 0.33 · 0.020: 0.018 |
| C2b detection (per sign): CP lower ≥ 0.95 | 500 | ≥ 485 detected (0.970) | see the C2b tables below |
| C2b B-sign (per sign): CP lower ≥ 0.98 **on n_B = n_detected** | 485 … 500 | **≤ 3 sign errors** for every n_B in range (485 → ≥ 482 correct; 500 → ≥ 497 correct) | see the C2b tables below |
| C2a detection (diagnostic) | 900 | ≥ 868 detected | 0.99: 1.000 · 0.97: 0.86 · 0.96: 0.28 · 0.95: 0.024 |
| C2a winning-deviation sign (diagnostic) | 868 … 900 | ≤ 9 sign errors | see the C2a table below |
| C3a correct: CP lower ≥ 0.80 | 400 | ≥ 336 | 0.95: 1.000 · 0.90: 1.000 · 0.85: 0.74 · 0.80: 0.024 |
| C3a wrong: CP upper ≤ 0.02 | 400 | ≤ 2 | 0: 1.000 · 0.002: 0.95 · 0.005: 0.68 · 0.01: 0.24 · 0.02: 0.013 |
| C3b correct: CP lower ≥ 0.50 | 600 | ≥ 325 | 0.65: 1.000 · 0.60: 0.998 · 0.55: 0.67 · 0.50: 0.023 |
| C3b wrong: CP upper ≤ 0.02 | 600 | ≤ 5 | 0.002: 0.999 · 0.005: 0.92 · 0.01: 0.44 · 0.02: 0.019 |
| C3c: CP upper ≤ 0.05 | 200 | ≤ 3 | 0: 1.000 · 0.01: 0.86 · 0.02: 0.43 · 0.03: 0.15 · 0.05: 0.009 |
| C5-sign: CP lower ≥ 0.90 | 100 | ≥ 96 | 0.999: 1.000 · 0.99: 0.997 · 0.98: 0.95 · 0.97: 0.82 · 0.95: 0.44 |
| C4 coverage, standard joint rule: each of 6 cells ≥ 85/100 **and** pooled in [552, 594]/600 | 6 × 100 | as stated | 0.99: 0.56 · 0.97: 0.9997 · **0.95: 0.9992** · 0.93: 0.85 · 0.92: 0.54 · 0.90: 0.055 · 0.85: 0.000 |
| C4-PURE coverage, V3 branch: each of 4 cells ≥ 85/100 **and** pooled in [368, 396]/400 | 4 × 100 | as stated | 0.99: 0.568 · 0.97: 0.998 · **0.95: 0.996** · 0.93: 0.812 · 0.92: 0.546 · 0.90: 0.103 · 0.85: 0.000 |

### C2b (primary gate)

**Model (per sign):**
- D ~ Binomial(N, d) detections among N analysed target cases;
- K | D ~ Binomial(D, s) correct B-signs;
- the sign passes if L_CP(D, N) ≥ 0.95 and L_CP(K, D) ≥ 0.98, where L_CP is the two-sided 95 % Clopper–Pearson lower bound;
- P(pass) = Σ_D P(D | N, d) · 1{L_CP(D, N) ≥ 0.95} · P(K ≥ k_min(D) | D, s).

**Thresholds for every target size (fixed now; N is the actual number of analysed target cases of that sign):**

| N (analysed, one sign) | target cells | detections needed | observed detection rate needed | sign errors allowed (for every D in the passing range) | smallest correct count |
|---|---|---|---|---|---|
| 500 | all 5 | ≥ 485 | 0.9700 | ≤ 3 | 482 at D = 485; 497 at D = 500 |
| 400 | 4 (one non-evaluable) | ≥ 389 | 0.9725 | ≤ 2 | 387 at D = 389; 398 at D = 400 |
| 300 | 3 | ≥ 293 | 0.9767 | ≤ 1 | 292 at D = 293; 299 at D = 300 |
| 200 | 2 | ≥ 197 | 0.9850 | 0 | 197 at D = 197; 200 at D = 200 |
| 100 | 1 | 100 detected for detection gate; sign gate cannot pass | 1.0000 | **unattainable**: even 100/100 correct gives L_CP < 0.98 | C2b fails by construction |
| < 183 detected | — | — | — | **unattainable**: even D/D correct gives L_CP < 0.98 | C2b fails by construction |

For an N that is not in the table (exclusions inside evaluable cells), the same two inequalities are applied to the actual counts. The table illustrates them; it does not replace them.

**Operating characteristics, N = 500 (one sign; "both" assumes equal d and s and independent sign groups):**

| true detection d | s = 0.999 | 0.998 | 0.997 | 0.995 | 0.99 | 0.985 |
|---|---|---|---|---|---|---|
| 0.999 | 0.998 | 0.981 | 0.935 | 0.758 | 0.264 | 0.058 |
| 0.99 | 0.998 (both 0.997) | 0.982 (0.964) | 0.936 (0.877) | 0.763 (0.582) | 0.271 (0.073) | 0.061 (0.004) |
| 0.98 | 0.951 | 0.936 | 0.894 | 0.732 | 0.264 | 0.061 |
| 0.97 | 0.567 | 0.558 | 0.534 | 0.438 | 0.160 | 0.037 |

**Operating characteristics, N = 400 (one sign):**

| true detection d | s = 0.999 | 0.998 | 0.997 | 0.995 | 0.99 | 0.985 |
|---|---|---|---|---|---|---|
| 0.999 | 0.992 | 0.953 | 0.880 | 0.677 | 0.237 | 0.061 |
| 0.99 | 0.991 | 0.953 | 0.882 | 0.682 | 0.242 | 0.063 |
| 0.98 | 0.884 | 0.850 | 0.787 | 0.611 | 0.220 | 0.058 |
| 0.97 | 0.456 | 0.439 | 0.407 | 0.317 | 0.115 | 0.031 |

**Detection part alone** (perfect sign): d = 0.97 passes with probability 0.568 (N = 500) and 0.460 (N = 400); d = 0.98 with 0.953 and 0.890.

**The difference criterion:** among passing outcomes the B-sign accuracy is ≥ 0.9938 (N = 500) and ≥ 0.9949 (N = 400). So \|acc₊ − acc₋\| ≤ 0.02 holds automatically whenever both signs pass, and it changes no probability.

**If the two signs have different N:** each sign's probability is read from its own row. The joint figure is reported only under the stated assumption of independent groups, as the product. Otherwise the two figures are reported separately.

The v0.3 figures for n = 900 do **not** apply to C2b and are not used to justify it.

### C2a (diagnostic only; the v0.3 rule, unchanged)

- **Model:** per sign, n_detected ~ Binomial(900, d); given n_detected, correct signs of a_{p*} ~ Binomial(n_detected, s).
- **Pass bars:** n_detected ≥ 868 and CP lower ≥ 0.98, which permits ≤ 9 errors.
- **Both signs:** P(both) = P(one)² under equal rates.
- **The difference criterion:** implied, since accuracy is ≥ 0.9896 among passing outcomes.

| true detection d | true sign accuracy s = 0.999 | 0.995 | 0.99 | 0.985 |
|---|---|---|---|---|
| 0.999 | one sign 1.000 · both 1.000 | 0.983 · 0.967 | 0.589 · 0.346 | 0.134 · 0.018 |
| 0.99 | 1.000 · 1.000 | 0.984 · 0.969 | 0.599 · 0.359 | 0.141 · 0.020 |
| 0.97 | 0.858 · 0.736 | 0.846 · 0.716 | 0.533 · 0.284 | 0.134 · 0.018 |
| 0.96 | 0.282 · 0.079 | 0.278 · 0.077 | 0.176 · 0.031 | 0.045 · 0.002 |

These C2a figures describe a diagnostic and do not enter the verdict.

**Reading:**
- **C4 standard:** a calibrated interval (true coverage 0.95) fails the six-cell **joint** coverage rule with probability ≈ 8.5e-4. Coverage of 0.90 fails with probability ≈ 0.95, and a strongly over-wide interval (0.99) fails with probability ≈ 0.44. This corrects the v0.1 text, which gave only per-condition bounds. (The v0.1 comparison is confirmed: a per-cell [0.90, 0.99] rule over six cells would fail a calibrated estimator with probability 0.10.)
- **C4-PURE (V3 only):** the exact joint-rule calculation is repeated for four independent Binomial(100, p) cells, each truncated to k ≥ 85, with pooled coverage 368–396/400. At true coverage 0.95 the pass probability is ≈ 0.996. This is a separate conditional rule, not a substitute for standard six-cell C4.
- **C1:** a scanner exactly at its nominal α passes with probability 0.95. At 1.5α it passes with probability 0.33.
- **C2b and C3a wrong:** kept strict by decision (§10). C2b's effective demand is a B-sign accuracy of about 0.998–0.999 per sign, with detection of about ≥ 0.98.

**Computation:** `scipy.stats.beta` and `scipy.stats.binom`. The C4 joint probability is the exact convolution of six truncated binomial(100, p) distributions (each truncated to k ≥ 85), summed over the pooled window.
