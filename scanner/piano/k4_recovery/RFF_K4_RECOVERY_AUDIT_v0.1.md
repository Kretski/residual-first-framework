# A: Can the existing residual-first scanners recover an unknown k⁴ coefficient? Audit of Module 0 at e7c7520

**Status:** draft for review. Document-only audit. No code from the audited repository was run; nothing was re-computed from the data. The only computations were: reading committed outputs, `git` history, SHA-256 of files, and re-deriving one config hash statically (AST literal of `CONFIG`, no execution of the scanner).
**Date:** 2026-11-22.
**Repository:** github.com/Kretski/residual-first-framework, **HEAD e7c7520d59299e0535fb101c4c13a8d6f377b9bd**, full history and tags fetched read-only.
**Scope:** Module 0 only (`scanner/piano`, `scanner/globex`, `scanner/waves`). Module 1 (GW) is out of scope; it has a fixed α = 4 template and no shape test, so it cannot answer the shape question in any case.

**Status labels (as requested):**
- **PASS:** shown by an existing, committed test.
- **PARTIAL:** shown only in a limited case (one sign, one config, synthetic only, or with an unverifiable step).
- **NOT TESTED:** no executed evidence in the repository.
- **FAIL:** an existing test shows a mismatch.
- **CODE:** a property of the code (for example "accepts a negative value"). This is **not evidence** of recovery and never counts as PASS.

---

## §0 The question, and the three sign claims kept apart

Main question: can the existing method detect a residual, establish that it has k⁴ dependence, and recover its coefficient with free sign and amplitude, without knowing in advance what it is looking for?

For sign, three different claims are scored separately:
- **S1:** the code accepts a coefficient of either sign (CODE);
- **S2:** the code has been shown to recover a coefficient of a given sign;
- **S3:** the method identifies the sign correctly **when the sign is unknown**.

Only S3 answers the main question. The same three-level split is used for shape (CODE family / identified in injections / identified against competitors and on real data) and for amplitude (estimate exists / bias measured / uncertainty and coverage measured).

---

## §1 What each environment measures, in canonical form

| environment | baseline subtracted | template family (fixed) | k⁴-equivalent template | canonical link |
|---|---|---|---|---|
| Piano (`piano_scanner.py`) | f_n = n·f0 (f0 fitted) | Δf_n ∝ n^p, p ∈ {0, 2, 3, 4} (l. 77) | p = 3 | stiff string f_n = n f0 √(1+Bn²) ⇒ ω² = c²k²(1+Λk²), k = nπ/L, Λ = B L²/π² |
| GLOBEX (`globex_scanner.py`) | ω = k√(gh)(1 + s_b), s_b free per block | r = s_b + a(kh)^q, q ∈ {1, 2, 3} (l. 99) | q = 2 | r ≈ ½Λk² ⇒ Λ = 2a h²; linear theory a ≈ −1/6 ⇒ Λ = −h²/3 |
| Black Sea (`waves_scanner.py`) | ω = √(g_eff k) + k·U | δω = a√(gk)(k/k_ref)^q, q ∈ {1, 2, 3} (l. 85) | q = 2 | ω² = gk(1 + Ck²): a relative k² term on a **gk** baseline, not c²k² |

Black Sea has the same relative k² structure but a different baseline. It is included for completeness and is the weakest analogue.

---

## §2 Summary matrix

P = piano, G = GLOBEX, W = Black Sea.

| # | claim | P | G | W |
|---|---|---|---|---|
| **1** | **Sign** | | | |
| S1 | estimator accepts either sign | CODE: yes | CODE: yes | CODE: yes |
| S2− | negative coefficient recovered | **NOT TESTED** | **PASS** (synthetic ×0.25–×2, real A3 both halves) | **NOT TESTED** |
| S2+ | positive coefficient recovered | **PARTIAL** (synthetic v0.3 config only; real: B > 0 in 39/39 n³ winners) | **NOT TESTED** (only positive outcomes are reflection false alarms) | **PARTIAL** (synthetic only; real 'not found') |
| S3 | sign identified when unknown | **NOT TESTED** | **NOT TESTED** | **NOT TESTED** |
| — | decision rule uses sign | ignores sign | hard-coded a < 0 | hard-coded a > 0 |
| **2** | **Shape** | | | |
| F1 | k⁴-equivalent in the family | CODE: yes, p = 3 of {0,2,3,4} | CODE: yes, q = 2 of {1,2,3} | CODE: yes, q = 2 of {1,2,3} |
| F2 | k⁴-equivalent identified in injections | **PARTIAL** (100 % at all amplitudes, but v0.3 config, + sign, no systematics) | **PARTIAL** (0.93–1.00 at ×0.25–×1; 0.00 at ×2; 0.17 / 0.00 with R = 0.03 / 0.10; − sign only) | **PARTIAL** (×1: 0.6; ×2: 1.0; + sign only) |
| F3 | winner must beat the competing forms (discrimination statistic) | **NOT TESTED** (no such statistic exists) | CODE: yes (partial F-tests, l. 306–316) | CODE: yes (partial F-tests, l. 258–276) |
| F4 | shape on real data | **PARTIAL** (n³ wins 21/21 and 18/19, but n² and n⁴ also significant in every note) | **FAIL** (undetermined in both halves; search winner q = 2, confirmation winner q = 3) | **FAIL** (non-physical baseline, winner q = 1) |
| F5 | residual outside the family or a mixture | **NOT TESTED** | **FAIL** (full tanh law at ×2 read as k³ in 30/30) | **NOT TESTED** |
| **3** | **Recovery** | | | |
| A1 | amplitude estimate produced | CODE: yes (B_est, l. 362) | CODE: yes (a, l. 297–300) | CODE: yes |
| A2 | bias measured | **PARTIAL** (median B_est/B 0.998 → 0.884 for B = 1e-5 → 1e-3; v0.3 config; + only; no spread) | **PARTIAL** (median a/a_exp 0.93–0.98; m = 30; − only; no spread) | **PARTIAL** (0.93–0.98; m = 10; + only) |
| A3 | uncertainty interval and coverage | **NOT TESTED** | **NOT TESTED** | **NOT TESTED** |
| A4 | amplitude on real data vs an independent reference | **NOT TESTED** (no reference value in the repo) | **FAIL** (ratio 1.47 vs linear theory, outside calibrated 0.90–0.99) | — (not found) |
| **4** | **Injections:** executed vs only possible by arguments | see §6 | see §6 | see §6 |
| **5** | **Null control** | **PARTIAL** (independent real null passed; data identity not verifiable; v0.4.1 synthetic null not committed) | **PARTIAL** (synthetic only, 0/200) | **PARTIAL** (synthetic only, 0/200) |
| **6** | **Blind confirmation** | **PARTIAL** | **PARTIAL** | **PARTIAL** |
| **7** | **Systematics** | listed, **NOT TESTED** by injection | listed, **partly tested** | listed, **NOT TESTED**; dominated real data |

---

## §3 Sign (question 1): evidence

**S1, CODE: all three estimators are unconstrained in sign.**
- P: `fit_templates` (l. 295–329) solves an unconstrained least squares for a. The significance flag uses |cents10| (l. 353–354).
- G: `_wls` (l. 273–279), unconstrained.
- W: `_fit` (l. 218–223), `least_squares` without bounds.

This shows only that the code *can* return a negative value.

**S2: recovered sign, per pipeline.**
- **P:**
  - All committed injections are positive: `results/synth_inject/injections.json` has cents ∈ {0.5, 1, 2, 5, 20} and B ∈ {1e-5, 1e-4, 4e-4, 1e-3}, and **no negative case**. The `B_ratio_median` values are positive (0.884–0.998), so the positive sign was recovered.
  - Real data (`results/real_v041/notes_piano.csv`): B_est > 0 for all 39 notes whose winner is p = 3. The one other note (A2, winner p = 4) has B_est = −1.04e-4. That is a coefficient of a template that did not win; it is not a recovery result.
- **G:**
  - All scale injections are generated from the negative linear-theory law. All detected outcomes carry "−" (`v010_inject/synth_inject.json`, counts).
  - Real A3: a < 0 for all three templates in both halves (`v010_real_A3/summary.json`).
  - The only "+" outcomes anywhere are false alarms from the null with reflection R = 0.10: 18/30 flagged, all "+".
- **W:**
  - Injections are positive only: capillary C ≥ 0 (l. 421), and wrong shapes use a_cap > 0 (l. 423).
  - The declared outcome string has **no sign** (`_declared`, l. 392–396), so a recovered sign is not recorded in the injection output.

**S3: sign under unknown sign: NOT TESTED in all three.**
- No pipeline was ever run on injections of both signs, let alone with the sign hidden.
- The real-data decision rules are sign-specific or sign-blind by construction:
  - P ignores sign: `group_summary` l. 428–450 counts winners regardless of sign, and `found`/`confirmed` (l. 691–697) do not look at sign.
  - G requires a < 0 (`ok`, l. 557–558).
  - W requires a > 0 (`b_ok`, l. 506–507).
- The one observation relevant to S3 is negative evidence. In G, a known systematic (reflection) produces significant deviations **of the opposite sign** to the physical effect.

---

## §4 Shape (question 2): evidence

**F1: the family is closed and small.**
- P: {0, 2, 3, 4}, with p = 1 absorbed by f0.
- G and W: {1, 2, 3}.
- The scanners choose among these powers. They do not reconstruct an arbitrary residual. All three have an "undetermined" outcome:
  - P: share < 0.6 among significant notes (l. 442–449);
  - G and W: winner not distinguishable from every competitor (G l. 316; W l. 258–277).
- None has an outcome meaning "significant, but not any template".

**F2: injections.**
- **P:** `correct_shape` = 1.00 for every p ∈ {0, 2, 3, 4} and every amplitude, and for the exact stiff-string law. But:
  - these results carry config hash **639993ecfbe5fc3b (v0.3.0)**;
  - the real-data analysis used **5bdf34b190a4b9ef (v0.4.1)**;
  - the two configs differ in the frequency estimator (`freq_method` phase_demod), the onset rule and the minimum effect (`min_effect_cents` 0.1);
  - **no injection result for v0.4.1 is committed**;
  - the output does not record seed or noise arguments (`injections.json` has no `args`).
- **G:** q = 2 identified 0.93 / 1.00 / 1.00 at ×0.25 / ×0.5 / ×1. Weaknesses:
  - ×2: 0.00 (read as k³);
  - dispersion + R = 0.03: 0.17;
  - dispersion + R = 0.10: 0.00;
  - wrong shapes (kh)¹ and (kh)³: correctly identified 30/30 each, negative sign only.
- **W (v0.5.0):** capillary ×1 gives 6/10 "k² distinguishable", ×2 gives 10/10, ×0.5 gives 3/10 (read from `counts`; the `form_identification_rate` field is NaN for ×0.25 and ×0.5 in the committed file). Wrong shape k¹: 5/30 distinguishable; k³: 30/30.

**F3: discrimination.**
- P chooses the winner by minimum RSS (l. 328) and has **no statistic** testing whether the winner is needed in the presence of the other templates.
- On the real piano, every template is significant in almost every note:
  - search: n² 21/21, n³ 21/21, n⁴ 21/21, n⁰ 15/21;
  - confirmation: n⁴ 19/19, n³ 18/19, n² 18/19.
- The confirmation rule (l. 695–697) requires only that the same winner holds and that n³'s binomial p is < 0.01. In the confirmation group, n⁴ has the smaller binomial p (1e-38 vs 1.9e-35).
- So the piano result is "n³ fits best", not "n³ is distinguished from n² and n⁴".
- G and W do have this statistic.

**F4: real data.**
- P: winner n³ in 21/21 (search) and 18/19 (confirmation), with the F3 caveat.
- G: the code verdict is 'found but not confirmed'. Search winner is q = 2 and confirmation winner is q = 3; both are undetermined. Partial-test padj: search, 0.028 against both q = 1 and q = 3 (above α = 0.01); confirmation, 0.0006 against q = 1 but 0.22 against q = 2.
- W: all three templates are significant, the winner is q = 1, and g_eff = −2.7 / −3.7 (non-physical).

**F5: outside the family.**
- Only G contains a relevant case: its scale injections use the full tanh law, which has (kh)⁴ and higher terms. At ×2 that mixture is misread as k³ in 30/30.
- No pipeline injects a non-power residual.

---

## §5 Recovery (question 3): evidence

- **A2, bias.**
  - P (v0.3 config, + sign): median B_est/B = 0.998, 0.986, 0.947, 0.884 at B = 1e-5, 1e-4, 4e-4, 1e-3. The underestimate grows with B (higher-order terms). Only the median is stored, not the spread.
  - G (− sign): median a/a_exp = 0.929, 0.966, 0.982, 0.973 at ×0.25, 0.5, 1, 2; 0.987 and 0.902 with R = 0.03 and 0.10.
  - W (+ sign): 0.931–0.984.
  - Each value is a median over m = 10–60 realizations, for one sign per pipeline.
- **A3, uncertainty and coverage: NOT TESTED.** No Module 0 script computes a standard error or interval for a or B. The fit covariances are not stored (P `fit_templates`, G `_wls`, W `_fit`). So coverage cannot be assessed from existing outputs and would need code changes.
- **A4, real amplitude against a reference.**
  - P: B_est 1.1e-4 … 4.9e-4 (median 1.4e-4). `scanner/README.md` calls this "consistent with published B values", but **no published value or citation is in the repository**.
  - G: a(q = 2) = −0.202 / −0.205 vs linear theory −0.138 / −0.139, a ratio of 1.47. That is outside every calibrated systematic, and the authors themselves do not interpret it as a measured coefficient.
  - W: no amplitude (not found).

---

## §6 Injections (question 4): executed vs only possible

| pipeline | executed and committed | config hash of those results | possible only via arguments (never run) | not possible without code change |
|---|---|---|---|---|
| P | synth-null m = 300 (3 flags, FPR 0.010); synth-inject: power p ∈ {0,2,3,4} × cents {0.5,1,2,5,20}, m = 60; stiff B ∈ {1e-5,1e-4,4e-4,1e-3} | **639993ecfbe5fc3b (v0.3.0)**; none for v0.4.1 | negative `--cents`, negative `--B` (argparse float, l. 734–735); `--jitter`, `--noise-db` | recording the recovered sign per realization (the table stores only `detect`, `correct_shape`, median ratio); the README's v0.4.1 "synthetic nulls with ±5–30 cent vibrato: FPR 0" has **no committed output** |
| G | synth-null m = 200 seed 11001 (0/200); synth-inject m = 30 seed 12001: scale ×{0,0.25,0.5,1,2}, wrong (kh)¹,(kh)³, reflection 0.03/0.10, bound ×1/×3, and combinations | 4b2d81580df7981b (= code at prereg and at HEAD) | negative `--mults` (positive physical sign) | correct scoring of a positive-sign injection: the truth label is hard-coded `"k^2-"` (l. 491). The `counts` with sign are saved, so re-scoring from the JSON is possible. Wrong shapes and systematics with the opposite sign: `a_ref` is fixed negative (l. 488, 494) |
| W | synth-null B m = 200 seed 5001 (0/200); capillary ×{0…2} m = 10 seed 6001; wrong k¹,k³ m = 30 seed 7001 | 069ba514c71704e5 (v0.5.0); real run under c66f827f9b3b8c56 (v0.5.1, a translation, stated to give identical output; not verified here) | negative `--mults` (negative C) | any sign record: the outcome label has no sign (l. 396) |

No script has a blind mode where the injected truth is hidden from the analysis step and revealed afterwards.

---

## §7 Null control (question 5)

| | definition of a false alarm | barrier | synthetic null | real or independent null |
|---|---|---|---|---|
| P | note with winner template significant: Holm over 4 templates, α = 0.01, and \|deviation at n = 10\| ≥ 0.1 cent (l. 353–358; `barrier` uses `significant` only, l. 418–425) | one-sided binomial test of the count vs α = 0.01; **stop only if** p < 0.05 | v0.3 only: 300 notes, 3 flags | violin (development) v0.4.1: 49 notes, 1 flag, p = 0.389; **cello (final, untouched per README)**: 89 notes, 2 flags, p = 0.224 (`real_v041/summary.json`) |
| G | `any_flag` over the 3 templates (stricter) | same binomial rule | 0/200, non-dispersive linear waves | none (no non-dispersive real record) |
| W | `any_flag` | same | 0/200 (v0.5.0) | none; the real record is declared "not a null test" |

**Notes:**
- The barrier stops only when the false-alarm count is *significantly above* α. Passing it does not show FPR ≤ α; with M = 49–200 its power against a modest excess is limited.
- P's independent control (cello) is the only real null in Module 0. Its audio is not in the repository and has no hash in `results/data_manifest.csv` (86 piano + 90 violin entries only). So its identity cannot be verified from the repository.
- The v0.4.1 changes were made after the v0.3 violin null failed (FPR 0.368, `results/real/summary.json`). That is development on null data, which the README reports openly.

---

## §8 Blind confirmation (question 6)

| | split | criteria fixed in | order (git) | separate opening of confirmation | data identity |
|---|---|---|---|---|---|
| P | MIDI parity: even = search, odd = confirmation (CONFIG l. 87; l. 681) | CONFIG and code; hash 5bdf34b190a4b9ef, re-derived statically at tag `v0.4.1-piano-prereg` and at e7c7520 | prereg 44b56ac 18:52:41 → violin null 6010bb2 18:56:16 → result e32051a 19:00:00 (same day). Code identical prereg → result; to HEAD only comments translated | **no**: both groups are computed in one `run_real` call | piano audio committed in the prereg commit 44b56ac; manifest from 7a78143 (18:33, before prereg); the v0.3 run stopped at the violin barrier before the piano step (`real/summary.json` "stopped": true) |
| G | first / second half in time (l. 540) | code at prereg df689a8 (seeds included) + **an interpretation rule added in the calibration commit message** dee3de1 (17:44), before the real run 0e01b6e (17:51) | code unchanged prereg → HEAD | **no**: one call | data not committed, no hash; prereg message says only A1 coherence magnitudes were inspected beforehand |
| W | first / second half | prereg 5cb3947 (v0.5.0) | code translated before the real run (v0.5.1) | **no** | data not committed |

**Consequence for G:** two verdicts coexist.
- The code says 'found but not confirmed'.
- The commit-message rule says 'negative sign detected and confirmed, shape undetermined'.

Both are reported in the README. Only the code rule was frozen at preregistration.

All three: **PARTIAL.** Criteria were fixed before the real data were analysed, as far as git can show. But the confirmation half was never opened separately, and two of three data sets have no recorded hash.

---

## §9 Systematics that can mimic the residual (question 7)

| | identified in the repo | tested by injection | effect |
|---|---|---|---|
| P | unison beating (upper register excluded), onset transients and vibrato (null instruments), higher-order stiff-string terms, detuning of the nominal f0 | only the higher-order terms (exact √ law) and nominal detuning (`detune`, always on) | B underestimated by 12 % at B = 1e-3. **Any mechanism with Δf ∝ n³ would be indistinguishable**; none was injected |
| G | reflection, bound harmonics, depth and set-down (constant per block, absorbed by s_b), nonlinear shoaling (diagnostic only, not registered), higher-order dispersion | reflection 0.03 / 0.10, bound ×1 / ×3, higher-order law (×2) | R = 0.10 null: 60 % false alarms, **all positive sign**. R = 0.03 keeps amplitude but kills shape distinguishability. Higher order turns k² into k³. Real 1.47× unexplained |
| W | currents and drift (≈ 6 % shift, 40–80× the capillary term), stereo calibration (g_eff), leakage amplified by prewhitening | random U in synthetic fields only | dominated the real record: baseline non-physical, 'not found' |

This matches the repository's own conclusion (`scanner/README.md`, "Lessons"). Known physics missing from the baseline appears in the residual with a size comparable to or larger than the target, and can share its shape.

---

## §10 Answer to the main question

| part | status | basis |
|---|---|---|
| detect a residual without knowing it in advance | **PARTIAL** | P and G real data: reproducible, known-sign deviation; W: FAIL |
| establish k⁴ dependence | **PARTIAL** | P only, by best fit, without a discrimination statistic; G real: undetermined; G synthetic: fails for a mixed-order law |
| recover the coefficient with free sign | **NOT TESTED** | S3 never run; each pipeline saw one sign; decision rules are sign-specific (G, W) or sign-blind (P) |
| recover the amplitude | **PARTIAL** (bias) / **NOT TESTED** (coverage) / **FAIL** (G real vs theory) | medians only; no intervals anywhere |

**Conclusion.** At e7c7520, Module 0 does **not** demonstrate that the method can extract an unknown k⁴ coefficient with free sign and amplitude. It demonstrates three narrower things:
1. detection of a known-sign dispersive deviation on real data in two systems;
2. selection of the k⁴-equivalent power within a fixed four-member family in one system (piano), by best fit;
3. synthetic amplitude bias within about 12 %, for one sign per pipeline.

**Missing pieces, each traceable to the matrix:**
- (i) sign identification under an unknown sign (S3);
- (ii) a discrimination statistic in the piano pipeline (F3);
- (iii) interval estimation and coverage for the amplitude (A3);
- (iv) validation of the configuration actually applied to real piano data (F2, A2: v0.4.1 injections not committed);
- (v) behaviour on mixtures and on residuals outside the family (F5);
- (vi) verifiable identity of the independent null and real data sets (§7, §8).

---

## §11 Facts for the later review of B (not a specification)

From the code alone, without running anything:

1. **P** can receive symmetric injections through arguments (negative `--cents`, negative `--B`; B below −1/n² removes partials, because NaN frequencies are dropped by the `keep` filter, l. 389). But:
   - the scoring does not record the recovered sign per realization;
   - a sign test therefore needs either a code change or a new wrapper that calls `analyze_signal`, which is new code;
   - the injections must be redone under v0.4.1 in any case (F2).
2. **G** can receive the opposite sign (negative `--mults`), but scores it against a hard-coded `"k^2-"`. Re-scoring from the saved `counts` (which carry the sign) is possible without a code change. Wrong-shape and systematic cases cannot be flipped without one.
3. **W** cannot report sign at all without a code change.
4. **No pipeline** produces an amplitude interval, so a coverage test (A3) is impossible without code changes in all three.
5. **Tuning risk:** P and G real data have already been analysed. A B test that uses only synthetic injections does not touch them. But its pass criteria must be justified from the synthetic calibration, not from the known real outcomes (n³ in 39/40, 1.47×).

Whether a symmetric test is admissible as "no change to the registered protocols" is therefore a decision for the B review. For P and W, it cannot be done purely by arguments.

---

## §12 Limits of this audit

- Read from a read-only clone. Statements about what was run rely on committed outputs and commit messages. An uncommitted run cannot be seen.
- Line numbers refer to e7c7520.
- The config hash was re-derived only for piano (prereg tag and HEAD). G is covered by "no diff between prereg tag and HEAD". The W v0.5.0 → v0.5.1 equivalence was not checked.
- Commit timestamps are author-set and are not proof of order. Tags and parent order are.

---

## Appendix: SHA-256 of audited files at e7c7520

```
65bcc9bece26d8ee28e4f2c76a9fbeb1fb151d6908c3af1b33bbcf5290c16750  scanner/README.md
ffb785e938e4af7913b1e43c39790b9b4b925a1a387e0d2756000c02f86a668c  scanner/piano/README.md
335ba8b8abaef4b424b36a516f55df568b9345f83cc8b6bc1c0dc1cc33c9e5d4  scanner/piano/piano_scanner.py
31be210f5c1672e42d8b06a3d2dc8bb727c274ce954dbaaf7a1e34fae1296f49  scanner/globex/globex_scanner.py
2cf79ac28a21658c3f1ef9bf40c252fa4dfbd5547639287893182d93794af44d  scanner/waves/waves_scanner.py
20dfba3795c2791f169a9f16eec55db4600a995bf5bbb54985d3664e911e2a8b  scanner/piano/results/synth_inject/injections.json
12879dc2eb6cd746655c7bb2962f2834f8fc950ab9bbdbccebcc7fad20ab5aba  scanner/piano/results/synth_null/summary.json
1a67ab0beb6e564b18419d2ba85843624934ab7cb90a0d397f038c1665182de2  scanner/piano/results/null_violin_v041/summary.json
f5a68ac3f2fd66c3bca00a7289b87ee2f9c6bfe749a4de22419fcc6669d5152d  scanner/piano/results/real_v041/summary.json
84e2f0dbad975ce07627a9f997a0c3b60afe787fb73407cce4c4ff876b6ac793  scanner/piano/results/real_v041/notes_piano.csv
7351777557a72717027d6c9aa79931d8d44d810fc811ccb39ec5f3c2e75472d2  scanner/piano/results/real/summary.json
6af5f0cc15b371f1fdaf17e8cf773f5c08ea7d70c9e03e60a1455ef7cb76237c  scanner/globex/results/v010_null/synth_null.json
eb64525c956397758f3c55583820674e5161f8620aee15b478a96583b9c5ab48  scanner/globex/results/v010_inject/synth_inject.json
83d9c89c9195f8926e1199115a13138b82fcae2724d93d2dbb22f3d741f5d93e  scanner/globex/results/v010_real_A3/summary.json
6b1b02e27432bf121f1af5938c2477ab3f00cea63d83028477d7d4e78f120665  scanner/waves/results/v050_null_B200/synth_null.json
ebd40ee00269bd9b16d72978da63bee7c7181746c0815f5494b66937259479a8  scanner/waves/results/v050_inject_cap/synth_inject.json
f265fde3b57d0a479f8f8377ac4f7cfdd04f0e1d79ca8f49496c31cd9dc50aa4  scanner/waves/results/v050_inject_wrong/synth_inject.json
4a1e8bea79a066ae4ce8ac68a3446595456696dc40b51f0a8b4ce575df4d3be4  scanner/waves/results/v051_real/summary.json
```

Tags used: `v0.3-piano-prereg` 078e546, `v0.3.1-piano-data` 7a78143, `v0.4.1-piano-prereg` 44b56ac, `v0.4.1-piano-result` e32051a, `v0.5.0-waves-prereg` 5cb3947, `v0.5.1-waves-english` 4a4a756, `v0.1.0-globex-prereg` df689a8, `v0-module0-closed` 6d77eaf.
