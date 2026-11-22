# B2 results: blind k⁴-equivalent recovery test, piano pipeline v0.4.1 + B1

**Status:** draft for review. Not committed. It reports the single B2 run defined by B0 v0.4 and does not change any criterion.
**Date:** 2026-11-22.
**Verdict:** **V4 — DETECTION, SIGN OF B AND SHAPE; NOT QUANTITATIVE**, set by the frozen ladder of B0 §7.

---

## 1. Frozen provenance

- **Repository:** github.com/Kretski/residual-first-framework.
- **Frozen commit:** `d40144d7f3ed61421ffc5d2b7a1e1e807c1c4ca7`, tag `k4-recovery-b0b1-frozen`, on `main`. Its parent is `e7c7520`. It adds only `scanner/piano/k4_recovery/` (7 files).

| frozen file | SHA-256 |
|---|---|
| `RFF_K4_RECOVERY_AUDIT_v0.1.md` (audit A) | `eb7cb961569d5becda3a277ffaa1b1c4b607f82d48315952084238307952fc97` |
| `RFF_K4_RECOVERY_TEST_SPEC_B0_v0.4_CORRECTED.md` (B0) | `662b2303aff00bb1000d8a0ed93be64fc5f7fdc72aae811567fb2f5810a4d042` |
| `k4_recovery_b1.py` (wrapper v0.3) | `75b8bbf6517887300be52fa403fcaa14523377a7c5e4402ed60deab14e72eb4a` |
| `k4_score_b1.py` (scorer) | `8415d1df7a89863b9c62386d3c2f2bf11d774a6358ada88a506e33f35f24e8dc` |
| `B1_IMPLEMENTATION_NOTES_v0.5.md` | `2f58b30e46f6120e26d50ad007729b9961f5999bb224619001690e63f6fae9b9` |
| `COPY_VERIFICATION_e7c7520.txt` | `ab32613f1fcd4fa848ed68eabcfdb5ae23c2a376c7808b30de3bec930e2ac252` |
| `K4_RECOVERY_FREEZE_RECORD.md` | `4fcabceedfa32bf2d2902fdf62dbff70279f1b8e6e47738527faa38ee750c7cb` |
| `scanner/piano/piano_scanner.py` (unchanged since e7c7520) | `335ba8b8abaef4b424b36a516f55df568b9345f83cc8b6bc1c0dc1cc33c9e5d4` (config hash `5bdf34b190a4b9ef`) |

**Chronology of the freeze (stated as agreed in review):**
- Commit `d40144d` was created on 2026-11-22 at 04:15 (UTC+02:00), **before the final written approval for the push**.
- Its content was then checked file by file against the approved set by SHA-256, both on GitHub and in the B2 clone, and found identical.
- The frozen commit has not been changed since.

---

## 2. Execution integrity

- **Script:** `run_k4_b2_lenovo.ps1` v2 (SHA-256 `84d23cc02a9f57ebc87f006d3bd8dccb02e9d98d27bcf6f90214309eddbf8555`). It runs outside the repository and was approved after review.
- **Machine and folder:** Lenovo, `D:\research\rff_k4_b2_20261122_042506`.
- **Clone:** fresh clone, detached at `d40144d`; the tag was verified and all eight frozen files were checked by SHA-256 in place.
- **Start:** `B2_STARTED.txt` was written at **2026-11-22T04:25:16+02:00**, immediately before the wrapper started. There is no other B2 folder.
- **Environment:** Python 3.13.12, numpy 2.5.2, scipy 1.18.1, recorded with `python -I`. This is the same environment as smoke run 2.
- **Run:** seed 20261122, **5900/5900 cases, 0 excluded**, wall time 4679 s (0.79 s/case).
  - The parity guard against `analyze_signal` never triggered, which it would have done on any mismatch.
  - `truth.jsonl` and `analysis.jsonl` have 5900 rows each.
- **Scorer:** verified the four run hashes before reading anything, then wrote `score.json`.

| output | SHA-256 |
|---|---|
| `b2_out/manifest.json` | `c4d5f4743c6bab479a07243cae213e6613a75277a45bbb72828ac4c940a7c295` |
| `b2_out/truth.jsonl` | `4c690be00df2588fd08610edcca809de2cc87551cda19df078371107c53e8412` |
| `b2_out/analysis.jsonl` | `eda9e9857ee214d0aab7f02a4996053ae1b95fffd64ed3103e9e4c2aa66c306b` |
| `b2_out/run.log` | `6af34c02d5fcac2a2c18f5a27eec8b5e8c602f66b4759e6adfdc62d1383a6ca1` |
| `b2_out/hashes.json` | `de5572725215227c3995a4270fcb8a3206a6d63ace63d0a8109f4251c3a0bb0d` |
| `b2_out/score.json` | `3fa92001293d472319080660eea2a80ff893382e518518bea80ca8bded9341c0` |
| `B2_STARTED.txt` | `9e128c5e376084f77c6010e4fbf74895debacbcfaa3442fe039e52a7554cec62` |
| `b2_environment.txt` (UTF-16, as written on Lenovo) | `e86295ad1bdf4b12bfa14c2159dfac620bb56f899e4b57af6cc6191eba213921` |
| `b2_wrapper_stdout.txt` (UTF-16) | `7ddcebcb338e8ba3fbef787e73e27fd30720e1c0083c124e365ea66be4b5e4a8` |
| `b2_scorer_stdout.txt` (UTF-16) | `a1e9d0ef3ad687a19adabd767b8be9991b63b43a5b81b0e634809202ed90da36` |

The uploaded copies of the data files (all except the three console files) were re-hashed on receipt and match exactly.

---

## 3. Results, criterion by criterion (from `score.json`; CP = Clopper–Pearson 95 %)

All 36 cells were evaluable: 100 % of cases were analysed in every cell.

| criterion | measured | rule (B0 §5) | result |
|---|---|---|---|
| **C1** null FPR | 0/2000 (CP upper 0.0018) | CP upper ≤ 0.02 | **PASS** |
| **C2b** sign of B, + | detection 500/500, sign 500/500 (CP lower 0.9926 each) | detection ≥ 0.95, sign ≥ 0.98 | **PASS** |
| **C2b** sign of B, − | detection 500/500, sign 500/500 (CP lower 0.9926 each) | same | **PASS** |
| C2a winning-deviation sign (diagnostic) | +: 900/900 and 900/900; −: 900/900 and 900/900 | — | would pass its bars |
| **C3a** pure shapes | n²: 400/400 correct, 0 wrong; n³: 400/400, 0; n⁴: 400/400, 0 (CP lower 0.9908, CP upper of wrong 0.0092) | correct ≥ 0.80; wrong ≤ 0.02 | **PASS** |
| **C3c** out of family | n⁵ claimed as n³: 0/200 (CP upper 0.018); declared `n^4± distinguishable` 200/200 | CP upper ≤ 0.05 | **PASS** |
| **C3b** stiff law | n³ distinguishable 600/600 (every cell 100/100); wrong 0/600 | n³ ≥ 0.50; wrong ≤ 0.02 | **PASS** |
| **C4** standard: the 6 decision cells P n³ 2 c ±, P n³ 5 c ±, S 1e-4 ± (see §4, group D) | the four pure n³ cells meet the per-cell limits; S ±1e-4 coverage 0/100 each; pooled 378/600 = 0.63 | each cell \|bias\| ≤ 0.10 and coverage ≥ 0.85; pooled ∈ [0.92, 0.99] | **FAIL** |
| C4-PURE (V3 branch only): the 4 cells P n³ 2 c ±, 5 c ± (group D, marked ᴾ) | pooled 378/400 = 0.945 | same limits | would pass; **not used**, because the ladder did not enter V3 |
| C4-range: the scored stiff levels S 1e-4, 4e-4, 1e-3 ± (groups D and R) | no level passes for both signs | — | none |
| **C5-sign** contamination | sign of B̂ correct 100/100 for each sign (CP lower 0.964) | CP lower ≥ 0.90 | **PASS** |

**Ladder:** C1 → C2b → C3a/C3c → C3b all pass; **the first failing level is C4**, which gives **V4**.

**Qualifier (B0 §7):** n³ was recognised in ≥ 50 % of pooled stiff-law cases (in fact 600/600). The per-|B| rates are reported, and no per-value claim is made from them.

---

## 4. C4 failure analysis (descriptive; no criterion changed)

**Per-cell amplitude results (bias = median(B̂/B) − 1; coverage of the 95 % W3 interval).** The table shows more cells than the C4 decision uses. The **group** column says which role each cell has:

- **D: decision cells** of standard C4 (B0 §5, C4-core): P n³ 2 c ±, P n³ 5 c ±, S 1e-4 ±. ᴾ marks the four of them that also form C4-PURE (used only in the V3 branch, which was not reached).
- **R: range cells**, which enter only the C4-range statement: S 4e-4 ±, S 1e-3 ±. S 1e-4 ± are decision cells **and** the first range level.
- **S: sensitivity cells** (|d₁₀| < 2 cents): **reported only**, not part of any C4 rule.

| group | cell | B | bias (+ / −) | coverage (+ / −) |
|---|---|---|---|---|
| S | P n³ 0.5 c ± | ±5.8e-6 | 9e-7 / 4e-6 | 94/100, 98/100 |
| D ᴾ | P n³ 2 c ± | ±2.3e-5 | 8e-7 / 2e-7 | 93/100, 96/100 |
| D ᴾ | P n³ 5 c ± | ±5.8e-5 | 5e-7 / −1.5e-7 | 98/100, 91/100 |
| S | S 1e-5 ± | ±1e-5 | −0.0014 / +0.0015 | 0/100, 0/100 |
| D, R | S 1e-4 ± | ±1e-4 | −0.0140 / +0.0147 | 0/100, 0/100 |
| R | S 4e-4 ± | ±4e-4 | −0.0525 / +0.0577 | 0/100, 0/100 |
| R | S 1e-3 ± | ±1e-3 | −0.1182 / +0.0862 | 0/100, 0/100 |

**The C4 decision (group D only):**
- the four ᴾ cells meet the per-cell limits;
- S 1e-4 ± fail on coverage (0/100 < 0.85; their bias of 1.4 % is within 0.10);
- the pooled coverage over the six D cells is 378/600 = 0.63, outside [0.92, 0.99].

So C4 fails. The S and R rows do not affect the verdict.

**Observations:**
1. **Pure power law: B̂ is unbiased and its interval is calibrated.** Bias is ≤ 4e-6 and coverage is 0.91–0.98 per cell, 0.945 pooled. This is where the p = 3 model is correctly specified (B0 W3).
2. **Exact stiff-string law: the bias is systematic and grows with |B|.** It is about 14 × |B| in relative terms up to 4e-4. Its sign does not follow the sign of B: B̂ is low for B > 0, and |B̂| is high for B < 0. This pattern is consistent with the term −B²n⁵/8 of the expansion f_n = n f0 (1 + Bn²/2 − B²n⁴/8 + …), which the p = 3 model omits. **This is an interpretation consistent with the numbers. It was not tested separately.**
3. **The interval does not cover the truth at any |B|.** With noise at −60 dB the statistical se is far smaller than even a 0.14 % model bias, so coverage is 0/100 in every stiff-law cell. At |B| = 1e-4 the bias (1.4 %) is well within the 10 % bias limit; the failure is entirely one of coverage.
4. **The scanner registers the omitted term.** In every stiff-law cell `mixture_like` = 0.92–1.00, meaning the reverse tests show a significant second contribution. In the pure-power cells it is ≤ 0.02, and in the null 0.001.
5. **This mis-specification was stated in advance.** B0 W3 says that the stiff-law coverage "is never presented as calibrated by construction". The C4 failure is therefore the pre-declared risk, now measured. It is not a code defect, and B0 §7 allows no revision for it.

---

## 5. Reported diagnostics (do not change the verdict)

- **R stratum (jitter 0.3 cent, noise −40 dB):**
  - null: 4/500 declared (0.8 %, CP upper 0.020);
  - pure n³ 5 c: detection and sign 100/100 for both signs, n³ distinguishable 100/100 and 96/100, coverage 0.92 / 0.89;
  - stiff 4e-4: detection and sign 100/100, n³ 100/100, coverage 0.01 / 0.02, bias −5.2 % / +6.1 %.
- **Block C (contamination, opposite-sign n²):** the sign is preserved in 100/100 cases for each sign. The bias in B̂ is −17 % for both signs, as expected for a mis-specified model; `mixture_like` is 0.91–0.94.
- **Block O (n⁵):** every case was declared `n^4± distinguishable`, with the correct sign. This is the closed-family limitation of B0 §1. The method does not claim k⁴ here, but it does assign the nearest family member with confidence.

---

## 6. Process notes

1. **Scorer exit code 1 (console only).**
   - The scorer wrote `score.json` completely and then failed while **printing** the first qualifier: the Windows console encoding cp1251 cannot represent "³" (`UnicodeEncodeError`).
   - The process exit code was therefore **1**. It is reported as such and **not as a successful exit**.
   - The verdict line had already been printed. No computation was affected, and `score.json` is the authoritative output.
   - A re-run of the scorer on the same files with `PYTHONIOENCODING=utf-8` is an optional, separate re-scoring of the same data, not a new B2. It has not been done.
2. **Encoding of the console files.** `Tee-Object` under Windows PowerShell 5.1 writes UTF-16. The uploaded copies of the three console files were converted to UTF-8, so their hashes differ. For `b2_environment.txt`, re-encoding the uploaded text to UTF-16 reproduces the Lenovo hash exactly.
3. **No second B2.** B2 is complete. No B0, B1 or threshold change is permitted under B0 §6.5 and §7.

---

## 7. Final verdict

**V4 — DETECTION, SIGN OF B AND SHAPE; NOT QUANTITATIVE**, for the piano pipeline v0.4.1 + B1 wrapper v0.3, synthetic, within the family {0, 2, 3, 4}.

This verdict is final for this pipeline under the B0 stopping rule. It enters the Λ map as the status of "independent extraction method, Module 0, piano pipeline".

---

## 8. Scope and limitations

- **Synthetic only.** These results describe the behaviour of the method on the `synth_note` generator distribution. Real recordings, real systematics and real dispersion laws were not tested.
- **Shown under these conditions:**
  - no false alarms;
  - blind recovery of the sign of B with no error in 1000 target cases;
  - recognition of n³ against n² and n⁴;
  - no k⁴ claim for an n⁵ residual;
  - a calibrated interval for B **only when the injected law is exactly the p = 3 template**.
- **Not shown:**
  - a calibrated interval for B under the physical stiff-string law, or under any law with additional contributions. That is the C4 failure;
  - recognition of residual shapes outside the family;
  - anything about Λ in any physical sector.
- **Consequence:** under this pipeline, B may be reported as **sign and shape only, not as a value with an interval**, whenever the true law may contain more than the p = 3 term. That includes the earlier real-piano B values.
- **Further development is outside this test and was not started.** For example, a model term for higher orders, or an interval that accounts for systematic model error. It would need its own specification and its own synthetic validation. It is not a retrospective repair of B2.

---

## 9. Results package (`scanner/piano/k4_recovery/results_b2/`)

All files are byte-identical to the files written on Lenovo. Their SHA-256 match the hashes recorded by the run (`hashes.json`, `b2_output_hashes.txt`).

| file | SHA-256 | note |
|---|---|---|
| `manifest.json` | `c4d5f4743c6bab479a07243cae213e6613a75277a45bbb72828ac4c940a7c295` |  |
| `hashes.json` | `de5572725215227c3995a4270fcb8a3206a6d63ace63d0a8109f4251c3a0bb0d` |  |
| `score.json` | `3fa92001293d472319080660eea2a80ff893382e518518bea80ca8bded9341c0` | authoritative scorer output |
| `run.log` | `6af34c02d5fcac2a2c18f5a27eec8b5e8c602f66b4759e6adfdc62d1383a6ca1` |  |
| `truth.jsonl` | `4c690be00df2588fd08610edcca809de2cc87551cda19df078371107c53e8412` | 5900 rows |
| `analysis.jsonl` | `eda9e9857ee214d0aab7f02a4996053ae1b95fffd64ed3103e9e4c2aa66c306b` | 5900 rows |
| `B2_STARTED.txt` | `9e128c5e376084f77c6010e4fbf74895debacbcfaa3442fe039e52a7554cec62` |  |
| `b2_environment.txt` | `e86295ad1bdf4b12bfa14c2159dfac620bb56f899e4b57af6cc6191eba213921` | UTF-16LE with BOM, as written by `Tee-Object` |
| `b2_wrapper_stdout.txt` | `7ddcebcb338e8ba3fbef787e73e27fd30720e1c0083c124e365ea66be4b5e4a8` | UTF-16LE with BOM |
| `b2_scorer_stdout.txt` | `a1e9d0ef3ad687a19adabd767b8be9991b63b43a5b81b0e634809202ed90da36` | UTF-16LE with BOM; contains the printing traceback |

**The three console files.**
- The uploaded copies had been converted to UTF-8 on the way.
- Re-encoding them as UTF-16LE with a byte-order mark reproduces the Lenovo SHA-256 of **all three** exactly. They are therefore committed in their original encoding.
- `b2_scorer_stdout.txt` shows the em dash of the verdict line as "Ч", which is how the cp1251 console rendered it. That is part of the original record and is not corrected.
