# B1 implementation notes v0.5 — wrapper and scorer for the blind k⁴ recovery test

**Status (as of writing, 2026-11-22):**
1. **Executed in smoke runs only.**
   - Wrapper v0.3 and the scorer ran in **smoke run 2** on Lenovo (seed 1, 108 cases, technical checks only, no rates computed).
   - Wrapper v0.2 ran in smoke run 1.
   - Details are in §8.
2. **Smoke run 2 passed:**
   - 108/108 cases analysed, 0 excluded;
   - parity guard clean;
   - W3 cross-check present in 108/108 cases, **0 failures**;
   - scorer exit 0.
3. **No commit** to residual-first-framework exists at the time this document is written. The freeze commit is proposed separately (`K4_RECOVERY_FREEZE_RECORD.md`).
4. **B2 (seed 20261122, 5900 cases) has not been run and is not authorised.**
5. A smoke pass is a technical check. It does not validate detection, sign recovery or the recovery of B.

**Changes from v0.4 (after smoke run 1, 2026-11-22):**
- **Wrapper v0.3, smoke path only:** `crosscheck_se` passes the analytic Jacobian [n, n³] to `curve_fit`.
- **Reason (smoke run 1, `D:\research\rff_k4_smoke_20261122_035543`):**
  - 21/108 cross-check failures;
  - B agreed to ≤ 7e-5, but the reference se differed by 1e-6 … 1.1e-1;
  - cause: the finite-difference Jacobian of the reference.
- The B2 analysis path, B0, every threshold and the tolerance 1e-6 are unchanged.
- **Smoke run 2** (`D:\research\rff_k4_smoke_20261122_040135`) passed (§8).

**Changes from v0.3 (second smoke-script review of 2026-11-22):**
- **Smoke script only:** the environment check now runs as `& $Py -I -c ...`, under the same isolated mode as the wrapper and scorer calls. It cannot pass under conditions different from the real smoke run.
- Wrapper and scorer are unchanged.

**Smoke script (outside the repo, run on Lenovo):** `run_k4_smoke_lenovo.ps1`, SHA-256 `50ea4315fc538be2f2fb901bf91a4ccca58ac880c88c5d8815fd4649b9581d48` (ASCII, CRLF; selects wrapper v0.3). Earlier versions: `fdfa8429…`, `46762d48…` and `2b3c7eb1…` (the last used for smoke run 1).

**Changes from v0.2 (smoke-script review of 2026-11-22):**
- **Scorer, smoke mode:** every analysed (non-excluded) case must carry the W3 cross-check. A missing cross-check, a count different from the number of analysed cases, or no analysed cases at all are structural problems, and the scorer then exits with code 1. The wrapper is unchanged.
- **The smoke script** (`run_k4_smoke_lenovo.ps1`, outside the repo):
  - the environment check runs in the same stderr-safe block as the wrapper and the scorer;
  - the scorer runs only if the wrapper exits with 0;
  - the script ends with exit code 2, 3 or 4 and a FAILED line if the environment, the wrapper or the scorer fails, and with 0 only when all structural checks pass.
- **No B0 criterion, threshold, seed or case count changed.**

**Changes from v0.1 (static review of 2026-11-22):**
- **A, wrapper:** `sign_B` and `sign_winner` are now three-valued, via `np.sign`. B̂ = 0 gives 0, which never counts as a correct sign for either injection. A declared outcome cannot carry sign 0, because a significant template has \|deviation at n = 10\| ≥ 0.1 cent; the label `0` exists only as a guard.
- **B, scorer:** `hashes.json` must list exactly `manifest.json`, `truth.jsonl`, `analysis.jsonl` and `run.log`. A missing or extra key aborts the scoring before anything is read.
- **C, scorer:** `verdict.required_qualifiers` collects every qualifier required by B0 §7:
  - C2a rates (diagnostic);
  - C5-sign;
  - C3b pooled qualifier and per-cell rates;
  - C3c declared shapes;
  - `mixture_like` per cell;
  - C4-range (none in V3);
  - non-evaluable cells;
  - the R stratum.

  The main ones are also printed. No criterion changed.
- **D, scorer:** the spawn indices must be exactly 0 … n_cases − 1 (complete, unique, in range).
- **Copy claims made checkable:** `COPY_VERIFICATION_e7c7520.txt` (§7) compares the copied lines with the e7c7520 sources.
- **Audit A** added to the folder unchanged (sha256 `eb7cb961…7fc97`).

Status text of v0.1 (historical; superseded by the status above): at v0.1 the only checks done were `python -m py_compile` (byte-compilation, no module code runs) and `pyflakes` (static lint, clean).
**Date:** 2026-11-22.
**Specification:** `RFF_K4_RECOVERY_TEST_SPEC_B0_v0.4_CORRECTED.md`, frozen, SHA-256 `662b2303aff00bb1000d8a0ed93be64fc5f7fdc72aae811567fb2f5810a4d042`. The hash was re-verified, and the diff against the draft v0.4 contains only the reviewer's C4-PURE and N = 100 additions. The C4-PURE operating characteristics were recomputed independently and match.
**Base:** residual-first-framework **e7c7520**. `scanner/piano/piano_scanner.py` is not modified.

---

## 1. Files (all new; nothing existing changes)

| path in residual-first-framework | lines | SHA-256 |
|---|---|---|
| `scanner/piano/k4_recovery/k4_recovery_b1.py` | 547 | `75b8bbf6517887300be52fa403fcaa14523377a7c5e4402ed60deab14e72eb4a` (v0.3; v0.2 was `b61e9b44…`) |
| `scanner/piano/k4_recovery/k4_score_b1.py` | 430 | `8415d1df7a89863b9c62386d3c2f2bf11d774a6358ada88a506e33f35f24e8dc` (v0.2 was `6c8305e7…`) |
| `scanner/piano/k4_recovery/COPY_VERIFICATION_e7c7520.txt` | — | `ab32613f1fcd4fa848ed68eabcfdb5ae23c2a376c7808b30de3bec930e2ac252` |
| `scanner/piano/k4_recovery/RFF_K4_RECOVERY_AUDIT_v0.1.md` | — | `eb7cb961569d5becda3a277ffaa1b1c4b607f82d48315952084238307952fc97` (audit A, unchanged) |
| `scanner/piano/k4_recovery/RFF_K4_RECOVERY_TEST_SPEC_B0_v0.4_CORRECTED.md` | — | `662b2303…a042` (frozen B0, copied unchanged) |
| `scanner/piano/k4_recovery/B1_IMPLEMENTATION_NOTES_v0.5.md` | — | this file (replaces v0.4) |
| `B1_wrapper_v0.2_to_v0.3.diff` (review aid) | — | `b6ddf1aacff200fc8e3209b57ad51795ae54fb2042dbf8bf47392124fc4f1a65` |
| `B1_v0.2_to_v0.3.diff` (review aid; not proposed for commit) | — | the scorer and smoke-script changes of v0.3 |
| `B1_v0.1_to_v0.2.diff` (review aid; not proposed for commit) | — | `cf25d8b193c8a6a80c290738353da5055e29256f0eab0e5b8509bcf4c8cc6b6d`. Exact: the v0.1 sides were reconstructed by reversing the edits and re-hash to `6277d9ef…` and `2813ac16…` |

**Traceability chain in the folder:**
- audit A: why B0 was needed;
- B0 v0.4 (frozen): the pre-fixed experiment;
- B1: the implementation;
- B2: the single run (later).

**Changes relative to e7c7520:** only the files above are added. The wrapper verifies at start-up:
- the scanner file hash `335ba8b8…e5d4`;
- `config_hash()` = `5bdf34b190a4b9ef`;
- templates `[0, 2, 3, 4]` and α = 0.01.

It aborts on any mismatch.

---

## 2. Mapping B0 → code

| B0 item | function (file) | notes |
|---|---|---|
| §2 frozen inputs | `verify_frozen_inputs` (wrapper) | abort on mismatch |
| W1 sign record | `analyse_case` → `sign_winner`, `sign_B` | three-valued (+1, 0, −1) |
| W2 discrimination | `discrimination`, `ftest_globex` | `ftest_globex` is a literal copy of `globex_scanner._ftest`; Holm is the scanner's own `ps.holm`; dof = N − 3 |
| W3 B̂ and interval | `b_hat_interval` | design built exactly as `fit_templates` (x = n/n_max, w = SNR/median, columns × w); delta method; t₀.₉₇₅,N−2 |
| W3 checks | `analyse_case` (parity of a₃, f0, B_est at 1e-12); `crosscheck_se` (smoke only, `curve_fit`, rtol 1e-6) | |
| W4 declared outcome | `analyse_case` | `none` / `n^p± distinguishable` / `n^p± undetermined` |
| W5 generators | `make_signal` (uses `ps.synth_note`, `ps.random_f0s`, `ps.detune` unchanged); `synth_mixture` (new) | `synth_mixture` copies `synth_note` l. 387–403 line by line after its own frequency law, keeping the same RNG order |
| W6 manifest and seeds | `cell_definitions`, `build_manifest`, `assert_manifest`, `run` | counts asserted before any signal exists; `SeedSequence(seed).spawn(n_cases + 1)`; the last child shuffles |
| parity guard | `analyse_case` | exclusion reason, winner, significant, per-template p and padj, W3 coefficients; `ParityError` aborts the run |
| §6.2 hidden truth | `run` | the analysis path gets `(x, sr, f_nom)` only; `truth.jsonl` and `analysis.jsonl` are separate; case ids follow the shuffled position; `hashes.json` is written after all files are closed |
| §6.4 smoke | mode `smoke` (seed 1, 3 per cell, 108 cases, cross-check on); scorer `--smoke` | the scorer prints counts, schema, hashes and cross-check failures; **no rates** |
| §5 C1–C5, §7 ladder | `score` (scorer) | §4 below lists the implementation choices |

---

## 3. Cell definitions (`cell_definitions`; asserted against B0 §4)

| block | cells | count |
|---|---|---|
| N | `N` (m = 2000) | 2000 |
| P | `P_p{2,3,4}_c{0.5,2,5}_{+,-}` | 1800 |
| S | `S_B{1e-05,0.0001,0.0004,0.001}_{+,-}` | 800 |
| O | `O_p5_c2_{+,-}` | 200 |
| C | `C_mix_{+,-}` | 200 |
| R | `R_N` (500), `R_P_p3_c5_{+,-}`, `R_S_B0.0004_{+,-}`; noise −40 dB, jitter 0.3 cent | 900 |

**Total:** 5900 cases in 36 cells. Truth B: for pure n³, B = 0.02(2^{c/1200} − 1), which does not depend on f0. For S it is the injected B. For C it is the B of the n³ part.

---

## 4. Implementation choices that the review should confirm

These are points where B0 fixes the rule but leaves a detail open. Each is implemented as stated and can be changed only before freezing B1.

1. **`mixture_like` (W2.6).** The forward family T(p*|q) and the reverse family T(q|p*) are each Holm-adjusted over the three competitors. A case is flagged if for some q both adjusted p < 0.01. The flag does not affect the verdict.
2. **Denominators.**
   - C1, C3a, C3b, C3c and C4 use the **analysed** cases of a cell, i.e. those not excluded by the scanner's rule.
   - C2a and C2b use n_analysed for detection and n_detected for the sign, as B0 states.
   - C5 uses the analysed cases of block C.
3. **Non-evaluable cells.**
   - They are dropped from pools (C2a, C2b, C3a, C3b, C3c).
   - If a criterion is left with no evaluable cell, it fails: CP lower = 0, CP upper = 1.
   - For C4, a non-evaluable **core** cell makes C4-standard or C4-PURE fail, and the cell is listed.
   - For C1, a non-evaluable `N` cell makes C1 fail.
4. **C4-range.** The largest scored |B| (1e-4, 4e-4, 1e-3) such that the per-cell conditions hold for both signs at that level **and every smaller scored level**, i.e. contiguous from the smallest. The per-level status is also reported, so a non-contiguous pattern is visible.
5. **C4 bias and coverage** are computed over all analysed cases of the cell, not only detected ones. B̂ exists for every analysed case (W3).
6. **The smoke cross-check** uses `curve_fit(method="lm", xtol = ftol = gtol = 1e-15)`, starting from median(f/n) rather than from the WLS solution. An exception is recorded as a failed check, not an abort.
7. **Negative B.** `synth_note` returns NaN for partials with 1 + Bn² < 0 and drops them through its own `keep` filter. The wrapper only suppresses the NumPy "invalid value" warning around the call.
8. **Output directory** must be empty; otherwise the run aborts.

---

## 4a. Review decisions on §4 (2026-11-22)

- Items 1–6 and 8 are accepted.
- Item 7 is accepted conditionally, pending an independent check of the copy claims.
- That check is now possible without trusting this note: `COPY_VERIFICATION_e7c7520.txt` lists, line by line, the e7c7520 source text next to the wrapper text, and the e7c7520 sources are public at that commit.

---

## 5. Run outputs (in `--out`)

- `manifest.json`: spec hash, mode, seed, cell definitions, wrapper hash, scanner hash, config hash, numpy and Python versions.
- `truth.jsonl`: one row per case (case id, spawn index, cell, block, sign, p_true, B_true, level, f0, noise, jitter).
- `analysis.jsonl`: one row per case, with no truth fields (exclusion, partial counts, winner, significant, distinguishable, mixture_like, declared, a_winner, signs, B̂ with se and interval, per-template a/p/padj/flag, W2 tests, timing; in smoke mode also the cross-check).
- `run.log`: progress every 100 cases, the s/case rate, exclusions; any parity failure.
- `hashes.json`: SHA-256 of the four files above, written after they are closed. The scorer verifies them before reading anything.
- `score.json` (scorer, full mode only): all criteria, per-cell information, the verdict, `qualifiers` (texts) and `required_qualifiers` (the full B0 §7 reporting block).

---

## 7. Copy verification (static)

`COPY_VERIFICATION_e7c7520.txt` was produced mechanically from `git show e7c7520:<path>`, comparing text only:
- **[1]** `piano_scanner.synth_note` l. 387–403 against `synth_mixture` (wrapper l. 199–215): **17/17 lines identical** after stripping indentation;
- **[2]** the body of `globex_scanner._ftest` l. 282–285 against `ftest_globex`: **4/4 lines identical**;
- **[3]** the `synth_note` frequency law (l. 377–386), listed only to show what was **not** copied; `synth_mixture` has its own law.

The reviewer can repeat this from the public repository at e7c7520 without running any project code.

---

## 6. Status of the steps

**Done:**
1. Static review of B1 (v0.1 → v0.3).
2. Smoke run 1 (wrapper v0.2): 21 cross-check failures. That led to wrapper v0.3.
3. Smoke run 2 (wrapper v0.3): passed (§8).

**Not done (each requires separate approval):**
4. **Freeze commit:** a commit to residual-first-framework, the tag `k4-recovery-b0b1-frozen` and the push, via `commit_rff_k4_freeze.ps1` on Lenovo with the files selected by SHA-256. **At the time of writing there is no commit.**
5. **B2:** the single run (seed 20261122, 5900 cases), with its own script and its own approval. **Not run, not authorised.**

---

## 8. Smoke runs (technical only; no rates computed)

| | run 1 | run 2 |
|---|---|---|
| folder (Lenovo) | `D:\research\rff_k4_smoke_20261122_035543` | `D:\research\rff_k4_smoke_20261122_040135` |
| wrapper | v0.2 `b61e9b44…` | v0.3 `75b8bbf6…` |
| scorer | `8415d1df…` | `8415d1df…` |
| environment | Python 3.13.12, numpy 2.5.2, scipy 1.18.1 | same |
| cases / excluded | 108 / 0 | 108 / 0 |
| parity guard | no mismatch | no mismatch |
| W3 cross-check | present 108/108, **21 failures** | present 108/108, **0 failures** |
| scorer structural problems | none (exit 1 because of the cross-check) | none, **exit 0** |
| speed | 0.73 s/case | 0.74 s/case |
| `truth.jsonl` | `601b3694…6eddef6` | `601b3694…6eddef6` (**identical**: generation is deterministic under the seed) |
| `analysis.jsonl` | `e2e9344c…5b9e7` | `508b7cdf…a4f25c` |
| `manifest.json` | `8798f465…0268` | `680ef76f…eb0b3` (records the wrapper hash) |
| `run.log` | `998016d6…848a` | `aa9c899b…9060` |
| `hashes.json` | `67b6225d…1f8f` | `f2ea4281…875a` |

**Runtime estimate for B2** (from run 2, which includes the smoke-only cross-check): 5900 × ≈ 0.74 s ≈ **73 min**.
