# Shammass R0 — synthetic gate run record

**Repository:** local `R0_IMPLEMENTATION` (Lenovo, `Desktop\Q6\Z\R0_IMPLEMENTATION`), no remote.  
**Spec:** `docs/SHAMMASS_R0_PREANALYSIS_SPEC_v1.0.3_DRAFT.md`, SHA-256 `24d487ab6799cd90868cc5fabe378101b0bd6d705f1b60e6de1022be95fd6914`.

## Timeline (2026-11-25, local time, from file LastWriteTime and git)

| Time | Event | Evidence |
|---|---|---|
| 09:15:05 | Commit `b70f1188…` with `shammass_r0_v1_0_3.py` (`18630f51…0a39`), spec v1.0.3, audit and curve-only results; tag `r0-v1.0.3-pre-gate` | git |
| 09:27:29 | `shammass_r0_v1_0_3_reviewed.py` (`007d56ad…06d9`) and `test_shammass_r0_v1_0_3_reviewed.py` written | LastWriteTime |
| 09:40:57 | `test_shammass_r0_gate_validation.py` (`43f694a1…dd23`) written | LastWriteTime |
| 09:52:15 | `R0_REVIEWED_RESULTS/gate/r0_gate_results.json` (`fd94f4b9…f4e3`) written; records `code_sha256` = `007d56ad…06d9` | LastWriteTime, file content |

File timestamps are evidence of order, not proof: they can change on copy or edit. The code hash recorded inside the results file is the stronger link.

## Deviations under spec §11

1. **Spec not formally approved at run time.** §13 checkboxes of v1.0.3 were unchecked when the gate ran. The spec text itself was not changed afterwards (hash above).
2. **Code not committed before the run.** The gate ran with `007d56ad…06d9`, written 12 minutes after the pre-gate tag and before the results. The tagged version `18630f51…0a39` was not used for the run.
3. **Changes between tagged and run code** (224 insertions, 65 deletions), all written before the results file:
   - gate loads only experimental `k` (the tagged version loaded experimental frequencies, unused);
   - G1 applied literally per §8.5 (the tagged version added `not_testable == 0`, not in the spec);
   - added hash, environment, conversion and jitter records; real-stage validator (not run).
4. **Implementation deviations found in audit** (`R0_GATE_CODE_AUDIT_v1.md`): D1 domain rule applied per replicate, D2 G3 reported `false` when not evaluated, D3 `wpd.json` hash recorded but not enforced.

## Result and its robustness

Verdict: `NOT_TESTABLE` for the large and the small series. G1 passes; G2 fails (largest CP lower bound of q = 2 detection 0.122 large, 0.034 small, required 0.80). The verdict does not depend on any deviation listed above (audit §4).

## Not done

No real-data stage was run. No experiment-minus-GPE residual was computed. No value of ΔΛ or Λ was estimated.
