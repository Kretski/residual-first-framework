# Shammass R0 — gate code audit v1

**Audited code:** `shammass_r0_v1_0_3_reviewed.py`, SHA-256 `007d56ad92cca87959cab349231c751c2532f0c921c64f4974a1ae40036c06d9` (matches `code_sha256` in `r0_gate_results.json`).  
**Against:** `SHAMMASS_R0_PREANALYSIS_SPEC_v1.0.3_DRAFT.md`, SHA-256 `24d487ab6799cd90868cc5fabe378101b0bd6d705f1b60e6de1022be95fd6914`.  
**Also inspected:** local commit `b70f1188b7eb7c24b6fc617b6a4461ff248ea8ac` (tag `r0-v1.0.3-pre-gate`, repository `Desktop\Q6\Z\R0_IMPLEMENTATION`, no remote), which contains `shammass_r0_v1_0_3.py`, SHA-256 `18630f51f12561f168ef755494ea5dfb225285e69c2fb3364a874a6933790a39`.  
**Method:** static reading of both scripts against the spec; independent recomputation of G1/G2 from `r0_gate_results.json`. The gate was not re-run.

## 1. Run history (facts)

| Item | Fact |
|---|---|
| Spec v1.0.3 approval | §13 checkboxes unchecked; status `DRAFT / NOT APPROVED`. The committed README states this. |
| Code at tag `r0-v1.0.3-pre-gate` | `18630f51…0a39` |
| Code that produced the gate result | `007d56ad…06d9` — **not committed, not tagged** |
| Difference | 224 insertions, 65 deletions |
| Gate results | not committed |
| Order of edits vs gate run | not established (file timestamps not yet recorded) |

## 2. Conformance of the audited code to v1.0.3

| Spec | Requirement | Finding |
|---|---|---|
| §1 | input hashes and counts, abort on mismatch | conforms |
| §1, blindness | gate reads experimental `k` only | conforms (`load_inputs(..., curve_only=True)`; CSV reader parses only requested columns) |
| §3.1–3.2 | sort, merge \|Δk\|<1e-9, record non-increasing | conforms |
| §3.4 | domain rule: a point outside the domain or in a segment > 0.5 µm⁻¹ is **not used** | **deviation D1**: in synthetic replicates the whole replicate is marked `NOT_TESTABLE` instead of dropping the point. Real-data path: all 25 points pass (curve-only record), so no effect there. |
| §4.1–4.2 | linear primary, PCHIP without extrapolation | conforms |
| §5.2 | \(v_i=(1+(2.0294 s_i)^2)(1+(1-t_i)^2+t_i^2)\) | conforms |
| §5.3 | extend δ grid if LOO RMS > 8 px | not implemented; not triggered (LOO 0.56 and 1.03 px) |
| §6 | column scaling, lstsq, cond > 1e8 → ill-conditioned, partial F (1, N−3) | conforms |
| §7 | Holm per series at 0.01; q* = smallest raw p among rejected, ties → smaller q; three-valued sign | conforms (adjusted-p form is equivalent to step-down) |
| §8.1 | truth \(F(k_i)+A\,T_q(k_i)\) on the observed \(k\), prepared curve | conforms |
| §8.2 | pixel-space noise on both axes; calibration jitter N(0,1) px on X1, X2, Y1, Y2; same jittered map for both; re-prepare | conforms. WPD `XYAxes` map reproduces all four exports to ≤ 1.8e-15. |
| §8.2 | draw order | conforms (calibration 8 values in stated order; then experimental points, then curve nodes; x,y interleaved per point) |
| §8.3 | \(\max_i\|A T_q\| = m\,\delta\,a_f\) | conforms |
| §8.4 | 10,000 / 500 replicates; `SeedSequence(20261125)`, spawn order series → δ → (null, q × m × sign) | conforms |
| §8.5 | G1, G2, G3, NOT_ESTIMABLE, verdict logic | conforms; **reporting nit D2**: when \(m^*\) does not exist, G3 is reported `false` instead of "not evaluated" |
| §8.6 | outputs written before printing | conforms |
| §1 (implicit) | `wpd.json` identity | **minor D3**: hash recorded, not enforced against an expected value; recorded value equals R0 §11.1 |

## 3. Differences from the committed version that bear on validity

1. **Blindness.** The committed `gate()` loaded the experimental CSVs in full, including frequencies (they were not used in the simulation). The audited version loads experimental `k` only. The audited version is the correct one.
2. **G1.** The committed version added a criterion not in the spec (`not_testable == 0`); the audited version applies §8.5 literally. With the committed rule, G1 would have failed at δ = 4 for both series.
3. Additional records (input/spec/code hashes, jitter summary, conversion record) and a validator for the real-data stage. Not run.

## 4. Robustness of the verdict

Recomputed independently from `r0_gate_results.json` with exact Clopper–Pearson intervals:

- G1 passes for both series, and still passes if every per-replicate `NOT_TESTABLE` were counted as a false detection.
- G2 fails by a wide margin: the largest CP lower bound of the q = 2 detection rate over all δ ∈ {1, 2, 4}, m and signs is 0.122 (large) and 0.034 (small), against the required 0.80 (at most 76/500 and 26/500 detections).

The verdict `NOT_TESTABLE` for both series therefore does not depend on D1, D2, D3, on the G1 difference between the two code versions, or on whether the code changes preceded the run.

## 5. Required records before packaging

1. Record that the gate ran with `007d56ad…06d9`, not with the tagged `18630f51…0a39`, and that spec v1.0.3 was not formally approved at run time (deviation under §11).
2. Record file timestamps of `shammass_r0_v1_0_3_reviewed.py` and `r0_gate_results.json` to establish order.
3. Record D1–D3 as implementation deviations that do not affect the verdict.
4. Commit the audited code, both test files and the gate results together, with this audit, in a new commit; do not move or rewrite the existing tag.
