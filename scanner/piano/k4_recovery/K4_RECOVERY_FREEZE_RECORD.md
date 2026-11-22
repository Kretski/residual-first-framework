# Freeze record: blind k⁴-equivalent recovery test (B0 + B1), before B2

**Status:** proposed freeze, for approval. When it is committed and tagged in residual-first-framework, this record becomes the pre-registration of B2. **B2 itself requires a separate, explicit approval.**
**Date:** 2026-11-22.
**Base:** residual-first-framework `e7c7520d59299e0535fb101c4c13a8d6f377b9bd` (branch `main`).

## Frozen set (folder `scanner/piano/k4_recovery/`)

| file | role | SHA-256 |
|---|---|---|
| `RFF_K4_RECOVERY_AUDIT_v0.1.md` | audit A: why the test is needed | `eb7cb961569d5becda3a277ffaa1b1c4b607f82d48315952084238307952fc97` |
| `RFF_K4_RECOVERY_TEST_SPEC_B0_v0.4_CORRECTED.md` | B0: the pre-registered test | `662b2303aff00bb1000d8a0ed93be64fc5f7fdc72aae811567fb2f5810a4d042` |
| `k4_recovery_b1.py` | B1 wrapper v0.3 | `75b8bbf6517887300be52fa403fcaa14523377a7c5e4402ed60deab14e72eb4a` |
| `k4_score_b1.py` | B1 scorer | `8415d1df7a89863b9c62386d3c2f2bf11d774a6358ada88a506e33f35f24e8dc` |
| `B1_IMPLEMENTATION_NOTES_v0.5.md` | implementation notes, review decisions, smoke record | `2f58b30e46f6120e26d50ad007729b9961f5999bb224619001690e63f6fae9b9` |
| `COPY_VERIFICATION_e7c7520.txt` | mechanical check of the copied lines | `ab32613f1fcd4fa848ed68eabcfdb5ae23c2a376c7808b30de3bec930e2ac252` |
| `K4_RECOVERY_FREEZE_RECORD.md` | this record | (hash of this file, reported at commit) |

**Unchanged and verified by the wrapper at start-up:**
- `scanner/piano/piano_scanner.py`: `335ba8b8abaef4b424b36a516f55df568b9345f83cc8b6bc1c0dc1cc33c9e5d4`;
- its `config_hash()`: `5bdf34b190a4b9ef`.

## Fixed for B2 (from B0; restated, not changed)

- **Command:** `python -I k4_recovery_b1.py b2 --out <empty dir>`, then `python -I k4_score_b1.py <dir>`.
- **Seed** 20261122; `SeedSequence(20261122).spawn(5901)`, with index 5900 used for the shuffle.
- **Cases:** 5900 in 36 cells, asserted before any signal is generated.
- **Criteria, verdict ladder and stopping rule:** B0 §5, §7 and Appendix B. The scorer's thresholds are the constants at the top of `k4_score_b1.py`.
- **One run only** (B0 §6.5). At most one revision after B2, and only for a demonstrated code defect.

## Technical readiness (smoke run 2, Lenovo, 2026-11-22)

- 108/108 cases, 0 excluded, parity guard clean;
- W3 cross-check 108/108 present with 0 failures;
- scorer exit 0;
- 0.74 s/case, so B2 is expected to take about 73 min.

Smoke run 1 (wrapper v0.2) had 21 cross-check failures. The cause was the reference's finite-difference Jacobian, fixed in wrapper v0.3. The full record is in the notes, §8.

## Not frozen here

- The Lenovo helper scripts (`run_k4_smoke_lenovo.ps1`, and the B2 script still to be written) are run tools outside the repository; their hashes are reported separately.
- The review diffs (`B1_v0.1_to_v0.2.diff`, `B1_v0.2_to_v0.3.diff`, `B1_wrapper_v0.2_to_v0.3.diff`) are review aids and are not committed.
