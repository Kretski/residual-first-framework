# Pre-registration: independent replication of the V2.1 HD test on EPTA DR2

Author: Dimitar Kretski
Status: **to be committed and tagged BEFORE any EPTA residual is computed or inspected.**
Frozen analysis code: `common_residual_search_v2.py` (V2.1), SHA-256 recorded at tag time:
`089a94e073c962eb332779e0cc50807ebbb30a6e3f02bb49bced0226497ec2b2`
(git tag `v2.1-epta-prereg`, commit ff05b3e)

## 1. Background (what was already seen)

On NANOGrav 15-yr (v2.1.0, epoch-averaged, 68 pulsars, 30-d bins) the V2.1 pipeline gave:

| test | p (NANOGrav) | status |
|---|---|---|
| cross_correlation (z²) | 0.045 | **invalid** — FPR 0.69 on real sampling (see diagnostics_shared_structure.py) |
| monopole | 0.035 | calibrated |
| hd_sky_scramble | 0.034 | calibrated; leaks under dipole (0.60) |
| hd_perp_mono_dipole:scramble | 0.046 | calibrated (0.04); robust to monopole (0.05); partial dipole leakage (0.10–0.20) |

The last test was added after seeing V2 results (exploratory on NANOGrav).
Its choice as the primary EPTA test is therefore informed by NANOGrav, but EPTA
data have NOT been seen.

## 2. Data

EPTA DR2, doi:10.5281/zenodo.8300645, file EPTA-DR2.zip, MD5 307b2b99cd352876b07409899448af8b.
Primary dataset: the full 25-pulsar release as distributed. If the release contains
separate variants (e.g. "DR2full" and "DR2new"), DR2full is primary and DR2new is
secondary; both are reported.

Residuals: computed with PINT from the released .par/.tim files, using the released
timing model unchanged (no re-fit beyond what PINT needs to evaluate residuals).
Pulsars that PINT cannot load are listed and excluded; no other exclusions.

## 3. Analysis (frozen)

`common_residual_search_v2.py` V2.1 with defaults: 30-d bins, min_overlap = 12,
mirror-shift temporal null, quadratic detrend, n_null = n_scramble = 2000, seed = 0.
No `--min-span`, no other options.

**Primary test (single, no multiplicity correction):**
`hd_perp_mono_dipole:scramble`, one-sided, α = 0.05.

Secondary (reported, not used for the decision): all other V2.1 tests.
`cross_correlation` is reported but declared uninterpretable in advance.

## 4. Calibration gate (run BEFORE the real-data run)

`--validate --real-gaps` on EPTA sampling/errors/positions, H₀ only, 200 trials.
If the H₀ rate of the primary test is > 0.09 (≈ 3σ above 0.05 for 200 trials),
the real-data run is NOT performed and the calibration failure is reported instead.

## 5. Decision rule and reporting

- p ≤ 0.05: "EPTA DR2 shows an HD-like correlation beyond monopole + linear dipole,
  consistent with the NANOGrav result" — with the dipole-leakage caveat.
- p > 0.05: "not replicated at α = 0.05". Reported with equal prominence.
- EPTA and NANOGrav share many pulsars and the same GWB realisation; agreement tests
  independence of *systematics*, not of the signal.
- All outputs (JSON provenance, validation logs) are published regardless of outcome.
- No further variants are run on EPTA before this result is published.

## Amendment 1 — technical, recorded BEFORE any EPTA residual was computed

Made after inspecting only the release's folder structure, README, one .par header
and one .tim header (output of `fetch_epta_dr2.py`); no TOAs were processed.

1. **Dataset names.** The release contains four variants: DR2full, DR2full+,
   DR2new, DR2new+ ("+" = combined with InPTA DR1 for 10 pulsars).
   Primary = `DR2full`; secondary = `DR2new`. The "+" variants are NOT analysed
   (they add a different PTA's data and would not be an EPTA-only replication).
2. **Numerical precision.** PINT on native Windows reports "platform does not
   support extended precision floating-point". Residuals will be computed only
   on a platform with 80-bit long double (Linux / WSL2), verified with
   `numpy.finfo(numpy.longdouble).eps < 1e-18` and recorded in the output.
3. **Model conversions required by PINT.** The released .par files use
   `UNITS TCB` and `BINARY T2`. PINT's own TCB→TDB and T2→equivalent-binary
   conversions are applied; these are re-expressions of the released model, not
   re-fits. Any pulsar for which conversion fails is listed and excluded.
4. **Clock files.** The corrected NRT clock file shipped in the release is used,
   as the release README instructs.
5. **Residual type.** PINT pre-fit residuals with the released (converted) model,
   epoch-averaged per observing day and backend before the frozen V2.1 binning.

## Amendment 2 — gate simulator span, recorded BEFORE the gate was run

The frozen module's `synthetic_array()` hard-codes a 16-yr time grid (NANOGrav).
EPTA DR2 spans up to ~25 yr, so beyond 16 yr its simulated red noise would be
constant. The calibration gate therefore uses `run_epta.synthetic_h0_array()`:
the same H₀ noise model (real TOA errors, 1 µs β = 3 red noise on every third
pulsar, weighted quadratic removal) on a grid spanning the actual data span.
The analysis itself (`common_residual_search_v2.py`) is unchanged; the gate uses
200 trials with n_null = n_scramble = 200 per trial. The real-data run uses the
frozen defaults (n_null = n_scramble = 2000, seed = 0).

## Amendment 3 — TOA reading, recorded BEFORE the gate and before any residual

Found by loading J1909-3744 TOAs only (no residuals): PINT treats an `END` line in
an INCLUDEd .tim file as the end of ALL input, whereas tempo2 (used by EPTA) stops
only that file. Result: 183 of ~2800 TOAs read, all NUPPI data missing.
Fix: `run_epta.merge_tim_tempo2()` flattens the INCLUDEs with tempo2 semantics
(END ends the current file only; commented 'C' lines dropped) before PINT reads them.
Epoch averaging (Amendment 1.5) uses the `-group` flag as "backend" (the NUPPI
`-sys` flag differs per frequency channel), falling back to `-sys`, then observatory.
The pre-fit vs. re-fit question (PINT warns its TCB→TDB conversion is approximate)
is NOT decided here; it will be recorded as a separate amendment before the
real-data run, and before any residual is computed.

## Amendment 4 — TOA reader bug fix, recorded BEFORE any residual

The first gate run (17/25 pulsars, primary H₀ rate 0.035, PASSED) excluded 8 pulsars
(J0613-0200, J0751+1807, J1012+5307, J1022+1001, J1640+2224, J1738+0333,
J1744-1134, J1911+1347) because of a bug in the Amendment-3 reader: merging all
INCLUDEd files under one global `FORMAT 1` header forced older fixed-column
(Princeton/Parkes-style) TOA files to be parsed as tempo2 format. These are not
"pulsars PINT cannot load" in the sense of Section 2, so they may not be excluded.
Fix (`split_tim_tempo2`, `read_leaf`): each INCLUDEd file is truncated at its own
END and read by PINT separately with its own format; the TOAs are then merged
(`pint.toa.merge_TOAs`). The gate is RE-RUN on all loadable pulsars and the re-run
is the one that counts; the 17-pulsar run is reported for completeness.

Note from the 17-pulsar gate: on EPTA sampling the shift-null tests are NOT
calibrated (monopole 0.145, hd_perp_*:shift 0.15–0.16, cross_correlation 0.95);
only the sky-scramble tests are (0.035–0.07). Secondary shift-null results on
EPTA will therefore be reported as uninterpretable.

### Amendment 4b (same day, still before any residual)

`check J0613-0200` with the Amendment-4 reader failed on EFF.EBPP.2639.tim: a line
whose first token is a lower-case `c` (a tempo2 comment — tempo2 commands are
case-insensitive, as the lower-case `end` lines in the same release show) was
parsed by PINT as a TOA. Fix: lines starting with `C` or `c` (or `#`) are dropped
from the per-file copies. On any remaining parse failure the offending line is
printed.
