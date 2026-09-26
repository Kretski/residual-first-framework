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

### Amendment 4c (same day, still before any residual) — supersedes 4 and 4b

The printed offending lines showed the Amendment-4 diagnosis (mixed TOA formats)
was WRONG. All EPTA files are tempo2 FORMAT 1:
- J0751+1807, JBO.DFB.1520.tim line 2 is a normal FORMAT 1 TOA line; the file has
  no FORMAT line of its own and inherits FORMAT 1 from the master .tim (tempo2
  behaviour). Read on its own, PINT fell back to the Princeton format and failed.
- J0613-0200, EFF.EBPP.2639.tim line 66 starts with `C??` — a TOA commented out by
  EPTA (not a bare `C`, so Amendment 4b missed it).
Reader rule now: FORMAT inherited from including files; a line is kept only if it
is a valid FORMAT 1 TOA line (numeric frequency, MJD, error) or a known tempo2
command; all other lines are dropped and the number dropped is printed per file.

### Amendment 4d (same day, still before any residual)

Remaining failures (J0751+1807 JBO.DFB.1520.tim, J0613-0200 LEAP.1396.tim) were
standard FORMAT 1 lines parsed by PINT as Princeton format. Cause: the 4c reader
passed FORMAT only from parent to child, whereas in tempo2 (and in PINT's own
INCLUDE handling) FORMAT is global state in reading order — once set by any
earlier file it applies to all following files. Fix: FORMAT tracked as global
state; FORMAT lines kept in place; a copy gets a leading 'FORMAT 1' only if that
state was already set when the file started. Dropped lines remain only EPTA's own
commented-out TOAs, printed per file.

### Amendment 4e (same day, still before any residual)

Root cause of the remaining J0751+1807 / J0613-0200 failures: the TOA lines begin
with a blank, and PINT's per-line format detection classifies a line that starts
with a blank and has '.' at column 41 as fixed-column Parkes format BEFORE applying
FORMAT 1 (e.g. ' J090811_044909.NEFTp 1520.00000000 55054.2…': column 41 is the
MJD decimal point). FORMAT 1 is whitespace-delimited, so leading blanks are
stripped from TOA lines read under FORMAT 1. No MJD, frequency or error value is
altered.

## Amendment 5 — invalid gate run, runner fix (before any residual)

Gate run #2 (after Amendment 4e) is INVALID and not counted: WSL filesystem errors
on /mnt/c (Errno 5 "Input/output error", Errno 14 "Bad address") made 22 of 25
pulsars fail to load; the runner wrongly treated these as pulsar exclusions and
completed the gate on 3 pulsars (primary H₀ rate 0.000, meaningless with 3 pairs).
The gate JSON was not written (same I/O error), so no real-data run was possible.
Fix: any OSError now aborts the run; only PINT loading errors exclude a pulsar.
All further runs are executed from the Linux filesystem (~/rff), outputs copied back.

## Gate result (valid run #3, from ~/rff, after Amendment 5)

25/25 pulsars loaded, 0 excluded; 200 H₀ trials. Primary
`hd_perp_mono_dipole:scramble` H₀ rate **0.045** ≤ 0.09 → **GATE PASSED**.
Other rates: hd_sky_scramble 0.035, hd_perp_mono:scramble 0.045 (calibrated);
monopole 0.135, hd_perp_mono:shift 0.120, hd_perp_mono_dipole:shift 0.120,
cross_correlation 1.000 (NOT calibrated on EPTA sampling → reported as
uninterpretable). Record: `epta_gate_DR2full.json`.

## Amendment 6 — pre-fit vs re-fit decided, BEFORE any residual is computed

Decision **A**: because PINT's TCB→TDB conversion is approximate (PINT warns the
model "should be re-fit"), each pulsar's released model is re-fitted once with
PINT's downhill WLS fitter (max 10 iterations) over the parameters marked free
in the released .par. Post-fit residuals are used (as for NANOGrav, whose
released residuals are also post-fit). A pulsar whose fit raises an error other
than a max-iteration warning is excluded and listed. Per-pulsar QC (number of free
parameters, pre-/post-fit RMS, reduced χ², fit status) is printed and saved;
QC is not used to select pulsars. Note: the gate simulator applied only a
quadratic removal, not the full timing-model fit; the sky-scramble null of the
primary test does not depend on the temporal model.
Programming errors (AttributeError, NameError, TypeError, ImportError) abort the
run instead of excluding pulsars; after fixing, the run is repeated in full.

## Method check M1 — nuisance leakage in EPTA geometry (specified BEFORE running)

Purpose: the NANOGrav validation showed the primary test leaks under an injected
dipole (0.10–0.20). This check measures leakage with EPTA sampling, errors and
positions. It uses synthetic data only (no EPTA residuals), so it can be run after
the real-data results without affecting them; its criterion is fixed here first.

Design (`run_epta.py dipole`, for DR2full and DR2new): H₀ noise as in the gate
plus an injected common red process (β = 13/3): H₀ (100 trials), dipole 0.5 µs and
1.0 µs (100 each, random direction per trial), monopole 0.5 µs (50), HD 1.0 µs (50,
power reference). n_null = n_scramble = 200 per trial, seed 1.
Criterion: the primary test is called ROBUST for this geometry if its p ≤ 0.05
rate is ≤ 0.10 under both dipole amplitudes and under the monopole; otherwise the
leakage rate is reported as a limitation of the EPTA results. The amplitudes are
stress-test values, not estimates of real ephemeris errors.

---

# RESULTS (reported in full, regardless of outcome)

## R1. Primary result — EPTA DR2full

Gate (run #3, valid): 25/25 pulsars, 200 H₀ trials, primary H₀ rate **0.045** → PASSED.
Real run (Amendment 6 WLS re-fit): 24 pulsars; **J0030+0451 excluded** (PINT
`MissingTOAs`: JUMP9 on -sys JBO.DFB.1520 has no TOAs). 301 × 30-d bins, 276 pairs,
median overlap 158 bins.

**Primary `hd_perp_mono_dipole:scramble`: p = 0.0135 → p ≤ 0.05.**
Pre-registered wording: "EPTA DR2 shows an HD-like correlation beyond monopole +
linear dipole, consistent with the NANOGrav result" — with the caveats below.

| test | null | calibrated on DR2full? (gate H₀ rate) | p |
|---|---|---|---|
| hd_perp_mono_dipole | scramble | yes (0.045) | **0.0135** (primary) |
| hd_perp_mono | scramble | yes (0.045) | 0.0080 |
| hd_sky_scramble | scramble | yes (0.035) | 0.0095 |
| monopole (projection −1.32) | shift | **no** (0.135) | 0.868 — uninterpretable |
| hd_perp_* | shift | **no** (0.120) | 0.005 / 0.0085 — uninterpretable |
| cross_correlation (z²) | shift | **no** (1.000) | 0.105 — uninterpretable |

QC caveat: several pulsars have implausibly large post-fit residuals / χ²_red
(J1744-1134 7208, J1857+0943 2067, J1713+0747 321, J1022+1001 254, J1730-2304 124),
most likely unmodelled DM / solar-wind (chromatic) noise and low-frequency WSRT data.
Per Section 3 QC was not used to select pulsars.

## R2. Secondary result — EPTA DR2new

Gate: 25/25 pulsars, primary H₀ rate **0.040** → PASSED. Real run: 25 pulsars,
0 excluded, 126 × 30-d bins, 300 pairs, median overlap 107. QC clean
(χ²_red 0.70–5.81; largest J1600-3053 5.81, J1843-1113 5.58, J1909-3744 5.40).

**Primary test on DR2new: p = 0.0065 → p ≤ 0.05.**

| test | calibrated on DR2new? (gate H₀) | p |
|---|---|---|
| hd_perp_mono_dipole:scramble | yes (0.040) | **0.0065** |
| hd_perp_mono:scramble | yes (0.065) | 0.0095 |
| hd_sky_scramble | yes (0.060) | 0.011 |
| monopole (+3.64) | **no** (0.150) | 0.014 — uninterpretable |
| hd_perp_*:shift | **no** (0.105–0.115) | 0.0005 — uninterpretable |
| cross_correlation (z²) | **no** (1.000) | 0.0005 — uninterpretable |

The signal is stronger in the cleaner DR2new subset, which argues against the
DR2full QC problems as its origin. DR2new is a subset of DR2full: not independent.

## R3. Method check M1 (nuisance leakage, EPTA geometry)

Rate of p ≤ 0.05 (H₀/dipole 100 trials, monopole/HD 50):

| injection | DR2full: primary / hd_sky / hd_perp_mono | DR2new: primary / hd_sky / hd_perp_mono |
|---|---|---|
| H₀ | 0.040 / 0.040 / 0.050 | 0.050 / 0.040 / 0.050 |
| dipole 0.5 µs | 0.060 / 0.090 / 0.070 | 0.020 / 0.110 / 0.080 |
| dipole 1.0 µs | 0.010 / 0.100 / 0.110 | 0.000 / 0.160 / 0.160 |
| monopole 0.5 µs | 0.000 / 0.020 / 0.020 | 0.020 / 0.020 / 0.040 |
| HD 1.0 µs (power) | 0.28 / 0.28 / 0.32 | 0.52 / 0.48 / 0.50 |

Criterion (≤ 0.10 under both dipoles and the monopole): **ROBUST for both variants.**
Removing the 6-dim dipole subspace is what makes the difference (0.10–0.16 → 0.00–0.01
at 1 µs). In NANOGrav geometry the same test leaked 0.10–0.20; robustness is
geometry-dependent.

## R4. Deviations and record-keeping notes

- Amendments 1–6 and 4b–4e were all committed before the corresponding runs.
- Commit 2b3d95f is labelled "Amendment 6" but contained only the gate JSON; the
  Amendment 6 files were committed in ca75402, still before any residual.
- Gate run #2 is invalid (Amendment 5) and not counted.
- The EPTA, DR2full and DR2new results are not independent of each other or of
  NANOGrav (shared pulsars, same GWB realisation, DR2new ⊂ DR2full). No combined
  p-value is reported.

## R5. Limitations

No DM / chromatic / solar-wind noise model; no EFAC/EQUAD from the EPTA noise files;
PINT TCB→TDB conversion is approximate (hence the re-fit); some Effelsberg/LEAP TOAs
fall outside the shipped clock files; power of the primary test is modest (0.28–0.52
at 1 µs HD); M1 amplitudes are stress values, not estimates of real ephemeris errors.
This analysis recovers the known PTA HD signal with a simple residual-level method;
it is not a new detection and has no bearing on modified-gravity models.
