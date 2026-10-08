# Module 1 design v2.0 — Neyman belts that include the posterior/reference mismatch

**Status: DRAFT, not frozen.** Written 2026-10-08, after Module 1 v1.0 stopped at rule 0a
(tag `v1.0-stopped-0a`) and after the pre-specified diagnostic of that failure
(`DIAG_0a_SPEC.md`, `DIAG_0a_REPORT.md`): no mechanism was identified, and the descriptive
sampler audit (`diag_0a_sampler_audit.txt`) shows that waveform model and sampler are fully
confounded in the public O1–O3 release. Nothing in v1.0 is changed by this document; v1.0
remains stopped and its records stand.

## 1. Motivation

v1.0 built its Neyman belts from Gaussian noise plus calibration errors, i.e. conditional on
the true signal being exactly the reference waveform h_ref. In real off-source noise with
injections of posterior samples (set C), the SEOBNRv4PHM search labels under-covered
(155/200; pooled EOB 300/360 = 0.833 < 0.859), while IMRPhenomXPHM passed (313/352 = 0.889).
The diagnostic could not attribute the failure to a mis-centred reference or to time offsets,
and the data cannot separate model from sampler.

v2.0 therefore does not select a model or a sampler after the failure. It changes the
construction so that it no longer assumes h_true = h_ref: the belts are **constructed to
include the posterior/reference mismatch as a nuisance uncertainty, and their coverage is then
checked empirically on new data.** Adding this term does not by itself guarantee coverage;
the validation of §5 does.

## 2. What is kept from v1.0 (unchanged)

- data, event lists, cohorts and every selection rule, including the pair-level exclusions of
  Addendum D (search IMRPhenomXPHM 24, SEOBNRv4PHM 25; confirmation 20 and 20);
- both waveform models in the primary test (no post-hoc removal of SEOBNRv4PHM);
- reference point (stored maximum-likelihood sample), subspaces A and B, calibration as a
  covariance, profile scan in Lambda, Neyman construction, union of the A and B sets;
- the cached geometry of v1.0 (no rebuild) and the stage-2 checks (rules 1, 1a, 2a), whose
  inputs do not depend on the belt construction;
- every stopping rule and threshold of v1.0 §11, applied to the v2.0 belts.

## 3. The single statistical change: belt trials include a mismatch term

For an event–model pair, every belt and expected-interval trial uses

    n_trial = n_Gaussian + n_cal + P [ h(theta_D) − h_ref ],

where theta_D is a posterior sample of set D, h(theta_D) is generated and projected exactly as
the off-source injections (own sky position, polarization and time; calibration mean factor),
and P is the same projection as for the noise. Trial j uses D sample (j mod 400).

**Set D (fixed before any v2.0 computation):**
- exactly 400 usable samples per event–model pair;
- drawn from the same fixed permutation of the published posterior samples (seed 20261001;
  the published samples are equally weighted draws, so the selection is uniform over them);
- D starts at position 2n + f_AB + 100 of that permutation (n = 800 per set, f_AB = generation
  failures of sets A and B in the v1.0 build), i.e. strictly after A, B and the set-C samples
  of v1.0 (at most 8 plus failures), with a margin of 100 positions;
- D is not chosen by event, coverage or any result; generation failures are skipped and
  recorded.

D enters the belts only. It is never used to validate them.

## 4. Catalog construction

As in v1.0 (module1_catalog: common grid, sums over events, profile-likelihood ordering),
with the per-event trials of §3. The catalog belt is the sum of the per-event mismatch-
including trials. The union of the A and B sets is kept, with both belts built from the
same noise and the same D samples.

## 5. Validation on new data (rules 0 and 0a of v1.0, unchanged thresholds)

- **Rule 0 (construction):** coverage of the v2.0 belts at the five true values, from
  independent trials that also include the mismatch term, with D samples taken in a
  different order (fixed shift of 200) so that the check trial and the belt trial never pair
  the same noise with the same D sample.
- **Rule 0a (off-source):** NEW segments and NEW injections, disjoint from everything used
  by v1.0 and its diagnostic:
  - segments: the first 8 valid segments with ladder index k >= 8 (alternating sides, same
    spacing and guard rule as v1.0), none overlapping any v1.0 off-source segment;
  - injections: set E = the 8 usable samples following set D in the same permutation;
  - same chain as v1.0 (time shift ±100 ms, complex amplitude, profile scan, q(0) against the
    event's own v2.0 critical value at Lambda = 0);
  - acceptance per model, both cohorts pooled, below the lower end of the 99% binomial
    interval around 0.90 -> stop before any on-source estimate. Reference-waveform
    injections remain a diagnostic.
- If SEOBNRv4PHM again fails rule 0a, v2.0 stops as well; this is reported as such.

## 6. Order and registration

1. this document (draft) and the v2.0 code;
2. synthetic test of the code: coverage with posterior-sample injections in Gaussian noise
   for two events (one per model), before freezing;
3. freeze: code hash in FREEZE_v2.0.json, tag v2.0-module1, before any v2.0 belt is built;
4. v2.0 belts and catalogs (geometry reused from the v1.0 cache);
5. rule 0, then rule 0a on the new segments and set E;
6. the reporting layer applied to all rules;
7. only if no stopping rule fired: the first on-source estimates, search cohort only; the
   confirmation cohort after the search result is committed.

## 7. Expected cost and consequences

- set D: 400 waveforms per pair (about 25 min per SEOBNRv4PHM pair, minutes for the others);
  belts rebuilt; 1–2 days of computation in total.
- Intervals are wider where the mismatch is large; for IMRPhenomXPHM, where the mismatch is
  small, they should stay close to v1.0. The width is reported with the result.

## 8. What v2.0 does not claim

It does not identify the cause of the v1.0 failure, and it does not make any model or sampler
preferred. It makes the interval construction valid under the posterior/reference mismatch
present in the published samples, and then tests that validity on independent data.
