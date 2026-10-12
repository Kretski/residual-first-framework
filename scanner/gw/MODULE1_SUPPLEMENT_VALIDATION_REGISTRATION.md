# Module 1 supplement — two post-result validation checks: registration

Date: 2026-10-12. Written after the supplementary IMRPhenomXPHM on-source results were committed
(tags supplement-search-result, supplement-confirm-result) and before either check is run.

## 1. Status

These are **validation checks of the method**, decided after the on-source results. They use
**off-source data only**, on segments and posterior samples never used before. They **cannot
change** the reported confidence sets. Their outcome is reported next to the sets as a statement
of how far the sets can be trusted, whatever that outcome is. They address two gaps:
(1) rule 0a tested coverage only at a true value Lambda = 0;
(2) no check tested whether a waveform-model difference can masquerade as Lambda, although the
EOB failures of v1.0 and v2.0 show that model-related mismatch is real.

## 2. Common setup (both checks)

- Model analysed: IMRPhenomXPHM, v2.0 construction, frozen supplement code
  (tag supplement-xphm-frozen-r2) for the chain, belts and seeds; pairs of Addendum D
  (search 24, confirmation 20).
- Segments: NEW off-source segments, the first 8 valid ones with ladder index k >= 16
  (alternating sides, same spacing and guard rule), none overlapping any segment used by the
  v1.0 or v2.0 off-source stages or the diagnostic.
- Injected samples: set F = the 8 usable posterior samples following set E in the fixed
  permutation (seed 20261001), disjoint from A, B, C, E. One set-F sample per segment.
- Chain: as in the off-source stage (common time shift ±100 ms and complex amplitude against
  the reference, profile scan), with each event's own grid and own Neyman belt (v2.0 belt
  construction, 2000 trials, seed SEED + 11 + 1000 k).
- Acceptance: per check, both cohorts pooled (about 352 segments), lower end of the 99% binomial
  interval around 0.90, as rule 0a.

## 3. Check V1 — coverage at a non-zero true Lambda in real noise

- Injection: the set-F waveform multiplied by the dispersion phase of the analysis template
  at a true value Lambda_t, exp(i K f^3 Lambda_t) with the event's K, i.e. exactly the signal
  model the estimator assumes.
- True values: Lambda_t = +2 sigma_lin and −2 sigma_lin of the event (grid points ±40), half of
  the segments each (alternating), so that the non-linear regime is included.
- Statistic: covered if q(Lambda_t) <= c(Lambda_t) of the event's own belt.
- Outcome: pooled coverage within or above the 99% interval -> "coverage established away
  from zero at ±2 sigma_lin"; below -> "coverage away from zero not established", reported next
  to the sets together with the measured coverage.

## 4. Check V2 — can a waveform-model difference look like Lambda?

- Injection: a waveform with the parameters of the set-F sample but generated with
  **SEOBNRv5PHM** (same sky position, polarization, time, reference frequency and sampling as
  the label; Lambda = 0). Analysis with the IMRPhenomXPHM reference and chain, unchanged.
- Statistic: covered if q(0) <= c(0); also the scan estimate of each segment in units of the
  event's sigma_lin.
- Outcome: pooled coverage of Lambda = 0 within or above the 99% interval -> "no leakage of
  the XPHM–SEOBNRv5PHM difference into Lambda detected at this sensitivity"; below -> "the model
  difference leaks into Lambda", reported with the median and spread of the estimates in
  sigma_lin, and stated next to the sets as an unresolved model systematic.
- Pairs whose SEOBNRv5PHM waveform cannot be generated with the label settings are skipped and
  listed; they are not replaced.

## 5. Order and what does not change

1. this registration (commit);
2. the code, tested first by reproducing one v2.0 off-source row (no injection change), then
   frozen with its own hash before the checks are run;
3. V1, then V2, each committed when finished.
No rule, threshold, pair, segment rule or construction is changed after either result. The
on-source sets stay as committed. Neither check is a new on-source computation.
