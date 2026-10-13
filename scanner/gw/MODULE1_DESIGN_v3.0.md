# Module 1 design v3.0 — full-chain Neyman belts

**Status: DRAFT, not frozen; choices of §6 settled.** Written 2026-10-12, after the secondary IMRPhenomXPHM results and
the two post-result validation checks (V1, V2), which showed that the committed bounds are not
validated. Nothing in v1.0, v2.0 or the secondary analysis is changed by this document.

## 1. What v3.0 fixes, and why

The analysis chain refines a common time shift (±100 ms) and complex amplitude against the
reference before it scans in Lambda. The v1.0/v2.0 belts did not contain that step: they assumed
that the residual is P[s(Lambda)] + noise. At Lambda = 0 this is harmless (rule 0a passed for
IMRPhenomXPHM), but at a non-zero true Lambda the refinement absorbs part of a large dispersion
phase. V1 measured the consequence: coverage 0.549 at Lambda = ±2 sigma_lin, with estimates pulled
toward zero, most strongly for weak events. V2 showed a second, separate problem: an
IMRPhenomXPHM–SEOBNRv5PHM waveform difference adds scatter to Lambda comparable to the
statistical width (coverage 0.845 at Lambda = 0).

v3.0 builds the belts from **the full chain**, refinement included, and **includes the waveform-
model difference as a nuisance**. The belts are then checked with pre-registered acceptance rules
that include V1- and V2-type tests, not only Lambda = 0.

## 2. Scope and status of the data

- Reference model: IMRPhenomXPHM (search C01, confirmation C00 IMRPhenomXPHM-SpinTaylor), Addendum D
  pairs (24 and 20). The EOB labels are not used as references; this choice follows the v1.0 and
  v2.0 failures and is stated as such.
- The O1–O4a on-source data have already been analysed (secondary analysis). A v3.0 on-source
  result on these data is therefore **not blind** and is reported as a re-analysis with a corrected
  construction. The blind test of v3.0 is a later application, frozen before the data are opened,
  to events not analysed here (for example GWTC-5.0), registered separately.

## 3. The construction

**Trial data.** For an event, a trial at true value Lambda_t is

    d = h_M(theta) exp(i K f^3 Lambda_t) + n_G + n_cal,

with theta a posterior sample not used for the subspace of the belt (cross-sets as in v2.0: belts
of subspace A use set-B samples and vice versa), h_M generated with waveform model M (below), n_G
Gaussian noise of the event's PSD, n_cal a calibration draw. The trial goes through **the exact
chain of the analysis**: refinement against h_ref, residual, projection, profile scan on the
event's grid, giving a chi2 curve over the grid.

**Model term.** M = IMRPhenomXPHM in half of the trials and M = SEOBNRv5PHM in the other half,
generated with the label settings, so that the belts are marginalised over the two models with
equal weight. Pairs for which SEOBNRv5PHM cannot be generated with the label settings (known:
GW230627_015337, GW230924_124453) use IMRPhenomXPHM only and are listed; no settings are changed to
make them work.

**Truth grid.** Belts are built at 11 true values per cohort, Lambda_t = 0, ±0.5, ±1, ±2, ±4, ±8
times sigma_cat (the linear catalog width), with 400 trials per value and subspace. Critical values
between them are interpolated linearly in Lambda; outside ±8 sigma_cat the value at the edge is used.

**Catalog.** For each truth and trial index, the per-event chi2 curves are summed over the events
of the cohort (independent noise and samples per event); q(t) = chi2(t) − min chi2; c_t = 90%
quantile at each truth. The reported set is the union over subspaces A and B of
{Lambda : q_obs(Lambda) <= c(Lambda)}, as before.

**Cost.** One trial costs one refinement (about 0.5–1 s) plus a scan; 11 truths × 400 trials × 2
subspaces × 44 pairs is about 3.9 × 10^5 trials, i.e. about 2–3 days on the present hardware, plus
the generation of the samples (cached per pair).

## 4. Acceptance rules (all must pass before any on-source computation)

All on NEW material: off-source segments with ladder index k >= 24, none used before; posterior
samples beyond set F; the threshold of rule 0a (lower end = 0.90 − 2.576 sqrt(0.09/n)).

| Rule | Data | Injection | Statistic |
| --- | --- | --- | --- |
| 0 | synthetic, independent trials | full chain, model mixture, at the 11 truths and at 4 truths between grid points | catalog coverage at each truth (interpolation check included) |
| 0a | real off-source noise | IMRPhenomXPHM posterior sample, Lambda = 0 | per-event coverage, pooled over both cohorts |
| V1' | real off-source noise | IMRPhenomXPHM sample with Lambda = ±2 sigma_lin | pooled coverage |
| V2' | real off-source noise | SEOBNRv5PHM sample, Lambda = 0 | pooled coverage |
| 3 | synthetic | calibration-only, noise-free | as in the secondary analysis |

If any rule fails, v3.0 stops before any on-source computation, and the failure is reported.

## 5. Order

1. this draft; review; amendments before freeze recorded with their reason;
2. code (a new driver; the v1.0/v2.0 files imported unchanged), tested by reproducing a stored
   off-source chain row and on two pairs synthetically;
3. freeze: code hash, tag v3.0-module1;
4. belts and catalogs; rules 0, 0a, V1', V2', 3;
5. only if all pass: the (non-blind) on-source re-analysis, search cohort first, each result
   committed before the next;
6. a separate registration for a blind application to new events.

## 6. Choices settled before freeze (2026-10-12)

- **Truth grid and trials:** 11 true values (0, ±0.5, ±1, ±2, ±4, ±8 sigma_cat), 400 trials per value
  and subspace. A 90% quantile from 400 trials is known to about ±0.03 in coverage; the
  interpolation between grid points is checked by rule 0 at four intermediate truths, and a
  failure there stops v3.0.
- **Model weights:** 50/50 between IMRPhenomXPHM and SEOBNRv5PHM. There is no independent reason
  to trust one model more than the other; any other mixture would need a justification that we
  do not have.
- **Per-event belts of rules 0a, V1' and V2':** built from the full chain as well, so that the
  checks test exactly the construction that is used on-source. The mismatch between the checked
  and the used construction is what V1 exposed in v2.0; it is not repeated. This adds about one
  day of computation.

## 7. Amendment before freeze (2026-10-12): failed SEOBNRv5PHM samples

One rule for every pair: a trial that calls for SEOBNRv5PHM uses the SEOBNRv5PHM waveform of its sample if it was generated with the label settings, and the IMRPhenomXPHM waveform of the same position otherwise; V2' rows are skipped where the set-G sample has no SEOBNRv5PHM waveform. This replaces the pair-level wording of section 3. Generated SEOBNRv5PHM samples (sets A, B of 200 each; set G of 8): 40 pairs complete; GW190706_222641 191/193/6; GW230811_032116 198/199/8; GW230924_124453 28/21/1; GW230627_015337 0/0/0 (IMRPhenomXPHM only). Counts recorded in v3_samples.log. No settings were changed to make generation succeed.

