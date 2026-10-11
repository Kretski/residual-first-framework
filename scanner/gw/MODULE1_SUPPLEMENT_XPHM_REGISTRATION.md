# Module 1 — supplementary IMRPhenomXPHM-only analysis: registration

Date: 2026-10-11. Written before any on-source quantity of Module 1 has been computed and
before the code that computes it exists.

## 1. Status and how this analysis must be read

The primary Module 1 analysis stopped twice at its own rule 0a, both times because of the EOB
models: v1.0 (tag `v1.0-stopped-0a`, EOB pooled off-source coverage 0.833) and v2.0 (tag
`v2.0-stopped-0a`, 0.849; threshold about 0.859). IMRPhenomXPHM passed rule 0a in both versions
(0.889 and 0.890).

This analysis is **selected after the EOB failures**. It is a **secondary, supplementary**
analysis. It does not replace the stopped primary analysis, it does not carry the two-model
robustness of the registered design, and its result is never presented as the main result of
the work (not in the title, not as the headline of the abstract). Every place that reports its
number states that the companion waveform model failed validation and that the selection was
made after that failure.

## 2. What is analysed

- Waveform model: IMRPhenomXPHM only (search labels C01 IMRPhenomXPHM; confirmation labels C00
  IMRPhenomXPHM-SpinTaylor).
- Construction: **v2.0 exactly as frozen** (tag `v2.0-module1`; belts including the
  posterior/reference mismatch with the cross-set rule). v1.0 is not computed for this
  supplement, so that there is no choice between constructions after seeing on-source numbers.
- Pairs: the primary lists of Addendum D, unchanged: `primary_search_xphm.txt` (24 events) and
  `primary_confirm_xphm.txt` (20 events).
- Geometry, calibration covariance, profile scan, Neyman construction, union of the A and B
  confidence sets, grid (±8 × the smallest sigma_lin of the cohort, step 0.05): as in the
  v2.0 catalogs `v2_catalog_search_xphm.json` and `v2_catalog_confirm_xphm.json`.

## 3. Conditions that must hold before any on-source computation

All v1.0 §11 stopping rules, evaluated for IMRPhenomXPHM alone, with the v2.0 belts where the
belts enter: rule 0 (catalog coverage; passed in both v2.0 XPHM catalogs), rule 0a (passed,
306/344 = 0.890), rules 1, 1a, 2a (pair-level, already applied through Addendum D), rule 3
(calibration). If the reporting layer finds any rule fired for IMRPhenomXPHM, no on-source
quantity is computed and this is reported.

## 4. The on-source procedure (to be implemented, then frozen, then run)

Per event: the on-source segment of the PE analysis; the residual chain exactly as in the
off-source stage (common time shift ±100 ms and complex amplitude against the reference,
residual r = d − a h_ref); projection and profile scan on the cohort's common grid with
subspace A and with subspace B. Catalog: chi2 summed over the events of the cohort; q_obs(t)
= chi2_cat(t) − min_g chi2_cat(g); the 90% confidence set is the union over A and B of
{t : q_obs(t) <= c_t}, with c_t the v2.0 catalog critical values, recomputed with the frozen
seeds of the v2.0 catalog run (reproduction of the stored critical-value range is checked
first; a mismatch stops the run).

Reported per cohort:
- the 90% confidence set (and, if it is not contiguous, its envelope from the smallest to the
  largest included Lambda), in m^2;
- whether Lambda = 0 is inside it;
- the leave-one-out set without the most informative event (search: GW191129_134029;
  confirmation: the event with the smallest sigma_lin), as a diagnostic;
- the comparison with the LVK GWTC-4.0 bound [−2.4e−11, +7.4e−12] m^2 (90%).

## 5. Order

1. this registration (commit);
2. the on-source code, tested on an off-source segment as if it were on-source (the chain must
   reproduce the off-source stage's q(0) for that segment), then frozen: hash in
   `FREEZE_supplement_xphm.json`, tag `supplement-xphm-frozen`;
3. the reporting layer for IMRPhenomXPHM (rules of §3);
4. **search cohort only**; its result committed;
5. then the confirmation cohort; its result committed.
The two cohorts are reported separately; no combined set is formed.

## 6. Interpretation fixed in advance

- A set containing Lambda = 0: reported as a bound, compared with LVK.
- A set excluding Lambda = 0: reported as a candidate deviation of the class alpha ≈ 2–4
  (the f^2, f^3, f^4 shapes are not separable, §2 of the design), not as evidence for
  Lambda, with the explicit caveats that (i) the EOB companion model failed validation,
  (ii) this analysis was selected after that failure, and (iii) a confirmation cohort result
  of the same sign is required before it is described as anything more than a candidate.
- No rule, threshold, pair or construction choice is changed after any on-source number is
  seen.
