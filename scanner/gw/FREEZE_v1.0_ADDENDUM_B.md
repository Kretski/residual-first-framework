# Addendum B to the Module 1 registration (v1.0-module1)

Date: 2026-10-05. Written before any on-source estimate of Lambda.

A count check (cached geometry 27 + 20 = 47 pairs of events, not 51) led to a
regeneration of the cohort lists that applied only the fit gate and the detector
rule. It reported four "missing" search events (GW170608, GW190707_093326,
GW190728_064510, GW190924_021846), and commit 462bf63 added them back.

That correction was wrong. Commit 6f46913 records that the four events are
excluded under the existing section 3 calibration rule: none of their labels has
recalibration priors (the calibration source is the bilby priors since v0.9.8).
Commit 462bf63 is reverted; the search cohort remains 27 events.

The four events are low-mass, long-inspiral events, i.e. among the most sensitive
to the f^3 dispersion phase. Their exclusion is therefore reported explicitly as
a limitation of the search cohort.

Status: no on-source estimate had been computed at any point of this check.

# Addendum C: pre-specified diagnostic for the calibration-excluded events

Date: 2026-10-05. Written before any on-source estimate of Lambda.

Events: exactly those excluded from the search cohort solely by the section 3
calibration rule (commit 6f46913): GW170608, GW190707_093326, GW190728_064510,
GW190924_021846. Their identity follows from data availability, not from any
analysis result.

Calibration term: J built from the calibration_envelope tables of the PE files
(same spline form and node placement as the primary analysis). If a detector has
no envelope, J is omitted for that detector and the resulting interval is marked
as not including calibration uncertainty.

Analysis: the frozen v1.0 estimator (profile scan, Neyman construction, union of
subspaces A and B) on these four events as a separate catalog, both models.

Reporting: the diagnostic 90% interval is reported next to the primary search
result, irrespective of its value. It is never combined with the primary result
and cannot change it.