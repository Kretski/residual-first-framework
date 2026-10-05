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