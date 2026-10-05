# Addendum D to the Module 1 registration (v1.0-module1)

Date: 2026-10-06. Written before any on-source estimate of Lambda.

The reporting layer (module1_report.py, FREEZE_v1.1_report.json) applied the
section 7 rules to the official stage-2 outputs. As registered, the A/B subspace
rule and the GR-leakage KS rule act on event-model PAIRS (section 7: "moves the
event-model pair from the primary test to the separately reported set"); only the
reference-fit gate acts on events. Result:

  search  IMRPhenomXPHM: 24 events (pairs out: GW190512_180714, GW190620_030421,
          GW190915_235702)
  search  SEOBNRv4PHM:   25 events (pairs out: GW170823, GW191109_010717)
  confirmation: 20 and 20 events, no pair out.

Primary catalogs are computed per model on these lists (primary_<cohort>_<model>.txt),
as registered. Consequence for the two-model comparison (section 10): the two search
catalogs contain different events. The comparison is therefore ALSO reported on the
intersection (22 search events, intersection_search.txt; confirmation unchanged, 20),
as a diagnostic. This adds a reporting item and changes no registered rule.

Queue note: four catalog commands queued in the build terminal used the full cohort
lists; they were disarmed by moving cohort_*.txt to hold/ so that they stop on a
missing argument without computing anything. Their failure lines in
official_catalog.log are expected.