# Addendum B to the Module 1 registration (v1.0-module1)

Date: 2026-10-05. Written before any on-source estimate of Lambda.

Finding: the committed cohort_search.txt (31 events expected by the frozen rules)
contained 27 events. Four search events that pass every registered rule were missing:
GW170608, GW190707_093326, GW190728_064510, GW190924_021846 (all pass the
reference-fit gate with both models and have at least two detectors).

How it was found: the cached geometry covered 27 + 20 = 47 events instead of 51;
the lists were regenerated from the frozen inputs (f5/f5_results.csv,
residual_check.csv) with the frozen rules (fit gate; detector rule) and compared.
The rule in module1_run.py select() is the registered one; the error was in the
committed list file only. Its origin is recorded in the git history of that file.

Correction: cohort_search.txt is replaced by the regenerated list (31 events);
cohort_confirm.txt is unchanged (20 events). The four events are low-mass,
long-inspiral events, i.e. among the most sensitive to the f^3 dispersion phase,
so their omission would have weakened the search cohort.

Status: no on-source estimate had been computed when this was found and corrected.
The geometry of the four events is built with the frozen code before stage 2.