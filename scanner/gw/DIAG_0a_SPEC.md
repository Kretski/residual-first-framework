# Pre-specified diagnostic of the rule 0a failure (Module 1 v1.0)

Date: 2026-10-07. Written before the diagnostic is run. It cannot change the v1.0
result (tag v1.0-stopped-0a); it only tests the two hypotheses of Addendum E.

Data: only the existing off-source segments and set-C injections of
primary_offsource_search_eob.csv (same segments, same injected samples, same seeds).
No on-source data are used.

Events, fixed by rule:
  failing: the four search SEOBNRv4PHM events with the lowest per-event coverage (4/8):
           GW150914, GW170814, GW190706_222641, GW200129_065458;
  controls: the four events of primary_search_eob.txt with the smallest |dt_ms| in
           residual_check.csv for their SEOBNRv4PHM label.
  Cross-model check: the IMRPhenomXPHM labels of the four failing events.

Measurement 1 (no change to the method): for every label, the time shift that best
aligns each of 200 set-B posterior samples with the reference waveform (network
overlap, ±100 ms); reported: median and standard deviation, and their relation to the
per-event coverage.

Variants, each on the same segments and injections:
  v1.0  baseline, recomputed; must reproduce the stored 'covered' flags exactly,
        otherwise the diagnostic stops;
  (a)   reference = the posterior medoid (the sample with the smallest mean whitened
        distance to the other samples) instead of the stored LALInference sample;
  (b)   posterior samples aligned in time (and phase) to the reference before the
        principal-component subspace is built.

Metric: pooled set-C coverage of the four failing events (baseline 16/32) and of the
four controls, per variant.

Reading, fixed in advance: a hypothesis is supported if its variant raises the pooled
set-C coverage of the four failing events to at least 24/32 (24/32 is the lower end of
the 99% interval of Binomial(32, 0.90): P(X <= 23) = 0.0033, P(X <= 24) = 0.0117), while
the pooled coverage of the four controls does not fall by more than 3/32 from its
recomputed baseline. Otherwise it is not supported. Both variants are reported
whatever the outcome. Neither variant is adopted on the basis of this diagnostic
alone; a fix is a new version with its own registration and full validation.
Clarification (2026-10-07, before Gate 2 is run): the medoid of variant (a) is taken among the same 200 set-B samples used in Measurement 1(ii), with the raw whitened distance (no alignment), as written above. The reading threshold remains 24/32.
