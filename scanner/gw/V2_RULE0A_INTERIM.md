# Module 1 v2.0: interim record before the confirmation off-source result

Date: 2026-10-11. Written while the SEOBNRv5PHM confirmation off-source stage is running,
before its result is known.

Search cohort, SEOBNRv4PHM, new segments and set E (v2.0):
  E coverage 167/200 = 0.835 (v1.0, old segments and set C: 155/200 = 0.775)
  reference-waveform injections (diagnostic): 183/200 = 0.915 (v1.0: 170/200 = 0.850)
IMRPhenomXPHM, both cohorts pooled: 306/344 = 0.890 -> passes rule 0a.

Reading fixed now, before the confirmation result:
- Rule 0a is applied as registered: per model, both cohorts pooled, lower end of the 99%
  binomial interval around 0.90. The per-cohort STOP line printed by the driver is a
  diagnostic.
- Whatever the pooled verdict, the search-cohort SEOBNRv4PHM coverage of 0.835 is reported
  explicitly and separately. If the pooled rule passes only because of the confirmation
  cohort, v2.0 has only partly repaired the v1.0 failure, and any on-source search result with
  SEOBNRv4PHM is reported with this limitation stated next to it.
- No rule, threshold or pairing is changed in response to either result.