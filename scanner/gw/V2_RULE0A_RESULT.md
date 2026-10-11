# Module 1 v2.0: rule 0a fired — v2.0 stops before any on-source estimate

Date: 2026-10-11. The reading was fixed in V2_RULE0A_INTERIM.md (68df69d) before the
confirmation result was known.

Rule 0a (per model, both cohorts pooled, new segments, set-E injections):
  IMRPhenomXPHM: 306/344 = 0.890 (99% lower end ~0.858) -> pass
  EOB: search SEOBNRv4PHM 167/200 = 0.835, confirmation SEOBNRv5PHM 132/152 = 0.868,
       pooled 299/352 = 0.849 (99% lower end ~0.859) -> STOP
  (exact binomial at p = 0.90: P(X <= 299 | n = 352) ~ 1e-3)
Reference-waveform injections (diagnostic): 183/200 = 0.915 and 133/152 = 0.875.

Comparison with v1.0 (old segments, set C): search SEOBNRv4PHM 0.775 -> 0.835;
confirmation SEOBNRv5PHM 0.906 -> 0.868; EOB pooled 0.833 -> 0.849; XPHM 0.889 -> 0.890.

As registered, v2.0 stops before any on-source estimate. No on-source estimate has been
computed in v1.0 or v2.0. Including the posterior/reference mismatch in the belts improved
the search-cohort SEOBNRv4PHM coverage but did not restore it, and the EOB pooled coverage
remains below the threshold: the mismatch term explains part of the v1.0 failure, not all of
it. The v2.0 result recorded here stands.