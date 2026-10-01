# Notes on event list v1

- "missing mass_2_source" (43 GWTC-4.0 candidates): these candidates have no public
  parameter-estimation results in GWTC-4.0 (GWOSC provides PE values only for the
  events of Table III of the catalog paper; PE was run only for candidates with
  FAR < 1/yr). The method requires a public PE reference point, so they cannot be
  analysed. Two of them have search SNR >= 12 (GW230728_083628: 13.1,
  GW240105_151143: 25.9); this SNR is from the search, not from PE.
- Three selected events have network_matched_filter_snr exactly 12.0 (GW170818,
  GW190924_021846, GW231020_142947). Catalog values are rounded to 0.1, so the
  inclusive threshold makes the selection sensitive to rounding at 12.0; recorded,
  not changed.
- O1 and O2 events are taken from GWTC-2.1-confident (reanalysed with
  IMRPhenomXPHM and SEOBNRv4PHM), not from the original GWTC-1 results.