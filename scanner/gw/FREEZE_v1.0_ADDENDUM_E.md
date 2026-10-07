\# Addendum E: rule 0a fired — Module 1 v1.0 stops before any on-source estimate



Date: 2026-10-07.



Rule 0a (catalog level, per model, set-C injections, both cohorts pooled):

&#x20; IMRPhenomXPHM: 313/352 = 0.889 (99% lower end \~0.859) -> pass

&#x20; EOB (SEOBNRv4PHM search + SEOBNRv5PHM confirmation): 300/360 = 0.833 -> STOP

&#x20; (search SEOBNRv4PHM alone 155/200 = 0.775; confirmation SEOBNRv5PHM 145/160 = 0.906)

Reference-waveform injections (diagnostic): search SEOBNRv4PHM 170/200 = 0.850.

The search/confirmation breakdown is diagnostic and is not used to bypass rule 0a.



As registered, the analysis stops before any on-source estimate. No on-source estimate

has been computed. The v1.0 result recorded here stands.



Per-event diagnostic (search SEOBNRv4PHM, 8 segments each): the four lowest

(4/8: GW150914, GW170814, GW190706\_222641, GW200129\_065458) are the events whose

LALInference reference needed the largest time shifts in the residual check

(-15.7, -15.1, -39.0, -23.0 ms); events needing < 2 ms reach 6-8/8.



Hypotheses (not yet tested):

&#x20; (a) the LALInference reference point may lie materially away from the likelihood

&#x20;     maximum, so that the linearised free subspace built around it is mis-centred;

&#x20; (b) the time convention of the LALInference samples may depend on the source

&#x20;     parameters, so that posterior samples differ from the reference by

&#x20;     parameter-dependent time offsets that the principal-component subspace cannot

&#x20;     represent linearly.



Any diagnosis uses only the existing off-source segments, is specified before it is

run, and is reported whatever its outcome. Any fix is a new code version with its own

registration, validated by all stages again, before any on-source estimate.

