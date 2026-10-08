\# Diagnostic report: rule 0a failure of Module 1 v1.0



Date: 2026-10-08. Pre-specified in DIAG\_0a\_SPEC.md (d96415a, clarified 534d243);

code diag\_0a\_m1.py and diag\_0a\_gate2.py (b75b5da), committed before running.

This report cannot change the v1.0 result (tag v1.0-stopped-0a).



Baseline: the v1.0 off-source flags were reproduced exactly for all 12 pairs

(0 mismatches).



Measurement 1: posterior-sample time offsets relative to the reference do not separate

the failing SEOBNRv4PHM events from the controls (large in failing GW190706\_222641 and

control GW200112\_155838; compact in GW150914, GW170814, GW200129\_065458). SEOBNRv4PHM

samples scatter by 2-18 ms in failing and control events alike; IMRPhenomXPHM samples

by 0.1-0.6 ms.



Gate 2 (set-C coverage, pooled, same segments, injections, belts and seeds):

&#x20; variant              failing EOB   controls EOB   failing events, XPHM

&#x20; v1.0                 16/32         27/32          30/32

&#x20; (a) medoid           20/32         28/32          29/32

&#x20; (b) time/phase       23/32         30/32          26/32

Reading (threshold 24/32, controls may not drop by more than 3/32):

&#x20; (a) NOT SUPPORTED;  (b) NOT SUPPORTED.



Descriptive observations (no effect on the reading): the stored reference lies 13-40

whitened units from the posterior medoid in the failing events, but re-centring does

not restore coverage consistently (GW190706 4->7/8, GW200129 4->2/8); alignment

improves three of four failing events but lowers the IMRPhenomXPHM cross-check from

30/32 to 26/32; both variants change sigma\_lin, i.e. the sensitivity of the

construction. The same four events pass with IMRPhenomXPHM (30/32): the failure is

specific to the SEOBNRv4PHM labels.



Conclusion: no mechanism identified. Neither a mis-centred reference nor

parameter-dependent time offsets, taken separately, explain the rule 0a failure.

Any v2.0 is a new, separately registered architecture with full validation before

any on-source estimate.

