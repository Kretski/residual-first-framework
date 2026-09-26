# Residual-level Hellings–Downs diagnostics for pulsar timing arrays

A simple, fully reproducible residual-level analysis of NANOGrav 15-yr and EPTA DR2
timing residuals, with a **pre-registered replication on EPTA** and a documented
record of which statistics are — and are not — calibrated on real PTA data.

📄 **Full write-up:** [`METHODS_NOTE.md`](METHODS_NOTE.md) ·
📋 **Pre-registration, amendments and all results:** [`PREREGISTRATION_EPTA.md`](PREREGISTRATION_EPTA.md)

## Main results

| Data | Test | p | Status |
|---|---|---|---|
| EPTA DR2full | HD ⟂ monopole + dipole, sky-scramble null | **0.0135** | pre-registered **primary** |
| EPTA DR2new | same test | **0.0065** | pre-registered **secondary** |
| NANOGrav 15-yr | same test | 0.046 | exploratory (test chosen after seeing these data) |

Both EPTA calibration gates passed (false-positive rate under H₀: 0.045 and 0.040).
In EPTA geometry the test is robust to injected monopole and dipole processes.

**Calibration failures (equally important):**

- An omnibus cross-correlation statistic, calibrated on idealised simulations (0.065),
  has a false-positive rate of **0.69 (NANOGrav) to 1.00 (EPTA)** on real sampling.
  Cause: variance that changes at the same dates for many pulsars, plus red-noise-dominated
  series (see `diagnostics_shared_structure.py`).
- Temporal-shift nulls fail for monopole and angular tests (0.10–0.15).
- Only sky-scramble tests stay calibrated in all three data sets.

## What this does *not* claim

This recovers the HD correlations already reported by the PTA collaborations with a much
simpler method. It is **not a new detection** and has **no bearing on modified-gravity
models** (see [`PTA_LAMBDA_SENSITIVITY.md`](PTA_LAMBDA_SENSITIVITY.md)).

## Reproduce

Run from a Linux filesystem (WSL is fine); PINT needs 80-bit long double.

```bash
pip install pint-pulsar numpy scipy
python fetch_nanograv15.py                     # NANOGrav 15-yr v2.1.0 (Zenodo, MD5-checked)
python common_residual_search_v2.py --real
python common_residual_search_v2.py --validate --real-gaps --trials 200
python diagnostics_shared_structure.py         # why the omnibus test fails
python fetch_epta_dr2.py                       # EPTA DR2 (Zenodo, MD5-checked)
python run_epta.py gate   [--variant DR2new]   # calibration gate
python run_epta.py real   [--variant DR2new]   # pre-registered real-data run
python run_epta.py dipole [--variant DR2new]   # nuisance-leakage check
```

## Key files

| File | Role |
|---|---|
| `common_residual_search_v2.py` | frozen analysis (tag `v2.1-epta-prereg`, SHA-256 `089a94e0…`) |
| `run_epta.py` | EPTA runner: PINT loading with tempo2 `.tim` semantics, gate, real run, checks |
| `nanograv_loader.py`, `fetch_*.py` | data download and reading |
| `diagnostics_shared_structure.py` | synthetic scenarios isolating the omnibus failure |
| `*_results.json`, `epta_gate_*.json`, `epta_dipole_check_*.json`, `*_log.txt` | provenance records |

Git history is part of the record: every decision was committed before the run it affects.

## Data

- NANOGrav 15-Year Data Set v2.1.0, [doi:10.5281/zenodo.16051178](https://doi.org/10.5281/zenodo.16051178)
- EPTA DR2, [doi:10.5281/zenodo.8300645](https://doi.org/10.5281/zenodo.8300645)

## Use of AI

Analysis code and text were developed with the assistance of an AI model (Claude, Anthropic).
All runs, decisions and pre-registration commits were made and checked by the author.

## Author

Dimitar Kretski — Institute of Metal Science, Equipment and Technologies "Acad. A. Balevski",
Center for Hydro- and Aerodynamics, Bulgarian Academy of Sciences, Varna.
