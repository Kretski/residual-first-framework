# realdata/shammass2012 — Shammass et al. (2012), R0

Application of the residual-first framework to Figure 2 of I. Shammass, S. Rinott, A. Berkovitz, R. Schley and J. Steinhauer, arXiv:1207.3440v2 (2012); journal version Phys. Rev. Lett. 109, 195301 (2012), DOI 10.1103/PhysRevLett.109.195301.

## Status

- **Analysis: closed.** The pre-registered synthetic power/identifiability gate (spec v1.0.3) gave `NOT_TESTABLE` for both series (large condensate, N = 16; small condensate, N = 9): G1 passed, G2 failed. Statement: *not testable with available data.*
- **No experiment-minus-GPE residual has been computed**, and none is authorised under this R0. No value of ΔΛ or Λ is estimated or implied.
- **Protocol approval: not recorded.** Spec v1.0.3 was not formally approved before the gate run. This and the other deviations are recorded in `audit/R0_GATE_RUN_RECORD.md`; the code audit is in `audit/R0_GATE_CODE_AUDIT_v1.md`. The verdict does not depend on them.
- Open: comparison of the journal Figure 2 with arXiv v2 (`R0/`, §10).

## Layout

```
R0/            provenance record v0.5.3; pre-analysis spec v1.0.3
unified_v1/    primary data: one WebPlotDigitizer 5.2 session, one calibration
  exports/     exp_large.csv (16), exp_small.csv (9), gpe_large.csv (38), gpe_small.csv (28)
  project/     wpd.json (calibration and pixel coordinates), info.json
legacy/        earlier digitization and LDA table; comparison and context only
source/        citation and hash of the source paper (paper not included)
S35_VISUAL_CHECK/  record of the 9-point small-series identity check
code/          shammass_r0_v1_0_3_reviewed.py = the code that ran the gate; its tests;
               pre_gate/ = the earlier tagged version, for the audit trail only
results/       gate/ (gate result), audit/ and curve_only/ (from the run code);
               pre_gate/ (from the earlier version)
audit/         code audit and run record
design/        post-gate design calculation v2.1 (accepted artifact; does not change the verdict)
MANIFEST.sha256  SHA-256 of every file in this folder except itself
```

CSV columns: `k_um_inv` (µm⁻¹), `omega_over_2pi_kHz` (kHz), `px`, `py` (pixel coordinates in the digitization image).

## Not included, and why

The arXiv record shows no open licence, so the default arXiv non-exclusive licence applies and does not grant redistribution. Not included here: the source PDF, the screenshot of Figure 2 used for digitization, the three WebPlotDigitizer `.tar` projects that embed it, and the two images of the S35 check derived from it. The S35 README lists the files of the original check sheet: `s35_checklist.csv` and its generator `make_s35_sheet.py` are included here; the two images are kept in the local archive. Its sentence that `decision` and `reason` are blank describes the sheet before the check; the completed record is at the end of that README. Hashes of all excluded files are in `R0/SHAMMASS_R0_DATA_PROVENANCE_v0.5.3.md` (§11). This is a conservative choice, not a legal opinion on the digitized values.

## Reproducing the non-gated stages

From the repository root (Python ≥ 3.10, numpy, scipy):

```
python realdata/shammass2012/code/shammass_r0_v1_0_3_reviewed.py . --stage audit --out <dir>
python realdata/shammass2012/code/shammass_r0_v1_0_3_reviewed.py . --stage curve-only --out <dir>
```

The audit stage checks the four input hashes and reproduces every exported value from its pixel coordinates. The code finds `wpd.json` by searching below the root, so only one `wpd.json` may exist in the repository. It looks for the spec beside itself or in `docs/` and will not find it in this layout; this affects only the record in the gated real-data stage, which is not run. The code is kept byte-identical (SHA-256 `007d56ad92cca87959cab349231c751c2532f0c921c64f4974a1ae40036c06d9`) because that hash links it to the gate result.

## History outside this repository

The implementation was developed in a separate local repository. Its commits: `b70f1188b7eb7c24b6fc617b6a4461ff248ea8ac` (tag `r0-v1.0.3-pre-gate`), `a8107b4` (tag `r0-v1.0.3-gate-result`), `065dcc0`. They are recorded here by hash only.
