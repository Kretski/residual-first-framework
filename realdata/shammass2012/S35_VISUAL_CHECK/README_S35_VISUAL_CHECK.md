# R0 §3.5 visual check — small-condensate points (exp_small)

**Purpose:** confirm or exclude each of the 9 digitized hollow-circle points in
`unified_v1/exports/exp_small.csv` by comparison with Figure 2 of the source PDF.
This is a visual identity check only. It is not an analysis.

**Inputs (unchanged, from `realdata/shammass2012/`):**
- `unified_v1/project/Screenshot_1.png` (sha `3c868d83…ed2e`) — the digitization image;
- `unified_v1/project/wpd.json` (sha `291f83d0…eb`) and `unified_v1/exports/exp_small.csv` (sha `a3923913…aa02`) — the 9 points;
- `source/1207.3440v2(1).pdf` (sha `ef5416e2…207c`), Figure 2, page 7 — the reference.

**Files in this folder:**
- `s35_overview_numbered.png` — `Screenshot_1.png` at 2× with the 9 points numbered 1–9 in click order (WPD points 45–53). Rings are drawn outside each marker (radius 11 px at 1×; markers are about 6 px), so no marker is covered.
- `s35_zoom_tiles.png` — one 44 × 44 px crop per point from the raw screenshot, enlarged 6× (nearest neighbour). Four tick marks point at the stored px/py with a gap, so the marker itself is untouched.
- `s35_checklist.csv` — one row per point: `no, px, py, k_um_inv, omega_over_2pi_kHz, decision, reason`. The coordinate columns are copied as text from `exp_small.csv`; `decision` and `reason` are blank.
- `make_s35_sheet.py` — the script that made the three files above (audit only).

The sheet shows only the digitized points. It contains no GPE values, no interpolation and no residuals.

## Procedure

1. **Count first, blind to the sheet.** Open the PDF, Figure 2, main panel, at 400% or more. Count the hollow circles in the main panel (not in the two insets). Write the count in the record below before opening the sheet.
2. For each point 1–9, find it in `s35_overview_numbered.png`, look at its tile in `s35_zoom_tiles.png`, and locate the same place in the PDF.
3. Apply the criterion below. Write `CONFIRMED` or `EXCLUDED` in `decision`. For `EXCLUDED`, `reason` is mandatory.
4. Do not change `px`, `py`, `k_um_inv`, `omega_over_2pi_kHz`; do not delete or add rows; do not recompute anything. If a point needs a new position, that is a new digitization (`unified_v2/`), not an edit.
5. Do not consult the legacy CSVs, the LDA table, or any distance to the GPE or LDA curves while deciding. The decision is about the marker only.
6. After all 9 rows are filled, compare the result with the legacy count of 7, for the record only.

## Criterion

**CONFIRMED** if all hold:
- the PDF main panel shows a hollow-circle marker at that location;
- the stored point lies inside that marker (within about one marker radius, ~6 px at screenshot scale, of its centre);
- it is not a filled circle, an inset element, an error bar, a curve dash, a tick, or a label;
- no other of the 9 points was placed on the same marker.

A hollow marker partly hidden behind a filled marker counts as present only if enough of its outline is visible in the PDF to identify it unambiguously as a separate hollow circle.

**EXCLUDED** otherwise. Use one of these reasons, plus a short free-text note if useful:
- `no_marker` — no hollow circle at that location in the PDF;
- `filled_marker` — the point sits on a filled (large-condensate) marker;
- `not_a_marker` — curve, error bar, tick, label or inset element;
- `duplicate_of_N` — a second placement on the same marker as point N;
- `offset` — a hollow marker exists, but the stored point is outside it;
- `ambiguous` — the PDF cannot settle it (state why).

## Record (fill in)

- Checker: ______________________   Date/time: ______________________
- PDF viewer and zoom: ______________________
- Hollow circles counted in the PDF main panel (step 1, before the sheet): ____
- CONFIRMED: ____   EXCLUDED: ____   (total must be 9)

Return the filled `s35_checklist.csv` and this record. The result goes into R0 v0.5 §3.5. Nothing is committed until the check is complete and R0 v0.5 is reviewed.


## Completed visual check record (2026-11-25)

- Checker: Dimitar Kretski.
- Reference: original arXiv:1207.3440v2 PDF, Figure 2, PDF page 7.
- Independent hollow-circle count in the main panel, before opening the numbered sheet: **9**.
- Point-by-point decisions: **9 CONFIRMED, 0 EXCLUDED**.
- Points 1–9 were reviewed sequentially against the original Figure 2 and the corresponding numbered overview/zoom tiles. The checker answered “yes” for each point.
- The decision was based on marker identity and position only. Legacy CSVs, LDA values, GPE curves/distances and residuals were not used.
- No coordinates, digitized values, or rows were changed; only `decision` and `reason` were populated in the checklist.
- This closes the visual identity/count check in R0 §3.5. It does not close the separate journal-version check, uncertainty model, synthetic power gate, interpolation specification, or formal R0 approval.
