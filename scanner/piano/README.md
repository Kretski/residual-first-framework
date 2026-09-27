# Module 0 (piano): validation of the residual-first anomaly scanner

This folder contains the first validation of a residual-first scanner: a tool that
looks for small, reproducible frequency-dependent structure in measurement
residuals **before** any physical interpretation. The scanner is meant for later
use on gravitational-wave data, so it is first tested on systems where the answer
is known.

## What is tested

- **Baseline theory (subtracted):** ideal string, `f_n = n · f0` (f0 fitted).
- **Residuals:** `Δf_n = f_n − n · f0_fit`.
- **Fixed templates:** `Δf_n ∝ n^p`, `p ∈ {0, 2, 3, 4}` (`p = 1` is absorbed by f0).
- **Known effect:** a stiff string gives `f_n = n f0 √(1 + B n²)`, so `Δf_n ∝ n³`.
  The scanner must find structure **and** pick `n³`, not `n²` or `n⁴`.
- **Null tests:** bowed strings (arco), whose partials are exactly harmonic.

## Pre-registered rules (config hash `5bdf34b190a4b9ef`, v0.4.1)

- Per-note F-tests of each template against the baseline, Holm correction across
  templates, α = 0.01, minimum effect size 0.1 cent at n = 10.
- Partial frequencies from phase demodulation (robust to vibrato), SNR weights.
- Fixed split: even MIDI numbers = search group, odd = confirmation group.
- Calibration barrier: stop if the false-alarm count on a null set is significantly
  above α (one-sided binomial test, level 0.05).
- Exclusion rule fixed in advance: notes with fewer than 8 usable partials.

## Data

University of Iowa Musical Instrument Samples (theremin.music.uiowa.edu), freely
available without restriction:
- Piano: Steinway model B, 2001, `mf`, 86 notes (A0 and Bb0 not available).
- Violin 2012, arco ff, 90 notes (development null).
- Cello 2012, arco ff, 95 notes (final null, untouched before the final run).

`download_iowa.py` downloads everything; `results/data_manifest.csv` lists file hashes.

## History (all results reported)

| Version | Step | Result |
|---|---|---|
| v0.3 | Synthetic null and injections | Passed after adding SNR weights |
| v0.3 | Real null, violin | **Failed**: FPR 0.368 (7/19), 71/90 notes excluded. Causes found: vibrato and onset detection. Piano not analysed. |
| v0.4.1 | Fixes: onset detection, phase-demodulation frequencies, minimum effect size | Synthetic nulls with ±5–30 cent vibrato: FPR 0 |
| v0.4.1 | Violin null (development) | Passed: 1/49, binomial p = 0.389 |
| v0.4.1 | **Cello null (final barrier)** | **Passed: 2/89, binomial p = 0.224** |
| v0.4.1 | **Piano, search group** | Structure in 21/21 notes, shape **n³** (100 %) |
| v0.4.1 | **Piano, confirmation group** | Same shape, n³ in 18/19 significant notes (95 %) → **confirmed** |

Secondary interpretation (separate from detection): stiff-string coefficient
B ≈ 1.1·10⁻⁴ … 4.9·10⁻⁴ (median 1.4·10⁻⁴). Synthetic tests show B is
underestimated by about 12 % at B = 10⁻³, due to higher-order terms.

## Limitation

Only 40 of 86 piano notes (roughly B0–G4) passed the pre-registered exclusion rule.
In the upper register too few partials are recovered, most likely because
unison strings beat and destabilise the partial phases. **The conclusion applies to
the bass and middle register only.** The upper register was not tuned post hoc; the
untouched `pp` and `ff` recordings of the same piano are reserved for any future
version.

## Reproduce

```
pip install numpy scipy soundfile matplotlib
python download_iowa.py
python piano_scanner.py synth-null   --out results/synth_null
python piano_scanner.py synth-inject --out results/synth_inject
python piano_scanner.py real --null-data data/Violin --null-only --out results/null_violin
python piano_scanner.py real --data data/Piano --null-data data/Cello --plots 5 --out results/real
```

Git tags: `v0.3-piano-prereg`, `v0.3.1-piano-data`, `v0.4.1-piano-prereg`,
`v0.4.1-piano-result`.

## Next step

Module 0 with the opposite sign: shallow-water gravity waves, where the finite-depth
term gives `ω² ≈ g h k² [1 − (h²/3) k²]`, a known negative coefficient that scales
with depth.
