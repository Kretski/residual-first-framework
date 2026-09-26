# PTA sensitivity to the Λ_GW dispersion (Paper 3) — closed, negative

**Conclusion.** For the currently constrained range of Λ_GW and its f³ dispersion,
PTAs have no practically measurable sensitivity to Λ_GW. The EPTA/NANOGrav HD results
(DR2full p = 0.0135, DR2new p = 0.0065) concern spatial residual correlations only and
must not be used as evidence for or against the Λ model.

## Assumptions
- ΔΨ(f) = −4π³ Λ_GW K(z) f³ / c³ with [Λ_GW] = m³/s, hence [K] = s; K ≈ D/c assumed.
  Cosmological factors in K(z) change results by factors of order a few.
- |Λ_GW| = 7.24×10⁻³ m³/s (upper edge of the published GW range).
- Equivalent dispersion ω² = c²k² + αk⁴, α = Λ_GW c, giving |v_ph/c − 1| ≈ 2π²|Λ_GW|f²/c³.

## Numbers (`pta_lambda_sensitivity.py`)

| f | \|ΔΨ\| at 1 Gpc | \|v_ph/c − 1\| | HD-curve effect (L = 1 kpc) |
|---|---|---|---|
| 1 nHz | 3×10⁻³⁶ rad | 5×10⁻⁴⁵ | 3×10⁻⁴² |
| 3 nHz | 9×10⁻³⁵ rad | 5×10⁻⁴⁴ | 9×10⁻⁴¹ |
| 10 nHz | 3×10⁻³³ rad | 5×10⁻⁴³ | 3×10⁻³⁹ |
| 30 nHz | 9×10⁻³² rad | 5×10⁻⁴² | 9×10⁻³⁸ |

For a stochastic background the source phase ΔΨ is unobservable anyway (random source
phases); the HD curve can change only through v_ph, where the effect is < 10⁻³⁷.

Reference points from the same formula (order of magnitude only — NOT statements about
what GWTC or any detector can measure, which depends on waveform degeneracies, noise,
K(z) and the full likelihood): 10 mHz (LISA band), 1 Gpc: ~3×10⁻¹⁵ rad;
100 Hz, 40 Mpc: ~1×10⁻⁴ rad; 1 kHz, 1 Gpc: ~3 rad.

## Implication
Searches for Λ_GW should focus on high-frequency GW observations (upper LIGO band,
Einstein Telescope, Cosmic Explorer). Lower-frequency detectors, including LISA, are
suppressed by (f/f_LIGO)³.
