# Residual-First Discovery Framework
## Design specification for independent morphology audits

**Status:** Framework specification for Lambda-model residual-seeking pipeline  
**Authors:** Dimitar Kretski  
**Date:** September 2026  
**Scope:** GW, pulsar timing, cosmological residuals — without Lambda assumptions  

---

## 1. Problem Statement

**Before (failed):** $\Lambda$ → seek evidence in all systems simultaneously  
**Now (honest):** DATA → baseline → residual → blind morphology → (replication) → *only then* Lambda  

Three sectors, three independent null hypotheses, one common statistical framework.

---

## 2. Core Null Hypotheses

### 2.1 Gravitational Waves (Primary pilot)

**System:** GWTC-4.0 events (LVK public release)

**Residual definition:**
$$
r_i(t) = d_i(t) - h_{\rm GR,i}(t)
$$

where $d_i$ is strain from detector $i$ and $h_{\rm GR,i}$ is GR prediction (IMR waveform from LIGO's own pipeline).

**Null hypothesis ($H_0$):**
> The residual $r(t)$ is consistent with stationary, Gaussian detector noise, independent of:
> - event properties (mass, spin, redshift)
> - detector orientation
> - gravitational-wave frequency content
> - propagation distance

**Rejection criterion ($H_1$):**
> Existence of structure in $r(t)$ that:
> - is reproducible across ≥2 independent events
> - is not explained by known detector systematics
> - exhibits coherence between separated detectors (H1↔L1)
> - has frequency or time-frequency structure not consistent with Gaussian noise

**Key:** Rejection of $H_0$ does **not** imply Lambda. It implies: "something is wrong with the GR baseline or detector model." Lambda is tested later, *if at all*.

---

### 2.2 Pulsar Timing (Independent cross-check)

**System:** NANOGrav DR15 (public posterior chains available)

**Residual definition:**
$$
R_p(t) = t_{\rm obs,p}(t) - t_{\rm model,p}(t)
$$

where model includes standard timing corrections (dispersion measure, parallax, proper motion, binary orbit if applicable).

**Null hypothesis:**
> Timing residuals are consistent with:
> - white noise (flat spectral density)
> - red noise from spin-down evolution (red noise with known index)
> - independent between pulsars
> - no correlation with observing frequency $\nu_{\rm rf}$ (achromatic)

**Rejection criterion:**
> - Common residual across multiple pulsars (correlated noise)
> - Frequency-dependent structure inconsistent with known astrophysics
> - Systematic deviation from red+white noise model
> - Reproducible in independent timing arrays

**Key:** This is a *fundamentally different* physical scale (years vs ms/s) — if structure emerges independently in both GW and pulsar data, that is *interesting* but does not automatically connect them to Lambda.

---

### 2.3 Cosmological Residuals (CMB/BAO — lower priority)

**System:** Planck DR4 + DESI DR1

**Residual definition (CMB):**
$$
r_\ell = C_\ell^{\rm obs} - C_\ell^{\Lambda \rm CDM, best-fit}
$$

**Null hypothesis:**
> Residuals are consistent with:
> - cosmic variance (expected for finite sky coverage)
> - measurement uncertainty
> - foreground/calibration systematics
> - *independent* of assumed base cosmology parameters

**Rejection criterion:**
> - Systematic trend across multiple $\ell$ ranges
> - Correlated with known systematics (dust, calibration) but with opposite sign
> - Reproducible in independent CMB datasets
> - Not explained by parameter degeneracies

**Note:** Cosmological residuals are *last* because instrumental systematics dominate. Start with GW and pulsar timing where instrumental models are better understood.

---

## 3. Statistical Framework

### 3.1 Baseline Subtraction and Noise Characterization

**Stage 0: Blind baseline**

For each event/system:

1. Compute $r(t) = d(t) - \text{model}(t)$
2. Do not inspect $r(t)$ yet
3. Compute power spectral density (PSD) of $r(t)$:
   $$S_r(f) = \frac{1}{T}\left|\tilde{r}(f)\right|^2$$
4. Fit PSD to standard noise model:
   - GW: $S_r(f) = S_0 f^{\alpha}$ (assume $\alpha \approx -3$ to $+3$)
   - Pulsar: $S_r(t) = S_0 + S_r(t) \sim f^{-\beta}$ with $\beta \in [0, 3]$ for red noise

### 3.2 Residual Characterization (Blind, no Lambda)

Test 1: **Spectral shape anomaly**
$$
\Delta S(f) = S_{\rm obs}(f) - S_{\rm noise\ model}(f)
$$

Reject $H_0$ if $|\Delta S(f)|$ exceeds $k\sigma$ in one or more frequency bands, with look-elsewhere correction (see 3.5).

Test 2: **Time-frequency concentration**

Compute Morlet wavelet transform:
$$
W(t, f) = \int r(t') \psi^*\left(\frac{t'-t}{a(f)}\right) dt'
$$

Test null hypothesis: $|W(t,f)|^2$ is consistent with Gaussian noise → exponential distribution in $\log$ space.

Reject if coherent time-frequency tracks emerge (not Gaussian scatter).

Test 3: **Autocorrelation anomaly**
$$
\rho_r(\tau) = \frac{\langle r(t) r(t+\tau) \rangle}{\sigma_r^2}
$$

Null: $\rho_r(\tau) = 0$ for $|\tau| > \tau_{\rm coherence}$ (instrument-dependent).  
Reject if excess correlation at unexpected lags.

Test 4: **Cross-detector coherence** (GW only)

For H1, L1 separated by $\Delta L \approx 3000$ km:
$$
\gamma^2(f) = \frac{|C_{H1,L1}(f)|^2}{S_{H1}(f) S_{L1}(f)}
$$

where $C_{H1,L1}$ is cross-spectrum.

Null: $\gamma^2(f) < \gamma^2_{\rm noise}$ (expected from independent noise).  
Reject if $\gamma^2(f)$ exceeds noise floor in frequency bands where $H_0$ would predict coherence.

**Critical:** Tests 1–4 are run *before* inspecting results qualitatively. Stopping rule is statistical, not subjective.

### 3.3 Anomaly Gating

**For each test $T_i$ (spectral, time-frequency, autocorr, coherence):**

Compute p-value:
$$
p_i = P(T_i \geq T_i^{\rm obs} | H_0)
$$

**Do not use raw $p_i$.** Apply look-elsewhere correction (§3.5).

### 3.4 Multiple Testing and Look-Elsewhere Effect

**Problem:** Testing $N_f$ frequency bins, $N_t$ time windows, $N_s$ statistics → multiple comparisons inflate false-positive rate.

**Solution:** Benjamini–Hochberg (BH) step-down procedure.

1. Sort all $p_i$ from all tests across all parameters: $p_1 \leq p_2 \leq \ldots \leq p_M$
2. Find largest $i$ such that $p_i \leq \frac{i}{M} \alpha_{\rm global}$
3. Reject $H_0$ only for $T_1, \ldots, T_i$

**Global $\alpha$:** Recommend $\alpha_{\rm global} = 0.05$ per independent dataset.

**Effective $N$ for look-elsewhere:** Not the number of raw bins, but the number of independent resolution elements:

For GW:
- Frequency resolution: $\Delta f \approx 1 / T_{\rm observation}$
- Time resolution: $\Delta t \approx 1 / \Delta f$
- Number of independent bins: $N_f \approx T \Delta f_{\rm max}$ (not $> 10^5$ even for full catalog)

For pulsar timing:
- Time span: $\sim 15$ years
- Typical TOA cadence: $\sim 2$ weeks → $N_t \approx 400$

### 3.5 Replication Criterion (Stage 3)

**Single-event anomaly = not publishable.**

**Replication rule:**

$$
\boxed{\text{Anomaly in event 1} \implies \text{Same anomaly in} \geq 2 \text{ independent events} \implies \text{Not consistent with} H_0}
$$

**"Same anomaly" definition:**

- Same frequency band (within $\Delta f < 2\Delta f_{\rm resolution}$)
- Same morphology (coherent track, peak, or spectral feature)
- Same statistical significance after independent BH correction
- Independent detector pair or pulsar system

**For GW:** At least one other GWTC-4.0 event (H1/L1 pair) must show matching structure.  
**For pulsar:** At least 3 pulsars in independent arrays (NANOGrav + IPTA baseline) must show common residual.

If no replication → $H_0$ not rejected, move to null-result documentation.

---

## 4. Architecture and Stopping Rules

### 4.1 Pipeline Stages

```
┌─────────────────────────────────────────┐
│ STAGE 0: BLIND BASELINE                 │
├─────────────────────────────────────────┤
│ r(t) = d(t) - model(t)                  │
│ No inspection, no plotting yet          │
└──────────────┬──────────────────────────┘
               │
┌──────────────▼──────────────────────────┐
│ STAGE 1: RESIDUAL CHARACTERIZATION      │
├─────────────────────────────────────────┤
│ • PSD fitting                           │
│ • Spectral anomaly test                 │
│ • Time-frequency concentration test     │
│ • Autocorrelation test                  │
│ • (GW: coherence test)                  │
│                                         │
│ Compute p-values (blind)                │
└──────────────┬──────────────────────────┘
               │
┌──────────────▼──────────────────────────┐
│ STAGE 2: ANOMALY GATING                 │
├─────────────────────────────────────────┤
│ Benjamini–Hochberg multiple testing     │
│ correction across all tests & bins      │
│                                         │
│ Decision: Reject H_0?                   │
└──────────────┬──────────────────────────┘
               │
        ┌──────┴──────┐
        │             │
    NO  │             │  YES
        │             │
   ┌────▼──┐      ┌───▼─────────────────────┐
   │ NULL  │      │ STAGE 3: REPLICATION    │
   │ result│      ├─────────────────────────┤
   └───────┘      │ Independent 2nd event   │
                  │ (or 2nd system/pulsar)  │
                  │                         │
                  │ Same anomaly?           │
                  └────┬──────┬─────────────┘
                       │      │
                   NO  │      │  YES
                       │      │
                  ┌────▼──┐ ┌─┴──────────────┐
                  │ NULL  │ │ STAGE 4:       │
                  │ result│ │ LAMBDA TEST    │
                  └───────┘ │ (if at all)    │
                            │                │
                            │ Compatible    │
                            │ with Lambda   │
                            │ prediction?   │
                            │                │
                            │ YES/NO/MAYBE  │
                            └────────────────┘
```

**Stopping rule:** 
- If anomaly not replicated: STOP. Report null result.
- If anomaly replicated: proceed to Stage 4 *only if* there is a testable Lambda prediction.
- If no Lambda prediction exists: STOP. Report "Anomaly in residuals, not explained by GR baseline."

### 4.2 Code Structure

```
residual_first/
│
├── common/
│   ├── __init__.py
│   ├── baseline_subtraction.py     # r(t) = d(t) - model(t)
│   ├── psd_fitting.py              # Fit S_r(f) to noise model
│   ├── spectral_anomaly.py         # Test 1
│   ├── timefreq_anomaly.py         # Test 2 (Morlet wavelet)
│   ├── autocorr_anomaly.py         # Test 3
│   ├── coherence_anomaly.py        # Test 4 (GW only)
│   ├── multiple_testing.py         # Benjamini–Hochberg
│   ├── replication_check.py        # Match anomalies across events
│   └── report_generator.py         # Null-hypothesis-safe reporting
│
├── gw/
│   ├── __init__.py
│   ├── gwtc4_residual_loader.py    # Load LVK GWTC-4 data
│   ├── gw_morphology_audit.py      # Run stages 0–3 on single event
│   ├── gw_catalog_scan.py          # Apply to all GWTC-4 events
│   └── gw_replication_framework.py # Implement replication rule
│
├── pulsar/
│   ├── __init__.py
│   ├── nanograv_loader.py          # Load NANOGrav DR15
│   ├── timing_residual_audit.py    # Stages 0–3 on pulsars
│   └── common_residual_search.py   # Test for correlated residuals
│
├── cosmology/
│   ├── cmb_residual_audit.py       # (Lower priority)
│   └── bao_residual_audit.py       # (Lower priority)
│
└── tests/
    ├── test_synthetic_anomalies.py # Generate false positives,
                                    # verify BH catches them
    ├── test_replication_logic.py   # Check replication rule
    └── test_look_elsewhere.py      # Verify multiple testing correction
```

---

## 5. Success Criteria (What Does "Finding Something" Mean?)

**Not enough:**
- Single event shows deviation from GR baseline
- Spectral feature appears in one frequency band
- Visual inspection suggests "interesting structure"

**Necessary for publication:**
1. ✅ Anomaly identified in event A (Stage 1–2)
2. ✅ Same morphology in event B, independent detector pair (Stage 3)
3. ✅ Does not appear in ≥3 other events (not universal)
4. ✅ Not explained by known systematics (instrumental, calibration)
5. ✅ Survives BH multiple-testing correction
6. ✅ Clear null hypothesis and stopping rule documented *before* opening data

**Claim credibility:**
- "GW residuals show repeatable structure" → publishable as methodology paper
- "Structure is inconsistent with Gaussian noise" → publishable as discovery
- "Structure matches Lambda prediction" → publishable only if Lambda model was *not* used to design the scan

**Not publishable (common trap):**
- "We looked for Lambda and found something; let's call it a Lambda signal"
  - This is exactly the look-elsewhere effect that BH corrects for.

---

## 6. Timeline and Priorities

### Phase 1 (Week 1–2): Framework Setup
- [ ] Finalize null hypotheses (document in YAML)
- [ ] Implement `common/` module (baseline, PSD, BH correction)
- [ ] Write synthetic test suite (verify false positive rate)
- [ ] Design documentation template (captures hypothesis before data)

### Phase 2 (Week 3–4): GW Pilot
- [ ] Load GWTC-4.0 public residual data (LVK RES.tar.gz)
- [ ] Run stages 0–1 on 5–10 events
- [ ] Characterize noise floors (establish baseline anomaly rate)
- [ ] Identify any candidates for replication check

### Phase 3 (Week 5–6): Replication and Interpretation
- [ ] Apply replication rule to candidates
- [ ] Generate null-result report (if no anomalies survive BH)
- [ ] Or document anomaly class (if replication succeeds)
- [ ] Do *not* mention Lambda at this stage

### Phase 4 (Week 7–8): Pulsar Timing (if GW shows signal)
- [ ] Load NANOGrav DR15
- [ ] Apply same stages 0–3
- [ ] Check for independent anomaly class

### Phase 5 (Conditional): Lambda Testing
- [ ] Only if anomaly class is robust and replicated
- [ ] And only if there is a concrete Lambda prediction to test
- [ ] Compare observed morphology to Lambda propagation equation
- [ ] Use independent parameter (not the one used for discovery)

---

## 7. Reporting and Reproducibility

### 7.1 Required Documentation

**Before opening any residual data:**

```yaml
# residual_audit_hypothesis.yaml
experiment: "GWTC-4.0 gravitational waves"
events: [GW150914, GW151226, ...]
null_hypothesis: "Residuals are stationary Gaussian noise"

baseline_model: "IMR waveform from LIGO pipeline"
residual_definition: "r(t) = d(t) - h_GR(t)"

tests:
  - name: "Spectral anomaly"
    method: "PSD deviation from power-law fit"
    correction: "Benjamini–Hochberg, alpha=0.05"
  
  - name: "Time-frequency concentration"
    method: "Morlet wavelet, excess coherence"
    correction: "BH, accounting for N_t * N_f bins"
  
  - name: "Coherence between H1 and L1"
    method: "Cross-spectrum correlation coefficient"
    correction: "BH"

replication_rule: "Same anomaly in ≥2 independent events"
stopping_rule: "If no replication → report null; if replication → Stage 4"

lambda_hypothesis: "NOT YET FORMULATED" (filled in only if Stage 3 succeeds)
```

**After Stage 3 (if anomaly found):**

```yaml
anomaly_class:
  morphology: "Coherent spectral peak at f = X Hz"
  events_with_anomaly: [GW150914, GW151226]
  p_value_after_BH: 0.0023
  not_explained_by: [calibration, detector glitch, template mismatch]

lambda_test_hypothesis: "IF APPLICABLE"
  prediction: "r(t) should satisfy Λ dispersion equation with parameter Λ = X"
  testable_difference: "Phase coherence over frequency range Y"
  expected_effect_size: "ΔΦ = Z radians"
```

### 7.2 Reproducibility

- [ ] All data loaded from public LVK/GWTC servers (not private)
- [ ] Code version control (GitHub)
- [ ] Synthetic false-positive tests (show BH catches them)
- [ ] Null hypothesis fixed before analysis
- [ ] Multiple-testing correction applied algorithmically (not by eye)
- [ ] Results tables include raw and corrected p-values

---

## 8. Known Risks and Mitigations

| Risk                                    | Mitigation                                                    |
| --------------------------------------- | ------------------------------------------------------------- |
| Visual inspection biases analyst        | BH correction is algorithmic, not subjective                  |
| False positives from look-elsewhere     | Replication rule (≥2 events) + BH correction                  |
| Confusing "structure" with "anomaly"    | Clear null hypothesis; tests are predefined                   |
| Jumping to Lambda prematurely           | Stage 4 is skipped if Stage 3 fails; stopping rule is explicit |
| Data-dependent hypothesis in Lambda     | If anomaly is found first, Lambda prediction must be independent |

---

## 9. Example: What a Null Result Looks Like

```
GWTC-4.0 Residual Morphology Audit — Null Result

Null hypothesis: GW residuals are consistent with GR baseline + 
                 stationary Gaussian detector noise.

Sample: GWTC-4.0 events with public residual data (N=42 events)

Stages 0–1 (Baseline and characterization):
- All 42 events show PSD consistent with instrument noise model
- No spectral anomalies exceed BH threshold (α=0.05)
- No time-frequency concentration in any event
- Cross-detector coherence consistent with expected noise floor

Stage 2 (Multiple testing correction):
- 12 frequency ranges flagged in individual events
- After Benjamini–Hochberg correction: 0 survive (all p > 0.05)

Stage 3 (Replication):
- No anomalies to replicate

Conclusion:
H_0 is not rejected. Residuals are consistent with GR + known noise.

No evidence for:
- Λ-type dispersion
- Systematic mismatch in GR waveforms
- Novel post-merger physics
- Unexpected frequency-dependent effects

This is scientifically valuable: it tightens the constraint on what 
the GR baseline does NOT include.
```

---

## 10. Next Concrete Step

**This week:** Implement `common/multiple_testing.py` and `synthetic_test_suite.py`.

Test the algorithm on *fake* data where you know the answer:
- Generate $r(t)$ from pure Gaussian noise → BH should not reject $H_0$
- Generate $r(t)$ with injected 5σ spectral peak → BH should reject at $p < 0.05$
- Verify that false positive rate is actually $\alpha$ under $H_0$

Once that works, GWTC-4 data can be opened.

---

**End of specification.**

This is not a paper. This is a methodological contract: hypothesis ↔ algorithm ↔ stopping rule.

