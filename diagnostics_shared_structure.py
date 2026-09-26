"""
Why does the V2 omnibus z² test (P1) fail on real PTA sampling?

Three synthetic 68-pulsar arrays under H_0 (no common signal), all with
independent white + intrinsic red noise and quadratic removal:
  A  independent gaps, constant per-pulsar errors   (idealised)
  B  SHARED observing gaps (common telescope outages)
  C  SHARED error step: errors 5× larger before year 7 for ALL pulsars
     (mimics a backend upgrade at the same date, e.g. ASP/GASP → PUPPI/GUPPI)

Result obtained when this was written (40 trials, 150 null draws):
  A: z² FPR 0.03 | B: 0.07 (n_ij changes 32.9 → 24.4 but Fisher z absorbs it)
  C: z² FPR 1.00  ← shared heteroscedasticity is the failure mechanism
  D  SHARED high-cadence campaigns (constant errors!) and D with one TOA
     per bin: see printed output — tests the TOA-density route to the same
     shared-variance mechanism.
  E  NO shared structure at all, one TOA per bin, constant errors, but small
     errors (0.05–0.2 µs) so the 1 µs intrinsic red noise DOMINATES:
     z² FPR ≈ 0.28. The mirror-shift null is itself inaccurate for short,
     red-dominated series — a second, non-synchronous mechanism.

Real-gaps ablation chain (user's machine, 100 trials, H_0):
  real errors & sampling 0.69 → constant errors 0.43 → + one TOA/bin 0.31
  ≈ E. Together: shared time-varying variance (errors, TOA density) +
  red-dominance account for the z² failure.
Under C the signed/sky-scramble tests stay calibrated (0.00–0.10).

  python diagnostics_shared_structure.py

Author: Dimitar Kretski
License: MIT
"""
import numpy as np
import common_residual_search_v2 as v


def template(rng, shared_gaps=False, shared_err_step=False, shared_campaigns=False, n=68,
             err_range=(0.3, 2.0)):
    T = 16 * 365.25
    windows = [(s, s + rng.uniform(60, 250)) for s in rng.uniform(0, T, 12)] if shared_gaps else []
    out = {}
    for i in range(n):
        start = rng.uniform(0, 13) * 365.25
        t = np.sort(rng.uniform(start, T, int((T - start) / 25)))
        for a, b in windows:
            t = t[(t < a) | (t > b)]
        if shared_campaigns:   # high-cadence campaigns at the SAME dates for all pulsars
            for c in (3.0, 9.0, 13.5):
                c0 = c * 365.25
                if c0 > start:
                    t = np.sort(np.concatenate([t, rng.uniform(c0, c0 + 365.25, 60)]))
        e = np.full(len(t), rng.uniform(*err_range) * 1e-6)
        if shared_err_step:
            e[t < 7 * 365.25] *= 5.0
        out[f"P{i:02d}"] = v._P(t, np.zeros_like(t), e, rng.uniform(0, 360),
                               np.degrees(np.arcsin(rng.uniform(-0.6, 1))))
    return out


def run(label, n_trials=40, one_per_bin=0.0, **kw):
    rng = np.random.default_rng(5)
    tmpl = template(rng, **kw)
    rej, n_real, n_null = [], [], []
    for k in range(n_trials):
        s = v.CommonResidualSearchV2(v.synthetic_array(rng, template=tmpl, one_per_bin=one_per_bin),
                                     n_null=150,
                                     n_scramble=50, seed=k, verbose=False)
        s.test_cross_correlation()
        rej.append(s.p_values[0] <= 0.05)
        n_real.append(s.n_overlap[s.keep].mean())
        n_null.append(np.mean([v.pairwise_correlations(
            v.detrend_rows(v.mirror_shift_within_span(s.X, s.spans, s.rng)), 12)[1].mean()
            for _ in range(5)]))
    print(f"{label:<36s} z² FPR = {np.mean(rej):.2f}   mean n_ij real {np.mean(n_real):5.1f} "
          f"vs null {np.mean(n_null):5.1f}", flush=True)


if __name__ == "__main__":
    run("A: independent gaps, const errors")
    run("B: SHARED gaps", shared_gaps=True)
    run("C: SHARED error step (5x early)", shared_err_step=True)
    run("D: SHARED high-cadence campaigns", shared_campaigns=True)
    run("D + one TOA per 30-d bin", shared_campaigns=True, one_per_bin=30.0)
    run("E: RED-dominated (errors 0.05-0.2 us)", err_range=(0.05, 0.2), one_per_bin=30.0)
