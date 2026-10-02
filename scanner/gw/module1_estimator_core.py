#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_estimator_core.py — Module 1 core estimator (MODULE1_DESIGN_v0.9 §6).

Pure linear algebra, no lalsuite dependency. Meant to be tested on synthetic
data (this file's __main__) before being wired to real strain/PE data.

THE BUG THIS FILE FIXES (found by synthetic injection, before any real
residual was touched):

  §6 as drafted treats the detector-calibration directions (amplitude and
  phase spline nodes, ~10 log-spaced nodes each) as FREE nuisance directions,
  exactly like the posterior-PCA directions, the coalescence-time direction,
  and the overall amplitude/phase directions. But a calibration phase error
  is a smooth function of f described by a handful of spline nodes, and so
  is the f^3 dispersion phase. Freely projecting out unconstrained
  calibration directions can absorb essentially all of T3: the self-test
  below shows survival -> ~0 and poor recovery of an injected Lambda.

  FIX: calibration uncertainty is not unknown-and-unbounded, it is
  known-and-bounded (a few % in amplitude, a few degrees in phase; the
  envelope is in the PE release). So calibration directions must not be
  projected out for free. They are treated as a constrained noise
  covariance  C = I + J J^T  (J = calibration basis columns scaled by the
  envelope sigma(f)), and the estimator becomes generalized least squares
  (GLS) with C^{-1} instead of an orthogonal projector. Only directions
  that really are free (posterior PCA, coalescence time, overall amplitude
  and phase) are still projected out exactly.

  Woodbury identity avoids ever forming the dense C or C^{-1}:
      C^{-1} v = v - J (I + J^T J)^{-1} J^T v

  All the expensive linear algebra (the free-subspace SVD, the Woodbury
  factor) depends only on the event/model/subspace choice, NOT on the data
  r. It is therefore computed ONCE per event in the EventSetup classes and
  reused for every null/injection trial -- this is also how the real
  pipeline must be structured, since §7-9 run many hundreds of trials per
  event.

Run:
    python module1_estimator_core.py
"""
import time

import numpy as np

RNG_SEED = 20261001


def stack(v, sw):
    """Complex frequency-domain vector -> real stacked vector, so the real
    dot product equals Re<a,b> under the noise weight sw."""
    return np.concatenate([v.real * sw, v.imag * sw])


def toy_reference_waveform(f):
    """Dimensionless, O(1)-amplitude toy waveform. Deliberately NOT in real
    strain units (~1e-21): this file tests subspace/projector algebra, and
    mixing a unit-amplitude toy waveform with a realistic-scale PSD would
    make the injected signal trivially dominate (or vanish under) unit-
    variance noise regardless of which estimator is used, hiding the very
    effect being tested. f and h_ref are consistently dimensionless; real
    use multiplies by the physical T_p prefactor of MODULE1_DESIGN §2."""
    amp = np.exp(-((f - 5.0) / 3.0) ** 2 / 2)
    phase = 0.3 * f ** 2
    return amp * np.exp(1j * phase)


def templates(f, h_ref):
    return {p: 1j * f ** p * h_ref for p in (2, 3, 4)}


def posterior_direction_samples(f, h_ref, n_samples, rng):
    """Toy stand-in for real PCA-from-PE-samples (F3): smooth, low-
    effective-rank perturbations, NOT full-rank white noise, as a real
    posterior gives."""
    n_modes = 6
    basis = np.stack([np.cos(k * np.pi * (f - f.min()) / (f.max() - f.min()))
                       for k in range(1, n_modes + 1)])
    scales = 0.08 * np.array([1.0, 0.6, 0.4, 0.25, 0.15, 0.08])
    coeffs = rng.normal(size=(n_samples, n_modes)) * scales
    delta_phase = coeffs @ basis
    delta_amp = 0.006 * (coeffs[:, :3] @ basis[:3])
    return (delta_amp + 1j * delta_phase) * h_ref[None, :]


def calibration_basis(f, n_nodes=10):
    """~n_nodes log-spaced raised-cosine bumps (stand-in for the real
    cubic-spline calibration nodes)."""
    logf = np.log(f)
    nodes = np.linspace(logf.min(), logf.max(), n_nodes)
    width = (nodes[1] - nodes[0]) * 1.1
    bumps = []
    for c in nodes:
        x = np.clip((logf - c) / width, -1, 1)
        bumps.append(0.5 * (1 + np.cos(np.pi * x)))
    return np.stack(bumps)


def calibration_envelope(f, amp_sigma=0.05, phase_sigma_rad=np.deg2rad(3.0)):
    """Toy envelope (flat in f here; the real envelope from the PE release
    is frequency-dependent and must replace this)."""
    return amp_sigma * np.ones_like(f), phase_sigma_rad * np.ones_like(f)


def calibration_directions(f, h_ref, n_nodes=10):
    bumps = calibration_basis(f, n_nodes)
    amp_dirs = bumps * h_ref[None, :]
    phase_dirs = 1j * bumps * h_ref[None, :]
    return np.concatenate([amp_dirs, phase_dirs], axis=0)


def free_subspace(f, h_ref, post_samples, sw, frac=0.99):
    D = np.stack([stack(d, sw) for d in post_samples]).T
    U, S, _ = np.linalg.svd(D, full_matrices=False)
    cum = np.cumsum(S ** 2) / np.sum(S ** 2)
    k = int(np.searchsorted(cum, frac) + 1)
    tdir = stack(-2j * np.pi * f * h_ref, sw)
    adir = stack(h_ref, sw)
    pdir = stack(1j * h_ref, sw)
    Qraw = np.column_stack([U[:, :k], tdir, adir, pdir])
    Q, _ = np.linalg.qr(Qraw)
    return Q, k


def project_out(v, Q):
    return v - Q @ (Q.T @ v)


class EventSetupOld:
    """OLD (buggy) treatment: calibration stacked into the SAME free
    subspace and projected out exactly. Kept only to reproduce the bug."""

    def __init__(self, f, h_ref, post_samples, cal_dirs, sw, frac=0.99):
        Qfree, k = free_subspace(f, h_ref, post_samples, sw, frac)
        Jcols = np.stack([stack(c, sw) for c in cal_dirs]).T
        self.Qall, _ = np.linalg.qr(np.column_stack([Qfree, Jcols]))
        self.k = k
        T3 = templates(f, h_ref)[3]
        T3_raw_norm = np.linalg.norm(stack(T3, sw))
        self.T3perp = project_out(stack(T3, sw), self.Qall)
        self.denom = float(self.T3perp @ self.T3perp)
        self.s3 = (np.linalg.norm(self.T3perp) / T3_raw_norm
                   if T3_raw_norm > 0 else np.nan)

    def estimate(self, r):
        rperp = project_out(r, self.Qall)
        lam = float(rperp @ self.T3perp / self.denom) if self.denom > 0 else np.nan
        sigma = 1.0 / np.sqrt(self.denom) if self.denom > 0 else np.inf
        return lam, sigma


class EventSetupGLS:
    """FIXED treatment: free directions projected out exactly; calibration
    enters as a constrained GLS covariance C = I + J J^T via Woodbury."""

    def __init__(self, f, h_ref, post_samples, cal_dirs_unscaled, sw,
                 amp_sigma, phase_sigma, frac=0.99):
        n_nodes = cal_dirs_unscaled.shape[0] // 2
        Qfree, k = free_subspace(f, h_ref, post_samples, sw, frac)
        self.k = k
        self.Qfree = Qfree

        scale = np.concatenate([np.full(n_nodes, amp_sigma.mean()),
                                 np.full(n_nodes, phase_sigma.mean())])
        Jraw = np.stack([stack(c, sw) for c in cal_dirs_unscaled]).T
        J = Jraw * scale[None, :]
        J1 = np.column_stack([project_out(J[:, i], Qfree) for i in range(J.shape[1])])
        self.J1 = J1

        M = np.eye(J1.shape[1]) + J1.T @ J1
        self.Minv_J1T = np.linalg.solve(M, J1.T)   # precomputed Woodbury factor

        T3 = templates(f, h_ref)[3]
        T3_raw_norm = np.linalg.norm(stack(T3, sw))
        self.T3_1 = project_out(stack(T3, sw), Qfree)
        self.Cinv_T3 = self._cinv(self.T3_1)
        self.denom = float(self.T3_1 @ self.Cinv_T3)
        self.s3 = np.sqrt(self.denom) / T3_raw_norm if T3_raw_norm > 0 else np.nan

    def _cinv(self, v):
        return v - self.J1 @ (self.Minv_J1T @ v)

    def estimate(self, r):
        r1 = project_out(r, self.Qfree)
        Cinv_r = self._cinv(r1)
        lam = float(self.T3_1 @ Cinv_r / self.denom) if self.denom > 0 else np.nan
        sigma = 1.0 / np.sqrt(self.denom) if self.denom > 0 else np.inf
        return lam, sigma


def draw_calibration_perturbation(h_ref, amp_sigma, phase_sigma, n_nodes, cal_basis, rng):
    a = rng.normal(scale=amp_sigma[0], size=n_nodes) @ cal_basis
    p = rng.normal(scale=phase_sigma[0], size=n_nodes) @ cal_basis
    return (a + 1j * p) * h_ref


def main():
    t0 = time.time()
    rng = np.random.default_rng(RNG_SEED)
    # dimensionless "frequency" axis, ~1 decade: enough range for the f^3
    # shape to vary substantially (by ~10^3) while staying numerically sane.
    f = np.linspace(1.0, 10.0, 2000)
    sw = np.ones_like(f)  # already-whitened units: unit-variance noise per
                           # stacked real component, by construction

    h_ref = toy_reference_waveform(f)
    T3w = stack(templates(f, h_ref)[3], sw)
    print(f"  ||T3|| (whitened, dimensionless) = {np.linalg.norm(T3w):.1f}")

    post_A = posterior_direction_samples(f, h_ref, 400, np.random.default_rng(RNG_SEED))
    n_nodes = 10
    cal_dirs = calibration_directions(f, h_ref, n_nodes=n_nodes)
    cal_basis = calibration_basis(f, n_nodes)
    amp_sigma, phase_sigma = calibration_envelope(f, amp_sigma=0.05, phase_sigma_rad=np.deg2rad(3.0))

    # chosen so that lambda_true * ||T3|| gives a modest, realistic-feeling
    # "injected SNR" of a few, comparable to the noise and calibration scale
    # -- NOT many orders of magnitude larger (which would make recovery
    # trivial regardless of which estimator is used, hiding the bug).
    lambda_true = 4.0 / np.linalg.norm(T3w)

    print("Setting up subspaces (one-time SVD + Woodbury factor)...")
    setup_old = EventSetupOld(f, h_ref, post_A, cal_dirs, sw)
    setup_gls = EventSetupGLS(f, h_ref, post_A, cal_dirs, sw, amp_sigma, phase_sigma)
    print(f"  k(free posterior components) = {setup_old.k}  ({time.time() - t0:.1f} s elapsed)\n")

    print("=== Injection recovery: 200 independent noise/calibration draws per "
          "estimator, same injected Lambda_true ===")
    print("    (a single draw is dominated by noise when survival is small; "
          "the honest comparison is mean and spread over many draws, not one number)\n")
    for label, setup in (("OLD (free calibration, the bug)", setup_old),
                          ("FIXED (GLS-constrained calibration)", setup_gls)):
        ratios = []
        for _ in range(200):
            cal_draw = draw_calibration_perturbation(h_ref, amp_sigma, phase_sigma,
                                                       n_nodes, cal_basis, rng)
            noise = rng.normal(size=2 * len(f))
            r = lambda_true * T3w + noise + stack(cal_draw, sw)
            lam, sigma = setup.estimate(r)
            if np.isfinite(lam):
                ratios.append(lam / lambda_true)
        ratios = np.array(ratios)
        print(f"  {label}")
        print(f"    survival(f3)={setup.s3:.4f}  analytic sigma/Lambda_true={setup.estimate(lambda_true*T3w)[1]/lambda_true:.2f}")
        print(f"    recovered Lambda_hat/Lambda_true: mean={ratios.mean():+.2f}  "
              f"std={ratios.std():.2f}  (n={len(ratios)})\n")

    print(f"\n({time.time() - t0:.1f} s elapsed)")
    print("\n=== Null check (Lambda=0, calibration-shaped noise only, FIXED estimator, n=300) ===")
    z_vals = []
    for _ in range(300):
        cal_draw = draw_calibration_perturbation(h_ref, amp_sigma, phase_sigma, n_nodes, cal_basis, rng)
        noise = rng.normal(size=2 * len(f))
        r0 = noise + stack(cal_draw, sw)
        lam0, sigma0 = setup_gls.estimate(r0)
        z_vals.append(lam0 / sigma0)
    z_vals = np.array(z_vals)
    print(f"  n={len(z_vals)}  mean(z)={z_vals.mean():.3f}  std(z)={z_vals.std():.3f}  "
          f"(expect ~0, ~1 for a correctly calibrated GLS estimator)")
    print(f"\nTotal time: {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
