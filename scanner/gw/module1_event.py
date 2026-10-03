#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_event.py — Module 1 estimator on the real geometry of an event
(MODULE1_DESIGN v0.9.8, §6; feasibility item "Estimator core" of §12).

This version builds everything that does NOT depend on the strain around the event
and tests it with synthetic noise only. No on-source strain is read, no real
residual is formed. Per event and model label:

  1. reference point = maximum-likelihood sample, the label's own waveform settings
     (f3_derivatives.settings), polarizations and detector projection exactly as in
     residual_check.py (verified on all 116 event-model pairs);
  2. calibration (module1_calibration, v0.9.8 protocol): reference multiplied by the
     bilby spline factor at the prior means; J from the prior sigmas;
  3. template prefactor K(z) = -4 pi^3 I4(z)/c^3 (module1_cosmo, Planck15), z from the
     reference luminosity distance; T_p = i K f^p h for p = 2, 3, 4;
  4. free subspace: PCA (frac >= 0.99) of the network differences
     h(theta_i) - h_ref for n posterior samples (sets A and B, disjoint, seed
     20261001; each sample with its own sky position, polarization and time), plus
     the common coalescence-time, amplitude and phase directions;
  5. generalized least squares with C = I + J1 J1^T (Woodbury), as in
     module1_estimator_core.py.

Synthetic check (whitened, real-stacked network vectors; noise N(0,1) per
component, i.e. Gaussian noise with the released PSD scaled by the window power):
  - null: noise + a calibration error drawn from the priors (exact bilby factor,
    not its linearization) -> z = Lambda_hat / sigma, expect mean 0, std 1;
  - f^3 injection at +5 sigma -> recovered / injected;
  - f^2 and f^4 injections at 5 sigma in their own template -> z of the f^3
    estimate and the fraction of trials in which the largest |z| is the injected
    shape;
  - set A vs set B estimate on the same data -> |Lambda_A - Lambda_B| / sigma;
  - linearity: max |dPsi| over the band at the injected Lambda.
  - response curve (noiseless): the exact phase factor h (exp(i Lambda K f^3) - 1)
    injected at several Lambda (in units of this event's sigma and at fixed values
    1e-12, 1e-11 m^2); recovered / injected and the signal-weighted phase
    ||dPsi h|| / ||h|| show where the linear estimator holds;
  - correlations between the projected templates T2, T3, T4 in the C^-1 metric;
  - z and sigma from the maximum-likelihood distance and from the posterior median
    distance (the choice is a v0.9.9 decision).
Diagnostics only: no thresholds here; the final numbers are fixed in v1.0.

Version 7: every sub-test of run_checks/scan_checks uses its own seeded generator, so
adding a test cannot change the numbers of the others.
Version 6: A/B systematic measured on the scan (decision level and median difference),
v0.9.13 §7.
Version 5: build() also returns the raw calibration matrix J (for the module1_run cache).
Version 4: template prefactor K(z) from the posterior median distance (v0.9.9, §5);
sigma with the maximum-likelihood z reported for comparison.
Version 3: profile scan over Lambda with the exact phase (scan_checks), Wilks coverage,
blind-region distance, and the GR-leakage test (set-B GR waveforms as data).
Version 2 (after the first run): differences stored in float32 and the Gram matrix
accumulated in blocks (the GW230627 SEOBNRv5PHM run was killed for lack of memory),
progress printed during sample generation, one CSV row written per label as soon as
it is done.

  python module1_event.py --events GW150914,GW200129_065458,GW230627_015337 --n 800
"""
import argparse
import csv
import os
import sys
import time

import h5py
import numpy as np

import f3_derivatives as f3
import module1_calibration as cal
import module1_cosmo as cosmo
import residual_check as rc

SEED = 20261001
FRAC = 0.99
N_LEAK = 60      # set-B GR differences kept for the GR-leakage test of the scan

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def posterior_table(fh, label):
    ps = fh[label]["posterior_samples"]
    if isinstance(ps, h5py.Dataset) and ps.dtype.names:
        arr = ps[()]
        return {n: np.asarray(arr[n]) for n in arr.dtype.names}
    names = [f3.text(v) for v in ps["parameter_names"][()]]
    smp = ps["samples"][()]
    return {n: smp[:, i] for i, n in enumerate(names)}


NEED = list(f3.STEPS) + ["ra", "dec", "psi", "geocent_time"]


class Network:
    """Whitened, real-stacked network vectors on the analysis band."""

    def __init__(self, f, band, ifos, Seff):
        self.f, self.band, self.ifos = f, band, ifos
        df = f[1] - f[0]
        self.sw = {i: np.sqrt(4 * df / Seff[i][band]) for i in ifos}
        self.nb = int(band.sum())

    def stack(self, per_ifo):
        out = []
        for i in self.ifos:
            v = per_ifo[i][self.band] if per_ifo[i].shape[0] == len(self.f) else per_ifo[i]
            out += [v.real * self.sw[i], v.imag * self.sw[i]]
        return np.concatenate(out)

    def single(self, ifo, v_band):
        """Stack a band vector of one detector, zeros for the others."""
        out = []
        for i in self.ifos:
            if i == ifo:
                out += [v_band.real * self.sw[i], v_band.imag * self.sw[i]]
            else:
                out += [np.zeros(self.nb), np.zeros(self.nb)]
        return np.concatenate(out)


def orthonormal(cols):
    Q, R = np.linalg.qr(np.column_stack(cols))
    keep = np.abs(np.diag(R)) > np.abs(R).max() * 1e-10
    return Q[:, keep]


def pca_basis(D, frac, block=20000):
    """Leading principal directions of the columns of D (float32, rows x n) via the
    Gram matrix accumulated in float64 blocks; memory ~ D + rows x k."""
    rows, n = D.shape
    G = np.zeros((n, n))
    for a in range(0, rows, block):
        Dc = D[a:a + block].astype(np.float64)
        G += Dc.T @ Dc
    w, V = np.linalg.eigh(G)
    w, V = w[::-1], V[:, ::-1]
    w = np.clip(w, 0, None)
    cum = np.cumsum(w) / np.sum(w)
    k = int(np.searchsorted(cum, frac) + 1)
    M = V[:, :k] / np.sqrt(w[:k])
    U = np.empty((rows, k))
    for a in range(0, rows, block):
        U[a:a + block] = D[a:a + block].astype(np.float64) @ M
    return U, k


class Estimator:
    """GLS estimator with an exactly projected free subspace Q and a constrained
    calibration covariance C = I + J1 J1^T (same algebra as EventSetupGLS in
    module1_estimator_core.py)."""

    def __init__(self, Q, J, T):
        self.Q = Q
        J1 = J - Q @ (Q.T @ J)
        self.J1 = J1
        self.MJ = np.linalg.solve(np.eye(J1.shape[1]) + J1.T @ J1, J1.T)
        self.Tp, self.CTp, self.den, self.s = {}, {}, {}, {}
        for p, t in T.items():
            tp = t - Q @ (Q.T @ t)
            ctp = self._cinv(tp)
            self.Tp[p], self.CTp[p] = tp, ctp
            self.den[p] = float(tp @ ctp)
            self.s[p] = np.sqrt(self.den[p]) / np.linalg.norm(t)

    def _cinv(self, v):
        return v - self.J1 @ (self.MJ @ v)

    def corr(self, p, q):
        return float(self.Tp[p] @ self.CTp[q]) / np.sqrt(self.den[p] * self.den[q])

    def sigma(self, p):
        return 1.0 / np.sqrt(self.den[p])

    def estimate(self, r, p):
        r1 = r - self.Q @ (self.Q.T @ r)
        return float(self.CTp[p] @ r1) / self.den[p]


def build(event, label, xphm_label, fh, ifos, n, frac, log):
    """Reference geometry, templates, calibration, and the set-A and set-B estimators."""
    sample, cfg, _ = f3.read_label(fh, label)
    st = f3.settings(cfg)
    psds = {i: np.asarray(fh[label]["psds"][i][()]) for i in ifos}
    dur, srate = st["duration"], st["srate"]
    seg_start = sample["geocent_time"] + rc.POST_TRIGGER - dur
    from scipy.signal.windows import tukey
    w2 = np.mean(tukey(int(round(dur * srate)), alpha=2 * rc.ROLL_OFF / dur) ** 2)

    f, hp, hc = rc.polarizations(label, sample, st)
    band = (f >= st["f_an"]) & (f <= st["f_high"])
    Seff = {i: np.interp(f, psds[i][:, 0], psds[i][:, 1], left=np.inf, right=np.inf) * w2 for i in ifos}
    net = Network(f, band, ifos, Seff)
    fb = f[band]

    pri, cal_source = cal.load(fh, label, xphm_label, ifos)
    sb = {i: cal.SplineBasis(pri[i]) for i in ifos}
    for i in ifos:
        if fb[0] < pri[i]["freqs"][0] * (1 - 1e-6) or fb[-1] > pri[i]["freqs"][-1] * (1 + 1e-6):
            log(f"  WARNING {i}: band {fb[0]:.1f}–{fb[-1]:.1f} Hz outside calibration nodes "
                f"{pri[i]['freqs'][0]:.1f}–{pri[i]['freqs'][-1]:.1f} Hz (spline extrapolated)")
    calf = {i: sb[i].mean_factor(fb) for i in ifos}

    def network_h(s):
        fp, sp, sc = rc.polarizations(label, s, st)
        if len(fp) != len(f) or not np.allclose(fp, f):
            raise ValueError("frequency grid mismatch")
        return {i: rc.project(f, sp, sc, i, s, seg_start)[band] * calf[i] for i in ifos}

    href = network_h(sample)
    tab = posterior_table(fh, label)
    z = cosmo.z_from_dl(sample["luminosity_distance"] * cosmo.MPC_SI)
    z_med = cosmo.z_from_dl(float(np.median(tab["luminosity_distance"])) * cosmo.MPC_SI)
    K = cosmo.K(z_med)          # v0.9.9: posterior median distance (§5)
    T = {p: net.stack({i: 1j * K * fb ** p * href[i] for i in ifos}) for p in (2, 3, 4)}
    tdir = net.stack({i: -2j * np.pi * fb * href[i] for i in ifos})
    adir = net.stack(href)
    pdir = net.stack({i: 1j * href[i] for i in ifos})
    J = np.column_stack([net.single(i, d) for i in ifos for d in sb[i].directions(fb, href[i])])

    nsamp = len(tab["log_likelihood"])
    idx = np.random.default_rng(SEED).permutation(nsamp)
    href_s = net.stack(href)
    pos, fails, est, ks = 0, 0, {}, {}
    rows = len(href_s)
    t_start = time.time()
    done = 0
    for name in ("A", "B"):
        D = np.empty((rows, n), dtype=np.float32)
        m = 0
        while m < n and pos < nsamp:
            j = idx[pos]
            pos += 1
            s = {k: float(tab[k][j]) for k in NEED}
            try:
                D[:, m] = net.stack(network_h(s)) - href_s
                m += 1
            except Exception:
                fails += 1
            done += 1
            if done % 100 == 0:
                el = time.time() - t_start
                log(f"    samples {done}/{2 * n}  ({el:.0f} s, about {el / done * (2 * n - done):.0f} s left)")
        if m < min(n, 50):
            raise ValueError(f"set {name}: only {m} usable samples (failures {fails})")
        U, k = pca_basis(D[:, :m], frac)
        if name == "B":
            gr_leak = D[:, :min(m, N_LEAK)].copy()
        del D
        Q = orthonormal([U, tdir, adir, pdir])
        del U
        est[name], ks[name] = Estimator(Q, J, T), k
    info = {"J": J, "gr_leak": gr_leak, "z": z, "z_med": z_med, "I4_ratio": cosmo.I4(z_med) / cosmo.I4(z), "K": K, "cal_source": cal_source, "fails": fails, "nb": net.nb,
            "snr_ref": float(np.linalg.norm(href_s)), "f_high": float(fb[-1]),
            "n_used": n, "nsamp": nsamp}
    return est, ks, info, sb, href, net, fb


def run_checks(est, info, sb, href, net, fb, T, trials, rng):
    """Per-test generators, as in scan_checks."""
    def gen(tag):
        return np.random.default_rng([SEED, int(np.frombuffer(tag.encode().ljust(8, b"\0")[:8],
                                                              dtype=np.uint64)[0] % (2 ** 32))])
    A, B = est["A"], est["B"]
    s3 = A.sigma(3)
    nvec = 2 * net.nb * len(net.ifos)

    def cal_error(g):
        out = {}
        for i in net.ifos:
            pr = sb[i].pri
            out[i] = (sb[i].factor(fb, g.normal(0, pr["sg_a"]), g.normal(0, pr["sg_p"])) - 1.0) * href[i]
        return net.stack(out)

    def noise(g):
        return g.normal(size=nvec) + cal_error(g)

    res = {}
    gnull = gen("null")
    z3, dab = [], []
    for _ in range(trials):
        r = noise(gnull)
        a3 = A.estimate(r, 3)
        z3.append(a3 / s3)
        dab.append(abs(a3 - B.estimate(r, 3)) / s3)
    res["null_mean_z"], res["null_std_z"] = float(np.mean(z3)), float(np.std(z3))
    res["ab_median_over_sigma"] = float(np.median(dab))

    lam = 5.0 * s3
    ratios = []
    g3 = gen("inj3")
    for _ in range(trials):
        ratios.append(A.estimate(lam * T[3] + noise(g3), 3) / lam)
    res["inj_lambda"] = lam
    res["inj_ratio_mean"], res["inj_ratio_std"] = float(np.mean(ratios)), float(np.std(ratios))
    res["inj_max_dpsi"] = float(abs(lam * info["K"]) * info["f_high"] ** 3)

    for q in (2, 4):
        aq = 5.0 * A.sigma(q)
        zq3, hit = [], 0
        gq = gen(f"inj{q}")
        for _ in range(trials):
            r = aq * T[q] + noise(gq)
            z = {p: A.estimate(r, p) / A.sigma(p) for p in (2, 3, 4)}
            zq3.append(z[3])
            hit += max(z, key=lambda p: abs(z[p])) == q
        res[f"f{q}_inj_mean_z3"] = float(np.mean(zq3))
        res[f"f{q}_inj_shape_correct"] = hit / trials
    hit = 0
    gs = gen("shape3")
    for _ in range(trials):
        r = lam * T[3] + noise(gs)
        z = {p: A.estimate(r, p) / A.sigma(p) for p in (2, 3, 4)}
        hit += max(z, key=lambda p: abs(z[p])) == 3
    res["f3_inj_shape_correct"] = hit / trials

    # noiseless response to the EXACT dispersive phase
    hnorm = np.linalg.norm(net.stack(href))
    for tag, L in (("0.1sig", 0.1 * s3), ("0.3sig", 0.3 * s3), ("1sig", s3), ("3sig", 3 * s3),
                   ("1e-12", 1e-12), ("1e-11", 1e-11)):
        dpsi = {i: L * info["K"] * fb ** 3 for i in net.ifos}
        sig = net.stack({i: href[i] * (np.exp(1j * dpsi[i]) - 1.0) for i in net.ifos})
        res[f"resp_{tag}"] = A.estimate(sig, 3) / L
        res[f"wphase_{tag}"] = float(np.linalg.norm(net.stack({i: dpsi[i] * href[i] for i in net.ifos})) / hnorm)
    res["rho23"], res["rho34"], res["rho24"] = A.corr(2, 3), A.corr(3, 4), A.corr(2, 4)
    return res


KEYS = ["event", "label", "ifos", "k_A", "k_B", "s2", "s3", "s4", "sigma3_m2", "sigma3_m2_zmaxL",
        "z", "z_med", "cos_AB", "rho23", "rho34", "rho24", "null_mean_z", "null_std_z",
        "ab_median_over_sigma", "inj_lambda", "inj_ratio_mean", "inj_ratio_std", "inj_max_dpsi",
        "f3_inj_shape_correct", "f2_inj_mean_z3", "f2_inj_shape_correct", "f4_inj_mean_z3",
        "f4_inj_shape_correct"] + [f"{a}_{t}" for t in ("0.1sig", "0.3sig", "1sig", "3sig", "1e-12", "1e-11")
                                   for a in ("resp", "wphase")] + ["scan_dchi2_at_pm1", "scan_min_dchi2_far_from_0"] + [
            f"scan_{l}_{q}" for l in ("+0sig", "+1sig", "+3sig", "-1sig") for q in ("mean", "std", "cover90", "width")] + [
            "noise_frac_abs_ge1", "ab_scan_disagree", "ab_scan_median_diff", "leak_n", "leak_mean", "leak_std", "leak_frac_abs_ge1", "cal_source", "fails", "seconds", "error"]


def exact_signal(L, K, href, net, fb):
    return net.stack({i: href[i] * (np.exp(1j * L * K * fb ** 3) - 1.0) for i in net.ifos})


class ProfileScan:
    """Profile scan over Lambda with the EXACT dispersive phase (v0.9.9 primary
    statistic). Model of the residual for a given Lambda: h_ref (exp(i Lambda K f^3) - 1);
    the free subspace Q is projected out and calibration enters through C^-1, as in
    Estimator. chi2(L) = || P (r - s(L)) ||^2_{C^-1}.

    For a grid L_g: a_g = P s(L_g), b_g = C^-1 a_g, G_gh = a_g . b_h. For data
    r = s(L_t) + n (L_t on the grid): chi2_g - chi2_t = 2 (Pn).(b_t - b_g)
    + (G_tt - 2 G_gt + G_gg), so a trial costs one matrix-vector product."""

    def __init__(self, est, K, href, net, fb, grid):
        self.grid = np.asarray(grid)
        self.est = est
        ng = len(self.grid)
        a0 = est.Q.T  # alias for speed
        b = None
        for g, L in enumerate(self.grid):
            a = exact_signal(L, K, href, net, fb)
            a = a - est.Q @ (a0 @ a)
            bg = est._cinv(a)
            if b is None:
                b = np.empty((ng, len(a)), dtype=np.float32)
                arows = np.empty((ng, len(a)), dtype=np.float32)
            b[g] = bg
            arows[g] = a
        self.b = b
        G = np.empty((ng, ng))
        for c0 in range(0, ng, 32):
            G[:, c0:c0 + 32] = b.astype(np.float64) @ arows[c0:c0 + 32].astype(np.float64).T
        self.G = 0.5 * (G + G.T)
        del arows

    def delta(self, t):
        """Noiseless chi2 difference between Lambda_t (truth) and every grid point."""
        G = self.G
        return G[t, t] - 2 * G[:, t] + np.diag(G)

    def chi2(self, n_proj, t):
        """chi2_g - chi2_t for data s(L_t) + n, n_proj = P n."""
        v = self.b.astype(np.float64) @ n_proj
        return 2 * (v[t] - v) + self.delta(t)


def scan_checks(est, info, sb, href, net, fb, trials, rng, half_width=8.0, step=0.05):
    """Each sub-test draws from its OWN generator, seeded from SEED and a fixed label, so
    that adding or removing a sub-test cannot change the numbers of the others (v0.9.13;
    before this, the GR-leakage numbers changed when the A/B test was inserted ahead of
    them in the same stream)."""
    def gen(tag):
        return np.random.default_rng([SEED, int(np.frombuffer(tag.encode().ljust(8, b"\0")[:8],
                                                              dtype=np.uint64)[0] % (2 ** 32))])
    A = est["A"]
    s3 = A.sigma(3)
    K = info["K"]
    units = np.round(np.arange(-half_width, half_width + 1e-9, step), 6)
    grid = units * s3
    sc = ProfileScan(A, K, href, net, fb, grid)
    nvec = len(A.Q)

    def cal_error(g):
        out = {}
        for i in net.ifos:
            pr = sb[i].pri
            out[i] = (sb[i].factor(fb, g.normal(0, pr["sg_a"]), g.normal(0, pr["sg_p"])) - 1.0) * href[i]
        return net.stack(out)

    res = {}
    # blind regions: smallest noiseless chi2 distance from Lambda = 0 to |Lambda| >= 2 sigma_lin
    t0 = int(np.argmin(np.abs(units)))
    d0 = sc.delta(t0)
    far = np.abs(units) >= 2.0
    res["scan_min_dchi2_far_from_0"] = float(d0[far].min())
    res["scan_dchi2_at_pm1"] = float(min(d0[np.argmin(np.abs(units - 1))], d0[np.argmin(np.abs(units + 1))]))
    for lev in (0.0, 1.0, 3.0, -1.0):
        t = int(np.argmin(np.abs(units - lev)))
        est_u, cover, width = [], 0, []
        rg = gen(f"level{lev:+g}")
        for _ in range(trials):
            n = rg.normal(size=nvec) + cal_error(rg)
            npj = n - A.Q @ (A.Q.T @ n)
            c = sc.chi2(npj, t)
            g = int(np.argmin(c))
            est_u.append(units[g])
            inside = c - c[g] <= 2.71
            cover += bool(inside[t])
            width.append(units[inside].max() - units[inside].min())
        if lev == 0.0:
            res["noise_frac_abs_ge1"] = float(np.mean(np.abs(np.array(est_u)) >= 1.0))
        tag = f"{lev:+g}sig"
        res[f"scan_{tag}_mean"] = float(np.mean(est_u))
        res[f"scan_{tag}_std"] = float(np.std(est_u))
        res[f"scan_{tag}_cover90"] = cover / trials
        res[f"scan_{tag}_width"] = float(np.median(width))
    # A/B systematic on the scan (v0.9.13): same noise, subspace A vs subspace B.
    # Decision level: do A and B disagree about whether Lambda = 0 is in the 90% interval?
    # Quantitative: median |Lambda_A - Lambda_B| of the scan estimates, in sigma_lin(A).
    B = est.get("B")
    if B is not None:
        scB = ProfileScan(B, K, href, net, fb, grid)
        dis, diff = 0, []
        g = gen("ab")
        for _ in range(trials):
            n = g.normal(size=nvec) + cal_error(g)
            inA = sc.chi2(n - A.Q @ (A.Q.T @ n), t0)
            inB = scB.chi2(n - B.Q @ (B.Q.T @ n), t0)
            gA, gB = int(np.argmin(inA)), int(np.argmin(inB))
            dis += bool((inA[t0] - inA[gA] <= 2.71) != (inB[t0] - inB[gB] <= 2.71))
            diff.append(abs(units[gA] - units[gB]))
        res["ab_scan_disagree"] = dis / trials
        res["ab_scan_median_diff"] = float(np.median(diff))
        del scB

    # GR leakage: data = a real GR waveform of set B (not used for the subspace of A)
    # minus h_ref, plus noise. A GR-only signal must give Lambda ~ 0 with the spread of
    # pure noise; a shift or wider spread means GR nonlinearity leaks into Lambda.
    leak = info.get("gr_leak")
    if leak is not None and leak.shape[1] > 0:
        est_u = []
        g = gen("leak")
        for j in range(leak.shape[1]):
            for _ in range(max(1, trials // leak.shape[1])):
                n = leak[:, j].astype(np.float64) + g.normal(size=nvec) + cal_error(g)
                npj = n - A.Q @ (A.Q.T @ n)
                c = sc.chi2(npj, t0)
                est_u.append(units[int(np.argmin(c))])
        est_u = np.array(est_u)
        res["leak_n"] = len(est_u)
        res["leak_mean"] = float(np.mean(est_u))
        res["leak_std"] = float(np.std(est_u))
        res["leak_frac_abs_ge1"] = float(np.mean(np.abs(est_u) >= 1.0))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="GW150914,GW200129_065458,GW230627_015337")
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    ap.add_argument("--n", type=int, default=800, help="posterior samples per set (design: 800)")
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--labels", default="both", choices=["both", "xphm", "eob"])
    ap.add_argument("--out", default="module1_geometry_check.csv")
    ap.add_argument("--scan-trials", type=int, default=100, help="noise trials per scan level")
    args = ap.parse_args()
    rows5 = {r["commonName"]: r for r in csv.DictReader(open(args.f5, encoding="utf-8"))}
    new_file = not os.path.exists(args.out)
    fout = open(args.out, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fout, fieldnames=KEYS)
    if new_file:
        w.writeheader()
    for ev in args.events.split(","):
        r5 = rows5[ev]
        ifos = [i for i in r5["ifos"].split(",") if i]
        labels = []
        if args.labels in ("both", "xphm"):
            labels.append(r5["xphm_label"])
        if args.labels in ("both", "eob") and r5.get("eob_label"):
            labels.append(r5["eob_label"])
        with h5py.File(os.path.join(args.pe_dir, r5["file"].split(" ")[0]), "r") as fh:
            for label in labels:
                t0 = time.time()
                print(f"\n=== {ev}  {label}  detectors {','.join(ifos)}", flush=True)
                try:
                    est, ks, info, sb, href, net, fb = build(ev, label, r5["xphm_label"], fh, ifos,
                                                             args.n, FRAC, lambda m: print(m, flush=True))
                    A, B = est["A"], est["B"]
                    K = info["K"]
                    # full (unprojected) templates: what the data would contain
                    T = {p: net.stack({i: 1j * K * fb ** p * href[i] for i in ifos}) for p in (2, 3, 4)}
                    cosAB = float(A.Tp[3] @ B.Tp[3] / (np.linalg.norm(A.Tp[3]) * np.linalg.norm(B.Tp[3])))
                    s3m = A.sigma(3) * info["I4_ratio"]          # what sigma would be with z maxL
                    print(f"  z maxL {info['z']:.4f}, z median {info['z_med']:.4f}   ref SNR {info['snr_ref']:.1f}   "
                          f"band to {info['f_high']:.0f} Hz   calibration from {info['cal_source']}   "
                          f"samples {info['nsamp']} (failures {info['fails']})")
                    print(f"  k A/B {ks['A']}/{ks['B']}   survival f2 {A.s[2]:.4f}  f3 {A.s[3]:.4f}  f4 {A.s[4]:.4f}   "
                          f"sigma(Lambda) {A.sigma(3):.3e} m^2 (z median, used; with z maxL {s3m:.3e})   cos(T3perp A,B) {cosAB:.4f}")
                    res = run_checks(est, info, sb, href, net, fb, T, args.trials,
                                     np.random.default_rng(SEED + 1))
                    print(f"  template correlations: rho(f2,f3) {res['rho23']:+.4f}  rho(f3,f4) {res['rho34']:+.4f}  "
                          f"rho(f2,f4) {res['rho24']:+.4f}")
                    print(f"  null: mean z {res['null_mean_z']:+.3f}  std z {res['null_std_z']:.3f}   "
                          f"median |Lambda_A - Lambda_B|/sigma {res['ab_median_over_sigma']:.3f}")
                    print(f"  linear f3 at +5 sigma: recovered/injected {res['inj_ratio_mean']:.3f} ± "
                          f"{res['inj_ratio_std']:.3f}; shape correct {res['f3_inj_shape_correct']:.2f}")
                    print(f"  f2 at 5 sigma: mean z(f3) {res['f2_inj_mean_z3']:+.2f}, shape correct "
                          f"{res['f2_inj_shape_correct']:.2f};  f4 at 5 sigma: mean z(f3) "
                          f"{res['f4_inj_mean_z3']:+.2f}, shape correct {res['f4_inj_shape_correct']:.2f}")
                    print("  exact-phase response (noiseless):  Lambda   recovered/injected   weighted |dPsi| [rad]")
                    for t in ("0.1sig", "0.3sig", "1sig", "3sig", "1e-12", "1e-11"):
                        print(f"      {t:>7s}   {res['resp_' + t]:8.4f}   {res['wphase_' + t]:.3e}")
                    ts = time.time()
                    res.update(scan_checks(est, info, sb, href, net, fb, args.scan_trials,
                                           np.random.default_rng(SEED + 2)))
                    print(f"  profile scan (exact phase, grid ±8 sigma_lin, step 0.05; {time.time() - ts:.0f} s):")
                    print(f"    noiseless dchi2 from Lambda=0: at ±1 sigma_lin {res['scan_dchi2_at_pm1']:.2f}; "
                          f"smallest at |Lambda| >= 2 sigma_lin {res['scan_min_dchi2_far_from_0']:.2f}  "
                          f"(small = blind region)")
                    print("    injected   mean estimate   std   90% coverage   median 90% width   [units of sigma_lin]")
                    for lev in ("+0sig", "+1sig", "+3sig", "-1sig"):
                        print(f"    {lev:>6s}    {res[f'scan_{lev}_mean']:+8.3f}    {res[f'scan_{lev}_std']:6.3f}"
                              f"      {res[f'scan_{lev}_cover90']:.2f}          {res[f'scan_{lev}_width']:.2f}")
                    if "ab_scan_disagree" in res:
                        print(f"    A/B on the scan: decision-level disagreement "
                              f"{res['ab_scan_disagree']:.3f} (rule <= 0.05); median |Lambda_A - Lambda_B| "
                              f"{res['ab_scan_median_diff']:.3f} sigma_lin (rule <= 0.5)")
                    if "leak_mean" in res:
                        print(f"    GR leakage (set-B GR waveforms + noise, {res['leak_n']} trials): mean "
                              f"{res['leak_mean']:+.3f}  std {res['leak_std']:.3f}  P(|est|>=1) "
                              f"{res['leak_frac_abs_ge1']:.3f}   vs pure noise: std {res['scan_+0sig_std']:.3f}  "
                              f"P(|est|>=1) {res['noise_frac_abs_ge1']:.3f}")
                    row = {"event": ev, "label": label, "ifos": ",".join(ifos), "k_A": ks["A"], "k_B": ks["B"],
                           "s2": A.s[2], "s3": A.s[3], "s4": A.s[4], "sigma3_m2": A.sigma(3),
                           "sigma3_m2_zmaxL": s3m, "z": info["z"], "z_med": info["z_med"], "cos_AB": cosAB,
                           "cal_source": info["cal_source"], "fails": info["fails"], **res, "error": ""}
                except Exception as e:
                    print(f"  ERROR {type(e).__name__}: {str(e)[:200]}", flush=True)
                    row = {"event": ev, "label": label, "error": f"{type(e).__name__}: {str(e)[:200]}"}
                row["seconds"] = round(time.time() - t0)
                print(f"  ({row['seconds']} s)", flush=True)
                w.writerow({k: row.get(k, "") for k in KEYS})
                fout.flush()
    fout.close()
    print(f"\nappended to {args.out}")


if __name__ == "__main__":
    main()
