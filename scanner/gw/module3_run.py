#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module3_run.py — Module 1 v3.0 (MODULE1_DESIGN_v3.0.md): full-chain Neyman belts, marginalised
50/50 over IMRPhenomXPHM and SEOBNRv5PHM. The v1.0/v2.0 files are imported unchanged.

All trial data live in the whitened complex band representation of each detector
(w_i = v * sw_i; module1_event.Network.stack writes [Re w_i, Im w_i] per detector). One trial at a
true value Lambda_t:
    d_w = h_M(theta)_w exp(i K f^3 Lambda_t) + n (unit complex Gaussian) + calibration draw,
followed by the exact chain: refinement of a common time shift (±100 ms, 0.05 ms grid) and complex
amplitude against h_ref (fast_refine: same result as module1_offsource.refine, one inverse FFT per
detector), residual, projection on the subspace complement, profile scan (chi2 curve on the grid).

Stages:
  samples   per pair: the first 400 samples of sets A and B (v1.0 build positions), each generated
            with IMRPhenomXPHM; the first 200 of each also with SEOBNRv5PHM (same parameters, label
            settings); set G (8 usable samples after C, E, F) with both models. Isolated workers.
  test      fast_refine vs the frozen refine on a stored v2.0 off-source row; short synthetic run
  belts     catalog belts per cohort and subspace at 11 truths x 400 trials
  rule0     independent catalog trials at 15 truths (11 + 4 intermediate): union coverage, grid edge
  checks    per-event full-chain belts at Lambda = 0, ±2 sigma_lin (subspace A) and the real-noise
            checks 0a, V1', V2' on new segments (ladder k >= 24) with set G
  rule3     calibration-only, noise-free catalog injections
  gate      collects rules 0, 0a, V1', V2', 3 into v3_gate.json
  onsource  the (non-blind) on-source sets, only frozen, gated and with --unblind
  hash      code hash; --write-manifest TAG writes FREEZE_v3.0.json
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time

import numpy as np

import module2_run as m2                     # installs the v1.0 on-source safeguard
import module1_run as run                    # noqa: E402
import module1_catalog as mc                 # noqa: E402
import module1_event as me                   # noqa: E402
import module1_offsource as mo               # noqa: E402
import supplement_onsource as so             # noqa: E402  (driver_args, strain_segment)

HERE = os.path.dirname(os.path.abspath(__file__))
N_SET, N_EOB, N_G, N_TR = 400, 200, 8, 400
TRUTHS = [0.0, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0, 4.0, -4.0, 8.0, -8.0]          # units of sigma_cat
CHECK_TRUTHS = TRUTHS + [0.75, -0.75, 3.0, -3.0]
K_MIN3 = 24
COHORTS = ("search", "confirm")
EOB_LABEL = "C00:SEOBNRv5PHM"


# ------------------------------------------------------------------ hashes / status ----
def v3_hash():
    h2, _, _ = m2.v2_hash()
    own = run.sha256_file(os.path.join(HERE, "module3_run.py"), normalise_eol=True)
    return hashlib.sha256((h2 + own).encode()).hexdigest(), h2, own


def v3_status(h3):
    man = os.path.join(HERE, "FREEZE_v3.0.json")
    if not os.path.exists(man):
        return "trial-v3", "no v3.0 freeze manifest"
    m = json.load(open(man, encoding="utf-8"))
    if m.get("code_sha256_v3") != h3:
        return "trial-v3", f"code hash differs from {m.get('tag', '?')}"
    return "official-v3", f"code hash matches {m.get('tag', '?')}"


# ------------------------------------------------------------------ chain pieces ---------
def fast_refine(dw, hw, fb, tmax=0.1, step=5e-5):
    """module1_offsource.refine on whitened complex band vectors, via one inverse FFT per
    detector: grid t = k * step, |k| <= tmax/step, first maximum of |z| in increasing t."""
    df = fb[1] - fb[0]
    N = int(round(1.0 / (df * step)))
    nb = len(fb)
    K = int(round(tmax / step))
    acc = np.zeros(N, complex)
    for i in dw:
        c = np.zeros(N, complex)
        c[:nb] = np.conj(hw[i]) * dw[i]
        acc += np.fft.ifft(c) * N
    k = np.arange(-K, K + 1)
    t = k * step
    z = acc[k % N] * np.exp(2j * np.pi * fb[0] * t)
    j = int(np.argmax(np.abs(z)))
    hh = sum(float(np.sum(np.abs(hw[i]) ** 2)) for i in hw)
    a = z[j] / hh
    ph = np.exp(-2j * np.pi * fb * t[j])
    return {i: dw[i] - a * hw[i] * ph for i in dw}, float(t[j]), a


def stack_w(ifos, w):
    return np.concatenate([x for i in ifos for x in (w[i].real, w[i].imag)])


def split_w(ifos, vec):
    nb = len(vec) // (2 * len(ifos))
    return {i: vec[2 * k * nb:(2 * k + 1) * nb] + 1j * vec[(2 * k + 1) * nb:(2 * k + 2) * nb]
            for k, i in enumerate(ifos)}


class Scan:
    """Projection and chi2 curve of one subspace on one grid (as so.chi2_curve)."""

    def __init__(self, e, wsp, grid):
        A = e[wsp]
        sc = me.ProfileScan(A, e["K"], e["href"], e["net"], e["fb"], grid)
        self.Q = A.Q.astype(np.float64)
        self.B = sc.b.astype(np.float64)
        self.g = np.diag(sc.G).copy()

    def chi(self, R):
        pr = R - self.Q @ (self.Q.T @ R)
        return -2.0 * (self.B @ pr) + self.g


class Pair:
    def __init__(self, args, k, ev, label, r5):
        self.ev, self.label, self.k = ev, label, k
        self.e = m2._ORIG_LOAD(args, ev, label, r5, run.code_hash()[0])
        e = self.e
        self.ifos, self.fb, self.K = e["meta"]["ifos"], e["fb"], e["K"]
        self.sw = e["net"].sw
        self.hw = {i: e["href"][i] * self.sw[i] for i in self.ifos}
        self.st, self.tab, self.idx, self.seg_on, self.band, self.network_h = m2.pair_setup(args, ev, label, r5, e)

    def chain_chi(self, dw, scans):
        rw, tsh, a = fast_refine(dw, self.hw, self.fb)
        R = stack_w(self.ifos, rw)
        return [s.chi(R) for s in scans], tsh, a


# ------------------------------------------------------------------ samples ---------------
def samples_path(args, ev, label):
    return os.path.join(args.cache3, f"{ev}__{label.replace(':', '_')}_samples.npz")


def eob_settings(st):
    st2 = json.loads(json.dumps(st))
    for k, v in list(st2.items()):
        if isinstance(v, str) and "IMRPhenomXPHM" in v:
            st2[k] = "SEOBNRv5PHM"
    if isinstance(st2.get("wf_args"), dict):
        st2["wf_args"] = {k: v for k, v in st2["wf_args"].items()
                          if not any(t in k for t in ("PhenomX", "SpinTaylor", "Prec", "MB"))}
    return st2


def project_w(p, s, pol):
    f, hp, hc = pol
    if len(f) != len(p.band) or not np.allclose(f[p.band], p.fb):
        raise ValueError("frequency grid differs from the cache")
    calf = {i: p.e["sb"][i].mean_factor(p.fb) for i in p.ifos}
    out = {i: me.rc.project(f, hp, hc, i, s, p.seg_on)[p.band] * calf[i] * p.sw[i] for i in p.ifos}
    if not all(np.all(np.isfinite(out[i])) for i in p.ifos):
        raise ValueError("nonfinite")
    return out


def gen_set(p, label, st, samples, timeout):
    res = m2.generate_isolated(label, json.loads(json.dumps(st)), samples, timeout)
    W = np.zeros((len(samples), len(p.ifos), len(p.fb)), np.complex64)
    ok = np.zeros(len(samples), bool)
    for m, (s, (good, pol, _)) in enumerate(zip(samples, res)):
        if good:
            try:
                w = project_w(p, s, pol)
                W[m] = np.array([w[i] for i in p.ifos])
                ok[m] = True
            except Exception:
                pass
    return W, ok


def stage_samples(args, h3, status):
    os.makedirs(args.cache3, exist_ok=True)
    for coh in COHORTS:
        a = so.driver_args(coh)
        a.cache3 = args.cache3
        for k, (ev, label, r5) in enumerate(run.select(a)):
            out = samples_path(args, ev, label)
            if os.path.exists(out):
                print(f"=== {ev} {label}: samples cached, skipped", flush=True)
                continue
            t0 = time.time()
            p = Pair(a, k, ev, label, r5)
            n, fails = p.e["meta"]["n"], p.e["meta"]["fails"]
            if fails != 0:
                print(f"  WARNING {ev}: v1.0 build had {fails} failures; A/B positions by index", flush=True)
            take = lambda lo, cnt: [{kk: float(p.tab[kk][int(p.idx[q])]) for kk in me.NEED} for q in range(lo, lo + cnt)]
            SA, SB = take(0, N_SET), take(n, N_SET)
            # set G: C, E, F are the first 24 usable (in-process) after A and B; G the next 8
            pos, seen, SG = 2 * n + fails, 0, []
            while len(SG) < N_G and pos < len(p.idx):
                s = {kk: float(p.tab[kk][int(p.idx[pos])]) for kk in me.NEED}
                pos += 1
                try:
                    p.network_h(s)
                except Exception:
                    continue
                if seen < 24:
                    seen += 1
                    continue
                SG.append(s)
            st2 = eob_settings(p.st)
            WA, okA = gen_set(p, label, p.st, SA, args.worker_timeout)
            WB, okB = gen_set(p, label, p.st, SB, args.worker_timeout)
            EA, okEA = gen_set(p, EOB_LABEL, st2, SA[:N_EOB], args.worker_timeout)
            EB, okEB = gen_set(p, EOB_LABEL, st2, SB[:N_EOB], args.worker_timeout)
            WG, okG = gen_set(p, label, p.st, SG, args.worker_timeout)
            EG, okEG = gen_set(p, EOB_LABEL, st2, SG, args.worker_timeout)
            np.savez(out, WA=WA, okA=okA, WB=WB, okB=okB, EA=EA, okEA=okEA, EB=EB, okEB=okEB,
                     WG=WG, okG=okG, EG=EG, okEG=okEG, code_sha256_v3=h3)
            print(f"=== {ev} {label}: XPHM A {okA.sum()}/{N_SET} B {okB.sum()}/{N_SET} G {okG.sum()}/{len(SG)}; "
                  f"SEOBNRv5PHM A {okEA.sum()}/{N_EOB} B {okEB.sum()}/{N_EOB} G {okEG.sum()}/{len(SG)}  "
                  f"({time.time() - t0:.0f} s)", flush=True)


def load_samples(args, p):
    z = np.load(samples_path(args, p.ev, p.label))
    to = lambda W: [{i: W[m, j].astype(complex) for j, i in enumerate(p.ifos)} for m in range(W.shape[0])]
    return {key: (to(z[key]), z["ok" + key[1:]] if key[0] == "W" else z["ok" + key]) for key in
            ("WA", "WB", "EA", "EB", "WG", "EG")}


def trial_sample(S, j, shift, subspace):
    """Belt of subspace A uses set-B samples and vice versa. Even trials: IMRPhenomXPHM,
    odd trials: SEOBNRv5PHM when that sample was generated, otherwise IMRPhenomXPHM."""
    X, okX = S["WB"] if subspace == "A" else S["WA"]
    E, okE = S["EB"] if subspace == "A" else S["EA"]
    m = (j + shift) // 2
    if j % 2 == 1 and okE[m % N_EOB]:
        return E[m % N_EOB], "EOB"
    q = m % N_SET
    while not okX[q]:
        q = (q + 1) % N_SET
    return X[q], "XPHM"


def trial_data(p, sample, lam, rng, cal):
    ph = np.exp(1j * p.K * p.fb ** 3 * lam)
    nz = {i: rng.normal(size=len(p.fb)) + 1j * rng.normal(size=len(p.fb)) for i in p.ifos}
    cw = split_w(p.ifos, cal())
    return {i: sample[i] * ph + nz[i] + cw[i] for i in p.ifos}


# ------------------------------------------------------------------ catalog belts ----------
def cohort_setup(coh, args):
    a = so.driver_args(coh)
    a.cache3, a.worker_timeout = args.cache3, args.worker_timeout
    pairs = run.select(a)
    metas = [json.load(open(run.cache_paths(a, ev, label)[1], encoding="utf-8")) for ev, label, _ in pairs]
    sig = np.array([m["sigma_lin"] for m in metas])
    grid = np.round(np.arange(-a.half_width, a.half_width + 1e-9, a.step), 6) * sig.min()
    sig_cat = 1.0 / np.sqrt(np.sum(1 / sig ** 2))
    return a, pairs, grid, sig_cat, sig


def truth_index(grid, values):
    return [int(np.argmin(np.abs(grid - v))) for v in values]


def catalog_trials(args, coh, truths_sig, n_tr, seed_base, shift):
    """Summed chi2 curves (subspace -> array [truth, trial, grid]) for the cohort."""
    a, pairs, grid, sig_cat, _ = cohort_setup(coh, args)
    tix = truth_index(grid, [t * sig_cat for t in truths_sig])
    tot = {s: np.zeros((len(tix), n_tr, len(grid))) for s in "AB"}
    for k, (ev, label, r5) in enumerate(pairs):
        t0 = time.time()
        p = Pair(a, k, ev, label, r5)
        S = load_samples(args, p)
        scans = {s: Scan(p.e, s, grid) for s in "AB"}
        for ti, t in enumerate(tix):
            for wsp in "AB":
                rng = np.random.default_rng(seed_base + 1000 * k + 10 * ti + (0 if wsp == "A" else 1))
                cal = mc.cal_error_fn(p.e, rng)
                for j in range(n_tr):
                    smp, _ = trial_sample(S, j, shift, wsp)
                    dw = trial_data(p, smp, grid[t], rng, cal)
                    (c,), _, _ = p.chain_chi(dw, [scans[wsp]])
                    tot[wsp][ti, j] += c
        print(f"  {coh} {ev:18s} {len(tix)} truths x {n_tr} trials x 2 subspaces ({time.time() - t0:.0f} s)", flush=True)
        del p, S
    return tot, grid, tix, sig_cat


def crit_interp(grid, tix, c_truth):
    order = np.argsort(tix)
    xs, ys = np.array(tix)[order], np.array(c_truth)[order]
    return np.interp(np.arange(len(grid)), xs, ys)


def stage_belts(args, h3, status):
    for coh in COHORTS:
        tot, grid, tix, sig_cat = catalog_trials(args, coh, TRUTHS, N_TR, me.SEED + 300, 0)
        out = {"grid": grid, "tix": np.array(tix), "sig_cat": sig_cat}
        for wsp in "AB":
            q = np.array([tot[wsp][ti, :, t] - tot[wsp][ti].min(axis=1) for ti, t in enumerate(tix)])
            c = np.quantile(q, 0.90, axis=1)
            out[f"c_{wsp}"] = c
            out[f"cgrid_{wsp}"] = crit_interp(grid, tix, c)
        np.savez(f"v3_belts_{coh}.npz", code_sha256_v3=h3, status=status, **out)
        print(f"{coh}: belts at {len(tix)} truths; c_A {out['c_A'].min():.2f}-{out['c_A'].max():.2f}, "
              f"c_B {out['c_B'].min():.2f}-{out['c_B'].max():.2f}; written v3_belts_{coh}.npz (status {status})", flush=True)


def stage_rule0(args, h3, status):
    res = {}
    for coh in COHORTS:
        b = np.load(f"v3_belts_{coh}.npz")
        tot, grid, tix, sig_cat = catalog_trials(args, coh, CHECK_TRUTHS, N_TR, me.SEED + 500, N_TR // 2)
        lo, _ = mo.binomial_interval(0.90, N_TR)
        rows, edge = [], 0
        for ti, t in enumerate(tix):
            cov = 0
            for j in range(N_TR):
                inside = np.zeros(len(grid), bool)
                for wsp in "AB":
                    chi = tot[wsp][ti, j]
                    inside |= (chi - chi.min()) <= b[f"cgrid_{wsp}"]
                cov += bool(inside[t])
                if CHECK_TRUTHS[ti] == 0.0:
                    edge += bool(inside[0] or inside[-1])
            rows.append({"truth_sigma_cat": CHECK_TRUTHS[ti], "coverage": cov / N_TR, "pass": bool(cov / N_TR >= lo)})
            print(f"  rule 0 {coh} truth {CHECK_TRUTHS[ti]:+5.2f} sigma_cat: coverage {cov}/{N_TR} = {cov / N_TR:.3f} "
                  f"(>= {lo:.3f})", flush=True)
        ok = all(r["pass"] for r in rows) and edge / N_TR <= 0.01
        res[coh] = {"truths": rows, "edge_fraction": edge / N_TR, "pass": bool(ok)}
        print(f"rule 0 {coh}: {'pass' if ok else 'FAIL'} (edge {edge / N_TR:.3f})", flush=True)
    json.dump({"status": status, "code_sha256_v3": h3, "rule0": res}, open("v3_rule0.json", "w"), indent=1)


# ------------------------------------------------------------------ per-event checks -------
def used_segments_all(ev, label):
    out = list(m2.v1_segments(ev, label))
    for fn in ("v2_offsource_search_xphm.csv", "v2_offsource_confirm_xphm.csv",
               "supplement_validation_v1.csv", "supplement_validation_v2.csv"):
        if os.path.exists(fn):
            for r in csv.DictReader(open(fn)):
                if r["event"] == ev and r["label"] == label:
                    out.append(float(r["seg_start"]))
    return sorted(set(out))


def stage_checks(args, h3, status):
    out = "v3_checks.csv"
    keys = ["cohort", "event", "check", "seg_start", "lambda_true", "q_true", "c_true", "covered", "estimate_sigma"]
    fo = open(out, "w", newline="", encoding="utf-8")
    w = csv.DictWriter(fo, fieldnames=keys)
    w.writeheader()
    for coh in COHORTS:
        a = so.driver_args(coh)
        for k, (ev, label, r5) in enumerate(run.select(a)):
            t0 = time.time()
            p = Pair(a, k, ev, label, r5)
            S = load_samples(args, p)
            sig = p.e["sigma_lin"]
            own = np.round(np.arange(-a.half_width, a.half_width + 1e-9, a.step), 6) * sig
            scan = Scan(p.e, "A", own)
            tz = int(np.argmin(np.abs(own)))
            tt = {0: tz, +2: tz + int(round(2 / a.step)), -2: tz - int(round(2 / a.step))}
            crit = {}
            for key, t in tt.items():                      # full-chain event belts, subspace A
                rng = np.random.default_rng(me.SEED + 700 + 1000 * k + (key + 2))
                cal = mc.cal_error_fn(p.e, rng)
                qs = []
                for j in range(N_TR):
                    smp, _ = trial_sample(S, j, 0, "A")
                    (c,), _, _ = p.chain_chi(trial_data(p, smp, own[t], rng, cal), [scan])
                    qs.append(c[t] - c.min())
                crit[key] = float(np.quantile(qs, 0.90))
            old = used_segments_all(ev, label)
            dur = p.st["duration"]
            windows = mo.forbidden_windows(a.events_list, dur)
            segs = []
            for kk in range(K_MIN3, K_MIN3 + 8 * N_G):
                for side in (-1, +1):
                    if len(segs) >= N_G:
                        break
                    s0 = p.seg_on + side * kk * (dur + mo.OFF_GAP)
                    if any(abs(s0 - o) < dur for o in old) or not mo.segment_allowed(s0, dur, windows)[0]:
                        continue
                    try:
                        segs.append((s0, so.strain_segment(a, ev, p.ifos, p.st, p.band, s0)))
                    except ValueError:
                        pass
            WG, okG = S["WG"]
            EG, okEG = S["EG"]
            for m, (s0, dseg) in enumerate(segs):
                if m >= len(WG):
                    break
                base = {i: dseg[i] * p.sw[i] for i in p.ifos}
                for check in ("0a", "V1'", "V2'"):
                    if check == "V2'":
                        if not okEG[m]:
                            continue
                        inj, key = EG[m], 0
                    else:
                        if not okG[m]:
                            continue
                        inj = WG[m]
                        key = 0 if check == "0a" else (+2 if m % 2 == 0 else -2)
                    t = tt[key]
                    ph = np.exp(1j * p.K * p.fb ** 3 * own[t])
                    (c,), _, _ = p.chain_chi({i: base[i] + inj[i] * ph for i in p.ifos}, [scan])
                    q = c - c.min()
                    w.writerow({"cohort": coh, "event": ev, "check": check, "seg_start": s0,
                                "lambda_true": own[t], "q_true": float(q[t]), "c_true": crit[key],
                                "covered": bool(q[t] <= crit[key]),
                                "estimate_sigma": float((own[int(np.argmin(c))] - own[t]) / sig)})
            fo.flush()
            print(f"  {coh} {ev:18s} segments {len(segs)}  c(0) {crit[0]:.2f}  ({time.time() - t0:.0f} s)", flush=True)
            del p, S
    fo.close()
    rows = list(csv.DictReader(open(out)))
    res = {}
    for check in ("0a", "V1'", "V2'"):
        r = [x for x in rows if x["check"] == check]
        nc, n = sum(x["covered"] == "True" for x in r), len(r)
        lo, _ = mo.binomial_interval(0.90, n)
        res[check] = {"covered": nc, "n": n, "lower_99": lo, "pass": bool(nc / n >= lo)}
        print(f"check {check:4s}: {nc}/{n} = {nc / n:.3f} (>= {lo:.3f}) -> {'pass' if nc / n >= lo else 'FAIL'}")
    json.dump({"status": status, "code_sha256_v3": h3, "checks": res}, open("v3_checks.json", "w"), indent=1)


# ------------------------------------------------------------------ rule 3, gate, on-source
def stage_rule3(args, h3, status):
    res = {}
    for coh in COHORTS:
        a, pairs, grid, sig_cat, _ = cohort_setup(coh, args)
        b = np.load(f"v3_belts_{coh}.npz")
        n_cal = 100
        tot = {s: np.zeros((n_cal, len(grid))) for s in "AB"}
        for k, (ev, label, r5) in enumerate(pairs):
            p = Pair(a, k, ev, label, r5)
            scans = [Scan(p.e, s, grid) for s in "AB"]
            rng = np.random.default_rng(me.SEED + 31 + 1000 * k)
            cal = mc.cal_error_fn(p.e, rng)
            for j in range(n_cal):
                cw = split_w(p.ifos, cal())
                (cA, cB), _, _ = p.chain_chi({i: p.hw[i] + cw[i] for i in p.ifos}, scans)
                tot["A"][j] += cA
                tot["B"][j] += cB
            del p
        tz = int(np.argmin(np.abs(grid)))
        est, excl = [], 0
        for j in range(n_cal):
            inside = np.zeros(len(grid), bool)
            for s in "AB":
                inside |= (tot[s][j] - tot[s][j].min()) <= b[f"cgrid_{s}"]
            excl += not inside[tz]
            est.append(grid[int(np.argmin(tot["A"][j]))])
        med = float(np.median(est))
        ok = abs(med) <= sig_cat and excl / n_cal <= 0.10
        res[coh] = {"median_estimate": med, "sigma_cat": float(sig_cat), "fraction_excluding_zero": excl / n_cal,
                    "pass": bool(ok)}
        print(f"rule 3 {coh}: median {med:+.2e} ({med / sig_cat:+.2f} sigma_cat), excluding 0: {excl}/{n_cal} -> "
              f"{'pass' if ok else 'FAIL'}", flush=True)
    json.dump({"status": status, "code_sha256_v3": h3, "rule3": res}, open("v3_rule3.json", "w"), indent=1)


def stage_gate(args, h3, status):
    r0 = json.load(open("v3_rule0.json"))
    ch = json.load(open("v3_checks.json"))
    r3 = json.load(open("v3_rule3.json"))
    same = all(x.get("code_sha256_v3") == h3 for x in (r0, ch, r3))
    ok = (same and all(v["pass"] for v in r0["rule0"].values()) and all(v["pass"] for v in ch["checks"].values())
          and all(v["pass"] for v in r3["rule3"].values()))
    json.dump({"status": status, "code_sha256_v3": h3, "same_code": same, "all_pass": bool(ok),
               "rule0": r0["rule0"], "checks": ch["checks"], "rule3": r3["rule3"]}, open("v3_gate.json", "w"), indent=1)
    print(f"gate: {'ALL RULES PASS' if ok else 'A RULE FIRED (or results from different code): no on-source computation'}")


def stage_onsource(args, h3, status):
    if not args.unblind or status != "official-v3":
        raise SystemExit("refused: needs --unblind and the frozen v3.0 code")
    g = json.load(open("v3_gate.json"))
    if not g["all_pass"] or g["code_sha256_v3"] != h3:
        raise SystemExit("refused: no passing v3_gate.json for this code")
    if args.cohort == "confirm" and not os.path.exists("v3_result_search.json"):
        raise SystemExit("refused: the search result must exist (and be committed) first")
    out = f"v3_result_{args.cohort}.json"
    if os.path.exists(out):
        raise SystemExit(f"refused: {out} exists")
    a, pairs, grid, sig_cat, _ = cohort_setup(args.cohort, args)
    b = np.load(f"v3_belts_{args.cohort}.npz")
    tot = {s: np.zeros(len(grid)) for s in "AB"}
    for k, (ev, label, r5) in enumerate(pairs):
        p = Pair(a, k, ev, label, r5)
        d = so.strain_segment(a, ev, p.ifos, p.st, p.band, p.seg_on)
        (cA, cB), tsh, amp = p.chain_chi({i: d[i] * p.sw[i] for i in p.ifos}, [Scan(p.e, s, grid) for s in "AB"])
        tot["A"] += cA
        tot["B"] += cB
        print(f"  {ev:18s} on-source: time shift {tsh * 1e3:+.2f} ms, |a| {abs(amp):.3f}", flush=True)
        del p
    inside = np.zeros(len(grid), bool)
    for s in "AB":
        inside |= (tot[s] - tot[s].min()) <= b[f"cgrid_{s}"]
    idx = np.where(inside)[0]
    res = {"cohort": args.cohort, "status": status, "code_sha256_v3": h3, "non_blind": True,
           "set": None if len(idx) == 0 else {
               "lower": float(grid[idx[0]]), "upper": float(grid[idx[-1]]),
               "contiguous": bool(idx[-1] - idx[0] + 1 == len(idx)),
               "contains_zero": bool(inside[int(np.argmin(np.abs(grid)))]),
               "touches_edge": bool(idx[0] == 0 or idx[-1] == len(grid) - 1)}}
    json.dump(res, open(out, "w"), indent=1)
    s = res["set"]
    print(f"\n90% set ({args.cohort}, v3.0, non-blind re-analysis): "
          + ("empty" if s is None else f"[{s['lower']:+.3e}, {s['upper']:+.3e}] m^2, contains 0: {s['contains_zero']}"))
    print(f"written {out}; commit it before any further on-source computation")


def stage_test(args, h3, status):
    a = so.driver_args("search")
    ev, label, r5 = run.select(a)[0]
    p = Pair(a, 0, ev, label, r5)
    row = next(r for r in csv.DictReader(open("v2_offsource_search_xphm.csv"))
               if r["event"] == ev and r["label"] == label and r["inj"] == "C")
    s = {kk: float(p.tab[kk][int(row["e_index"])]) for kk in me.NEED}
    hE = p.network_h(s)
    dseg = so.strain_segment(a, ev, p.ifos, p.st, p.band, float(row["seg_start"]))
    d = {i: dseg[i] + hE[i] for i in p.ifos}
    Seff = {i: 4 * (p.fb[1] - p.fb[0]) / p.sw[i] ** 2 for i in p.ifos}
    r1, t1, a1 = mo.refine(d, p.e["href"], Seff, p.fb[1] - p.fb[0], fb=p.fb)
    rw, t2, a2 = fast_refine({i: d[i] * p.sw[i] for i in p.ifos}, p.hw, p.fb)
    rel = max(float(np.max(np.abs(r1[i] * p.sw[i] - rw[i])) / np.max(np.abs(rw[i]))) for i in p.ifos)
    own = np.round(np.arange(-a.half_width, a.half_width + 1e-9, a.step), 6) * p.e["sigma_lin"]
    tz = int(np.argmin(np.abs(own)))
    c = Scan(p.e, "A", own).chi(stack_w(p.ifos, rw))
    q0 = c[tz] - c.min()
    print(f"fast_refine vs frozen refine ({ev}): t {t1:+.5f} vs {t2:+.5f}; residual rel diff {rel:.1e}; "
          f"q0 {q0:.6f} vs stored {float(row['q0']):.6f} -> "
          f"{'PASS' if rel < 1e-9 and abs(q0 - float(row['q0'])) < 1e-6 * max(1, abs(q0)) else 'FAIL: do not freeze'}",
          flush=True)
    if os.path.exists(samples_path(args, ev, label)):
        S = load_samples(args, p)
        scan = Scan(p.e, "A", own)
        rng = np.random.default_rng(1)
        cal = mc.cal_error_fn(p.e, rng)
        t0 = time.time()
        for j in range(20):
            smp, model = trial_sample(S, j, 0, "A")
            p.chain_chi(trial_data(p, smp, own[tz], rng, cal), [scan])
        print(f"synthetic smoke test: 20 full-chain trials in {time.time() - t0:.1f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["samples", "test", "belts", "rule0", "checks", "rule3", "gate", "onsource", "hash"])
    ap.add_argument("--cohort", choices=list(COHORTS), default="search")
    ap.add_argument("--unblind", action="store_true")
    ap.add_argument("--cache3", default=os.path.expanduser("~/gw/cache_module3"))
    ap.add_argument("--worker-timeout", type=int, default=180)
    ap.add_argument("--write-manifest", default="", metavar="TAG")
    args = ap.parse_args()
    h3, h2, own = v3_hash()
    status, why = v3_status(h3)
    print(f"v2.0 code SHA-256 {h2}\nmodule3_run.py SHA-256 {own}\nv3.0 code SHA-256 {h3}\nstatus: {status} ({why})",
          flush=True)
    if args.stage == "hash":
        if args.write_manifest:
            man = os.path.join(HERE, "FREEZE_v3.0.json")
            if os.path.exists(man):
                raise SystemExit(f"{man} exists; refusing to overwrite")
            json.dump({"tag": args.write_manifest, "code_sha256_v3": h3, "code_sha256_v2": h2,
                       "module3_run_sha256": own, "design": "MODULE1_DESIGN_v3.0.md",
                       "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
                      open(man, "w", encoding="utf-8"), indent=1)
            print(f"written {man}")
        return
    if args.stage not in ("samples", "test") and status != "official-v3":
        raise SystemExit("refused: belts, rules and on-source need the frozen v3.0 code")
    {"samples": stage_samples, "test": stage_test, "belts": stage_belts, "rule0": stage_rule0,
     "checks": stage_checks, "rule3": stage_rule3, "gate": stage_gate,
     "onsource": stage_onsource}[args.stage](args, h3, status)


if __name__ == "__main__":
    main()
