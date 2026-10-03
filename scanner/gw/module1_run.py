#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_run.py — driver for the Module 1 pre-freeze computations, with a verified cache
(MODULE1_DESIGN v0.9.11, §6, §7, §10, §12).

Stages (all synthetic; NO on-source data):
  build    per event and model: geometry of module1_event.build (reference, calibration,
           K(z), PCA subspaces of sets A and B, calibration directions, GR-leakage
           waveforms) -> cache/<event>__<label>.npz + .json (float32 arrays)
  checks   per event and model, from the cache: synthetic checks of module1_event
           (null, injections, exact-phase response, profile scan, GR leakage)
  catalog  from the cache: common absolute grid, Neyman belts, expected interval,
           coverage, N_eff, leave-one-out, joint f^2/f^3/f^4 conditioning

Safeguards required before v1.0:
  1. On-source data cannot be read: the strain readers of residual_check are replaced,
     before anything else is imported, by functions that stop the program. No stage
     projects a real residual.
  2. Every cache entry records: SHA-256 of the code (all modules listed in CODE_FILES,
     line endings normalised), seeds, n, PCA rule, MD5 and SHA-256 of the PE file,
     SHA-256 of the input tables, package versions. A cache entry is used only if all of
     these match the current run; otherwise it is rebuilt.
  Leave-one-out (catalog stage) uses its own grid, built from the smallest sigma of the
  remaining events, whenever the left-out event set the common grid.
  Off-source stage: the only stage that reads strain; it reads off-source segments only,
  refusing any segment that overlaps an event of the event list (module1_offsource.py).
  3. Every output is marked "official" only if the code hash equals the hash in the
     freeze manifest (default FREEZE_v1.0.json, written at the v1.0 tag); otherwise
     "trial".

  python module1_run.py build   --events GW150914,GW230627_015337 --labels both
  python module1_run.py checks  --events GW150914,GW230627_015337 --labels both
  python module1_run.py catalog --events GW150914,GW230627_015337 --labels xphm
  python module1_run.py build   --events all-primary --labels both      (official run)
"""
import argparse
import csv
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# ---- safeguard 1: on-source strain cannot be read (installed before module1_event) ----
import residual_check as _rc


def _forbidden(*a, **k):
    raise RuntimeError("module1_run: reading detector strain is forbidden before the v1.0 tag "
                       "(MODULE1_DESIGN §12); this driver works on synthetic noise only")


_ORIG_READ_STRAIN = _rc.read_strain        # kept only for the off-source stage
_ORIG_SEGMENT = _rc.segment
_rc.read_strain = _forbidden
_rc.segment = _forbidden
_rc.check = _forbidden

import h5py                                  # noqa: E402
import module1_calibration as cal            # noqa: E402
import module1_catalog as mc                 # noqa: E402
import module1_event as me                   # noqa: E402
import module1_offsource as mo               # noqa: E402

CODE_FILES = ["module1_run.py", "module1_event.py", "module1_catalog.py", "module1_offsource.py", "module1_calibration.py",
              "module1_cosmo.py", "f3_derivatives.py", "residual_check.py"]
CACHE_VERSION = 1

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


# ------------------------------------------------------------------ hashes ----------
def sha256_file(path, normalise_eol=False):
    h = hashlib.sha256()
    if normalise_eol:
        h.update(open(path, "rb").read().replace(b"\r\n", b"\n"))
    else:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 22), b""):
                h.update(chunk)
    return h.hexdigest()


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def code_hash():
    h = hashlib.sha256()
    per = {}
    for name in CODE_FILES:
        d = sha256_file(os.path.join(HERE, name), normalise_eol=True)
        per[name] = d
        h.update(name.encode() + b"\0" + d.encode() + b"\n")
    return h.hexdigest(), per


def versions():
    out = {"python": platform.python_version(), "numpy": np.__version__}
    for mod in ("lal", "lalsimulation", "bilby", "h5py", "scipy"):
        try:
            out[mod] = __import__(mod).__version__
        except Exception:
            out[mod] = "?"
    return out


def official_status(chash, manifest, per=None):
    """A run is official only if the combined code hash AND every per-file hash match the
    freeze manifest; a per-file mismatch is named, so a silent edit cannot pass."""
    if not os.path.exists(manifest):
        return "trial", "no freeze manifest"
    m = json.load(open(manifest, encoding="utf-8"))
    tag = m.get("tag", "?")
    if m.get("code_sha256") != chash:
        return "trial", f"code hash differs from {tag}"
    want = m.get("files", {})
    if per and want:
        bad = sorted(k for k in set(want) | set(per) if want.get(k) != per.get(k))
        if bad:
            return "trial", f"per-file hash differs from {tag}: {', '.join(bad)}"
    return "official", f"code hash matches {tag}"


# ------------------------------------------------------------------ selection -------
def select(args):
    rows5 = {r["commonName"]: r for r in csv.DictReader(open(args.f5, encoding="utf-8"))}
    if args.events != "all-primary":
        evs = args.events.split(",")
    else:
        bad = set()
        for r in csv.DictReader(open(args.residual_check, encoding="utf-8")):
            if r.get("error") or not (0.8 <= float(r["snr_ratio"]) <= 1.2):
                bad.add(r["event"])
        evs = [e for e, r in rows5.items() if r["status"] == "BOTH_OK" and e not in bad]
        n_s = sum(rows5[e]["xphm_label"].startswith("C01") for e in evs)
        print(f"primary test: {n_s} search (C01) + {len(evs) - n_s} confirmation (C00) events; "
              f"excluded by the reference-fit gate: {sorted(bad)}")
    pairs = []
    for e in evs:
        r5 = rows5[e]
        labs = []
        if args.labels in ("both", "xphm"):
            labs.append(r5["xphm_label"])
        if args.labels in ("both", "eob") and r5.get("eob_label"):
            labs.append(r5["eob_label"])
        for lab in labs:
            pairs.append((e, lab, r5))
    return pairs


def cache_paths(args, ev, label):
    base = os.path.join(args.cache, f"{ev}__{label.replace(':', '_')}")
    return base + ".npz", base + ".json"


def expected_meta(args, ev, label, r5, chash, pe_path):
    return {"cache_version": CACHE_VERSION, "event": ev, "label": label, "code_sha256": chash,
            "seed": me.SEED, "n": args.n, "frac": me.FRAC, "n_leak": me.N_LEAK,
            "pe_file": os.path.basename(pe_path), "pe_md5": md5_file(pe_path),
            "pe_sha256": sha256_file(pe_path), "f5_sha256": sha256_file(args.f5)}


KEY_FIELDS = ["cache_version", "code_sha256", "seed", "n", "frac", "n_leak", "pe_md5", "pe_sha256", "f5_sha256"]


def cache_valid(meta_path, exp):
    if not os.path.exists(meta_path):
        return False, "absent"
    m = json.load(open(meta_path, encoding="utf-8"))
    diff = [k for k in KEY_FIELDS if m.get(k) != exp[k]]
    return (not diff), ("ok" if not diff else "differs in " + ",".join(diff))


# ------------------------------------------------------------------ build -----------
def stage_build(args, chash, status):
    for ev, label, r5 in select(args):
        pe_path = os.path.join(args.pe_dir, r5["file"].split(" ")[0])
        npz, js = cache_paths(args, ev, label)
        exp = expected_meta(args, ev, label, r5, chash, pe_path)
        ok, why = cache_valid(js, exp)
        if ok and os.path.exists(npz):
            print(f"\n=== {ev} {label}: cache valid, skipped", flush=True)
            continue
        print(f"\n=== {ev} {label}: building (cache {why})", flush=True)
        t0 = time.time()
        ifos = [i for i in r5["ifos"].split(",") if i]
        with h5py.File(pe_path, "r") as fh:
            est, ks, info, sb, href, net, fb = me.build(ev, label, r5["xphm_label"], fh, ifos, args.n,
                                                        me.FRAC, lambda m: print(m, flush=True))
        arrays = {"QA": est["A"].Q.astype(np.float32), "QB": est["B"].Q.astype(np.float32),
                  "J": info["J"].astype(np.float32), "gr_leak": info["gr_leak"].astype(np.float32),
                  "fb": fb, "K": np.array(info["K"])}
        for i in ifos:
            arrays[f"href_{i}"] = href[i].astype(np.complex64)
            arrays[f"sw_{i}"] = net.sw[i]
            for k in ("freqs", "mu_a", "sg_a", "mu_p", "sg_p"):
                arrays[f"pri_{i}_{k}"] = sb[i].pri[k]
        os.makedirs(args.cache, exist_ok=True)
        tmp = npz + ".tmp.npz"
        np.savez(tmp, **arrays)
        os.replace(tmp, npz)
        meta = dict(exp)
        meta.update({"ifos": ifos, "k_A": ks["A"], "k_B": ks["B"], "sigma_lin": est["A"].sigma(3),
                     "s3": est["A"].s[3], "z": info["z"], "z_med": info["z_med"], "K": info["K"],
                     "f_high": info["f_high"], "cal_source": info["cal_source"], "fails": info["fails"],
                     "nsamp": info["nsamp"], "snr_ref": info["snr_ref"], "npz_sha256": sha256_file(npz),
                     "versions": versions(), "status": status, "seconds": round(time.time() - t0),
                     "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
        json.dump(meta, open(js, "w", encoding="utf-8"), indent=1)
        print(f"  cached: k {ks['A']}/{ks['B']}, sigma_lin {meta['sigma_lin']:.3e} m^2, "
              f"{os.path.getsize(npz) / 1e6:.0f} MB, {meta['seconds']} s", flush=True)


# ------------------------------------------------------------------ load ------------
def load(args, ev, label, r5, chash):
    npz, js = cache_paths(args, ev, label)
    pe_path = os.path.join(args.pe_dir, r5["file"].split(" ")[0])
    ok, why = cache_valid(js, expected_meta(args, ev, label, r5, chash, pe_path))
    if not ok or not os.path.exists(npz):
        raise RuntimeError(f"{ev} {label}: cache not usable ({why}); run the build stage first")
    meta = json.load(open(js, encoding="utf-8"))
    if sha256_file(npz) != meta["npz_sha256"]:
        raise RuntimeError(f"{ev} {label}: cache file altered after it was written")
    d = np.load(npz)
    ifos = meta["ifos"]
    fb = d["fb"]
    df = fb[1] - fb[0]
    Seff = {i: 4 * df / d[f"sw_{i}"].astype(np.float64) ** 2 for i in ifos}
    net = me.Network(fb, np.ones(len(fb), bool), ifos, Seff)
    href = {i: d[f"href_{i}"].astype(np.complex128) for i in ifos}
    sb = {i: cal.SplineBasis({k: d[f"pri_{i}_{k}"] for k in ("freqs", "mu_a", "sg_a", "mu_p", "sg_p")})
          for i in ifos}
    K = float(d["K"])
    T = {p: net.stack({i: 1j * K * fb ** p * href[i] for i in ifos}) for p in (2, 3, 4)}
    J = d["J"].astype(np.float64)
    est = {}
    for name in ("A", "B"):
        Q, _ = np.linalg.qr(d["Q" + name].astype(np.float64))
        est[name] = me.Estimator(Q, J, T)
    info = {"K": K, "f_high": float(fb[-1]), "gr_leak": d["gr_leak"], "z": meta["z"], "z_med": meta["z_med"]}
    return {"name": ev, "label": label, "A": est["A"], "B": est["B"], "est": est, "K": K, "href": href,
            "net": net, "fb": fb, "sb": sb, "T": T, "info": info, "meta": meta,
            "sigma_lin": est["A"].sigma(3)}


# ------------------------------------------------------------------ checks ----------
def stage_checks(args, chash, status):
    out = args.out or "module1_checks.csv"
    new = not os.path.exists(out)
    fout = open(out, "a", newline="", encoding="utf-8")
    keys = ["event", "label", "status", "code_sha256"] + me.KEYS
    w = csv.DictWriter(fout, fieldnames=keys)
    if new:
        w.writeheader()
    for ev, label, r5 in select(args):
        t0 = time.time()
        e = load(args, ev, label, r5, chash)
        A, B = e["A"], e["B"]
        res = me.run_checks(e["est"], e["info"], e["sb"], e["href"], e["net"], e["fb"], e["T"], args.trials,
                            np.random.default_rng(me.SEED + 1))
        res.update(me.scan_checks(e["est"], e["info"], e["sb"], e["href"], e["net"], e["fb"], args.scan_trials,
                                  np.random.default_rng(me.SEED + 2)))
        m = e["meta"]
        cosAB = float(A.Tp[3] @ B.Tp[3] / (np.linalg.norm(A.Tp[3]) * np.linalg.norm(B.Tp[3])))
        row = {"event": ev, "label": label, "status": status, "code_sha256": chash, "ifos": ",".join(m["ifos"]),
               "k_A": m["k_A"], "k_B": m["k_B"], "s2": A.s[2], "s3": A.s[3], "s4": A.s[4],
               "sigma3_m2": A.sigma(3), "z": m["z"], "z_med": m["z_med"], "cos_AB": cosAB,
               "cal_source": m["cal_source"], "fails": m["fails"], **res, "seconds": round(time.time() - t0),
               "error": ""}
        w.writerow({k: row.get(k, "") for k in keys})
        fout.flush()
        print(f"{ev:18s} {label:30s} s3 {A.s[3]:.3f} sigma_lin {A.sigma(3):.2e}  null z {res['null_mean_z']:+.2f}/"
              f"{res['null_std_z']:.2f}  scan std {res['scan_+0sig_std']:.2f}  Wilks cov "
              f"{res['scan_+0sig_cover90']:.2f}  leak std {res.get('leak_std', float('nan')):.2f}  "
              f"({row['seconds']} s)", flush=True)
    fout.close()


# ------------------------------------------------------------------ catalog ---------
def combine(parts, grid):
    """parts: list over subspaces of (list of Gram matrices, of belt V, of trial V), one
    entry per event. The reported set is the UNION over subspaces (v0.9.14, §7); the same
    noise realisations are used in both, so the sets of a trial come from the same data."""
    Deltas, cs, Ves = [], [], []
    for Gs, Vbs, Ves_ in parts:
        D = mc.delta_matrix(0.5 * (sum(Gs) + sum(Gs).T))
        Deltas.append(D)
        cs.append(mc.neyman(sum(Vbs), D))
        Ves.append(sum(Ves_))
    t0 = int(np.argmin(np.abs(grid)))
    n = Ves[0].shape[0]
    lo, hi, cont, edge = [], [], 0, 0
    for j in range(n):
        qs = [mc.q_obs_all(V[j], D, t0) for V, D in zip(Ves, Deltas)]
        a, b, ok = mc.union_interval(qs, cs, grid)
        lo.append(a)
        hi.append(b)
        cont += ok
        edge += (a <= grid[0]) or (b >= grid[-1])
    cov = {}
    for frac in (-0.5, -0.2, 0.0, 0.2, 0.5):
        t = int(np.argmin(np.abs(grid - frac * grid[-1])))
        covered = np.zeros(n, bool)
        wilks = np.zeros(n, bool)
        for V, D, c in zip(Ves, Deltas, cs):
            q = np.array([mc.q_obs_all(V[j], D, t)[t] for j in range(n)])
            covered |= q <= c[t]
            wilks |= q <= 2.71
        cov[f"{grid[t]:+.3e}"] = {"neyman": float(covered.mean()), "wilks": float(wilks.mean())}
    return {"exp_lower_median": float(np.nanmedian(lo)), "exp_upper_median": float(np.nanmedian(hi)),
            "contiguous": cont / n, "edge": edge / n, "coverage": cov,
            "subspaces": len(parts),
            "critical_values": [float(min(c.min() for c in cs)), float(max(c.max() for c in cs))]}


def stage_catalog(args, chash, status):
    pairs = select(args)
    metas = []
    for ev, label, r5 in pairs:
        _, js = cache_paths(args, ev, label)
        metas.append(json.load(open(js, encoding="utf-8")))
    smin = min(m["sigma_lin"] for m in metas)
    grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * smin
    w = np.array([1 / m["sigma_lin"] ** 2 for m in metas])
    neff = w.sum() ** 2 / (w ** 2).sum()
    print(f"catalog: {len(pairs)} pairs, grid ±{args.half_width:g} x {smin:.3e} m^2 ({len(grid)} points), "
          f"N_eff {neff:.2f}", flush=True)
    subs = ("A", "B")
    acc = {w: ([], [], []) for w in subs}
    F = np.zeros((3, 3))
    for k, (ev, label, r5) in enumerate(pairs):
        t0 = time.time()
        e = load(args, ev, label, r5, chash)
        A = e["A"]
        for a, p in enumerate((2, 3, 4)):
            for b, q in enumerate((2, 3, 4)):
                F[a, b] += float(A.Tp[p] @ A.CTp[q])
        for wsp in subs:
            G, vb, ve = mc.event_trials(e, grid, args.belt_trials, args.exp_trials,
                                        me.SEED + 7 + 1000 * k, which=wsp)
            for store, val in zip(acc[wsp], (G, vb, ve)):
                store.append(val)
        del e
        print(f"  {ev} {label}: {time.time() - t0:.0f} s (subspaces {', '.join(subs)})", flush=True)
    parts = [acc[w] for w in subs]
    F = 0.5 * (F + F.T)
    R = F / np.sqrt(np.outer(np.diag(F), np.diag(F)))
    wv = np.linalg.eigvalsh(R)
    res = {"status": status, "code_sha256": chash, "pairs": [(m["event"], m["label"], m["sigma_lin"]) for m in metas],
           "N_eff": neff, "grid": [float(grid[0]), float(grid[-1]), float(grid[1] - grid[0])],
           "joint_condition": float(wv[-1] / wv[0]), "joint_inflation_f3": float(np.sqrt(np.linalg.inv(R)[1, 1])),
           "belt_trials": args.belt_trials, "exp_trials": args.exp_trials}
    res["all"] = combine(parts, grid)
    top = int(np.argmax(w))
    if len(pairs) > 1:
        keep = [i for i in range(len(pairs)) if i != top]
        smin_loo = min(metas[i]["sigma_lin"] for i in keep)
        if smin_loo <= smin * (1 + 1e-9):
            loo = combine([tuple([x[i] for i in keep] for x in part) for part in parts], grid)
            loo_grid = grid
        else:
            # own grid from the smallest sigma of the remaining events, same seeds per event
            loo_grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * smin_loo
            print(f"\nleave-one-out without {metas[top]['event']}: own grid ±{args.half_width:g} x "
                  f"{smin_loo:.3e} m^2", flush=True)
            acc2 = {w: ([], [], []) for w in subs}
            for k in keep:
                ev, label, r5 = pairs[k]
                e = load(args, ev, label, r5, chash)
                for wsp in subs:
                    G, vb, ve = mc.event_trials(e, loo_grid, args.belt_trials, args.exp_trials,
                                                me.SEED + 7 + 1000 * k, which=wsp)
                    for store, val in zip(acc2[wsp], (G, vb, ve)):
                        store.append(val)
                del e
            loo = combine([acc2[w] for w in subs], loo_grid)
        res["leave_one_out"] = {"left_out": metas[top]["event"],
                                "grid": [float(loo_grid[0]), float(loo_grid[-1]), float(loo_grid[1] - loo_grid[0])],
                                **loo}
    for key in ("all", "leave_one_out"):
        if key not in res:
            continue
        r = res[key]
        print(f"\n[{key}{' without ' + r['left_out'] if key == 'leave_one_out' else ''}] expected 90% interval "
              f"[{r['exp_lower_median']:+.2e}, {r['exp_upper_median']:+.2e}] m^2; contiguous {r['contiguous']:.2f}; "
              f"edge {r['edge']:.3f}; critical values {r['critical_values'][0]:.2f}–{r['critical_values'][1]:.2f}")
        if r.get("subspaces", 1) > 1:
            print(f"   union of {r['subspaces']} subspaces (v0.9.14): coverage above 0.90 is expected "
                  f"and reported as conservative; only below the interval stops the analysis")
        if r["edge"] > 0.01:
            print("   WARNING: more than 1% of the intervals touch the grid edge (v0.9.11 rule): "
                  "widen the grid before any real estimate")
        for kk, v in r["coverage"].items():
            print(f"   true {kk}: Neyman {v['neyman']:.3f}  Wilks {v['wilks']:.3f}")
    print(f"\njoint f2/f3/f4 condition {res['joint_condition']:.0f}, f3 inflation x{res['joint_inflation_f3']:.1f}; "
          f"LVK GWTC-4.0 [{mc.LVK[0]:+.1e}, {mc.LVK[1]:+.1e}]; status {status}")
    out = args.out or f"module1_catalog_{args.labels}.json"
    json.dump(res, open(out, "w"), indent=1)
    print(f"written {out}")



# ------------------------------------------------------------------ off-source ------
def stage_offsource(args, chash, status):
    from scipy.signal import resample_poly
    from scipy.signal.windows import tukey
    import f3_derivatives as f3
    out = args.out or f"module1_offsource_{args.labels}.csv"
    keys = ["event", "label", "status", "code_sha256", "seg_start", "inj", "q0", "c0", "covered", "lambda_hat",
            "dt_ms", "abs_a", "c_index"]
    new = not os.path.exists(out)
    fout = open(out, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fout, fieldnames=keys)
    if new:
        w.writeheader()
    manifest = {(r["event"], r["ifo"]): r for r in
                csv.DictReader(open(os.path.join(args.strain_dir, "strain_manifest.csv"), encoding="utf-8"))}
    summary = {"C": [0, 0], "ref": [0, 0]}
    for k, (ev, label, r5) in enumerate(select(args)):
        t0 = time.time()
        e = load(args, ev, label, r5, chash)
        meta = e["meta"]
        ifos, fb = meta["ifos"], e["fb"]
        pe_path = os.path.join(args.pe_dir, r5["file"].split(" ")[0])
        with h5py.File(pe_path, "r") as fh:
            sample, cfg, _ = f3.read_label(fh, label)
            st = f3.settings(cfg)
            tab = me.posterior_table(fh, label)
        dur, srate = st["duration"], st["srate"]
        seg_on = sample["geocent_time"] + me.rc.POST_TRIGGER - dur
        f, hp, hc = me.rc.polarizations(label, sample, st)
        band = (f >= st["f_an"]) & (f <= st["f_high"])
        if not np.allclose(f[band], fb):
            raise RuntimeError(f"{ev} {label}: frequency band differs from the cache")
        calf = {i: e["sb"][i].mean_factor(fb) for i in ifos}
        # set C: indices after those consumed by sets A and B (2n attempts + failures)
        idx = np.random.default_rng(me.SEED).permutation(len(tab["log_likelihood"]))
        pos = 2 * meta["n"] + meta["fails"]
        hC = []
        while len(hC) < mo.N_C and pos < len(idx):
            j = int(idx[pos])
            pos += 1
            sC = {kk: float(tab[kk][j]) for kk in me.NEED}
            try:
                fp, sp, sc_ = me.rc.polarizations(label, sC, st)
                hC.append((j, {i: me.rc.project(f, sp, sc_, i, sC, seg_on)[band] * calf[i] for i in ifos}))
            except Exception:
                continue
        if len(hC) < mo.N_C:
            print(f"  {ev} {label}: only {len(hC)} set-C waveforms", flush=True)
        windows = mo.forbidden_windows(args.events_list, dur)
        win = tukey(int(round(dur * srate)), alpha=2 * me.rc.ROLL_OFF / dur)
        strain = {}
        for i in ifos:
            x, dt, t0f = _ORIG_READ_STRAIN(os.path.join(args.strain_dir, manifest[(ev, i)]["file"]))
            if abs(1 / dt - srate) > 1e-6:
                x = resample_poly(x, int(srate), int(round(1 / dt)))
                dt = 1.0 / srate
            strain[i] = (x, dt, t0f)
        segs = []
        for kk in range(1, 4 * mo.N_OFF):
            for side in (-1, +1):
                if len(segs) >= mo.N_OFF:
                    break
                s0 = seg_on + side * kk * (dur + mo.OFF_GAP)
                ok, _ = mo.segment_allowed(s0, dur, windows)
                if not ok:
                    continue
                try:
                    dseg = {}
                    for i in ifos:
                        x, dt, t0f = strain[i]
                        dseg[i] = (np.fft.rfft(_ORIG_SEGMENT(x, dt, t0f, s0, dur) * win) * dt)[band]
                    segs.append((s0, dseg))
                except ValueError:
                    pass
        del strain
        own_grid = np.round(np.arange(-args.half_width, args.half_width + 1e-9, args.step), 6) * e["sigma_lin"]
        c0, scan, tz = mo.event_c0(e, own_grid, args.belt_trials, me.SEED + 11 + 1000 * k)
        Seff = {i: 4 * (fb[1] - fb[0]) / e["net"].sw[i] ** 2 for i in ifos}
        for m, (s0, dseg) in enumerate(segs):
            for inj in ("C", "ref"):
                if inj == "C":
                    if m >= len(hC):
                        continue
                    cj, hinj = hC[m]
                else:
                    cj, hinj = -1, e["href"]
                d = {i: dseg[i] + hinj[i] for i in ifos}
                r, tsh, a = mo.refine(d, e["href"], Seff, fb[1] - fb[0], fb=fb)
                q0, g = mo.q0_of_residual(e["net"].stack(r), e, scan, tz)
                cov = q0 <= c0
                summary[inj][0] += cov
                summary[inj][1] += 1
                w.writerow({"event": ev, "label": label, "status": status, "code_sha256": chash, "seg_start": s0,
                            "inj": inj, "q0": q0, "c0": c0, "covered": cov, "lambda_hat": own_grid[g],
                            "dt_ms": tsh * 1e3, "abs_a": abs(a), "c_index": cj})
        fout.flush()
        print(f"{ev:18s} {label:30s} segments {len(segs)}  c0 {c0:.2f}  ({time.time() - t0:.0f} s)", flush=True)
    fout.close()
    for inj in ("C", "ref"):
        nc, n = summary[inj]
        if n:
            p = nc / n
            lo, hi = mo.binomial_interval(0.90, n)
            verdict = ("STOP: below the 99% binomial interval around 0.90" if p < lo else
                       "conservative (above the interval)" if p > hi else "within the interval")
            print(f"{inj:4s} coverage {nc}/{n} = {p:.3f}; 99% interval around 0.90: {lo:.3f}–{hi:.3f} -> "
                  f"{verdict if inj == 'C' else 'diagnostic'}")

# ------------------------------------------------------------------ main ------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["build", "checks", "catalog", "offsource", "hash"])
    ap.add_argument("--events", default="GW150914,GW230627_015337")
    ap.add_argument("--labels", default="both", choices=["both", "xphm", "eob"])
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    ap.add_argument("--residual-check", default="residual_check.csv")
    ap.add_argument("--cache", default=os.path.expanduser("~/gw/cache_module1"))
    ap.add_argument("--n", type=int, default=800)
    ap.add_argument("--trials", type=int, default=200)
    ap.add_argument("--scan-trials", type=int, default=100)
    ap.add_argument("--belt-trials", type=int, default=5000)
    ap.add_argument("--exp-trials", type=int, default=2000)
    ap.add_argument("--half-width", type=float, default=8.0)
    ap.add_argument("--step", type=float, default=0.05)
    ap.add_argument("--strain-dir", default=os.path.expanduser("~/gwdata/strain"))
    ap.add_argument("--events-list", default="event_selection/event_list_v1.csv")
    ap.add_argument("--freeze-manifest", default=os.path.join(HERE, "FREEZE_v1.0.json"))
    ap.add_argument("--write-manifest", default="", metavar="TAG",
                    help="with stage hash: write the freeze manifest for TAG (refuses to overwrite)")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    chash, per = code_hash()
    status, why = official_status(chash, args.freeze_manifest, per)
    print(f"code SHA-256 {chash}\nstatus: {status} ({why})", flush=True)
    if args.stage == "hash":
        for k, v in per.items():
            print(f"  {v}  {k}")
        if args.write_manifest:
            if os.path.exists(args.freeze_manifest):
                raise SystemExit(f"{args.freeze_manifest} exists; refusing to overwrite a freeze manifest")
            json.dump({"tag": args.write_manifest, "code_sha256": chash, "files": per,
                       "code_files": CODE_FILES, "cache_version": CACHE_VERSION, "seed": me.SEED,
                       "frac": me.FRAC, "n_leak": me.N_LEAK, "versions": versions(),
                       "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
                      open(args.freeze_manifest, "w", encoding="utf-8"), indent=1)
            print(f"\nwritten {args.freeze_manifest} for tag {args.write_manifest}")
        return
    if args.stage == "offsource" and args.belt_trials == 5000:
        args.belt_trials = 2000          # event-only belt at Lambda = 0
    {"build": stage_build, "checks": stage_checks, "catalog": stage_catalog,
     "offsource": stage_offsource}[args.stage](args, chash, status)


if __name__ == "__main__":
    main()
