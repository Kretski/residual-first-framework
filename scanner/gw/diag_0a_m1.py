#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diag_0a_m1.py — Measurement 1(ii) of DIAG_0a_SPEC.md (pre-specified diagnostic of the
rule 0a failure of Module 1 v1.0). It cannot change the v1.0 result.

For every pair of the specification (four failing SEOBNRv4PHM search events, four
controls with the smallest |dt_ms| in residual_check.csv, and the IMRPhenomXPHM labels
of the four failing events): 200 set-B posterior samples are generated and projected
exactly as the set-C injections of the off-source stage (own sky position, polarization
and time; calibration mean factor), and each is aligned to the reference waveform with
the same refine() of module1_offsource (common time shift ±100 ms and complex
amplitude), noise-free. Reported per pair: median and standard deviation of the time
shift, fraction with |shift| > 10 ms, median |a|, and the v1.0 off-source coverage.

No strain is read and no coverage is computed. The cache is loaded with module1_run's
own loader and must be valid for the frozen code.

  python diag_0a_m1.py --out diag_0a_m1.csv
"""
import argparse
import csv
import os
import time

import h5py
import numpy as np

import f3_derivatives as f3
import module1_event as me
import module1_offsource as mo
import module1_run as run

FAILING = ["GW150914", "GW170814", "GW190706_222641", "GW200129_065458"]
N_B = 200


def driver_args():
    """The module1_run defaults used to build and validate the cache."""
    return argparse.Namespace(
        pe_dir=os.path.expanduser("~/gwdata/pe"), f5="f5/f5_results.csv",
        residual_check="residual_check.csv", cache=os.path.expanduser("~/gw/cache_module1"),
        n=800, strain_dir=os.path.expanduser("~/gwdata/strain"),
        events_list="event_selection/event_list_v1.csv")


def controls(args):
    prim = [e for e in open("primary_search_eob.txt").read().split(",") if e]
    rows = [r for r in csv.DictReader(open(args.residual_check, encoding="utf-8"))
            if "SEOBNRv4PHM" in r["label"] and r["event"] in prim and r["event"] not in FAILING]
    rows.sort(key=lambda r: abs(float(r["dt_ms"])))
    return [r["event"] for r in rows[:4]]


def coverage(ev, label):
    fn = "primary_offsource_search_eob.csv" if "SEOBNR" in label else "primary_offsource_search_xphm.csv"
    rows = [r for r in csv.DictReader(open(fn)) if r["event"] == ev and r["label"] == label and r["inj"] == "C"]
    return sum(r["covered"] in ("True", "1") for r in rows), len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="diag_0a_m1.csv")
    ap.add_argument("--n-b", type=int, default=N_B)
    cli = ap.parse_args()
    args = driver_args()
    chash, _ = run.code_hash()
    rows5 = {r["commonName"]: r for r in csv.DictReader(open(args.f5, encoding="utf-8"))}
    ctrl = controls(args)
    pairs = [(e, rows5[e]["eob_label"], "failing") for e in FAILING] + \
            [(e, rows5[e]["eob_label"], "control") for e in ctrl] + \
            [(e, rows5[e]["xphm_label"], "failing-XPHM") for e in FAILING]
    print(f"code SHA-256 {chash}\ncontrols (smallest |dt_ms|, SEOBNRv4PHM): {ctrl}", flush=True)

    fo = open(cli.out, "w", newline="", encoding="utf-8")
    w = csv.DictWriter(fo, fieldnames=["event", "label", "role", "sample_index", "shift_ms", "abs_a"])
    w.writeheader()
    summary = []
    for ev, label, role in pairs:
        t0 = time.time()
        r5 = rows5[ev]
        e = run.load(args, ev, label, r5, chash)
        meta, ifos, fb = e["meta"], e["meta"]["ifos"], e["fb"]
        with h5py.File(os.path.join(args.pe_dir, r5["file"].split(" ")[0]), "r") as fh:
            sample, cfg, _ = f3.read_label(fh, label)
            st = f3.settings(cfg)
            tab = me.posterior_table(fh, label)
        seg_on = sample["geocent_time"] + me.rc.POST_TRIGGER - st["duration"]
        f, _, _ = me.rc.polarizations(label, sample, st)
        band = (f >= st["f_an"]) & (f <= st["f_high"])
        if not np.allclose(f[band], fb):
            raise RuntimeError(f"{ev} {label}: frequency band differs from the cache")
        calf = {i: e["sb"][i].mean_factor(fb) for i in ifos}
        Seff = {i: 4 * (fb[1] - fb[0]) / e["net"].sw[i] ** 2 for i in ifos}
        idx = np.random.default_rng(me.SEED).permutation(len(tab["log_likelihood"]))
        pos = meta["n"] + (0 if meta.get("fails", 0) == 0 else meta["fails"])   # set B starts after set A
        shifts, amps = [], []
        while len(shifts) < cli.n_b and pos < len(idx):
            j = int(idx[pos])
            pos += 1
            s = {k: float(tab[k][j]) for k in me.NEED}
            try:
                _, sp, sc = me.rc.polarizations(label, s, st)
                hs = {i: me.rc.project(f, sp, sc, i, s, seg_on)[band] * calf[i] for i in ifos}
            except Exception:
                continue
            _, t, a = mo.refine(hs, e["href"], Seff, fb[1] - fb[0], fb=fb)
            shifts.append(t * 1e3)
            amps.append(abs(a))
            w.writerow({"event": ev, "label": label, "role": role, "sample_index": j,
                        "shift_ms": t * 1e3, "abs_a": abs(a)})
        fo.flush()
        sh = np.array(shifts)
        nc, n = coverage(ev, label)
        summary.append((ev, label, role, np.median(sh), np.std(sh), np.mean(np.abs(sh) > 10), np.median(amps), nc, n))
        print(f"  {ev:18s} {label:22s} {role:13s} shift median {np.median(sh):+7.2f} ms  sd {np.std(sh):6.2f}  "
              f"|shift|>10ms {np.mean(np.abs(sh) > 10):.2f}  |a| {np.median(amps):.3f}  coverage {nc}/{n}  "
              f"({len(sh)} samples, {time.time() - t0:.0f} s)", flush=True)
    fo.close()
    print("\nMeasurement 1(ii) summary (reported whatever the outcome):")
    for s in summary:
        print(f"  {s[0]:18s} {s[2]:13s} median {s[3]:+7.2f} ms  sd {s[4]:6.2f}  frac>10ms {s[5]:.2f}  "
              f"|a| {s[6]:.3f}  coverage {s[7]}/{s[8]}")


if __name__ == "__main__":
    main()
