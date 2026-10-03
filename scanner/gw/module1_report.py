#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
module1_report.py — reporting layer over the FROZEN Module 1 code (tag v1.0-module1).
Applies the decision rules of MODULE1_DESIGN §11 to the official outputs and says which
event-model pairs leave the primary test and whether any stopping rule has fired.

This layer does NOT recompute anything of its own and does not read strain. The estimates
it needs for the Kolmogorov-Smirnov test of rule 2a (which the frozen CSV summarises but
does not store) are produced by calling the FROZEN functions of module1_event through the
module1_run cache, with the frozen seeds — the same code path that produced the CSV.

Two guarantees, both enforced at run time:
  1. the frozen modules must match FREEZE_v1.0.json (combined and per-file hashes);
  2. the summaries recomputed from the samples (leak_n, leak_mean, leak_std,
     leak_frac_abs_ge1, scan_+0sig_std, noise_frac_abs_ge1) must equal those in the frozen
     CSV to 1e-12; otherwise the report stops. This proves the KS samples are the ones the
     recorded numbers came from.

  python module1_report.py --checks official_checks.csv \
      --catalog official_catalog_search_xphm.json,official_catalog_confirm_xphm.json \
      --offsource official_offsource.csv --cohort search --labels xphm
"""
import argparse
import csv
import json
import os
import sys

import numpy as np
from scipy.stats import ks_2samp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import module1_run as mr          # noqa: E402  (installs the strain guard on import)
import module1_event as me        # noqa: E402

TOL = 1e-12
RULES = {
    "1": "null bias |mean Lambda| > 0.3 sigma_lin",
    "1a": "subspace systematic median |Lambda_A - Lambda_B| > 0.5 sigma_lin",
    "2a": "GR leakage: two-sided KS p < 0.01",
}
KS_ALPHA = 0.01
NULL_BIAS = 0.3
AB_MEDIAN = 0.5

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def leakage_samples(e, trials, half_width=8.0, step=0.05):
    """Scan estimates for the GR-leakage and the pure-noise suites, produced by the frozen
    code path: the same generators (gen("leak"), gen("level+0")), the same grid, the same
    ProfileScan. Returns (leak, noise) in units of sigma_lin."""
    A = e["A"]
    K, href, net, fb, sb = e["K"], e["href"], e["net"], e["fb"], e["sb"]
    units = np.round(np.arange(-half_width, half_width + 1e-9, step), 6)
    sc = me.ProfileScan(A, K, href, net, fb, units * A.sigma(3))
    nvec = len(A.Q)

    def gen(tag):
        return np.random.default_rng([me.SEED, int(np.frombuffer(tag.encode().ljust(8, b"\0")[:8],
                                                                 dtype=np.uint64)[0] % (2 ** 32))])

    def cal_error(g):
        out = {}
        for i in net.ifos:
            pr = sb[i].pri
            out[i] = (sb[i].factor(fb, g.normal(0, pr["sg_a"]), g.normal(0, pr["sg_p"])) - 1.0) * href[i]
        return net.stack(out)

    t0 = int(np.argmin(np.abs(units)))
    noise = []
    g = gen("level+0")
    for _ in range(trials):
        n = g.normal(size=nvec) + cal_error(g)
        noise.append(units[int(np.argmin(sc.chi2(n - A.Q @ (A.Q.T @ n), t0)))])
    leak = []
    L = e["info"]["gr_leak"]
    g = gen("leak")
    for j in range(L.shape[1]):
        for _ in range(max(1, trials // L.shape[1])):
            n = L[:, j].astype(np.float64) + g.normal(size=nvec) + cal_error(g)
            leak.append(units[int(np.argmin(sc.chi2(n - A.Q @ (A.Q.T @ n), t0)))])
    return np.array(leak), np.array(noise)


def verify(row, leak, noise):
    """Recompute the recorded summaries from the samples; any mismatch stops the report."""
    got = {"leak_n": float(len(leak)), "leak_mean": float(leak.mean()), "leak_std": float(leak.std()),
           "leak_frac_abs_ge1": float(np.mean(np.abs(leak) >= 1.0)),
           "scan_+0sig_std": float(noise.std()),
           "noise_frac_abs_ge1": float(np.mean(np.abs(noise) >= 1.0))}
    bad = []
    for k, v in got.items():
        rec = row.get(k, "")
        if rec in ("", None):
            continue
        if abs(float(rec) - v) > TOL * max(1.0, abs(v)):
            bad.append(f"{k}: CSV {float(rec):.12g} vs recomputed {v:.12g}")
    return bad


def binomial(p, n, z=2.576):
    s = np.sqrt(p * (1 - p) / n)
    return p - z * s, p + z * s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--checks", default="official_checks.csv")
    ap.add_argument("--catalog", default="", help="comma-separated catalog JSON files")
    ap.add_argument("--offsource", default="")
    ap.add_argument("--cohort", default="search", choices=["search", "confirm"])
    ap.add_argument("--labels", default="xphm", choices=["xphm", "eob", "both"])
    ap.add_argument("--pe-dir", default=os.path.expanduser("~/gwdata/pe"))
    ap.add_argument("--f5", default="f5/f5_results.csv")
    ap.add_argument("--cache", default=os.path.expanduser("~/gw/cache_module1"))
    ap.add_argument("--n", type=int, default=800)
    ap.add_argument("--scan-trials", type=int, default=100)
    ap.add_argument("--freeze-manifest", default=os.path.join(HERE, "FREEZE_v1.0.json"))
    ap.add_argument("--out", default="module1_report.json")
    args = ap.parse_args()

    chash, per = mr.code_hash()
    status, why = mr.official_status(chash, args.freeze_manifest, per)
    print(f"frozen code SHA-256 {chash}\nstatus: {status} ({why})")
    if status != "official":
        raise SystemExit("the report runs only against the frozen code")

    events = [e for e in open(f"cohort_{args.cohort}.txt").read().split(",") if e]
    rows5 = {r["commonName"]: r for r in csv.DictReader(open(args.f5, encoding="utf-8"))}
    checks = {(r["event"], r["label"]): r for r in csv.DictReader(open(args.checks, encoding="utf-8"))}

    print(f"\ncohort {args.cohort}: {len(events)} events, labels {args.labels}")
    out = {"cohort": args.cohort, "labels": args.labels, "code_sha256": chash, "pairs": {}, "fired": []}
    dropped = []
    for ev in events:
        r5 = rows5[ev]
        labels = []
        if args.labels in ("both", "xphm"):
            labels.append(r5["xphm_label"])
        if args.labels in ("both", "eob"):
            labels.append(r5["eob_label"])
        for label in labels:
            row = checks.get((ev, label))
            if row is None:
                print(f"  {ev} {label}: MISSING from {args.checks}")
                out["pairs"][f"{ev}|{label}"] = {"error": "missing"}
                continue
            e = mr.load(args, ev, label, r5, chash)
            leak, noise = leakage_samples(e, args.scan_trials)
            bad = verify(row, leak, noise)
            if bad:
                raise SystemExit(f"{ev} {label}: recomputed summaries differ from the frozen CSV:\n  "
                                 + "\n  ".join(bad))
            ks = ks_2samp(leak, noise)
            fails = []
            if abs(float(row["null_mean_z"])) * 1.0 > NULL_BIAS:
                fails.append("1")
            if float(row["ab_scan_median_diff"]) > AB_MEDIAN:
                fails.append("1a")
            if ks.pvalue < KS_ALPHA:
                fails.append("2a")
            out["pairs"][f"{ev}|{label}"] = {
                "ks_p": float(ks.pvalue), "ks_stat": float(ks.statistic),
                "leak_std": float(leak.std()), "noise_std": float(noise.std()),
                "null_mean": float(row["null_mean_z"]), "ab_median": float(row["ab_scan_median_diff"]),
                "fails": fails}
            if fails:
                dropped.append((ev, label, fails))
            print(f"  {ev:20s} {label:30s} KS p {ks.pvalue:7.4f}  leak/noise std "
                  f"{leak.std():.2f}/{noise.std():.2f}  null {float(row['null_mean_z']):+.2f}  "
                  f"A/B {float(row['ab_scan_median_diff']):.2f}  "
                  f"{'FAILS ' + ','.join(fails) if fails else 'ok'}")

    npairs = len([k for k, v in out["pairs"].items() if "error" not in v])
    print(f"\nleaving the primary test: {len(dropped)} of {npairs}")
    for ev, label, f in dropped:
        print(f"   {ev} {label}: rule(s) " + ", ".join(RULES[x] for x in f))
    if npairs and len([d for d in dropped if "2a" in d[2]]) > npairs / 3:
        out["fired"].append("2a: more than one third of the pairs fail the GR-leakage test "
                            "-> the method is not interpreted")
    if npairs and len([d for d in dropped if "1" in d[2]]) > npairs / 3:
        out["fired"].append("1: more than one third of the pairs show a null bias "
                            "-> the catalog result is not interpreted")

    for path in [p for p in args.catalog.split(",") if p]:
        c = json.load(open(path, encoding="utf-8"))
        r = c["all"]
        lo, hi = binomial(0.90, c["exp_trials"])
        low = [k for k, v in r["coverage"].items() if v["neyman"] < lo]
        print(f"\n{os.path.basename(path)}: N_eff {c['N_eff']:.2f}; expected 90% interval "
              f"[{r['exp_lower_median']:+.2e}, {r['exp_upper_median']:+.2e}] m^2; contiguous "
              f"{r['contiguous']:.2f}; edge {r['edge']:.3f}; coverage "
              f"{min(v['neyman'] for v in r['coverage'].values()):.3f}-"
              f"{max(v['neyman'] for v in r['coverage'].values()):.3f} (rule: not below {lo:.3f})")
        if low:
            out["fired"].append(f"0 ({os.path.basename(path)}): coverage below {lo:.3f} at {low}")
        if r["edge"] > 0.01:
            out["fired"].append(f"0 ({os.path.basename(path)}): {r['edge']:.3f} of the intervals at the grid edge")

    if args.offsource and os.path.exists(args.offsource):
        rows = [r for r in csv.DictReader(open(args.offsource, encoding="utf-8")) if r["inj"] == "C"]
        n = len(rows)
        cov = sum(r["covered"] == "True" for r in rows) / n if n else float("nan")
        lo, hi = binomial(0.90, n)
        print(f"\noff-source (set C, {n} segments): coverage {cov:.3f}; 99% interval {lo:.3f}-{hi:.3f}")
        if n and cov < lo:
            out["fired"].append(f"0a: off-source coverage {cov:.3f} below {lo:.3f}")

    print("\n" + ("STOPPING RULES FIRED:" if out["fired"] else "no stopping rule fired"))
    for f in out["fired"]:
        print("   " + f)
    print("\nReal-data estimates may be made only if no stopping rule fired (MODULE1_DESIGN §11).")
    json.dump(out, open(args.out, "w"), indent=1)
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
