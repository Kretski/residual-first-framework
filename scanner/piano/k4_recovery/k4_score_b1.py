#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
k4_score_b1.py — B1 scorer for the blind k⁴-equivalent recovery test.

Specification: RFF_K4_RECOVERY_TEST_SPEC_B0_v0.4_CORRECTED.md
               SHA-256 662b2303aff00bb1000d8a0ed93be64fc5f7fdc72aae811567fb2f5810a4d042

The only code that reads truth.jsonl. It
  1. verifies the hashes written by the run (hashes.json) before reading anything,
  2. joins truth and analysis by case id, asserts unique spawn indices and the
     per-cell counts of the manifest,
  3. --smoke: prints structural checks only (counts, schema, hashes, W3 cross-check
     failures) and NO rates (B0 §6.4),
  4. otherwise: evaluates C1, C2a (diagnostic), C2b, C3a, C3b, C3c, C4 / C4-PURE,
     C5 and the verdict ladder of B0 §7, and writes score.json.

Usage:  python k4_score_b1.py <run_dir> [--smoke]
"""

import argparse
import hashlib
import json
import os
import sys

import numpy as np
from scipy import stats

SPEC_SHA256 = "662b2303aff00bb1000d8a0ed93be64fc5f7fdc72aae811567fb2f5810a4d042"

# Thresholds of B0 v0.4 (§5); none may be changed after freezing
EVALUABLE_MIN = 0.80
SCORED_MIN_CENTS = 2.0
C1_UPPER = 0.02
C2_DET_LOWER = 0.95
C2_SIGN_LOWER = 0.98
C2_DIFF = 0.02
C3A_CORRECT_LOWER = 0.80
C3A_WRONG_UPPER = 0.02
C3B_WRONG_UPPER = 0.02
C3B_N3_LOWER = 0.50
C3C_UPPER = 0.05
C4_BIAS = 0.10
C4_CELL_COV = 0.85
C4_POOLED = (0.92, 0.99)
C5_LOWER = 0.90

ANALYSIS_KEYS = {"case_id", "excluded", "reason", "f_nom", "seconds"}
ANALYSED_KEYS = {"winner", "significant", "distinguishable", "mixture_like", "declared",
                 "a_winner", "sign_winner", "B_hat", "B_se", "B_lo", "B_hi", "sign_B",
                 "per", "discrimination", "n_partials", "n_fit"}
TRUTH_KEYS = {"case_id", "spawn_index", "cell", "block", "sign", "p_true", "B_true",
              "level", "f0", "noise_db", "jitter_cents"}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cp_lower(k, n):
    if n == 0 or k == 0:
        return 0.0
    return float(stats.beta.ppf(0.025, k, n - k + 1))


def cp_upper(k, n):
    if n == 0 or k == n:
        return 1.0
    return float(stats.beta.ppf(0.975, k + 1, n - k))


def rate(k, n):
    return {"k": int(k), "n": int(n), "rate": (k / n if n else None),
            "cp_lower": cp_lower(k, n), "cp_upper": cp_upper(k, n)}


def d10_cents(cell):
    """True deviation of partial 10 in cents (B0 §5 'scored' rule)."""
    g = cell["gen"]
    if g["model"] == "power":
        return abs(g["cents_at_10"])
    if g["model"] == "stiff":
        return abs(1200.0 * np.log2(np.sqrt(1.0 + g["B"] * 100.0)))
    if g["model"] == "mixture":
        return abs(g["cents3"])
    return 0.0


# ----------------------------------------------------------------------------
REQUIRED_HASHED = {"manifest.json", "truth.jsonl", "analysis.jsonl", "run.log"}


def load(run_dir):
    hashes = json.load(open(os.path.join(run_dir, "hashes.json"), encoding="utf-8"))
    # B1 review B: exactly these four files must be listed, nothing missing, nothing extra
    if set(hashes) != REQUIRED_HASHED:
        sys.exit(f"ABORT: hashes.json lists {sorted(hashes)}; required exactly {sorted(REQUIRED_HASHED)}")
    for name, h in hashes.items():
        got = sha256_file(os.path.join(run_dir, name))
        if got != h:
            sys.exit(f"ABORT: {name} sha256 {got} != recorded {h}")
    manifest = json.load(open(os.path.join(run_dir, "manifest.json"), encoding="utf-8"))
    if manifest["spec_sha256"] != SPEC_SHA256:
        sys.exit("ABORT: the run was made under a different B0")
    truth = [json.loads(l) for l in open(os.path.join(run_dir, "truth.jsonl"), encoding="utf-8")]
    ana = [json.loads(l) for l in open(os.path.join(run_dir, "analysis.jsonl"), encoding="utf-8")]
    return hashes, manifest, truth, ana


def join_and_check(manifest, truth, ana):
    problems = []
    tmap = {t["case_id"]: t for t in truth}
    amap = {a["case_id"]: a for a in ana}
    if len(tmap) != len(truth) or len(amap) != len(ana):
        problems.append("duplicate case ids")
    if set(tmap) != set(amap):
        problems.append("truth and analysis case ids differ")
    spawn = [t["spawn_index"] for t in truth]
    if len(set(spawn)) != len(spawn):
        problems.append("spawn indices not unique")
    # B1 review D: the indices must be exactly 0 … n_cases − 1 (complete and in range)
    if set(spawn) != set(range(manifest["n_cases"])):
        out_of_range = sorted(i for i in set(spawn) if not (0 <= i < manifest["n_cases"]))
        problems.append(f"spawn indices are not exactly 0..{manifest['n_cases'] - 1}"
                        + (f" (out of range: {out_of_range[:5]})" if out_of_range else ""))
    if manifest["n_cases"] != len(truth):
        problems.append(f"manifest n_cases {manifest['n_cases']} != {len(truth)} truth rows")
    if manifest["shuffle_spawn_index"] in spawn:
        problems.append("shuffle spawn index used by a case")
    counts = {}
    for t in truth:
        counts[t["cell"]] = counts.get(t["cell"], 0) + 1
    for c in manifest["cells"]:
        if counts.get(c["cell"], 0) != c["m"]:
            problems.append(f"cell {c['cell']}: {counts.get(c['cell'], 0)} cases != {c['m']}")
    for t in truth:
        if set(t) != TRUTH_KEYS:
            problems.append(f"truth schema at {t['case_id']}")
            break
    for a in ana:
        need = ANALYSIS_KEYS | (set() if a.get("excluded") else ANALYSED_KEYS)
        if not need <= set(a):
            problems.append(f"analysis schema at {a.get('case_id')}: missing {sorted(need - set(a))}")
            break
    joined = [(tmap[k], amap[k]) for k in sorted(tmap) if k in amap]
    return joined, counts, problems


# ----------------------------------------------------------------------------
def score(manifest, joined):
    cells = {c["cell"]: c for c in manifest["cells"]}
    by_cell = {c: [] for c in cells}
    for t, a in joined:
        by_cell[t["cell"]].append((t, a))

    def analysed(cell_id):
        return [(t, a) for t, a in by_cell[cell_id] if not a["excluded"]]

    cell_info = {}
    for cid, c in cells.items():
        n_an = len(analysed(cid))
        cell_info[cid] = {"m": c["m"], "n_analysed": n_an,
                          "evaluable": (n_an / c["m"] >= EVALUABLE_MIN) if c["m"] else False,
                          "d10_cents": d10_cents(c),
                          "scored": d10_cents(c) >= SCORED_MIN_CENTS,
                          "mixture_like_rate": (sum(a["mixture_like"] for _, a in analysed(cid)) / n_an
                                                if n_an else None)}

    def cells_where(pred):
        return [cid for cid, c in cells.items() if pred(c) and cell_info[cid]["evaluable"]]

    res = {"cells": cell_info}

    # ---- C1 null --------------------------------------------------------------
    rows = analysed("N")
    k = sum(a["declared"] != "none" for _, a in rows)
    r1 = rate(k, len(rows))
    c1 = cell_info["N"]["evaluable"] and r1["cp_upper"] <= C1_UPPER
    res["C1"] = {"fpr": r1, "evaluable": cell_info["N"]["evaluable"], "pass": bool(c1)}

    # ---- C2a (diagnostic) and C2b (primary) ----------------------------------
    def sign_gate(cell_ids, sign_key):
        out = {}
        for s, lab in ((+1, "+"), (-1, "-")):
            rows = [(t, a) for cid in cell_ids if cells[cid]["sign"] == s for t, a in analysed(cid)]
            det = [(t, a) for t, a in rows if a["declared"] != "none"]
            corr = sum(a[sign_key] == t["sign"] for t, a in det)
            rd, rs = rate(len(det), len(rows)), rate(corr, len(det))
            out[lab] = {"cells": sorted({t["cell"] for t, _ in rows}), "detection": rd, "sign": rs,
                        "pass": bool(rd["cp_lower"] >= C2_DET_LOWER and rs["cp_lower"] >= C2_SIGN_LOWER)}
        accs = [out[l]["sign"]["rate"] for l in ("+", "-")]
        diff = abs(accs[0] - accs[1]) if None not in accs else None
        out["acc_diff"] = diff
        out["pass"] = bool(out["+"]["pass"] and out["-"]["pass"])
        return out

    scored_PS = cells_where(lambda c: c["block"] in ("P", "S") and d10_cents(c) >= SCORED_MIN_CENTS)
    c2a = sign_gate(scored_PS, "sign_winner")
    c2a["pass"] = bool(c2a["pass"] and c2a["acc_diff"] is not None and c2a["acc_diff"] <= C2_DIFF)
    res["C2a_diagnostic"] = c2a

    target = cells_where(lambda c: ((c["block"] == "P" and c["p_true"] == 3) or c["block"] == "S")
                         and d10_cents(c) >= SCORED_MIN_CENTS)
    res["C2b"] = sign_gate(target, "sign_B")
    res["C2b"]["target_cells_evaluable"] = target

    # ---- C3a pure shapes -------------------------------------------------------
    c3a = {}
    for p in (2, 3, 4):
        ids = cells_where(lambda c, p=p: c["block"] == "P" and c["p_true"] == p
                          and d10_cents(c) >= SCORED_MIN_CENTS)
        rows = [(t, a) for cid in ids for t, a in analysed(cid)]
        corr = sum(a["significant"] and a["distinguishable"] and a["winner"] == p for _, a in rows)
        wrong = sum(a["significant"] and a["distinguishable"] and a["winner"] != p for _, a in rows)
        rc, rw = rate(corr, len(rows)), rate(wrong, len(rows))
        c3a[str(p)] = {"cells": ids, "correct": rc, "wrong": rw,
                       "pass": bool(rc["cp_lower"] >= C3A_CORRECT_LOWER and rw["cp_upper"] <= C3A_WRONG_UPPER)}
    c3a["pass"] = all(c3a[str(p)]["pass"] for p in (2, 3, 4))
    res["C3a"] = c3a

    # ---- C3b stiff law ---------------------------------------------------------
    ids = cells_where(lambda c: c["block"] == "S" and d10_cents(c) >= SCORED_MIN_CENTS)
    rows = [(t, a) for cid in ids for t, a in analysed(cid)]
    n3 = sum(a["significant"] and a["distinguishable"] and a["winner"] == 3 for _, a in rows)
    wrong = sum(a["significant"] and a["distinguishable"] and a["winner"] != 3 for _, a in rows)
    per_level = {}
    for cid in ids:
        rr = analysed(cid)
        per_level[cid] = rate(sum(a["significant"] and a["distinguishable"] and a["winner"] == 3
                                  for _, a in rr), len(rr))
    rn3, rw = rate(n3, len(rows)), rate(wrong, len(rows))
    res["C3b"] = {"cells": ids, "n3_distinguishable": rn3, "wrong": rw, "per_cell": per_level,
                  "pass": bool(rn3["cp_lower"] >= C3B_N3_LOWER and rw["cp_upper"] <= C3B_WRONG_UPPER)}

    # ---- C3c out of family -----------------------------------------------------
    ids = cells_where(lambda c: c["block"] == "O")
    rows = [(t, a) for cid in ids for t, a in analysed(cid)]
    k = sum(a["significant"] and a["distinguishable"] and a["winner"] == 3 for _, a in rows)
    decl = {}
    for _, a in rows:
        decl[a["declared"]] = decl.get(a["declared"], 0) + 1
    r3c = rate(k, len(rows))
    res["C3c"] = {"cells": ids, "n3_claimed": r3c, "declared": decl,
                  "pass": bool(r3c["cp_upper"] <= C3C_UPPER)}

    # ---- C4 amplitude ------------------------------------------------------------
    c4cells = {}
    for cid, c in cells.items():
        if c["B_true"] is None or c["block"] not in ("P", "S"):
            continue
        rr = analysed(cid)
        ratios = [a["B_hat"] / t["B_true"] for t, a in rr]
        cov = sum(a["B_lo"] <= t["B_true"] <= a["B_hi"] for t, a in rr)
        bias = float(np.median(ratios) - 1.0) if ratios else None
        ok = (cell_info[cid]["evaluable"] and bias is not None and abs(bias) <= C4_BIAS
              and len(rr) > 0 and cov / len(rr) >= C4_CELL_COV)
        c4cells[cid] = {"bias": bias, "coverage": rate(cov, len(rr)),
                        "evaluable": cell_info[cid]["evaluable"], "cell_conditions": bool(ok)}

    def c4_rule(core_ids):
        missing = [cid for cid in core_ids if cid not in c4cells or not c4cells[cid]["evaluable"]]
        kk = sum(c4cells[cid]["coverage"]["k"] for cid in core_ids if cid in c4cells)
        nn = sum(c4cells[cid]["coverage"]["n"] for cid in core_ids if cid in c4cells)
        pooled = kk / nn if nn else None
        ok = (not missing and all(c4cells[cid]["cell_conditions"] for cid in core_ids)
              and pooled is not None and C4_POOLED[0] <= pooled <= C4_POOLED[1])
        return {"cells": core_ids, "non_evaluable": missing, "pooled_coverage": {"k": kk, "n": nn, "rate": pooled},
                "pass": bool(ok)}

    pure_core = [f"P_p3_c{c:g}_{s}" for c in (2.0, 5.0) for s in ("+", "-")]
    std_core = pure_core + ["S_B0.0001_+", "S_B0.0001_-"]
    c4 = {"per_cell": c4cells, "standard": c4_rule(std_core), "pure": c4_rule(pure_core)}
    # C4-range: largest scored |B| such that the per-cell conditions hold for both
    # signs at that level and at every smaller scored level (contiguous from the smallest)
    levels = sorted({c["level"] for c in cells.values() if c["block"] == "S"
                     and d10_cents(c) >= SCORED_MIN_CENTS})
    rng_max, per_level_ok = None, {}
    for L in levels:
        ok = all(c4cells.get(f"S_B{L:g}_{s}", {}).get("cell_conditions", False) for s in ("+", "-"))
        per_level_ok[f"{L:g}"] = ok
    for L in levels:
        if per_level_ok[f"{L:g}"]:
            rng_max = L
        else:
            break
    c4["range"] = {"per_level_both_signs": per_level_ok, "contiguous_max_absB": rng_max}
    res["C4"] = c4

    # ---- C5 robustness -------------------------------------------------------------
    c5 = {}
    for s, lab in ((+1, "+"), (-1, "-")):
        cid = f"C_mix_{lab}"
        rr = analysed(cid) if cell_info[cid]["evaluable"] else []
        k = sum(a["sign_B"] == s for _, a in rr)
        r5 = rate(k, len(rr))
        bias = float(np.median([a["B_hat"] / t["B_true"] for t, a in rr]) - 1.0) if rr else None
        c5[lab] = {"evaluable": cell_info[cid]["evaluable"], "sign_B": r5, "bias_reported": bias,
                   "pass": bool(r5["cp_lower"] >= C5_LOWER)}
    c5["pass"] = bool(c5["+"]["pass"] and c5["-"]["pass"])
    rs = {}
    for cid in [c for c in cells if c.startswith("R_")]:
        rr = analysed(cid)
        n_an = len(rr)
        det = [(t, a) for t, a in rr if a["declared"] != "none"]
        d = {"n_analysed": n_an, "evaluable": cell_info[cid]["evaluable"],
             "declared_not_none": rate(len(det), n_an)}
        if cells[cid]["B_true"] is not None:
            d["sign_B_among_detected"] = rate(sum(a["sign_B"] == t["sign"] for t, a in det), len(det))
            d["bias"] = (float(np.median([a["B_hat"] / t["B_true"] for t, a in rr]) - 1.0) if rr else None)
            d["coverage"] = rate(sum(a["B_lo"] <= t["B_true"] <= a["B_hi"] for t, a in rr), n_an)
            d["n3_distinguishable"] = rate(sum(a["significant"] and a["distinguishable"] and a["winner"] == 3
                                               for _, a in rr), n_an)
        rs[cid] = d
    c5["R_stratum_reported"] = rs
    res["C5"] = c5

    # ---- Verdict ladder (B0 §7) --------------------------------------------------
    if not res["C1"]["pass"]:
        verdict = "V0 — NOT VALIDATED (null control)"
    elif not res["C2b"]["pass"]:
        verdict = "V1 — NOT VALIDATED (target coefficient sign)"
    elif not (res["C3a"]["pass"] and res["C3c"]["pass"]):
        verdict = "V2 — SIGN OF B ONLY"
    elif not res["C3b"]["pass"]:
        verdict = ("V3 — PURE-POWER SHAPE ONLY"
                   + (" + quantitative (pure n³)" if c4["pure"]["pass"] else ", not quantitative"))
    elif not c4["standard"]["pass"]:
        verdict = "V4 — DETECTION, SIGN OF B AND SHAPE; NOT QUANTITATIVE"
    else:
        verdict = "V5 — EXTRACTION VALIDATED (synthetic, within the family {0, 2, 3, 4})"
    qualifiers = []
    if not res["C5"]["pass"]:
        qualifiers.append("an opposite-sign second mechanism can flip the inferred sign")
    if res["C3b"]["pass"]:
        qualifiers.append("n³ recognised in ≥ 50 % of pooled stiff-law cases; per-|B| rates reported; "
                          "no per-value claim")
    non_eval = [cid for cid, ci in cell_info.items() if not ci["evaluable"]]
    c4_range = None if verdict.startswith("V3") else rng_max   # B0 §7: no stiff-law range in V3
    # B1 review C: every qualifier required by B0 §7, in one reporting block
    # (diagnostic; nothing here changes the verdict)
    required_qualifiers = {
        "C2a_winning_deviation_sign": {
            "+": {"detection": c2a["+"]["detection"], "sign": c2a["+"]["sign"]},
            "-": {"detection": c2a["-"]["detection"], "sign": c2a["-"]["sign"]},
            "acc_diff": c2a["acc_diff"], "would_pass_bars": c2a["pass"]},
        "C5_sign": {"pass": res["C5"]["pass"],
                    "+": res["C5"]["+"]["sign_B"], "-": res["C5"]["-"]["sign_B"]},
        "C3b_pooled_qualifier": res["C3b"]["pass"],
        "C3b_per_cell_n3_distinguishable": res["C3b"]["per_cell"],
        "C3c_declared_shapes_for_n5": res["C3c"]["declared"],
        "mixture_like_rate_per_cell": {cid: ci["mixture_like_rate"] for cid, ci in cell_info.items()},
        "C4_range": {"contiguous_max_absB": c4_range,
                     "per_level_both_signs": (None if verdict.startswith("V3")
                                              else per_level_ok)},
        "non_evaluable_cells": non_eval,
        "R_stratum": res["C5"]["R_stratum_reported"],
    }
    res["verdict"] = {"verdict": verdict, "qualifiers": qualifiers,
                      "required_qualifiers": required_qualifiers}
    return res


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="B1 scorer: blind k⁴-equivalent recovery test")
    ap.add_argument("run_dir")
    ap.add_argument("--smoke", action="store_true", help="structural checks only, no rates")
    args = ap.parse_args()

    hashes, manifest, truth, ana = load(args.run_dir)
    joined, counts, problems = join_and_check(manifest, truth, ana)
    if args.smoke != (manifest["mode"] == "smoke"):
        problems.append(f"--smoke flag does not match run mode '{manifest['mode']}'")

    print("hashes verified:")
    for k, v in hashes.items():
        print(f"  {v}  {k}")
    print(f"mode {manifest['mode']}, cases {len(truth)}, cells {len(counts)}")

    if args.smoke:
        n_ex = sum(a["excluded"] for a in ana)
        analysed = [a for a in ana if not a["excluded"]]
        # every analysed smoke case must carry the W3 cross-check (B1 review, smoke item 1)
        missing_xc = [a["case_id"] for a in analysed if "crosscheck" not in a]
        xc = [a["crosscheck"]["ok"] for a in analysed if "crosscheck" in a]
        if missing_xc:
            problems.append(f"W3 cross-check missing in {len(missing_xc)} analysed cases "
                            f"(first: {missing_xc[:3]})")
        if len(xc) != len(analysed):
            problems.append(f"W3 cross-check count {len(xc)} != analysed cases {len(analysed)}")
        if not analysed:
            problems.append("no analysed cases: the W3 cross-check was not exercised")
        print(f"structural: truth rows {len(truth)}, analysis rows {len(ana)}, excluded {n_ex}")
        print(f"structural: W3 cross-check present in {len(xc)}/{len(analysed)} analysed cases, "
              f"failures {xc.count(False)}")
        print("structural problems: " + ("none" if not problems else "; ".join(problems)))
        print("(smoke mode: no rates are computed or printed)")
        return 1 if problems or (False in xc) else 0

    if problems:
        sys.exit("ABORT: " + "; ".join(problems))
    res = score(manifest, joined)
    res["inputs"] = hashes
    out = os.path.join(args.run_dir, "score.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, default=float, ensure_ascii=False)
    v = res["verdict"]
    print(f"verdict: {v['verdict']}")
    for q in v["qualifiers"]:
        print(f"  qualifier: {q}")
    rq = v["required_qualifiers"]
    for lab in ("+", "-"):
        c = rq["C2a_winning_deviation_sign"][lab]
        print(f"  C2a ({lab}, diagnostic): detection {c['detection']['k']}/{c['detection']['n']}, "
              f"sign {c['sign']['k']}/{c['sign']['n']}")
    print(f"  C3c declared shapes for n^5: {rq['C3c_declared_shapes_for_n5']}")
    print(f"  C4-range (contiguous max |B|): {rq['C4_range']['contiguous_max_absB']}")
    print(f"  non-evaluable cells: {rq['non_evaluable_cells']}")
    print("  (mixture_like rates and the R stratum are in score.json)")
    print(f"score.json sha256 {sha256_file(out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
