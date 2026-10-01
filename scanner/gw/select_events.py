#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
select_events.py — Module 1 event selection (frozen BEFORE any residual is computed).

Selection rule (MODULE1_DESIGN §3), fixed here:
  - Catalogs (GWOSC event API, /eventapi/json/<catalog>/):
        discovery    : GWTC-2.1-confident (O1, O2, O3a), GWTC-3-confident (O3b)
        confirmation : GWTC-4.0 (O4a)
  - BBH criterion  : mass_2_source (catalog median, source frame) > 3 M_sun
  - SNR criterion  : network_matched_filter_snr of the catalog listing
                     (the parameter-estimation median, not a search-pipeline
                     value) >= 12.0
  - Missing values : the event is kept in the list as EXCLUDED with the reason.

Provenance: the raw JSON response of every catalog is saved, with its SHA-256,
the URL and the UTC access time, so that later catalog updates cannot change
the frozen list.

Outputs (in --out):
  catalog_snapshots/<catalog>.json     raw responses
  provenance.json                      URLs, access times, SHA-256, rule
  event_list_v1.csv                    every event with status and reason

  python select_events.py --out event_selection
  python select_events.py --out test_sel --offline-dir <folder with <catalog>.json>
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import sys
import urllib.request

API = "https://gwosc.org/eventapi/json/{catalog}/"
COHORTS = {
    "discovery": ["GWTC-2.1-confident", "GWTC-3-confident"],
    "confirmation": ["GWTC-4.0"],
}
RULE = {
    "bbh_field": "mass_2_source",
    "bbh_min_msun": 3.0,
    "snr_field": "network_matched_filter_snr",
    "snr_min": 12.0,
    "snr_inclusive": True,
}
FIELDS = ["cohort", "catalog", "commonName", "version", "GPS", "mass_1_source",
          "mass_2_source", "network_matched_filter_snr", "luminosity_distance",
          "redshift", "jsonurl", "status", "reason"]

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def fetch(catalog, out_dir, offline_dir=None):
    snap_dir = os.path.join(out_dir, "catalog_snapshots")
    os.makedirs(snap_dir, exist_ok=True)
    if offline_dir:
        url = os.path.join(offline_dir, f"{catalog}.json")
        with open(url, "rb") as fh:
            raw = fh.read()
    else:
        url = API.format(catalog=catalog)
        req = urllib.request.Request(url, headers={"User-Agent": "residual-first-framework/module1"})
        with urllib.request.urlopen(req, timeout=120) as r:
            raw = r.read()
    path = os.path.join(snap_dir, f"{catalog}.json")
    with open(path, "wb") as fh:
        fh.write(raw)
    return {"catalog": catalog, "url": url, "accessed_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}, json.loads(raw)


def classify(ev):
    m2 = ev.get(RULE["bbh_field"])
    snr = ev.get(RULE["snr_field"])
    if m2 is None:
        return "EXCLUDED", f"missing {RULE['bbh_field']}"
    if not m2 > RULE["bbh_min_msun"]:
        return "EXCLUDED", f"{RULE['bbh_field']} = {m2} <= {RULE['bbh_min_msun']} (not BBH)"
    if snr is None:
        return "EXCLUDED", f"missing {RULE['snr_field']}"
    if not snr >= RULE["snr_min"]:
        return "EXCLUDED", f"{RULE['snr_field']} = {snr} < {RULE['snr_min']}"
    return "SELECTED", ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="event_selection")
    ap.add_argument("--offline-dir", default=None, help="read <catalog>.json from here (testing)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    prov = {"rule": RULE, "cohorts": COHORTS, "catalogs": []}
    rows = []
    for cohort, catalogs in COHORTS.items():
        for cat in catalogs:
            try:
                meta, data = fetch(cat, args.out, args.offline_dir)
            except Exception as e:
                print(f"ERROR: catalog {cat} could not be fetched: {type(e).__name__}: {e}")
                print("Stop: the list must not be frozen with a missing catalog.")
                sys.exit(2)
            events = data.get("events", {})
            meta["n_events"] = len(events)
            prov["catalogs"].append(meta)
            for key, ev in sorted(events.items(), key=lambda kv: kv[1].get("GPS") or 0):
                status, reason = classify(ev)
                rows.append({
                    "cohort": cohort, "catalog": cat,
                    "commonName": ev.get("commonName", key), "version": ev.get("version"),
                    "GPS": ev.get("GPS"), "mass_1_source": ev.get("mass_1_source"),
                    "mass_2_source": ev.get("mass_2_source"),
                    "network_matched_filter_snr": ev.get("network_matched_filter_snr"),
                    "luminosity_distance": ev.get("luminosity_distance"),
                    "redshift": ev.get("redshift"), "jsonurl": ev.get("jsonurl"),
                    "status": status, "reason": reason})

    # duplicate names across catalogs within the run are reported, not resolved silently
    seen = {}
    for r in rows:
        seen.setdefault(r["commonName"], []).append(r["catalog"])
    dups = {k: v for k, v in seen.items() if len(v) > 1}
    prov["duplicates_across_catalogs"] = dups

    with open(os.path.join(args.out, "event_list_v1.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(args.out, "provenance.json"), "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=2)

    print("catalog snapshots:")
    for m in prov["catalogs"]:
        print(f"  {m['catalog']:22s} events {m['n_events']:4d}  sha256 {m['sha256'][:16]}…  {m['accessed_utc']}")
    for cohort in COHORTS:
        sel = [r for r in rows if r["cohort"] == cohort and r["status"] == "SELECTED"]
        exc = [r for r in rows if r["cohort"] == cohort and r["status"] == "EXCLUDED"]
        print(f"\n{cohort}: {len(sel)} selected, {len(exc)} excluded")
        for r in sel:
            print(f"   {r['commonName']:18s} {r['catalog']:20s} m2 = {r['mass_2_source']:>6}  "
                  f"SNR = {r['network_matched_filter_snr']}")
        reasons = {}
        for r in exc:
            reasons[r["reason"].split(" =")[0]] = reasons.get(r["reason"].split(" =")[0], 0) + 1
        print("   exclusion reasons:", reasons)
    if dups:
        print("\nWARNING: names appearing in more than one catalog:", dups)


if __name__ == "__main__":
    main()
