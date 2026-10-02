#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download_strain.py — Module 1: download GWOSC strain around every event that
passed F5 (status BOTH_OK or SINGLE_MODEL; the separately analysed events are
included, the F5 FAIL event is not).

For each event and each detector used in its PE (detector list from the F5
output), the GWOSC event file of 4096 s duration at 4096 Hz is downloaded
(the event version listed in the frozen event list). One file serves both the
on-source analysis segment and the off-source segments of the null suite.

Polite downloading: one file at a time, pause between files, resume, retries on
HTTP 429/5xx honouring Retry-After. A manifest records URL, size and SHA-256 of
every file (GWOSC does not publish checksums through this API).

  python download_strain.py --events event_selection/event_list_v1.csv \
         --f5 f5/f5_results.csv --dest ~/gwdata/strain
"""
import argparse
import csv
import hashlib
import os
import sys
import time

import requests

DURATION = 4096
SAMPLE_RATE = 4096
PAUSE = 5
MAX_RETRIES = 8

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def sha256(path, chunk=2 ** 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def download(url, path):
    tmp = path + ".part"
    for attempt in range(1, MAX_RETRIES + 1):
        have = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        headers = {"User-Agent": "residual-first-framework/module1"}
        if have:
            headers["Range"] = f"bytes={have}-"
        try:
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                if r.status_code in (429, 500, 502, 503, 504):
                    wait = int(r.headers.get("Retry-After", 0)) or min(600, 30 * 2 ** (attempt - 1))
                    print(f"    HTTP {r.status_code}; waiting {wait} s ({attempt}/{MAX_RETRIES})", flush=True)
                    time.sleep(wait)
                    continue
                if r.status_code == 200 and have:
                    have = 0
                r.raise_for_status()
                with open(tmp, "ab" if have else "wb") as fh:
                    for chunk in r.iter_content(chunk_size=2 ** 20):
                        fh.write(chunk)
            os.replace(tmp, path)
            return True
        except (requests.RequestException, OSError) as e:
            wait = min(600, 30 * 2 ** (attempt - 1))
            print(f"    {type(e).__name__}: {e}; waiting {wait} s ({attempt}/{MAX_RETRIES})", flush=True)
            time.sleep(wait)
    return False


def event_urls(name, version, ifo):
    from gwosc.locate import get_event_urls
    kw = dict(duration=DURATION, sample_rate=SAMPLE_RATE, detector=ifo, format="hdf5")
    try:
        return get_event_urls(name, version=int(version) if version else None, **kw)
    except Exception as e:
        print(f"    GWOSC lookup with version {version} failed ({type(e).__name__}: {e}); trying latest")
        return get_event_urls(name, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="event_selection/event_list_v1.csv")
    ap.add_argument("--f5", default="f5/f5_results.csv")
    ap.add_argument("--dest", default=os.path.expanduser("~/gwdata/strain"))
    args = ap.parse_args()
    os.makedirs(args.dest, exist_ok=True)

    ev = {r["commonName"]: r for r in csv.DictReader(open(args.events, encoding="utf-8"))
          if r["status"] == "SELECTED"}
    f5 = [r for r in csv.DictReader(open(args.f5, encoding="utf-8"))
          if r["status"] in ("BOTH_OK", "SINGLE_MODEL")]
    plan = []
    for r in f5:
        name = r["commonName"]
        ifos = [i for i in r["ifos"].split(",") if i]
        for ifo in ifos:
            plan.append((name, ev[name]["version"], ifo))
    print(f"{len(f5)} events, {len(plan)} detector files (4096 s, 4096 Hz) → {args.dest}")

    man_path = os.path.join(args.dest, "strain_manifest.csv")
    done = {}
    if os.path.exists(man_path):
        for m in csv.DictReader(open(man_path, encoding="utf-8")):
            done[(m["event"], m["ifo"])] = m
    rows = list(done.values())
    failed = []
    for i, (name, version, ifo) in enumerate(plan, 1):
        if (name, ifo) in done and os.path.exists(os.path.join(args.dest, done[(name, ifo)]["file"])):
            print(f"[{i:3d}/{len(plan)}] {name:18s} {ifo}  already present")
            continue
        try:
            urls = event_urls(name, version, ifo)
        except Exception as e:
            print(f"[{i:3d}/{len(plan)}] {name:18s} {ifo}  GWOSC lookup FAILED: {type(e).__name__}: {e}")
            failed.append((name, ifo, "lookup"))
            continue
        if len(urls) != 1:
            print(f"[{i:3d}/{len(plan)}] {name:18s} {ifo}  expected 1 URL, got {len(urls)}: {urls[:3]}")
            failed.append((name, ifo, f"{len(urls)} urls"))
            continue
        url = urls[0]
        fname = url.rstrip("/").split("/")[-1]
        path = os.path.join(args.dest, fname)
        print(f"[{i:3d}/{len(plan)}] {name:18s} {ifo}  {fname}", flush=True)
        if not download(url, path):
            failed.append((name, ifo, "download"))
            continue
        rows.append({"event": name, "version": version, "ifo": ifo, "url": url, "file": fname,
                     "bytes": os.path.getsize(path), "sha256": sha256(path)})
        with open(man_path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["event", "version", "ifo", "url", "file", "bytes", "sha256"])
            w.writeheader()
            w.writerows(rows)
        time.sleep(PAUSE)
    total = sum(int(r["bytes"]) for r in rows)
    print(f"\nfiles in manifest: {len(rows)}, total {total / 1e9:.1f} GB; failed: {failed if failed else 'none'}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
