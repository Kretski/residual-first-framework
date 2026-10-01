#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download_pe.py — download the LVK PE release files of the SELECTED events of the
frozen event list, slowly and politely (Zenodo rate limits), with resume,
retries on HTTP 429/5xx (honouring Retry-After) and MD5 verification against
the checksums published by Zenodo.

File resolution uses exactly the same rules as f5_batch.py (latest version of
each release record; full event names from GPS; cosmological-prior files).

  python download_pe.py --events event_selection/event_list_v1.csv --dest ~/gwdata/pe
"""
import argparse
import csv
import hashlib
import os
import sys
import time

import requests

import f5_batch as f5

PAUSE_BETWEEN_FILES = 10      # s
MAX_RETRIES = 8


def md5sum(path, chunk=2 ** 22):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def download(url, path, size=None):
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
                    print(f"    HTTP {r.status_code}; waiting {wait} s (attempt {attempt}/{MAX_RETRIES})", flush=True)
                    time.sleep(wait)
                    continue
                if r.status_code == 200 and have:
                    have = 0                                   # server ignored Range: restart
                r.raise_for_status()
                mode = "ab" if have else "wb"
                done = have
                t0 = time.time()
                with open(tmp, mode) as fh:
                    for chunk in r.iter_content(chunk_size=2 ** 20):
                        fh.write(chunk)
                        done += len(chunk)
                        if size and time.time() - t0 > 5:
                            print(f"    {done / 1e6:8.1f} / {size / 1e6:.1f} MB", end="\r", flush=True)
                            t0 = time.time()
            os.replace(tmp, path)
            return True
        except (requests.RequestException, OSError) as e:
            wait = min(600, 30 * 2 ** (attempt - 1))
            print(f"    {type(e).__name__}: {e}; waiting {wait} s (attempt {attempt}/{MAX_RETRIES})", flush=True)
            time.sleep(wait)
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="event_selection/event_list_v1.csv")
    ap.add_argument("--dest", default=os.path.expanduser("~/gwdata/pe"))
    args = ap.parse_args()
    os.makedirs(args.dest, exist_ok=True)
    events = [r for r in csv.DictReader(open(args.events, encoding="utf-8")) if r["status"] == "SELECTED"]

    files_by_cat = {}
    for cat, rec in f5.RECORDS.items():
        files_by_cat[cat] = f5.latest_record(rec).get("files", [])
        time.sleep(2)

    plan, total = [], 0
    for ev in events:
        try:
            fm, why = f5.pick_file(files_by_cat[ev["catalog"]], ev["commonName"], ev.get("GPS"))
        except ValueError as e:
            fm, why = None, str(e)
        if fm is None:
            print(f"  {ev['commonName']}: cannot resolve file ({why})")
            continue
        plan.append((ev["commonName"], fm))
        total += fm.get("size") or 0
    print(f"{len(plan)} files, total {total / 1e9:.1f} GB → {args.dest}")

    failed = []
    for i, (name, fm) in enumerate(plan, 1):
        path = os.path.join(args.dest, fm["key"])
        expected = (fm.get("checksum") or "").replace("md5:", "")
        if os.path.exists(path) and expected and md5sum(path) == expected:
            print(f"[{i:2d}/{len(plan)}] {name:18s} already present, MD5 ok")
            continue
        url = fm.get("links", {}).get("self") or fm.get("links", {}).get("content")
        print(f"[{i:2d}/{len(plan)}] {name:18s} {fm.get('size', 0) / 1e6:.0f} MB", flush=True)
        ok = download(url, path, fm.get("size"))
        if ok and expected:
            got = md5sum(path)
            if got != expected:
                print(f"    MD5 MISMATCH ({got} ≠ {expected}); file removed")
                os.remove(path)
                ok = False
            else:
                print("    MD5 ok" + " " * 30)
        if not ok:
            failed.append(name)
        time.sleep(PAUSE_BETWEEN_FILES)
    print(f"\ndone; failed: {failed if failed else 'none'}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
