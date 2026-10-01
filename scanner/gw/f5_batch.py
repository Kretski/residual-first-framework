#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
f5_batch.py — Module 1 feasibility check F5 for every SELECTED event of the
frozen event list (structure only; no residual, no waveform, no analysis).

For each event the LVK parameter-estimation (PE) release file is located through
the Zenodo API (latest version of the release record, resolved and recorded) and
read REMOTELY (HTTP range requests; only metadata and small datasets are fetched).

Pre-registered label rules (MODULE1_DESIGN §5):
  discovery    : XPHM label   "C01:IMRPhenomXPHM"
                 EOB  label   "C01:SEOBNRv4PHM"
  confirmation : XPHM label   first of ["C00:IMRPhenomXPHM-SpinTaylor", "C00:IMRPhenomXPHM"]
                 EOB  label   "C00:SEOBNRv5PHM"
  "Mixed" labels are never used. A label with a different name is NOT picked ad hoc:
  the event is marked FAIL with the list of available labels, for review.

Per label, F5 requires:
  posterior samples readable, with log_likelihood and the parameters needed for
  waveform generation; PSDs for every analysed detector; a calibration envelope
  (from the label itself or, if absent, from the XPHM label of the same event —
  the source is recorded); waveform settings readable from the configuration
  (bilby or LALInference keys).

Status per event:
  BOTH_OK         both models pass   → primary analysis
  SINGLE_MODEL    exactly one passes → analysed separately, not in the primary test
  FAIL            neither passes, or the file/labels could not be resolved

Outputs: f5_results.csv, f5_provenance.json (records, files, checksums).

  pip install fsspec aiohttp requests      (once, in the gw2 environment)
  python download_pe.py --events event_selection/event_list_v1.csv --dest ~/gwdata/pe
  python f5_batch.py --events event_selection/event_list_v1.csv --out f5 --local-dir ~/gwdata/pe
"""
import argparse
import csv
import json
import os
import re
import sys
import time

import h5py
import numpy as np
import requests

RECORDS = {                          # Zenodo release records (resolved to latest version)
    "GWTC-2.1-confident": "6513631",
    "GWTC-3-confident": "8177023",
    "GWTC-4.0": "17602505",
}
LABELS = {
    "discovery": {"XPHM": ["C01:IMRPhenomXPHM"], "EOB": ["C01:SEOBNRv4PHM"]},
    "confirmation": {"XPHM": ["C00:IMRPhenomXPHM-SpinTaylor", "C00:IMRPhenomXPHM"],
                     "EOB": ["C00:SEOBNRv5PHM"]},
}
NEEDED = ("mass_1", "mass_2", "luminosity_distance", "geocent_time", "a_1", "a_2",
          "tilt_1", "tilt_2", "phi_12", "phi_jl", "theta_jn", "psi", "ra", "dec",
          "phase", "log_likelihood")
CONFIG_MAP = {   # canonical name → (bilby key, LALInference key); keys are normalised
                 # so that '_' and '-' are equivalent (both spellings occur in releases)
    "f_low": ("minimum-frequency", "flow"),
    "f_ref": ("reference-frequency", "fref"),
    "f_high": ("maximum-frequency", "fhigh"),
    "srate": ("sampling-frequency", "srate"),
    "duration": ("duration", "seglen"),
    "fmin_template": (None, "fmin-template"),
    "waveform_arguments": ("waveform-arguments-dict", None),
    "mode_array": ("mode-array", None),
}
FIELDS = ["cohort", "commonName", "record_id", "file", "status", "xphm_label", "xphm_ok",
          "xphm_reason", "eob_label", "eob_ok", "eob_reason", "calib_source", "ifos",
          "available_labels", "settings_xphm", "settings_eob"]

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass


def text(x):
    if isinstance(x, (bytes, np.bytes_)):
        return x.decode(errors="replace")
    if isinstance(x, np.ndarray):
        return ", ".join(text(v) for v in x.ravel()[:20])
    return str(x)


def latest_record(rec_id):
    r = requests.get(f"https://zenodo.org/api/records/{rec_id}/versions/latest", timeout=60)
    if r.status_code != 200:
        r = requests.get(f"https://zenodo.org/api/records/{rec_id}", timeout=60)
    r.raise_for_status()
    return r.json()


def full_name(name, gps):
    """Short catalog names (e.g. GW150914, GW190521) are completed to the
    GWYYMMDD_HHMMSS form from the GPS time (UTC, leap seconds 17 before 2017,
    18 from 2017). Several events can share a date, e.g. GW190521_030229 and
    GW190521_074359, so the short name alone is ambiguous."""
    if "_" in name or gps in (None, ""):
        return name
    import datetime as dt
    g = float(gps)
    leap = 17 if g < 1167264018 else 18
    t = dt.datetime(1980, 1, 6) + dt.timedelta(seconds=g - leap)
    full = "GW" + t.strftime("%y%m%d_%H%M%S")
    if full[:8] != name[:8]:
        raise ValueError(f"GPS {gps} gives {full}, inconsistent with {name}")
    return full


def pick_file(files, name, gps=None):
    """Select the PE file of an event: exact full-name match, .h5/.hdf5, prefer the
    cosmological distance prior ('cosmo' but not 'nocosmo') when several variants exist."""
    name = full_name(name, gps)
    cands = [f for f in files
             if f["key"].lower().endswith((".h5", ".hdf5"))
             and ("PEDataRelease" in f["key"] or "PE" in f["key"])
             and (f"{name}_" in f["key"] or f"{name}." in f["key"] or f"{name}-" in f["key"])]
    if not cands:
        return None, "no PE file found"
    cosmo = [f for f in cands if "cosmo" in f["key"] and "nocosmo" not in f["key"]]
    if len(cosmo) == 1:
        return cosmo[0], ""
    if len(cands) == 1:
        return cands[0], ""
    return None, "ambiguous PE files: " + "; ".join(f["key"] for f in cands[:6])


def open_remote(url):
    import fsspec
    fobj = fsspec.open(url, mode="rb", block_size=2 ** 20, cache_type="blockcache").open()
    return h5py.File(fobj, "r")


def posterior_names(g):
    ps = g["posterior_samples"]
    if isinstance(ps, h5py.Dataset) and ps.dtype.names:
        return set(ps.dtype.names)
    if isinstance(ps, h5py.Group):
        keys = list(ps.keys())
        if "parameter_names" in keys:
            return {text(v) for v in ps["parameter_names"][()]}
        return set(keys)
    return set()


def config(g):
    flat = {}
    if "config_file" not in g:
        return flat

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            flat[name.split("/")[-1].replace("_", "-")] = text(obj[()])
    g["config_file"].visititems(visit)
    out = {}
    for canon, (kb, kl) in CONFIG_MAP.items():
        for k in (kb, kl):
            if k and k in flat:
                out[canon] = flat[k]
                break
    return out


def check_label(f, label):
    """Return (ok, reason, ifos, has_calib, settings)."""
    if label not in f:
        return False, "label missing", [], False, {}
    g = f[label]
    try:
        names = posterior_names(g)
    except Exception as e:
        return False, f"posterior unreadable ({type(e).__name__})", [], False, {}
    missing = [p for p in NEEDED if p not in names]
    if missing:
        return False, "missing parameters: " + ",".join(missing), [], False, {}
    ifos = sorted(g["psds"].keys()) if "psds" in g else []
    if not ifos:
        return False, "no PSDs", [], False, {}
    calib = "calibration_envelope" in g and len(g["calibration_envelope"].keys()) > 0
    if calib:
        missing_cal = [i for i in ifos if i not in g["calibration_envelope"]]
        if missing_cal:
            calib = False
    cfg = config(g)
    for need in ("f_low", "f_ref", "srate", "duration"):
        if need not in cfg:
            return False, f"config lacks {need}", ifos, calib, cfg
    return True, "", ifos, calib, cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", default="event_selection/event_list_v1.csv")
    ap.add_argument("--out", default="f5")
    ap.add_argument("--local-dir", default=os.path.expanduser("~/gwdata/pe"),
                    help="read PE files from here when present (download_pe.py); "
                         "otherwise read remotely")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    events = [r for r in csv.DictReader(open(args.events, encoding="utf-8")) if r["status"] == "SELECTED"]
    print(f"selected events: {len(events)}")

    prov = {"records": {}, "label_rules": LABELS, "needed_parameters": NEEDED, "files": {}}
    files_by_cat = {}
    for cat, rec in RECORDS.items():
        meta = latest_record(rec)
        files_by_cat[cat] = meta.get("files", [])
        prov["records"][cat] = {"requested": rec, "resolved_id": meta.get("id"),
                                "version": meta.get("metadata", {}).get("version"),
                                "n_files": len(files_by_cat[cat])}
        print(f"  {cat}: record {rec} → resolved {meta.get('id')} "
              f"(version {prov['records'][cat]['version']}), {len(files_by_cat[cat])} files")

    rows = []
    for i, ev in enumerate(events, 1):
        name, cat, cohort = ev["commonName"], ev["catalog"], ev["cohort"]
        row = {k: "" for k in FIELDS}
        row.update(cohort=cohort, commonName=name,
                   record_id=prov["records"][cat]["resolved_id"])
        t0 = time.time()
        try:
            fmeta, why = pick_file(files_by_cat[cat], name, ev.get("GPS"))
        except ValueError as e:
            fmeta, why = None, str(e)
        if fmeta is None:
            row.update(status="FAIL", xphm_reason=why, eob_reason=why)
            rows.append(row)
            print(f"[{i:2d}/{len(events)}] {name:18s} FAIL: {why}")
            continue
        url = fmeta.get("links", {}).get("self") or fmeta.get("links", {}).get("content")
        row["file"] = fmeta["key"]
        prov["files"][name] = {"key": fmeta["key"], "size": fmeta.get("size"),
                               "checksum": fmeta.get("checksum"), "url": url}
        local = os.path.join(args.local_dir, fmeta["key"])
        row["file"] = fmeta["key"] + (" (local)" if os.path.exists(local) else " (remote)")
        try:
            with (h5py.File(local, "r") if os.path.exists(local) else open_remote(url)) as f:
                labels = [k for k in f.keys() if k not in ("version", "history")]
                row["available_labels"] = " | ".join(labels)
                res = {}
                for fam in ("XPHM", "EOB"):
                    chosen = next((l for l in LABELS[cohort][fam] if l in labels), None)
                    if chosen is None:
                        res[fam] = (None, False, "preferred label not present", [], False, {})
                    else:
                        res[fam] = (chosen,) + check_label(f, chosen)
        except Exception as e:
            row.update(status="FAIL", xphm_reason=f"file unreadable ({type(e).__name__}: {e})")
            rows.append(row)
            print(f"[{i:2d}/{len(events)}] {name:18s} FAIL: unreadable")
            continue
        (lx, okx, rx, ifx, cx, sx), (le, oke, re_, ife, ce, se) = res["XPHM"], res["EOB"]
        calib_ok_x = okx and cx
        calib_src = "own" if cx else ""
        if oke and not ce and calib_ok_x:
            calib_src = "own (XPHM); EOB borrows XPHM envelope"
        elif oke and ce:
            calib_src = "own (both)"
        if okx and not cx:
            okx, rx = False, "no calibration envelope"
        if oke and not (ce or cx):
            oke, re_ = False, "no calibration envelope (none to borrow)"
        status = "BOTH_OK" if (okx and oke) else ("SINGLE_MODEL" if (okx or oke) else "FAIL")
        row.update(status=status, xphm_label=lx or "", xphm_ok=okx, xphm_reason=rx,
                   eob_label=le or "", eob_ok=oke, eob_reason=re_, calib_source=calib_src,
                   ifos=",".join(ifx or ife), settings_xphm=json.dumps(sx),
                   settings_eob=json.dumps(se))
        rows.append(row)
        print(f"[{i:2d}/{len(events)}] {name:18s} {status:12s} XPHM {'ok' if okx else rx[:40]:40s} "
              f"EOB {'ok' if oke else re_[:40]:40s} ({time.time() - t0:.0f} s)", flush=True)
        with open(os.path.join(args.out, "f5_results.csv"), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            w.writerows(rows)

    with open(os.path.join(args.out, "f5_results.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(args.out, "f5_provenance.json"), "w", encoding="utf-8") as fh:
        json.dump(prov, fh, indent=2, default=str)
    for cohort in LABELS:
        sub = [r for r in rows if r["cohort"] == cohort]
        counts = {s: sum(r["status"] == s for r in sub) for s in ("BOTH_OK", "SINGLE_MODEL", "FAIL")}
        print(f"\n{cohort}: {len(sub)} events → {counts}")


if __name__ == "__main__":
    main()
