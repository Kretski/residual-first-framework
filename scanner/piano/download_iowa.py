#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
download_iowa.py — downloads the Module 0 (piano) data from University of Iowa MIS.

  data/Piano   ← Piano.mf.<note>.aiff  (Steinway B, 2001, 16 bit / 44.1 kHz, stereo)
  data/Violin  ← Violin.arco.ff.sul*.stereo.zip (2012, anechoic chamber) — development null test
  data/Cello   ← Cello.arco.ff.sul*.stereo.zip  (2012, anechoic chamber) — final null test

Run from the folder scanner\\piano:
  python download_iowa.py
Files already downloaded are skipped, so it can be run repeatedly.
"""
import os
import sys
import urllib.parse
import urllib.request
import zipfile

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

BASE = "https://theremin.music.uiowa.edu/"
PIANO_DIR = "sound files/MIS/Piano_Other/piano/"
STRINGS_DIR = "sound files/MIS Pitches - 2014/Strings/"
STRING_SETS = {"Violin": ["sulG", "sulD", "sulA", "sulE"],
               "Cello": ["sulC", "sulG", "sulD", "sulA"]}
NAMES = ["C", "Db", "D", "Eb", "E", "F", "Gb", "G", "Ab", "A", "Bb", "B"]
DYNAMIC = "mf"


def url(path):
    return BASE + urllib.parse.quote(path)


def fetch(u, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "already present"
    try:
        with urllib.request.urlopen(u, timeout=60) as r:
            data = r.read()
        if data[:4] not in (b"FORM", b"PK\x03\x04"):      # AIFF or ZIP
            return "not audio/zip (skipped)"
        with open(dest, "wb") as fh:
            fh.write(data)
        return f"ok ({len(data) / 1e6:.1f} MB)"
    except Exception as e:
        return f"error: {e}"


def main():
    os.makedirs("data/Piano", exist_ok=True)
    
    print("=== Piano, dynamic", DYNAMIC, "===")
    ok = 0
    for midi in range(21, 109):                       # A0 … C8
        name = f"{NAMES[midi % 12]}{midi // 12 - 1}"
        fn = f"Piano.{DYNAMIC}.{name}.aiff"
        status = fetch(url(PIANO_DIR + fn), os.path.join("data/Piano", fn))
        ok += status.startswith(("ok", "already"))
        print(f"  {fn:22s} {status}")
    print(f"piano: {ok} files")

    for inst, strings in STRING_SETS.items():
        folder = os.path.join("data", inst)
        os.makedirs(folder, exist_ok=True)
        print(f"\n=== {inst} arco (null test) ===")
        for st in strings:
            fn = f"{inst}.arco.ff.{st}.stereo.zip"
            dest = os.path.join(folder, fn)
            print(f"  {fn:34s} {fetch(url(STRINGS_DIR + inst + '/' + fn), dest)}")
            if os.path.exists(dest):
                with zipfile.ZipFile(dest) as z:
                    for m in z.namelist():
                        base = os.path.basename(m)
                        # skip macOS service files (__MACOSX, ._*)
                        if "__MACOSX" in m or base.startswith("._") or not base:
                            continue
                        if base.lower().endswith((".aif", ".aiff")):
                            out = os.path.join(folder, base)
                            if not os.path.exists(out):
                                with open(out, "wb") as fh:
                                    fh.write(z.read(m))
                os.remove(dest)
        n_f = len([f for f in os.listdir(folder) if f.lower().endswith((".aif", ".aiff"))])
        print(f"{inst}: {n_f} files")


if __name__ == "__main__":
    main()
