"""
Download + extract the NANOGrav 15-yr data set (v2.1.0) from Zenodo.

  python fetch_nanograv15.py                 # download, verify MD5, extract, inspect
  python fetch_nanograv15.py --inspect-only  # only print what is on disk

* Resumes an interrupted download (HTTP Range).
* Verifies the MD5 published on Zenodo.
* Extracts READMEs, narrowband .par files and the narrowband residual
  tables (not the templates / MCMC chains — saves disk space).
  Use --all to extract everything.
* Prints the residual files it found and the first lines of one of them,
  so the column format can be checked.

Data: The NANOGrav Collaboration (2025), The NANOGrav 15-Year Data Set
(v2.1.0), Zenodo, doi:10.5281/zenodo.16051178, CC-BY 4.0.
Cite Agazie et al. 2023, ApJL 951, L9 when using it.

Author: Dimitar Kretski
License: MIT
"""

import hashlib
import os
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

URL = ("https://zenodo.org/records/16051178/files/"
       "NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz?download=1")
ARCHIVE = Path("NANOGrav15yr_PulsarTiming_v2.1.0.tar.gz")
MD5 = "557d42dd8486a5f8272d90dec9b228a8"
EXPECTED_SIZE = 638.7e6
DATA_DIR = Path("NANOGrav15yr")

KEEP_EXT = (".par", ".txt", ".dat", ".res", ".resid", ".csv", ".tsv", ".md")


def download(url=URL, dest=ARCHIVE, chunk=1 << 20):
    have = dest.stat().st_size if dest.exists() else 0
    req = urllib.request.Request(url, headers={"User-Agent": "residual-first-framework"})
    if have:
        req.add_header("Range", f"bytes={have}-")
        print(f"Resuming from {have/1e6:.1f} MB")
    try:
        resp = urllib.request.urlopen(req, timeout=60)
    except urllib.error.HTTPError as e:
        if e.code == 416:          # already complete
            print("Archive already fully downloaded.")
            return dest
        raise
    if have and resp.status != 206:   # server ignored Range → start over
        print("Server does not support resume — restarting download.")
        have = 0
    total = have + int(resp.headers.get("Content-Length", 0))
    mode = "ab" if have else "wb"
    t0, done = time.time(), have
    with open(dest, mode) as f:
        while True:
            block = resp.read(chunk)
            if not block:
                break
            f.write(block)
            done += len(block)
            rate = (done - have) / max(time.time() - t0, 1e-6) / 1e6
            pct = 100 * done / total if total else 0
            print(f"\r  {done/1e6:7.1f} / {total/1e6:.1f} MB  ({pct:5.1f}%)  {rate:5.1f} MB/s",
                  end="", flush=True)
    print()
    return dest


def md5sum(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _wanted(name: str, extract_all: bool) -> bool:
    if extract_all:
        return True
    low = name.lower()
    base = os.path.basename(low)
    if base.startswith("readme") or base == "description.txt":
        return True
    if "/residuals/" in low.replace("\\", "/") and base.endswith(".res"):
        return True                       # *_nb.avg.res / *_nb.full.res (v2.1.0)
    if "narrowband" not in low:
        return False
    if base.endswith(".par"):
        return True
    return "resid" in low and base.endswith(KEEP_EXT)


def extract(archive=ARCHIVE, out=DATA_DIR, extract_all=False):
    out.mkdir(exist_ok=True)
    n = 0
    with tarfile.open(archive, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.isfile() and _wanted(m.name, extract_all)]
        for m in members:
            target = (out / m.name).resolve()
            if not str(target).startswith(str(out.resolve())):   # path-traversal guard
                continue
            try:
                tar.extract(m, out, filter="data")      # Python ≥ 3.12 (safe extraction)
            except TypeError:
                tar.extract(m, out)
            n += 1
    print(f"Extracted {n} files to {out.resolve()}")


def inspect(out=DATA_DIR, n_show=8):
    if not out.exists():
        print(f"{out} does not exist yet.")
        return
    files = [p for p in out.rglob("*") if p.is_file()]
    res = sorted(p for p in files if p.suffix.lower() == ".res")
    pars = sorted(p for p in files if p.suffix.lower() == ".par")
    print(f"\n.par files: {len(pars)}   .res files: {len(res)} "
          f"(avg: {sum('.avg.' in p.name for p in res)}, full: {sum('.full.' in p.name for p in res)})")
    for p in res[:n_show]:
        print(f"   {p.relative_to(out)}   ({p.stat().st_size/1e3:.0f} kB)")

    readme = next((p for p in files if p.name == "README.residuals"), None)
    if readme:
        print(f"\n--- {readme.name} (full) ---")
        print(readme.read_text(errors="replace").rstrip())
        print("--- end ---")

    if res:
        sample = next((p for p in res if "J1909" in p.name and ".avg." in p.name), res[0])
        print(f"\n--- first 15 lines of {sample.name} ---")
        with open(sample, errors="replace") as f:
            for i, line in enumerate(f):
                if i >= 15:
                    break
                print(line.rstrip())
        print("--- end ---")
    else:
        print("\nNo .res files found. Run:  python fetch_nanograv15.py --all")


if __name__ == "__main__":
    if "--inspect-only" not in sys.argv:
        if not (ARCHIVE.exists() and ARCHIVE.stat().st_size >= EXPECTED_SIZE * 0.999):
            print(f"Downloading {ARCHIVE.name} (~{EXPECTED_SIZE/1e6:.0f} MB) from Zenodo…")
            download()
        print("Verifying MD5…", end=" ", flush=True)
        got = md5sum(ARCHIVE)
        if got != MD5:
            print(f"MISMATCH ({got}). Delete {ARCHIVE} and run again.")
            sys.exit(1)
        print("OK")
        extract(extract_all="--all" in sys.argv)
    inspect()
