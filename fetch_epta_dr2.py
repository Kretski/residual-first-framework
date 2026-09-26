"""
Download + extract EPTA DR2 (Zenodo 10.5281/zenodo.8300645) and INSPECT it.

  python fetch_epta_dr2.py              # download (39 MB), verify MD5, extract, inspect
  python fetch_epta_dr2.py --inspect-only

This step only shows what is in the release (folders, .par/.tim files, the
first lines of one of each) and whether PINT is installed. No residuals are
computed and no statistics are run — the EPTA analysis must be frozen
(PREREGISTRATION_EPTA.md committed + tagged) before any residual is looked at.

Data: EPTA Collaboration (2023), A&A 678, A48; doi:10.5281/zenodo.8300645, CC-BY 4.0.

Author: Dimitar Kretski
License: MIT
"""
import sys
import zipfile
from pathlib import Path

from fetch_nanograv15 import download, md5sum

URL = "https://zenodo.org/records/8300645/files/EPTA-DR2.zip?download=1"
ARCHIVE = Path("EPTA-DR2.zip")
MD5 = "307b2b99cd352876b07409899448af8b"
EXPECTED_SIZE = 39.3e6
OUT = Path("EPTA_DR2")


def extract(archive=ARCHIVE, out=OUT):
    out.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        for m in z.infolist():
            target = (out / m.filename).resolve()
            if str(target).startswith(str(out.resolve())):      # path-traversal guard
                z.extract(m, out)
    print(f"Extracted to {out.resolve()}")


def head(path, n):
    print(f"\n--- first {n} lines of {path.name} ---")
    with open(path, errors="replace") as f:
        for i, line in enumerate(f):
            if i >= n:
                break
            print(line.rstrip())
    print("--- end ---")


def inspect(out=OUT):
    files = [p for p in out.rglob("*") if p.is_file()]
    dirs = sorted({str(p.parent.relative_to(out)) for p in files})
    pars = sorted(p for p in files if p.suffix == ".par")
    tims = sorted(p for p in files if p.suffix == ".tim")
    print(f"\nFiles: {len(files)}   .par: {len(pars)}   .tim: {len(tims)}")
    print(f"Folders ({len(dirs)}), first 25:")
    for d in dirs[:25]:
        print("  ", d)
    exts = {}
    for p in files:
        exts[p.suffix or "(none)"] = exts.get(p.suffix or "(none)", 0) + 1
    print("By extension:", dict(sorted(exts.items(), key=lambda x: -x[1])))
    for p in sorted(p for p in files if p.name.lower().startswith("readme"))[:3]:
        head(p, 40)
    par = next((p for p in pars if "J1909" in p.name), pars[0] if pars else None)
    tim = next((p for p in tims if "J1909" in p.name), tims[0] if tims else None)
    if par:
        head(par, 40)
    if tim:
        head(tim, 15)
    try:
        import pint
        print(f"\nPINT installed: version {pint.__version__}")
    except ImportError:
        print("\nPINT NOT installed. Install with:   pip install pint-pulsar")


if __name__ == "__main__":
    if "--inspect-only" not in sys.argv:
        if not (ARCHIVE.exists() and ARCHIVE.stat().st_size >= EXPECTED_SIZE * 0.99):
            print(f"Downloading {ARCHIVE.name} (~39 MB) from Zenodo…")
            download(URL, ARCHIVE)
        print("Verifying MD5…", end=" ", flush=True)
        got = md5sum(ARCHIVE)
        if got != MD5:
            print(f"MISMATCH ({got}). Delete {ARCHIVE} and run again.")
            sys.exit(1)
        print("OK")
        extract()
    inspect()
