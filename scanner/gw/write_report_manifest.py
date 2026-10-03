#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
write_report_manifest.py — freeze manifest for the REPORTING layer (module1_report.py),
separate from FREEZE_v1.0.json, which covers the computational core.

  python write_report_manifest.py v1.1-report
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import module1_run as mr        # noqa: E402

OUT = os.path.join(HERE, "FREEZE_v1.1_report.json")
FILES = ["module1_report.py", "write_report_manifest.py"]

if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "v1.1-report"
    if os.path.exists(OUT):
        raise SystemExit(f"{OUT} exists; refusing to overwrite a freeze manifest")
    core, core_per = mr.code_hash()
    per = {f: mr.sha256_file(os.path.join(HERE, f), normalise_eol=True) for f in FILES}
    import hashlib
    h = hashlib.sha256()
    for f in FILES:
        h.update(f.encode() + b"\0" + per[f].encode() + b"\n")
    json.dump({"tag": tag, "report_sha256": h.hexdigest(), "files": per, "report_files": FILES,
               "core_sha256": core, "core_files": core_per,
               "rules": {"null_bias_sigma": 0.3, "ab_median_sigma": 0.5, "ks_alpha": 0.01,
                         "verify_tolerance": 1e-12},
               "versions": mr.versions(),
               "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
              open(OUT, "w", encoding="utf-8"), indent=1)
    print(f"report SHA-256 {h.hexdigest()}")
    for f in FILES:
        print(f"  {per[f]}  {f}")
    print(f"core SHA-256 {core}")
    print(f"written {OUT} for tag {tag}")
