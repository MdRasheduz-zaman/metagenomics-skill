"""Snakemake script: provision one module DB (db.provision) via its canonical downloader.

Idempotent — dbprovision.provision() skips the download when the DB is already present. No
`from __future__ import annotations` (Snakemake prepends a preamble; it must stay first).
"""
import os
import sys

from metagx import dbprovision

tool = snakemake.params.tool          # noqa: F821
db_dir = snakemake.params.db_dir      # noqa: F821

result = dbprovision.provision(tool, db_dir, run=True)

os.makedirs(os.path.dirname(snakemake.output.sentinel), exist_ok=True)  # noqa: F821
with open(snakemake.output.sentinel, "w") as fh:                        # noqa: F821
    fh.write((result.get("skipped") or ("ok" if result.get("ok") else "failed")) + "\n")

if not result.get("ok"):
    sys.exit(f"db.provision '{tool}' failed: {(result.get('tail') or result.get('note') or '')[-500:]}")
