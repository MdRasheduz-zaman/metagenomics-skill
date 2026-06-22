# Build the validation BLAST DB from the SAME genomes as the classifier (in-scope benchmark).
# No `from __future__` import (Snakemake prepends a preamble — it must be the first statement).
import os
import sys

from metagx import validation as v

snk = snakemake  # noqa: F821  (injected by Snakemake)
source = snk.params.source
prefix = snk.params.prefix
os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)

res = v.build_blast_db(source, prefix, run=True)
if not res.get("ok"):
    sys.exit(f"makeblastdb failed for {source!r}: {res.get('note') or res.get('tail') or res}")
# sentinel so the DAG can depend on a ready DB
open(snk.output.sentinel, "w").write(res.get("skipped", "built"))
