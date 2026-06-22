"""Snakemake script: build the kraken2 + Bracken DB per config['db']['build'].

Invoked by rules/dbbuild.smk. No `from __future__ import annotations` — Snakemake prepends a
preamble, so a __future__ import would not be the first statement (SyntaxError at run time).
"""
import sys

from metagx import dbbuild

build = dict(snakemake.config["db"]["build"])  # noqa: F821 (snakemake injects this)
db_dir = snakemake.params.db_dir               # noqa: F821

result = dbbuild.build_database(
    db_dir=db_dir,
    strategy=build.get("strategy", "standard"),
    taxonomy=build.get("taxonomy", "real"),
    libraries=build.get("libraries"),
    source=build.get("source"),
    read_lengths=build.get("read_lengths", [150]),
    threads=snakemake.threads,                 # noqa: F821
    kmer_len=build.get("kmer_len", 35),
    minimizer_len=build.get("minimizer_len", 31),
    minimizer_spaces=build.get("minimizer_spaces"),
    max_db_size=build.get("max_db_size"),
    no_masking=build.get("no_masking"),
    use_ftp=build.get("use_ftp", True),
    build_blast=bool(build.get("blast", False)),  # aligned BLAST validation DB, built together
    run=True,
)

# Provenance manifest doubles as the rule's output sentinel (so the build is idempotent).
dbbuild.write_manifest(db_dir, result)

if not result.get("ok"):
    sys.exit(f"db.build failed at step '{result.get('failed_step')}' — see logs under {db_dir}")
# The aligned BLAST DB is opt-in; if it was requested but failed, fail loudly (don't silently
# leave the user with an out-of-scope validation later).
_b = result.get("blast")
if _b and not _b.get("ok"):
    sys.exit(f"db.build: kraken2 DB built but the aligned BLAST DB failed: {_b.get('error') or _b}")
