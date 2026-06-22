# Build the validation BLAST DB IN SCOPE with the classifier.
#
#   build_from: <FASTA|folder>  -> makeblastdb on those genomes (organism via subject title).
#   build_from: "classifier"    -> resolve the kraken2 DB dir's own genomes (custom_library.fasta
#                                  or library/) + seqid2taxid.map, normalize the map to bare
#                                  accession->taxid, and tag the BLAST subjects with kraken2's
#                                  EXACT taxids so validation is taxid-vs-taxid (no NCBI taxdb).
# No `from __future__` import (Snakemake prepends a preamble — it must be the first statement).
import os
import sys

from metagx import validation as v

snk = snakemake  # noqa: F821  (injected by Snakemake)
build_from = snk.params.build_from
kraken_db = snk.params.kraken_db
prefix = snk.params.prefix
os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)

taxid_map = None
if build_from == "classifier":
    try:
        src = v.kraken2_db_sources(kraken_db)
    except ValueError as e:
        sys.exit(str(e))
    fastas = src["fastas"]
    source = fastas[0] if len(fastas) == 1 else os.path.dirname(fastas[0])
    if src["seqid2taxid"]:
        taxid_map = prefix + ".acc2taxid.tsv"
        mapping = v.normalize_seqid2taxid(src["seqid2taxid"])
        with open(taxid_map, "w") as fh:
            for acc, tx in mapping.items():
                fh.write(f"{acc}\t{tx}\n")
else:
    source = build_from

res = v.build_blast_db(source, prefix, run=True, taxid_map=taxid_map)
if not res.get("ok"):
    sys.exit(f"makeblastdb failed for {source!r}: {res.get('note') or res.get('tail') or res}")
open(snk.output.sentinel, "w").write(res.get("skipped", "built"))
