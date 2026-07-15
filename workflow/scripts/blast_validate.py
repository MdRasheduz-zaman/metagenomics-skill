# BLAST validation of read-level classifier calls (the `validate` module).
#
# For a sample: take the top taxa from the read kraken2 report, pull a seeded subsample of the
# reads kraken2 assigned to each, BLAST them against db.blast (or NCBI -remote), and check
# whether the best alignment's organism agrees with the classifier's call. Emits a per-taxon
# agreement table + an overall concordance verdict. No `from __future__` import (Snakemake
# prepends a preamble — it must stay the first statement).
import json
import os
import random
import subprocess

from metagx import registry
from metagx import validation as validate

snk = snakemake  # noqa: F821  (injected by Snakemake)
cfg = dict(snk.params.validate)
blastn_cfg = dict(snk.params.blastn or {})
sample = snk.params.sample
threads = int(snk.threads)
db = snk.params.blast_db          # local BLAST db path/prefix (or remote db name when remote)
names_dmp = getattr(snk.params, "names_dmp", "")  # kraken2 names.dmp -> in-sync taxid->name
remote = bool(cfg.get("remote"))
seed = int(cfg.get("seed", 42))
rng = random.Random(seed)

os.makedirs(os.path.dirname(snk.output.json), exist_ok=True)
workdir = os.path.dirname(snk.output.json)

# 1) top taxa + per-read assignments
level = cfg["level"]                   # 'genus' | 'species' — the rank NAME for nodes.dmp roll-up
taxa = validate.top_taxa(snk.input.kreport, level=cfg["rank"], top_n=int(cfg["top_n"]))
with open(snk.input.kraken) as fh:
    assignments = validate.parse_kraken_assignments(fh.read())

# kraken2's own taxonomy tree (nodes.dmp) lets us roll each read's LEAF assignment up to the
# validation rank, so a genus clade collects the reads assigned at species beneath it (and the
# agreement can be taxid-vs-taxid, robust to non-binomial names). Absent -> group by the leaf,
# which is correct when reads are assigned at this rank (e.g. species-level validation).
nodes_dmp = getattr(snk.params, "nodes_dmp", "")
tree = validate.parse_nodes_dmp(nodes_dmp) if nodes_dmp and os.path.isfile(nodes_dmp) else {}

by_taxid = {}
for rid, tid in assignments.items():
    key = (validate.ancestor_at_rank(tid, level, tree) or tid) if tree else tid
    by_taxid.setdefault(key, []).append(rid)

query_taxon = {}      # read id -> classifier taxon NAME (name-fallback agreement)
query_taxid = {}      # read id -> classifier clade TAXID at this rank (taxid agreement)
n_per = int(cfg["reads_per_taxon"])
for t in taxa:
    ids = by_taxid.get(t["taxid"], [])
    if len(ids) > n_per:
        ids = rng.sample(ids, n_per)
    for rid in ids:
        query_taxon[rid] = t["name"]
        query_taxid[rid] = t["taxid"]

# 2) pull the sequences and write the query FASTA
seqs = validate.extract_sequences([str(r) for r in snk.input.reads], set(query_taxon))
query_fa = os.path.join(workdir, f"{sample}.validate_query.fasta")
n_written = validate.write_fasta({rid: seqs[rid] for rid in query_taxon if rid in seqs}, query_fa)
query_taxon = {rid: name for rid, name in query_taxon.items() if rid in seqs}  # only BLASTed ones

blast_out = os.path.join(workdir, f"{sample}.validate_blast.tsv")
outfmt = "6 " + " ".join(validate.BLAST6_FIELDS)

# 3) run blastn (subprocess argv — keeps the multi-token -outfmt intact, no shell quoting)
if n_written == 0:
    open(blast_out, "w").close()
    if not taxa:
        # No taxa at the requested rank — the usual cause is a synthetic-taxonomy DB (flat,
        # species-only) validated at genus. Name the mismatch so it's fixable, not opaque.
        present = validate.ranks_present(snk.input.kreport)
        blast_note = (f"no taxa at rank '{cfg['rank']}' in the kreport (ranks present: "
                      f"{', '.join(present) or 'none'}). If this is a synthetic-taxonomy DB, set "
                      f"validate.level: species (rank S); otherwise check validate.rank.")
    else:
        blast_note = "no classified reads to validate (top taxa had no assigned reads)"
else:
    base = dict(blastn_cfg)
    managed = {"query": query_fa, "out": blast_out, "outfmt": outfmt}
    if remote:
        base["remote"] = True
        managed["db"] = db or "nt"          # remote db name; -num_threads is illegal with -remote
    else:
        managed["db"] = db
        managed["num_threads"] = threads
    argv = ["blastn"] + registry.render_args("blastn", base, managed=managed)
    proc = subprocess.run(argv, capture_output=True, text=True)
    blast_note = ("ok" if proc.returncode == 0
                  else f"blastn exited {proc.returncode}: {(proc.stderr or '')[-500:]}")

# 4) parse hits, assess agreement per taxon, summarize
with open(blast_out) as fh:
    hits = validate.best_hit_per_query(validate.parse_blast6(fh.read()))

# in-sync name resolution: when the BLAST DB carries kraken2's taxids, resolve them through
# kraken2's OWN names.dmp so the comparison uses the same taxonomy (no NCBI taxdb needed).
name_resolver = None
if names_dmp and os.path.isfile(names_dmp):
    _names = validate.parse_names_dmp(names_dmp)
    name_resolver = _names.get

per_taxon = []
rows_tsv = ["taxon\ttaxid\tn_queries\tn_with_hits\tn_agree\thit_rate\tagreement_rate\tverdict"]
for t in taxa:
    qt = {rid: name for rid, name in query_taxon.items() if name == t["name"]}
    if not qt:
        continue
    a = validate.assess(qt, {q: hits[q] for q in qt if q in hits}, level=level,
                        name_resolver=name_resolver,
                        query_taxids={rid: query_taxid[rid] for rid in qt}, tree=tree)
    a.update(taxon=t["name"], taxid=t["taxid"], reads_classified=t["reads"],
             verdict=validate.verdict(a["agreement_rate"], a["hit_rate"]))
    per_taxon.append(a)
    rows_tsv.append(f"{t['name']}\t{t['taxid']}\t{a['n_queries']}\t{a['n_with_hits']}\t"
                    f"{a['n_agree']}\t{a['hit_rate']}\t{a['agreement_rate']}\t{a['verdict']}")

summary = validate.summarize(sample, per_taxon, level=cfg["level"],
                             db=("NCBI-remote:" + (db or "nt")) if remote else db)
summary["blast_note"] = blast_note
summary["n_query_reads"] = n_written

with open(snk.output.json, "w") as fh:
    json.dump(summary, fh, indent=2)
with open(snk.output.tsv, "w") as fh:
    fh.write("\n".join(rows_tsv) + "\n")
