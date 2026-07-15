# BLAST validation of classifier calls (the `validate` module).
#
# A kraken2/Bracken assignment is a k-mer call, not proof. For each WGS sample this BLASTs a
# seeded subsample of the reads assigned to the top taxa against a reference and reports whether
# the best alignment's organism agrees with the classifier — a per-taxon agreement rate + an
# overall verdict. The reference is one of:
#   * db.blast            — a local BLAST+ nucleotide DB the user provides;
#   * validate.build_from — build the BLAST DB from the SAME genomes as the classifier (FASTA/
#                           folder, or "classifier" => reuse the db.build source) so the
#                           benchmark is IN SCOPE (the right design — see SKILL.md); or
#   * validate.remote     — search NCBI remotely (only for a few sequences).

VALIDATE = config.get("validate", {})
_BUILD_FROM = VALIDATE.get("build_from")
_INSYNC_PREFIX = f"{OUT}/validation/blastdb/insync"


if _BUILD_FROM:
    rule build_validate_blast_db:
        # Build an in-scope BLAST DB from the classifier's own genomes (idempotent via makeblastdb).
        # For build_from: classifier, also tag subjects with kraken2's taxids (taxid-vs-taxid).
        output:
            sentinel=f"{_INSYNC_PREFIX}.done",
        params:
            build_from=_BUILD_FROM,
            kraken_db=DB.get("kraken2", ""),
            prefix=_INSYNC_PREFIX,
        script:
            "../scripts/build_blast_db.py"


def _blast_db_for_validate():
    return _INSYNC_PREFIX if _BUILD_FROM else DB.get("blast", "")


def _names_dmp():
    # kraken2's own names.dmp -> in-sync taxid->name resolution for the agreement check
    cand = os.path.join(DB.get("kraken2", ""), "taxonomy", "names.dmp") if DB.get("kraken2") else ""
    return cand if cand and os.path.isfile(cand) else ""


def _nodes_dmp():
    # kraken2's own nodes.dmp -> roll a read's leaf assignment up to the validation rank (so a
    # genus clade collects its species-assigned reads, and agreement is taxid-vs-taxid)
    cand = os.path.join(DB.get("kraken2", ""), "taxonomy", "nodes.dmp") if DB.get("kraken2") else ""
    return cand if cand and os.path.isfile(cand) else ""


rule blast_validate:
    input:
        kreport=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kreport",
        kraken=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kraken",
        reads=lambda wc: reads_for_classify(wc.sample),
        # when building an in-sync DB, wait for it; otherwise no extra dep
        blastdb_ready=lambda wc: ([f"{_INSYNC_PREFIX}.done"] if _BUILD_FROM else []),
    output:
        json=f"{OUT}/validation/{{sample}}.blast_validation.json",
        tsv=f"{OUT}/validation/{{sample}}.blast_validation.tsv",
    threads: THREADS
    params:
        validate=VALIDATE,
        blastn=config.get("blastn", {}),
        blast_db=_blast_db_for_validate(),
        names_dmp=_names_dmp(),
        nodes_dmp=_nodes_dmp(),
        sample=lambda wc: wc.sample,
    script:
        "../scripts/blast_validate.py"
