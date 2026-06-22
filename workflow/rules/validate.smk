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


def _blast_source():
    if _BUILD_FROM == "classifier":
        return config.get("db", {}).get("build", {}).get("source", "")
    return _BUILD_FROM or ""


if _BUILD_FROM:
    rule build_validate_blast_db:
        # Build an in-scope BLAST DB from the classifier's own genomes (idempotent via makeblastdb).
        output:
            sentinel=f"{_INSYNC_PREFIX}.done",
        params:
            source=lambda wc: _blast_source(),
            prefix=_INSYNC_PREFIX,
        script:
            "../scripts/build_blast_db.py"


def _blast_db_for_validate():
    return _INSYNC_PREFIX if _BUILD_FROM else DB.get("blast", "")


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
        sample=lambda wc: wc.sample,
    script:
        "../scripts/blast_validate.py"
