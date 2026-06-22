# BLAST validation of classifier calls (the `validate` module).
#
# A kraken2/Bracken assignment is a k-mer call, not proof. For each WGS sample this BLASTs a
# seeded subsample of the reads assigned to the top taxa against db.blast (or NCBI -remote)
# and reports whether the best alignment's organism agrees with the classifier — a per-taxon
# agreement rate + an overall verdict. Depends on the read-level kraken2 outputs (report +
# per-read calls at the canonical sweep label) and the reads themselves.

VALIDATE = config.get("validate", {})


rule blast_validate:
    input:
        kreport=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kreport",
        kraken=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kraken",
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        json=f"{OUT}/validation/{{sample}}.blast_validation.json",
        tsv=f"{OUT}/validation/{{sample}}.blast_validation.tsv",
    threads: THREADS
    params:
        validate=VALIDATE,
        blastn=config.get("blastn", {}),
        blast_db=DB.get("blast", ""),
        sample=lambda wc: wc.sample,
    script:
        "../scripts/blast_validate.py"
