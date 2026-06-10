# Run-level reporting/visualization (Tier 2), gated by modules.aggregate (requires
# modules.classify). Two run-level artifacts that summarize the whole study:
#   MultiQC -> one interactive HTML aggregating fastp + kraken2 (read-flow accounting)
#   Krona   -> one interactive zoomable taxonomy chart across all samples
# Composition barplot + PCoA ordination already come from modules.stats. Run with --use-conda.

rule multiqc:
    input:
        # depend on per-sample classification matrices so MultiQC runs after the reports exist
        matrices=expand(f"{OUT}/summary/{{sample}}.matrix.json", sample=list(SAMPLES)),
    output:
        html=f"{OUT}/report/multiqc/multiqc_report.html",
    conda:
        "../envs/aggregate.yaml"
    params:
        scandir=OUT,
        args=lambda wc: " ".join(registry.render_args(
            "multiqc", config.get("multiqc", {}),
            managed={"outdir": f"{OUT}/report/multiqc", "filename": "multiqc_report",
                     "force": True})),
    shell:
        "multiqc {params.scandir} {params.args}"


rule krona:
    input:
        reports=expand(f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kreport", sample=list(SAMPLES)),
    output:
        html=f"{OUT}/report/krona.html",
    conda:
        "../envs/aggregate.yaml"
    params:
        # taxid in column 5, magnitude (reads) in column 3 of the kraken2 report
        args=lambda wc: " ".join(registry.render_args(
            "krona", config.get("krona", {}),
            managed={"tax_field": 5, "magnitude_field": 3,
                     "output": f"{OUT}/report/krona.html", "taxonomy": DB.get("krona", "")})),
    shell:
        "ktImportTaxonomy {params.args} {input.reports}"
