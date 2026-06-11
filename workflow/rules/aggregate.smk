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


# Convert each kraken2 report to Krona text (lineage from the report's own indentation) — needs
# no taxonomy DB, so it works with custom kraken2 DBs whose synthetic taxids aren't in NCBI.
rule kreport2krona_text:
    input:
        kreport=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kreport",
    output:
        txt=f"{OUT}/report/krona/{{sample}}.krona.txt",
    script:
        "../scripts/kreport2krona.py"


rule krona:
    input:
        texts=expand(f"{OUT}/report/krona/{{sample}}.krona.txt", sample=list(SAMPLES)),
    output:
        html=f"{OUT}/report/krona.html",
    conda:
        "../envs/aggregate.yaml"
    params:
        args=lambda wc: " ".join(registry.render_args(
            "krona", config.get("krona", {}), managed={"output": f"{OUT}/report/krona.html"})),
        # "<file>,<sample>" makes each sample a labelled, switchable dataset in the chart
        inputs=lambda wc, input: " ".join(
            f"{t},{os.path.basename(t)[:-len('.krona.txt')]}" for t in input.texts),
    shell:
        "ktImportText {params.args} {params.inputs}"
