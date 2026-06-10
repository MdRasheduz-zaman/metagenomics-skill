# Domain-aware contig/genome taxonomy & quality, gated by modules.domain_taxonomy + the
# `domains` list. Real metagenomes are mixed: classify reads broadly, then route the
# assembly to domain-appropriate tools.
#   viral       -> geNomad (identify + ICTV taxonomy) + CheckV (completeness)
#   prokaryote  -> GTDB-Tk (taxonomy) + CheckM2 (quality) on MetaBAT2 bins
#   eukaryote   -> EukRep (separate euk contigs) + EukCC (completeness)
# Each rule declares an isolated conda env (workflow/envs/*.yaml); run with --use-conda to
# have Snakemake provision the tools automatically.

def contigs_of(wc):
    return f"{OUT}/assembly/{wc.sample}/final.contigs.fa"


# ----------------------- viral -----------------------
rule genomad:
    input:
        contigs=lambda wc: contigs_of(wc),
    output:
        done=f"{OUT}/viral/{{sample}}/genomad.done",
    threads: THREADS
    conda:
        "../envs/viral.yaml"
    params:
        outdir=lambda wc: f"{OUT}/viral/{wc.sample}/genomad",
        db=DB.get("genomad", ""),
        args=lambda wc, threads: " ".join(
            registry.render_args("genomad", config.get("genomad", {}), managed={"threads": threads})
        ),
    shell:
        "genomad end-to-end {params.args} {input.contigs} {params.outdir} {params.db} "
        "&& touch {output.done}"


rule checkv:
    input:
        contigs=lambda wc: contigs_of(wc),
    output:
        summary=f"{OUT}/viral/{{sample}}/checkv/quality_summary.tsv",
    threads: THREADS
    conda:
        "../envs/viral.yaml"
    params:
        outdir=lambda wc: f"{OUT}/viral/{wc.sample}/checkv",
        db=DB.get("checkv", ""),
    shell:
        "checkv end_to_end {input.contigs} {params.outdir} -t {threads} -d {params.db}"


# ----------------------- prokaryote (needs bins) -----------------------
rule gtdbtk:
    input:
        bins=f"{OUT}/binning/{{sample}}/bins.done",
    output:
        done=f"{OUT}/prok/{{sample}}/gtdbtk.done",
    threads: THREADS
    conda:
        "../envs/prok.yaml"
    params:
        bindir=lambda wc: f"{OUT}/binning/{wc.sample}/bins",
        outdir=lambda wc: f"{OUT}/prok/{wc.sample}/gtdbtk",
        db=DB.get("gtdbtk", ""),
    shell:
        "GTDBTK_DATA_PATH={params.db} gtdbtk classify_wf --genome_dir {params.bindir} "
        "--out_dir {params.outdir} --cpus {threads} --extension fa --skip_ani_screen "
        "&& touch {output.done}"


rule checkm2:
    input:
        bins=f"{OUT}/binning/{{sample}}/bins.done",
    output:
        report=f"{OUT}/prok/{{sample}}/checkm2/quality_report.tsv",
    threads: THREADS
    conda:
        "../envs/prok.yaml"
    params:
        bindir=lambda wc: f"{OUT}/binning/{wc.sample}/bins",
        outdir=lambda wc: f"{OUT}/prok/{wc.sample}/checkm2",
        db=DB.get("checkm2", ""),
    shell:
        "checkm2 predict --threads {threads} --input {params.bindir} "
        "--output-directory {params.outdir} -x fa --database_path {params.db} --force"


# ----------------------- eukaryote -----------------------
rule eukrep:
    input:
        contigs=lambda wc: contigs_of(wc),
    output:
        euk=f"{OUT}/euk/{{sample}}/euk_contigs.fa",
    conda:
        "../envs/euk.yaml"
    params:
        args=lambda wc: " ".join(registry.render_args("eukrep", config.get("eukrep", {}))),
    shell:
        "EukRep -i {input.contigs} -o {output.euk} {params.args}"


rule eukcc:
    input:
        euk=f"{OUT}/euk/{{sample}}/euk_contigs.fa",
    output:
        csv=f"{OUT}/euk/{{sample}}/eukcc/eukcc.csv",
    threads: THREADS
    conda:
        "../envs/euk.yaml"
    params:
        outdir=lambda wc: f"{OUT}/euk/{wc.sample}/eukcc",
        db=DB.get("eukcc", ""),
    shell:
        "eukcc single --out {params.outdir} --threads {threads} --db {params.db} {input.euk}"
