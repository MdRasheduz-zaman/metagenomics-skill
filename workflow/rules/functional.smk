# Functional layer (Tier 2), gated by modules.functional. Data-driven sub-steps run only
# when their inputs exist, so one flag scales with how much of the pipeline is enabled:
#   pathways    -> HUMAnN on QC'd shotgun reads        (read-based; no assembly needed)
#   amr         -> AMRFinderPlus + ABRicate on contigs (needs modules.assembly)
#   annotation  -> Bakta gene-calling + eggNOG-mapper on bins (needs modules.binning)
# WGS-only (functional profiling does not apply to amplicon/marker-gene data). Each sub-step
# declares its own isolated conda env; DBs come from config["db"]. Run with --use-conda.

# ----------------------- read-based pathways: HUMAnN -----------------------
rule humann_concat_reads:
    # HUMAnN consumes a single sequence file; merge pe/interleaved QC'd reads into one.
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        merged=temp(f"{OUT}/functional/{{sample}}/humann/input.fastq.gz"),
    shell:
        "cat {input.reads} > {output.merged}"


rule humann:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        reads=f"{OUT}/functional/{{sample}}/humann/input.fastq.gz",
    output:
        path=f"{OUT}/functional/{{sample}}/humann/input_pathabundance.tsv",
    threads: THREADS
    conda:
        "../envs/functional.yaml"
    params:
        args=lambda wc, input, threads: " ".join(registry.render_args(
            "humann", config.get("humann", {}),
            managed={"input": input.reads,
                     "output": f"{OUT}/functional/{wc.sample}/humann",
                     "nucleotide_database": DB.get("humann_nucleotide", ""),
                     "protein_database": DB.get("humann_protein", ""),
                     "threads": threads})),
    shell:
        "humann {params.args}"


# ----------------------- contig AMR: AMRFinderPlus + ABRicate -----------------------
rule amrfinder:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    output:
        tsv=f"{OUT}/functional/{{sample}}/amr/amrfinder.tsv",
    threads: THREADS
    conda:
        "../envs/amr.yaml"
    params:
        args=lambda wc, input, output, threads: " ".join(registry.render_args(
            "amrfinderplus", config.get("amrfinderplus", {}),
            managed={"nucleotide": input.contigs, "output": output.tsv,
                     "database": DB.get("amrfinderplus", ""), "threads": threads})),
    shell:
        "amrfinder {params.args}"


rule abricate:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    output:
        tsv=f"{OUT}/functional/{{sample}}/amr/abricate.tsv",
    threads: THREADS
    conda:
        "../envs/amr.yaml"
    params:
        args=lambda wc, threads: " ".join(registry.render_args(
            "abricate", config.get("abricate", {}), managed={"threads": threads})),
    shell:
        "abricate {params.args} {input.contigs} > {output.tsv}"


# ----------------------- MAG annotation: Bakta + eggNOG-mapper -----------------------
rule bakta_bins:
    # Annotate the recovered bins together (concatenated) — gene calling + features.
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        bins=f"{OUT}/binning/{{sample}}/bins.done",
    output:
        faa=f"{OUT}/functional/{{sample}}/annotate/bakta/bins.faa",
    threads: THREADS
    conda:
        "../envs/annotate.yaml"
    params:
        bindir=lambda wc: f"{OUT}/binning/{wc.sample}/bins",
        outdir=lambda wc: f"{OUT}/functional/{wc.sample}/annotate/bakta",
        catfa=lambda wc: f"{OUT}/functional/{wc.sample}/annotate/bakta/all_bins.fa",
        args=lambda wc, threads: " ".join(registry.render_args(
            "bakta", config.get("bakta", {}),
            managed={"db": DB.get("bakta", ""),
                     "output": f"{OUT}/functional/{wc.sample}/annotate/bakta",
                     "prefix": "bins", "threads": threads})),
    shell:
        "mkdir -p {params.outdir} && cat {params.bindir}/*.fa > {params.catfa} && "
        "bakta {params.args} {params.catfa}"


rule eggnog:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        faa=f"{OUT}/functional/{{sample}}/annotate/bakta/bins.faa",
    output:
        annot=f"{OUT}/functional/{{sample}}/annotate/eggnog/bins.emapper.annotations",
    threads: THREADS
    conda:
        "../envs/annotate.yaml"
    params:
        outdir=lambda wc: f"{OUT}/functional/{wc.sample}/annotate/eggnog",
        args=lambda wc, input, threads: " ".join(registry.render_args(
            "eggnog", config.get("eggnog", {}),
            managed={"input": input.faa, "output": "bins",
                     "output_dir": f"{OUT}/functional/{wc.sample}/annotate/eggnog",
                     "data_dir": DB.get("eggnog", ""), "cpu": threads})),
    shell:
        "mkdir -p {params.outdir} && emapper.py {params.args}"
