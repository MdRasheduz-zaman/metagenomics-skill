# Binning refinement (Tier 2), gated by modules.bin_refinement (requires modules.binning).
# Runs two more binners on the same contigs+coverage, reconciles them into a consensus MAG
# set with DAS_Tool, then dereplicates across all samples with dRep:
#   MetaBAT2 (from binning) + MaxBin2 + CONCOCT  --DAS_Tool-->  refined per-sample bins
#   refined bins of every sample                  --dRep-->     study-level genome catalog
# WGS-only. Isolated conda envs (binrefine.yaml / drep.yaml); run with --use-conda.

# ----------------------- second binner: MaxBin2 -----------------------
rule maxbin2:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
        depth=f"{OUT}/binning/{{sample}}/depth.txt",
    output:
        done=f"{OUT}/binning/{{sample}}/maxbin/maxbin.done",
    threads: THREADS
    conda:
        "../envs/binrefine.yaml"
    params:
        outdir=lambda wc: f"{OUT}/binning/{wc.sample}/maxbin",
        abund=lambda wc: f"{OUT}/binning/{wc.sample}/maxbin/abund.txt",
        args=lambda wc, input, threads: " ".join(registry.render_args(
            "maxbin2", config.get("maxbin2", {}),
            managed={"contig": input.contigs,
                     "abund": f"{OUT}/binning/{wc.sample}/maxbin/abund.txt",
                     "out": f"{OUT}/binning/{wc.sample}/maxbin/bin",
                     "thread": threads})),
    shell:
        "mkdir -p {params.outdir} && "
        "awk 'NR>1{{print $1\"\\t\"$3}}' {input.depth} > {params.abund} && "
        "run_MaxBin.pl {params.args} && touch {output.done}"


# ----------------------- third binner: CONCOCT -----------------------
rule concoct:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
        bam=f"{OUT}/binning/{{sample}}/aln.sorted.bam",
    output:
        done=f"{OUT}/binning/{{sample}}/concoct/concoct.done",
    threads: THREADS
    conda:
        "../envs/binrefine.yaml"
    params:
        outdir=lambda wc: f"{OUT}/binning/{wc.sample}/concoct",
        args=lambda wc, threads: " ".join(registry.render_args(
            "concoct", config.get("concoct", {}), managed={"threads": threads})),
    shell:
        "mkdir -p {params.outdir}/bins && cd {params.outdir} && "
        "cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b contigs_10K.bed "
        "> contigs_10K.fa && "
        "concoct_coverage_table.py contigs_10K.bed {input.bam} > coverage_table.tsv && "
        "concoct {params.args} --composition_file contigs_10K.fa "
        "--coverage_file coverage_table.tsv -b ./ && "
        "merge_cutup_clustering.py clustering_gt1000.csv > clustering_merged.csv && "
        "extract_fasta_bins.py {input.contigs} clustering_merged.csv --output_path bins && "
        "touch {output.done}"


# ----------------------- consensus selection: DAS_Tool -----------------------
def _dastool_cmd(wc, threads):
    base = f"{OUT}/binning/{wc.sample}"
    refined = f"{base}/refined"
    args = registry.render_args(
        "das_tool", config.get("das_tool", {}),
        managed={"bins": f"{refined}/metabat.tsv,{refined}/maxbin.tsv,{refined}/concoct.tsv",
                 "labels": "metabat,maxbin,concoct",
                 "contigs": f"{OUT}/assembly/{wc.sample}/final.contigs.fa",
                 "outputbasename": f"{refined}/dastool",
                 "write_bins": True, "threads": threads},
    )
    return (
        f"mkdir -p {refined} && "
        f"Fasta_to_Contig2Bin.sh -i {base}/bins -e fa > {refined}/metabat.tsv && "
        f"Fasta_to_Contig2Bin.sh -i {base}/maxbin -e fasta > {refined}/maxbin.tsv && "
        f"Fasta_to_Contig2Bin.sh -i {base}/concoct/bins -e fa > {refined}/concoct.tsv && "
        "DAS_Tool " + " ".join(args)
    )


rule das_tool:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
        metabat=f"{OUT}/binning/{{sample}}/bins.done",
        maxbin=f"{OUT}/binning/{{sample}}/maxbin/maxbin.done",
        concoct=f"{OUT}/binning/{{sample}}/concoct/concoct.done",
    output:
        done=f"{OUT}/binning/{{sample}}/refined/dastool.done",
    threads: THREADS
    conda:
        "../envs/binrefine.yaml"
    params:
        cmd=lambda wc, threads: _dastool_cmd(wc, threads),
    shell:
        "{params.cmd} && touch {output.done}"


# ----------------------- cross-sample dereplication: dRep -----------------------
rule drep:
    # Run-level: dereplicate the refined bins of every sample into one representative set.
    input:
        refined=expand(f"{OUT}/binning/{{sample}}/refined/dastool.done", sample=WGS_SAMPLES),
    output:
        done=f"{OUT}/binning/drep/drep.done",
    threads: THREADS
    conda:
        "../envs/drep.yaml"
    params:
        outdir=f"{OUT}/binning/drep/dereplicated",
        bins=" ".join(
            f"{OUT}/binning/{s}/refined/dastool_DASTool_bins/*.fa" for s in WGS_SAMPLES
        ),
        args=lambda wc, threads: " ".join(registry.render_args(
            "drep", config.get("drep", {}), managed={"processors": threads})),
    shell:
        "dRep dereplicate {params.outdir} {params.args} -g {params.bins} && "
        "touch {output.done}"
