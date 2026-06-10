# Amplicon (marker-gene) analysis branch. Assembly does NOT apply to amplicon data, so these
# samples are routed here instead: primer removal happens in QC (cutadapt, qc.smk); then
#   short reads (Illumina/MGI) -> VSEARCH OTU table (merge/filter/derep/cluster/chimera)
#   long reads  (ONT/PacBio)   -> Emu relative-abundance (EM against a 16S DB; needs db.emu)
# Read-level kraken2/Bracken classification still runs on amplicon samples (allowed, with a
# warning that a marker-gene DB + these methods are more appropriate). Conda env: amplicon.yaml.


def _vsearch_cmd(wc, output, threads):
    reads = reads_for_classify(wc.sample)              # cutadapt-trimmed reads
    d = os.path.dirname(output.table)
    v = config.get("vsearch", {})
    idv = v.get("cluster_id", 0.97)
    minsize = v.get("minsize", 2)
    merged = f"{d}/merged.fastq"
    filt = f"{d}/filtered.fasta"
    derep = f"{d}/derep.fasta"
    raw = f"{d}/otus_raw.fasta"
    if is_paired(wc.sample):
        merge = f"vsearch --fastq_mergepairs {reads[0]} --reverse {reads[1]} --fastqout {merged} --threads {threads}"
    else:
        merge = f"gunzip -c {reads[0]} > {merged} 2>/dev/null || cp {reads[0]} {merged}"
    return (
        f"mkdir -p {d} && {merge} && "
        f"vsearch --fastq_filter {merged} --fastaout {filt} --fastq_maxee 1.0 && "
        f"vsearch --derep_fulllength {filt} --output {derep} --sizeout --minuniquesize {minsize} && "
        f"vsearch --cluster_size {derep} --id {idv} --centroids {raw} --sizein --sizeout --threads {threads} && "
        f"vsearch --uchime3_denovo {raw} --nonchimeras {output.otus} && "
        f"vsearch --usearch_global {filt} --db {output.otus} --id {idv} --otutabout {output.table} --threads {threads}"
    )


rule vsearch_otus:
    wildcard_constraints:
        sample=_alt(SHORT_AMPLICON_OTU),
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        otus=f"{OUT}/amplicon/{{sample}}/otus.fasta",
        table=f"{OUT}/amplicon/{{sample}}/otu_table.txt",
    threads: THREADS
    conda:
        "../envs/amplicon.yaml"
    params:
        cmd=lambda wc, output, threads: _vsearch_cmd(wc, output, threads),
    shell:
        "{params.cmd}"


# DADA2 ASV inference (amplicon.method: asv) — exact sequence variants instead of 97% OTUs.
# R/Bioconductor; user params come from the dada2 registry, I/O is injected here.
# workflow.basedir is the main Snakefile's dir (workflow/); scripts live in workflow/scripts/.
_DADA2_R = os.path.join(workflow.basedir, "scripts", "dada2_asv.R")


rule dada2_asv:
    wildcard_constraints:
        sample=_alt(SHORT_AMPLICON_ASV),
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        table=f"{OUT}/amplicon/{{sample}}/asv_table.tsv",
        seqs=f"{OUT}/amplicon/{{sample}}/asv_seqs.fasta",
    threads: THREADS
    conda:
        "../envs/dada2.yaml"
    params:
        script=_DADA2_R,
        r2=lambda wc: f"--r2 {reads_for_classify(wc.sample)[1]}" if is_paired(wc.sample) else "",
        args=lambda wc: " ".join(registry.render_args("dada2", config.get("dada2", {}))),
    shell:
        "Rscript {params.script} --r1 {input.reads[0]} {params.r2} "
        "--out_table {output.table} --out_seqs {output.seqs} "
        "--sample {wildcards.sample} --threads {threads} {params.args}"


def _emu_type(sample):
    return "map-ont" if platform_of(sample) in ONT_PLATFORMS else "map-pb"


rule emu_abundance:
    wildcard_constraints:
        sample=_alt(LONG_AMPLICON),
    input:
        reads=lambda wc: reads_for_classify(wc.sample)[0],
    output:
        tsv=f"{OUT}/amplicon/{{sample}}/emu_rel-abundance.tsv",
    threads: THREADS
    conda:
        "../envs/amplicon.yaml"
    params:
        d=lambda wc: f"{OUT}/amplicon/{wc.sample}",
        db=DB.get("emu", ""),
        type=lambda wc: _emu_type(wc.sample),
        args=lambda wc, threads: " ".join(
            registry.render_args("emu", config.get("emu", {}), managed={"threads": threads})),
    shell:
        "emu abundance --type {params.type} --db {params.db} {params.args} "
        "--output-dir {params.d} {input.reads} && "
        "mv {params.d}/*rel-abundance.tsv {output.tsv}"
