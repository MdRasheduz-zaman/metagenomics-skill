# Optional de novo assembly + genome binning, dispatched by platform.
#   short reads -> MEGAHIT ;  long reads (ONT/PacBio) -> Flye/metaFlye
# Both emit {OUT}/assembly/{sample}/final.contigs.fa, so binning is platform-agnostic.
# Read mapping uses a platform-appropriate minimap2 preset (sr / map-ont / map-hifi / map-pb).

# ---------- user-provided contigs: skip assembly, stage the FASTA ----------
# When a sample supplies pre-assembled contigs (an isolate genome, a previous assembly, or
# downloaded references), copy them to the assembler's output path so every downstream module
# (AMR/functional, BGC, binning, domain taxonomy, reconcile) consumes them unchanged.
rule stage_provided_contigs:
    wildcard_constraints:
        sample=_alt(PROVIDED_CONTIGS),
    input:
        contigs=lambda wc: provided_contigs(wc.sample),
    output:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    shell:
        "mkdir -p $(dirname {output.contigs}) && cp {input.contigs} {output.contigs}"


# ---------- short-read assembly: MEGAHIT ----------
def _megahit_cmd(wc, output, threads):
    reads = reads_for_classify(wc.sample)
    args = registry.render_args("megahit", dict(config.get("megahit", {})),
                                managed={"threads": threads})
    ins = f"-1 {reads[0]} -2 {reads[1]}" if is_paired(wc.sample) else f"-r {reads[0]}"
    tmp = f"{output.contigs}.megahit_tmp"
    return (f"rm -rf {tmp} && megahit {ins} -o {tmp} " + " ".join(args) +
            f" && cp {tmp}/final.contigs.fa {output.contigs} && rm -rf {tmp}")


rule megahit:
    wildcard_constraints:
        sample=_alt(SHORT_ASM),
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    threads: THREADS
    params:
        cmd=lambda wc, output, threads: _megahit_cmd(wc, output, threads),
    shell:
        "{params.cmd}"


# ---------- short-read assembly: metaSPAdes (optional, paired-end / hybrid) ----------
def _metaspades_cmd(wc, output, threads):
    reads = reads_for_classify(wc.sample)  # paired-end [R1, R2] (validated in common.smk)
    tmp = f"{output.contigs}.spades_tmp"
    managed = {"meta": True, "pe1": reads[0], "pe2": reads[1], "out": tmp, "threads": threads}
    lr = long_reads_of(wc.sample)
    if lr:  # hybrid: fold long reads into the SPAdes assembly
        lp = SAMPLES[wc.sample].get("long_platform", "ont")
        managed["pacbio" if lp.startswith("pacbio") else "nanopore"] = lr
    args = registry.render_args("metaspades", dict(config.get("metaspades", {})), managed=managed)
    return (f"rm -rf {tmp} && spades.py " + " ".join(args) +
            f" && cp {tmp}/contigs.fasta {output.contigs} && rm -rf {tmp}")


rule metaspades:
    wildcard_constraints:
        sample=_alt(SPADES_ASM),
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    threads: THREADS
    conda:
        "../envs/spades.yaml"
    params:
        cmd=lambda wc, output, threads: _metaspades_cmd(wc, output, threads),
    shell:
        "{params.cmd}"


# ---------- long-read assembly: Flye / metaFlye ----------
def _flye_cmd(wc, output, threads):
    reads = reads_for_classify(wc.sample)[0]
    tmp = f"{output.contigs}.flye_tmp"
    args = registry.render_args("flye", dict(config.get("flye", {})),
                                managed={"threads": threads, "out_dir": tmp})
    return (f"rm -rf {tmp} && flye {flye_read_flag(wc.sample)} {reads} " + " ".join(args) +
            f" && cp {tmp}/assembly.fasta {output.contigs} && rm -rf {tmp}")


rule flye:
    wildcard_constraints:
        sample=_alt(LONG_ASM),
    input:
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    threads: THREADS
    params:
        cmd=lambda wc, output, threads: _flye_cmd(wc, output, threads),
    shell:
        "{params.cmd}"


# ---------- mapping + binning (platform-agnostic; preset varies) ----------
rule map_to_contigs:
    input:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
        reads=lambda wc: reads_for_classify(wc.sample),
    output:
        bam=f"{OUT}/binning/{{sample}}/aln.sorted.bam",
    threads: THREADS
    params:
        preset=lambda wc: minimap2_preset(wc.sample),
    shell:
        "minimap2 -ax {params.preset} -t {threads} {input.contigs} {input.reads} "
        "| samtools sort -@ {threads} -o {output.bam} - && samtools index {output.bam}"


rule contig_depth:
    input:
        bam=f"{OUT}/binning/{{sample}}/aln.sorted.bam",
    output:
        depth=f"{OUT}/binning/{{sample}}/depth.txt",
    shell:
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}"


rule metabat2:
    input:
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
        depth=f"{OUT}/binning/{{sample}}/depth.txt",
    output:
        marker=f"{OUT}/binning/{{sample}}/bins.done",
    threads: THREADS
    params:
        bindir=lambda wc: f"{OUT}/binning/{wc.sample}/bins",
        prefix=lambda wc: f"{OUT}/binning/{wc.sample}/bins/bin",
        log=lambda wc: f"{OUT}/binning/{wc.sample}/metabat2.log",
        args=lambda wc, threads: " ".join(
            registry.render_args("metabat2", config.get("metabat2", {}), managed={"threads": threads})
        ),
    shell:
        # Recovering zero MAGs is a legitimate outcome for low-coverage / low-diversity /
        # single-sample inputs (MetaBAT2 exits non-zero with "no large target contigs").
        # Don't let that crash the whole run: capture MetaBAT2's output, count the bins
        # actually written, and record the count in the marker so zero is visible, not hidden.
        # `find` (not a glob) so the count is 0 cleanly under `set -euo pipefail` when no
        # bins are written — a bare `ls *.fa` would error on no-match and abort the rule.
        "mkdir -p {params.bindir} && "
        "( metabat2 -i {input.contigs} -a {input.depth} -o {params.prefix} {params.args} "
        "  > {params.log} 2>&1 || true ) && "
        "n=$(find {params.bindir} -maxdepth 1 -name '*.fa' | wc -l | tr -d ' ') && "
        "echo \"bins=$n\" > {output.marker} && "
        "echo \"MetaBAT2: recovered $n bin(s) — see {params.log}\""
