# Read QC, dispatched by platform. All paths emit gzipped FASTQ under {OUT}/qc/.
#   short reads (illumina/mgi):  fastp           (pe / se / interleaved)
#   ONT:                         porechop_abi -> chopper
#   PacBio:                      chopper         (adapters assumed handled upstream)
# Command lines are built from the per-tool registries.

# ---------- host/contaminant removal (pre-classification) ----------
# Map reads to a host genome and keep the UNmapped reads (-f 4). First-class QC step so the
# whole pipeline (classify, assembly, ...) runs on host-depleted reads. SE for now.
rule host_remove:
    wildcard_constraints:
        sample=_alt([s for s in SAMPLES if layout(s) == "se"]) if HOST_GENOME else r"(?!)",
    input:
        reads=lambda wc: raw_reads(wc.sample)[0],
    output:
        fq=f"{OUT}/hostclean/{{sample}}.fastq",
    threads: THREADS
    params:
        host=HOST_GENOME or "",
        preset=lambda wc: minimap2_preset(wc.sample),
    shell:
        "minimap2 -ax {params.preset} -t {threads} {params.host} {input.reads} "
        "| samtools fastq -n -f 4 - > {output.fq}"


# ---------- short reads: fastp ----------
rule fastp_pe:
    wildcard_constraints:
        sample=_alt(SHORT_PE_FASTQ),
    input:
        r1=lambda wc: staged_reads(wc.sample)[0],
        r2=lambda wc: staged_reads(wc.sample)[1],
    output:
        r1=f"{OUT}/qc/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/qc/{{sample}}_R2.fastq.gz",
        json=f"{OUT}/qc/{{sample}}.fastp.json",
        html=f"{OUT}/qc/{{sample}}.fastp.html",
    threads: THREADS
    params:
        args=lambda wc, input, output, threads: " ".join(
            registry.render_args("fastp", config.get("fastp", {}), managed={
                "thread": threads, "in1": input.r1, "in2": input.r2,
                "out1": output.r1, "out2": output.r2, "json": output.json, "html": output.html,
            })
        ),
    shell:
        "fastp {params.args}"


rule fastp_se:
    wildcard_constraints:
        sample=_alt(SHORT_SE_FASTQ),
    input:
        r1=lambda wc: staged_reads(wc.sample)[0],
    output:
        r1=f"{OUT}/qc/{{sample}}_R1.fastq.gz",
        json=f"{OUT}/qc/{{sample}}.fastp.json",
        html=f"{OUT}/qc/{{sample}}.fastp.html",
    threads: THREADS
    params:
        args=lambda wc, input, output, threads: " ".join(
            registry.render_args("fastp", config.get("fastp", {}), managed={
                "thread": threads, "in1": input.r1,
                "out1": output.r1, "json": output.json, "html": output.html,
            })
        ),
    shell:
        "fastp {params.args}"


rule fastp_interleaved:
    wildcard_constraints:
        sample=_alt(SHORT_IL_FASTQ),
    input:
        r1=lambda wc: staged_reads(wc.sample)[0],
    output:
        r1=f"{OUT}/qc/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/qc/{{sample}}_R2.fastq.gz",
        json=f"{OUT}/qc/{{sample}}.fastp.json",
        html=f"{OUT}/qc/{{sample}}.fastp.html",
    threads: THREADS
    params:
        args=lambda wc, input, output, threads: " ".join(
            registry.render_args("fastp", config.get("fastp", {}), managed={
                "thread": threads, "interleaved_in": True, "in1": input.r1,
                "out1": output.r1, "out2": output.r2, "json": output.json, "html": output.html,
            })
        ),
    shell:
        "fastp {params.args}"


# ---------- ancient DNA: fastp read-merging (collapse overlapping pairs -> single reads) ----------
rule fastp_ancient:
    wildcard_constraints:
        sample=_alt(ANCIENT_MERGE),
    input:
        r1=lambda wc: staged_reads(wc.sample)[0],
        r2=lambda wc: staged_reads(wc.sample)[1],
    output:
        r1=f"{OUT}/qc/{{sample}}_R1.fastq.gz",
        json=f"{OUT}/qc/{{sample}}.fastp.json",
        html=f"{OUT}/qc/{{sample}}.fastp.html",
    threads: THREADS
    params:
        args=lambda wc, input, output, threads: " ".join(
            registry.render_args("fastp", config.get("fastp", {}), managed={
                "thread": threads, "in1": input.r1, "in2": input.r2,
                "merge": True, "merged_out": output.r1,
                "json": output.json, "html": output.html,
            })
        ),
    shell:
        "fastp {params.args}"


# ---------- ONT: porechop_abi -> chopper ----------
rule ont_qc:
    wildcard_constraints:
        sample=_alt(ONT_FASTQ),
    input:
        r1=lambda wc: staged_reads(wc.sample)[0],
    output:
        r1=f"{OUT}/qc/{{sample}}_R1.fastq.gz",
    threads: THREADS
    params:
        chop=lambda wc, input, threads: " ".join(
            registry.render_args("porechop_abi", config.get("porechop_abi", {}),
                                  managed={"input": input.r1, "output": f"{OUT}/qc/{wc.sample}.porechop.fastq",
                                           "threads": threads})
        ),
        chopper=lambda wc, threads: " ".join(
            registry.render_args("chopper", config.get("chopper", {}), managed={"threads": threads})
        ),
        tmp=lambda wc: f"{OUT}/qc/{wc.sample}.porechop.fastq",
    shell:
        "porechop_abi {params.chop} && "
        "chopper {params.chopper} < {params.tmp} | gzip > {output.r1} && "
        "rm -f {params.tmp}"


# ---------- amplicon: cutadapt primer removal (replaces fastp for marker-gene reads) ----------
_AMP = config.get("amplicon", {})


rule cutadapt_pe:
    wildcard_constraints:
        sample=_alt(AMPLICON_PE),
    input:
        r1=lambda wc: staged_reads(wc.sample)[0],
        r2=lambda wc: staged_reads(wc.sample)[1],
    output:
        r1=f"{OUT}/qc/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/qc/{{sample}}_R2.fastq.gz",
    threads: THREADS
    conda:
        "../envs/amplicon.yaml"
    params:
        args=lambda wc, input, output, threads: " ".join(
            registry.render_args("cutadapt", config.get("cutadapt", {}), managed={
                "threads": threads, "fwd_primer": _AMP.get("fwd_primer", ""),
                "rev_primer": _AMP.get("rev_primer", "")})),
    shell:
        "cutadapt {params.args} -o {output.r1} -p {output.r2} {input.r1} {input.r2}"


rule cutadapt_se:
    wildcard_constraints:
        sample=_alt(AMPLICON_SE),
    input:
        r1=lambda wc: staged_reads(wc.sample)[0],
    output:
        r1=f"{OUT}/qc/{{sample}}_R1.fastq.gz",
    threads: THREADS
    conda:
        "../envs/amplicon.yaml"
    params:
        args=lambda wc, input, output, threads: " ".join(
            registry.render_args("cutadapt", config.get("cutadapt", {}), managed={
                "threads": threads, "fwd_primer": _AMP.get("fwd_primer", "")})),
    shell:
        "cutadapt {params.args} -o {output.r1} {input.r1}"


# ---------- PacBio: chopper only ----------
rule pacbio_qc:
    wildcard_constraints:
        sample=_alt(PB_FASTQ),
    input:
        r1=lambda wc: staged_reads(wc.sample)[0],
    output:
        r1=f"{OUT}/qc/{{sample}}_R1.fastq.gz",
    threads: THREADS
    params:
        chopper=lambda wc, threads: " ".join(
            registry.render_args("chopper", config.get("chopper", {}), managed={"threads": threads})
        ),
        decompress=lambda wc, input: cat_cmd(input.r1),
    shell:
        "{params.decompress} | chopper {params.chopper} | gzip > {output.r1}"
