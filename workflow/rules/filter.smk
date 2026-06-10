# Taxonomic read filtering before assembly + filtered-vs-unfiltered comparison.
# Default is DEPLETION (remove host/contaminant taxids, keep everything else incl.
# unclassified). Optional mapping-based host depletion (minimap2, gold standard) runs first
# when read_filter.host_genome is set. Single-end only for now. Needs modules.assembly
# (unfiltered baseline) + modules.classify (per-read kraken calls).

from metagx import readfilter as _readfilter

_RF = config.get("read_filter", {})
_NODES = os.path.join(DB["kraken2"], "taxonomy", "nodes.dmp")
_SE_FASTA = [s for s in FASTA_SAMPLES if layout(s) == "se"]
_SE_FASTQ = [s for s in FASTQ_SAMPLES if layout(s) == "se"]


def prefilter_reads(sample):
    """Reads entering the taxid filter: host-depleted (if a host genome is set) else staged."""
    if _RF.get("host_genome"):
        return f"{OUT}/filtered/{sample}.hostdepleted.fastq"
    return staged_reads(sample)[0]


def filtered_reads(sample):
    return f"{OUT}/filtered/{sample}.filtered{formats.canonical_ext(SAMPLES[sample]['r1'])}"


# Optional gold-standard host depletion: map to host genome, keep the UNmapped reads.
rule host_deplete:
    input:
        reads=lambda wc: staged_reads(wc.sample)[0],
    output:
        fq=f"{OUT}/filtered/{{sample}}.hostdepleted.fastq",
    threads: THREADS
    params:
        host=_RF.get("host_genome", ""),
        preset=lambda wc: minimap2_preset(wc.sample),
    shell:
        "minimap2 -ax {params.preset} -t {threads} {params.host} {input.reads} "
        "| samtools fastq -n -f 4 - > {output.fq}"


def _do_filter(wc, input, output):
    return _readfilter.filter_reads(
        input.reads, input.kraken, output[0], _RF.get("taxids", []),
        nodes_path=_NODES, include_children=_RF.get("include_children", True),
        mode=_RF.get("mode", "exclude"), keep_unclassified=_RF.get("keep_unclassified", True))


rule kraken_filter_fasta:
    wildcard_constraints:
        sample=_alt(_SE_FASTA),
    input:
        reads=lambda wc: prefilter_reads(wc.sample),
        kraken=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kraken",
    output:
        f"{OUT}/filtered/{{sample}}.filtered.fasta",
    run:
        stats = _do_filter(wildcards, input, output)
        print(f"[read_filter] {wildcards.sample}: kept {stats['kept']}/{stats['total']} "
              f"({stats['mode']}, removed {stats['removed']})")


rule kraken_filter_fastq:
    wildcard_constraints:
        sample=_alt(_SE_FASTQ),
    input:
        reads=lambda wc: prefilter_reads(wc.sample),
        kraken=f"{OUT}/kraken2/{{sample}}.{READ_LABEL}.kraken",
    output:
        f"{OUT}/filtered/{{sample}}.filtered.fastq",
    run:
        stats = _do_filter(wildcards, input, output)
        print(f"[read_filter] {wildcards.sample}: kept {stats['kept']}/{stats['total']} "
              f"({stats['mode']}, removed {stats['removed']})")


def _filtered_asm_cmd(wc, output, threads):
    reads = filtered_reads(wc.sample)
    tmp = f"{output.contigs}.tmp"
    if assembler(wc.sample) == "flye":
        args = registry.render_args("flye", dict(config.get("flye", {})),
                                    managed={"threads": threads, "out_dir": tmp})
        return (f"rm -rf {tmp} && flye {flye_read_flag(wc.sample)} {reads} " + " ".join(args) +
                f" && cp {tmp}/assembly.fasta {output.contigs} && rm -rf {tmp}")
    args = registry.render_args("megahit", dict(config.get("megahit", {})),
                                managed={"threads": threads})
    return (f"rm -rf {tmp} && megahit -r {reads} -o {tmp} " + " ".join(args) +
            f" && cp {tmp}/final.contigs.fa {output.contigs} && rm -rf {tmp}")


rule assemble_filtered:
    input:
        reads=lambda wc: filtered_reads(wc.sample),
    output:
        contigs=f"{OUT}/filtered_assembly/{{sample}}/final.contigs.fa",
    threads: THREADS
    params:
        cmd=lambda wc, output, threads: _filtered_asm_cmd(wc, output, threads),
    shell:
        "{params.cmd}"


rule compare_assemblies:
    input:
        unfiltered=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
        filtered=f"{OUT}/filtered_assembly/{{sample}}/final.contigs.fa",
    output:
        tsv=f"{OUT}/filtered_assembly/{{sample}}.assembly_comparison.tsv",
        json=f"{OUT}/filtered_assembly/{{sample}}.assembly_comparison.json",
    params:
        sample=lambda wc: wc.sample,
    script:
        "../scripts/compare_assemblies.py"
