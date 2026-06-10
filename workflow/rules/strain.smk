# Strain-level microdiversity (Tier 3), gated by modules.strain (requires assembly). Profiles
# within-population SNVs / nucleotide diversity (π) from the reads-vs-contigs mapping with
# inStrain — resolving mixed strains that abundance profiling cannot. WGS-only. --use-conda.

rule instrain:
    wildcard_constraints:
        sample=_alt(WGS_SAMPLES),
    input:
        bam=f"{OUT}/binning/{{sample}}/aln.sorted.bam",
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    output:
        done=f"{OUT}/strain/{{sample}}/instrain.done",
    threads: THREADS
    conda:
        "../envs/strain.yaml"
    params:
        args=lambda wc, threads: " ".join(registry.render_args(
            "instrain", config.get("instrain", {}),
            managed={"output": f"{OUT}/strain/{wc.sample}/IS", "processors": threads})),
    shell:
        "inStrain profile {input.bam} {input.contigs} {params.args} && touch {output.done}"
