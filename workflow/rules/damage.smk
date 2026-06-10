# Ancient-DNA damage authentication (Tier 3), gated by modules.damage (requires assembly +
# ≥1 ancient sample). Maps the merged ancient reads to the assembly (reuses map_to_contigs)
# and quantifies post-mortem cytosine deamination with mapDamage2, then emits a verdict:
# authentic aDNA shows elevated C→T at 5' ends and G→A at 3' ends. Run with --use-conda.

rule mapdamage:
    wildcard_constraints:
        sample=_alt(ANCIENT_SAMPLES),
    input:
        bam=f"{OUT}/binning/{{sample}}/aln.sorted.bam",
        contigs=f"{OUT}/assembly/{{sample}}/final.contigs.fa",
    output:
        ct5=f"{OUT}/ancient/{{sample}}/mapdamage/5pCtoT_freq.txt",
        ga3=f"{OUT}/ancient/{{sample}}/mapdamage/3pGtoA_freq.txt",
    threads: THREADS
    conda:
        "../envs/ancient.yaml"
    params:
        args=lambda wc, input: " ".join(registry.render_args(
            "mapdamage", config.get("mapdamage", {}),
            managed={"input": input.bam, "reference": input.contigs,
                     "folder": f"{OUT}/ancient/{wc.sample}/mapdamage"})),
    shell:
        "mapDamage {params.args}"


rule damage_authenticate:
    wildcard_constraints:
        sample=_alt(ANCIENT_SAMPLES),
    input:
        ct5=f"{OUT}/ancient/{{sample}}/mapdamage/5pCtoT_freq.txt",
        ga3=f"{OUT}/ancient/{{sample}}/mapdamage/3pGtoA_freq.txt",
    output:
        json=f"{OUT}/ancient/{{sample}}/authentication.json",
    params:
        sample=lambda wc: wc.sample,
        threshold=config.get("damage_ct_threshold", 0.05),
    script:
        "../scripts/damage_authenticate.py"
