# Decontamination via negative/blank controls (Tier 3), gated by modules.decontam (requires
# abundance + ≥1 control sample, marked `control: true` in the sample sheet). Prevalence-based:
# taxa as/more prevalent in controls than in real samples are flagged as contaminants and
# removed from a cleaned abundance table. Pure-Python — no external tool.

rule decontam:
    input:
        bracken=f"{OUT}/summary/bracken_combined.tsv",
    output:
        flagged=f"{OUT}/stats/decontam_flagged.tsv",
        cleaned=f"{OUT}/stats/abundance_decontaminated.tsv",
    params:
        controls=CONTROL_SAMPLES,
        label=READ_LABEL,
    script:
        "../scripts/decontam.py"
