# Differential abundance between two sample groups (gated by modules.differential; requires
# abundance + a `group` column with two labels, >=2 samples each). CLR + permutation test +
# BH FDR over the combined Bracken table. Pure-Python (numpy) — no external tool.

rule differential:
    input:
        bracken=f"{OUT}/summary/bracken_combined.tsv",
    output:
        tsv=f"{OUT}/stats/differential_abundance.tsv",
        json=f"{OUT}/stats/differential_abundance.json",
        png=f"{OUT}/stats/differential_volcano.png",
    params:
        label=READ_LABEL,
        groups=GROUPS,
        differential=config.get("differential", {}),
    script:
        "../scripts/differential.py"
