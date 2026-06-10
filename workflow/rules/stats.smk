# Cross-sample community statistics (Tier 1): normalization, alpha/beta diversity, PCoA.
# Pure-Python (numpy) — no external tools. Needs modules.abundance (the combined Bracken
# table) and >=2 samples for meaningful beta/ordination (alpha works with one).

rule diversity:
    input:
        bracken=f"{OUT}/summary/bracken_combined.tsv",
    output:
        alpha=f"{OUT}/stats/alpha_diversity.tsv",
        beta=f"{OUT}/stats/beta_braycurtis.tsv",
        pcoa=f"{OUT}/stats/pcoa.tsv",
        json=f"{OUT}/stats/diversity.json",
        barplot=f"{OUT}/stats/composition_barplot.png",
        pcoa_png=f"{OUT}/stats/pcoa.png",
    params:
        label=READ_LABEL,
        outdir=f"{OUT}/stats",
    script:
        "../scripts/diversity.py"
