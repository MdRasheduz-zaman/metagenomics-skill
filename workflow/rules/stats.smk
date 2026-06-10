# Cross-sample community statistics (Tier 1): normalization, alpha/beta diversity, PCoA.
# Pure-Python (numpy) — no external tools. Needs modules.abundance (the combined Bracken
# table) and >=2 samples for meaningful beta/ordination (alpha works with one).

rule diversity:
    input:
        bracken=f"{OUT}/summary/bracken_combined.tsv",
    output:
        alpha=f"{OUT}/stats/alpha_diversity.tsv",
        beta=f"{OUT}/stats/beta_braycurtis.tsv",
        jaccard=f"{OUT}/stats/beta_jaccard.tsv",
        pcoa=f"{OUT}/stats/pcoa.tsv",
        json=f"{OUT}/stats/diversity.json",
        barplot=f"{OUT}/stats/composition_barplot.png",
        pcoa_png=f"{OUT}/stats/pcoa.png",
        rarefaction=f"{OUT}/stats/rarefaction.tsv",
        rarefaction_png=f"{OUT}/stats/rarefaction.png",
        core=f"{OUT}/stats/core_taxa.tsv",
    params:
        label=READ_LABEL,
        outdir=f"{OUT}/stats",
        core_prevalence=config.get("stats", {}).get("core_prevalence", 0.8),
    script:
        "../scripts/diversity.py"
