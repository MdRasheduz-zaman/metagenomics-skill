# Phylogenetics (Tier 2): MAFFT → optional TrimAl → IQ-TREE 2 or FastTree.
# Gated by modules.phylogenetics. Requires phylogenetics.input (FASTA) or
# phylogenetics.aligned_input (skip alignment). Outputs under {OUT}/phylogenetics/.
# Run with --use-conda (workflow/envs/phylogenetics.yaml).

_PHYLO = config.get("phylogenetics", {})
_PHYLO_INPUT = _PHYLO.get("input")
if MODULES.get("phylogenetics") and not (_PHYLO.get("input") or _PHYLO.get("aligned_input")):
    raise ValueError(
        "modules.phylogenetics is on but phylogenetics.input (or aligned_input) is missing"
    )


rule phylogenetics:
    input:
        fasta=_PHYLO_INPUT or _PHYLO.get("aligned_input"),
    output:
        aligned=f"{OUT}/phylogenetics/aligned.fasta",
        tree=f"{OUT}/phylogenetics/tree.nwk",
        json=f"{OUT}/phylogenetics/phylogenetics.json",
        plot=f"{OUT}/phylogenetics/tree.png",
    threads: THREADS
    conda:
        "../envs/phylogenetics.yaml"
    params:
        phylo=_PHYLO,
        mafft=config.get("mafft", {}),
        iqtree=config.get("iqtree", {}),
        fasttree=config.get("fasttree", {}),
        threads=THREADS,
    script:
        "../scripts/phylogenetics.py"
