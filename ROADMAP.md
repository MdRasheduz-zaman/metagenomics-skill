# metagx roadmap

## Phylogenetics module (shipped)

`modules.phylogenetics` follows the same pattern as other modules:

1. **Input** — `phylogenetics.input` (unaligned FASTA) or `phylogenetics.aligned_input`.
2. **Align** — MAFFT (`mafft` registry; method presets in workflow script).
3. **Trim** — optional TrimAl (`phylogenetics.trim: true`).
4. **Tree** — IQ-TREE 2 (default) or FastTree (`phylogenetics.method: iqtree|fasttree|auto`).
5. **Outputs** — `results/<project>/phylogenetics/` (aligned FASTA, Newick, JSON, PNG).

Evidence: `phylogenetics_mafft.yaml`, `phylogenetics_iqtree.yaml`. Advisor includes
`mafft`, `iqtree`, `fasttree` when the module is enabled.

## Advisor layer (shipped)

- `metagx/evidence/` — kraken2 sweeps, bracken read length, platform routing, secondary kraken2 params, cutadapt/metaspades/kaiju/metaphlan guidance
- `metagx/tool_advisor.py` — multi-tool `recommend_config()` (QC routing, optional modules)
- `metagx recommend --config` — all enabled tools; `--tool` for single-param mode
- Registry `recommend` / `warn_if` on kraken2, bracken, fastp, chopper, megahit, metabat2, porechop_abi, cutadapt, metaspades, kaiju, metaphlan, mafft, iqtree, …
- `metagx advise` — post-run metrics + recommendations → `next_config.suggested.yaml`
- `metagx history`, `sync-help`, `catalog`

**Incomplete (extend with new validation runs):** per-platform MAFFT method from real alignment benchmarks, binning depth thresholds, functional module defaults.
