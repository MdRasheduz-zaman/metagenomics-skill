# metagx evidence files

YAML summaries of **validation runs** (platform sweeps, mock communities, etc.).
Consumed by `metagx recommend` and `metagx advise` — they do **not** change Snakemake
execution. Registries remain the source of truth for which flags the workflow accepts.

Add a new file when you complete an experiment whose metrics should guide future
interviews (e.g. `bracken_threshold.yaml` after a large-DB benchmark).

| File | Purpose |
|------|---------|
| `kraken2_confidence.yaml` | Platform confidence sweeps (validation runs) |
| `kraken2_secondary.yaml` | minimum_hit_groups, minimum_base_quality, quick |
| `bracken_read_length.yaml` | Bracken `-r` vs platform / median read length |
| `platform_routing.yaml` | QC/assembly/classify routing + install alternatives |
| `cutadapt_amplicon.yaml` | Amplicon minimum length after primer trim |
| `metaspades_memory.yaml` | metaSPAdes RAM cap (-m) |
| `kaiju_consensus.yaml` | Kaiju second-classifier sensitivity |
| `metaphlan_consensus.yaml` | MetaPhlAn stat_q for consensus |
| `phylogenetics_mafft.yaml` | MAFFT method by alignment set size |
| `phylogenetics_iqtree.yaml` | IQ-TREE bootstrap guidance |

Registry `recommend` / `warn_if` blocks in `metagx/parameters/*.yaml` are merged by
`metagx recommend` and `tool_advisor.recommend_config()`.
