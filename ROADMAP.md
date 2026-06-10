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
- `metagx history`, `sync-help`, `catalog`, `schedulers`

**Incomplete (extend with new validation runs):** per-platform MAFFT method from real alignment benchmarks, binning depth thresholds, functional module defaults.

## Verified by real execution (2026-06-10, run 3)

Beyond dry-run, the following were executed against real tools on the dev machine and
caught real bugs that dry-run cannot:

- **Phylogenetics** — MAFFT (linsi) → TrimAl → FastTree on a 6-taxon fixture using the
  bundled `.conda/phylogenetics` env (real alignment + Newick tree + plot; A/B clades
  recovered). Fixed: `newick_stats` counted FastTree internal-node support values as leaves.
  Gated real-execution test added (skips in CI).
- **Consensus classifier (Kaiju)** — added `metagx build-kaiju-db` (custom Kaiju protein DB
  from genomes via prodigal + kaiju-mkbwt/mkfmi, synthetic taxids, no NCBI download), built it
  from the bundled genomes, ran Kaiju on the real ont_sim reads (9.7k classified), and the
  classifier_consensus cross-check vs kraken2: 30 shared species, Jaccard 0.94, full top-10
  overlap. Converts `modules.classify_consensus` (kaiju path) from dry-run-only to verified.
- **Functional / AMR (ABRicate)** — ran `abricate --db card` (real tool + bundled CARD/NCBI/
  ResFinder/VFDB DBs, installed osx-64) on the real ont_sim contigs (valid output, 0 hits =
  correct for viral) AND on an ErmA positive-control fixture (correctly detected the resistance
  gene). Gated real-run test added (`tests/test_functional_amr.py`, skips without ABRicate).
- **Core profiling + cross-sample stats** — fastp → kraken2 → Bracken → `modules.stats` on
  the viral DB + 4 samples (`metagx-bio` env). The new diversity outputs (Chao1/ACE/Good's
  coverage, analytic rarefaction, Jaccard, core microbiome) and the diversity-aware advisor
  were verified on real Bracken output. Fixed: the post-run advisor crashed (`KeyError:
  'qc_key'`) for any **TSV sample sheet** — the common case — because `sample_contexts`
  returned an incomplete fallback; it now parses the sheet and uses real per-sample platforms.

## HPC schedulers (2026-06-10)

`metagx run --executor {local,slurm,lsf,sge,pbs,generic}` via bundled `workflow/profiles/`,
one registry in `metagx/schedulers.py`, exposed on CLI + MCP + HTTP. slurm/lsf use native
Snakemake plugins; sge/pbs/generic use cluster-generic with site-editable submit commands.
