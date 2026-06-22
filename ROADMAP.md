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

## Verified by real execution (2026-06-11, run 4)

- **Decontam** — ran the real `workflow/scripts/decontam.py` (`run()`) on a real multi-sample
  Bracken table and on a realistic blank-control scenario: control samples are dropped and
  flagged taxa removed while real biology is retained. Confirmed `run()`'s `n_flagged` matches
  the written `decontam_flagged.tsv`. Converts decontam from dry-run-only **and adds the first
  unit tests** (`tests/test_decontam.py`, 3) — it previously had none.
- **Consensus classifier (Kaiju), true end-to-end** — kraken2 (`metagx-bio`, viral_custom DB) +
  kaiju/kaiju2table (`metagx-kaiju`, kaiju_custom DB) on the same 20k real reads → cross-checked
  with `classifier_consensus.py`. **Found + fixed a real parser bug** the prior run's recorded
  0.94 had masked: `parse_kaiju_table` leaked two non-species rows — kaiju's DB-dependent
  `"cannot be assigned to a (non-viral) species"` (the literal guard only matched
  `"...to a species"`) and the 0-read `"Viruses"` (taxid 10239) higher-rank catch-all. Fix:
  prefix-match the "cannot be assigned" variants and drop zero-abundance rows. After the fix,
  kaiju species 32→30, **Jaccard 0.94→1.0**, no spurious `kaiju_only` (correct: both DBs are
  built from the same genomes). Regression test added (`tests/test_classifier_consensus.py`).
- **Reconcile (incl. CAT cross-check), true end-to-end** — ran the real `reconcile.py` `main()`
  on the bundled ont_sim assembly experiment (51 contigs, real kraken2 contig calls + jgi depth +
  samtools breadth + read kreport + **CAT `add_names`**). Exercised the CAT cross-check branch the
  plain ont_sim runs never had: 25 contigs classified, kraken2-vs-CAT 25 agree / 0 conflict;
  concordance both=19 / reads_only=11 / contigs_only=0; 13 flags; all output tables well-formed and
  internally consistent. No bug — but reconcile had **no unit tests**; added 4
  (`tests/test_reconcile.py`: contig-call/CAT/read-report parsers + a main() concordance/flag/
  cross-check integration).
- **Aggregate (MultiQC + Krona) — provisioned + run for real.** Built the `metagx-aggregate`
  conda env (osx-64; conda solved cleanly — recipe valid). **MultiQC 1.35** run on the real
  `results/multi` dir: found 6 kraken reports, wrote a 2.1 MB HTML. **Found + fixed a real defect
  in the Krona rule:** it used `ktImportTaxonomy -tax <db>`, which needs a Krona/NCBI taxonomy DB
  the custom kraken2 DBs' **synthetic taxids are absent from** — so a custom-DB run would render an
  empty chart (and it forced a large taxonomy download). Reworked to **kreport2krona →
  `ktImportText`** (new pure-Python `workflow/scripts/kreport2krona.py`, no KrakenTools dep; krona
  registry → `ktImportText`): the kreport's own indented lineage carries the taxonomy, so it needs
  no taxonomy DB and works with custom DBs. Verified end-to-end on the real custom-DB kreports
  (per-sample labelled datasets: sampA/B/C, 30/31/31 taxa, Powassan virus etc. in the chart) and
  the rewritten 2-rule DAG dry-runs clean. Tests updated/added (`tests/test_aggregate_args.py`:
  MultiQC + ktImportText tokens + kreport2krona converter).
- **`metagx run --use-conda` verified end-to-end (Snakemake's own provisioning).** Ran the real CLI
  on a 3-sample custom-DB config (`project: aggrun`): Snakemake built its own per-rule conda env
  (`.snakemake/conda/09a03e…`) and ran classify → abundance → aggregate, producing
  `results/aggrun/report/multiqc/multiqc_report.html` (2.1 MB, 71 kraken refs) + `report/krona.html`
  (custom-DB datasets embedded) + per-sample `report/krona/*.krona.txt`; re-invocation reports
  "Nothing to be done". Worked around the disk/conda limits with a dedicated conda-26 env on PATH +
  `CONDA_PKGS_DIRS` pointed at the warm package cache (hardlink, no re-download). Disk (100%-full
  volume) was the recurring obstacle throughout — freed room by removing two unrelated user conda
  envs (with approval) + cleaning caches.
- **DX fix shipped: conda-frontend preflight (`metagx/runner.py`).** The finding above — modern
  Snakemake (8/9) hard-requires **conda ≥ 24.7.1**, so on the box's 23.10.0 `--use-conda` died with
  a cryptic `CreateCondaEnvironmentException` — is now handled. The runner **prefers mamba** when
  present (`pick_conda_frontend`, no version gate, faster) and otherwise **preflights the conda
  version** (`conda_preflight`), raising `CondaFrontendError` *before* launching Snakemake with an
  actionable message ("conda >= 24.7.1 … update it or install mamba"). `cmd_run` prints it cleanly
  and exits 1; a plain `--dry-run` is never blocked (it provisions nothing). Verified on the stock
  machine (conda 23.10, no mamba → clean error, exit 1) + 7 unit tests (`tests/test_runner_preflight.py`).
- **Read↔contig accuracy, true end-to-end** — ran the real `read_contig_accuracy.py` `main()`
  (samtools from `metagx-bio`) on the bundled ont_sim BAM (14,806 aligned reads): **93.3%
  concordant**, 0.04% discordant, 6.66% read-unclassified on classified contigs; the nodes.dmp
  lineage branch fired. Extracted the read/contig bucketing into a pure `bucket()` helper
  (behavior-preserving — real numbers identical pre/post) so the scientific core is testable
  without a BAM; added 4 unit tests (`tests/test_read_contig_accuracy.py`: parsers, lineage walk,
  all five bucket categories, nodes.dmp).
- **Ancient-DNA damage authentication** — pure-Python; ran the real `damage_authenticate.py`
  `run()` on real-format mapDamage frequency fixtures: authentic case (terminal C→T 0.28 / G→A
  0.24) → `damage_present: true`; modern case (≤0.015) → false. Logic requires **both** ends
  elevated (conservative double-stranded signature). Added 3 unit tests (`tests/test_damage.py`).
  Real mapDamage execution still needs the ancient/damage conda env on real aDNA.
- Suite 148→166 (decontam ×3, kaiju-parser ×1, reconcile ×4, aggregate/kreport2krona ×3,
  read-accuracy ×4, damage ×3). Seven modules went from no/low-tests to tested; **two real bugs
  fixed** (kaiju consensus parser; Krona custom-DB taxonomy).

## HPC schedulers (2026-06-10)

`metagx run --executor {local,slurm,lsf,sge,pbs,generic}` via bundled `workflow/profiles/`,
one registry in `metagx/schedulers.py`, exposed on CLI + MCP + HTTP. slurm/lsf use native
Snakemake plugins; sge/pbs/generic use cluster-generic with site-editable submit commands.

## Open / next pass

- **In-scope BLAST validation from a *standard/prebuilt* kraken2 index.** `validate.build_from`
  already builds the BLAST DB from the classifier's genomes for `db.build` custom-fasta/folder/
  spike-in (and any user FASTA). For a *standard* build or a `fetch-db` prebuilt index, auto-derive
  the in-scope BLAST DB from the index's `library/*/library.fna` (when retained). Until then, point
  `validate.build_from` at the genome FASTA(s) you classified against. See
  `docs/ARCHITECTURE-WIRING.md` Part 2.
