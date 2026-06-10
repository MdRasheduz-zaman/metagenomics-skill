---
name: metagx
description: >
  Interview-driven metagenomics pipeline. Use when the user wants to run taxonomic
  classification or a metagenomics workflow — kraken2, Bracken, fastp QC, a kraken2
  confidence sweep/matrix, or MAG assembly/binning (MEGAHIT/MetaBAT2). Triggers:
  "classify reads", "run kraken2", "metagenomics", "taxonomic profiling", "confidence
  sweep", ".fastq", ".fastq.gz", ".fasta", ".kreport", "Bracken abundance". The skill interviews the user,
  generates a validated config.yaml from per-tool parameter registries, and runs a
  Snakemake workflow.
---

# metagx — metagenomics pipeline skill

A schema-driven pipeline. **Per-tool parameter registries** (`metagx/parameters/*.yaml`)
are the single source of truth that drives the interview, config validation, the MCP
tool schemas, the CLI, and the actual command-line construction in the Snakefile. Add a
flag to a registry once and it is available everywhere — no other file needs editing.

## When to use
The user has sequencing reads and wants taxonomic profiling, a kraken2 confidence
comparison ("which threshold should I use?"), abundance estimates, or assembly/binning.

**Inputs handled:** FASTA or FASTQ (±gzip); single-end, paired-end (R1/R2), or
interleaved; **short reads** (Illumina/MGI) and **long reads** (ONT, PacBio). Format is
auto-detected; the per-sample `platform` drives QC and assembler choice. Ask the user
their platform — it changes the tools (see below).

## How to run it (Claude Code / CLI available)
Do **not** guess kraken2 flags. Drive the interview from the registries.

1. **One-time setup** (if not done): `bash setup.sh`, then `uv pip install -e .`. The bio
   tools install via conda/mamba (on Apple-Silicon macOS use
   `bash scripts/install_bio_macos_arm64.sh`).
2. **Custom database (optional).** If the user has reference genomes rather than a prebuilt
   index: `metagx build-db --genomes <genomes.fasta> --db <dir> --read-length <len>` builds
   a kraken2 + Bracken db with a synthetic taxonomy (no NCBI download).
3. **Offer a preset.** Run `metagx presets` and let the user pick a starting point
   (`pathogen-detection`, `gut-profiling`, `soil-deep-assembly`, `quick-screen`,
   `amr-surveillance`, `ancient-dna`). A preset
   pre-fills modules + parameters; the user then only adjusts what they care about. This
   replaces most of the interview. Pass the name as `preset` to `build-config`.
3. **Interview / refine.** For each enabled module, run `metagx interview <tool>` (e.g.
   `kraken2`, `fastp`, `bracken`). Ask the returned questions in plain language, one
   cluster at a time, showing defaults. Use `metagx params <tool>` for the full advanced
   flag list. Offer the **confidence sweep** as the headline feature.
   **Skip what the request already answers, but say so.** The interview is a funnel
   (goal/preset → platform → data → per-module params); if the user's request already fixes
   the goal, platform, or dataset, start at the first *unresolved* decision instead of
   re-asking. When you skip steps, name them in one line ("goal, platform, and dataset are set
   by your request — starting at the differential parameters") so the entry point is
   transparent.
4. **Write the config.** Assemble answers (plus any `preset`) into the `build_config`
   kwargs and either:
   - write an `answers.json` and run `metagx build-config answers.json`, or
   - construct `config.yaml` directly following `config/config.example.yaml`,
   then verify with `metagx validate config.yaml`.
5. **Dry run, then run:** `metagx run --config config.yaml --dry-run`, then drop
   `--dry-run`. (Database not present? `snakemake --snakefile workflow/Snakefile
   fetch_kraken_db` — only when the user agrees, the index is large.)
   After a successful run, metagx auto-writes `results/<project>/advisor/` and appends
   to `.metagx/history.jsonl` (use `--no-advisor` / `--no-history` to skip).
6. **Recommend before sweeping:** `metagx recommend --platform pacbio_clr` (or
   `--config config.yaml` to infer platform) returns evidence-based confidence grids
   and warnings from `metagx/evidence/`.
7. **Post-run advisor:** `metagx advise --config config.yaml --write` inspects finished
   results and writes `advisor.json`, `trial_log.md`, and optional
   `next_config.suggested.yaml`. Use `metagx history` to list prior trials.
8. **Report results:** `metagx results --config config.yaml` for the JSON summaries
   (`results/<project>/summary/`: `*.matrix.json`, `*.heatmap.png`,
   `bracken_combined.tsv`).
9. **Generate provenance + Methods.** `metagx report --config config.yaml` writes
   `results/<project>/report/`: `manifest.json` (tool versions, exact commands, DB
   identity, QC + % classified metrics), `methods.md` (a paste-into-paper Methods
   paragraph + citations — the "Copy as Methods" feature), and `report.md` (full report
   with parameter table, figures, abundance table). Add `--format latex|pdf` (needs
   pandoc) for those outputs.
10. **Write the full paper.** `metagx paper --config config.yaml` elaborates the whole run
   into a structured **IMRaD manuscript** (Introduction / Methods / Results / Discussion +
   abstract, tables, figures, references) and compiles it to PDF with **pdflatex** —
   `results/<project>/report/paper.{tex,pdf}`. Every number is read back from the result
   files; the design comes from the interview and the methods/citations from the registries,
   so it is a real first draft to refine, not fabricated. Use `--no-pdf` to emit only the
   `.tex` if no LaTeX engine is installed.

## Key concepts
- **Advisor layer (evidence, not codegen):** `metagx recommend --config config.yaml`
  returns **multi-tool** guidance: platform QC routing (fastp vs porechop/chopper vs
  cutadapt), Bracken `read_length` from median read length, kraken2 confidence +
  `minimum_hit_groups`, assembly defaults, optional modules to enable, and **alternatives
  to install** (e.g. trimmomatic, filtlong) when not wired in metagx. Single-param mode:
  `metagx recommend --tool bracken --param read_length --platform ont`. `metagx advise`
  merges post-run metrics with the same engine — including **diversity-aware advice** when
  `modules.stats` ran: it reads `stats/diversity.json` and flags under-sampling (mean Good's
  coverage < 0.95 → rarefaction not saturated, sequence deeper) and an empty core microbiome
  (high heterogeneity → check grouping or lower `stats.core_prevalence`). History:
  `.metagx/history.jsonl`. Does **not** replace Snakemake — it informs the next interview pass.
- **Phylogenetics:** `modules.phylogenetics` runs MAFFT → optional TrimAl → IQ-TREE 2 or
  FastTree on `phylogenetics.input` (or skip alignment with `aligned_input`). Registries:
  `mafft`, `iqtree`, `fasttree`. Outputs: `results/<project>/phylogenetics/` (aligned FASTA,
  Newick, JSON stats, tree figure). Use `method: auto` to pick FastTree when >500 sequences.
- **Confidence sweep ("k-dense" matrix):** `sweep: {param: confidence, values: [...]}`
  runs kraken2 at each value and produces a per-sample matrix + line plot showing how
  each organism's read count changes with threshold. Never also pin the swept param in
  the `kraken2:` section.
- **Modules:** `qc` (fastp/porechop+chopper) → `classify` (kraken2) → `abundance`
  (Bracken); optional `assembly` (MEGAHIT/metaSPAdes/Flye), `binning` (MetaBAT2),
  `bin_refinement`, `reconcile`, `domain_taxonomy`, `filtered_assembly`, `stats`,
  `differential`, `classify_consensus`, `functional`, `bgc`, `aggregate`, `strain`,
  `damage`, `decontam`, and `phylogenetics`. Toggle in `modules:`.
- **Reconcile** (`modules.reconcile`, needs `assembly`+`classify`): classifies the
  assembled contigs with kraken2 and joins them to per-contig coverage and the read-level
  calls. Outputs under `results/<project>/reconcile/`: a per-contig taxonomy table, a
  per-taxon reconciliation (read abundance vs coverage-weighted contig abundance +
  concordance: both / reads-only / contigs-only), and discordance/chimera flags. Use it to
  separate high-confidence calls (seen in both) from read-only (low-coverage/unassembled)
  and to spot misassembly/novel contigs. Do **not** Bracken contigs — contig classification
  is presence/identity, not abundance.
- **CAT cross-check** (optional): set `db.cat` and reconcile also classifies contigs with
  **CAT** (per-ORF voting, more chimera-robust than whole-contig LCA) and reports
  kraken2-vs-CAT agreement per contig. Build it with
  `metagx build-cat-db --genomes <fa> --db <dir> --taxonomy <kraken2_db>/taxonomy`.
  Bin-level taxonomy (BAT/GTDB-Tk) needs prokaryotic MAGs + ~100 GB refs — N/A to
  viral/unbinnable assemblies.
- **Domain taxonomy** (`modules.domain_taxonomy` + `domains: [...]`): real metagenomes are
  mixed, so reads are classified broadly once, then the assembly is routed per domain —
  **viral** → geNomad (identify+ICTV taxonomy) + CheckV (completeness); **prokaryote** →
  GTDB-Tk + CheckM2 on bins (needs `binning`); **eukaryote** → EukRep (separate euk contigs)
  + EukCC. Each needs its reference DB under `db.{genomad,checkv,gtdbtk,checkm2,eukcc}`. These
  tools live in isolated conda envs (`workflow/envs/`) — run with `metagx run --use-conda`
  to auto-provision them; the core tools come from `environment.yml`.
- **Filtered assembly** (`modules.filtered_assembly` + `read_filter`): taxonomically filter
  reads before assembly, then compare filtered vs unfiltered. Default is **depletion** —
  `mode: exclude` removes host/contaminant `taxids` and **keeps unclassified** (so novel
  reads aren't lost). `mode: include` is targeted recovery (still keep_unclassified to avoid
  losing divergent target reads + fragmenting the assembly). Optional `host_genome` adds
  gold-standard minimap2 host depletion. Abundance stays on UNfiltered reads. Output:
  `filtered_assembly/<sample>.assembly_comparison.{tsv,json}` (contigs/N50/bp delta) — so the
  benefit (or harm) of filtering is measured, not assumed. Uses our own read extractor (no
  KrakenTools dep). Single-end for now.
- **Library strategy (WGS vs amplicon).** Each sample has `library: wgs|amplicon` (default
  wgs). **Assembly only applies to WGS shotgun data** — amplicon (16S/18S/ITS) is many copies
  of one locus, so assembly/binning/reconcile/filtered_assembly/domain_taxonomy are *skipped*
  for amplicon samples. Instead amplicon is routed to: **cutadapt** primer removal (set
  `amplicon.fwd_primer`/`rev_primer`) → **VSEARCH** OTU table (short reads) or **Emu**
  relative abundance (long reads, needs `db.emu`). Read-level kraken2/Bracken still runs on
  amplicon (allowed) but the skill **warns** it's rough — prefer a marker-gene DB
  (SILVA/GreenGenes/UNITE) + the OTU/Emu output. An all-amplicon run with assembly modules on
  is rejected with a clear error.
- **Cross-sample stats** (`modules.stats`, needs `abundance` + ≥2 samples): α-diversity
  (Shannon/Simpson/richness/Pielou **+ Chao1 & ACE asymptotic richness estimators + Good's
  coverage**), β-diversity (Bray–Curtis abundance **+ Jaccard presence/absence**), PCoA
  ordination, TSS + CLR matrices, **analytic rarefaction curves** (`rarefaction.{tsv,png}` —
  "did I sequence deeply enough?", Hurlbert's expected richness, no random subsampling), and the
  **core microbiome** (`core_taxa.tsv` — taxa shared across ≥`stats.core_prevalence` of samples,
  default 0.8), plus composition barplot + PCoA plot under `results/<project>/stats/`. Pure-Python
  (numpy only), so it is fully unit-tested. Chao1/ACE/rarefaction are count-based (est. reads
  rounded to integers); the curve flattening toward observed richness means the depth was enough.
- **Differential abundance** (`modules.differential`, needs `abundance` + a `group` column):
  answers *which taxa differ between conditions?* Tests CLR-transformed Bracken counts between
  two groups (e.g. case/control) with a two-sided **permutation test** + Benjamini-Hochberg
  FDR (an ALDEx2-lite; pure-Python, no scipy/R). Outputs `stats/differential_abundance.tsv`
  (per-taxon group means, CLR diff, effect size, p, q, significant), a JSON summary, and a
  volcano plot. Mark samples with a `group` label in the sheet; needs ≥2 groups and **≥2
  samples per group to run** — but a 2-vs-2 permutation p-value floors at ~0.33, so use more
  replicates for real power. Optional `differential.reference_group` sets the contrast direction.
- **BGC mining** (`modules.bgc`, needs `assembly`; WGS-only): runs **antiSMASH** on the
  assembled contigs to find biosynthetic gene clusters (antibiotics, NRPS/PKS, siderophores).
  Tune `antismash: {taxon, cb_general, cb_knownclusters}`; reference DBs via `db.antismash`
  (`download-antismash-databases`). Output tree under `results/<project>/bgc/<sample>/`.
- **Amplicon ASV vs OTU** (`amplicon.method`): short-read marker-gene samples default to
  VSEARCH 97% **OTUs** (`method: otu`); set `method: asv` to denoise into exact **amplicon
  sequence variants with DADA2** (`dada2: {trunc_len_f, trunc_len_r, max_ee_f, ...}`) — ASVs
  are reproducible across studies and resolve single-base differences OTUs cannot. Long-read
  amplicon still routes to Emu.
- **Functional layer** (`modules.functional`, WGS-only): HUMAnN gene-family + MetaCyc pathway
  profiling on reads (always); AMRFinderPlus + ABRicate AMR/virulence on contigs (when
  `assembly` on); Bakta + eggNOG-mapper MAG annotation (when `binning` on). DBs:
  `db.{humann_nucleotide,humann_protein,amrfinderplus,bakta,eggnog}`. Tune e.g.
  `abricate: {db: card}`. Preset: `amr-surveillance`.
- **Assembler choice / hybrid** (`assembly: {assembler: metaspades}`): use metaSPAdes instead
  of MEGAHIT for short reads (paired-end only). Give a short-read sample a `long_reads`
  (+ `long_platform`) column for hybrid short+long assembly. Long-read samples always use Flye.
- **Bin refinement** (`modules.bin_refinement`, needs `binning`): MaxBin2 + CONCOCT alongside
  MetaBAT2 → DAS_Tool consensus per sample → dRep dereplication across samples into one
  representative genome catalog. Tune `das_tool.score_threshold`, `drep.{completeness,s_ani}`.
- **Consensus classifier** (`modules.classify_consensus`, needs `classify`, WGS-only): runs a
  second, independent classifier — `consensus: {classifier: metaphlan|kaiju}` — and writes a
  per-sample species concordance JSON vs kraken2 (`results/<project>/consensus/`). Agreement =
  confidence; kraken2-only taxa flag DB-completeness false positives. DB: `db.{metaphlan|kaiju}`.
- **Aggregate report** (`modules.aggregate`, needs `classify`): run-level MultiQC report +
  interactive Krona taxonomy chart under `results/<project>/report/`.
- **Ancient DNA** (`library: ancient` + `modules.damage`, needs `assembly`): short-read PE
  ancient samples are read-merged (fastp `--merge`, collapse overlapping fragments → single
  reads, treated as SE downstream), classified/assembled, then mapped back and run through
  mapDamage2. `ancient/<sample>/authentication.json` reports 5' C→T / 3' G→A deamination +
  a `damage_present` verdict — the test that separates authentic aDNA from modern
  contamination. Preset: `ancient-dna`. Tune `damage_ct_threshold` (default 0.05).
- **Decontam** (`modules.decontam`, needs `abundance` + a `control: true` sample): prevalence
  test over the combined Bracken table — taxa as/more prevalent in negative/blank controls
  than in real samples are flagged and removed (`stats/decontam_flagged.tsv`,
  `stats/abundance_decontaminated.tsv`). Pure-Python. For low-biomass samples.
- **Strain level** (`modules.strain`, needs `assembly`): inStrain SNV/microdiversity profiling
  over the reads-vs-contigs mapping, resolving mixed strains abundance can't. **Short-read
  (Illumina) only in practice** — inStrain's default `--min_read_ani 0.95` rejects most
  ONT/PacBio reads (~5–15% error), so the advisor warns when `strain` is on with long-read
  samples. Provision inStrain via `--use-conda` (bioconda), not pip — it pins an ancient
  biopython that won't build on modern arm64/Python.
- **Per-sample Bracken length**: kmer_distrib is length-specific. Priority for `-r`: sample
  sheet `bracken_read_length` column > `bracken_read_length_by_platform: {illumina: 150,
  ont: 1000}` > global `bracken.read_length`. The DB must have that `databaseXmers` built.
- **Reproducibility**: `Dockerfile` (core image) + `containers/README.md` (conda-lock for
  bit-exact pins, Apptainer/Singularity). Per-rule conda envs via `metagx run --use-conda`.
- **HPC / schedulers** (`metagx run --executor <name>`): the same workflow submits to a cluster
  via a bundled Snakemake profile. `metagx schedulers` lists the backends — **local** (fat node,
  no scheduler), **slurm** and **lsf** (native Snakemake v8 plugins), **sge** (SGE/UGE/OGS),
  **pbs** (PBS Pro/TORQUE), and **generic** (any other: HTCondor/Moab/OAR/Flux). slurm/lsf use
  native plugins; sge/pbs/generic drive the `cluster-generic` executor with an explicit, **site-
  editable `qsub`/`bsub` submit command**. Each bundled profile under `workflow/profiles/<name>/`
  must be edited once for your site (partition/account/queue, parallel-environment, memory
  resource) — they ship with `CHANGE_ME` placeholders and per-rule thread/memory sizing.
  `--slurm` is a back-compat alias for `--executor slurm`; `--profile <dir>` points at a custom
  external profile. The required plugin (`snakemake-executor-plugin-{slurm,lsf,cluster-generic}`)
  must be in the env. `--executor` is exposed on every surface (CLI, MCP `run_pipeline`,
  HTTP `build-and-run`), and the backends come from one registry (`metagx/schedulers.py`).
- **Host removal** (`host_removal: {genome: <fasta>}`): a first-class pre-classification step
  — maps reads to the host and keeps the unmapped, so the *whole* pipeline runs host-depleted
  (accuracy + PHI). Single-end for now.
- **Managed params** (db/threads/io/`--paired`/`--gzip-compressed`) are injected by the
  workflow — never interview about them.
- **FASTA input** (reads with no quality scores) is auto-detected: QC (fastp) is skipped
  for those samples and `--minimum-base-quality` is dropped from the kraken2 command. No
  user action needed.
- **Subsampling:** add `subsample: {fraction: 0.2, seed: 42}` to classify a reproducible
  random subset (faster/cheaper). Single-end only for now.
- **Cross-platform comparison** (`scripts/compare_platforms.py`): when you have the **same
  biological sample sequenced on multiple platforms** (a mock community, a benchmark, a
  re-sequenced isolate), this puts them side by side. Run metagx per platform (or just point
  at existing kreports/reads/contigs), describe each in a TSV manifest
  (`config/cross_platform.manifest.tsv`: `label, platform, platform_class, kreport, reads,
  contigs, reference`), then `… .venv/bin/python scripts/compare_platforms.py [manifest] [outdir]`.
  It emits an integrated table + figures comparing **classification** (classified %, species
  recovered, diversity, genome-length bias) and **assembly** (contigs, N50, longest, reference
  breadth via minimap2, read-to-contig concordance). Any blank manifest cell skips that block,
  so it works for classify-only or assembly-only sets. Add `--paper` to also emit an IMRaD
  comparison manuscript (`comparison_paper.{tex,pdf}`, pdflatex — same layout as `metagx paper`).
  See experiments 07–08 and
  `BENCHMARKING-DATASETS.md` for valid comparison design (match the abundance model + depth, or
  differences are confounded).
- **Report scope = one project.** A `metagx report`/`paper` covers exactly one project (one
  config = one sample set + the modules enabled on it). All enabled modules combine into that
  single report — adding a module to the *same* data means re-running the *same* config/project,
  and the report regenerates over everything (it grows; you don't get a second paper). Switching
  the **data or the question** is a new project → a new report, and that is correct, not
  fragmentation — stapling unrelated analyses into one paper would misrepresent them. The one
  case where you legitimately combine *across* projects is a **comparison** (same community,
  different platforms/params): that is a deliberate meta-analysis via `metagx compare`
  (`scripts/compare_platforms.py`), not automatic concatenation of per-project reports.

## Input handling (auto-dispatched by `platform` + `layout`)
| Platform | QC | Assembly | minimap2 preset |
|---|---|---|---|
| illumina / mgi (short) | fastp (pe / se / interleaved) | MEGAHIT | `sr` |
| ont | porechop_abi → chopper | Flye `--nano-hq` (`--meta`) | `map-ont` |
| pacbio_hifi | chopper | Flye `--pacbio-hifi` | `map-hifi` |
| pacbio_clr | chopper | Flye `--pacbio-raw` | `map-pb` |
| any, FASTA input | (skipped — no quality) | as above | as above |

Layout: `se` / `pe` / `interleaved` (short only; split + classified as pairs). Long-read
platforms are single-end. Subsampling is single-end only for now.

## Other clients
- **MCP** (Claude Desktop, Cursor, any MCP client): `python mcp_server.py`. Tools:
  `list_pipeline_tools`, `list_presets`, `get_preset`, `list_schedulers`, `get_parameters`,
  `get_interview`, `build_config`, `run_pipeline` (with `executor`), `get_results`,
  `generate_report`, `generate_paper`, `compare_platforms`.
- **Web agents** (ChatGPT Actions, Gemini, Perplexity): `uvicorn mcp_server:app` exposes
  `/api/v1/tools`, `/api/v1/interview/{tool}`, `/api/v1/build-and-run`.
- **No tools at all** (plain chat, Ollama): paste `prompts/INTERVIEW.md`.
