# metagx

A flexible, **interview-driven metagenomics pipeline** that runs the same way whether
it's driven by Claude, Cursor, ChatGPT/Codex, Gemini, Perplexity, Ollama, or a plain
terminal. An LLM interviews the user about what they want, then generates a validated
`config.yaml` for a Snakemake workflow that runs read QC, taxonomic classification,
abundance estimation, and optional assembly + binning.

## The idea: one registry drives everything

kraken2, fastp, and Bracken each expose dozens of flags. Hard-coding them everywhere is
brittle. Instead, each tool has a **parameter registry** in `metagx/parameters/*.yaml`
describing every flag — its type, range, default, whether it's worth asking about, and a
ready-made interview question. That single source of truth drives:

```
metagx/parameters/kraken2.yaml ─┬─> the interview questions   (what the LLM asks)
                                ├─> config validation         (types, ranges, choices)
                                ├─> MCP tool schemas / CLI    (every surface)
                                └─> command-line construction (the Snakefile)
```

Add a flag once → it's available in the interview, the config, and the workflow. No other
file changes.

## Pipeline

`qc` (fastp) → `classify` (kraken2) → `abundance` (Bracken), plus optional `assembly`
(MEGAHIT **or metaSPAdes**) → `binning` (MetaBAT2, with optional `bin_refinement` →
DAS_Tool consensus + dRep). Single-end and paired-end reads are both supported, per
sample. The headline feature is the **confidence sweep**: run kraken2 at several
`--confidence` thresholds and get a per-sample matrix + plot showing how each organism's
read count responds — so you can pick a threshold with evidence.

On top of the core profiling/assembly engine, optional modules add a **functional layer**
(pathways + AMR + MAG annotation), a **second-classifier consensus** cross-check, and a
**run-level aggregate report** (MultiQC + Krona) — see "Functional, consensus & refinement
layers" below.

## Layout

```
metagx/                  Python core (the engine)
  parameters/*.yaml      per-tool parameter registries  <- single source of truth
  presets/*.yaml         named workflow templates (pathogen, gut, soil, quick)
  registry.py            load registries, render args, validate
  config_builder.py      interview answers (+ preset) -> validated config.yaml
  report.py              provenance manifest + Methods paragraph + report
  cli.py                 the `metagx` command
  runner.py              invoke Snakemake
workflow/
  Snakefile              wires modules from config; builds rule all
  rules/*.smk            common, qc, classify, abundance, assembly, binning_refine,
                         consensus, functional, bgc, aggregate, reconcile, domains,
                         stats, differential, amplicon (OTU/ASV)
  envs/*.yaml            isolated conda envs per heavy tool group (--use-conda)
  scripts/*.py, *.R      matrix, combined-abundance, reconcile, consensus, differential, DADA2
mcp_server.py            MCP tools + FastAPI HTTP wrapper for web agents
prompts/INTERVIEW.md     paste-in playbook for tool-less clients
config/config.example.yaml, samples.example.tsv
```

## Read vs contig reconciliation (`modules.reconcile`)

Assembled contigs answer a *different* question than reads: reads give **abundance** (and
the unassembled fraction); contigs give **confident identity** (a long, error-corrected
consensus → deeper, more reliable taxonomy) and recoverable genomes. Enable `reconcile`
(needs `assembly` + `classify`) to classify contigs with kraken2 and join them to
per-contig coverage and the read-level calls. Outputs under `results/<project>/reconcile/`:

- `*.contig_taxonomy.tsv` — per contig: length, coverage, taxon, classified Y/N
- `*.reconciliation.tsv` — per taxon: read % vs **coverage-weighted contig %** (length×depth,
  taxonomy from the consensus) + concordance (`both` / `reads_only` / `contigs_only`)
- `*.flags.tsv` — discordances: reads-only taxa that didn't assemble (low coverage), and
  high-coverage **unclassified** contigs (candidate novel/divergent or chimeric)

Interpretation: `both` = high-confidence core; `reads_only` = low-abundance/fragmented or a
read-level false positive; `contigs_only` = resolved by assembly or a misassembly. Don't run
Bracken on contigs — that's presence/identity, not abundance.

**CAT cross-check (optional).** Build a custom CAT db from your genomes
(`metagx build-cat-db --genomes genomes.fasta --db DIR --taxonomy <kraken2_db>/taxonomy`)
and set `db.cat: DIR/catdb` in the config. Reconcile then classifies contigs a second way —
**CAT**, which votes per predicted ORF (prodigal→DIAMOND→LCA) and so is more robust to
chimeras than whole-contig k-mer LCA — and reports kraken2-vs-CAT agreement per contig
(`methods_agree`) plus a summary count. On the viral test set the two methods agreed on
25/25 classified contigs. **Bin-level taxonomy (BAT / GTDB-Tk)** needs prokaryotic MAGs and
~100 GB references; it's not applicable to viral or unbinnable assemblies.

## Multi-domain taxonomy (`modules.domain_taxonomy`)

Real metagenomes mix viruses, prokaryotes, and eukaryotes — which need *different* tools
after assembly. metagx keeps one broad read-classification pass, then routes the assembly
per domain (set `domains: [viral, prokaryote, eukaryote]`):

| domain | tools | needs |
|---|---|---|
| viral | geNomad (identify + ICTV taxonomy) → CheckV (completeness) | assembly + `db.genomad`, `db.checkv` |
| prokaryote | GTDB-Tk (taxonomy) + CheckM2 (quality) on bins | binning + `db.gtdbtk`, `db.checkm2` |
| eukaryote | EukRep (separate euk contigs) → EukCC (completeness) | assembly + `db.eukcc` |

It also adds **horizontal coverage (breadth)** per contig to `reconcile` — at large reference
scale, breadth is a better "is it really present?" signal than depth. See
`DESIGN-multidomain-and-db-scaling.md` for the full rationale + references.

## Filtered assembly: deplete/target reads, then compare (`modules.filtered_assembly`)

Filter reads by taxon before assembly and measure whether it helped. **Default is
depletion** (`read_filter.mode: exclude`) — remove host/contaminant taxids, keep everything
else *including unclassified* (so novel reads aren't discarded). `mode: include` is targeted
recovery of a clade/species (`include_children`), and you should keep `keep_unclassified:
true` — positively selecting only *classified* reads drops novel/divergent target reads and
can fragment the very assembly you want. Optional `host_genome` adds gold-standard
minimap2-based host depletion before the taxid filter.

Reads are filtered using our own extractor over the persisted per-read kraken2 calls (no
KrakenTools dependency). The unfiltered and filtered assemblies are both produced and
compared in `filtered_assembly/<sample>.assembly_comparison.{tsv,json}` (contig count, total
bp, N50, longest + deltas). Abundance (Bracken) always runs on the **unfiltered** reads.
Caveat: positive selection is reference-biased toward known taxa — prefer depletion for
discovery. (Single-end for now.)

## Library strategy: WGS vs amplicon

Assembly assumes reads are **random genome-wide fragments** — true for WGS shotgun, false for
**amplicon** (16S/18S/ITS), which is many copies of one marker locus. So each sample carries
`library: wgs|amplicon` (default `wgs`), and amplicon samples are routed away from assembly:

| | WGS (shotgun) | amplicon (marker gene) |
|---|---|---|
| QC | fastp / porechop+chopper | **cutadapt** (primer removal) |
| profiling | kraken2/Bracken + assembly/binning/reconcile/domain-taxonomy | **VSEARCH** OTUs or **DADA2** ASVs (short) / **Emu** abundance (long) |
| assembly | ✅ | ❌ skipped (not applicable) |

Set primers under `amplicon: {fwd_primer, rev_primer}`; choose `amplicon.method: otu` (VSEARCH
97% OTUs, default) or `asv` (DADA2 exact sequence variants) for short reads. Emu needs `db.emu` (a 16S DB).
Read-level kraken2/Bracken still runs on amplicon reads (with a warning to prefer a
marker-gene DB + OTU/Emu). An all-amplicon run with assembly modules enabled errors out.
Mixed WGS+amplicon runs are fine — each sample is routed by its `library`. A third value,
`library: ancient`, routes degraded/ancient samples through read-merging + damage
authentication (see "Ancient DNA, decontamination & strain-level" above). See
`DESIGN-multidomain-and-db-scaling.md` for the rationale + references (DADA2/QIIME2, Emu, SILVA).

## Cross-sample statistics & host removal

- **`modules.stats`** (needs `abundance` + ≥2 samples) computes α-diversity
  (Shannon/Simpson/richness/Pielou), β-diversity (Bray–Curtis), **PCoA** ordination, and
  TSS + CLR normalized matrices, with a composition barplot and PCoA plot under
  `results/<project>/stats/`. Pure-Python (numpy) — no external tools.
- **`modules.differential`** (needs `abundance` + a `group` column) answers *which taxa differ
  between conditions?* — CLR-transformed Bracken counts compared between two groups with a
  two-sided **permutation test** + Benjamini-Hochberg FDR (ALDEx2-lite; pure-Python). Writes
  `stats/differential_abundance.{tsv,json}` + a volcano plot. Mark samples `group: case|control`
  in the sheet (≥2 per group to run; more for power — a 2-vs-2 permutation can't reach
  significance). Optional `differential.reference_group` sets the contrast direction.
- **`host_removal: {genome: host.fasta}`** runs a first-class host/contaminant depletion
  (minimap2 → keep unmapped) **before** classification, so the whole pipeline operates on
  host-depleted reads (accuracy + PHI compliance).

## Functional, consensus & refinement layers

These optional modules extend metagx from a "what's here + can I assemble it" engine toward a
study-level platform. Each is a registry + rule + isolated conda env + routing flag, so they
provision their own tools under `--use-conda` and appear in the Methods/citations automatically.

| module | what it does | tools | needs |
|---|---|---|---|
| `functional` | read-based pathways; AMR/virulence on contigs; functional annotation of MAGs | HUMAnN (+MetaPhlAn), AMRFinderPlus + ABRicate, Bakta + eggNOG-mapper | reads (pathways); `assembly` adds AMR; `binning` adds annotation. WGS-only |
| `bgc` | biosynthetic gene clusters (antibiotics, NRPS/PKS, siderophores) on contigs | antiSMASH | `assembly`. WGS-only; DB via `db.antismash` |
| `bin_refinement` | two more binners → consensus bins → cross-sample dereplication | MaxBin2 + CONCOCT → DAS_Tool → dRep | `binning` |
| `classify_consensus` | cross-check kraken2 against an independent classifier (species concordance JSON) | MetaPhlAn (markers) **or** Kaiju (protein) | `classify`. WGS-only |
| `aggregate` | one run-level QC/read-flow report + interactive taxonomy chart | MultiQC, Krona | `classify` |

- **Functional layer** (`modules.functional`) is data-driven: HUMAnN pathway profiling always
  runs on the reads; **AMR/virulence** screening (AMRFinderPlus + ABRicate) runs when
  `assembly` is on (it needs contigs); **Bakta + eggNOG** MAG annotation runs when `binning`
  is on (it needs bins). DBs: `db.humann_nucleotide`/`db.humann_protein`, `db.amrfinderplus`,
  `db.bakta`, `db.eggnog`. Tune e.g. `abricate: {db: card}`.
- **Assembler choice** — set `assembly: {assembler: metaspades}` to use metaSPAdes instead of
  MEGAHIT for short reads (paired-end only). **Hybrid assembly**: give a short-read sample a
  `long_reads` (+ optional `long_platform`) column and metaSPAdes folds the long reads in.
- **Bin refinement** (`modules.bin_refinement`) runs MaxBin2 and CONCOCT alongside MetaBAT2,
  reconciles all three into a consensus per sample with DAS_Tool, then dereplicates the refined
  MAGs across all samples with dRep into a study-level genome catalog.
- **Consensus classifier** (`modules.classify_consensus`) runs MetaPhlAn or Kaiju (pick with
  `consensus: {classifier: kaiju}`) and reports species-level concordance with kraken2 —
  agreement is a confidence signal; kraken2-only taxa flag DB-completeness false positives.
- **Aggregate** (`modules.aggregate`) writes `report/multiqc/multiqc_report.html` and
  `report/krona.html`.

## Ancient DNA, decontamination & strain-level

| module / field | what it does | tools | needs |
|---|---|---|---|
| `library: ancient` + `modules.damage` | collapse overlapping fragments, then authenticate post-mortem damage | fastp `--merge`, mapDamage2 | `assembly` + an ancient sample |
| `modules.decontam` + `control: true` | flag & remove reagent/lab contaminants from negative controls | pure-Python (prevalence) | `abundance` + a control sample |
| `modules.strain` | within-population SNV / microdiversity profiling | inStrain | `assembly` |

- **Ancient DNA** — mark degraded samples with `library: ancient` (short-read, paired-end).
  fastp collapses the overlapping pairs into single fragments; the reads are classified and
  assembled, then `modules.damage` maps them to the assembly and runs mapDamage2 to measure
  the C→T (5') / G→A (3') deamination signature. `ancient/<sample>/authentication.json` gives a
  verdict (`damage_present`) — the test that separates authentic ancient DNA from modern
  contamination. Preset: `ancient-dna`.
- **Decontamination** — mark negative/blank samples with `control: true`. `modules.decontam`
  applies the decontam *prevalence* test over the combined Bracken table: taxa that are as (or
  more) prevalent in the controls than in real samples are flagged as contaminants and removed,
  producing `stats/decontam_flagged.tsv` + `stats/abundance_decontaminated.tsv`. Essential for
  low-biomass samples. Pure-Python — no external tool.
- **Strain level** — `modules.strain` runs inStrain on the reads-vs-contigs mapping to profile
  SNVs and nucleotide diversity, resolving mixed strains that abundance profiling cannot.
- **Per-sample Bracken length** — Bracken's k-mer distribution is length-specific, so mixed
  short+long runs need the right `-r` per sample. Priority: a sample-sheet `bracken_read_length`
  column > `bracken_read_length_by_platform: {illumina: 150, ont: 1000}` > the global
  `bracken.read_length`. The DB must have that `databaseXmers` built (`metagx build-db` accepts
  a comma list of lengths).
- **Reproducibility** — `Dockerfile` builds a core-tool image; `containers/README.md` documents
  conda-lock (bit-exact pins) and Apptainer/Singularity conversion. See the three levels there.

## Tests

`pip install -e ".[test]" && pytest` — unit tests for the registry, config builder,
diversity, formats, subsampling, and read-filter logic. CI runs them on every push
(`.github/workflows/ci.yml`). See `CRITIQUE.md` for the full gap analysis + roadmap.

## Self-provisioning tools (conda)

Core tools come from `environment.yml` (`conda env create -f environment.yml`). The heavy
domain tools live in isolated envs under `workflow/envs/`; run `metagx run --use-conda` and
Snakemake creates them per-rule on first use — so the pipeline installs what it needs.

## Presets, provenance & reports (inspired by K-Dense BYOK)

- **Presets** — `metagx presets` lists ready-to-run templates (`pathogen-detection`,
  `gut-profiling`, `soil-deep-assembly`, `quick-screen`, `amr-surveillance`, `ancient-dna`).
  Pick one to pre-fill modules and parameters, then override only what you care about. They're
  built on the same registries, so they can't drift from the real flags.
- **Provenance + "Copy as Methods"** — after a run, `metagx report` captures a
  `manifest.json` (tool versions, exact command lines, database identity/size, QC read
  counts, % classified) and writes a `methods.md` paragraph with citations, ready to paste
  into a paper.
- **Report** — `metagx report --format md|latex|pdf` bundles methods + a registry-annotated
  parameter table + figures + top-abundance table into one document (LaTeX/PDF need pandoc).
- **Full IMRaD paper** — `metagx paper` elaborates the run into a complete manuscript
  (Introduction / Methods / Results / Discussion + abstract, tables, figures, references) and
  compiles it to **PDF with pdflatex** → `results/<project>/report/paper.{tex,pdf}`. The design
  comes from the interview, methods + citations from the registries, and every result number is
  read back from the output files — a publishable first draft to refine, not a fabricated one.
  `--no-pdf` writes only the `.tex`. (MCP: `generate_paper`.)

## Quick start

```bash
bash setup.sh                 # uv venv + deps + dirs
uv pip install -e .           # install the `metagx` CLI and package

metagx tools                  # list pipeline steps
metagx presets                # list workflow presets (pick a starting point)
metagx interview kraken2      # questions an LLM should ask (JSON)
metagx build-db --genomes data/genomes.fasta --db local_databases/custom   # custom DB
metagx build-config answers.json   # validate answers (+preset) -> config.yaml
metagx run --config config.yaml --dry-run
metagx run --config config.yaml
metagx results --config config.yaml
metagx report --config config.yaml # provenance manifest + Methods + report
```

Don't have a database yet? `snakemake --snakefile workflow/Snakefile fetch_kraken_db`
downloads a standard index (large — only run when you mean it).

## Driving it from each client

| Client | How |
| --- | --- |
| Claude Code / Claude.ai | the `SKILL.md` skill |
| Claude Desktop, Cursor (MCP) | `python mcp_server.py` |
| ChatGPT Actions, Gemini, Perplexity | `uvicorn mcp_server:app --port 8000` (HTTP) |
| Codex, terminal agents | the `metagx` CLI |
| Ollama / plain chat | paste `prompts/INTERVIEW.md` |

## Flexible inputs: FASTA, FASTQ, custom DBs, subsampling

- **Custom databases** — `metagx build-db --genomes genomes.fasta --db DIR` builds a
  kraken2 **and** Bracken database from your own reference genomes. It assigns each genome
  a synthetic taxid and a minimal taxonomy (no multi-GB NCBI download), tags sequences with
  `kraken:taxid|` headers, and runs `kraken2-build` + `bracken-build`. Degrades gracefully:
  if `bracken-build` is absent, the kraken2 db is still built.
- **FASTA reads (no quality scores)** — handled automatically. Format is detected per
  sample; for FASTA the QC step (fastp, a FASTQ-only tool) is **skipped**, the FASTQ-only
  `--minimum-base-quality` flag is **dropped** from the kraken2 command, and
  `--gzip-compressed` is only added for `.gz` inputs.
- **Layout & platform** — each sample carries `layout` (`se` / `pe` / `interleaved`) and
  `platform`. The platform dispatches the right tools, since adapter/QC chemistry and
  assemblers differ by technology:

  | platform | QC | assembly | minimap2 |
  |---|---|---|---|
  | `illumina` / `mgi` (short) | fastp | MEGAHIT | `sr` |
  | `ont` | porechop_abi → chopper | Flye `--nano-hq` | `map-ont` |
  | `pacbio_hifi` | chopper | Flye `--pacbio-hifi` | `map-hifi` |
  | `pacbio_clr` | chopper | Flye `--pacbio-raw` | `map-pb` |

  Interleaved short reads are split by fastp and classified as pairs (kraken2 `--paired`).
  Long-read platforms are single-end. All driven from the registries — no rule edits to
  retune a flag.
- **Subsampling** — add `subsample: {fraction: 0.2, seed: 42}` to classify a random,
  reproducible 20% of reads (faster/cheaper). Format-aware and dependency-free
  (single-end for now). Subsampling feeds QC/classification automatically.

### Worked example (bundled `data/`)

```bash
# 30 viral reference genomes -> custom kraken2 + Bracken DB
metagx build-db --genomes data/genomes.fasta --db local_databases/viral_custom --read-length 150

# Classify 20% of 20,000 FASTA reads across a confidence sweep, with abundance + report
metagx build-config config/viral-test.answers.json --out config.yaml
metagx run    --config config.yaml --cores 4
metagx report --config config.yaml
```
This produces the confidence sensitivity curve (here: ~49% → ~47% → ~0.6% classified at
confidence 0.0 / 0.1 / 0.5), Bracken abundances, and a paste-ready Methods paragraph.

## Requirements

Python ≥3.10 and these bioinformatics tools on `PATH` for the steps you enable:
`kraken2`, `bracken`, `fastp`, and (for assembly/binning) `megahit`, `minimap2`,
`samtools`, `metabat2`. Install via conda/mamba, e.g.
`mamba install -c bioconda kraken2 bracken fastp megahit minimap2 samtools metabat2`.

The heavier optional layers — `metaspades`; `humann`/`amrfinderplus`/`abricate`/`bakta`/
`eggnog-mapper`; `maxbin2`/`concoct`/`das_tool`/`drep`; `metaphlan`/`kaiju`; `multiqc`/`krona` —
each have an isolated env under `workflow/envs/` and are provisioned automatically by
`metagx run --use-conda` (Snakemake creates each env on first use), so you only install what
the modules you enable actually need. They also need their own reference DBs (see
`db.*` in `config/config.example.yaml`).

**Apple Silicon (arm64) macOS:** bioconda has no native arm64 builds for kraken2/bracken,
and its osx-64 `bracken` package is a broken placeholder. Run
`bash scripts/install_bio_macos_arm64.sh`, which creates an osx-64 (Rosetta) conda env for
kraken2/fastp/seqtk and builds Bracken from source (with Homebrew `libomp`), including a
fix for an upstream `bracken` wrapper bug. Then put the env on `PATH` before running metagx.
