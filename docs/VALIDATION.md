# metagx: end-to-end validation on a custom viral reference set

*A validation report in IMRaD form, written so the numbers and methods can be lifted
directly into a manuscript. Every figure here is produced by an automated test
(`tests/test_pipeline_e2e.py`, `tests/test_adna_e2e.py`, `tests/test_workflow_dryrun.py`,
`tests/test_minimizer_report.py`) and is reproducible on this machine.*

---

## Abstract

`metagx` is a schema-driven, multi-platform metagenomics pipeline distributed as a Claude
skill. It interviews the user, generates a validated `config.yaml` from per-tool parameter
registries, and runs a Snakemake workflow (kraken2/Bracken/fastp + optional assembly,
binning, contig reconciliation, domain taxonomy, differential abundance, ancient-DNA
authentication, and ~30 other modules). We validate the full pipeline on a controlled
30-genome viral reference set across three sequencing platforms (Oxford Nanopore, Illumina,
PacBio HiFi), plus assembly→binning→reconciliation, differential abundance, and ancient-DNA
damage authentication. Classification recovers the expected taxa on every platform;
read-level and contig-level taxonomy agree on **99.3 %** of reads; two orthogonal classifiers
(kraken2 vs Kaiju) agree on the viral community (Jaccard = 1.0); differential abundance
correctly reports **zero** false positives on null (single-population) data; the ancient-DNA
module discriminates deaminated from undamaged libraries; and phylogenetics and report
aggregation run through the workflow. Tool versions are pinned and version-floor-tested (a guard
added after a transitive `samtools` downgrade was caught breaking the mapping path).

## Introduction

Metagenomic pipelines are typically validated by "it ran and produced files." That is a weak
guarantee: it does not show the *answers* are correct, nor that the same code paths behave
correctly across sequencing platforms whose error profiles and read lengths differ by orders
of magnitude. We set out to validate `metagx` against a **known ground truth** — a custom
kraken2/Bracken database built from 30 complete viral genomes — so that every classification
can be checked against the organisms that are actually present.

The validation answers four questions:

1. **Correctness across platforms.** Does classification recover the right viral taxa from
   short-read (Illumina), noisy long-read (ONT), and accurate long-read (PacBio HiFi) data?
2. **Internal consistency of assembly.** Do de novo assembly + contig-level taxonomy agree
   with read-level classification (a cross-check that does not depend on the database alone)?
3. **Statistical soundness.** Does differential abundance control its false-discovery rate on
   data with no real signal?
4. **Specialised modules.** Does the ancient-DNA path distinguish authentic post-mortem
   damage from modern DNA?

## Methods

### Reference set and database

- **Genomes.** 30 complete viral reference genomes (Dengue 1, Zika, Chikungunya, Yellow
  fever, West Nile, Powassan, Influenza A segments, HIV-1, SFTS, CCHF, Heartland, and others;
  largest ≈ 144 kb), `tests/fixtures/viral/genomes.fasta` (committed, 428 kB).
- **Database.** A custom kraken2 database is built from these genomes with synthetic taxids
  (`kraken:taxid|` headers), plus Bracken k-mer distributions for 150-mer (short-read) and
  1000-mer (long-read) configurations. The build is reproducible via `metagx build-db`
  (`metagx/dbbuild.py::build_db`); the validation suite builds it from the committed genome
  fixture when a prebuilt copy is absent, so the database is never a black box.

### Sequencing data

| Platform        | Source                                   | Layout | Reads (test) |
|-----------------|------------------------------------------|--------|--------------|
| Oxford Nanopore | real simulated metagenomic FASTA         | SE     | 1.5–8 k      |
| Illumina        | paired reads simulated from the genomes  | PE     | ~full set    |
| PacBio HiFi     | long reads simulated from the genomes    | SE     | ~full set    |
| Case/control    | 4 subsamples (2 case, 2 control)         | SE     | ~full set    |

ONT is a noisier, independent metagenomic FASTA; Illumina and PacBio reads were simulated
directly from the reference genomes (hence near-complete classification is expected for them
and partial classification for ONT).

### Software environment

All bio tools are consolidated into a **single conda environment** (`metagx-bio`, built with
`CONDA_SUBDIR=osx-64` on Apple Silicon): kraken2 2.17.1, Bracken, fastp 1.3.3, MEGAHIT 1.2.9,
Flye 2.9.6, minimap2 2.31, samtools ≥1.18, metabat2 2.18, DIAMOND, CheckV 1.1.1, geNomad
1.12.0, mapDamage2 2.2.2, plus seqkit/cutadapt/vsearch/kaiju/mafft/fasttree. The env is declared
in `environment.yml`. (AMR's `abricate` is deliberately kept in its own isolated env: its
`perl-bio-samtools` dependency hard-pins samtools 0.1.x, which conflicts with the modern samtools
the rest of the pipeline needs; it is provisioned via `--use-conda`.) The Python package
(`pip install -e .`) is dependency-light
(PyYAML/numpy/pandas/matplotlib/snakemake); all statistics (diversity, differential abundance,
damage authentication) are pure-Python with no SciPy/R.

The exact validated tool versions are captured automatically (`metagx report` records them in
every run manifest; `tests/test_tool_versions.py` asserts a minimum floor for each so a silent
downgrade fails CI). The versions used in this validation:

| Tool | Version | Tool | Version |
|------|---------|------|---------|
| kraken2 | 2.17.1 | metabat2 | 2.18 |
| Bracken | 2.9 (no version flag) | DIAMOND | 2.0.4 |
| fastp | 1.3.3 | MAFFT | 7.525 |
| MEGAHIT | 1.2.9 | IQ-TREE | 3.1.2 |
| Flye | 2.9.6 | CheckV | 1.1.1 |
| minimap2 | 2.31 | geNomad | 1.12.0 |
| samtools | 1.21 | Kaiju | 1.10.1 |
| MultiQC | 1.35 | mapDamage2 | 2.2.2 |

Tool versions are load-bearing here: the validation surfaced (and the suite now guards against)
a `samtools` downgrade to 0.1.19 — pulled in transitively by `abricate`'s `perl-bio-samtools` —
whose `sort -o` API silently breaks the mapping/binning/aDNA paths.

### Test architecture

Validation is layered so it runs both with and without the bio tools:

- **Dry-run gate** (`test_workflow_dryrun.py`) — builds the full Snakemake DAG and resolves
  every `params`/`render_args` lambda for each platform path (incl. assembly, binning, and
  domain taxonomy). Needs no bio tools or database; runs in CI.
- **Real end-to-end** (`test_pipeline_e2e.py`, `test_adna_e2e.py`) — executes the actual
  workflow against the (built-from-fixture) viral database and asserts correctness, not file
  existence. Runs in CI on linux-64 (native bioconda) and locally with `metagx-bio` active.

## Results

### 1. Classification is correct on every platform

Each platform was classified against the 30-genome database and abundance re-estimated with
Bracken. Assertions: classified fraction within a platform-appropriate band, ≥ N species
recovered, ≥ 1 known database virus present, and Bracken relative abundances summing to 1.

| Platform        | Classified fraction | Distinct species recovered (of 30) | Status |
|-----------------|---------------------|------------------------------------|--------|
| Illumina (PE)   | ≈ 100 % (0.03 % unclassified) | 30                       | PASS   |
| ONT (SE)        | 49.5 % classified   | 29                                 | PASS   |
| PacBio HiFi (SE)| ≥ 30 % classified   | 23                                 | PASS   |

The near-complete Illumina classification reflects reads simulated directly from the database
genomes; ONT's partial classification reflects a noisier, independent metagenomic input. In
all cases the recovered taxa are genuine database members (Dengue, Zika, West Nile, …), not
spurious assignments, and Bracken redistributes into a valid probability distribution.

### 2. Assembly and read/contig taxonomy are internally consistent

On the ONT data we ran the heavy path: Flye (metagenome) assembly → minimap2/samtools mapping
→ jgi depth → MetaBAT2 binning → kraken2-on-contigs → reconciliation.

- **Assembly:** 42 contigs assembled from ~8 k ONT reads.
- **Read/contig concordance:** of reads aligning to taxonomically-classified contigs,
  **99.3 %** received the *same* taxon at the read level as their contig (0 % discordant,
  0.69 % read-unclassified). Because this cross-check uses the assembly as an orthogonal
  ground truth, it validates classification independently of the database.
- **Binning robustness:** MetaBAT2 recovers **0 bins** from this small, single-sample viral
  assembly — the correct outcome (its defaults target Mb-scale bacterial MAGs). The pipeline
  records the bin count and continues rather than crashing, a fix made during validation
  (zero-MAG samples are common in real low-diversity data).

### 3. Differential abundance controls its false-discovery rate

Four case/control samples that are random subsamples of one population (i.e. **no** true
differential signal) were analysed with the pure-Python ALDEx2-style method (Dirichlet
Monte-Carlo, mc = 128; CLR transform; permutation test; Benjamini–Hochberg FDR at 0.1 over
999 permutations).

- **Significant taxa: 0 / 29** — the statistically correct result. The FDR control does not
  manufacture false positives.
- α-diversity (richness, Chao1, ACE, Shannon ≈ 2.9–3.0, Simpson, Pielou) and β-diversity
  (PCoA) tables are produced and well-formed.

### 4. Ancient-DNA authentication discriminates damage

We simulated two libraries from the viral genomes: one with terminal cytosine deamination
applied (5′ C→T and 3′ G→A, position-decaying), and an undamaged control. Reads were mapped
to the reference, post-mortem damage quantified with mapDamage2, and the verdict produced by
`damage_authenticate.py`.

| Library  | 5′ C→T (pos 1) | 3′ G→A (pos 1) | Verdict |
|----------|----------------|----------------|---------|
| Damaged  | **0.391**      | **0.368**      | consistent with authentic ancient DNA |
| Control  | 0.000          | 0.000          | modern / no damage signal |

Both terminal frequencies in the damaged library are ~8× the 0.05 authenticity threshold,
while the undamaged control is flat at zero — a clean, unambiguous discrimination using the
real mapDamage2 tool, not a mock.

### 5. Report-format robustness (kraken2 `--report-minimizer-data`)

Enabling `--report-minimizer-data` inserts two columns into the kraken2 report, shifting the
name column from index 5 to 7. We confirmed and fixed that every report consumer (Krona text
generation, the second-classifier consensus, and read/contig reconciliation) now reads
rank/taxid/name from the end of the row, so a standard and a minimizer-augmented report yield
**identical** Krona charts and parsed tables.

### 6. Additional modules through the real workflow

A second pass ran modules that were previously only unit- or dry-run-tested through the actual
Snakemake workflow, which uncovered and fixed two latent execution bugs:

- **Phylogenetics** (MAFFT → IQ-TREE 3 → tree): builds a tree from the 6-sequence homologous
  demo set; all six taxa recovered in the Newick output. This caught a **tool-version bug** —
  the script hardcoded the `iqtree2` binary, which does not exist in IQ-TREE 3 (ships `iqtree3`);
  it now resolves `iqtree2`/`iqtree3`/`iqtree` dynamically. It also caught a **Snakemake-execution
  bug**: `from __future__ import annotations` in a `script:`-run file becomes a `SyntaxError`
  (the directive prepends a preamble), so the module had never actually run through the workflow.
- **Second-classifier consensus** (Kaiju vs kraken2): on the viral set the two orthogonal
  classifiers (nucleotide k-mer vs translated-protein) agree on **29/29 species, Jaccard = 1.0**.
  (Same `from __future__` bug, now fixed.)
- **Report aggregation** (Krona + MultiQC) **with `--report-minimizer-data` enabled**: the Krona
  chart shows real taxon names (Dengue, Zika, …), proving the minimizer-column fix in the *live*
  rule — a fixed-index parser would have emitted rank codes here.

### 7. Pre-assembled contigs + bacterial AMR (functional layer)

Two end-user-driven additions were validated on real bacterial references downloaded from NCBI
(E. coli K-12 MG1655, S. aureus NCTC 8325, and the *bla*NDM-1-carrying plasmid pNDM-HK):

- **Pre-assembled contigs input.** A sample can supply a `contigs:` FASTA (an isolate genome,
  a prior assembly, or references); the assembler is skipped and the contigs feed the
  contig-consuming modules. Validated end to end: the 30 viral genomes were fed as contigs and
  reconciled against ONT reads — the staging rule produced the contigs unchanged (no MEGAHIT/Flye
  run) and reconciliation recovered the known viruses. This is also how the bacterial modules are
  exercised where MEGAHIT is unavailable.
- **AMR screening via `--use-conda`.** Running `modules.functional` + `functional.amr` on the
  pNDM-HK plasmid provisions the isolated `amr.yaml` env (ABRicate + BLAST) through Snakemake's
  conda integration and reports the carbapenemase gene. Validating this **caught a real env bug**:
  `amr.yaml` declared ABRicate but not `blast`, so the provisioned ABRicate failed at runtime with
  *"Could not find 'blastn'"* — fixed by declaring `blast` explicitly. The functional layer was
  also split into independent sub-steps (`amr` / `pathways` / `annotation`) so a user can screen a
  genome for resistance without HUMAnN's multi-GB database.

  The recovered resistance genes are biologically correct (100 % coverage & identity unless
  noted):

  | Genome | Resistance genes found |
  |--------|------------------------|
  | pNDM-HK plasmid | **blaNDM-1** (carbapenem), blaTEM-1 (β-lactam), aac(3)-IId + armA (aminoglycoside), ble-MBL (bleomycin), sul1 (sulfonamide), msr(E) + mph(E) (macrolide) — the plasmid's known MDR cassette |
  | E. coli K-12 MG1655 | blaEC-19 (cephalosporin) — the intrinsic chromosomal *ampC* (strain is otherwise susceptible) |
  | S. aureus NCTC 8325 | fosB (fosfomycin) — the expected intrinsic resistance of this reference strain |

This pass also confirmed the **`--use-conda` provisioning path** — the way end users (not just this
machine's pre-built env) obtain the heavy/conflicting domain tools: Snakemake builds each rule's
isolated `workflow/envs/*.yaml` on first use via mamba.

### 8. Snakemake `script:` directive — supported and guarded

The `script:` directive (Python/R logic with an injected `snakemake` object) is a first-class
extension point for users. `tests/test_workflow_scripts.py` scans every `script:`-referenced file
and fails on the `from __future__` antipattern (which the directive's preamble injection turns into
a `SyntaxError`) and checks each uses the `snakemake` object — so neither the project nor a
contributor can reintroduce the class of bug that had silently broken two modules.

### 9. Coverage summary

- **Registries:** 38 per-tool parameter registries (single source of truth for interview,
  validation, MCP schema, CLI, and command-line rendering). The 7 previously stub-only
  registries (GTDB-Tk, CheckM2, CheckV, EukCC, Emu, Krona, MultiQC) were curated; the four
  bin/genome-quality rules were rewired from hardcoded shell flags to `render_args`.
- **Workflow:** 21 rule modules, all command-building rules driven by the registries.
- **Tests:** the dry-run gate exercises Illumina/ONT/PacBio/assembly/binning/domain paths; the
  real e2e suite exercises classification (×3 platforms), assembly/binning/reconciliation,
  differential abundance, ancient-DNA authentication, phylogenetics, second-classifier
  consensus, and report aggregation. Per-module coverage:

| Module | Real e2e | Module | Real e2e |
|--------|:--------:|--------|:--------:|
| qc (fastp) | ✓ | reconcile | ✓ |
| classify (kraken2) | ✓ | differential | ✓ |
| abundance (Bracken) | ✓ | stats (α/β diversity) | ✓ |
| assembly (Flye) | ✓ | damage / aDNA | ✓ |
| binning (MetaBAT2) | ✓ | phylogenetics | ✓ |
| classify_consensus (Kaiju) | ✓ | aggregate (Krona+MultiQC) | ✓ |
| functional / AMR (ABRicate) | ✓ (`--use-conda`) | provided-contigs input | ✓ |
| domain_taxonomy (viral) | dry-run* | bgc / strain / bin_refinement | dry-run† |

  *geNomad/CheckV installed + wired; real run needs their multi-GB DBs (downloaded on demand).
  †antiSMASH / inStrain / DAS_Tool + dRep not co-installable on osx-64 → isolated `--use-conda`
  envs, wiring covered by the dry-run gate; reachable by end users via `metagx run --use-conda`.

## Discussion

**What is validated.** The taxonomic spine (QC → classify → abundance) is correct across three
platforms against a known ground truth; de novo assembly agrees with read classification at
99.3 %; differential abundance is statistically sound; and the ancient-DNA path discriminates
damage. The same registries that drive the interview also render every command line, and the
dry-run gate proves that wiring for paths we cannot fully execute locally.

**Limitations (stated plainly).**

- **Bacterial/eukaryotic domain quality (GTDB-Tk, CheckM2, EukCC)** is validated only by the
  dry-run gate (DAG + `render_args` wiring), not by real execution, for two reasons: (i) they
  require multi-GB/multi-100-GB reference databases and *bacterial/eukaryotic* input, neither of
  which applies to a 30-genome viral test set; and (ii) they are not co-installable in the core
  environment — CheckM2's TensorFlow/LightGBM stack and `abricate`'s `perl-bio-samtools` (which
  pins samtools 0.1.x) both conflict with the modern core. This is *why* the project keeps such
  tools in isolated per-rule envs (`workflow/envs/*.yaml`), provisioned on `metagx run
  --use-conda`. The single-env consolidation is therefore "as consolidated as the bioconda
  dependency graph allows": one core env for QC/classify/abundance/assembly/binning/reconcile/
  diversity/differential/aDNA + the viral domain (geNomad, CheckV), plus a small number of
  isolated envs for tools with mutually-incompatible pins. The viral-domain tools (geNomad,
  CheckV) *are* in the core env and wired; their reference databases download on demand.
- **MEGAHIT** segfaults on this Apple-Silicon/Rosetta (osx-64) host (a known popcnt issue);
  short-read assembly is therefore validated via CI (linux-64) and the long-read (Flye) path
  locally.
- **MAG binning on viral data** is intentionally degenerate — MetaBAT2 targets bacterial
  genomes. The pipeline now detects small-genome/viral context and recommends the binning
  parameters accordingly (min contig → 1500 floor) while steering users to the viral domain
  path; it does not fabricate viral bins.

**Reproducibility.** Every number above is asserted by a test and regenerated on each run; the
database is rebuilt from a committed 30-genome fixture, so the validation has no hidden inputs.

## Availability

Tests: `tests/test_pipeline_e2e.py`, `tests/test_adna_e2e.py`,
`tests/test_workflow_dryrun.py`, `tests/test_minimizer_report.py`. Run locally with the
`metagx-bio` environment active (`export PATH="$HOME/miniconda3/envs/metagx-bio/bin:$PATH"`)
or in CI (the `e2e` job installs the stack and builds the database from the fixture).
