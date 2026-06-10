# Test datasets for metagx (by sequencing platform)

Public metagenomic datasets to exercise each pipeline direction. **Accessions are
best-effort — verify before bulk download.** Fetch from SRA/ENA with
`prefetch` + `fasterq-dump` (SRA) or the ENA "Download files" links.

> **Tip:** the single most useful resource is a **mock community sequenced on multiple
> platforms** — the *same sample* with known composition lets you test every direction and
> actually measure accuracy (precision/recall, breadth, read-vs-contig concordance).

## Cross-platform gold standard (recommended starting point)

| Resource | Platforms | Where |
|---|---|---|
| **ZymoBIOMICS Microbial Community Standards** (D6300 even, D6310/D6311 log, D6331 gut) | Illumina, ONT, PacBio | product info: zymoresearch.com; exact read accessions below |
| **Loman Lab mock community** (Nicholls et al. 2019) | ONT (GridION, PromethION) + Illumina | ENA **PRJEB29504**; github.com/LomanLab/mockcommunity — *GigaScience* 8:giz043, doi:10.1093/gigascience/giz043 |
| **CAMI / CAMI II** (simulated + truth, all-domain) | Illumina-like (configurable) | data.cami-challenge.org — Meyer et al. 2022 *Nat Methods* 19:429, doi:10.1038/s41592-022-01431-4 |

### ZymoBIOMICS — exact download links per platform

The ZymoBIOMICS standards are the de-facto multi-platform mock community. **Caveat for a
*matched* cross-platform comparison:** the public ONT + Illumina reads are of the **D6300
(even) / D6310 (log)** standards (Loman Lab), while the public PacBio HiFi reads are of the
**D6331 gut** standard — *different communities*. The cleanest *same-sample* pair is the Loman
ONT + Illumina (identical physical standards); treat the PacBio HiFi set as a separate,
HiFi-specific benchmark unless you sequence one standard yourself on all three. See
`BENCHMARKING-DATASETS.md` for why matching the community matters.

**Products (Zymo Research):**
- Microbial Community Standard D6300 (even): https://www.zymoresearch.com/products/zymobiomics-microbial-community-standard
- Microbial Community DNA Standard D6305 (even, DNA): https://www.zymoresearch.com/products/zymobiomics-microbial-community-dna-standard
- Gut Microbiome Standard D6331: https://www.zymoresearch.com/products/zymobiomics-gut-microbiome-standard

**Illumina + Oxford Nanopore — Loman Lab (D6300 even / D6310 log), ENA study `PRJEB29504`**
Portal (download index + S3 links): https://lomanlab.github.io/mockcommunity/

| Platform | Even (D6300) | Log (D6310) | ENA browse |
|---|---|---|---|
| ONT GridION (FASTQ) | `ERR3152364` | `ERR3152366` | https://www.ebi.ac.uk/ena/browser/view/PRJEB29504 |
| ONT PromethION (FASTQ) | `ERR3152365` | `ERR3152367` | (same study) |
| Illumina (isolates, for truth refs) | `ERR2935848`–`ERR2935857` | — | one accession per isolate |

Direct ENA record, e.g. GridION even FASTQ: https://www.ebi.ac.uk/ena/browser/view/ERR3152364 ·
example S3 asset: https://nanopore.s3.climb.ac.uk/mockcommunity/v2/Zymo-Isolates-SPAdes-Illumina.fasta

**PacBio HiFi — ZymoBIOMICS D6331 gut, NCBI project `PRJNA680590`** (PacBio data catalog:
https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/master/docs/PacBio-Data.md)

| Library | SRA run | Reads / yield |
|---|---|---|
| Standard | `SRR13128014` | 1.98 M reads, 9.1 kb mean, 17.9 Gb (Q39) |
| Low | `SRR13128013` | 2.77 M reads, 9.3 kb mean, 25.8 Gb |
| Ultra-Low | `SRR13128012` | 2.48 M reads, 8.6 kb mean, 21.3 Gb |

Browse: https://www.ncbi.nlm.nih.gov/sra/?term=SRR13128014 — fetch with
`prefetch SRR13128014 && fasterq-dump SRR13128014`. Other PacBio HiFi mocks in the same catalog:
ATCC MSA-1003 `SRR9328980` (PRJNA546278); 71-species mock `ERR9765783` (PRJEB52977).

## Illumina (short read)

| Dataset | Notes | Link |
|---|---|---|
| Human Microbiome Project (HMP1/HMP2) WGS | Large real human metagenomes | hmpdacc.org → SRA |
| Zymo mock (Illumina) | Ground-truth composition | within PRJEB29504 / Zymo studies |
| CAMI II toy/marine/strain-madness | Gold-standard assemblies + profiles | data.cami-challenge.org |

## Oxford Nanopore (ONT, long read)

| Dataset | Notes | Link |
|---|---|---|
| Loman mock (GridION/PromethION) | Even + log, ground truth — ideal for our ONT path | ENA **PRJEB29504** |
| Long-read profiling benchmark (Portik et al. 2022) | Mock ONT + HiFi, with profiling truth | *BMC Bioinformatics* 23:541, doi:10.1186/s12859-022-05103-0 (data links in paper) |

## PacBio (HiFi / CLR, long read)

| Dataset | Notes | Link |
|---|---|---|
| Zymo D6331 gut standard (HiFi) | Mock, HiFi | PacBio DevNet / SRA — search "Zymo D6331 HiFi" |
| Bickhart et al. 2022 (sheep gut, HiFi) | Real complex community → complete MAGs | *Nat Biotechnol* 40:711, doi:10.1038/s41587-021-01130-z (SRA in paper) |
| Portik et al. 2022 (HiFi mock) | Profiling benchmark with truth | doi:10.1186/s12859-022-05103-0 |

## Ion Torrent

Shotgun metagenomics on Ion Torrent is uncommon (it's mostly 16S amplicon), so curated
public *shotgun* sets are sparse. Best approach — filter SRA directly:

- NCBI SRA Advanced Search: `Platform = ION_TORRENT` AND `Strategy = WGS` AND
  `Source = METAGENOMIC`.
- For amplicon/16S Ion Torrent test data: `Platform = ION_TORRENT` AND `Strategy = AMPLICON`.

## Bundled Illumina simulation (wgsim — test-only)

> **Removed 2026-06-10.** `wgsim` (and the PacBio simulators below) were uninstalled from
> `metagx-bio` to reclaim space (env 1.8 GB → 1.6 GB; `numpy`/`scipy`/`pysam` went with
> them as orphaned deps — the pipeline does not use those from conda, it uses `.venv`).
> The simulated FASTQ under `data/illumina_sim/` and `data/pacbio_sim/` is **kept** and the
> tests still run. To regenerate from scratch, reinstall first:
> `conda install -n metagx-bio -c bioconda wgsim pbsim3 pbccs simlord`.

When you need **real FASTQ** (not FASTA) to exercise fastp QC and paired-end kraken2 on the
30-virus reference set without downloading public data:

```bash
export PATH="$HOME/miniconda3/envs/metagx-bio/bin:$PATH"
bash scripts/simulate_illumina_wgsim.sh
metagx run --config config/illumina-sim.yaml --cores 4
metagx run --config config/illumina-sim-diff.yaml --cores 4   # + subsample + differential
```

`wgsim` is a **test-only** dependency — remove it on production machines; see
`scripts/REMOVE-wgsim-test-only.md`.

## Bundled PacBio simulation (PBSIM3 / SimLoRD / CCS — test-only)

Synthetic PacBio CLR and HiFi FASTQ from `data/genomes.fasta` for chopper QC and
long-read kraken2/Bracken (requires `database1000mers.kmer_distrib` on the viral DB):

```bash
export PATH="$HOME/miniconda3/envs/metagx-bio/bin:$PATH"
# one-time: bracken-build -d local_databases/viral_custom -t 4 -k 35 -l 1000
bash scripts/simulate_pacbio.sh
metagx run --config config/pacbio-hifi-sim.yaml --cores 4
metagx run --config config/pacbio-clr-sim.yaml --cores 4
```

- **CLR**: PBSIM3 single-pass (`pacbio_clr.fastq.gz`) or SimLoRD (`pacbio_clr_simlord.fastq.gz`)
- **HiFi**: SimLoRD multipass on macOS; PBSIM3 multipass + `ccs` on Linux (bioconda `pbccs`)

Remove test tools with `scripts/REMOVE-pacbio-sim-test-only.md`.

## Generic SRA/ENA recipe (any platform)

SRA Advanced Search → set **Platform** (`ILLUMINA` / `OXFORD_NANOPORE` / `PACBIO_SMRT` /
`ION_TORRENT`) + **Strategy** `WGS` + **Source** `METAGENOMIC`, then:

```bash
prefetch SRRXXXXXXX
fasterq-dump --split-files SRRXXXXXXX        # paired Illumina
fasterq-dump SRRXXXXXXX                       # single-end / long reads
```

In the metagx sample sheet, set `platform` accordingly (`illumina|ont|pacbio_hifi|
pacbio_clr`) so QC + assembler routing is correct.

## Datasets for specific analysis modules (the gap-closing features)

These exercise the modules that go beyond plain profiling. Each is a *reproduction target*:
run the dataset, then check our output against the published / expected result. Start small —
**subsample** big datasets (`seqtk sample -s<seed> reads.fastq <N>`) so a smoke test is minutes,
not hours. When no public set is handy, subsample our bundled `data/simulated_metagenomic_reads.fasta`.

### Ancient DNA — damage authentication (`library: ancient` + `modules.damage`)
The key test is recovering the C→T / G→A post-mortem deamination signal and the authenticity
verdict. **Reproduction target: the i2B group (University of Bergen, Norway)** ancient
metagenome(s) — run their sample through `modules.damage` (mapDamage2) and confirm we reproduce
their reported damage profile / authentication call. *(Fill in the exact ENA/SRA accession from
their publication before bulk download.)*

| Dataset | Notes | Where |
|---|---|---|
| **i2B group, Bergen** ancient metagenome | reproduction target — match their damage/auth result | ENA/SRA accession from their paper (verify) |
| **nf-core/eager test data** | tiny, purpose-built aDNA pipeline test reads (fast smoke test) | github.com/nf-core/test-datasets (eager branch) |
| **aMeta test data** | small ancient metagenomic test set with damage | github.com/NBISweden/aMeta |
| Weyrich et al. 2017 (Neanderthal calculus) | real ancient oral metagenomes, strong damage | ENA **PRJEB14706** (verify) — *Nature* 544:357, doi:10.1038/nature21674 |
| Warinner et al. 2014 (dental calculus) | foundational ancient microbiome | SRA in paper — *Nat Genet* 46:336, doi:10.1038/ng.2906 |

To exercise aDNA QC end-to-end (read-merging), the input must be **FASTQ paired-end** (FASTA
skips QC). `library: ancient` collapses overlapping pairs (fastp `--merge`) before classify/assembly.

### Differential abundance + grouping (`modules.differential`, needs a `group` column)
Needs a study with ≥2 groups and **replication** (≥2 samples/group to run; more for power — a
2-vs-2 permutation p-value floors at ~0.33 and can never reach significance).

| Dataset | Design | Where |
|---|---|---|
| **iHMP IBDMDB** (Lloyd-Price 2019) | IBD vs non-IBD gut, many replicates | SRA **PRJNA398089** — *Nature* 569:655, doi:10.1038/s41586-019-1237-9 |
| **curatedMetagenomicData** CRC cohorts | colorectal-cancer case/control, harmonized | bioconductor.org/packages/curatedMetagenomicData |
| Our bundled reads (quick smoke test) | subsample into ≥2 samples/group with `seqtk` | `data/simulated_metagenomic_reads.fasta` (see below) |

Smoke test we actually ran: 4 pseudo-samples (`seqtk sample` seeds, 2 case / 2 control) →
`classify → abundance → stats → differential` produced `stats/differential_abundance.{tsv,json}`
+ `differential_volcano.png`; with random subsamples (no real signal) it correctly reported 0
significant taxa (FDR control working).

### Amplicon ASV vs OTU (`amplicon.method: asv` → DADA2 | `otu` → VSEARCH)
Use a **16S/18S/ITS mock community** so the expected member list is known.

| Dataset | Notes | Where |
|---|---|---|
| **mockrobiota** (Bokulich 2016) | curated marker-gene mock communities w/ expected taxonomy | github.com/caporaso-lab/mockrobiota |
| **DADA2 tutorial data** (mothur MiSeq SOP) | small paired 16S V4, the canonical ASV test | benjjneb.github.io/dada2/tutorial.html |
| QIIME 2 "Atacama soil" / "moving pictures" | small 16S tutorial sets | docs.qiime2.org |
| **Hymos lab mock dataset** (user-referenced) | smaller mock to subsample for a quick ASV/OTU test | obtain from the Hymos lab; subsample with `seqtk` (verify accession) |
| Zymo 16S (within PRJEB29504) | even/log mock, known composition | ENA **PRJEB29504** |

Set per-sample `library: amplicon`, primers in `amplicon.{fwd_primer,rev_primer}`, and
`amplicon.method: asv` (DADA2, short-read) — long-read amplicon still routes to Emu.

### Biosynthetic gene clusters (`modules.bgc` → antiSMASH; needs assembly)
BGC discovery is easiest to validate on a genome with **known** clusters.

| Dataset | Notes | Where |
|---|---|---|
| *Streptomyces coelicolor* A3(2) genome | ~20+ well-characterized BGCs — gold standard | NCBI assembly **GCF_000203835** |
| **MIBiG** reference BGCs | the known-cluster DB antiSMASH compares against (`--cb-knownclusters`) | mibig.secondarymetabolites.org |
| Soil / marine WGS metagenome | BGC-rich communities for metagenomic BGC mining | SRA (Source=METAGENOMIC, soil/marine) |

antiSMASH needs its reference DBs once: `download-antismash-databases --database-dir <db.antismash>`.

## Bundled validation experiments (no interview — run directly)

Configs live under `config/experiments/`. Batch run + IMRaD papers:

```bash
export PATH="$HOME/miniconda3/envs/metagx-bio/bin:/Library/TeX/texbin:$PATH"
bash scripts/run_experiments_and_papers.sh          # run + report + paper
SKIP_RUN=1 bash scripts/run_experiments_and_papers.sh   # papers only
```

| # | Config | Results dir | Manuscript |
|---|---|---|---|
| 01 | `01-illumina-kraken-confidence-sweep.yaml` | `results/experiments/illumina_kraken2_confidence_sweep/` | `report/manuscript_kraken2_confidence_sensitivity_illumina.pdf` |
| 02 | `02-viral-fasta-confidence-sweep.yaml` | `results/experiments/viral_fasta_confidence_sweep/` | `report/manuscript_kraken2_confidence_sweep_viral_fasta_subsample.pdf` |
| 03 | `03-ont-assembly-confidence-sweep.yaml` | `results/experiments/ont_assembly_reconcile_confidence_sweep/` | `report/manuscript_ont_assembly_reconcile_confidence_sweep.pdf` |
| 04 | `04-pacbio-hifi-confidence-sweep.yaml` | `results/experiments/pacbio_hifi_confidence_sweep/` | `report/manuscript_pacbio_hifi_confidence_sweep.pdf` |
| 05 | `05-pacbio-clr-confidence-sweep.yaml` | `results/experiments/pacbio_clr_confidence_sweep/` | `report/manuscript_pacbio_clr_confidence_sweep.pdf` |
| 06 | `06-multisample-alpha-beta-diversity.yaml` | `results/experiments/multisample_alpha_beta_diversity/` | `report/manuscript_multisample_alpha_beta_diversity.pdf` |
| 07 | *(standalone)* `scripts/compare_technologies.py` | `results/experiments/cross_technology_comparison/` | `manuscript_cross_technology_comparison.pdf` |
| 08 | *(standalone)* `scripts/run_platform_assemblies.sh` + `scripts/compare_platforms.py` | `results/experiments/cross_platform_comparison/` | `manuscript_assembly_comparison.pdf` |

Profiling / differential baselines (also produce papers):

| Config | Results | Manuscript |
|---|---|---|
| `illumina-sim.yaml` | `results/illumina_sim/` | `manuscript_illumina_pe_qc_profiling.pdf` |
| `illumina-sim-diff.yaml` | `results/illumina_sim_diff/` | `manuscript_illumina_case_control_differential.pdf` |
| `diff-demo.yaml` | `results/diff_demo/` | `manuscript_four_sample_differential_abundance.pdf` |
| `pacbio-hifi-sim.yaml` | `results/pacbio_hifi_sim/` | `manuscript_pacbio_hifi_profiling.pdf` |
| `pacbio-clr-sim.yaml` | `results/pacbio_clr_sim/` | `manuscript_pacbio_clr_profiling.pdf` |

**What each experiment tests**

- **01 Illumina confidence sweep** — wgsim PE reads; Kraken2 `--confidence` ∈ {0, 0.05, 0.1, 0.2, 0.5}; stable until 0.5 where ~25% reads drop out of top species.
- **02 Viral FASTA subsample sweep** — 20% seqtk subsample of bundled FASTA; confidence ∈ {0, 0.1, 0.5}; cheaper smoke test of sweep machinery.
- **03 ONT assembly + reconcile** — nanopore sim reads + MEGAHIT + CAT/Binning reconciliation; confidence ∈ {0, 0.05, 0.1, 0.2}.
- **04 PacBio HiFi sweep** — SimLoRD HiFi; confidence ∈ {0, 0.05, 0.1, 0.2}; classified fraction falls sharply as confidence rises.
- **05 PacBio CLR sweep** — noisy CLR; use low thresholds {0, 0.01, 0.02} — at ≥0.05 almost all reads become unclassified.
- **06 Multi-sample diversity** — four subsampled FASTA pseudo-samples (`data/multi/samp{A-D}.fasta`); alpha/beta diversity + PCoA (`stats/`).
- **07 Cross-technology comparison** — *standalone, not a `metagx run` config.* Compares the
  ONT / Illumina / PacBio-CLR / PacBio-HiFi profiles already produced for the same 30-genome
  reference. Run `.venv/bin/python scripts/compare_technologies.py`. Headline result: the
  uniform-coverage simulators (wgsim, SimLoRD) impose a strong **genome-length bias**
  (Spearman ρ 0.90–0.96 between genome length and assigned reads; the 190 kb *Glossina*
  genome captures 44–60% of reads), while the ONT set (≈equal reads/genome) shows none
  (ρ = 0.07). Treat the ONT vs uniform-coverage sets as **different ground-truth communities**
  when benchmarking accuracy.
- **08 Cross-platform assembly comparison** — *standalone.* Assembles each platform
  (MEGAHIT short / Flye long) and compares contiguity, reference breadth, and read concordance.
  Run `bash scripts/run_platform_assemblies.sh` then `scripts/compare_platforms.py`. Findings:
  long reads win **contiguity** (PacBio CLR N50 30.6 kb, one 143 kb near-complete genome vs
  Illumina N50 1.5 kb / 373 contigs); deep even short reads win **breadth** (Illumina 90%,
  30/30 genomes); PacBio HiFi at ~5× (415 reads) **failed to assemble at all** — assembly has a
  coverage floor classification does not. The comparison is **depth-confounded** (20k/8k/415/415
  reads); read it qualitatively. `compare_platforms.py` is the **generalized, manifest-driven
  successor to `compare_technologies.py`** — use it for any same-sample multi-platform set.

Multi-sample data is **not** subsampled from a single file when you already have distinct samples — here we use four independent `seqtk` draws (seeds 1–4) from `data/simulated_metagenomic_reads.fasta` to mimic biological replicate spread while keeping ground-truth taxa comparable.

## Key references (methods/benchmarks)
- Nicholls SM, et al. *Ultra-deep, long-read nanopore sequencing of mock microbial community standards.* GigaScience. 2019;8:giz043.
- Meyer F, et al. *Critical Assessment of Metagenome Interpretation (CAMI II).* Nature Methods. 2022;19:429-440.
- Portik DM, Brown CT, Pierce-Ward NT. *Evaluation of taxonomic classification and profiling methods for long-read shotgun metagenomic sequencing datasets.* BMC Bioinformatics. 2022;23:541.
- Bickhart DM, et al. *Generating lineage-resolved, complete metagenome-assembled genomes from complex microbial communities.* Nature Biotechnology. 2022;40:711-719.
- Lloyd-Price J, et al. *Multi-omics of the gut microbial ecosystem in inflammatory bowel diseases.* Nature. 2019;569:655-662. (differential-abundance case/control reference cohort)
- Bokulich NA, et al. *mockrobiota: a public resource for microbiome bioinformatics benchmarking.* mSystems. 2016;1(5):e00062-16. (amplicon mock communities)
- Callahan BJ, et al. *DADA2: high-resolution sample inference from Illumina amplicon data.* Nature Methods. 2016;13:581-583. (ASV inference)
- Blin K, et al. *antiSMASH 7.0.* Nucleic Acids Research. 2023;51(W1):W46-W50. (BGC discovery)
- Fellows Yates JA, et al. *Reproducible, portable, and efficient ancient genome reconstruction with nf-core/eager.* PeerJ. 2021;9:e10947. (aDNA pipeline + test data)
- Jónsson H, et al. *mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters.* Bioinformatics. 2013;29(13):1682-1684. (damage authentication)
