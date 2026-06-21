# metagx — design notes: multi-domain awareness & reference-database scaling

> Internal design/reasoning document. Covers (1) making the pipeline aware of viruses,
> prokaryotes, and eukaryotes, and (2) what changes — statistically and methodologically —
> when the classification database grows from a handful of references to comprehensive
> (e.g. ~14k+ RefSeq viral genomes) for real samples where only a few taxa are present.
> Inline markers `[n]` reference the list at the end (papers + the tools/algorithms we use).

---

## 0. Context

The current pipeline: read QC (fastp / porechop_abi+chopper) `[12,17,18]` → taxonomic
classification (kraken2 `[1]` + Bracken `[2]`) → optional assembly (MEGAHIT `[4]` /
Flye–metaFlye `[5,6]`) → binning (MetaBAT2 `[7]`) → a **reconcile** step that classifies
contigs (kraken2 + CAT per-ORF voting `[10]`) and joins them to per-contig coverage
(minimap2 `[8]` / SAMtools `[9]`) and the read-level calls. Orchestrated with Snakemake
`[11]`. Validated on 30 unrelated viral genomes + simulated ONT reads.

Two questions arise when moving from that toy setup toward real work.

---

## 1. Domain awareness (virus / prokaryote / eukaryote)

### 1.1 Principle

**Domain barely changes *read classification*, but it changes almost everything *after*
assembly** — and real metagenomes are *mixed*. So the design is not "pick a domain" but
"classify broadly, then **route each contig to domain-appropriate tools**."

- Read-level kraken2/Bracken `[1,2]` are domain-agnostic *provided the reference DB and its
  NCBI taxonomy `[20]` span the target domains*. A single broad index (e.g. the standard
  PlusPF set, or a custom DB built with `kraken2-build --download-library viral|bacteria|
  fungi|protozoa`) classifies bacteria, archaea, viruses, fungi, and protozoa through one
  LCA pass. Domain intelligence therefore is **not** needed at the read step.
- Domain intelligence **is** needed at the contig/genome layer: gene calling, binnability,
  completeness/quality estimation, and genome-level taxonomy all differ by domain.

### 1.2 Where domain actually matters

| Stage | Virus | Prokaryote (bac/arch) | Eukaryote |
|---|---|---|---|
| Read classify | kraken2/Bracken `[1,2]` (DB must include viral) | same | same; or EukDetect `[16]` for sensitive euk read detection |
| Gene prediction | Prodigal `[14]` (ok) | Prodigal `[14]` | ❌ Prodigal is prokaryotic — needs a eukaryote-aware caller |
| Binning | ✗ viruses don't bin by tetranucleotide+coverage | ✓ MetaBAT2 `[7]` | hard; separate euk contigs first (EukRep `[15]`) |
| Contig/genome ID | geNomad `[21]` (identify + ICTV taxonomy) | GTDB-Tk `[13]` | EukRep `[15]` |
| Completeness/quality | CheckV `[22]` | CheckM2 `[23]` | EukCC `[24]` |

Two facts behind the table:
- **Binning is prokaryote-centric.** MetaBAT2 `[7]` groups contigs by tetranucleotide
  frequency + coverage covariance across samples; viruses (small, often one contig =
  one genome) and eukaryotes (large, repeat-rich, intron-containing) do not behave that
  way. This is why our viral test produced **0 bins** and why BAT/GTDB-Tk were not
  applicable — contig-level is the correct granularity for viruses.
- **Gene calling diverges.** Prodigal `[14]` (used by CAT `[10]`) is built for prokaryotic
  ORFs and works acceptably on viruses, but is wrong for eukaryotes — a euk path needs a
  dedicated gene model.

### 1.3 Mixed-community routing (the real design)

A shotgun metagenome contains all domains at once. The robust pattern:

```
reads ──kraken2/Bracken (broad DB)──> community profile (all domains)
contigs ──┬── geNomad [21] ───────────> viral/plasmid contigs ─> CheckV [22]
          ├── EukRep  [15] ───────────> eukaryotic contigs   ─> EukCC [24]
          └── remainder (prokaryotic) ─> MetaBAT2 [7] ─> GTDB-Tk [13] + CheckM2 [23]
```

geNomad `[21]` and EukRep `[15]` act as **contig sorters**; each fraction then flows to its
domain's completeness + taxonomy tools. Reads stay a single broad pass.

### 1.4 Incorporation into metgx (fits the existing pattern)

`domain` becomes a routing dimension exactly like `platform` already is for QC/assembler:

1. A `domains` (or `target`) config field — e.g. `domains: [viral, prokaryote, eukaryote]`.
   It controls (a) which reference libraries `build-db` includes, and (b) which
   post-assembly modules run.
2. A **contig-router** step (geNomad `[21]` + EukRep `[15]`) after assembly.
3. New per-tool registries — `genomad`, `checkv`, `checkm2`, `gtdbtk`, `eukrep`, `eukcc` —
   each adds its flags once (the single-source-of-truth pattern already used for
   kraken2/fastp/flye/etc.).
4. A euk-aware gene-calling branch so CAT/annotation isn't limited to Prodigal `[14]`.

Recommended first slice (matches current viral testing): **geNomad `[21]` + CheckV `[22]`**
attached to `reconcile`, then add the prokaryote (GTDB-Tk `[13]` + CheckM2 `[23]`) and
eukaryote (EukRep `[15]` + EukCC `[24]`) branches by the same shape.

---

## 2. Scaling the reference DB: 30 → 14,000+ viral genomes

### 2.1 Computational cost: negligible

Viral genomes are tiny (~10–200 kb; ~14k × ~50 kb ≈ 0.7 Gb of sequence → a few-GB kraken2
index, modest RAM). One would not hand-build it: `kraken2-build --download-library viral`
pulls RefSeq viral `[20]` (~15k genomes) **with real NCBI taxonomy**, then `bracken-build`
at the read length. (Large RAM is a concern for *bacterial* indices, not viral.)

### 2.2 Statistical consequences (the real change)

The clean species-level results on 30 *unrelated* genomes are partly an artifact of that
setup. With a comprehensive, internally-related reference set:

- **LCA push-up.** Related viruses share k-mers (conserved domains, strains, species
  complexes). kraken2 assigns each read to the **lowest common ancestor** of every genome
  a k-mer hits `[1]`, so reads migrate from species toward genus/family. Nasko et al.
  showed directly that **growing the reference database shifts k-mer LCA identification to
  higher ranks / reduces species-level accuracy** `[25]` — exactly this effect.
- **% classified rises, precision falls.** More references ⇒ more reads find *some* match,
  including spurious ones; raw "% classified" becomes a poor quality proxy. Benchmarks of
  metagenomic classifiers document this sensitivity/precision trade-off and the resulting
  false-positive burden `[26,27]`.
- **A false-positive long tail.** When "not all viruses are present" (always true), a
  near-complete DB recruits stray reads (chance k-mers, conserved regions, contaminants)
  to absent taxa, producing many low-abundance phantom species `[26,27]`.
- **Bracken starts doing real work.** With 30 unrelated genomes it was nearly idle; at
  scale it is essential — it redistributes ambiguous/higher-rank assignments back down to
  species using the per-genome k-mer distribution `[2]`.

### 2.3 Controls (mostly already exposed)

- **Confidence sweep** (the "k-dense" matrix) — kraken2 `--confidence` `[1]`; sweeping it is
  precisely how to find where the FP tail collapses. Tune *up* vs the toy run.
- **`--minimum-hit-groups`** `[1]` — require independent evidence before a call.
- **Bracken `-t` threshold** `[2]` — drop taxa below a minimum read count.

### 2.4 The key interpretation upgrade: presence ≠ read count

At 14k references, depth/read-count is misleading — a single conserved gene can recruit
reads to a virus that is absent. The defensible "is it really present?" signals:

1. **Breadth / horizontal coverage** — what *fraction of the reference genome* is covered,
   not depth at one locus. Breadth-of-coverage thresholds are the standard way to separate
   true presence from spurious recruitment (e.g. inStrain uses breadth for population
   detection `[28]`; MIUViG reporting standards for uncultivated viral genomes formalize
   completeness/coverage expectations `[29]`).
2. **Assembly + CheckV** — a virus that assembles into a substantial, high-completeness
   contig (CheckV `[22]`) is far stronger evidence than scattered read hits.
3. **Read ↔ contig concordance** — the existing `reconcile` output.

### 2.5 Why this makes `reconcile` *more* valuable at scale

The expected pattern flips. On 30 unrelated genomes, reads and contigs agreed 25/25.
At 14k related references, **reads become ambiguous (LCA push-up) while long contigs — and
CAT's per-ORF voting `[10]` — still resolve species.** So the signal becomes:

> reads → a noisy abundance tail; contigs → the trustworthy identities.

That read-vs-contig split is exactly what `reconcile` surfaces. Adding a
**breadth-of-coverage column** (from the existing minimap2/SAMtools mapping `[8,9]`) and
**CheckV completeness** `[22]` would turn a depth-ranked "list of maybes" into a defensible
virome call.

---

## 3. Library strategy: assembly applies to WGS, not amplicon

De novo assembly assumes reads are **random fragments sampled across whole genomes** (WGS
shotgun). **Amplicon** sequencing (16S/18S rRNA, ITS) targets a single marker locus, so the
data is millions of near-identical copies of one short region — assembly, binning, contig
taxonomy, and MAG recovery are meaningless on it. The pipeline treats `library`
(`wgs|amplicon`) as a per-sample routing dimension:

- **WGS** → the shotgun path: QC → kraken2/Bracken → assembly → binning → reconcile →
  domain taxonomy.
- **Amplicon** → a dedicated path: primer removal (cutadapt `[30]`) → OTU clustering for
  short reads (VSEARCH `[31]`; or ASV denoising à la DADA2 `[32]` / QIIME 2 `[33]`) or
  full-length 16S abundance for long reads (Emu `[34]`), against a marker-gene database such
  as SILVA `[35]`. Assembly-based modules are skipped; an all-amplicon run that enables them
  is rejected with a clear error.

Read-level kraken2/Bracken *can* still be run on amplicon reads, but it is a rough
substitute: classifying a single conserved region against a whole-genome database inflates
ambiguity, so a marker-gene DB + ASV/OTU methods are preferred (the tool warns accordingly).
Per-sample routing means mixed WGS + amplicon runs are handled correctly.

Other non-shotgun library types (targeted/hybrid-capture, metatranscriptomics) share the
same caveat — the random-fragment assumption behind assembly may not hold; `library` is the
extension point for handling them.

## 4. Summary of design implications

1. Keep **one broad read-classification pass**; put domain logic in the contig/genome layer.
2. Add a `domains` routing dimension + contig sorters (geNomad `[21]`, EukRep `[15]`) →
   domain-specific completeness/taxonomy (CheckV `[22]`; GTDB-Tk `[13]`/CheckM2 `[23]`;
   EukCC `[24]`); fix the Prodigal-only `[14]` gene-calling assumption for eukaryotes.
3. For large viral DBs: trivial to build (RefSeq viral `[20]`), but tighten FP controls
   (confidence sweep `[1]`, min-hit-groups `[1]`, Bracken `-t` `[2]`) and **shift the
   presence criterion to breadth + assembly/CheckV + read-contig concordance** `[22,28,29]`.
4. `reconcile` is the right home for this and gains importance as references scale.
5. Gate assembly-based analysis on `library: wgs`; route amplicon to primer-trim + OTU/Emu.
   Assembly's random-fragment assumption is violated by marker-gene data.

---

## References

### Pipeline tools we use
1. Wood DE, Lu J, Langmead B. *Improved metagenomic analysis with Kraken 2.* Genome Biology. 2019;20:257. doi:10.1186/s13059-019-1891-0
2. Lu J, Breitwieser FP, Thielen P, Salzberg SL. *Bracken: estimating species abundance in metagenomics data.* PeerJ Computer Science. 2017;3:e104. doi:10.7717/peerj-cs.104
4. Li D, Liu CM, Luo R, Sadakane K, Lam TW. *MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph.* Bioinformatics. 2015;31(10):1674-1676. doi:10.1093/bioinformatics/btv033
5. Kolmogorov M, Yuan J, Lin Y, Pevzner PA. *Assembly of long, error-prone reads using repeat graphs (Flye).* Nature Biotechnology. 2019;37:540-546. doi:10.1038/s41587-019-0072-8
6. Kolmogorov M, et al. *metaFlye: scalable long-read metagenome assembly using repeat graphs.* Nature Methods. 2020;17:1103-1110. doi:10.1038/s41592-020-00971-x
7. Kang DD, et al. *MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies.* PeerJ. 2019;7:e7359. doi:10.7717/peerj.7359
8. Li H. *Minimap2: pairwise alignment for nucleotide sequences.* Bioinformatics. 2018;34(18):3094-3100. doi:10.1093/bioinformatics/bty191
9. Danecek P, et al. *Twelve years of SAMtools and BCFtools.* GigaScience. 2021;10(2):giab008. doi:10.1093/gigascience/giab008
10. von Meijenfeldt FAB, Arkhipova K, Cambuy DD, Coutinho FH, Dutilh BE. *Robust taxonomic classification of uncharted microbial sequences and bins with CAT and BAT.* Genome Biology. 2019;20:217. doi:10.1186/s13059-019-1817-x
11. Mölder F, et al. *Sustainable data analysis with Snakemake.* F1000Research. 2021;10:33. doi:10.12688/f1000research.29032.2
12. Chen S, Zhou Y, Chen Y, Gu J. *fastp: an ultra-fast all-in-one FASTQ preprocessor.* Bioinformatics. 2018;34(17):i884-i890. doi:10.1093/bioinformatics/bty560
14. Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. *Prodigal: prokaryotic gene recognition and translation initiation site identification.* BMC Bioinformatics. 2010;11:119. doi:10.1186/1471-2105-11-119
17. Bonenfant Q, Noé L, Touzet H. *Porechop_ABI: discovering unknown adapters in Oxford Nanopore Technology sequencing reads for downstream analysis.* Bioinformatics Advances. 2023;3(1):vbac085. doi:10.1093/bioadv/vbac085
18. De Coster W, Rademakers R. *NanoPack2: population-scale evaluation of long-read sequencing data (chopper).* Bioinformatics. 2023;39(5):btad311. doi:10.1093/bioinformatics/btad311
19. Buchfink B, Reuter K, Drost HG. *Sensitive protein alignments at tree-of-life scale using DIAMOND.* Nature Methods. 2021;18:366-368. doi:10.1038/s41592-021-01101-x
20. O'Leary NA, et al. *Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation.* Nucleic Acids Research. 2016;44(D1):D733-745. doi:10.1093/nar/gkv1189

### Domain-specific genome/contig tools (proposed)
13. Chaumeil PA, Mussig AJ, Hugenholtz P, Parks DH. *GTDB-Tk v2: memory friendly classification with the Genome Taxonomy Database.* Bioinformatics. 2022;38(23):5315-5316. doi:10.1093/bioinformatics/btac672
15. West PT, Probst AJ, Grigoriev IV, Thomas BC, Banfield JF. *Genome-reconstruction for eukaryotes from complex natural microbial communities (EukRep).* Genome Research. 2018;28:569-580. doi:10.1101/gr.228429.117
16. Lind AL, Pollard KS. *Accurate and sensitive detection of microbial eukaryotes from whole metagenome shotgun sequencing (EukDetect).* Microbiome. 2021;9:58. doi:10.1186/s40168-021-01015-y
21. Camargo AP, et al. *Identification of mobile genetic elements with geNomad.* Nature Biotechnology. 2024;42:1303-1312. doi:10.1038/s41587-023-01953-y
22. Nayfach S, Camargo AP, Schulz F, Eloe-Fadrosh E, Roux S, Kyrpides NC. *CheckV assesses the quality and completeness of metagenome-assembled viral genomes.* Nature Biotechnology. 2021;39:578-585. doi:10.1038/s41587-020-00774-7
23. Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. *CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning.* Nature Methods. 2023;20:1203-1212. doi:10.1038/s41592-023-01940-w
24. Saary P, Mitchell AL, Finn RD. *Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC.* Genome Biology. 2020;21:244. doi:10.1186/s13059-020-02155-4

### Methodological / conceptual basis for the reasoning
25. Nasko DJ, Koren S, Phillippy AM, Treangen TJ. *RefSeq database growth influences the accuracy of k-mer-based lowest common ancestor species identification.* Genome Biology. 2018;19:165. doi:10.1186/s13059-018-1554-6  — *basis for §2.2 "LCA push-up".*
26. Ye SH, Siddle KJ, Park DJ, Sabeti PC. *Benchmarking metagenomics tools for taxonomic classification.* Cell. 2019;178(4):779-794. doi:10.1016/j.cell.2019.07.010  — *sensitivity/precision & false-positive trade-offs.*
27. McIntyre ABR, et al. *Comprehensive benchmarking and ensemble approaches for metagenomic classifiers.* Genome Biology. 2017;18:182. doi:10.1186/s13059-017-1299-7  — *false-positive control, multi-tool agreement.*
28. Olm MR, et al. *inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains.* Nature Biotechnology. 2021;39:727-736. doi:10.1038/s41587-020-00797-0  — *breadth-of-coverage as a presence/strain criterion (§2.4).*
29. Roux S, et al. *Minimum Information about an Uncultivated Virus Genome (MIUViG).* Nature Biotechnology. 2019;37:29-37. doi:10.1038/nbt.4306  — *reporting standards: completeness/coverage for viral genome calls (§2.4).*

### Amplicon-branch tools (§3)
30. Martin M. *Cutadapt removes adapter sequences from high-throughput sequencing reads.* EMBnet.journal. 2011;17(1):10-12. doi:10.14806/ej.17.1.200
31. Rognes T, Flouri T, Nichols B, Quince C, Mahé F. *VSEARCH: a versatile open source tool for metagenomics.* PeerJ. 2016;4:e2584. doi:10.7717/peerj.2584
32. Callahan BJ, et al. *DADA2: high-resolution sample inference from Illumina amplicon data.* Nature Methods. 2016;13:581-583. doi:10.1038/nmeth.3869
33. Bolyen E, et al. *Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2.* Nature Biotechnology. 2019;37:852-857. doi:10.1038/s41587-019-0209-9
34. Curry KD, et al. *Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data.* Nature Methods. 2022;19:845-853. doi:10.1038/s41592-022-01520-4
35. Quast C, et al. *The SILVA ribosomal RNA gene database project.* Nucleic Acids Research. 2013;41(D1):D590-D596. doi:10.1093/nar/gks1219
