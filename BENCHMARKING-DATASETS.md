# How to build valid benchmarking datasets for cross-platform comparison

A comparison is only as good as the dataset behind it. If you want to see **real** comparative
features between sequencing platforms (ONT vs Illumina vs PacBio), the dataset has to be built
so that the *only* thing that differs is the variable you are studying. Everything else must
be held constant. This note is the checklist for getting that right — written from the
mistakes our own bundled datasets make (see the self-assessment at the end).

> **The one rule:** *change one variable at a time.* To attribute a difference to "platform",
> every other input — the genomes present, their relative abundance, and the sequencing depth —
> must be identical across the datasets you compare. Otherwise your "platform effect" is
> actually a depth effect or an abundance-model effect wearing a platform costume.

## The three things that must match

### 1. Shared reference and known truth
All datasets must derive from the **same reference genomes**, and you must record the ground
truth so accuracy is measurable:

- the exact set of genomes (accessions) and their **lengths**,
- the **intended relative abundance** of each genome (the truth vector),
- provenance (which simulator/run produced each file).

Without a truth vector you can measure *agreement between platforms* but not *accuracy*
(precision/recall, breadth, abundance error). Keep a `truth.tsv`: `accession, length_bp,
rel_abundance`.

### 2. Matched abundance model — the subtle one
This is the trap that silently breaks most home-made benchmarks. Read simulators offer two
fundamentally different notions of "abundance", and they produce **different communities**:

| Model | What it fixes | Effect | Tools |
|---|---|---|---|
| **Equal coverage** (depth ×) | reads per *base* per genome | read count ∝ **genome length** → long genomes dominate | wgsim, SimLoRD, ART (per-genome ×) |
| **Equal reads** (count) | reads per *genome* | even read share regardless of length | per-genome read-count sims |

In our own data this produced a 12× distortion: under equal-coverage simulation the single
190 kb *Glossina* genome (44.5 % of total reference length) captured 44–60 % of reads, while
under equal-reads simulation it captured ~8 %. **Same genomes, completely different profile.**

For a valid cross-platform comparison you must pick **one** abundance model and apply the
**identical per-genome abundance vector** to every platform. If platform A is simulated
equal-coverage and platform B equal-reads, you are comparing two different communities, not two
platforms. (Either model is fine — equal-reads is more even and easier to reason about;
equal-coverage mimics what real shotgun sequencing of equal-mass DNA actually does. Just be
consistent and write down which.)

### 3. Matched sequencing depth
Depth must be comparable across platforms, normalized **the same way**:

- **Equal total bases** (Σ read length) — the fair choice for **assembly** comparisons, since
  assembly contiguity and breadth scale with bases of coverage, not read count.
- **Equal read count** — acceptable for **classification** comparisons (each read is one
  classification event), but unfair to short reads for assembly.

Mismatched depth is the most common confound. Our assembly experiment compared 20,000 / 8,000 /
415 / 415 reads and consequently measured *depth × platform*, not platform: PacBio HiFi at ~5×
failed to assemble entirely, which says more about depth than about HiFi. Decide on a target
(e.g. "30× per genome" or "20,000 reads") and hit it on every platform.

## Recipe: a clean simulated benchmark

```bash
# 0. Reference + truth ----------------------------------------------------------
#    genomes.fasta = the shared community; write truth.tsv (accession,length,rel_abundance)

# 1. Pick ONE abundance vector and turn it into per-genome read counts ----------
#    e.g. equal-reads: N_total / n_genomes reads each, OR a deliberate skew you record.
#    For equal-coverage, per-genome reads = depth * genome_len / mean_read_len.

# 2. Simulate EACH platform from the SAME per-genome read counts ----------------
#    Illumina  : wgsim / ART     (short, ~150 bp, low error)
#    ONT       : pbsim3 / Badread (long, high error)
#    PacBio HiFi: pbsim3 multipass + ccs / SimLoRD multipass (long, low error)
#    PacBio CLR : pbsim3 1-pass / SimLoRD 1-pass (long, high error)
#    -> drive every simulator from the same per-genome counts so abundance is identical.

# 3. Normalize depth across platforms -------------------------------------------
seqtk sample -s42 platform.fastq <N_READS>     # equal read count, OR
#   downsample to equal total bases with a short script.

# 4. Record provenance ----------------------------------------------------------
#    seed, simulator version, command, target depth, abundance model -> README per file.
```

For a quick smoke test you can subsample one deep FASTA into pseudo-replicates with
`seqtk sample -s<seed>`, but pseudo-replicates share the same reads — fine for plumbing tests,
not for measuring real biological or platform variance.

## Better: a real multi-platform mock community
The strongest benchmark is one physical DNA sample sequenced on every platform — the abundance
is then **intrinsically matched** (it is the same molecules), and any per-platform skew (GC
bias, length bias, chimera rate) is a *genuine* platform feature you want to measure rather than
a simulation artifact.

### ZymoBIOMICS — exact links (the de-facto multi-platform mock)
Product info: ZymoBIOMICS Microbial Community Standard **D6300** (even,
https://www.zymoresearch.com/products/zymobiomics-microbial-community-standard) and Gut
Microbiome Standard **D6331** (https://www.zymoresearch.com/products/zymobiomics-gut-microbiome-standard).

| Platform | Community | Accession(s) | Where |
|---|---|---|---|
| Illumina + ONT (GridION/PromethION) | D6300 even / D6310 log | ENA study **PRJEB29504**; ONT FASTQ even `ERR3152364` (GridION) / `ERR3152365` (PromethION), log `ERR3152366` / `ERR3152367` | https://lomanlab.github.io/mockcommunity/ · https://www.ebi.ac.uk/ena/browser/view/PRJEB29504 |
| PacBio HiFi | D6331 gut | NCBI **PRJNA680590**: standard `SRR13128014`, low `SRR13128013`, ultra-low `SRR13128012` | https://github.com/PacificBiosciences/pb-metagenomics-tools/blob/master/docs/PacBio-Data.md |

> **Matched-community caveat:** the public ONT + Illumina reads are the **D6300/D6310**
> standards; the public PacBio HiFi reads are the **D6331 gut** standard — *different
> communities*. So ZymoBIOMICS public data gives a clean **ONT-vs-Illumina** same-sample
> comparison (identical standards, Loman Lab), but a head-to-head that *includes* PacBio
> requires either sequencing one standard yourself on all three platforms, or treating the
> HiFi D6331 set as its own benchmark. This is exactly the "match the community" rule applied
> to real data.

Other options:
- **Loman Lab mock** (Nicholls 2019, ENA `PRJEB29504`) — ONT + Illumina of the same standards,
  with ground truth. The recommended *matched* pair.
- **CAMI / CAMI II** (data.cami-challenge.org) — simulated with gold-standard truth profiles and
  assemblies; abundance and depth are controlled by design.
- **Portik et al. 2022** (*BMC Bioinformatics* 23:541) — mock ONT + PacBio HiFi with profiling
  truth, purpose-built for long-read method benchmarking.

## Validity checklist

- [ ] Same reference genomes across all platforms (record accessions + lengths).
- [ ] A written truth vector (`truth.tsv`) of intended per-genome relative abundance.
- [ ] **One** abundance model (equal-reads OR equal-coverage), applied identically to all.
- [ ] Depth normalized the **same way** for all platforms (equal bases for assembly; equal reads
      for classification) and at a level that actually supports the analysis (assembly needs
      ≥ ~10–30× per genome; classification tolerates far less).
- [ ] Provenance recorded per file: simulator + version, seed, command, target depth.
- [ ] Analysis run identically across platforms (same DB, same params, same pipeline), varying
      only the platform-routed steps (QC, assembler, minimap2 preset).
- [ ] Confounds named in the writeup. If something is not matched, say so and read results
      qualitatively.

## Self-assessment of the bundled metagx datasets

Honest accounting of what our shipped test data does and doesn't support:

| Property | ONT (`simulated_metagenomic_reads.fasta`) | Illumina (wgsim) | PacBio (SimLoRD) |
|---|---|---|---|
| Shared reference (30 genomes) | ✅ | ✅ | ✅ |
| Abundance model | equal-reads (even) | equal-coverage (length-biased) | equal-coverage |
| Depth | ~20,000 reads | ~8,000 reads | ~415 reads (~5×) |

**Good for:** exercising each platform's pipeline path end-to-end (QC routing, assembler
routing, classification, the comparison machinery). This is their purpose and they do it well.

**Not valid for, as-is:** a head-to-head *platform accuracy* ranking — the ONT set uses a
different abundance model than the uniform-coverage sets, and depths differ by ~48×. Experiments
07–08 therefore report results **qualitatively** and label the confounds explicitly.

**To turn them into a valid benchmark:** re-simulate all four platforms from one per-genome
abundance vector at one matched depth target (e.g. 20,000 reads each, equal-reads model), then
re-run `metagx compare`. That single change converts the suite from "pipeline smoke tests" into
"platform comparison with measurable accuracy."

## See also
- `DATASETS.md` — dataset catalog + experiments 07 (classification) and 08 (assembly).
- `scripts/compare_platforms.py` / `metagx compare` — the comparison engine.
- `results/experiments/cross_technology_comparison/` and `…/cross_platform_comparison/` —
  the comparisons that surfaced these lessons.
