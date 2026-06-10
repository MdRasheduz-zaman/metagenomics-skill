# metagx — comprehensive critique & roadmap (internal)

Honest assessment of the tool's coverage of metagenomics, with a prioritized roadmap.
Status tags below: ✅ done & run · 🟡 implemented, dry-run only · ⬜ planned.

## Executive verdict

metagx is an unusually well-architected profiling + assembly **engine** with excellent DX
(registry single-source-of-truth, multi-client, provenance/Methods, reconcile layer). With
Tier 1 (cross-sample stats, normalization, host removal, tests/CI), Tier 2 (functional/AMR,
bin refinement→DAS_Tool/dRep, metaSPAdes/hybrid, second-classifier consensus, MultiQC/Krona)
and Tier 3 (ancient-DNA damage authentication, decontam/controls, strain-level, per-sample
Bracken length, containers/conda-lock) now implemented, it spans a **study-level analysis
platform**. Remaining gaps are adjacent subfields explicitly out of scope here
(metatranscriptomics, viral-host linkage, long-read-specific polishing). The architecture is
ready for all of it — the gap is breadth of coverage, not design. (Tiers 2-3 are implemented +
dry-run-validated; full execution needs external tools/DBs that aren't installed here.)

## Strengths
- Registry-as-truth keeps interview, validation, CLI/MCP, and command construction in sync.
- Per-sample routing dimensions: platform, layout, format, **library (wgs/amplicon)**, domain.
- Reconcile (read↔contig concordance, breadth, read-accuracy, CAT cross-check) — original.
- Provenance + Copy-as-Methods, presets, conda self-provisioning, SLURM, multi-client surface.

## CRITICAL gaps (limit scientific conclusions)
1. ✅ **Cross-sample statistics & diversity** — α (Shannon/Simpson/evenness), β (Bray–Curtis),
   ordination (PCoA) [Tier 1], **and differential abundance** (`modules.differential`: CLR +
   permutation + BH FDR, pure-Python) [2026-06-10, RUN on subsampled data]. Comparative result
   now available.
2. **Compositional normalization** — TSS/CLR/rarefaction; raw counts are compositional.
3. **Host/contaminant removal as a first-class pre-classification QC step** (human T2T;
   accuracy + PHI). Currently only feeds filtered_assembly.
4. **Decontamination / controls** — negative/blank-aware filtering (decontam) for low-biomass.
5. **Single classifier** — kraken2-only inherits DB-completeness FP behavior; add a marker
   (MetaPhlAn) or protein (Kaiju) classifier for consensus.

## IMPORTANT gaps
6. **Functional / AMR / resistome** — HUMAnN (pathways), Bakta/Prokka + eggNOG (MAG
   annotation), AMRFinderPlus / RGI-CARD / ABRicate (resistance) [Tier 2], ✅ **antiSMASH
   (BGCs)** (`modules.bgc`) [2026-06-10, dry-run].
7. **Binning** — single binner, single sample. Add MaxBin2+CONCOCT → **DAS_Tool**, **dRep**
   dereplication, multi-sample co-binning; add **metaSPAdes** (Illumina) + hybrid assembly.
8. **Strain-level** — inStrain / StrainPhlAn (SNV/microdiversity).
9. **Aggregate QC + read-flow accounting** (MultiQC-style) and **visualization** (Krona,
   composition barplots, ordination plots).
10. **Duplicate removal** as a deliberate, quantification-aware step.

## Ancient DNA / paleogenomics (distinct subfield)
- **Damage authentication is mandatory**: C→T/G→A deamination at fragment ends distinguishes
  ancient from modern contamination — metaDMG / mapDamage2 / PMDtools.
- Ultra-short degraded reads (~30–70 bp) → short-read params; **read merging/collapsing** of
  overlapping pairs (fastp --merge / leeHom).
- Authentication metrics: read-length dist, edit distance, **breadth** (already have) + damage.
- Decontamination paramount. Interop: nf-core/eager, aMeta, HOPS/MALT.
- Fits a new `library: ancient` routing value.

## Engineering / robustness
- **No tests, no CI** — biggest engineering risk; pure-Python core is highly unit-testable
  (already hit Bracken `echo` + `platform`-shadow bugs). Add pytest + GitHub Actions.
- **No containers / no pinned envs** (conda-lock) — not bit-reproducible. (Singularity deferred.)
- Per-read `.kraken` persistence (for read-accuracy) should be opt-in (disk at scale).
- Config-time file-existence + DB↔intent compatibility checks would fail faster.
- Known caveats: Bracken read-length global (not per-sample); porechop --ab_initio aborts on
  IUPAC; most domain/amplicon paths dry-run-only.

## Prioritized roadmap & status

**Tier 1 (study tool) — ✅ DONE & RUN (2026-06-09):**
- ✅ Cross-sample stats: `modules.stats` → α (Shannon/Simpson/richness/Pielou), β
  (Bray–Curtis), PCoA, TSS + CLR matrices, composition barplot + PCoA plot
  (`metagx/diversity.py`, `rules/stats.smk`). Run on 3 pseudo-samples.
- ✅ Compositional normalization (TSS + CLR) — part of stats.
- ✅ Host removal as first-class pre-classification QC step (`host_removal.genome`,
  minimap2 + keep-unmapped; routes the whole pipeline). Run (764 host reads depleted).
- ✅ pytest suite (20 tests) + GitHub Actions CI (`tests/`, `.github/workflows/ci.yml`).

**Tier 2 — 🟡 IMPLEMENTED, dry-run only (2026-06-09; tools/DBs not installed on this machine):**
- 🟡 functional + AMR — `modules.functional`: HUMAnN pathways (read-based) + AMRFinderPlus &
  ABRicate (contig AMR/virulence) + Bakta & eggNOG-mapper (MAG annotation). Data-driven:
  AMR runs with `assembly`, annotation with `binning` (`rules/functional.smk`).
- 🟡 binning refinement — `modules.bin_refinement`: MaxBin2 + CONCOCT → DAS_Tool consensus →
  dRep cross-sample dereplication (`rules/binning_refine.smk`).
- 🟡 metaSPAdes + hybrid assembly — `assembly.assembler: megahit|metaspades`; hybrid via a
  per-sample `long_reads` column (`rules/assembly.smk`, `common.smk`).
- 🟡 second classifier consensus — `modules.classify_consensus`: MetaPhlAn or Kaiju vs kraken2,
  species concordance JSON (`rules/consensus.smk`, `scripts/classifier_consensus.py`).
- 🟡 MultiQC aggregate + Krona viz — `modules.aggregate` (`rules/aggregate.smk`). Composition
  + PCoA plots already in `stats`.

All wired through the same registry+rule+env+routing+preset+docs pattern, validated with
`snakemake -n` (full 52-job DAG builds across both metaSPAdes and Kaiju branches) and a
pytest suite (12 added tests). Real execution needs the external tools/DBs (conda envs under
`workflow/envs/`, `--use-conda`). New preset: `amr-surveillance`.

**Tier 3 — 🟡 IMPLEMENTED, dry-run only (2026-06-09; tools/DBs not installed on this machine):**
- 🟡 aDNA branch — `library: ancient` routing (short-read PE; fastp `--merge` collapses
  overlapping fragments → single reads) + `modules.damage`: mapDamage2 5' C→T / 3' G→A
  deamination → `authentication.json` verdict (`rules/damage.smk`,
  `scripts/damage_authenticate.py`). Preset: `ancient-dna`.
- 🟡 decontam/controls — `modules.decontam` + per-sample `control: true`: prevalence test over
  the Bracken table flags & removes reagent contaminants (`rules/decontam.smk`,
  `scripts/decontam.py`, pure-Python).
- 🟡 strain-level — `modules.strain`: inStrain SNV/microdiversity over reads-vs-contigs
  (`rules/strain.smk`).
- 🟡 per-sample Bracken length — sample-sheet `bracken_read_length` > per-platform map >
  global, wired in `common.smk` + `abundance.smk`.
- 🟡 containers/conda-lock — `Dockerfile` (core image) + `containers/README.md` (conda-lock
  pins + Apptainer/Singularity). Per-rule `container:` directives deferred per request.

Validated with `snakemake -n` (28-job DAG across ancient/damage/strain/decontam; mapDamage +
damage only target ancient samples) and pytest (41 tests total, all green). Registries now 33;
presets 6.

**Gap-closing modules — (2026-06-10), from the independent audit (ASSESSMENT-2026-06-10.md):**
- ✅ **differential abundance** — `modules.differential`: CLR transform + two-sided permutation
  test + Benjamini-Hochberg FDR over the Bracken table between two sample `group` labels
  (an ALDEx2-lite; pure-numpy, `metagx/differential.py` + `rules/differential.smk`). Adds a
  `group` sample-sheet column. **RUN** on 4 subsampled pseudo-samples (correct FDR control: 0
  false positives on null data). Closes the last CRITICAL gap (comparative significance).
- 🟡 **antiSMASH BGC mining** — `modules.bgc` (needs assembly, WGS-only): biosynthetic gene
  clusters on contigs (`antismash.yaml` registry, `rules/bgc.smk`, env `bgc.yaml`, `db.antismash`).
- 🟡 **DADA2 ASV amplicon** — `amplicon.method: asv` denoises short-read marker-gene reads into
  exact sequence variants instead of VSEARCH 97% OTUs (`dada2.yaml` registry,
  `scripts/dada2_asv.R`, env `dada2.yaml`). OTU path unchanged; long-read still Emu.
- **Deduplication**: kept as the opt-in fastp `--dedup` flag (`fastp: {dedup: true}`), not a
  separate module — for shotgun metagenomics deliberate dedup is usually discouraged (duplicates
  can be biological), so flag-level granularity is the correct scope.
- ✅ **Full IMRaD paper generation** (K-Dense-style) — `metagx paper` / MCP `generate_paper`
  (`metagx/paper.py`): elaborates the interview design + registry methods + run results into a
  structured Introduction/Methods/Results/Discussion manuscript (abstract, tables, figures,
  references) and compiles it to **PDF via pdflatex**. **RUN** (real PDF) on the diff-demo.

Registries now 35; pytest 89. The architecture is now broad enough for study-level metagenomics;
what remains is real execution against the external tools/DBs (conda envs under `workflow/envs/`,
`--use-conda`) and subfields explicitly out of scope (metatranscriptomics, viral-host linkage, etc.).

_Each remaining item is "add a registry + a rule + an env + a routing flag + docs" — the
established pattern. They are deferred because they need external tools/DBs that can't be
installed/run in the current environment, not because of design gaps._
