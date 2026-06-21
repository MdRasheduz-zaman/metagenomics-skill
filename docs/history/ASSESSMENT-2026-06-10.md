# metagx — independent audit (2026-06-10)

> **Follow-up (same day): the audit's gaps were then implemented.** Added (registry+rule+env+
> routing+report+tests+docs, the standard pattern): **`modules.differential`** — differential
> abundance between sample `group`s (CLR + permutation test + BH FDR, pure-numpy
> `metagx/differential.py`); **`modules.bgc`** — antiSMASH BGC mining on contigs;
> **`amplicon.method: asv`** — DADA2 ASVs as an alternative to VSEARCH OTUs. New sample-sheet
> `group` column. Registries 33→35; pytest 44→89. **Differential abundance was RUN for real**
> (not just dry-run) on subsampled bundled data — see `config/diff-demo.yaml`
> (`metagx run --config config/diff-demo.yaml`) and `results/diff_demo/`; BGC and ASV are
> dry-run-validated (need antiSMASH / R-DADA2 + DBs). Dataset sources for every gap (incl. the
> i2B-Bergen aDNA reproduction target and amplicon mock communities) are in `DATASETS.md`.
> Dedup was intentionally left as the fastp `--dedup` flag (shotgun dedup is usually discouraged).
> The remaining items below (P2 especially) still stand.
>
> **Also added (K-Dense-style): full IMRaD paper generation** — `metagx paper` /
> MCP `generate_paper` (`metagx/paper.py`) elaborates the interview design + registry methods +
> run results into a structured Introduction/Methods/Results/Discussion manuscript (abstract,
> tables, figures, references) and compiles it to **PDF with pdflatex**. RUN for real on the
> diff-demo (`results/diff_demo/report/paper.pdf`). pytest 89→93.


Goal of this pass: *is the interview-based, BYOK metagenomics tool fine, does it cover
the breadth of metagenomics, and is it flexible enough?* Below is what I **verified by
running**, the **bugs I found and fixed**, and a prioritized list of **what's still worth
doing**. Status tags: ✅ verified · 🟡 implemented but never executed (dry-run only) · ⬜ absent.

## Verdict

The architecture and DX are excellent and the breadth is genuinely wide — it behaves like a
study-level platform, not a toy. The registry-as-single-source-of-truth principle largely
holds. **But two things qualify "fine":** (1) ~20 of the modules (all of Tier 2/3) have only
ever been *dry-run*, never executed against real tools/DBs — dry-run proves DAG wiring and
command *rendering*, not correctness; (2) I found three real drift bugs where downstream
surfaces had fallen out of sync with the builder (now fixed). One true scientific gap remains
(**differential abundance + sample grouping**), plus a few smaller coverage gaps.

## What I verified by running (✅)

- **Unit tests: 41 → 44 pass** (added 3 regression tests for the validate bug below).
- **Dry-run DAGs build for every major branch**, with real input files:
  - maximal short-read (megahit) — **75-job** DAG: classify+sweep, abundance, assembly,
    binning, bin_refinement (CONCOCT/DAS_Tool/dRep), reconcile (+read-accuracy/breadth),
    functional (HUMAnN/AMRFinder/ABRicate), MAG annotation (Bakta/eggNOG), consensus (Kaiju),
    strain (inStrain), stats (diversity), aggregate (MultiQC/Krona).
  - metaSPAdes + **hybrid** (short+long), amplicon (VSEARCH short / Emu long), **ancient-DNA**
    (fastp `--merge` → mapDamage), **decontam** (control sample), **3-domain** taxonomy
    (geNomad/CheckV/GTDB-Tk/CheckM2/EukRep/EukCC) — all validate + build.
- **Interview + params work for all 33 registries** (the core interview engine).
- **Real run intact**: the bundled ONT run (`config.yaml`) is complete; `metagx report`
  regenerates `manifest.json` + `methods.md` + `report.md` cleanly.

## Bugs found and fixed this pass

1. **`metagx validate` had drifted out of sync with the builder** (correctness + the project's
   core principle). `cmd_validate` round-trips the YAML through `build_config` "to reuse all
   validation" but forwarded only ~29 of ~50 kwargs — every Tier 2/3 section
   (`assembly`, `metaspades`, `consensus`, `metaphlan`, `kaiju`, `humann`, `amrfinderplus`,
   `abricate`, `bakta`, `eggnog`, `maxbin2`, `concoct`, `das_tool`, `drep`, `multiqc`, `krona`,
   `mapdamage`, `instrain`, `bracken_read_length_by_platform`) was dropped. Effects: a valid
   hybrid `assembler: metaspades` config was **falsely rejected**, and bad params in any of
   those sections **passed validation silently**. `metagx run` does *not* validate, so this was
   the only CLI gate. **Fix:** forward every `build_config` parameter generically via
   `inspect.signature` (`metagx/cli.py`) — cannot drift again. Added `tests/test_cli_validate.py`.
2. **MCP server wouldn't import without FastAPI** — `mcp_server.py` imported `fastapi` at module
   top, so a Claude-Desktop/Cursor user who installs only the MCP stdio surface (the documented
   *primary* surface) couldn't load the server at all. **Fix:** import FastAPI lazily behind
   `_HAVE_FASTAPI`; the HTTP block is built only when it's present.
3. **HTTP `BuildRequest` exposed only ~13 params** — web agents (ChatGPT Actions/Gemini/
   Perplexity) couldn't drive any Tier 2/3 feature. **Fix:** `model_config = {extra: allow}` +
   filter to the builder's signature, so the HTTP surface reaches the full feature set in
   lockstep with `build_config`.

(All three were the same class of bug: a downstream surface re-listing builder params by hand.
The interview, CLI `build-config`, MCP `build_config`, and the Snakefile were already generic.)

## Coverage map

| Aspect | Status |
|---|---|
| QC (short/ONT/PacBio/amplicon), host removal, subsample | ✅ run |
| Read taxonomy: kraken2 + Bracken + **confidence sweep**; 2nd classifier (MetaPhlAn/Kaiju) | ✅ kraken2/Bracken run; 🟡 consensus |
| Assembly: MEGAHIT / metaSPAdes / Flye / **hybrid** | ✅ MEGAHIT+Flye run; 🟡 metaSPAdes/hybrid |
| Binning + refinement: MetaBAT2 / MaxBin2 / CONCOCT → DAS_Tool → dRep | 🟡 |
| Bin/contig QC + taxonomy: CheckM2/CheckV/EukCC, GTDB-Tk/geNomad, CAT cross-check | ✅ CAT run; 🟡 domains |
| Functional / AMR: HUMAnN, Bakta, eggNOG, AMRFinderPlus, ABRicate | 🟡 |
| Reconcile (read↔contig concordance, breadth, read-accuracy) | ✅ run |
| Strain (inStrain), Ancient-DNA damage (mapDamage), Decontam (controls) | 🟡 |
| Diversity: α (Shannon/Simpson/richness/Pielou), β (Bray–Curtis), PCoA, TSS/CLR | ✅ run |
| **Differential abundance / significance across groups** | ⬜ **absent** |
| **Sample metadata / grouping (case-vs-control study design)** | ⬜ **absent** |
| Amplicon: primer trim + **OTU** (VSEARCH) / Emu; **ASV (DADA2)** | 🟡 OTU; ⬜ ASV |
| BGC discovery (antiSMASH); deliberate dedup step | ⬜ (dedup only as a fastp flag) |
| Reporting/provenance, MultiQC, Krona, presets, SLURM, Docker/conda | ✅ report run; 🟡 MultiQC/Krona |
| Metatranscriptomics, viral-host linkage, long-read polishing | ⬜ out of scope (documented) |

## What still needs doing (prioritized)

**P1 — Differential abundance + sample grouping (only true scientific gap). ✅ DONE & RUN.**
Implemented as `modules.differential` (CLR + permutation + BH FDR, pure-numpy) + a `group`
sample column; run on subsampled data (`config/diff-demo.yaml`). The tool can now answer "which
taxa differ between conditions?" — closing the headline comparative-result gap. (One caveat to
teach users: ≥2 samples/group is the *minimum to run*; a 2-vs-2 permutation p-value floors at
~0.33, so real power needs more replicates — use a real case/control cohort, see DATASETS.md.)

**P2 — Smoke-test the dry-run-only modules at least once (biggest validation risk).**
~20 modules (all 🟡 above) have never actually executed; dry-run does not catch wrong flags for
the installed tool version, parsers that choke on real output, or conda envs that don't solve.
On a machine with the tools/DBs (or `--use-conda` on HPC), run a minimal real input through
each Tier 2/3 branch once and capture the output. This is the single highest-value next step.

**P3 — ASV amplicon (DADA2). ✅ DONE (dry-run).** Added `amplicon.method: asv` → DADA2
(`scripts/dada2_asv.R`, env `dada2.yaml`) alongside the unchanged VSEARCH OTU path. Needs
R/Bioconductor to execute (smoke-test when available; see DATASETS.md mock communities).

**P4 — Smaller coverage. ✅ antiSMASH/BGC DONE (dry-run)** as `modules.bgc` (registry+rule+env+
`db.antismash`). **Dedup:** kept as the opt-in fastp `--dedup` flag, not a separate module —
deliberate dedup is usually discouraged for shotgun metagenomics (duplicates can be biological),
so flag granularity is correct.

**Out of scope (document, don't build now):** metatranscriptomics, viral-host linkage,
long-read polishing — correctly deferred.

## Flexibility (BYOK) — strong

Five client surfaces (Claude SKILL, CLI, MCP stdio, HTTP for web agents, paste-in
`INTERVIEW.md`); user supplies their own DB paths; custom kraken2/Bracken/CAT DB builders need
no NCBI download. With the three fixes above, the MCP-only install works and the HTTP surface
reaches the full feature set — so all five surfaces now expose the same registry-driven
capabilities. No flexibility gaps remain beyond the coverage items listed.
