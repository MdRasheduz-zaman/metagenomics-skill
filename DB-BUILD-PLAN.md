# Plan: first-class database building as a pipeline step

Status: **implemented** (2026-06-22). Turns DB construction from a detached CLI pre-step
(`metagx build-db` / `fetch-db`) into a registry-driven, interview-asked, Snakemake-wired
**`db.build`** step — while reusing the existing `dbbuild.py` synthetic path.

Done: registries (kraken2-build/bracken-build), `db.build` config schema + invariants,
`dbbuild.build_database` (4 strategies), `rules/dbbuild.smk` + `scripts/build_db.py` (idempotent
manifest), Bracken read-length coupling, `--use-ftp` default (NCBI rsync is dead) + prefer-prebuilt
doctor warning, air-gap advisory, `metagx build-db --strategy` CLI, `db.build.auto:false` escape
hatch, high-mem `build_kraken2_db` SLURM resources, report DB-provenance manifest, real-taxonomy
header→taxid precheck, and a network-gated viral fetch+classify e2e. Validated end-to-end: a real
`metagx run` built a viral DB and classified noisy ONT reads (29/30 species, real ICTV names).

## Goal

A user with no DB can answer a few interview questions and get the *right* kraken2 +
Bracken (+ optional Kaiju/CAT) database built as part of the run, with the read-lengths,
taxonomy regime, and resources their downstream config actually needs.

## Strategy matrix (what the interview resolves)

| Strategy | Source | Taxonomy | Status |
|----------|--------|----------|--------|
| `have-it` | user-supplied path | n/a | exists (just validate via `doctor`) |
| `prebuilt` | `fetch-db` standard index | real | exists |
| `standard` | NCBI libraries `{bacteria,viral,archaea,fungi,protozoa,human,plasmid,UniVec_Core,nt}` | real | **new** |
| `custom-fasta` | one multifasta | synthetic **or** real | partial (synthetic only today) |
| `custom-folder` | folder of per-genome FASTAs | synthetic **or** real | **new** |
| `spike-in` | custom genomes **+** a standard library | real (custom genomes need real taxids) | **new** |

### Taxonomy fork (asked only for custom/spike-in)
- **synthetic** — flat species-under-root, fabricated taxids, zero NCBI dependency.
  Good for closed-set presence/abundance + simulations. *Loses genus/family rollups.*
- **real** — NCBI `download-taxonomy` + header→taxid resolution (accession2taxid or
  `kraken:taxid|<real>`). Needed for rank rollups, real names, or spike-in mixing.

## Design considerations baked in

1. **Read-lengths come from the sample sheet.** The build rule derives the set of Bracken
   `-l` lengths from the platforms present (ONT→1000, Illumina→150, …) so every
   `databaseLmers.kmer_distrib` the run needs exists. (This is the core reason it must be a
   rule, not a detached call.)
2. **RAM, not disk, is the constraint.** Separate high-mem/long-walltime resources for the
   build job (see scheduler profile below); offer `--max-db-size` cap and `--memory-mapping`
   for classify; `--clean` after build to drop the ~2–3× intermediate library.
3. **Downloads run wherever the rule runs**, honoring `http_proxy/https_proxy`. An **opt-in**
   `db.build.download_on: login` split is offered for the air-gapped minority — NOT default.
4. **Build once, idempotent.** `db.build` is an opt-in target, skipped when a valid
   `.metagx_db.json` manifest + `hash.k2d` exist. Atomic temp-dir + move (mirror `flye`/
   `dbbuild`). Never a hard `input:` that rebuilds on config churn.
5. **Masking dep.** Standard/real builds either declare BLAST+ (`dustmasker`) in the env or
   default `--no-masking` with a noted small false-positive cost (custom path already does).
6. **Scope honesty.** Unify only the kraken2 family (kraken2/bracken/kaiju/CAT, shared
   taxonomy). Domain-module DBs (genomad, checkv, gtdbtk, checkm2, antismash, humann,
   amrfinder, bakta, eggnog) stay as thin per-tool `fetch` wrappers — each has a bespoke
   downloader.
7. **Thread clamp** reused from `dbbuild.py` (online-CPU cap; the kmer2read_distr saga).
8. **Provenance.** `.metagx_db.json` manifest (strategy, libraries, NCBI snapshot date,
   kraken2 version, k/l/s, read-lengths) consumed by `metagx report`/`paper`.

## Work breakdown

### A. Registries (source of truth)
- `metagx/parameters/kraken2-build.yaml` — `--kmer-len`, `--minimizer-len`,
  `--minimizer-spaces`, `--max-db-size`, `--no-masking`, `--load-factor`, library list
  (interpreted), managed `--threads`/`--db`. `interview_spec` carries the strategy +
  taxonomy questions with `promote_when` on goal=simulation / has-genomes-of-interest.
- `metagx/parameters/bracken-build.yaml` — `-k` (must equal kraken2 `--kmer-len`; enforced
  in `validate`), `-l` (managed, from sample sheet), managed `-d`/`-t`.

### B. Config schema
- `db.build:` block: `{strategy, source, libraries[], taxonomy: synthetic|real,
  read_lengths|auto, kraken2_build:{}, bracken_build:{}, download_on: rule|login,
  also:[kaiju,cat]}`. `config_builder.build_config` validates it; `db.kraken2` may be an
  *output* path the build writes to.

### C. `dbbuild.py` extensions
- `--download-library/--download-taxonomy/--build` standard path (real taxonomy).
- folder-of-FASTAs input (iterate `--add-to-library`).
- real-taxonomy option for custom (header→taxid resolution + accession2taxid hook).
- spike-in: add custom genomes to a downloaded standard library before `--build`.
- keep the existing synthetic path as the `taxonomy: synthetic` branch.

### D. Workflow rule
- `workflow/rules/dbbuild.smk`: `rule build_kraken2_db` (idempotent, manifest sentinel,
  atomic move) + `rule build_bracken_db` per derived read-length. `Snakefile` includes it
  when `db.build` is set. classify/abundance gain an order-only dep on the DB manifest.

### E. Runner / CLI
- `metagx build-db` grows `--strategy/--libraries/--taxonomy/--from-folder/--spike-into`.
- `metagx run`: if `db.build` present and DB absent, build first (or require `--build-db`).
- `metagx doctor`: fail-fast validation of `db.build` inputs before a multi-hour burn.

### F. Scheduler resources
- Per-rule resource override (`mem_mb`, long `runtime`, optional `partition`) for
  `build_kraken2_db` in the bundled profiles (`workflow/profiles/*`), since the build job's
  envelope differs sharply from classify.

### G. Advisor/evidence
- Lower default `confidence` for big `standard` DBs vs small custom ones (DB-completeness
  false-positive note already in SKILL.md).

### H. Tests
- registry well-formedness; `validate` enforces `-k == --kmer-len`; read-length derivation
  from a synthesized sample sheet; dry-run DAG gate includes `db.build`; a small real build
  from the committed viral fixture in `test_pipeline_e2e` (synthetic + folder branches).

### I. Docs
- SKILL.md step 2 rewritten around the strategy tree; README install/DB section; the
  handoff `run.sh` gets the DB step as a real (uncommented) stage + the `doctor` gate.

## Resolved decisions (+ guardrails they require)
1. **Default taxonomy = real NCBI** when the user is unsure. Implication: real raises the
   failure surface (network + headers must resolve to taxids). Guardrail: `doctor`/`validate`
   pre-checks header→taxid resolvability **before** the build and, if unresolvable, fails fast
   with a clear message offering `taxonomy: synthetic` or a supplied accession2taxid map.
   Synthetic remains a first-class opt-in, not a hidden fallback.
2. **`metagx run` auto-builds a missing DB** when `db.build` is configured. Implication: a
   plain run can launch a multi-hour/multi-GB job. Guardrails: (a) preflight prints a
   size/time **estimate** from the chosen libraries; (b) idempotent — skipped when a valid
   manifest exists; (c) escape hatch `db.build.auto: false` (then require `--build-db`) for
   batch users who don't want surprise jobs; (d) the build is its own high-mem rule, so it
   slots into the scheduler cleanly rather than blocking the head process.
3. **Kaiju/CAT built only when the consensus/reconcile module is enabled** — kraken2+bracken
   always; the protein/contig DBs follow the module toggle, from the same taxids.
