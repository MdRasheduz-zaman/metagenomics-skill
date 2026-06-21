# metagx — Critical Evaluation (2026-06-22)

**Reviewer stance:** brutal, evidence-based, end-user-first. The bar is not "does it run on
the author's Mac" — it is "can a stranger who is not the author run a correct analysis on
their own data, on their own machine, without reverse-engineering tribal knowledge."

**Method:** reproduced the test suite locally, reproduced the custom-DB build path, read the
CI config + the actual CI failure log (`e2e_error.txt`), inspected packaging, onboarding,
and the gap between what CI *claims* to test and what it *actually* runs.

**One-line verdict:** the engine is genuinely good and well-tested *where the author's
machine is the environment*. The project's stated goal — "work for end-users, not only
locally" — is **not yet met**: the green path is narrow, depends on undocumented
environment state, and the headline "CI runs the real pipeline" was masking a real,
user-facing bug (now fixed).

---

## 0. The CI e2e failure — ROOT CAUSE FOUND & FIXED ✅

**Symptom (from `e2e_error.txt`):** 7 errors, all in the session-scoped `viral_db` fixture:

```
build_db: OMP only wants you to use 2 threads
xargs: cat: terminated by signal 13      # SIGPIPE
kraken2-build --build ... -> returncode 64
```

**Root cause:** kraken2-build's internal `cat … | build_db` pipe. On a *small* DB,
`build_db` finishes consuming the library and closes the pipe before `cat` is done, so
`cat` dies with SIGPIPE (signal 13), `xargs` reports failure, and the bash wrapper exits
64 — **even though `hash.k2d`/`opts.k2d`/`taxo.k2d` were written correctly.** This is a
well-known kraken2 quirk.

**Why it only failed in CI:** locally `local_databases/viral_custom` is **prebuilt** (and
gitignored), so the fixture returns it and *never runs the build*. CI has no prebuilt DB,
so it builds fresh from the fixture and hits the SIGPIPE. This is the textbook
"works-only-locally" trap: the failing code path had **zero local coverage** because local
state hid it.

**Why this matters beyond CI:** `dbbuild.build_db` is the *same function* behind the
`metagx build-db` end-user command. **Any user building a small custom kraken2 DB on Linux
hits this exact false failure.** Fixing CI and fixing the user are the same fix.

**Fix applied (`metagx/dbbuild.py`):** trust artifacts over exit codes. After a build step
returns non-zero, check whether its expected outputs (`*.k2d` for the kraken2 build,
`database{L}mers.kmer_distrib` for Bracken) exist and are non-empty; if so, record a
`recovered` note and continue. A non-zero exit with *missing* artifacts still fails hard.
Regression tests added (`tests/test_dbbuild.py`):
`test_build_db_recovers_from_kraken2_build_sigpipe` (recovers) and
`test_build_db_real_failure_still_fails` (does not mask genuine failure). Full suite:
**286 passed, 16 skipped.**

> Note: the e2e test itself is a backstop here — it *uses* the DB downstream (asserts
> species recovered, Bracken sums to 1), so a truncated DB would still fail later. The fix
> is safe.

---

## 1. Critical end-user blockers (ranked)

### P0 — "Works only locally" is structural, not incidental
- **The failing CI path had no local coverage** (§0). The general disease: tests and
  daily dev run against gitignored local state (`local_databases/`, `data/`) that the
  author has and nobody else does. The green checkmark on the author's machine is not
  evidence the pipeline works for anyone else. *Mitigation started:* the DB-build fix means
  the fixture path now exercises the real build. But the pattern (local state masking
  bugs) should be hunted elsewhere too.

### P1 — Packaging: `workflow/` and `mcp_server.py` are not in the wheel
- `pyproject.toml` packages only `metagx/` (plus `parameters/`, `presets/`, `evidence/`
  force-included). The **Snakefile, all rules, scripts, and `workflow/envs/` are not
  installed.** `runner.workflow_path()` looks next to the package, then falls back to
  `cwd`. Consequence: a real `pip install metagx` user (non-editable, outside the repo)
  gets `FileNotFoundError: Could not locate workflow/Snakefile`. **The tool only works
  from a git clone, run inside the clone.** That is the definition of "works only locally."
  - *Mitigated by:* the Docker image (`pip install -e /opt/metagx` with the full repo) and
    the "git clone + editable" skill-distribution model. So severity is P1, not P0 — but
    the README/SKILL should be explicit that this is **not** a `pip install`-and-go tool,
    and `mcp_server.py` likewise lives only in the source tree.
  - *Fix options:* (a) move `workflow/` under `metagx/workflow/` and ship it as package
    data, resolving via `importlib.resources`; or (b) document loudly that metagx is a
    repo-checkout tool and make `workflow_path()` raise a *actionable* error naming the fix.

### P1 — Apple-Silicon / macOS path is a minefield (and it's the author's own platform)
  From `CLAUDE.md` + `environment.yml`, the macOS-arm64 story requires the user to know:
  `CONDA_SUBDIR=osx-64` (Rosetta) but **never leak it to base** or you corrupt base;
  Bracken's osx-64 conda build is **broken**; MEGAHIT **segfaults** under Rosetta (tests
  dodge it via provided-contigs); abricate hard-pins **samtools 0.1.x** and must be kept
  out of the core env; `PATH` must be **appended** not prepended. Every one of these is a
  landmine an end-user will step on. This is encoded as tribal knowledge in CLAUDE.md, not
  as guardrails the tool enforces or a preflight that detects-and-explains.
  - *Fix:* a `metagx doctor` preflight that detects arch/env hazards and prints the exact
    remedy; and steer Mac users to Docker/Linux as the supported path.

### P1 — The real-world blocker is databases, and onboarding underplays it
  Real metagenomics needs large reference DBs (kraken2 standard ≈ 50–100 GB+). The custom
  small-DB build path was broken on Linux until §0. There is no first-class, tested,
  documented "get a usable DB" flow for a new user (download standard? build custom? how
  big? where?). For most users this is step 1 and the highest abandonment point.

### P2 — Reproducibility is "floors," not pins
  `environment.yml` and the version-floor test use `>=`. Bioconda drift can silently change
  results between two installs a month apart — unacceptable for a *scientific* pipeline that
  also ships a `paper`/Methods generator claiming provenance. The Dockerfile comment admits
  this ("Pin to … a conda-lock file") but no lockfile is committed.
  - *Fix:* commit a `conda-lock`/explicit-spec lockfile; pin the Docker base by digest.

---

## 2. CI / test-coverage reality check

The "CI is no longer dry — it runs the real pipeline" claim is **half true and oversold**:

- The `e2e` job runs, but most e2e tests **skip in CI** because they need gitignored
  `data/`: illumina, pacbio, differential/diversity, kaiju-consensus, AMR, and **the
  entire aDNA test** (`data/genomes.fasta` absent). What actually executes in CI is a thin
  slice: ONT classify, assembly/bin/reconcile, phylogenetics, aggregate, provided-contigs.
- **`environment.yml` ≠ CI `create-args`.** CI installs a hand-picked subset and omits
  porechop_abi, chopper, seqtk, CAT, prodigal, checkv, genomad, cutadapt, vsearch, seqkit.
  So long-read QC, CAT contig taxonomy, viral domain (geNomad/CheckV), amplicon, and the
  consensus DB build are **never exercised in CI** — and the `--use-conda` heavyweights
  (GTDB-Tk, CheckM2, antiSMASH, DAS_Tool, dRep, inStrain) are tested **nowhere**, locally
  or in CI. "~20 modules" is the marketing number; the verified number is far smaller.
- **Action:** either commit tiny fixtures so the skipped paths actually run in CI, or stop
  advertising coverage the suite doesn't provide. A coverage matrix (module × CI-real /
  local-real / dry-only / untested) belongs in the README so claims are honest.

---

## 3. Documentation sprawl (contributor-facing, but it signals churn)

Root holds 11 markdown files + a PDF: `ASSESSMENT-2026-06-10.md`, `CRITIQUE.md`,
`LEARNING.md`, `ROADMAP.md`, `CO-SCIENTIST-DIRECTION.md`, `DESIGN-multidomain-and-db-scaling.md`,
`DATASETS.md`, `BENCHMARKING-DATASETS.md(+.pdf)`, plus `README/SKILL/CLAUDE` (18 tracked
`.md` total), now joined by `e2e_error.txt`. Several are overlapping, dated self-assessments
that age into noise. New contributors can't tell which is current truth.
- *Fix:* move historical assessments to `docs/history/`; keep one living `ROADMAP.md` and
  this file. Add `e2e_error.txt` to `.gitignore` (or delete after triage).

---

## 4. What's genuinely good (credit where due)

- **The registry-as-single-source-of-truth** design is excellent and the discipline holds:
  one YAML drives interview + validation + MCP schema + CLI + the actual command line. This
  is the project's real moat.
- **`test_workflow_dryrun.py`** (builds the full DAG and resolves every `render_args`
  lambda with no bio tools) is a smart, cheap gate that catches registry/rule drift.
- **The e2e tests, where they run, assert *answers*** (classified-fraction bands, species
  recovered, Bracken sums to 1, FDR yields zero false positives on null data) — not just
  "a file appeared." That is the right standard.
- **Pure-Python stats** (diversity, ALDEx2-lite differential, decontam) keep the core
  dependency-light and testable.
- **Hard-won pitfalls are documented** (samtools 0.1.x pin, `__future__`-in-`script:`,
  iqtree2→iqtree3, minimizer-report column shift). The knowledge exists — it just needs to
  become *enforcement*, not lore.

---

## 5. Tracked action checklist

- [x] **P0** Fix kraken2-build SIGPIPE false-failure in `dbbuild.build_db` (artifact-based
      success) + regression tests. *(done 2026-06-22; unblocks CI e2e and `metagx build-db`)*
- [x] **P0** Audit for other "local-state-masks-bugs" paths; ensure every code path with a
      gitignored-state shortcut also runs from-scratch in CI. *(done: the Kaiju consensus DB
      now builds from the committed fixture (`kaiju_db` session fixture) and the aDNA e2e falls
      back to the committed viral fixture instead of gitignored `data/genomes.fasta` — both now
      run in CI instead of skipping. Added always-on dry-run DAG scenarios for amplicon /
      strain / bgc / stats, whose `render_args` call sites were previously resolved by no test.)*
- [x] **P1** Packaging: ship `workflow/` (+ `mcp_server.py`) as package data *or* document
      loudly that metagx is a repo-checkout tool; make `workflow_path()` error actionable.
      *(done: hatchling force-includes `workflow/`→`metagx/workflow/`, `mcp_server.py`,
      `environment.yml`; `workflow_path()` resolves inside the installed package and raises an
      actionable error; verified a clean wheel install outside the repo resolves the Snakefile;
      `tests/test_packaging.py`.)*
- [x] **P1** `metagx doctor` preflight: detect arch/env hazards (arch mismatch, broken
      Bracken build, samtools downgrade, missing DB) and print exact remedies. *(done:
      `metagx/doctor.py` + `metagx doctor [--config --json --strict]` + MCP `check_environment`;
      `tests/test_doctor.py`. Version floors are now single-sourced in `doctor.VERSION_FLOORS`.)*
- [x] **P1** First-class, tested, documented "get a usable database" onboarding flow. *(done:
      `metagx/dbfetch.py` + `metagx fetch-db [--list]` + MCP `fetch_database`; curated,
      HEAD-verified prebuilt indices (the old Snakefile default URL was a 404 — fixed and routed
      through `dbfetch.index_url`); `tests/test_dbfetch.py` incl. a network-gated live-URL check;
      documented in README + SKILL.)*
- [x] **P1** Make CI `create-args` track `environment.yml` (or commit fixtures) so advertised
      modules are actually exercised; publish an honest coverage matrix. *(done: CI e2e builds
      from `environment.yml` (`environment-file:`), parity enforced by
      `tests/test_ci_env_parity.py`; consensus + aDNA now run in CI via committed fixtures; an
      honest module × verification-level coverage matrix is in the README.)*
- [x] **P2** Commit a conda-lock lockfile + pin Docker base by digest (reproducibility).
      *(done: committed `conda-lock.yml` (linux-64, 411 hash-pinned packages — samtools 1.23,
      no 0.1.x regression); `scripts/lock-env.sh` regenerates it; Dockerfile installs from the
      lock when present and pins `mambaorg/micromamba` by manifest-list digest. osx-64 is
      intentionally unlocked — the Rosetta stack doesn't solve as a unit.)*
- [x] **P2** Consolidate root docs into `docs/history/`; gitignore/delete `e2e_error.txt`.
      *(done: moved ASSESSMENT/CRITIQUE/LEARNING/CO-SCIENTIST-DIRECTION to `docs/history/`
      (with an index README) and the design doc to `docs/`; root now holds only living docs;
      `e2e_error.txt` deleted and already gitignored.)*

---

*Bottom line:* the architecture is strong and the author clearly knows the domain. The gap
is entirely between "runs on my machine, with my prebuilt DBs and my conda env" and "a
stranger can install and run it correctly." The CI failure was the first, concrete proof of
that gap — and it's now fixed at the source, for the user too, not just for the test.
