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
environment state, and the headline "CI runs the real pipeline" is *still red* — the
first-round "fix" was declared done in this very document without the one piece of evidence
that mattered: a green CI run.

---

## 0b. ROUND 2 (the "fix" did not hold) — re-opened ⛔

**This document said §0 was "FOUND & FIXED ✅" and checklist P0 was `[x]`. CI failed again
with the byte-identical error.** That is the single most important finding of this round and
it is a process failure, not just a code one: *victory was declared in the tracker before the
gate went green.* Everything below is written to not repeat that.

**New evidence (ground truth, this round):**
- Re-ran the exact fixture build locally (`dbbuild.build_db(tests/fixtures/viral/genomes.fasta,
  read_length=[150,1000], threads=4)`): **succeeds** — rc 0, `hash.k2d`/`opts.k2d`/`taxo.k2d`
  all written, "Database construction complete." **No SIGPIPE at all.**
- CI (2-core Linux runner) prints **`build_db: OMP only wants you to use 2 threads`** and then
  `xargs: cat: terminated by signal 13`, and **writes no `*.k2d`** — the build aborts *inside*
  step 3, not after it.

**Why the round-1 fix was structurally wrong:** the round-1 fix "trust the artifacts over the
exit code." But in the *actual* CI failure **there are no artifacts to trust** — the build
aborts before writing them. The recovery branch (`_artifacts_present`) is never taken, so it
falls straight through to `failed_step="build"`. The fix addressed a failure mode (build
succeeds, wrapper exits non-zero) that **was not the one CI was hitting** (build aborts,
writes nothing). It was never reproduced against the red environment — it was reasoned about
from the log and shipped. The log even *showed* "Building database files (step 3)..." with no
"complete" line, i.e. an abort, not a clean-build-with-noisy-exit.

**Real root cause (high confidence):** multithreaded `build_db` on a low-core runner. With
≥4 cores (macOS dev box) the build finishes cleanly; with 2 cores the OMP path caps threads,
the internal `cat | build_db` pipe races, `cat` dies with SIGPIPE, and step 3 aborts with no
database. It is environment- and core-count-specific, which is exactly why it is invisible on
the author's machine and only ever shows up in CI.

**Fix applied this round (`metagx/dbbuild.py`):** when `--build` exits non-zero **and** the
`*.k2d` artifacts are missing, retry the build **once with `--threads 1`** before failing
hard (single-threaded build has no pipe race). Recorded as a `recovered`/`retry_threads1`
note. Regression test added: `test_build_db_recovers_via_single_threaded_retry`
(multithreaded abort → single-threaded retry succeeds). Local suite green (6/6 in
`test_dbbuild.py`).

**This round the failure WAS reproduced locally** (the round-1 sin, not repeated). The bug is
core-count-gated, so it was forced by capping OMP threads:
`OMP_NUM_THREADS=2 metagx … build_db(threads=4)` on macOS reproduces CI **exactly**:

```
build_db: OMP only wants you to use 2 threads
xargs: cat: terminated with signal 13     # SIGPIPE
kraken2-build --build → rc 64, and ZERO *.k2d written   ("no matches found")
```

That is a *faithful* reproduction (same message, same no-artifacts outcome) — and against it:
the round-1 artifact-trust path **does not fire** (correct: there are no artifacts), and the
round-2 single-threaded retry **does** recover (`retry_threads1: 0`, full DB built, `ok:True`).
So the fix is verified against the actual failure mode, not merely reasoned from a log.

> **Honesty caveat (still do not check the box until CI is green):** the reproduction is
> faithful but it is still macOS-with-forced-OMP-cap, not the literal GitHub runner. Confidence
> is now high, not certain. **P0 stays OPEN until an actual green CI e2e run exists** — that is
> the only evidence that round 1 lacked, and its absence is what made the "FIXED ✅" claim
> false. On the next CI run, confirm the `retry_threads1` note appears in the build log (proof
> the retry path is what saved it). If `--threads 1` *still* aborts, escalate: force
> `--threads 1` for the whole small-DB build path unconditionally, or bypass the flaky
> `kraken2-build` bash wrapper with a direct `build_db` call.

---

## 0. The CI e2e failure — round-1 attempt (SUPERSEDED by §0b) ⚠️

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

- [ ] **P0 — RE-OPENED ⛔** CI e2e still red; the artifact-based fix did not address the
      actual failure (multithreaded `build_db` aborts on the 2-core runner, writing no
      artifacts). Round-2 mitigation applied: single-threaded retry on missing-artifact build
      failure (`dbbuild.build_db`) + regression test. **STAYS OPEN UNTIL A GREEN CI E2E RUN
      EXISTS** — the round-1 lesson is that this is unverifiable on macOS and must not be
      checked off from local evidence alone. (see §0b)
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

## 6. ROUND 2 — brutal critique addendum (the meta-problems)

The round-1 critique (§1–§5) was about the *product*. Round 2 is about the *process that keeps
shipping the product red*, because that is now the binding constraint.

### 6.1 The tracker declares victory before the gate is green — the cardinal sin
Every item in §5 is `[x] done`, written in confident past tense ("verified", "enforced",
"unblocks CI"). Yet the headline gate (CI e2e) is **red**, and was red when those boxes were
ticked. "Done" in this repo has meant *"reasoned about and shipped on a green-locally machine,"*
not *"observed working in the target environment."* Round 1 proved that standard fails: §0
called the exact bug "FOUND & FIXED ✅" and CI failed again, identically. **Until proven
otherwise, treat every other `[x]` with the same suspicion** (§6.2).

### 6.2 None of the "end-user hardening" claims have a single green Linux signal
Commit `8779fcd` ("Make metagx end-user ready: doctor, fetch-db, packaging, lockfile, honest
CI") is the foundation of the end-user story — and **CI e2e has never gone green since**, so:
- "verified a clean wheel install outside the repo resolves the Snakefile" (P1 packaging) —
  on which OS? The only automated Linux check is the same job that's failing at the DB-build
  step *before* it ever reaches packaging.
- "CI e2e builds from `environment.yml`, parity enforced" (P1 CI parity) — the env may *parse*
  and the parity *test* may pass, but we have **no evidence the env actually solves and the
  pipeline runs on Linux**, because the run dies at `build_db`.
- doctor / fetch-db / lockfile — all plausibly fine, but each is asserted from local runs.
- **Action:** re-audit every §5 `[x]` against one standard — "is there a green CI job that
  exercises this on Linux?" Demote any that fail that test back to `[ ]` with a note. An
  honest tracker beats a green-looking one.

### 6.3 "Fix by reading the log" is now twice-burned — make the red env reproducible
Round 1 reasoned from the log and shipped; it was wrong. This round it was only caught because
the failure was *forced* locally (`OMP_NUM_THREADS=2`). That trick should be **institutional,
not a one-off**: the entire class of "only fails on the small/low-core CI runner" bugs is
invisible until someone constrains the environment. **Action:** add a `make repro-ci` (or a
pytest marker) that runs the fixture DB build under `OMP_NUM_THREADS=2` (and ideally
`taskset`/cgroup core-pinning), so the low-core path is exercised on *every* dev machine, not
just discovered in CI weeks later.

### 6.4 Correctness depends on parsing a flaky third-party bash wrapper
The small-DB build path's success/failure hinges on interpreting `kraken2-build`'s noisy,
non-deterministic exit behavior (SIGPIPE → 64, sometimes with artifacts, sometimes without).
We are now on our **second** layer of heuristics (artifact-trust, then thread-retry) wrapped
around someone else's bash script. Each heuristic is a guess about a wrapper we don't control.
**Action (longer-term):** for the synthetic small-DB path, call `build_db` directly with known
args instead of going through `kraken2-build`, removing the wrapper and its pipe race entirely.
That converts two heuristics into one deterministic call.

### 6.5 The local prebuilt DB is a structural blind spot, not just one bug
The root cause of *both* CI rounds is the same shape as P0 in §1: locally the prebuilt,
gitignored `viral_custom` DB means **the build code never runs**, so every build-path bug is
shipped blind. The DB-build fix narrows it, but the pattern recurs anywhere a gitignored
artifact short-circuits a code path. **Action:** add a CI-only env flag (e.g.
`METAGX_FORCE_DB_BUILD=1`) the fixture honors to *always* build from the committed fixture even
when a prebuilt DB is present, so the build path has coverage on every CI run by construction.

### 6.7 ROUND 3 — the round-2 fix did not turn CI green, and the log was blind
The single-threaded retry shipped; CI failed again, **identically**. But this round exposed a
second discipline failure: **the e2e log was uninformative.** The fixture asserted
`f"DB build failed: {res}"`, and pytest's repr abbreviates the large `res` dict with `...`,
so `failed_step`, `logs`, the build's stderr, and *whether the retry even ran* were all hidden.
Three rounds of guessing partly because the test threw away the evidence.

**Round-3 actions:**
- **Observability first (`tests/test_pipeline_e2e.py`):** on build failure the fixture now
  raises with the *untruncated* `failed_step` + returncode + tool `tail` + the `retry_threads1`
  result. The next CI log will finally show what `kraken2-build` actually printed and whether
  the retry path executed — so the next fix is evidence-driven, not another guess.
- **Strengthened retry (`metagx/dbbuild.py`):** the single-threaded retry now also forces
  `OMP_NUM_THREADS=1` in the subprocess env (—threads 1 caps build_db's own threads but its
  libgomp regions could still spawn the racing reader).

**Honest status:** this push may *still* be red. That is acceptable this round only because its
purpose is to make the failure legible. **Do not attempt a fourth blind fix — read the new log
first.** If the retry ran and still produced no DB, escalate straight to §6.4 (bypass the
`kraken2-build` wrapper; call `build_db` directly) rather than adding a fourth heuristic.

### 6.8 ROUND 4 — real root cause found (the observability paid off)
The round-3 diagnostic dump immediately earned its keep. The next CI log showed the kraken2
build step **passing** and the failure moving to a *new* step with a *clear* message:

```
failed_step : bracken-build-150
returncode  : 1
kmer2read_distr:   thread count exceeds number of processors
```

So the SIGPIPE saga was real but secondary. The actual blocker on the 2-core runner is
**Bracken**: its `kmer2read_distr` is called with `-t 4` and **hard-aborts** when the thread
count exceeds online CPUs (kraken2 only warns and reduces — which is why kraken2 looked like
the culprit for three rounds). Same environmental constraint (2 cores), stricter tool.

**Fix (`metagx/dbbuild.py`): clamp `threads` to the online CPU count** (`_usable_cpus()`,
affinity-aware) before building any command. This deterministically fixes bracken **and**
removes the kraken2 thread>core mismatch that fed the SIGPIPE race in the first place — one
root cause, one fix, no heuristic. `result["threads"]`/`note_threads` record any clamp.

**Verified against the real constraint (not reasoned):** forcing `_usable_cpus()→2` and running
the **real** kraken2-build + bracken-build end-to-end builds the DB cleanly (`ok:True`,
threads=2, no failed step). Regression: `test_build_db_clamps_threads_to_cpu_count`. Suite 7/7.

**Residual follow-up:** the same thread>core hard-fail can bite an end-user with a small
machine anywhere Bracken/kmer2read_distr runs in the *workflow* (not just `dbbuild`). Audit the
bracken rule's thread param (`render_args`) and apply the same clamp at the workflow level.

### 6.6 Tracked checklist (round 2)
- [ ] **P0** Get CI e2e **actually green** — the only acceptance criterion. (round-2 retry +
      round-3 observability/`OMP_NUM_THREADS=1` applied; still awaiting a green run. **Next step
      is to READ the now-untruncated build log, not to guess again.**)
- [ ] **P1** Re-audit every §5 `[x]` against "is there a green Linux CI job proving it?"; demote
      the unproven ones to `[ ]`. (§6.2)
- [ ] **P1** `make repro-ci` / pytest marker: run the fixture DB build under `OMP_NUM_THREADS=2`
      (+ core-pinning) so low-core failures surface on dev machines. (§6.3)
- [ ] **P2** Replace the `kraken2-build` wrapper call with a direct `build_db` invocation for the
      synthetic small-DB path, eliminating the SIGPIPE race and both heuristic layers. (§6.4)
- [ ] **P2** `METAGX_FORCE_DB_BUILD=1` so the e2e fixture builds from the committed fixture even
      when a prebuilt DB exists — give the build path CI coverage by construction. (§6.5)

---

*Bottom line (round 2):* the architecture is strong and the author knows the domain — that has
not changed. What this round exposes is a **discipline gap, not a knowledge gap**: the project
keeps declaring fixes "done" in its own tracker on the strength of a green-locally machine,
while the one Linux gate that represents an actual end-user stays red. The DB-build bug is now
fixed *and reproduced* — but the rule that should outlive this bug is simple: **nothing is
"done" until the target environment says so. The box stays unchecked until CI is green.**
