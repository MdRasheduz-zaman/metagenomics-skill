# Learning `metagx` — an outsider's guide

You just landed in this repo and want to understand it well enough to change it with
confidence. This guide is the path. It assumes you can read Python and have seen YAML; it
does **not** assume you know bioinformatics or Snakemake. Work through it top to bottom —
each section builds on the last, and the "Try it" boxes are the point, not decoration.

> TL;DR of the whole codebase: **one YAML file per bioinformatics tool is the single source
> of truth.** That file drives the interview questions, the config validation, the
> command-line that actually runs, the MCP/CLI surfaces, and more. Learn that one idea and
> everything else falls into place.

---

## 0. What this project *is* (the 2-minute version)

`metagx` runs **metagenomics** pipelines — you give it sequencing reads (DNA fragments from
a sample, e.g. gut microbiome or soil) and it tells you which organisms are present, and
optionally assembles genomes, finds antibiotic-resistance genes, builds phylogenies, etc.

But the *interesting* part isn't the biology — it's the **architecture**. Most pipelines
hard-code each tool's flags in a script. `metagx` instead describes every tool declaratively
in a registry, and derives everything from that description. You can contribute meaningfully
without knowing what "Bracken" does, as long as you understand the registry pattern.

Three audiences use it, all through the same core:
- a **human** on the CLI (`metagx ...`),
- an **LLM agent** via MCP (`mcp_server.py`) or the bundled skill (`SKILL.md`),
- the **Snakemake workflow** that does the real work (`workflow/`).

---

## 1. Set up so you can poke at it

You don't need the heavy bioinformatics tools installed to *learn* the architecture — the
Python core (registries, interview, validation, command construction) runs standalone. The
actual pipeline run needs the bio tools, but skip that for now.

```bash
bash setup.sh                  # uv venv + editable install
uv pip install -e ".[test]"    # dev/test deps
pytest -q                      # should be all green — your safety net
```

> **Try it:** run `pytest -q`. Note the number that passes. After every change you make
> while learning, re-run it. A green suite is how this codebase tells you that the five
> things each registry drives are still in sync.

---

## 2. The keystone: read **one** registry slowly

Open `metagx/parameters/kraken2.yaml`. This is *the* file to understand. kraken2 is a tool
that classifies reads to organisms; its registry is the best-documented one (the field
schema is commented at the top).

Each entry under `params:` is one command-line flag, described by fields like:

| field | meaning |
|---|---|
| `flag` | the literal CLI token, e.g. `--confidence` |
| `type` | `bool` / `int` / `float` / `str` / `path` / `enum` (drives coercion + validation) |
| `default`, `min`, `max`, `choices` | validation bounds |
| `tier` | 1 core / 2 common / 3 advanced — how eagerly to ask about it |
| `ask` | should the interview raise it at all? |
| `question` | the natural-language prompt an LLM uses to ask the user |
| `managed` | injected by the workflow (db/threads/io); user may **not** set it |
| `recommend` / `warn_if` | evidence-based guidance (per platform) |

Three special markers are worth learning now because they're the spine of recent work:

- **`managed: true`** — the workflow owns this flag (e.g. `--db`, `--threads`). Never
  interviewed; the validator *rejects* a user who tries to set it.
- **`passthrough: true`** (the `extra_args` param) — a raw escape hatch. Whatever string the
  user puts here is split and appended to the command verbatim. No flag of its own. This is
  how a power user reaches a flag we didn't model.
- **`interpreted: true`** — user-facing and validated, but consumed by a workflow *script*
  rather than emitted as a CLI flag (e.g. a "method" selector a script maps to tool calls).

> **Try it:** in `kraken2.yaml`, find `confidence` (a `float`, tier 1, with `recommend:` per
> platform and `warn_if:` guardrails) and `db` (a `managed` path). Then find `extra_args` at
> the bottom — the passthrough valve. These three params represent the three "buckets" every
> registry sorts its flags into.

---

## 3. The core engine: `metagx/registry.py` (only ~230 lines)

This is the whole brain. It's intentionally dependency-light (only PyYAML) so it imports
cleanly in both the CLI and inside Snakemake. Read it in this order — these are the five
consumers the registry drives:

1. **`load_registry(tool)`** — loads `parameters/<tool>.yaml`.
2. **`interview_spec(tool, max_tier, context)`** — returns the questions to ask. Filters to
   `ask: true` params at/below the tier. **`context` is the goal-conditional promotion**: if
   you pass `{"goal": "strain_resolved"}`, a normally-quiet param whose `promote_when` matches
   gets pulled into the funnel with a `promoted` note explaining why. (More in §6.)
3. **`validate(tool, values)`** — coerces + bounds-checks user values; rejects unknown and
   `managed` keys.
4. **`render_args(tool, values, managed)`** — turns validated values into the actual list of
   CLI tokens. Bool flags emit only when true; `passthrough` values are split and appended
   last; `interpreted` params are skipped.
5. **`when_matches(when, context)`** — the small clause evaluator shared by `promote_when`
   (bare key = equality; `_gte`/`_lte`/`_gt`/`_lt` suffix = numeric compare).

> **Try it (this is the "aha"):** run these and watch one YAML drive three different outputs.
> ```bash
> python3 -c "from metagx import registry; print(registry.render_args('kraken2', {'confidence':0.1}, managed={'db':'DB','threads':8,'paired':True}))"
> python3 -c "from metagx import registry; import json; print(json.dumps(registry.interview_spec('kraken2', max_tier=2), indent=2))"
> python3 -c "from metagx import registry; registry.validate('kraken2', {'confidence': 5})"   # watch it raise: above max 1.0
> ```
> The first builds a command line, the second builds interview questions, the third validates
> — all from `kraken2.yaml`, no other file touched.

**The payoff to internalize:** to add or change a tool flag, you edit *only* the registry
YAML. The Snakefile rules, the interview, the validation, and the MCP schema all update
themselves because they call back into `registry.py`.

---

## 4. Follow one run end-to-end (the data flow)

Trace how reads become results. Don't run it — just read the files in this order and you'll
see the registry pattern reappear at each hop.

```
user answers ─▶ config_builder.build_config() ─▶ config.yaml ─▶ runner.py ─▶ Snakemake
                       (validates via registry)                    (shells out)
                                                                        │
                                              workflow/Snakefile includes rules/*.smk
                                              based on the `modules:` toggles
                                                                        │
                                  rule kraken2 calls registry.render_args("kraken2", ...)
                                              to build the real command line
```

Concretely:

1. **`config_builder.py` → `build_config(...)`** (line ~261) turns interview answers (+ an
   optional preset) into a validated `config.yaml`. Open `config/illumina-sim.yaml` to see
   what one looks like: a `project`, a `samples` sheet, `modules:` toggles, and a per-tool
   block (`kraken2:`, `bracken:`, `megahit:`) holding the user's chosen values.
2. **`runner.py`** (`metagx run`) shells out to Snakemake with that config.
3. **`workflow/Snakefile`** conditionally `include:`s a `rules/*.smk` file *per enabled
   module*. Read its top: `if MODULES.get("classify", True): include: "rules/classify.smk"`.
   This is why there's no per-experiment code generation — the workflow is fixed and
   reviewed; the *config* decides which rules run.
4. **`workflow/rules/classify.smk`** is the concrete example to read. The `_kraken_cmd`
   function pulls the user's `config["kraken2"]`, adds the workflow-`managed` flags
   (`db`, `threads`, `report`, `paired`...), and calls
   `registry.render_args("kraken2", base, managed=managed)`. That's the registry pattern
   closing the loop at runtime.

> **Try it:** in `workflow/rules/classify.smk`, find the line `args =
> registry.render_args("kraken2", base, managed=managed)`. Now you've seen the *same*
> function you ran by hand in §3 being used to build the actual pipeline command. The CLI,
> the tests, and the live workflow all go through this one door.

---

## 5. The map of the package (where things live)

You don't need all of these on day one, but here's the territory so nothing is a surprise:

- `metagx/registry.py` — **the core** (loader + interview + validate + render_args). Start here.
- `metagx/parameters/*.yaml` — **the source of truth**, one per tool (~38 of them).
- `metagx/cli.py` — the `metagx` command; each subcommand is a `cmd_*` function wired in
  `build_parser()`. The most readable index of "what can this thing do".
- `metagx/config_builder.py` — answers → validated `config.yaml`.
- `metagx/runner.py` — runs Snakemake.
- `metagx/scaffold.py` — generates a registry *stub* from a tool's `--help` (see §7).
- `metagx/sync_help.py` — diffs a tool's live `--help` against its registry (drift detection).
- `metagx/schedulers.py` — HPC backends (local/slurm/lsf/sge/pbs) for `run --executor`.
- `metagx/advise.py`, `tool_advisor.py`, `history.py`, `evidence/*.yaml` — the advisor layer
  (pre-run `recommend`, post-run `advise`, trial history). Reads evidence YAML, not codegen.
- `metagx/report.py`, `paper.py` — provenance manifest + Methods + a full IMRaD manuscript.
- pure-Python stats: `diversity.py`, `differential.py`, `formats.py`, `subsample.py`,
  `readfilter.py`, `dbbuild.py` — **no scipy/R** (a deliberate constraint; match it).
- `mcp_server.py` — MCP server (for agents) **and** a FastAPI app, both thin wrappers over
  the same core.
- `workflow/` — the Snakemake side: `Snakefile` + `rules/*.smk` + `scripts/*.py` + `envs/*.yaml`.
- `tests/` — 24 test files; pytest. CI runs `pytest -q`.

A good first read order: `registry.py` → one registry → `classify.smk` → `cli.py` (skim the
`cmd_*` list) → `config/illumina-sim.yaml`.

---

## 6. The "assume less" idea (recent design, worth understanding)

A design principle runs through the newest code: **the builder should assume as little as
possible about what data a user has or what they want to do.** It's implemented as three
buckets per registry plus a promotion mechanism:

- **Capability-complete:** a registry should model *every* flag a tool exposes, so nothing is
  unreachable. Flags split into `managed` (workflow owns it), curated user params (a small
  `ask:true` funnel + the rest `ask:false`/`tier:3` but reachable), and one **`extra_args`
  passthrough** valve for the genuine long tail. Every one of the ~38 tools now has the valve
  (enforced by `test_every_tool_has_extra_args_valve`).
- **Opinionated but minimal defaults:** the interview still asks only a handful of questions.
  The median case is a *default*, not a ceiling.
- **`promote_when` — the funnel widens on evidence.** Quiet params carry rules like:
  ```yaml
  promote_when:
    - when: { goal: strain_resolved }
      to_tier: 1
      reason: Strain-resolved goal — collapsing haplotypes would erase the signal.
  ```
  When `interview_spec(..., context={"goal": "strain_resolved"})` sees a match, that param is
  pulled into the funnel with its reason shown. Flye and kraken2 are the reference tools.

> **Try it:** watch the funnel widen.
> ```bash
> metagx interview flye --tier 2                          # just: meta
> metagx interview flye --tier 2 --goal strain_resolved   # meta + keep_haplotypes (with a reason)
> metagx interview flye --context '{"estimated_bases": 6e10}'  # meta + asm_coverage
> ```

---

## 7. How to make a change (your first contributions)

Pick the smallest real task and let the tests guide you.

**A. Add or tweak a flag on an existing tool.** Edit only the registry YAML
(`metagx/parameters/<tool>.yaml`). Add the param with `flag`, `type`, bounds, `tier`, `ask`,
and a `question`. Run `pytest tests/test_registry.py -q`. Nothing else needs touching — the
interview, validation, and command-line pick it up. (This is the architecture's whole promise;
prove it to yourself once.)

**B. Onboard a new tool.** Run `metagx scaffold <command>` (e.g. `metagx scaffold flye`). It
runs the tool's `--help`, parses the full flag set, and prints a registry *stub* — every flag
`ask:false`/`tier:3` with a guessed type, plus the `extra_args` valve. Then **curate up from
complete**: mark `managed` flags, promote the few worth interviewing to `ask:true` with real
questions, fix guessed types/bounds, add `recommend`/`warn_if`/`promote_when`. Save as
`metagx/parameters/<tool>.yaml`. (The scaffolder needs the tool installed; it makes a registry
*capability-complete*, not *correct* — judgment is still yours.)

**C. Add goal-conditional promotion.** Add a `promote_when:` block to a quiet param (see §6).
The clause semantics are in `registry.when_matches`. Add a test mirroring
`test_promote_when_*` in `tests/test_registry.py`.

The repo convention (in `CLAUDE.md`): keep stats pure-Python (no scipy/R); the interview is a
funnel (preserve tiers); reports cover exactly one project. Read `CLAUDE.md` for the codebase
contract and `SKILL.md` for how the pipeline is *operated*.

---

## 8. Where to read next, in order

1. `CLAUDE.md` — the codebase contract (architecture, conventions). **Read this next.**
2. `SKILL.md` — the operational playbook: the interview funnel, every module, the
   platform→tool dispatch table. Read it as "how an agent runs the pipeline".
3. `metagx/parameters/kraken2.yaml` and `flye.yaml` — the two reference registries.
4. `tests/test_registry.py` — the tests *are* the spec for the core. If you want to know what
   a function guarantees, read its test.
5. `workflow/rules/classify.smk` then `workflow/Snakefile` — the runtime side.

When in doubt, follow a single flag: pick `--confidence`, and find it in the registry, the
interview output, the validator's bounds check, and the rendered kraken2 command. That one
thread touches every layer of the system — once you can trace it, you understand `metagx`.
