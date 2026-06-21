# metagx — co-scientist direction & architecture critique

> Consolidated record of two evaluation sessions (2026-06-11). Part I is the architecture
> critique (breadth-vs-depth, agent-above-engine). Part II answers the follow-up question:
> *do we need an agent framework (LangChain/LangGraph) or a mathematical technique (MCMC/Bayes)
> to add the "intelligence" Part I called for?* Nothing here is lost — it supersedes nothing in
> CRITIQUE.md / ASSESSMENT-2026-06-10.md / ROADMAP.md; it sits on top of them.

---

## The one-sentence thesis

metagx is not a competitor to a K-dense–style autonomous co-scientist; it is the **reproducible
execution substrate that such an agent should stand on.** Every design decision should be read
through that lens: keep the engine deterministic and auditable, and put *judgment* in the agent
layer above it — not baked into static files, not delegated to a framework, and not faked with
heavyweight math where domain reasoning belongs.

---

# Part I — Architecture critique (session 1)

## Verdict

A genuinely well-built **execution substrate** that has been optimizing the wrong axis: it chases
**breadth** (38 tool registries) while the property that actually makes it trustworthy — **depth
of real execution** — lags. The honest one-line scorecard: **broad design, narrowly verified.**

## Strong bones — do not redesign

- **Registry-as-truth** + the `inspect.signature` generic-forwarding fix. One YAML drives the
  interview, validation, MCP/HTTP schemas, CLI, and the actual command tokens. The forwarding fix
  removed an entire *class* of drift bug rather than three instances of it.
- **Reconcile / CAT cross-check** (read↔contig concordance) — original and scientifically real;
  separates "seen in both" from read-only artifacts.
- **Provenance manifest with numbers read back from result files** — exactly right.
- **Pure-Python stats (no scipy/R)** — defensible for dependency-lightness and *actually unit-tested*.
- **One-report-per-project** editorial stance.
- Tests currently: **144 pass / 2 skip.**

## The core weakness

The project encodes **scientific judgment** in static files (`evidence/*.yaml`, `recommend` /
`warn_if`, hand-written rule logic) — but judgment does not fit there:

- The advisor is a **lookup** against grids benchmarked on **simulated reads** (wgsim, nanopore-sim)
  against a perfectly-matched DB. That is the *easy* case; real samples with DB gaps behave
  differently, and a static table cannot reason about *why* one sample is an outlier.
- `advise → next_config.suggested.yaml` is **if/else** ("Good's coverage < 0.95 → sequence deeper").
  A thermostat, not a scientist.
- The interview asks **pre-registered** questions; it never forms a hypothesis about the data.

Validation honesty: **~20 modules are dry-run-only.** Dry-run proves the DAG wires and the command
string renders; it does *not* prove the flag is valid for the installed tool version, the parser
survives real output, or the conda env solves. Every "verified by real execution" entry in
ROADMAP.md found a real bug dry-run could not (newick leaf miscount, `KeyError: qc_key` on the
common TSV path, the Bracken `echo`). That is empirical proof that dry-run hides bugs at a high rate.

There are really **two sources of truth**, and CLAUDE.md slightly oversells the single-source claim:
registries for *leaf flags*, and the evidence/routing layer for *cross-tool, data-dependent* decisions
(assembler constrains binner; DADA2 truncation depends on a runtime quality profile).

## The reframe and the seam

Right division of labor:

- **Judgment → the agent** (Claude + skills): per-sample, hypothesis-driven, improvised.
- **Execution → metagx**: every agent decision becomes a validated config + a logged Snakemake run —
  reproducible, provenance-tracked, auditable.
- **The seam is the `advise → suggest → history.jsonl` loop**, which already exists but is wired as
  if/else. Promote it from rule-based to agent-driven: **propose → run → observe → revise.**

## Part I priorities

1. **Freeze breadth; verify depth.** Run one minimal real input through every dry-run-only branch.
   Highest-value step; the value of a pipeline is in the *joints* (parsers, format handoffs, env
   solves), which dry-run skips.
2. **Turn static priors into posteriors** updated from the user's own `history.jsonl` sweeps.
3. **Gate the auto-generated Discussion.** Read-back prevents fake *numbers*, not fake
   *interpretation*. Default `metagx paper` to a Methods+Results skeleton; gate narrative prose.

---

# Part II — Frameworks vs. math (session 2)

The follow-up question: *to add the intelligence Part I called for, do we reach for an agent
framework (LangChain / LangGraph) or a mathematical technique (MCMC / Bayes)?*

**You are half right and half pointed at the wrong layer.** The confusion is that "intelligence"
means two different things here, and they want opposite tools:

| Need | What it is | Right tool | Wrong tool |
|---|---|---|---|
| **Orchestration intelligence** | deciding *what to do next* | LLM judgment (Claude) + the validated config action space | LangChain, RL, a learned policy |
| **Inferential intelligence** | quantifying *uncertainty in the numbers* | lightweight Monte-Carlo / empirical-Bayes, pure-numpy, **inside the stats modules** | full MCMC engines (PyMC/Stan), Bayesian optimization of cheap grids |

Keep those two columns straight and the whole question resolves cleanly.

## 2.1 — LangChain / LangGraph: mostly no, and here is why

**LangChain: no.** It is a heavy dependency with leaky abstractions and high churn, and it is
philosophically the opposite of this codebase's entire ethos ("`registry.py` depends only on PyYAML;
no hidden magic; pure-Python"). Adopting it would *regress* the property that makes metagx good.

**LangGraph: only if you go headless, and even then probably not.** LangGraph is for cyclic,
checkpointed, multi-agent state graphs across processes. The loop you actually need —
interview → config → `metagx run` → read `history.jsonl` → revise — is **a 4-node loop with one
back-edge.** You do not need a graph library to run that; you need a `while` loop and a critic.

The deeper reason the framework question is a red herring: **the hard part of agent design is
defining a safe, typed, validated action space, and you already built it** — the schema-checked
`config.yaml`. That is exactly what LangChain/LangGraph *don't* give you; they give orchestration
plumbing on top of an action space you'd still have to design. You're being offered scaffolding for
the part that's already done, in exchange for the dependency discipline that's your main asset.

**Who is the orchestrator, then?** Today: **Claude Code + `SKILL.md`** already *is* the agent
runtime — it has the conversational context, the tools, and the reasoning. Adding LangGraph means
running a *second*, context-blind agent runtime beside it. Don't.

**When a framework finally earns its place:** only when you want metagx to run the propose→run→
observe→revise loop **autonomously, headless, with no Claude session attached** (e.g. overnight on a
cluster). At that point you need *an* agent runtime — and I would reach for the **Claude Agent SDK**
(purpose-built to keep Claude as the reasoning engine, far less abstraction tax) or a ~100-line
hand-rolled state machine, **before** LangGraph. The principle holds either way: the framework is
plumbing for a loop; it never improves the *quality of the judgment* flowing through the loop, which
is the actual weakness.

> **Bottom line:** orchestration is not your weakness. No framework makes a better scientist; it
> only reschedules one. Spend the complexity budget on the judgment signal, not the pipe it flows in.

## 2.2 — MCMC / Bayes: yes, but inside the statistics, and as Monte-Carlo not MCMC

This is where your instinct is genuinely onto something — you've just aimed it at the agent brain
when it belongs in the leaf computations. Metagenomic counts are **compositional, sparse,
overdispersed, and low-replicate** — precisely the regime where frequentist point estimates are
weakest and Bayesian uncertainty propagation pays off. Ranked by value-per-effort:

### (A) HIGH VALUE — Make "ALDEx2-lite" actually ALDEx2 (Dirichlet Monte-Carlo) ✅ DONE 2026-06-11

> **Implemented.** `metagx/differential.py` now draws `mc_samples` (default 128) Monte-Carlo
> instances from each sample's `Dir(counts+0.5)` posterior, CLR-transforms each, runs the
> permutation test per instance, and BH-corrects the **expected** p-value (effect = median across
> instances). `mc_samples: 1` recovers the legacy point estimate. Plumbed through `config_builder`
> (`differential.mc_samples`) and the workflow script; summary JSON records `method`/`mc_samples`.
> Tests added (Dirichlet shape/centering, MC planted-signal, MC null FDR control, mc=1 fallback
> reproduces the legacy numbers); suite 144→148. **Verified on real data**: the bundled diff-demo
> Bracken table run through the production script at mc=128 gives 0 significant taxa on null
> subsamples (correct FDR control) with valid TSV/JSON/volcano outputs.

Original rationale (kept for the record): `metagx/differential.py` was a permutation test on a
CLR with a **fixed pseudocount of 0.5**
(`diversity.clr`). But the thing that makes real ALDEx2 *Bayesian* is exactly what's missing: ALDEx2
draws Monte-Carlo instances from the **Dirichlet posterior** `Dir(counts + 0.5)` per sample, CLR-
transforms each draw, and propagates that sampling-depth uncertainty into the test. Your version
collapses that posterior to a single point estimate — so it inherits the name without the mechanism.

The fix is **a few lines of pure numpy** (`np.random.dirichlet`), no new dependency, fully testable,
perfectly on-ethos: for each sample draw _k_ Monte-Carlo instances from the Dirichlet posterior,
run the existing permutation test per instance, and aggregate (expected p, expected effect). This
turns a point estimate into honest uncertainty propagation and **directly upgrades the headline
comparative result** the whole tool exists to produce. This is the single most defensible "insert
Bayes" move in the codebase, and it's small.

### (B) HIGH VALUE — Static evidence priors → posteriors from `history.jsonl`

This is Part I's recommendation #2, stated as the math it actually is. The per-platform confidence
recommendations are point estimates from *simulated* data. Model `percent_classified` vs `confidence`
as a **Beta-Binomial** (conjugate, closed-form) or a tiny **Gaussian Process** over the curve
(closed-form, numpy), seed it with the simulation grid as the **prior**, and update with each real
run logged in `history.jsonl`. Fixes the "benchmarked on sim data" weakness, and is the smallest
real instance of "agent above a deterministic engine." **No MCMC needed — conjugate updates.**

### (C) MEDIUM VALUE — Empirical-Bayes shrinkage for low-replicate differential abundance

The documented 2-vs-2 problem (permutation p-value floors at ~0.33) is the classic low-replicate
variance-estimation failure. The principled fix is a **hierarchical model that borrows strength
across taxa** — but full hierarchical inference normally means PyMC/Stan (heavy, off-ethos). The
**empirical-Bayes shrinkage** version (moderate per-taxon variance toward a global prior, à la
limma/voom's moderated _t_) is **closed-form and pure-numpy**, and recovers most of the benefit. It
degrades gracefully at low _n_ where the permutation test simply dies. Sweet spot.

### (D) LOW VALUE / SKIP — Bayesian optimization of the confidence sweep

Tempting ("active-learning which threshold to try next"), but the sweep is a **cheap, ~5-point, 1-D,
near-monotone** search. Just compute the grid exhaustively. BO earns its keep in *expensive,
high-dimensional* spaces; applying it here is sophistication that doesn't pay. Revisit only if sweeps
become joint/multi-param and expensive.

### (E) DON'T — RL / bandits / a learned policy for the agent loop

The propose→run→observe→revise loop is formally a sequential decision problem, so RL/POMDP framing is
*available*. Don't. The action space is huge, the reward ("good analysis") is sparse and ill-defined,
and you get a handful of runs per project — RL needs thousands of episodes. **An LLM brings the
domain prior that an RL agent would need millions of samples to learn.** The orchestration
intelligence must be LLM judgment, not a learned policy. (This is the mirror image of 2.1's verdict:
math belongs in the leaves, reasoning belongs in the loop.)

### (F) DON'T (yet) — full MCMC / PyMC / Stan

Everything valuable above needs **Monte-Carlo the *technique*** (Dirichlet draws, conjugate posterior
sampling) — all pure-numpy. It does **not** need **MCMC the *heavy machinery*** (Metropolis/HMC,
PyMC, Stan), which would import a large dependency tree and break the dependency-light discipline
that is the project's main asset. If a genuinely non-conjugate hierarchical model is ever required,
isolate it behind an optional extra (`pip install metagx[bayes]`) and a per-rule conda env, never in
`registry.py`'s dependency path.

### Honorable mention — Conformal prediction for calibrated call confidence

A lighter, distribution-free alternative to full Bayes for "how much do I trust this taxon call?":
**conformal prediction** gives calibrated intervals with almost no assumptions and a small numpy
footprint. Worth a look as an uncertainty layer on classification/abundance if (A)/(C) feel too
method-specific.

## 2.3 — The principle to carry forward

> **Math earns its place when it quantifies uncertainty in a *number you're already computing*
> (differential abundance, a sweep curve, a diversity estimate). It does not earn its place as the
> "brain" that decides what to do next — that is LLM judgment. Frameworks never add intelligence;
> they only reorganize where it runs. Reach for the lightest tool that closes the specific gap:
> Dirichlet draws before MCMC, conjugate updates before PyMC, a `while` loop before LangGraph, the
> Claude Agent SDK before LangChain.**

---

# Consolidated priority list (both sessions)

| # | Action | Source | Cost | Value |
|---|---|---|---|---|
| 1 | **Freeze breadth; real-execute every dry-run-only module once** — IN PROGRESS 2026-06-11: decontam, Kaiju-consensus, reconcile+CAT, read↔contig accuracy (93.3% concordant), damage-auth, aggregate (MultiQC run; Krona reworked for custom DBs). **2 real bugs fixed** (kaiju parser; Krona custom-DB taxonomy), 7 modules no-tests→tested. Suite 144→166. See ROADMAP run 4 | I | med | ★★★★★ |
| 2 | ✅ **Dirichlet Monte-Carlo in `differential.py`** (make ALDEx2-lite real) — DONE 2026-06-11 | II-A | low | ★★★★★ |
| 3 | **Evidence priors → posteriors from `history.jsonl`** (Beta-Binomial / GP) | I-2 / II-B | low–med | ★★★★ |
| 4 | **Agent-driven `advise→suggest` loop** (replace if/else with LLM critic) | I / seam | med | ★★★★ |
| 5 | **Empirical-Bayes shrinkage** for low-replicate differential abundance | II-C | med | ★★★ |
| 6 | **Gate auto-generated paper Discussion** behind a flag + caveat banner | I-3 | low | ★★★ |
| 7 | Correct CLAUDE.md: two sources of truth (flags vs routing/judgment) | I | low | ★★ |
| — | Conformal-prediction confidence layer | II | med | ★★ (optional) |
| ✗ | LangChain / LangGraph / RL / PyMC-MCMC / BO-of-sweep | II | — | avoid (see text) |

**If you do exactly two things:** (1) real-execute the dry-run modules, and (2) the Dirichlet
Monte-Carlo. The first makes the tool *honest*; the second makes its headline result *correct*. Both
are pure-Python and fit the existing ethos. Everything else is the agent layer, which is where the
K-dense partnership actually lives — and which a `while` loop with Claude as the critic, not a
framework, is enough to build.
