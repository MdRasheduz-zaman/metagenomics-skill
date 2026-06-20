# Spec: `metagx probe` — measured pre-flight context

**Status:** phase 1 (MVP) shipped · **Scope:** the "local probe" tier only (advisory and
external/closed-loop tiers are out of scope here) · **Owner:** TBD

## 1. Purpose

Replace *asserted* interview context ("I think my data is deep / HQ") with *measured* context
from the user's actual reads, so goal/data-conditional promotion (`promote_when`) and
platform/QC routing fire on facts instead of guesses.

The mechanism already exists: `registry.interview_spec(tool, max_tier, context=...)` consumes a
`context` dict. Today the LLM fills it from the interview. `metagx probe` fills it from a cheap,
local, non-reconstructive subsample of every sample in the sheet. **No new promotion machinery —
just a better source for `context`.**

## 2. Non-goals

- No sending data or summaries off the machine (that is the separate external/closed-loop tier).
- No silent per-dataset auto-tuning of parameter *values* (reproducibility hazard). Probe
  *promotes questions* and *routes platform/QC*; it does not secretly pick numbers. Any value it
  suggests is surfaced + logged, never applied invisibly.
- No dependency on bio tools being on `PATH`. Pure-Python, reusing `formats.py` / `subsample.py`.
- No scipy/R (project convention).

## 3. Consent model (decided)

The trust boundary is **"does anything leave the machine?"** — and the probe never crosses it.
A full `metagx run` already reads 100% of every read file locally; the probe reads a subsample
and stops, emitting only aggregates. So the probe is strictly *within* the footprint a user who
agreed to run the pipeline has already accepted.

Rules:

1. **Opt-in, asked once, remembered.** First use prompts for consent; the choice is stored in
   `.metagx/consent.json` (`{"probe": "local"|"off", "ts": ...}`) so we never re-nag.
2. **Re-confirm on escalation.** Any future capability that would move data or summaries off-box
   must ask again — the stored `local` consent never implies external consent.
3. **Decline = advisory, not failure.** `off`, non-interactive with no stored consent, or unreadable
   files → return an empty/partial context and fall back to today's a-priori behavior. Nothing breaks.
4. **Non-reconstructive by construction.** Output and provenance carry only aggregate statistics —
   never read sequences, never read IDs. This is the technical guarantee behind the promise and is
   enforced by a test (§9).
5. **Transparent + auditable.** The probe report states what was read (file, N reads sampled), what
   was computed, and what was written. Provenance logs *that* a probe ran and its benign summary.

Default-on policy when execution is available vs always-opt-in is a deployment choice (set in
`SKILL.md` guidance / a config flag); the code path is identical either way.

## 4. Scope: all samples, per-sample + reconciled

A metagx project is a **sample sheet** (see `config/*.samples.tsv`: `sample, r1[, r2], platform,
layout`). Samples legitimately differ (platform, depth, host load, quality), and the workflow
already dispatches the assembler per-sample by platform. So the probe is **per-sample**, then
reconciled to a project view:

- **Per-sample profile** for every sample (R1 representative; R2 spot-checked for pairing sanity).
- **Reconciled project summary**: min/median/max of each metric across samples.
- **Heterogeneity + sheet-error warnings**, e.g. "samples 1–3 measure ONT-HQ, sample 4 looks short
  & clean (Illumina?) — check the sheet"; "host fraction ranges 5–60%". This sheet-validation is
  valuable independent of promotion — it catches mislabeled inputs *before* a multi-hour run.
- **Bounded cost**: cap reads/file (`--max-reads`, default 100_000) and total samples scanned
  (`--max-samples`, default all; ceiling for very large sheets). Streaming, single pass, O(sampled
  reads) memory.

## 5. What it measures (per sample)

All derivable from raw reads with no reference except an optional host index:

| metric | how | feeds |
|---|---|---|
| `n_sampled`, `format` (fasta/fastq) | `formats.read_format`, streamed count | provenance, FASTQ-only logic |
| read-length: `min/median/p90/max`, dist | length of each sampled read | platform inference, `estimated_bases` |
| `gc_fraction` | base counts | contamination/oddity flags |
| quality: `mean_q`, `q20_frac`, `est_error` | Phred decode (FASTQ only) | HQ vs raw read-type routing, QC strictness |
| `dup_fraction` | hash of first k bp of sampled reads | library-complexity / over-amplification flag |
| `inferred_platform_class` | length + error heuristic → `illumina`/`ont`/`pacbio_hifi`/`pacbio_clr` | compared against sheet's declared platform |
| `host_fraction` (optional) | minimap2/mash vs a host index **iff** provided + on PATH; else skipped | promote host removal |
| `estimated_bases` | `n_in_file_est × median_len` (file-size extrapolation, not full read) | `--asm-coverage` promotion |

`host_fraction` is the only metric needing an external tool; it degrades to `null` (skipped) when
unavailable — everything else is pure-Python.

## 6. Output: the context contract

Two artifacts:

### 6.1 Human/agent report (`results/<project>/probe/probe.json` + a short markdown)

```json
{
  "consent": "local",
  "samples": {
    "gut1": {"format":"fastq","n_sampled":100000,"read_length":{"median":151,"p90":151,"max":151},
             "mean_q":34.2,"q20_frac":0.98,"est_error":0.004,"gc_fraction":0.43,
             "dup_fraction":0.07,"inferred_platform_class":"illumina","declared_platform":"illumina",
             "host_fraction":null,"estimated_bases":4.5e9},
    "gut2": {"...": "..."}
  },
  "project": {
    "n_samples":2,
    "read_length_median":{"min":151,"median":151,"max":151},
    "platform_consensus":"illumina",
    "warnings":[]
  }
}
```

### 6.2 The `context` dict for `interview_spec(context=)`

A flat, promotion-ready reduction of the project view. Reducer per fact is chosen for safety:

```python
context = {
    "goal": <from interview, untouched>,
    "estimated_bases": max(per-sample estimated_bases),   # ANY sample deep -> surface asm-coverage
    "platform_class": project.platform_consensus,          # None/"mixed" if heterogeneous
    "max_est_error": max(per-sample est_error),
    "max_host_fraction": max(host_fraction or 0),
    "any_sample_low_q": any(q20_frac < 0.9),
    "platform_mismatch": any(inferred != declared),        # routes a warning, not a silent change
    "measured": True,                                       # marks context as data-derived vs asserted
}
```

`promote_when` clauses then match these keys with the existing `when_matches` semantics
(`estimated_bases_gte`, `max_host_fraction_gte`, `platform_class` equality, ...). Per-sample facts
that drive per-sample routing (host removal, read-type flag) are read from `samples` by the config
builder, not flattened into the project context.

## 7. Integration points

- **New module** `metagx/probe.py` (pure-Python). Reuses `formats.read_format`, `formats.is_gzipped`,
  and the streaming readers in `subsample.py` (refactor `_iter_fastq`/`_iter_fasta` into a shared
  internal iterator if needed). No new heavy deps.
- **New consent module** `metagx/consent.py`: `get(scope) -> str|None`, `set(scope, value)`,
  storage in `.metagx/consent.json`. Tiny.
- **CLI** `metagx probe`:
  ```
  metagx probe --samples sheet.tsv [--max-reads 100000] [--max-samples N]
               [--host-index PATH] [--out results/<proj>/probe] [--yes] [--json]
  ```
  `--yes` records local consent non-interactively (for scripted/agent use). Prints the report;
  with `--json` prints the `context` dict for piping into the interview.
- **CLI** `metagx interview <tool> --probe results/<proj>/probe/probe.json` — load measured context
  instead of (or merged over) `--context`/`--goal`.
- **MCP** `run_probe(samples, max_reads, host_index, consent)` returning the report + context; and
  `get_interview(..., context=)` already accepts it (done).
- **Provenance** (`report.py`): record that a probe ran, the consent mode, and the non-reconstructive
  summary, so a Methods section can say "library stats were profiled on N-read subsamples (Table Sx)."
- **config_builder**: optional `probe=<report>` kwarg; uses per-sample facts to set per-sample
  routing defaults (host_removal on where `host_fraction` high, read-type from `inferred_platform_class`)
  — still surfaced to the user, never silently final.

## 8. Control flow

```
metagx probe
  ├─ consent.get("probe")  ──"off"/none(non-interactive)──▶  return {"measured": False}  (advisory)
  ├─ load sample sheet (reuse config_builder._validate_samples shape)
  ├─ for each sample (≤ max_samples):
  │     stream ≤ max_reads from R1 (spot-check R2 pairing)
  │     compute per-sample metrics  (pure-Python; host_fraction iff index+tool)
  ├─ reconcile → project summary + heterogeneity/sheet-error warnings
  ├─ write probe.json + probe.md ; log to provenance
  └─ return context dict  ──▶  interview_spec(context=)  ──▶  promotion fires on measured facts
```

## 9. Testing plan (no bio tools, no scipy)

- `test_probe_metrics_on_synthetic_fastq` — write tiny FASTQ/FASTA fixtures (short clean vs long
  noisy) and assert measured read-length / mean_q / inferred_platform_class.
- `test_probe_subsamples_all_samples` — multi-sample sheet → a profile per sample.
- `test_probe_flags_sheet_mismatch` — a sample declared `ont` but written as short clean reads →
  `platform_mismatch` warning.
- `test_probe_context_drives_promotion` — feed the probe's context into `interview_spec("flye", ...)`
  and assert `asm_coverage` / host-removal promotion fires (closes the loop with §6).
- `test_probe_output_is_non_reconstructive` — **privacy guarantee**: assert no read sequence or read
  ID substring from the input appears anywhere in `probe.json` / context. Enforces §3.4.
- `test_probe_consent_off_degrades` — consent `off` / non-interactive → `{"measured": False}`, no files read.
- `test_probe_bounded` — `--max-reads` honored; a huge file is not fully read (cap on bytes/records).

## 10. Phasing

1. ✅ **MVP (shipped)**: pure-Python metrics (everything except `host_fraction`), per-sample +
   reconciled, consent gate, `probe.json` + `probe.md` + context, CLI (`metagx probe`,
   `interview --probe`), 8 tests. Modules: `metagx/probe.py`, `metagx/consent.py`.
2. **+host_fraction** when a host index + minimap2/mash are available (optional, degrades to null).
3. **+config_builder** per-sample routing defaults from the probe.
4. **MCP `run_probe`** surface (CLI exists; agent tool still TODO).
5. (Later, separate consent) external/closed-loop tier.

## 11. Open decisions

- Default consent posture when execution is available: ask-once-recommend-on vs always-opt-in
  (deployment choice; set in `SKILL.md`).
- `inferred_platform_class` thresholds (read-length/error cutoffs) — start conservative, tune from
  the existing `evidence/*.yaml`.
- Whether `dup_fraction` (first-k-bp hashing) is "reconstructive enough" to drop under the §3.4 rule
  — current judgment: a 32-bp prefix hash count is an aggregate, not a stored sequence; keep, but it
  must never store the prefixes themselves.
