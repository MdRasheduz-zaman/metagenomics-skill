"""Intent-first planning: turn a research goal into modules + a database checklist.

This is the backend for the *top of the interview funnel*. The LLM classifies the user's
free-text "what are you trying to do?" into a preset (or an explicit module set — its job,
since that's fuzzy), then calls this to get the deterministic consequences:

  * which modules/domains the goal implies, and
  * **every reference DB those modules need**, each with its size and the three ways to
    satisfy it — already have it (a `db.<tool>` path), fetch/build it (`db.provision` /
    `fetch-db` / `db.build`), or (a few tools) download it manually.

The point (per the project TODO): DB needs should *fall out* of the research-question
conversation, not be a separate manual step. The LLM asks the returned `questions`, and each
answer routes straight into `build_config` kwargs — a path goes to `db.<tool>`; "fetch it"
adds the tool to `db.provision` (module DBs) or configures `db.build` (the kraken2 classifier).

The module-DB → needed-by mapping lives in `dbprovision` (single source of truth); this layer
only adds the kraken2/Bracken classifier DB (which `dbprovision.needed_dbs` deliberately omits,
since it's the `db.build`/`fetch-db` layer, not a per-tool module DB) and the routing prose.

Dependency-light (stdlib + sibling modules), so it imports cleanly in the CLI and MCP server.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

from . import dbprovision, presets, schedulers


def intake_prompt() -> Dict[str, Any]:
    """The grey-text *hint* for the opening "what are you trying to do?" question.

    Like a placeholder under a search box: it nudges the user to mention the few things that
    actually *route* the funnel, so the LLM can plan in one turn instead of a dozen round-trips.
    Every dimension is tied to where it lands in the design (and the example lists are pulled
    from the real registries — presets, platforms, schedulers — not hardcoded prose), so this
    can't drift from what the pipeline actually supports.
    """
    from .config_builder import KNOWN_PLATFORMS  # local import: avoid import cost at module load
    return {
        "prompt": "Describe your research question in plain words. It helps to mention:",
        "include": [
            {"field": "goal / research question", "routes": "preset + which modules run",
             "examples": [p["name"] for p in presets.describe_presets()]},
            {"field": "sequencing platform", "routes": "QC tool + assembler dispatch",
             # friendly labels; canonical sample-sheet values are in config_builder.KNOWN_PLATFORMS
             "examples": ["ONT", "Illumina", "PacBio HiFi", "PacBio CLR", "MGI/BGI",
                          "amplicon (16S/ITS)"],
             "config_values": sorted(KNOWN_PLATFORMS)},
            {"field": "number of samples & groups", "routes": "cross-sample comparison + "
             "differential abundance (needs >=2 groups)",
             "examples": ["12 samples", "treated vs control", "single sample"]},
            {"field": "databases you have or need", "routes": "db.<tool> path vs "
             "db.provision / db.build (I'll ask per-DB)",
             "examples": ["have a kraken2 DB at /path", "need a viral DB", "build from my genomes"]},
            {"field": "where to run", "routes": "executor / scheduler",
             "examples": schedulers.list_schedulers()},
            {"field": "run scope", "routes": "emit config.yaml + Snakemake to run yourself "
             "(BYOK handoff) vs run it end-to-end here",
             "examples": ["just generate the config + Snakemake", "run it for me"]},
        ],
        "example": ("e.g. \"12 paired-end Illumina gut samples in two groups (treated vs "
                    "control) — species profiling + differential abundance; I don't have a "
                    "kraken2 DB; just generate the config + Snakemake to run on our SLURM "
                    "cluster.\""),
    }

# Modules that consume the kraken2 + Bracken classifier DB (the db.build / fetch-db layer).
_CLASSIFIER_MODULES = ("classify", "abundance")


def _classifier_db_entry() -> Dict[str, Any]:
    """The kraken2/Bracken classifier DB — not a per-tool module DB, so dbprovision omits it."""
    return {
        "name": "kraken2",
        "kind": "classifier",
        "needed_by": "classify / abundance (Bracken reuses the kraken2 DB)",
        "size": "varies (prebuilt standard-8 ~6 GB; must fit in RAM)",
        "config_key": "db.kraken2",
        "have_it": "set db.kraken2: <path> (Bracken defaults to the same path)",
        "fetch": "metagx fetch-db standard-8 --dir <dir>  (prebuilt index, fast — preferred)",
        "build": ("add a db.build block (strategies: standard/custom-fasta/custom-folder/"
                  "spike-in) — best for viral/custom/spike-in; NCBI rsync is dead so standard "
                  "downloads are slow"),
        "manual": False,
    }


def _module_db_entry(tool: str, key: str) -> Dict[str, Any]:
    spec = dbprovision.SPECS.get(tool, {})
    manual = bool(spec.get("manual"))
    entry: Dict[str, Any] = {
        "name": tool,
        "kind": "module",
        "needed_by": spec.get("needed_by", ""),
        "size": spec.get("size", "unknown"),
        "config_key": f"db.{key}",
        "have_it": f"set db.{key}: <path>",
        "manual": manual,
    }
    if manual:
        entry["fetch"] = f"manual download — {spec.get('docs', 'see the tool docs')}"
    else:
        entry["fetch"] = (
            f"add '{tool}' to db.provision (auto-fetch at run time, idempotent) "
            f"OR metagx fetch-db --tool {tool} --dir <dir>"
        )
    if spec.get("self_gates"):
        entry["note"] = "the rule self-gates on this DB, so a missing one is a skip (not a failure)"
    # GTDB-Tk is the one DB you must never auto-fetch on a small machine.
    if tool == "gtdbtk":
        entry["warn"] = "~110 GB and no small build — confirm disk before fetching; never auto-fetch on a laptop"
    return entry


def _question_for(entry: Dict[str, Any]) -> str:
    name = entry["name"]
    nb = entry.get("needed_by") or "your enabled modules"
    size = entry.get("size", "")
    if entry["kind"] == "classifier":
        return (f"Do you already have a kraken2 database (for {nb})? "
                f"If yes, give the path → {entry['config_key']}. "
                f"If no, I can fetch a prebuilt index ({entry['fetch']}) "
                f"or build one ({entry['build']}).")
    if entry.get("manual"):
        return (f"The {name} DB ({size}, needed by {nb}) has no automatic downloader. "
                f"Do you have it? If yes, path → {entry['config_key']}; "
                f"if no, {entry['fetch']}.")
    q = (f"Do you already have a {name} DB ({size}, needed by {nb})? "
         f"If yes, path → {entry['config_key']}; if no, I can fetch it "
         f"(add '{name}' to db.provision).")
    if entry.get("warn"):
        q += f"  ⚠ {entry['warn']}"
    return q


def plan(
    *,
    preset: Optional[str] = None,
    modules: Optional[Dict[str, bool]] = None,
    domains: Optional[List[str]] = None,
    functional: Optional[List[str]] = None,
    have: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Resolve a goal (a preset and/or explicit modules) into modules + a DB checklist.

    Args:
      preset:     a preset name (its modules/params are the base).
      modules:    explicit module toggles, deep-merged over the preset (user wins).
      domains:    domain list for the domain_taxonomy module (viral/prokaryote/eukaryote).
      functional: enabled functional sub-steps (annotation/amr/pathways) — drives which
                  functional DBs are needed.
      have:       tool names the user already has a DB path for; marked resolved (not asked).

    Returns a dict with the effective modules/domains, a `databases` checklist (each entry
    self-describing how to satisfy it), and ready-to-ask `questions` for the unresolved DBs.
    """
    base_modules: Dict[str, bool] = {}
    summary = ""
    if preset:
        cfg = presets.get_preset_config(preset)  # raises KeyError on unknown preset
        base_modules = dict(cfg.get("modules", {}))
        if not domains and cfg.get("domains"):
            domains = list(cfg["domains"])
        for desc in presets.describe_presets():
            if desc["name"] == preset:
                summary = desc.get("when_to_use") or desc.get("description") or ""
                break
    mods = {**base_modules, **(modules or {})}
    mods = {k: v for k, v in mods.items() if v}  # keep only enabled, for a clean echo
    doms = [str(d).lower() for d in (domains or [])]
    func = {str(f).lower(): True for f in (functional or [])}

    have_set = {str(h).lower() for h in (have or [])}

    databases: List[Dict[str, Any]] = []
    # 1) the kraken2/Bracken classifier DB (if any classifying module is on)
    if any(mods.get(m) for m in _CLASSIFIER_MODULES):
        databases.append(_classifier_db_entry())
    # 2) per-tool module DBs — dbprovision owns the module→DB mapping (single source of truth)
    planned_cfg = {"modules": mods, "domains": doms, "functional": func}
    for tool, key in dbprovision.needed_dbs(planned_cfg).items():
        databases.append(_module_db_entry(tool, key))

    for entry in databases:
        entry["resolved"] = entry["name"] in have_set

    questions = [_question_for(e) for e in databases if not e["resolved"]]

    result: Dict[str, Any] = {
        "preset": preset,
        "summary": summary,
        "modules": mods,
        "domains": doms,
        "functional": sorted(func),
        "databases": databases,
        "questions": questions,
    }
    # Opening state (no goal resolved yet): include the intake hint so the LLM can ask a
    # well-framed "what are you trying to do?" with the prompts that route the funnel.
    if not mods:
        result["intake"] = intake_prompt()
    return result
