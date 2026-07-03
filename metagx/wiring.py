"""Wiring-integrity audit — the cross-part "DAG" that catches a half-wired tool or module.

Adding a tool or a module touches many independent moving parts:

    registry (parameters/*.yaml)
      → config_builder (DEFAULT_MODULES, a build_config kwarg section, db.<key> acceptance)
      → workflow/Snakefile (an include guard + final targets)
      → dbprovision (SPECS) + doctor (needed_dbs presence gate)
      → tool_advisor (MODULE_TOOLS, so recommend/advise know its tools)
      → report (active_tools version capture + a CITATIONS entry)
      → mcp_server (the same build_config section, for web/MCP agents)
      → prompts/INTERVIEW.md + SKILL.md (the interview surface)

Each part is defined on its own, so it is easy to wire one and forget another — and unit tests
that exercise a single part stay green. This module introspects the parts and cross-checks them,
returning a flat list of `gaps`. `metagx wiring` prints the report (exit 1 on any gap) and
`tests/test_wiring.py` asserts there are none — so integrity is *enforced*, not remembered.

Dependency-light: stdlib + the metagx modules it audits. It parses mcp_server.py with `ast`
(no import) so it works even when the optional `mcp` serve-extra isn't installed.
"""
from __future__ import annotations

import ast
import inspect
import os
from typing import Any, Dict, List

from . import config_builder, dbprovision, registry, report, tool_advisor

# kraken2-build / bracken-build are DB-construction tools, not user-interviewed tool sections.
BUILD_ONLY_TOOLS = {"kraken2-build", "bracken-build"}

# Modules whose tool set is resolved by platform/option routing (so a static superset check
# would mis-fire); their version capture is exercised by the dry-run + e2e tests instead.
ROUTING_MODULES = {"qc", "assembly", "classify_consensus", "phylogenetics", "filtered_assembly"}

# Modules implemented in pure Python with no external tool to version (no MODULE_TOOLS entry needed).
PURE_PYTHON_MODULES = {"stats", "differential", "decontam"}

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(path: str) -> str:
    try:
        with open(os.path.join(_ROOT, path)) as fh:
            return fh.read()
    except OSError:
        return ""


def _mcp_build_config_params() -> List[str]:
    """build_config's parameter names parsed from mcp_server.py source (no import needed)."""
    src = _read("mcp_server.py")
    if not src:
        return []
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "build_config":
            a = node.args
            return [arg.arg for arg in (a.posonlyargs + a.args + a.kwonlyargs)]
    return []


def _kitchen_sink_cfg() -> Dict[str, Any]:
    """A config enabling every module (short + long reads, all domains) so report.active_tools
    is exercised across the board for the version-capture coverage check."""
    return {
        "modules": {m: True for m in config_builder.DEFAULT_MODULES},
        "samples": [{"sample": "s", "platform": "illumina", "r1": "x"},
                    {"sample": "o", "platform": "ont", "r1": "y"}],
        "domains": ["viral", "prokaryote", "eukaryote"],
        "consensus": {"classifier": "metaphlan"},
        "assembly": {"assembler": "megahit"},
        "phylogenetics": {"method": "iqtree", "input": "x"},
        "validate": {"remote": True},
    }


def audit() -> Dict[str, Any]:
    """Cross-check every moving part; return {parts: {...}, gaps: [...]}. Empty gaps == in sync."""
    gaps: List[str] = []
    tools = set(registry.list_tools())
    user_tools = tools - BUILD_ONLY_TOOLS
    modules = set(config_builder.DEFAULT_MODULES)

    # A. registry user tool  <->  config_builder.build_config section (CLI surface)
    cli_params = set(inspect.signature(config_builder.build_config).parameters)
    for t in sorted(user_tools):
        if t not in cli_params:
            gaps.append(f"registry tool '{t}' has no build_config (CLI) section kwarg")

    # B. registry user tool  <->  mcp_server.build_config section (MCP surface)
    mcp_params = set(_mcp_build_config_params())
    if mcp_params:  # only when source is present
        for t in sorted(user_tools):
            if t not in mcp_params:
                gaps.append(f"registry tool '{t}' missing from mcp_server.build_config (MCP surface)")

    # C. dbprovision.SPECS  ⊆  config_builder.DB_EXTRA_KEYS (so db.<key> is accepted)
    for s in sorted(dbprovision.SPECS):
        if s not in config_builder.DB_EXTRA_KEYS:
            gaps.append(f"dbprovision SPEC '{s}' is not an accepted db.<key> (config_builder.DB_EXTRA_KEYS)")

    # D. every module has a tool_advisor.MODULE_TOOLS entry (recommend/advise coverage)
    for m in sorted(modules):
        if m not in tool_advisor.MODULE_TOOLS:
            gaps.append(f"module '{m}' missing from tool_advisor.MODULE_TOOLS")

    # E. every module is referenced in the Snakefile (an include guard)
    snake = _read("workflow/Snakefile")
    for m in sorted(modules):
        if f'"{m}"' not in snake and f"'{m}'" not in snake:
            gaps.append(f"module '{m}' not referenced in workflow/Snakefile")

    # F. every module is named in the interview surfaces
    interview = _read("prompts/INTERVIEW.md").lower()
    skill = _read("SKILL.md").lower()
    for m in sorted(modules):
        if m.lower() not in interview:
            gaps.append(f"module '{m}' not documented in prompts/INTERVIEW.md")
        if m.lower() not in skill:
            gaps.append(f"module '{m}' not documented in SKILL.md")

    # G. version capture: enabling a non-routing module surfaces its tools in report.active_tools
    actives = set(report.active_tools(_kitchen_sink_cfg()))
    for m in sorted(modules - ROUTING_MODULES - PURE_PYTHON_MODULES):
        for t in tool_advisor.MODULE_TOOLS.get(m, []):
            if t in tools and t not in actives:
                gaps.append(f"module '{m}' tool '{t}' is not captured by report.active_tools "
                            f"(its version won't appear in the provenance manifest)")

    # H. every captured tool has a citation (so report/paper can cite it)
    for t in sorted(actives):
        if t not in report.CITATIONS:
            gaps.append(f"active tool '{t}' has no report.CITATIONS entry")

    # I. every tool ANY config can activate is provisioned by SOME conda env — the core
    # environment.yml (tools that coexist) OR a per-rule workflow/envs/*.yaml (isolated for
    # size/collision, provisioned by `metagx run --use-conda`). A tool in neither means the
    # end-user's generated scripts would die with "command not found" for that module. This is
    # the collision-/dependency-safety guarantee: every config's tools have a home.
    provisioned = _conda_provisioned_packages()
    for t in sorted(actives):
        if t == "snakemake":
            continue
        pkg = report.conda_package(t)
        if pkg not in provisioned:
            gaps.append(f"active tool '{t}' (conda pkg '{pkg}') is provisioned by NO env — not in "
                        f"environment.yml and not in any workflow/envs/*.yaml; its module can't run")

    parts = {
        "registry_tools": sorted(tools),
        "user_tools": sorted(user_tools),
        "modules": sorted(modules),
        "db_extra_keys": list(config_builder.DB_EXTRA_KEYS),
        "dbprovision_specs": sorted(dbprovision.SPECS),
        "active_tools_kitchen_sink": sorted(actives),
        "conda_provisioned_packages": sorted(provisioned),
        "mcp_build_config_known": bool(mcp_params),
    }
    return {"parts": parts, "gaps": gaps, "ok": not gaps}


def _conda_provisioned_packages() -> set:
    """Package bases across the core environment.yml + every per-rule workflow/envs/*.yaml —
    i.e. everything the pipeline can put on PATH (core env) or provision via --use-conda."""
    import glob
    import yaml as _yaml

    def bases(path: str) -> set:
        out: set = set()
        try:
            with open(path) as fh:
                data = _yaml.safe_load(fh) or {}
        except OSError:
            return out
        for dep in data.get("dependencies", []) or []:
            if isinstance(dep, str):
                out.add(dep.split()[0].split("=")[0].split(">")[0].split("<")[0].strip())
            elif isinstance(dep, dict):  # e.g. a nested pip: list
                for sub in dep.get("pip", []) or []:
                    out.add(str(sub).split("=")[0].split(">")[0].strip())
        return out

    found = bases(os.path.join(_ROOT, "environment.yml"))
    for env in glob.glob(os.path.join(_ROOT, "workflow", "envs", "*.yaml")):
        found |= bases(env)
    return found
