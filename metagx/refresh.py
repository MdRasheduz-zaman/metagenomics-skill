"""`metagx refresh <tool>` — offline authoring aid that proposes registry updates from the
*installed* binary, for a human to review and commit.

Keeping a parameter registry current is version-dependent toil: a tool gains/renames/removes
flags between versions, and nothing tells the maintainer the curated params drifted from what is
installed. This orchestrates the pieces that already exist — it does **not** reimplement them:

    sync_help.capture_help / diff_registry  → installed version + flag diff (binary = truth)
    scaffold.parse_params                   → registry stubs for the new (unmodeled) flags
    registry.version_info                    → the version the registry was curated against

and emits a **reviewable proposal**: new flags scaffolded as ``_status: proposed`` drafts (inert
until a human curates them — see ``registry._is_proposed``), registry flags that vanished from the
binary flagged as drift, and a proposed ``tested_version`` bump.

Design constraints (deliberate):
  * **Pure core, no writes.** ``plan_refresh`` is a pure function over an injected capture, so it
    unit-tests with no tool installed (mirrors ``toollock``).
  * **Never mutates the curated registry.** The curated YAML carries hand-written comments a
    round-trip would destroy; ``write_proposal`` writes a *sidecar* (``refresh/<tool>.proposed
    .yaml``) the human pastes from. The tool is advisory — the command path stays committed + reviewed.
  * **Binary is ground truth for the flag contract; official docs (``docs_url``) are for prose.**
    The LLM-drafting of ``question``/``recommend``/``conflicts`` from the docs is a separate,
    human-gated curation step (the skill fills the ``<LLM draft — REVIEW>`` placeholders), not a
    hidden call here.
"""
from __future__ import annotations

import os
from typing import Any, Callable, Dict, Optional

import yaml

from . import registry, scaffold, sync_help

_DRAFT = "<LLM draft — REVIEW>"


def version_token(version_string: Optional[str]) -> Optional[str]:
    """Pull a dotted version token (e.g. ``2.17.1``) out of a tool's raw ``--version`` line.
    Shared with ``doctor`` so the drift check parses versions the same way ``refresh`` does."""
    if not version_string:
        return None
    import re
    m = re.search(r"\d+(?:\.\d+)+", version_string)
    return m.group(0) if m else None


def plan_refresh(tool: str, capture: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Compute a proposed registry update for ``tool`` — pure, no IO when ``capture`` is given.

    ``capture`` is a ``sync_help.capture_help``-shaped dict ({ok, version, help_text, flags});
    inject it in tests. When ``None`` the installed binary is probed live.

    Returns::

        {tool, capture_ok, error, version: {installed, tested, differs},
         new_params: {name: proposed-stub, ...}, removed_flags: [...], summary: str}
    """
    reg = registry.load_registry(tool)
    command = registry.resolve_command(tool)  # find iqtree3 when the registry says iqtree2, etc.
    if capture is None:
        capture = sync_help.capture_help(command, version_cmd=reg.get("version_probe"))

    tested = registry.version_info(tool).get("tested_version")
    if not capture.get("ok"):
        return {
            "tool": tool, "command": command, "capture_ok": False,
            "error": capture.get("error") or "tool not installed / help not parsed",
            "version": {"installed": None, "tested": tested, "differs": False},
            "new_params": {}, "removed_flags": [],
            "summary": f"{tool}: cannot probe installed binary — nothing to propose.",
        }

    diff = sync_help.diff_registry(tool, help_capture=capture)
    # diff_registry excludes *managed* flags from its registry set (its job is user-flag
    # reachability), so a managed flag like --db/--threads shows up as "missing". For refresh the
    # question is "what is unmodeled *at all*", so subtract every flag already in the registry,
    # managed included — otherwise we'd propose re-adding a managed flag as a new param.
    modeled = {s.get("flag") for s in reg.get("params", {}).values() if s.get("flag")}
    missing = set(diff.get("missing_in_registry") or []) - modeled
    removed = list(diff.get("in_registry_not_in_help") or [])

    # scaffold stubs for the whole help surface, then keep only the genuinely-new flags and mark
    # them proposed so they stay inert until curated.
    new_params: Dict[str, Any] = {}
    for name, spec in scaffold.parse_params(capture.get("help_text", "")).items():
        if spec.get("flag") in missing:
            new_params[name] = {**spec, "_status": "proposed", "question": _DRAFT}

    installed_raw = capture.get("version")
    inst_tok, test_tok = version_token(installed_raw), version_token(tested)
    differs = bool(test_tok) and inst_tok is not None and inst_tok != test_tok

    parts = []
    if new_params:
        parts.append(f"{len(new_params)} new flag(s) to curate")
    if removed:
        parts.append(f"{len(removed)} registry flag(s) not in installed --help")
    if differs:
        parts.append(f"version drift {test_tok} → {inst_tok}")
    summary = f"{tool}: " + ("; ".join(parts) if parts else "in sync with installed binary.")

    return {
        "tool": tool, "command": command, "capture_ok": True, "error": None,
        "version": {"installed": installed_raw, "tested": tested, "differs": differs},
        "new_params": new_params, "removed_flags": removed, "summary": summary,
    }


def render_sidecar(proposal: Dict[str, Any]) -> str:
    """The paste-ready YAML block of proposed param stubs (with a review banner)."""
    tool = proposal["tool"]
    header = (f"# PROPOSED registry additions for {tool} — REVIEW before pasting into\n"
              f"# metagx/parameters/{tool}.yaml. Each entry is `_status: proposed` (inert) until\n"
              f"# you curate its type/bounds/tier, replace the `{_DRAFT}` question, add any\n"
              f"# recommend/warn_if/conflicts, and remove the `_status` marker.\n")
    body = yaml.safe_dump({"params": proposal.get("new_params", {})},
                          sort_keys=False, default_flow_style=False)
    return header + body


def write_proposal(proposal: Dict[str, Any], root: str,
                   writer: Optional[Callable[[str, str], None]] = None) -> Optional[str]:
    """Write the paste-ready sidecar to ``<root>/refresh/<tool>.proposed.yaml``; return its path.

    Never touches the curated registry (see module docstring). Returns ``None`` when there is
    nothing new to propose. ``writer`` is injectable for tests (default writes to disk)."""
    if not proposal.get("new_params"):
        return None
    outdir = os.path.join(root, "refresh")
    path = os.path.join(outdir, f"{proposal['tool']}.proposed.yaml")
    content = render_sidecar(proposal)
    if writer is not None:
        writer(path, content)
    else:
        os.makedirs(outdir, exist_ok=True)
        with open(path, "w") as fh:
            fh.write(content)
    return path
