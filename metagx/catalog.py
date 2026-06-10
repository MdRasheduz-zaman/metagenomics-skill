"""Index of registries, evidence files, and workflow scripts for LLM discovery."""

from __future__ import annotations

import glob
import os
from typing import Any, Dict

from . import evidence_pack, registry


def workflow_scripts(root: str | None = None) -> list:
    here = root or os.path.join(os.path.dirname(os.path.dirname(__file__)), "workflow", "scripts")
    return sorted(os.path.basename(p) for p in glob.glob(os.path.join(here, "*.py")))


def build_catalog() -> Dict[str, Any]:
    tools = {}
    for t in registry.list_tools():
        reg = registry.load_registry(t)
        tools[t] = {
            "description": (reg.get("description") or "").strip(),
            "command": reg.get("command"),
            "metadata": registry.tool_metadata(t),
            "interview_params": len(registry.interview_spec(t, max_tier=2)),
        }
    return {
        "registry_tools": tools,
        "evidence_files": evidence_pack.list_evidence(),
        "workflow_scripts": workflow_scripts(),
        "planned_modules": [],
        "advisor_commands": ["recommend", "advise", "history", "sync-help", "catalog"],
    }
