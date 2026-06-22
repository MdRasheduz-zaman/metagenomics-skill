"""Capture live ``tool --help`` and diff against parameter registries."""

from __future__ import annotations

import re
import shutil
import subprocess
from typing import Any, Dict, List, Optional

from . import registry

FLAG_RE = re.compile(r"^\s+(-{1,2}[\w-]+)(?:[,\s]+(-{1,2}[\w-]+))?\s+(.*)$")


def capture_help(command: str, timeout: int = 30) -> Dict[str, Any]:
    """Run ``command --help`` and ``command --version`` when available."""
    exe = command.split()[0]
    if not shutil.which(exe):
        return {"command": command, "ok": False, "error": f"{exe} not found on PATH"}
    try:
        proc = subprocess.run(
            [exe, "--help"],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        help_text = proc.stdout or proc.stderr or ""
    except (subprocess.TimeoutExpired, OSError) as e:
        return {"command": command, "ok": False, "error": str(e)}

    version = None
    try:
        vproc = subprocess.run(
            [exe, "--version"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        version = (vproc.stdout or vproc.stderr or "").strip().split("\n")[0]
    except (subprocess.TimeoutExpired, OSError):
        version = None

    return {
        "command": command,
        "ok": True,
        "version": version,
        "help_text": help_text,
        "flags": parse_help_flags(help_text),
    }


def parse_help_flags(help_text: str) -> List[Dict[str, str]]:
    """Best-effort parse of CLI flags from help output."""
    flags: List[Dict[str, str]] = []
    seen = set()
    for line in help_text.splitlines():
        m = FLAG_RE.match(line)
        if not m:
            continue
        for g in (m.group(1), m.group(2)):
            if not g or g in seen:
                continue
            seen.add(g)
            flags.append({"flag": g, "description": (m.group(3) or "").strip()})
    return flags


def diff_registry(tool: str, help_capture: Dict[str, Any] | None = None) -> Dict[str, Any]:
    """Compare registry flags to live --help (when capture succeeded)."""
    reg = registry.load_registry(tool)
    command = reg.get("command", tool)
    cap = help_capture if help_capture is not None else capture_help(command)

    reg_flags = {
        spec["flag"]
        for spec in reg.get("params", {}).values()
        if spec.get("flag") and not spec.get("managed")
    }
    help_flags = {f["flag"] for f in cap.get("flags", [])} if cap.get("ok") else set()
    # Registries are deliberately curated, not exhaustive: a small typed funnel plus one
    # `passthrough: true` valve (extra_args) for the long tail. So a help flag absent from the
    # typed params is NOT unreachable when a passthrough valve exists — it's just unmodeled.
    has_passthrough = any(
        spec.get("passthrough") for spec in reg.get("params", {}).values()
    )

    missing_in_registry = sorted(help_flags - reg_flags) if help_flags else []
    registry_only = sorted(reg_flags - help_flags) if help_flags else []

    return {
        "tool": tool,
        "command": command,
        "capture_ok": cap.get("ok", False),
        "version": cap.get("version"),
        "error": cap.get("error"),
        "registry_user_flags": sorted(reg_flags),
        "help_flags_count": len(help_flags),
        "missing_in_registry": missing_in_registry,
        "in_registry_not_in_help": registry_only,
        "has_passthrough": has_passthrough,
        # With a passthrough valve every flag is still reachable (via extra_args), so the
        # unmodeled flags are informational — curate them up into typed params only if users
        # need bounds/recommendations on them. Without it, they are genuine reachability gaps.
        "unmodeled_reachable_via_passthrough": has_passthrough and bool(missing_in_registry),
        "metadata": registry.tool_metadata(tool),
    }


def sync_all(tools: Optional[List[str]] = None) -> Dict[str, Any]:
    """Diff every registry tool against live --help."""
    names = tools or registry.list_tools()
    results = {t: diff_registry(t) for t in names}
    n_ok = sum(1 for r in results.values() if r.get("capture_ok"))
    n_drift = sum(1 for r in results.values() if r.get("missing_in_registry"))
    # Genuine reachability gaps: unmodeled help flags on a tool with NO passthrough valve.
    # (Every metagx registry currently has one, so this should be empty — that's the invariant.)
    n_unreachable = sum(
        1 for r in results.values()
        if r.get("capture_ok") and r.get("missing_in_registry") and not r.get("has_passthrough")
    )
    return {
        "tools_checked": len(names),
        "capture_ok": n_ok,
        "tools_with_help_flags_missing_from_registry": n_drift,
        "tools_with_unreachable_flags": n_unreachable,
        "tools": results,
    }
