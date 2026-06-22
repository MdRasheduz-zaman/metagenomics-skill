"""Tool version + valid-argument locking.

Bioinformatics tools drift in ways that silently corrupt a pipeline: a flag is renamed or
removed (`iqtree2` → `iqtree3`), an option changes meaning, a subcommand moves, or a tool
disappears from PATH. The registries model flags, but nothing pins them to the *installed*
version. This module closes that gap with two capabilities:

  * **lock** — snapshot, per tool, the resolved binary, its reported version, and the set of
    flags its live `--help` accepts. Written to a lockfile so a run records exactly what it
    ran against (provenance) and a later environment can be checked for drift.
  * **verify / config-flag check** — diff the current install against a lock (version changed,
    flags that vanished), and — config-aware — confirm every flag a given config actually sets
    is present in the installed tool's `--help`, so a removed/renamed flag fails fast with a
    clear message instead of a cryptic tool error mid-run.

The diff/check logic is pure (operates on dicts), so it is unit-tested without any tool
installed; the probing reuses `sync_help.capture_help` (live `--help` + `--version`).
"""
from __future__ import annotations

import shutil
from typing import Any, Callable, Dict, List, Optional

from . import registry, sync_help


def probe_tool(tool: str, capture: Optional[Callable[[str], Dict]] = None) -> Dict[str, Any]:
    """Resolve + version + accepted flags for one registry tool's command."""
    reg = registry.load_registry(tool)
    command = reg.get("command", tool)
    exe = command.split()[0]
    resolved = shutil.which(exe)
    if capture:
        cap = capture(command)
    else:
        cap = sync_help.capture_help(command, version_cmd=reg.get("version_probe"))
    return {
        "tool": tool,
        "command": command,
        "resolved": resolved,
        "installed": bool(resolved),
        "version": cap.get("version"),
        "flags": sorted({f["flag"] for f in cap.get("flags", [])}) if cap.get("ok") else [],
        "help_ok": bool(cap.get("ok")),
    }


def snapshot(tools: Optional[List[str]] = None,
             capture: Optional[Callable[[str], Dict]] = None) -> Dict[str, Any]:
    """A lockable snapshot of every (or the given) registry tool's version + accepted flags."""
    names = tools or registry.list_tools()
    return {"tools": {t: probe_tool(t, capture=capture) for t in names}}


def diff_lock(locked: Dict[str, Any], current: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Drift findings between a saved lock and a current snapshot. Pure (no IO).

    Severity: ``error`` = the lock can no longer be honored (tool gone, or a flag the lock
    recorded has vanished — likely a rename/removal that breaks rendered commands);
    ``warn`` = version changed but flags still cover the lock (probably fine, worth noting).
    """
    out: List[Dict[str, Any]] = []
    lt, ct = locked.get("tools", {}), current.get("tools", {})
    for tool, lk in lt.items():
        if not lk.get("installed"):
            continue  # wasn't installed when locked — nothing was pinned, so nothing to drift
        cur = ct.get(tool)
        if not cur or not cur.get("installed"):
            out.append({"tool": tool, "severity": "error", "kind": "missing",
                        "message": f"{tool} ({lk.get('command')}) was locked but is not on PATH now"})
            continue
        if lk.get("version") and cur.get("version") and lk["version"] != cur["version"]:
            out.append({"tool": tool, "severity": "warn", "kind": "version",
                        "message": f"{tool} version changed: locked {lk['version']!r} -> "
                                   f"now {cur['version']!r}"})
        if lk.get("flags") and cur.get("help_ok"):
            gone = sorted(set(lk["flags"]) - set(cur.get("flags", [])))
            if gone:
                out.append({"tool": tool, "severity": "error", "kind": "flags_removed",
                            "message": f"{tool}: flags in the lock are gone from this version's "
                                       f"--help (renamed/removed?): {gone}"})
    return out


def config_flag_check(cfg: Dict[str, Any],
                      capture: Optional[Callable[[str], Dict]] = None) -> List[Dict[str, Any]]:
    """For every tool section in a config, confirm the flags it would emit exist in the live
    tool's --help. Catches "you set --foo but this installed version dropped/renamed it" before
    the run. Skips tools not on PATH (nothing to validate) and passthrough/extra_args (unmodeled).
    """
    out: List[Dict[str, Any]] = []
    cap = capture or sync_help.capture_help
    for tool in registry.list_tools():
        section = cfg.get(tool)
        if not section:
            continue
        reg = registry.load_registry(tool)
        command = reg.get("command", tool)
        if not shutil.which(command.split()[0]):
            continue
        try:
            tokens = registry.render_args(tool, section)
        except registry.ValidationError:
            continue  # a bad value is config-validation's job, not ours
        used = [t for t in tokens if t.startswith("-")]
        help_cap = cap(command)
        if not help_cap.get("ok"):
            continue
        have = {f["flag"] for f in help_cap.get("flags", [])}
        missing = sorted(f for f in used if f not in have)
        if missing:
            out.append({"tool": tool, "severity": "error", "kind": "flag_not_in_help",
                        "version": help_cap.get("version"),
                        "message": f"{tool}: config sets {missing} but the installed version "
                                   f"({help_cap.get('version')}) does not list them in --help "
                                   f"— flag renamed/removed, or wrong tool version"})
    return out
