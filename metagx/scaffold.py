"""Scaffold a parameter-registry stub from a tool's live ``--help``.

The output is a STARTING POINT a human must curate, never a finished registry. It makes a
tool *capability-complete* cheaply: every flag becomes an ``ask: false`` / ``tier: 3`` entry
with a guessed type, plus an ``extra_args`` passthrough valve for the long tail. The curator
then: promotes the few flags worth interviewing (``ask: true`` + a real ``question``), marks
workflow-owned flags ``managed: true``, fixes guessed types/defaults/bounds, and adds
``recommend`` / ``warn_if`` evidence. See ``flye.yaml`` for a curated result.

The parsing (``from_help_text``) is pure and unit-tested; ``scaffold`` wraps it with a live
``capture_help`` so it also works against an installed tool.
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple

import yaml

from . import registry, sync_help

# Flags that are never analysis parameters.
SKIP_FLAGS = {"-h", "--help", "-v", "--version", "--usage", "--citation"}

# Note: SIZE is intentionally absent — tools use it for human sizes like "5m"/"2.6g" (str).
_INT_METAVARS = {"INT", "N", "NUM", "NUMBER", "COUNT", "LEN", "LENGTH",
                 "K", "BP", "SEED", "MIN", "MAX"}
_FLOAT_METAVARS = {"FLOAT", "FRAC", "FRACTION", "RATE", "PROB", "PROBABILITY",
                   "CUTOFF", "THRESHOLD", "PERCENT"}
_PATH_METAVARS = {"PATH", "FILE", "DIR", "DIRECTORY", "FASTA", "FASTQ", "FQ", "FA",
                  "DB", "BAM", "SAM", "BED", "GFF", "GTF", "VCF", "OUTPUT", "OUT", "PREFIX"}
_PATH_LOWER = {"path", "file", "dir", "directory", "fasta", "fastq", "prefix"}

_HEADER = """\
# AUTO-SCAFFOLDED stub from `{command} --help` — CURATE BEFORE USE.
#
# Every flag below is ask:false / tier:3 with a *guessed* type. This is the
# capability-complete surface, not the interview funnel. Before shipping:
#   - set ask:true + a real `question` on the few flags worth interviewing,
#   - mark workflow-owned flags `managed: true` (db/threads/io/mutually-exclusive inputs),
#   - fix guessed types, add min/max bounds, defaults, choices,
#   - add `recommend:` / `warn_if:` evidence where you have it.
# See metagx/parameters/kraken2.yaml and flye.yaml for the documented schema + a curated result.
"""


def _type_from_metavar(metavar: str) -> str:
    token = metavar.strip("<>[](). ")
    up = token.upper()
    if up in _INT_METAVARS:
        return "int"
    if up in _FLOAT_METAVARS:
        return "float"
    if up in _PATH_METAVARS or token in _PATH_LOWER:
        return "path"
    return "str"


def _extract_choices(metavar: Optional[str], desc: str) -> Optional[List[str]]:
    for text in (metavar or "", desc):
        m = re.search(r"\{([^}]+)\}", text)
        if m:
            items = [c.strip() for c in m.group(1).split(",") if c.strip()]
            if len(items) >= 2:
                return items
    return None


def _numeric_default(desc: str, ptype: str) -> Optional[Any]:
    cands: List[str] = []
    for pat in (r"\(default[:=]?\s*([^)]+)\)", r"\[default[:=]?\s*([^\]]+)\]", r"\[([^\]]+)\]\s*$"):
        m = re.search(pat, desc, re.I)
        if m:
            cands.append(m.group(1).strip())
    for c in cands:
        try:
            return int(c) if ptype == "int" else float(c)
        except ValueError:
            continue  # "[auto]", "[not set]", etc. -> leave unset
    return None


def _parse_optspec(optpart: str) -> Optional[Tuple[str, Optional[str]]]:
    """From an option spec like ``-t INT, --threads INT`` -> (canonical_flag, metavar)."""
    flags: List[str] = []
    metavar: Optional[str] = None
    # split flags on commas, but not commas inside an enum metavar like {fast,normal}
    for tok in re.split(r",(?![^{]*\})", optpart):
        parts = tok.split()
        if not parts or not parts[0].startswith("-"):
            continue
        flags.append(parts[0])
        if metavar is None and len(parts) > 1:
            metavar = " ".join(parts[1:])
    if not flags:
        return None
    long_flags = [f for f in flags if f.startswith("--")]
    flag = max(long_flags, key=len) if long_flags else flags[0]
    return flag, metavar


def _make_param(flag: str, metavar: Optional[str], desc: str) -> Dict[str, Any]:
    choices = _extract_choices(metavar, desc)
    if choices:
        ptype = "enum"
    elif metavar:
        ptype = _type_from_metavar(metavar.split()[0])
    else:
        ptype = "bool"

    spec: Dict[str, Any] = {"flag": flag, "type": ptype}
    if ptype == "bool":
        spec["default"] = False
    elif ptype in ("str", "path"):
        spec["default"] = ""
    else:  # int / float
        dflt = _numeric_default(desc, ptype)
        if dflt is not None:
            spec["default"] = dflt
    if choices:
        spec["choices"] = choices
    spec["tier"] = 3
    spec["ask"] = False
    spec["question"] = " ".join(desc.split()) or f"(scaffolded) {flag}"
    return spec


def parse_params(help_text: str) -> Dict[str, Dict[str, Any]]:
    """Best-effort parse of every option in ``--help`` into registry param stubs."""
    lines = help_text.splitlines()
    out: Dict[str, Dict[str, Any]] = {}
    i, n = 0, len(lines)
    while i < n:
        raw = lines[i]
        i += 1
        if not raw.strip().startswith("-"):
            continue
        indent = len(raw) - len(raw.lstrip())
        segs = re.split(r"\s{2,}", raw.strip(), maxsplit=1)
        optpart, desc = segs[0], (segs[1].strip() if len(segs) > 1 else "")
        # description wrapped onto the next, more-indented, non-option line
        if not desc and i < n:
            nxt = lines[i]
            nindent = len(nxt) - len(nxt.lstrip())
            if nxt.strip() and not nxt.strip().startswith("-") and nindent > indent:
                desc = nxt.strip()
                i += 1
        parsed = _parse_optspec(optpart)
        if not parsed:
            continue
        flag, metavar = parsed
        if flag in SKIP_FLAGS:
            continue
        name = flag.lstrip("-").replace("-", "_")
        if not name or name in out:
            continue
        out[name] = _make_param(flag, metavar, desc)
    return out


def from_help_text(help_text: str, command: str, name: Optional[str] = None,
                   version: Optional[str] = None) -> str:
    """Render a registry-stub YAML string from raw help text (pure / testable)."""
    tool = name or command.split()[0]
    params = parse_params(help_text)
    # the raw escape hatch for everything we did not model (run-control, future flags)
    params["extra_args"] = {
        "type": "str",
        "default": "",
        "tier": 3,
        "ask": False,
        "passthrough": True,
        "question": f"Additional raw {tool} command-line tokens, passed verbatim "
                    "(escape hatch; not validated, logged in provenance).",
    }
    reg: Dict[str, Any] = {
        "tool": tool,
        "command": command,
        "version_probe": f"{command} --version",
        "description": f"TODO: one-line description of {tool} (scaffolded).",
        "params": params,
    }
    body = yaml.safe_dump(reg, sort_keys=False, default_flow_style=False, width=100)
    return _HEADER.format(command=command) + body


def scaffold(command: str, name: Optional[str] = None) -> Dict[str, Any]:
    """Capture ``command --help`` live and render a stub. Returns a result dict."""
    cap = sync_help.capture_help(command)
    if not cap.get("ok"):
        return {"ok": False, "error": cap.get("error"), "command": command}
    tool = name or command.split()[0]
    existing = tool in registry.list_tools()
    yaml_text = from_help_text(cap["help_text"], command, name=tool, version=cap.get("version"))
    return {
        "ok": True,
        "command": command,
        "tool": tool,
        "version": cap.get("version"),
        "registry_exists": existing,
        "yaml": yaml_text,
    }
