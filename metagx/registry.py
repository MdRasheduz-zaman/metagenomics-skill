"""Load the per-tool parameter registries and derive everything from them.

This module is intentionally dependency-light (only PyYAML) so it can be imported
both inside the Snakemake workflow and by the MCP server / CLI.
"""

from __future__ import annotations

import functools
import shlex
from importlib import resources
from typing import Any, Dict, List, Tuple

import yaml


class ValidationError(ValueError):
    """Raised when supplied parameter values violate the registry schema."""


@functools.lru_cache(maxsize=None)
def _load_yaml(tool: str) -> Dict[str, Any]:
    pkg = resources.files("metagx.parameters")
    path = pkg / f"{tool}.yaml"
    if not path.is_file():
        raise ValidationError(f"Unknown tool '{tool}'. Available: {', '.join(list_tools())}")
    return yaml.safe_load(path.read_text())


def list_tools() -> List[str]:
    """Names of every tool that has a parameter registry."""
    pkg = resources.files("metagx.parameters")
    return sorted(p.name[:-5] for p in pkg.iterdir() if p.name.endswith(".yaml"))


def load_registry(tool: str) -> Dict[str, Any]:
    """Full registry dict for a tool (description + params)."""
    return _load_yaml(tool)


def tool_metadata(tool: str) -> Dict[str, Any]:
    """Upstream docs / version probe fields (optional per registry)."""
    reg = load_registry(tool)
    return {
        k: reg[k]
        for k in ("source_repo", "docs_url", "version_probe")
        if reg.get(k)
    }


def _params(tool: str) -> Dict[str, Any]:
    return load_registry(tool)["params"]


# --------------------------------------------------------------------------- #
# Interview                                                                    #
# --------------------------------------------------------------------------- #
_CMP = {
    "_gte": lambda a, b: a >= b,
    "_lte": lambda a, b: a <= b,
    "_gt": lambda a, b: a > b,
    "_lt": lambda a, b: a < b,
}


def when_matches(when: Dict[str, Any], context: Dict[str, Any]) -> bool:
    """True when every clause in ``when`` holds against ``context`` (AND).

    Key conventions: a bare key (``goal``) means equality (case-insensitive for
    strings); a key with a ``_gte``/``_lte``/``_gt``/``_lt`` suffix (``estimated_bases_gte``)
    compares the numeric context value for the base key. A missing or ``None`` context
    value never matches — we don't promote on unknown facts.
    """
    if not when:
        return False
    for key, expected in when.items():
        suffix = next((s for s in _CMP if key.endswith(s)), None)
        if suffix:
            base = key[: -len(suffix)]
            actual = context.get(base)
            if actual is None:
                return False
            try:
                if not _CMP[suffix](float(actual), float(expected)):
                    return False
            except (TypeError, ValueError):
                return False
        else:
            actual = context.get(key)
            if actual is None:
                return False
            if isinstance(expected, str) and isinstance(actual, str):
                if actual.strip().lower() != expected.strip().lower():
                    return False
            elif actual != expected:
                return False
    return True


def _promotion(spec: Dict[str, Any], context: Dict[str, Any] | None) -> Dict[str, Any] | None:
    """First matching ``promote_when`` rule for a param, or None."""
    if not context:
        return None
    for rule in spec.get("promote_when") or []:
        if when_matches(rule.get("when") or {}, context):
            return {"to_tier": int(rule.get("to_tier", 1)), "reason": (rule.get("reason") or "").strip()}
    return None


def interview_spec(tool: str, max_tier: int = 2,
                   context: Dict[str, Any] | None = None) -> List[Dict[str, Any]]:
    """User-facing parameters an LLM should consider asking about.

    Returns params with ``ask: true`` at or below ``max_tier`` (1=core, 2=common,
    3=advanced). When ``context`` is given (e.g. ``{"goal": "strain_resolved",
    "estimated_bases": 6e10}``), a normally-quiet param whose ``promote_when`` matches
    is pulled into the funnel at the rule's ``to_tier`` and carries a ``promoted`` note
    explaining why. Managed params (db/threads/io) are never returned. Each entry carries
    everything the LLM needs to phrase a good question and validate the answer.
    """
    out: List[Dict[str, Any]] = []
    for name, spec in _params(tool).items():
        if spec.get("managed"):
            continue
        promo = _promotion(spec, context)
        ask = spec.get("ask", False) or promo is not None
        if not ask:
            continue
        tier = promo["to_tier"] if promo else spec.get("tier", 3)
        if tier > max_tier:
            continue
        entry: Dict[str, Any] = {
            "name": name,
            "type": spec["type"],
            "default": spec.get("default"),
            "min": spec.get("min"),
            "max": spec.get("max"),
            "choices": spec.get("choices"),
            "tier": tier,
            "sweepable": bool(spec.get("sweepable", False)),
            "question": " ".join(str(spec.get("question", "")).split()),
        }
        if spec.get("recommend"):
            entry["recommend"] = spec["recommend"]
        if spec.get("warn_if"):
            entry["warn_if"] = spec["warn_if"]
        if promo:
            entry["promoted"] = promo  # surfaced only because context matched; carries the reason
        out.append(entry)
    out.sort(key=lambda p: p["tier"])
    return out


# --------------------------------------------------------------------------- #
# Validation                                                                   #
# --------------------------------------------------------------------------- #
def _coerce(name: str, spec: Dict[str, Any], value: Any) -> Any:
    t = spec["type"]
    try:
        if t == "bool":
            if isinstance(value, bool):
                return value
            return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}
        if t == "int":
            return int(value)
        if t == "float":
            return float(value)
        # str / path / enum
        return str(value)
    except (TypeError, ValueError):
        raise ValidationError(f"'{name}': expected {t}, got {value!r}")


def _check_one(tool: str, name: str, spec: Dict[str, Any], value: Any) -> Any:
    value = _coerce(name, spec, value)
    if spec["type"] in ("int", "float"):
        lo, hi = spec.get("min"), spec.get("max")
        if lo is not None and value < lo:
            raise ValidationError(f"'{name}'={value} is below minimum {lo}")
        if hi is not None and value > hi:
            raise ValidationError(f"'{name}'={value} is above maximum {hi}")
    if spec["type"] == "enum":
        choices = spec.get("choices", [])
        if value and value not in choices:
            raise ValidationError(f"'{name}'={value!r} not in choices {choices}")
    return value


def validate(tool: str, values: Dict[str, Any]) -> Dict[str, Any]:
    """Validate & coerce a dict of user values. Returns the cleaned dict.

    Sweepable params may be given as a list; every element is validated. Unknown
    or managed keys raise (managed values are injected by the workflow, not the user).
    """
    params = _params(tool)
    cleaned: Dict[str, Any] = {}
    for name, value in values.items():
        if name not in params:
            raise ValidationError(
                f"'{name}' is not a {tool} parameter. "
                f"See `metagx params {tool}` for valid names."
            )
        spec = params[name]
        if spec.get("managed"):
            raise ValidationError(
                f"'{name}' is managed by the workflow and cannot be set manually."
            )
        if isinstance(value, list):
            if not spec.get("sweepable"):
                raise ValidationError(f"'{name}' is not sweepable; give a single value.")
            cleaned[name] = [_check_one(tool, name, spec, v) for v in value]
        else:
            cleaned[name] = _check_one(tool, name, spec, value)
    return cleaned


# --------------------------------------------------------------------------- #
# Command-line construction                                                    #
# --------------------------------------------------------------------------- #
def render_args(tool: str, values: Dict[str, Any], managed: Dict[str, Any] | None = None) -> List[str]:
    """Build a flat list of CLI tokens from user + managed values.

    ``values``  : validated user parameters (no lists — pick one sweep value first).
    ``managed`` : workflow-supplied values for managed params (db, threads, io, paired...).
    bool flags are emitted only when truthy; empty strings are skipped.

    A ``managed`` key that is not a registry param raises: this is the workflow-owned
    contract between a rule and the registry, so a typo (``reprot``) or a renamed param
    must fail loudly here rather than silently drop the flag and corrupt the command.
    """
    params = _params(tool)
    managed = managed or {}
    unknown = [k for k in managed if k not in params]
    if unknown:
        known = sorted(k for k, s in params.items() if s.get("managed"))
        raise ValidationError(
            f"render_args({tool!r}): managed key(s) {unknown} are not parameters of the "
            f"{tool} registry (typo, or the registry renamed the param). "
            f"Registry-declared managed params: {known}."
        )
    merged: Dict[str, Any] = {}
    merged.update(values or {})
    merged.update(managed)

    args: List[str] = []
    trailing: List[str] = []  # passthrough tokens, appended verbatim after all flags
    for name, value in merged.items():
        if name not in params:
            continue
        spec = params[name]
        # `interpreted` params are user-facing (asked + validated) but consumed by a
        # workflow script, not emitted as a CLI flag — never render them here.
        if spec.get("interpreted"):
            continue
        if value is None or value == "":
            continue
        # `passthrough` params are the raw escape hatch: the value is a string of
        # one or more whole CLI tokens, split shell-style and appended verbatim with
        # no flag of their own. They cannot be bounds-checked — provenance logs them.
        if spec.get("passthrough"):
            trailing.extend(shlex.split(str(value)))
            continue
        flag = spec.get("flag")
        if spec["type"] == "bool":
            if bool(value) and flag:
                args.append(flag)
        elif flag:
            args.extend([flag, str(value)])
    return args + trailing
