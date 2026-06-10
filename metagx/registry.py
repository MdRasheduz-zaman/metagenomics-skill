"""Load the per-tool parameter registries and derive everything from them.

This module is intentionally dependency-light (only PyYAML) so it can be imported
both inside the Snakemake workflow and by the MCP server / CLI.
"""

from __future__ import annotations

import functools
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


def _params(tool: str) -> Dict[str, Any]:
    return load_registry(tool)["params"]


# --------------------------------------------------------------------------- #
# Interview                                                                    #
# --------------------------------------------------------------------------- #
def interview_spec(tool: str, max_tier: int = 2) -> List[Dict[str, Any]]:
    """User-facing parameters an LLM should consider asking about.

    Returns only params with ``ask: true`` at or below ``max_tier`` (1=core,
    2=common, 3=advanced). Managed params (db/threads/io) are never returned.
    Each entry carries everything the LLM needs to phrase a good question and
    validate the answer.
    """
    out: List[Dict[str, Any]] = []
    for name, spec in _params(tool).items():
        if spec.get("managed"):
            continue
        if not spec.get("ask", False):
            continue
        if spec.get("tier", 3) > max_tier:
            continue
        out.append(
            {
                "name": name,
                "type": spec["type"],
                "default": spec.get("default"),
                "min": spec.get("min"),
                "max": spec.get("max"),
                "choices": spec.get("choices"),
                "tier": spec.get("tier", 3),
                "sweepable": bool(spec.get("sweepable", False)),
                "question": " ".join(str(spec.get("question", "")).split()),
            }
        )
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
    """
    params = _params(tool)
    merged: Dict[str, Any] = {}
    merged.update(values or {})
    merged.update(managed or {})

    args: List[str] = []
    for name, value in merged.items():
        if name not in params:
            continue
        spec = params[name]
        flag = spec.get("flag")
        if value is None or value == "":
            continue
        if spec["type"] == "bool":
            if bool(value) and flag:
                args.append(flag)
        elif flag:
            args.extend([flag, str(value)])
    return args
