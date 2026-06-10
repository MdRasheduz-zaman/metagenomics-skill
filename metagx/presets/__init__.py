"""Workflow presets: named partial-config templates the interview can offer.

A preset supplies module toggles and sensible tool parameters; the user always supplies
samples + db, and any answer they give overrides the preset. Presets are validated through
the same registries as everything else, so they cannot drift from the real flags.

The ``*.yaml`` files in this package directory are the preset definitions.
"""

from __future__ import annotations

import functools
from importlib import resources
from typing import Any, Dict, List

import yaml

_PKG = "metagx.presets"


@functools.lru_cache(maxsize=None)
def _load(name: str) -> Dict[str, Any]:
    path = resources.files(_PKG) / f"{name}.yaml"
    if not path.is_file():
        raise KeyError(f"Unknown preset '{name}'. Available: {', '.join(list_presets())}")
    return yaml.safe_load(path.read_text())


def list_presets() -> List[str]:
    pkg = resources.files(_PKG)
    return sorted(p.name[:-5] for p in pkg.iterdir() if p.name.endswith(".yaml"))


def describe_presets() -> List[Dict[str, str]]:
    """Name / description / when_to_use for each preset — for the interview to present."""
    out = []
    for name in list_presets():
        p = _load(name)
        out.append(
            {
                "name": p.get("name", name),
                "description": " ".join(str(p.get("description", "")).split()),
                "when_to_use": " ".join(str(p.get("when_to_use", "")).split()),
            }
        )
    return out


def get_preset_config(name: str) -> Dict[str, Any]:
    """The partial config dict a preset contributes (no samples/db)."""
    return dict(_load(name).get("config", {}))


def deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge ``override`` onto ``base``. Dicts merge; scalars/lists replace."""
    out = dict(base)
    for key, val in (override or {}).items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], val)
        else:
            out[key] = val
    return out
