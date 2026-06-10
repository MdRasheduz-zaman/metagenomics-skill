"""Append-only run trial log for personalization and advisor context."""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
from typing import Any, Dict, List, Optional

import yaml

DEFAULT_HISTORY = ".metagx/history.jsonl"


def history_path(path: str | None = None) -> str:
    return os.path.abspath(path or DEFAULT_HISTORY)


def _config_hash(cfg: Dict[str, Any]) -> str:
    blob = yaml.dump(cfg, sort_keys=True, default_flow_style=False)
    return hashlib.sha256(blob.encode()).hexdigest()[:16]


def append_entry(entry: Dict[str, Any], path: str | None = None) -> str:
    """Append one JSON line; returns the file path."""
    fp = history_path(path)
    os.makedirs(os.path.dirname(fp) or ".", exist_ok=True)
    with open(fp, "a") as fh:
        fh.write(json.dumps(entry, default=str) + "\n")
    return fp


def read_entries(path: str | None = None, limit: int = 50) -> List[Dict[str, Any]]:
    fp = history_path(path)
    if not os.path.isfile(fp):
        return []
    rows: List[Dict[str, Any]] = []
    with open(fp) as fh:
        for line in fh:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows[-limit:]


def record_from_run(
    cfg: Dict[str, Any],
    config_path: str,
    analysis: Dict[str, Any],
    *,
    success: bool,
    returncode: int = 0,
    path: str | None = None,
    verdict: str | None = None,
) -> Dict[str, Any]:
    """Build and append a history row after a pipeline run or advise pass."""
    metrics = analysis.get("metrics") or {}
    entry = {
        "ts": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "project": cfg.get("project"),
        "config_path": os.path.abspath(config_path),
        "config_hash": _config_hash(cfg),
        "platforms": analysis.get("platforms") or _platforms_from_config(cfg),
        "success": success,
        "returncode": returncode,
        "verdict": verdict or analysis.get("verdict"),
        "metrics": metrics,
        "suggestions": analysis.get("suggestions", []),
        "warnings": analysis.get("warnings", []),
        "best_sweep": analysis.get("suggested_sweep"),
    }
    append_entry(entry, path=path)
    return entry


def _platforms_from_config(cfg: Dict[str, Any]) -> List[str]:
    samples = cfg.get("samples")
    if isinstance(samples, list):
        return sorted({str(s.get("platform", "illumina")).lower() for s in samples})
    return ["illumina"]


def best_trial(
    path: str | None = None,
    metric: str = "mean_percent_classified",
) -> Optional[Dict[str, Any]]:
    """Return the trial with the highest value for a top-level metrics field."""
    rows = read_entries(path=path, limit=500)
    scored = []
    for r in rows:
        m = r.get("metrics") or {}
        val = m.get(metric)
        if val is not None:
            scored.append((float(val), r))
    if not scored:
        return None
    scored.sort(key=lambda x: x[0], reverse=True)
    return scored[0][1]
