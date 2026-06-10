"""Load bundled evidence YAML and merge with registry guidance for recommendations."""

from __future__ import annotations

import functools
from importlib import resources
from typing import Any, Dict, List, Optional

import yaml

from . import registry

PLATFORM_ALIASES = {
    "illumina": "illumina",
    "mgi": "mgi",
    "bgi": "bgi",
    "ont": "ont",
    "nanopore": "ont",
    "pacbio_hifi": "pacbio_hifi",
    "hifi": "pacbio_hifi",
    "pacbio_clr": "pacbio_clr",
    "clr": "pacbio_clr",
    "pacbio": "pacbio_clr",
}


def normalize_platform(platform: str) -> str:
    key = str(platform or "illumina").strip().lower()
    return PLATFORM_ALIASES.get(key, key)


@functools.lru_cache(maxsize=None)
def list_evidence() -> List[str]:
    pkg = resources.files("metagx.evidence")
    return sorted(
        p.name[:-5]
        for p in pkg.iterdir()
        if p.name.endswith(".yaml") and p.name != "README.md"
    )


def load_evidence(name: str) -> Dict[str, Any]:
    pkg = resources.files("metagx.evidence")
    path = pkg / f"{name}.yaml"
    if not path.is_file():
        raise KeyError(f"Unknown evidence '{name}'. Available: {', '.join(list_evidence())}")
    return yaml.safe_load(path.read_text())


EVIDENCE_BY_TOOL_PARAM = {
    ("kraken2", "confidence"): "kraken2_confidence",
    ("bracken", "read_length"): "bracken_read_length",
    ("cutadapt", "minimum_length"): "cutadapt_amplicon",
    ("metaspades", "memory_gb"): "metaspades_memory",
    ("kaiju", "min_match_length"): "kaiju_consensus",
    ("kaiju", "mismatches"): "kaiju_consensus",
    ("metaphlan", "stat_q"): "metaphlan_consensus",
    ("mafft", "method"): "phylogenetics_mafft",
    ("iqtree", "bootstrap"): "phylogenetics_iqtree",
}


def _registry_recommend(tool: str, param: str, platform: str) -> Any:
    reg = registry.load_registry(tool)
    spec = reg.get("params", {}).get(param, {})
    rec = spec.get("recommend") or {}
    plat = normalize_platform(platform)
    vals = rec.get(plat)
    if vals is None:
        vals = rec.get("default")
    if vals is None:
        return None
    if isinstance(vals, list):
        return [float(v) if isinstance(v, (int, float)) else v for v in vals]
    return vals


def _evidence_file_for(tool: str, param: str, evidence_name: str | None) -> str | None:
    if evidence_name:
        return evidence_name
    return EVIDENCE_BY_TOOL_PARAM.get((tool, param))


def recommend(
    tool: str,
    platform: str,
    param: str = "confidence",
    evidence_name: str | None = None,
) -> Dict[str, Any]:
    """Merge registry ``recommend`` blocks with bundled evidence for a platform."""
    plat = normalize_platform(platform)
    reg_val = _registry_recommend(tool, param, plat)
    out: Dict[str, Any] = {
        "tool": tool,
        "platform": plat,
        "param": param,
        "registry_value": reg_val,
        "evidence_value": None,
        "sweep_suggest": None,
        "value_suggest": None,
        "warnings": [],
        "notes": [],
    }

    ev_file = _evidence_file_for(tool, param, evidence_name)
    if ev_file:
        ev = load_evidence(ev_file)
        plat_ev = (ev.get("platforms") or {}).get(plat, {})
        if not plat_ev:
            plat_ev = (ev.get("platforms") or {}).get("default", {})
        if plat_ev:
            if param == "read_length":
                out["evidence_value"] = plat_ev.get("default")
                out["value_suggest"] = plat_ev.get("default")
                if plat_ev.get("notes"):
                    out["notes"].append(str(plat_ev["notes"]).strip())
                out["build_db_lengths"] = plat_ev.get("build_db_lengths")
            elif param in plat_ev:
                out["evidence_value"] = plat_ev[param]
                out["value_suggest"] = plat_ev[param]
                if plat_ev.get("notes"):
                    out["notes"].append(str(plat_ev["notes"]).strip())
            else:
                out["evidence_value"] = plat_ev.get("sweep_default") or plat_ev.get("sweep_recommend")
                if plat_ev.get("default") is not None and out["evidence_value"] is None:
                    out["evidence_value"] = plat_ev.get("default")
                    out["value_suggest"] = plat_ev.get("default")
                if plat_ev.get("notes"):
                    out["notes"].append(str(plat_ev["notes"]).strip())
                if plat_ev.get("avoid_above") is not None:
                    out["warnings"].append(
                        f"Avoid {param} above {plat_ev['avoid_above']} on {plat} (validation runs)."
                    )
                out["observations"] = plat_ev.get("observations")

    reg = registry.load_registry(tool)
    for rule in (reg.get("params", {}).get(param, {}).get("warn_if") or []):
        if _warn_rule_matches(rule, plat, {}):
            out["warnings"].append(rule.get("message", "").strip())

    if isinstance(out["evidence_value"], list) or isinstance(reg_val, list):
        out["sweep_suggest"] = (
            out["evidence_value"] if isinstance(out["evidence_value"], list) else None
        ) or (reg_val if isinstance(reg_val, list) else None) or [0.0, 0.1]
        out["sweep_config"] = {"param": param, "values": out["sweep_suggest"]}
    else:
        out["value_suggest"] = out["evidence_value"] if out["evidence_value"] is not None else reg_val
    return out


def _warn_rule_matches(rule: Dict[str, Any], platform: str, ctx: Dict[str, Any]) -> bool:
    when = rule.get("when") or {}
    if "platform" in when and normalize_platform(when["platform"]) != platform:
        return False
    if "confidence_gte" in when or "confidence_gt" in when:
        conf = ctx.get("confidence")
        if conf is None:
            return False
        if "confidence_gte" in when and float(conf) < float(when["confidence_gte"]):
            return False
        if "confidence_gt" in when and float(conf) <= float(when["confidence_gt"]):
            return False
    return True


def registry_warnings(tool: str, platform: str, param_values: Dict[str, Any]) -> List[str]:
    """Evaluate registry ``warn_if`` rules for a tool given user param values."""
    plat = normalize_platform(platform)
    ctx = dict(param_values or {})
    msgs: List[str] = []
    reg = registry.load_registry(tool)
    for pname, pspec in (reg.get("params") or {}).items():
        if pname in ctx:
            ctx_for_param = {**ctx, "confidence": ctx.get(pname) if pname == "confidence" else ctx.get("confidence")}
        else:
            ctx_for_param = ctx
        for rule in pspec.get("warn_if") or []:
            if _warn_rule_matches(rule, plat, ctx_for_param):
                msg = (rule.get("message") or "").strip()
                if msg:
                    msgs.append(msg)
    return msgs
