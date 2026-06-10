"""Post-run advisor: metrics → multi-tool recommendations → next config hints."""

from __future__ import annotations

import copy
import glob
import json
import os
from typing import Any, Dict, List, Optional

import yaml

from . import evidence_pack, tool_advisor
from .report import classification_metrics, _outdir  # noqa: PLC2701


def platforms_from_config(cfg: Dict[str, Any]) -> List[str]:
    return sorted({c["platform"] for c in tool_advisor.sample_contexts(cfg)})


def _mean_classified(metrics: Dict[str, Any]) -> Optional[float]:
    vals = [m.get("percent_classified") for m in metrics.values() if m.get("percent_classified") is not None]
    return round(sum(vals) / len(vals), 3) if vals else None


def _worst_classified(metrics: Dict[str, Any]) -> Optional[float]:
    vals = [m.get("percent_classified") for m in metrics.values() if m.get("percent_classified") is not None]
    return min(vals) if vals else None


def analyze(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Inspect finished results and return advisor payload (rules-first, all active tools)."""
    outdir = _outdir(cfg)
    plats = platforms_from_config(cfg)
    primary = plats[0] if plats else "illumina"

    class_metrics = classification_metrics(outdir) if os.path.isdir(outdir) else {}
    kraken_cfg = cfg.get("kraken2") or {}
    sweep = cfg.get("sweep") or {}

    multi = tool_advisor.recommend_config(cfg)
    warnings: List[str] = list(multi.get("warnings") or [])
    suggestions: List[str] = list(multi.get("suggestions") or [])
    config_patches: Dict[str, Any] = dict(multi.get("config_patches") or {})
    suggested_sweep: Optional[List[float]] = None

    for plat in plats:
        warnings.extend(evidence_pack.registry_warnings("kraken2", plat, kraken_cfg))

    mean_pc = _mean_classified(class_metrics)
    worst_pc = _worst_classified(class_metrics)

    if primary in ("pacbio_clr", "pacbio_hifi") and worst_pc is not None and worst_pc < 50:
        rec = evidence_pack.recommend("kraken2", primary)
        suggested_sweep = rec.get("sweep_suggest")
        suggestions.append(
            f"Low classification ({worst_pc}%) on {primary}; try confidence sweep "
            f"{suggested_sweep} (evidence-based)."
        )
        if not sweep and rec.get("sweep_config"):
            config_patches["sweep"] = rec["sweep_config"]

    if primary == "illumina" and worst_pc is not None and worst_pc < 90:
        conf = kraken_cfg.get("confidence")
        if conf is not None and float(conf) >= 0.5:
            suggestions.append("Lower kraken2 confidence below 0.5 to recover classified reads.")
            config_patches.setdefault("kraken2", {})["confidence"] = 0.2

    for opt in multi.get("optional_modules") or []:
        if opt.get("ready"):
            suggestions.append(
                f"Consider enabling modules.{opt['module']}: {opt.get('reason', '')}"
            )

    matrix_files = glob.glob(os.path.join(outdir, "summary", "*.matrix.json"))
    if matrix_files and not sweep:
        suggestions.append(
            "Consider a confidence sweep — run `metagx recommend --config <yaml>` for "
            "platform-aware grids."
        )

    verdict = "ok"
    if worst_pc is not None and worst_pc < 30:
        verdict = "poor"
    elif worst_pc is not None and worst_pc < 70:
        verdict = "marginal"
    if not suggestions and mean_pc is not None and mean_pc >= 90:
        verdict = "good"

    return {
        "project": cfg.get("project"),
        "outdir": outdir,
        "platforms": plats,
        "metrics": {
            "classification": class_metrics,
            "mean_percent_classified": mean_pc,
            "worst_percent_classified": worst_pc,
        },
        "recommendations": multi,
        "warnings": list(dict.fromkeys(w for w in warnings if w)),
        "suggestions": list(dict.fromkeys(suggestions)),
        "suggested_sweep": suggested_sweep,
        "config_patches": config_patches,
        "verdict": verdict,
    }


def write_advisor_outputs(
    cfg: Dict[str, Any],
    analysis: Dict[str, Any],
    advisor_dir: str | None = None,
) -> Dict[str, str]:
    """Write advisor JSON + optional patched config YAML for the next trial."""
    outdir = _outdir(cfg)
    adir = advisor_dir or os.path.join(outdir, "advisor")
    os.makedirs(adir, exist_ok=True)

    paths: Dict[str, str] = {}
    report_path = os.path.join(adir, "advisor.json")
    with open(report_path, "w") as fh:
        json.dump(analysis, fh, indent=2)
    paths["advisor_json"] = report_path

    if analysis.get("config_patches"):
        merged = copy.deepcopy(cfg)
        _deep_merge(merged, analysis["config_patches"])
        next_path = os.path.join(adir, "next_config.suggested.yaml")
        with open(next_path, "w") as fh:
            yaml.dump(merged, fh, default_flow_style=False, sort_keys=False)
        paths["next_config"] = next_path

    log_path = os.path.join(adir, "trial_log.md")
    with open(log_path, "w") as fh:
        fh.write(f"# Advisor — {cfg.get('project', 'run')}\n\n")
        fh.write(f"**Verdict:** {analysis.get('verdict')}\n\n")
        rec = analysis.get("recommendations") or {}
        if rec.get("qc_routing"):
            fh.write("## QC routing\n\n")
            for route in rec["qc_routing"]:
                fh.write(f"- **{route.get('platform')}**: {route.get('pipeline') or route.get('default')}\n")
            fh.write("\n")
        if analysis.get("warnings"):
            fh.write("## Warnings\n\n")
            for w in analysis["warnings"]:
                fh.write(f"- {w}\n")
            fh.write("\n")
        if analysis.get("suggestions"):
            fh.write("## Suggestions\n\n")
            for s in analysis["suggestions"]:
                fh.write(f"- {s}\n")
            fh.write("\n")
        m = analysis.get("metrics") or {}
        if m.get("mean_percent_classified") is not None:
            fh.write(f"Mean % classified: {m['mean_percent_classified']}\n")
    paths["trial_log"] = log_path
    return paths


def _deep_merge(base: Dict[str, Any], patch: Dict[str, Any]) -> None:
    for k, v in patch.items():
        if k in base and isinstance(base[k], dict) and isinstance(v, dict):
            _deep_merge(base[k], v)
        else:
            base[k] = v
