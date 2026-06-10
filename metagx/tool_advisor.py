"""Platform-aware recommendations across all workflow tools (rules + evidence + routing)."""

from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Tuple

from . import evidence_pack, formats, registry

# Module flag -> primary config key / tool(s) when enabled
MODULE_TOOLS: Dict[str, List[str]] = {
    "qc": [],  # resolved per-sample via routing
    "classify": ["kraken2"],
    "classify_consensus": [],  # metaphlan or kaiju from consensus.classifier
    "abundance": ["bracken"],
    "assembly": [],  # megahit | flye | metaspades via assembler routing
    "binning": ["metabat2"],
    "bin_refinement": ["maxbin2", "concoct", "das_tool", "drep"],
    "strain": ["instrain"],
    "damage": ["mapdamage"],
    "reconcile": ["kraken2"],
    "domain_taxonomy": ["genomad", "checkv", "gtdbtk", "checkm2", "eukrep", "eukcc"],
    "stats": [],
    "differential": [],
    "decontam": [],
    "phylogenetics": ["mafft", "iqtree", "fasttree"],
    "bgc": ["antismash"],
    "aggregate": ["multiqc", "krona"],
    "functional": ["humann", "amrfinderplus", "abricate", "bakta", "eggnog"],
    "filtered_assembly": ["megahit", "flye"],
}

SHORT_PLATFORMS = {"illumina", "mgi", "bgi"}
LONG_PLATFORMS = {"ont", "pacbio_hifi", "pacbio_clr", "pacbio"}


def sample_contexts(cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Per-sample platform / library context from the sample sheet."""
    samples = cfg.get("samples") or []
    if not isinstance(samples, list):
        return [{"platform": "illumina", "library": "wgs", "layout": "se"}]
    out = []
    for s in samples:
        plat = evidence_pack.normalize_platform(s.get("platform", "illumina"))
        lib = str(s.get("library", "wgs")).lower()
        layout = str(s.get("layout", "se")).lower()
        if lib == "amplicon":
            qc_key = "amplicon"
        elif lib == "ancient" and layout == "pe":
            qc_key = "ancient_pe"
        elif plat in SHORT_PLATFORMS:
            qc_key = plat
        elif plat in LONG_PLATFORMS:
            qc_key = plat
        else:
            qc_key = plat
        out.append({
            "sample": s.get("sample"),
            "platform": plat,
            "library": lib,
            "layout": layout,
            "qc_key": qc_key,
            "reads": s.get("r1") or s.get("reads"),
            "bracken_read_length": s.get("bracken_read_length"),
        })
    return out


def _routing() -> Dict[str, Any]:
    return evidence_pack.load_evidence("platform_routing")


def qc_routing_for(ctx: Dict[str, Any]) -> Dict[str, Any]:
    qc = _routing().get("qc") or {}
    block = qc.get(ctx["qc_key"]) or qc.get(ctx["platform"]) or {}
    return {
        "platform": ctx["platform"],
        "qc_key": ctx["qc_key"],
        "default": block.get("default"),
        "pipeline": block.get("pipeline", []),
        "alternatives": block.get("alternatives", []),
        "note": block.get("note"),
    }


def assembly_routing_for(platform: str) -> Dict[str, Any]:
    asm = (_routing().get("assembly") or {}).get(platform, {})
    return {
        "platform": platform,
        "default": asm.get("default"),
        "wired_tools": asm.get("wired_tools", []),
        "alternatives": asm.get("alternatives", []),
    }


def _nearest_standard(length: int, standards: List[int]) -> int:
    if not standards:
        return length
    return min(standards, key=lambda x: abs(x - length))


def suggest_bracken_read_length(
    platform: str,
    median_read_len: Optional[int] = None,
    current: Optional[int] = None,
    library: str = "wgs",
) -> Dict[str, Any]:
    """Pick Bracken -r from platform defaults + optional median read length."""
    plat = evidence_pack.normalize_platform(platform)
    if library == "ancient":
        plat = "ancient"
    ev = evidence_pack.load_evidence("bracken_read_length")
    standards = ev.get("standard_lengths") or [50, 100, 150, 250, 500, 1000]
    plat_ev = (ev.get("platforms") or {}).get(plat, {})
    default = int(plat_ev.get("default", 150))
    suggested = default
    if median_read_len and median_read_len > 0:
        if plat in SHORT_PLATFORMS:
            suggested = _nearest_standard(median_read_len, [100, 150, 200, 250])
        else:
            suggested = _nearest_standard(median_read_len, standards)
    out: Dict[str, Any] = {
        "platform": plat,
        "suggested_read_length": suggested,
        "platform_default": default,
        "median_input": median_read_len,
        "current": current,
        "build_db_lengths": plat_ev.get("build_db_lengths", standards[:6]),
        "notes": (plat_ev.get("notes") or "").strip(),
    }
    if current is not None and int(current) != suggested:
        out["warning"] = (
            f"Configured Bracken read_length {current} differs from suggested {suggested} "
            f"for {plat}; abundance may be miscalibrated unless the DB was built for {current}."
        )
    return out


def _median_from_reads(paths: List[str]) -> Optional[int]:
    medians = []
    for p in paths:
        if p and os.path.isfile(p):
            stats = formats.estimate_read_length(p)
            if stats.get("median"):
                medians.append(int(stats["median"]))
    if not medians:
        return None
    medians.sort()
    return medians[len(medians) // 2]


def bracken_recommendations(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Per-platform Bracken read_length + level/threshold from registry/evidence."""
    modules = cfg.get("modules") or {}
    if not modules.get("abundance", True):
        return {}
    br_cfg = cfg.get("bracken") or {}
    by_plat = cfg.get("bracken_read_length_by_platform") or {}
    contexts = sample_contexts(cfg)
    per_platform: Dict[str, Any] = {}
    warnings: List[str] = []

    paths_by_plat: Dict[str, List[str]] = {}
    for ctx in contexts:
        paths_by_plat.setdefault(ctx["platform"], []).append(ctx.get("reads") or "")

    for plat in sorted({c["platform"] for c in contexts}):
        median = _median_from_reads(paths_by_plat.get(plat, []))
        current = by_plat.get(plat)
        if current is None:
            current = br_cfg.get("read_length")
        for ctx in contexts:
            if ctx["platform"] == plat and ctx.get("bracken_read_length"):
                current = ctx["bracken_read_length"]
                break
        rec = suggest_bracken_read_length(
            plat, median_read_len=median, current=current, library=contexts[0].get("library", "wgs")
        )
        per_platform[plat] = rec
        if rec.get("warning"):
            warnings.append(rec["warning"])

    level_rec = evidence_pack.recommend("bracken", contexts[0]["platform"], param="level")
    threshold_rec = evidence_pack.recommend("bracken", contexts[0]["platform"], param="threshold")

    return {
        "read_length_by_platform": per_platform,
        "level": level_rec,
        "threshold": threshold_rec,
        "warnings": warnings,
        "config_patch": {
            "bracken_read_length_by_platform": {
                p: r["suggested_read_length"] for p, r in per_platform.items()
            }
        },
    }


def kraken2_secondary_recommendations(cfg: Dict[str, Any], platform: str) -> Dict[str, Any]:
    ev = evidence_pack.load_evidence("kraken2_secondary")
    k_cfg = cfg.get("kraken2") or {}
    plat = evidence_pack.normalize_platform(platform)
    out: Dict[str, Any] = {"platform": plat, "params": {}, "warnings": []}
    for param in ("minimum_hit_groups", "minimum_base_quality", "quick"):
        block = ev.get(param) or {}
        suggested = block.get(plat)
        if suggested is None:
            suggested = block.get("default")
        if suggested is None and param == "minimum_hit_groups":
            suggested = block.get("pathogen_screening") if "pathogen" in str(cfg.get("preset", "")) else 2
        current = k_cfg.get(param)
        out["params"][param] = {
            "suggested": suggested,
            "current": current,
            "notes": block.get("notes"),
        }
        if current is not None and suggested is not None and current != suggested:
            if param == "minimum_hit_groups" and int(current) < int(suggested):
                out["warnings"].append(
                    f"kraken2 {param}={current} is below suggested {suggested} for {plat}."
                )
    return out


def active_module_tools(cfg: Dict[str, Any]) -> List[Tuple[str, str]]:
    """Return (module, tool) pairs for enabled modules."""
    modules = cfg.get("modules") or {}
    pairs: List[Tuple[str, str]] = []
    for mod, tools in MODULE_TOOLS.items():
        if not _module_enabled(modules, mod):
            continue
        if mod == "classify_consensus":
            clf = (cfg.get("consensus") or {}).get("classifier", "metaphlan")
            pairs.append((mod, clf))
        elif mod == "assembly":
            for plat in {c["platform"] for c in sample_contexts(cfg)}:
                asm = assembly_routing_for(plat)
                default = asm.get("default")
                if default:
                    pairs.append((mod, default))
        elif mod == "qc":
            for ctx in sample_contexts(cfg):
                route = qc_routing_for(ctx)
                for t in route.get("pipeline") or [route.get("default")]:
                    if t:
                        pairs.append((mod, t))
        else:
            for t in tools:
                pairs.append((mod, t))
    # dedupe
    seen = set()
    unique = []
    for p in pairs:
        if p not in seen:
            seen.add(p)
            unique.append(p)
    return unique


def params_with_guidance(tool: str) -> List[str]:
    """Parameter names that have recommend or warn_if in the registry."""
    reg = registry.load_registry(tool)
    names = []
    for pname, pspec in (reg.get("params") or {}).items():
        if pspec.get("recommend") or pspec.get("warn_if"):
            names.append(pname)
    return names


def recommend_tool_params(
    tool: str,
    platform: str,
    cfg: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    """All parameter recommendations for one tool on one platform."""
    results = []
    ctx = (cfg or {}).get(tool) or {}
    for param in params_with_guidance(tool) or _default_params_for_tool(tool):
        rec = evidence_pack.recommend(tool, platform, param=param)
        rec["registry_warnings"] = evidence_pack.registry_warnings(tool, platform, ctx)
        results.append(rec)
    if tool == "kraken2" and cfg:
        sec = kraken2_secondary_recommendations(cfg, platform)
        for pname, pdata in sec.get("params", {}).items():
            results.append({
                "tool": "kraken2",
                "platform": platform,
                "param": pname,
                "value_suggest": pdata.get("suggested"),
                "current": pdata.get("current"),
                "source": "kraken2_secondary.yaml",
            })
    return results


def _default_params_for_tool(tool: str) -> List[str]:
    defaults = {
        "bracken": ["read_length", "level", "threshold"],
        "kraken2": ["confidence"],
        "fastp": ["qualified_quality_phred", "length_required"],
        "chopper": ["quality", "minlength"],
        "porechop_abi": ["ab_initio", "extra_end_trim"],
        "flye": ["meta"],
        "megahit": ["presets"],
        "metaspades": ["memory_gb"],
        "metabat2": ["min_contig"],
        "cutadapt": ["minimum_length"],
        "kaiju": ["min_match_length", "mismatches"],
        "metaphlan": ["stat_q"],
        "mafft": ["method"],
        "iqtree": ["bootstrap", "model"],
        "fasttree": ["model", "sequence_type"],
    }
    return defaults.get(tool, [])


def _module_enabled(modules: Dict[str, Any], name: str) -> bool:
    defaults = {"qc": True, "classify": True, "abundance": True}
    return bool(modules.get(name, defaults.get(name, False)))


def optional_modules(cfg: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Suggest modules/tools not yet enabled."""
    modules = cfg.get("modules") or {}
    samples = cfg.get("samples") or []
    n = len(samples) if isinstance(samples, list) else 0
    suggestions = []
    for entry in _routing().get("optional_modules") or []:
        when = entry.get("when") or {}
        needs = when.get("needs") or []
        mod = entry["module"]
        if _module_enabled(modules, mod):
            continue
        prereqs_met = all(_module_enabled(modules, m) for m in needs) if needs else True
        if when.get("min_samples") and n < int(when["min_samples"]):
            prereqs_met = False
        if not prereqs_met:
            suggestions.append({
                "module": mod,
                "tools": entry.get("tools", []),
                "reason": entry.get("reason"),
                "ready": False,
                "prerequisite": f"Enable modules.{', modules.'.join(needs)} first" if needs else None,
            })
        else:
            suggestions.append({
                "module": mod,
                "tools": entry.get("tools", []),
                "reason": entry.get("reason"),
                "ready": True,
                "config_patch": {"modules": {mod: True}},
            })
    return suggestions


def classify_alternatives() -> List[Dict[str, Any]]:
    block = _routing().get("classify") or {}
    return block.get("alternatives", [])


def recommend_config(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Full multi-tool recommendation payload for a config (interview / advise / CLI)."""
    contexts = sample_contexts(cfg)
    platforms = sorted({c["platform"] for c in contexts})
    primary = platforms[0] if platforms else "illumina"

    qc_routes = [qc_routing_for(c) for c in contexts]
    asm_routes = {p: assembly_routing_for(p) for p in platforms}

    tool_params: Dict[str, List[Dict[str, Any]]] = {}
    for _mod, tool in active_module_tools(cfg):
        for plat in platforms:
            key = f"{tool}@{plat}"
            tool_params[key] = recommend_tool_params(tool, plat, cfg)

    bracken = bracken_recommendations(cfg) if (cfg.get("modules") or {}).get("abundance", True) else {}

    warnings: List[str] = []
    suggestions: List[str] = []
    config_patches: Dict[str, Any] = {}

    for route in qc_routes:
        for alt in route.get("alternatives") or []:
            if not alt.get("wired"):
                suggestions.append(
                    f"QC alternative for {route['platform']}: {alt['tool']} — {alt.get('reason', '')}"
                )

    if bracken.get("warnings"):
        warnings.extend(bracken["warnings"])
    if bracken.get("config_patch"):
        existing = cfg.get("bracken_read_length_by_platform") or {}
        patch = bracken["config_patch"]["bracken_read_length_by_platform"]
        if any(existing.get(p) != v for p, v in patch.items()):
            config_patches.setdefault("bracken_read_length_by_platform", {}).update(patch)
            suggestions.append(
                f"Set bracken_read_length_by_platform from medians: {patch}"
            )

    conf_rec = evidence_pack.recommend("kraken2", primary, param="confidence")
    if conf_rec.get("warnings"):
        warnings.extend(conf_rec["warnings"])

    return {
        "platforms": platforms,
        "primary_platform": primary,
        "qc_routing": qc_routes,
        "assembly_routing": asm_routes,
        "tool_parameters": tool_params,
        "bracken": bracken,
        "classify_alternatives": classify_alternatives(),
        "optional_modules": optional_modules(cfg),
        "warnings": list(dict.fromkeys(w for w in warnings if w)),
        "suggestions": suggestions,
        "config_patches": config_patches,
        "confidence_sweep": conf_rec.get("sweep_config"),
    }
