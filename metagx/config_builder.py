"""Turn interview answers into a validated config.yaml for the Snakemake workflow.

The MCP server, the CLI, and any paste-in LLM all converge here: the LLM gathers
answers in natural language, then calls ``build_config`` which validates every section
against the registries and writes a deterministic config file.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List

import yaml

from . import presets, registry

MODULE_TOOLS = {
    "qc": "fastp",
    "classify": "kraken2",
    "abundance": "bracken",
    "assembly": "megahit",
    "binning": "metabat2",
    "bin_refinement": "das_tool",
    "classify_consensus": "metaphlan",
    "functional": "humann",
    "aggregate": "multiqc",
    "damage": "mapdamage",
    "strain": "instrain",
    "bgc": "antismash",
}

DEFAULT_MODULES = {
    "qc": True,
    "classify": True,
    "abundance": True,
    "assembly": False,
    "binning": False,
    "bin_refinement": False,
    "classify_consensus": False,
    "reconcile": False,
    "domain_taxonomy": False,
    "filtered_assembly": False,
    "stats": False,
    "differential": False,
    "functional": False,
    "bgc": False,
    "aggregate": False,
    "damage": False,
    "decontam": False,
    "strain": False,
}

KNOWN_DOMAINS = {"viral", "prokaryote", "eukaryote"}


KNOWN_PLATFORMS = {"illumina", "mgi", "bgi", "ont", "nanopore",
                   "pacbio_hifi", "pacbio_clr", "pacbio"}
LONG_PLATFORMS = {"ont", "nanopore", "pacbio_hifi", "pacbio_clr", "pacbio"}


def _validate_samples(samples: Any) -> Any:
    """Accept a path to a TSV sample sheet or an inline list of records.

    Inline records may carry optional 'platform' (illumina/mgi/ont/pacbio_hifi/pacbio_clr)
    and 'layout' (se/pe/interleaved). Combinations are checked so bad sheets fail clearly.
    (TSV paths are validated by the workflow at runtime.)
    """
    if isinstance(samples, str):
        return samples
    if not isinstance(samples, list):
        raise registry.ValidationError(
            "samples must be a path to a TSV or a list of {sample, r1, r2?, platform?, layout?}"
        )
    for i, rec in enumerate(samples):
        if "sample" not in rec or "r1" not in rec:
            raise registry.ValidationError(f"sample #{i} needs at least 'sample' and 'r1'")
        plat = str(rec.get("platform", "illumina")).lower()
        if plat not in KNOWN_PLATFORMS:
            raise registry.ValidationError(
                f"sample '{rec['sample']}': unknown platform '{plat}'"
            )
        lay = str(rec.get("layout", "")).lower()
        if lay and lay not in {"se", "pe", "interleaved"}:
            raise registry.ValidationError(
                f"sample '{rec['sample']}': layout must be se|pe|interleaved"
            )
        if plat in LONG_PLATFORMS and (rec.get("r2") or lay in {"pe", "interleaved"}):
            raise registry.ValidationError(
                f"sample '{rec['sample']}': long-read platforms are single-end (no r2, layout=se)"
            )
        if lay == "interleaved" and plat not in {"illumina", "mgi", "bgi"}:
            raise registry.ValidationError(
                f"sample '{rec['sample']}': interleaved layout is for short reads only"
            )
        lib = str(rec.get("library", "wgs")).lower()
        if lib not in {"wgs", "amplicon", "ancient"}:
            raise registry.ValidationError(
                f"sample '{rec['sample']}': library must be 'wgs', 'amplicon', or 'ancient'"
            )
        if lib == "ancient" and plat not in {"illumina", "mgi", "bgi"}:
            raise registry.ValidationError(
                f"sample '{rec['sample']}': ancient library is short-read shotgun; "
                f"platform must be illumina/mgi/bgi, got '{plat}'"
            )
        if rec.get("bracken_read_length") is not None:
            try:
                int(rec["bracken_read_length"])
            except (TypeError, ValueError):
                raise registry.ValidationError(
                    f"sample '{rec['sample']}': bracken_read_length must be an integer"
                )
        if rec.get("group") is not None and not str(rec["group"]).strip():
            raise registry.ValidationError(
                f"sample '{rec['sample']}': group, if given, must be a non-empty label"
            )
        if rec.get("long_reads"):
            lp = str(rec.get("long_platform", "ont")).lower()
            if lp not in LONG_PLATFORMS:
                raise registry.ValidationError(
                    f"sample '{rec['sample']}': long_platform must be a long-read platform "
                    f"({sorted(LONG_PLATFORMS)})"
                )
            if plat in LONG_PLATFORMS:
                raise registry.ValidationError(
                    f"sample '{rec['sample']}': long_reads is for hybrid assembly of a "
                    "short-read sample; this sample is already long-read"
                )
    return samples


def _all_amplicon(samples: Any) -> bool:
    """True only if samples is an inline list and every record is amplicon."""
    return (isinstance(samples, list) and len(samples) > 0
            and all(str(r.get("library", "wgs")).lower() == "amplicon" for r in samples))


def _validate_sweep(sweep: Dict[str, Any] | None, kraken2: Dict[str, Any]) -> Dict[str, Any] | None:
    if not sweep:
        return None
    param = sweep.get("param", "confidence")
    values = sweep.get("values")
    if not isinstance(values, list) or len(values) < 1:
        raise registry.ValidationError("sweep.values must be a non-empty list")
    spec = registry.load_registry("kraken2")["params"].get(param)
    if not spec or not spec.get("sweepable"):
        raise registry.ValidationError(f"kraken2 param '{param}' is not sweepable")
    # validate each sweep value through the registry
    registry.validate("kraken2", {param: values})
    # the swept param must not also be pinned in the base kraken2 section
    if param in kraken2:
        raise registry.ValidationError(
            f"'{param}' is both swept and pinned; remove it from the kraken2 section"
        )
    return {"param": param, "values": values}


def _validate_read_filter(rf: Dict[str, Any] | None) -> Dict[str, Any] | None:
    if not rf:
        return None
    mode = str(rf.get("mode", "exclude")).lower()
    if mode not in {"include", "exclude"}:
        raise registry.ValidationError("read_filter.mode must be 'include' or 'exclude'")
    taxids = rf.get("taxids", [])
    if not isinstance(taxids, list) or not all(isinstance(t, int) for t in taxids):
        raise registry.ValidationError("read_filter.taxids must be a list of integer taxids")
    if mode == "exclude" and not (taxids or rf.get("host_genome")):
        raise registry.ValidationError(
            "exclude mode needs taxids and/or host_genome to remove something")
    out = {
        "mode": mode,
        "taxids": taxids,
        "include_children": bool(rf.get("include_children", True)),
        "keep_unclassified": bool(rf.get("keep_unclassified", True)),
    }
    if rf.get("host_genome"):
        out["host_genome"] = rf["host_genome"]
    return out


def _validate_subsample(sub: Dict[str, Any] | None) -> Dict[str, Any] | None:
    if not sub:
        return None
    try:
        fraction = float(sub.get("fraction", 1.0))
    except (TypeError, ValueError):
        raise registry.ValidationError("subsample.fraction must be a number")
    if not (0 < fraction <= 1):
        raise registry.ValidationError("subsample.fraction must be in (0, 1]")
    seed = int(sub.get("seed", 42))
    return {"fraction": fraction, "seed": seed}


def _validate_differential(diff: Dict[str, Any] | None) -> Dict[str, Any]:
    """Tuning block for the differential-abundance module (CLR + permutation test)."""
    diff = diff or {}
    try:
        n_perm = int(diff.get("n_permutations", 999))
        fdr = float(diff.get("fdr", 0.05))
        min_count = int(diff.get("min_count", 1))
    except (TypeError, ValueError):
        raise registry.ValidationError(
            "differential: n_permutations/min_count must be ints and fdr a number")
    if n_perm < 99:
        raise registry.ValidationError("differential.n_permutations must be >= 99")
    if not (0 < fdr < 1):
        raise registry.ValidationError("differential.fdr must be in (0, 1)")
    out = {"group_column": str(diff.get("group_column", "group")),
           "n_permutations": n_perm, "fdr": fdr, "min_count": min_count,
           "seed": int(diff.get("seed", 42))}
    if diff.get("reference_group"):
        out["reference_group"] = str(diff["reference_group"])
    return out


def _sample_groups(samples: Any, group_column: str = "group") -> Dict[str, int]:
    """Count samples per group label (inline lists only; TSV checked at runtime)."""
    if not isinstance(samples, list):
        return {}
    counts: Dict[str, int] = {}
    for rec in samples:
        g = rec.get(group_column) or rec.get("group")
        if g:
            counts[str(g)] = counts.get(str(g), 0) + 1
    return counts


def build_config(
    *,
    project: str = "run",
    outdir: str = "results",
    threads: int = 8,
    samples: Any,
    db: Dict[str, str],
    preset: str | None = None,
    modules: Dict[str, bool] | None = None,
    sweep: Dict[str, Any] | None = None,
    subsample: Dict[str, Any] | None = None,
    read_filter: Dict[str, Any] | None = None,
    host_removal: Dict[str, Any] | None = None,
    fastp: Dict[str, Any] | None = None,
    kraken2: Dict[str, Any] | None = None,
    bracken: Dict[str, Any] | None = None,
    megahit: Dict[str, Any] | None = None,
    metabat2: Dict[str, Any] | None = None,
    porechop_abi: Dict[str, Any] | None = None,
    chopper: Dict[str, Any] | None = None,
    flye: Dict[str, Any] | None = None,
    domains: List[str] | None = None,
    genomad: Dict[str, Any] | None = None,
    checkv: Dict[str, Any] | None = None,
    checkm2: Dict[str, Any] | None = None,
    gtdbtk: Dict[str, Any] | None = None,
    eukrep: Dict[str, Any] | None = None,
    eukcc: Dict[str, Any] | None = None,
    amplicon: Dict[str, Any] | None = None,
    cutadapt: Dict[str, Any] | None = None,
    vsearch: Dict[str, Any] | None = None,
    emu: Dict[str, Any] | None = None,
    humann: Dict[str, Any] | None = None,
    amrfinderplus: Dict[str, Any] | None = None,
    abricate: Dict[str, Any] | None = None,
    bakta: Dict[str, Any] | None = None,
    eggnog: Dict[str, Any] | None = None,
    maxbin2: Dict[str, Any] | None = None,
    concoct: Dict[str, Any] | None = None,
    das_tool: Dict[str, Any] | None = None,
    drep: Dict[str, Any] | None = None,
    assembly: Dict[str, Any] | None = None,
    metaspades: Dict[str, Any] | None = None,
    consensus: Dict[str, Any] | None = None,
    metaphlan: Dict[str, Any] | None = None,
    kaiju: Dict[str, Any] | None = None,
    multiqc: Dict[str, Any] | None = None,
    krona: Dict[str, Any] | None = None,
    mapdamage: Dict[str, Any] | None = None,
    instrain: Dict[str, Any] | None = None,
    antismash: Dict[str, Any] | None = None,
    dada2: Dict[str, Any] | None = None,
    differential: Dict[str, Any] | None = None,
    bracken_read_length_by_platform: Dict[str, int] | None = None,
) -> Dict[str, Any]:
    """Validate every section and return a clean config dict.

    If ``preset`` is given, its template config is the base and the user's section
    arguments are deep-merged on top (user wins). Raises ``ValidationError`` with a
    human-readable message on any problem, so the LLM can relay it and re-ask.
    """
    base: Dict[str, Any] = {}
    if preset:
        try:
            base = presets.get_preset_config(preset)
        except KeyError as e:
            raise registry.ValidationError(str(e))
        modules = presets.deep_merge(base.get("modules", {}), modules or {})
        sweep = sweep if sweep is not None else base.get("sweep")
        # Tool param sections are merged generically below (preset base under user overrides),
        # so adding a new tool needs no special-casing here.

    mods = {**DEFAULT_MODULES, **(modules or {})}

    if mods.get("classify") and not db.get("kraken2"):
        raise registry.ValidationError("classify is enabled but db.kraken2 is missing")
    if mods.get("abundance") and not db.get("bracken", db.get("kraken2")):
        raise registry.ValidationError("abundance is enabled but no Bracken db is set")
    if mods.get("binning") and not mods.get("assembly"):
        raise registry.ValidationError("binning requires assembly to be enabled too")
    if mods.get("bin_refinement") and not mods.get("binning"):
        raise registry.ValidationError(
            "bin_refinement (MaxBin2+CONCOCT→DAS_Tool→dRep) requires binning to be enabled"
        )
    if mods.get("reconcile") and not (mods.get("assembly") and mods.get("classify")):
        raise registry.ValidationError(
            "reconcile requires both assembly (for contigs) and classify (for read calls)"
        )
    if mods.get("filtered_assembly") and not (mods.get("assembly") and mods.get("classify")):
        raise registry.ValidationError(
            "filtered_assembly requires assembly (unfiltered baseline) and classify (read calls)"
        )
    if mods.get("stats") and not mods.get("abundance"):
        raise registry.ValidationError("stats (diversity) requires abundance (Bracken table)")
    if mods.get("differential"):
        if not mods.get("abundance"):
            raise registry.ValidationError(
                "differential abundance requires abundance (it tests the Bracken table)")
        gcol = str((differential or {}).get("group_column", "group"))
        counts = _sample_groups(samples, gcol)
        if isinstance(samples, list):  # inline sheets are checked now; TSV at runtime
            if len(counts) < 2:
                raise registry.ValidationError(
                    "differential needs >=2 sample groups; add a `group` column with "
                    "two labels (e.g. case/control)")
            small = [g for g, n in counts.items() if n < 2]
            if small:
                raise registry.ValidationError(
                    f"differential needs >=2 samples per group for the permutation test; "
                    f"under-replicated group(s): {small}")
    if mods.get("bgc") and not mods.get("assembly"):
        raise registry.ValidationError(
            "bgc (antiSMASH biosynthetic gene clusters) requires assembly — it mines contigs")
    if mods.get("classify_consensus") and not mods.get("classify"):
        raise registry.ValidationError(
            "classify_consensus needs classify (it cross-checks kraken2 against a 2nd classifier)"
        )
    if mods.get("aggregate") and not mods.get("classify"):
        raise registry.ValidationError(
            "aggregate (MultiQC + Krona) needs classify (it summarizes the kraken2 reports)"
        )
    if mods.get("damage"):
        if not mods.get("assembly"):
            raise registry.ValidationError(
                "damage (aDNA authentication) requires assembly — reads are mapped to contigs"
            )
        if isinstance(samples, list) and not any(
                str(r.get("library", "wgs")).lower() == "ancient" for r in samples):
            raise registry.ValidationError(
                "damage needs at least one sample with library=ancient to authenticate"
            )
    if mods.get("strain") and not mods.get("assembly"):
        raise registry.ValidationError(
            "strain (inStrain) requires assembly — it profiles SNVs over the contigs+mapping"
        )
    if mods.get("decontam"):
        if not mods.get("abundance"):
            raise registry.ValidationError(
                "decontam requires abundance (it operates on the combined Bracken table)"
            )
        if isinstance(samples, list) and not any(
                str(r.get("control", "")).lower() in {"1", "true", "yes", "y"}
                for r in samples):
            raise registry.ValidationError(
                "decontam needs at least one sample marked `control: true` (negative/blank)"
            )
    doms = [d.lower() for d in (domains or [])]
    bad = [d for d in doms if d not in KNOWN_DOMAINS]
    if bad:
        raise registry.ValidationError(
            f"unknown domain(s) {bad}; choose from {sorted(KNOWN_DOMAINS)}"
        )
    if mods.get("domain_taxonomy"):
        if not doms:
            raise registry.ValidationError("domain_taxonomy needs a non-empty `domains` list")
        if ("viral" in doms or "eukaryote" in doms) and not mods.get("assembly"):
            raise registry.ValidationError("viral/eukaryote domain taxonomy requires assembly")
        if "prokaryote" in doms and not mods.get("binning"):
            raise registry.ValidationError("prokaryote domain taxonomy requires binning (bins)")

    # Assembly is WGS-only: block assembly-based modules for an all-amplicon run.
    # `functional` (HUMAnN pathways + AMR + MAG annotation) is shotgun-only too.
    if _all_amplicon(samples):
        blocked = [m for m in ("assembly", "binning", "bin_refinement", "reconcile",
                               "filtered_assembly", "domain_taxonomy", "functional",
                               "classify_consensus", "damage", "strain", "bgc")
                   if mods.get(m)]
        if blocked:
            raise registry.ValidationError(
                f"all samples are amplicon — these WGS-only modules do not apply: "
                f"{blocked}. Disable them (amplicon uses the OTU/Emu branch), or mark WGS "
                "samples with library=wgs."
            )

    sections = {
        "fastp": fastp, "kraken2": kraken2, "bracken": bracken,
        "megahit": megahit, "metabat2": metabat2,
        "porechop_abi": porechop_abi, "chopper": chopper, "flye": flye,
        "genomad": genomad, "checkv": checkv, "checkm2": checkm2,
        "gtdbtk": gtdbtk, "eukrep": eukrep, "eukcc": eukcc,
        "cutadapt": cutadapt, "vsearch": vsearch, "emu": emu,
        "humann": humann, "amrfinderplus": amrfinderplus, "abricate": abricate,
        "bakta": bakta, "eggnog": eggnog,
        "maxbin2": maxbin2, "concoct": concoct, "das_tool": das_tool, "drep": drep,
        "metaspades": metaspades, "metaphlan": metaphlan, "kaiju": kaiju,
        "multiqc": multiqc, "krona": krona, "mapdamage": mapdamage, "instrain": instrain,
        "antismash": antismash, "dada2": dada2,
    }
    # Merge any preset-provided tool params under the user's overrides (user wins),
    # for every registry tool — so presets can tune new tools with no code changes here.
    if base:
        for tool in sections:
            merged = presets.deep_merge(base.get(tool, {}), sections.get(tool) or {})
            sections[tool] = merged or None
    cleaned_sections: Dict[str, Any] = {}
    for tool, values in sections.items():
        if values:
            cleaned_sections[tool] = registry.validate(tool, values)

    cfg: Dict[str, Any] = {
        "project": project,
        "outdir": outdir,
        "threads": int(threads),
        "samples": _validate_samples(samples),
        "db": {"kraken2": db.get("kraken2"), "bracken": db.get("bracken", db.get("kraken2"))},
        "modules": mods,
    }
    for extra in ("cat", "genomad", "checkv", "gtdbtk", "checkm2", "eukcc", "emu",
                  "humann_nucleotide", "humann_protein", "amrfinderplus", "bakta", "eggnog",
                  "metaphlan", "kaiju", "antismash"):
        if db.get(extra):
            cfg["db"][extra] = db[extra]
    if doms:
        cfg["domains"] = doms
    if amplicon:
        cfg["amplicon"] = {k: amplicon[k] for k in ("fwd_primer", "rev_primer")
                           if amplicon.get(k)}
        method = str(amplicon.get("method", "otu")).lower()
        if method not in {"otu", "asv"}:
            raise registry.ValidationError(
                "amplicon.method must be 'otu' (VSEARCH, short reads) or 'asv' (DADA2)")
        cfg["amplicon"]["method"] = method
    if host_removal:
        if not host_removal.get("genome"):
            raise registry.ValidationError("host_removal needs a 'genome' path")
        cfg["host_removal"] = {"genome": host_removal["genome"]}
    # Short-read assembler choice (megahit default | metaspades). Long reads always use Flye.
    asm_choice = str((assembly or {}).get("assembler", "megahit")).lower()
    if assembly:
        if asm_choice not in {"megahit", "metaspades"}:
            raise registry.ValidationError(
                "assembly.assembler must be 'megahit' or 'metaspades'")
        cfg["assembly"] = {"assembler": asm_choice}
    if (isinstance(samples, list) and any(r.get("long_reads") for r in samples)
            and asm_choice != "metaspades"):
        raise registry.ValidationError(
            "a sample carries hybrid long_reads, which requires assembly.assembler=metaspades")
    # Second-classifier choice for the consensus cross-check (metaphlan default | kaiju).
    if consensus:
        clf = str(consensus.get("classifier", "metaphlan")).lower()
        if clf not in {"metaphlan", "kaiju"}:
            raise registry.ValidationError(
                "consensus.classifier must be 'metaphlan' or 'kaiju'")
        cfg["consensus"] = {"classifier": clf}
    # Per-platform Bracken read length (a routing hint, not a Bracken flag — kept top-level so
    # registry validation of the bracken section doesn't reject it). Per-sample sheet field wins.
    if bracken_read_length_by_platform:
        bad = [k for k in bracken_read_length_by_platform if k not in KNOWN_PLATFORMS]
        if bad:
            raise registry.ValidationError(
                f"bracken_read_length_by_platform: unknown platform(s) {bad}")
        cfg["bracken_read_length_by_platform"] = {
            k: int(v) for k, v in bracken_read_length_by_platform.items()}
    if preset:
        cfg["preset"] = preset
    sweep_clean = _validate_sweep(sweep, cleaned_sections.get("kraken2", {}))
    if sweep_clean:
        cfg["sweep"] = sweep_clean
    sub_clean = _validate_subsample(subsample)
    if sub_clean:
        cfg["subsample"] = sub_clean
    rf_clean = _validate_read_filter(read_filter)
    if rf_clean:
        cfg["read_filter"] = rf_clean
    if mods.get("differential") or differential:
        cfg["differential"] = _validate_differential(differential)
    cfg.update(cleaned_sections)
    return cfg


def write_config(cfg: Dict[str, Any], path: str = "config.yaml") -> str:
    """Write the config dict to YAML and return the absolute path."""
    with open(path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, default_flow_style=False)
    return os.path.abspath(path)
