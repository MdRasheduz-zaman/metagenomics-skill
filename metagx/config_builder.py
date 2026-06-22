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

# (the canonical module->tools map lives in tool_advisor.MODULE_TOOLS; a stale duplicate here
# and in report.py was dead code and was removed to avoid drift — see metagx/wiring.py)

# Optional `db.<key>` paths accepted beyond the always-present kraken2/bracken. Each is a
# per-tool module/validation DB; keep in sync with dbprovision.SPECS (the wiring audit checks
# this) — "cat" is the only one with no provisioner (built via the custom CAT helper).
DB_EXTRA_KEYS = ("cat", "genomad", "checkv", "gtdbtk", "checkm2", "eukcc", "emu",
                 "humann_nucleotide", "humann_protein", "amrfinderplus", "bakta", "eggnog",
                 "metaphlan", "kaiju", "antismash", "blast")

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
    "phylogenetics": False,
    "validate": False,
}

KNOWN_DOMAINS = {"viral", "prokaryote", "eukaryote"}


KNOWN_PLATFORMS = {"illumina", "mgi", "bgi", "ont", "nanopore",
                   "pacbio_hifi", "pacbio_clr", "pacbio"}
LONG_PLATFORMS = {"ont", "nanopore", "pacbio_hifi", "pacbio_clr", "pacbio"}


def _load_probe(probe: Any) -> Dict[str, Any] | None:
    """Accept a probe report dict or a path to probe.json (JSON is valid YAML)."""
    if probe is None:
        return None
    if isinstance(probe, str):
        with open(probe) as fh:
            return yaml.safe_load(fh)
    return probe


def _backfill_platforms(samples: Any, probe_report: Dict[str, Any] | None) -> List[tuple]:
    """Fill a MISSING per-sample platform from the probe's inferred class. Never overrides a
    declared platform. Returns the (sample, inferred) pairs applied, for provenance."""
    applied: List[tuple] = []
    if not isinstance(samples, list) or not probe_report:
        return applied
    profiles = probe_report.get("samples", {})
    for rec in samples:
        name = rec.get("sample") or rec.get("name")
        inferred = profiles.get(name, {}).get("inferred_platform_class")
        if inferred and not str(rec.get("platform", "")).strip():
            rec["platform"] = inferred
            applied.append((name, inferred))
    return applied


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
        mc_samples = int(diff.get("mc_samples", 128))
    except (TypeError, ValueError):
        raise registry.ValidationError(
            "differential: n_permutations/min_count/mc_samples must be ints and fdr a number")
    if n_perm < 99:
        raise registry.ValidationError("differential.n_permutations must be >= 99")
    if not (0 < fdr < 1):
        raise registry.ValidationError("differential.fdr must be in (0, 1)")
    if mc_samples < 1:
        raise registry.ValidationError(
            "differential.mc_samples must be >= 1 (1 = single point estimate; ALDEx2 uses 128)")
    out = {"group_column": str(diff.get("group_column", "group")),
           "n_permutations": n_perm, "fdr": fdr, "min_count": min_count,
           "mc_samples": mc_samples, "seed": int(diff.get("seed", 42))}
    if diff.get("reference_group"):
        out["reference_group"] = str(diff["reference_group"])
    return out


def _validate_validate(val: Dict[str, Any] | None, db: Dict[str, Any]) -> Dict[str, Any]:
    """Settings for the BLAST validation module (cross-check classifier calls vs alignment).

    Needs either a local BLAST DB (db.blast) or remote NCBI search (validate.remote). The
    nt DB is ~200 GB, so remote is the dev-box / few-sequences escape hatch.
    """
    val = val or {}
    target = str(val.get("target", "reads")).lower()
    if target not in {"reads"}:
        raise registry.ValidationError(
            "validate.target must be 'reads' (contig-level BLAST validation is not yet wired)")
    level = str(val.get("level", "genus")).lower()
    if level not in {"genus", "species"}:
        raise registry.ValidationError("validate.level must be genus or species")
    try:
        top_n = int(val.get("top_n", 10))
        reads_per_taxon = int(val.get("reads_per_taxon", 50))
    except (TypeError, ValueError):
        raise registry.ValidationError("validate.top_n and reads_per_taxon must be ints")
    if top_n < 1 or reads_per_taxon < 1:
        raise registry.ValidationError("validate.top_n and reads_per_taxon must be >= 1")
    remote = bool(val.get("remote", False))
    # build_from: keep the validation reference IN SCOPE with the classifier by building the
    # BLAST DB from the SAME genomes. Either a FASTA/folder path, or "classifier" => reuse the
    # db.build source (only meaningful for custom-fasta/custom-folder builds, where a local
    # source FASTA exists). Validating against a broader DB is a different benchmark.
    build_from = val.get("build_from")
    if build_from:
        build_from = str(build_from)
        if build_from == "classifier":
            # Resolve at run time from the kraken2 DB dir (its custom_library.fasta / library/
            # genomes + seqid2taxid.map, so the BLAST DB carries kraken2's exact taxids). Needs
            # the genomes ON DISK — a prebuilt/fetched or --clean'd index has only *.k2d and the
            # build rule will fail fast (validation.kraken2_db_sources) telling the user to pass a
            # FASTA. So require a kraken2 DB path here.
            if not db.get("kraken2"):
                raise registry.ValidationError(
                    "validate.build_from: 'classifier' needs db.kraken2 (the classifier DB dir) "
                    "so the BLAST DB can be built from its genomes + seqid2taxid.map. A "
                    "prebuilt/fetched index ships no genomes — then point validate.build_from at "
                    "the genome FASTA(s) you used, or set db.blast / validate.remote.")
    elif not remote and not db.get("blast"):
        raise registry.ValidationError(
            "validate needs a BLAST database: set db.blast to a local BLAST+ db, or "
            "validate.build_from to a FASTA/folder of the SAME genomes as your classifier DB "
            "(keeps the benchmark in scope), or validate.remote: true to search NCBI remotely "
            "(only practical for a few sequences). nt is ~200 GB — never auto-fetched.")
    out: Dict[str, Any] = {"target": target, "level": level, "top_n": top_n,
                           "reads_per_taxon": reads_per_taxon, "remote": remote,
                           "rank": "S" if level == "species" else "G",
                           "seed": int(val.get("seed", 42))}
    if build_from:
        out["build_from"] = build_from
    return out


def _validate_phylogenetics(phylo: Dict[str, Any]) -> Dict[str, Any]:
    if not phylo.get("input") and not phylo.get("aligned_input"):
        raise registry.ValidationError(
            "phylogenetics needs 'input' (FASTA) or 'aligned_input' (pre-aligned FASTA)"
        )
    method = str(phylo.get("method", "iqtree")).lower()
    if method not in {"iqtree", "fasttree", "auto"}:
        raise registry.ValidationError(
            "phylogenetics.method must be iqtree, fasttree, or auto"
        )
    seq_type = str(phylo.get("sequence_type", "nt")).lower()
    if seq_type not in {"nt", "aa"}:
        raise registry.ValidationError("phylogenetics.sequence_type must be nt or aa")
    out: Dict[str, Any] = {
        "method": method,
        "sequence_type": seq_type,
        "trim": bool(phylo.get("trim", True)),
    }
    if phylo.get("input"):
        out["input"] = str(phylo["input"])
    if phylo.get("aligned_input"):
        out["aligned_input"] = str(phylo["aligned_input"])
    if phylo.get("trimal_method"):
        out["trimal_method"] = str(phylo["trimal_method"])
    if phylo.get("fasttree_threshold") is not None:
        out["fasttree_threshold"] = int(phylo["fasttree_threshold"])
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


def _any_provided_contigs(samples: Any) -> bool:
    """True if any sample supplies pre-assembled contigs (list records or a TSV path)."""
    if isinstance(samples, list):
        return any(str(r.get("contigs", "")).strip() for r in samples)
    if isinstance(samples, str) and os.path.isfile(samples):
        import csv
        with open(samples) as fh:
            return any(str(row.get("contigs", "")).strip()
                       for row in csv.DictReader(fh, delimiter="\t"))
    return False


# NCBI libraries that `kraken2-build --download-library` accepts.
_KRAKEN2_LIBRARIES = {"bacteria", "viral", "archaea", "fungi", "protozoa", "human",
                      "plasmid", "UniVec_Core", "nt"}
# Map sample platform -> the Bracken read length the DB should carry (mirrors the
# bracken registry `recommend:` table). Used to derive read_lengths: auto.
_PLATFORM_READ_LENGTH = {"illumina": 150, "mgi": 150, "bgi": 150,
                         "ont": 1000, "pacbio_hifi": 1000, "pacbio_clr": 1000}


def _sheet_platforms(path: str) -> set:
    """The set of platform labels declared in a sample-sheet TSV."""
    if not (isinstance(path, str) and os.path.isfile(path)):
        return set()
    import csv
    with open(path) as fh:
        return {str(row.get("platform", "")).lower() for row in csv.DictReader(fh, delimiter="\t")}


def _platform_read_lengths(samples: Any) -> List[int]:
    """The sorted set of Bracken read lengths implied by the samples' platforms."""
    if isinstance(samples, list):
        plats = {str(r.get("platform", "")).lower() for r in samples}
    else:
        plats = _sheet_platforms(samples)
    lengths = {_PLATFORM_READ_LENGTH.get(p, 150) for p in plats if p}
    return sorted(lengths) or [150]


def _validate_db_build(build: Dict[str, Any] | None, samples: Any) -> Dict[str, Any] | None:
    """Validate a ``db.build`` block: how to construct the kraken2 + Bracken DB.

    The registry (kraken2-build) validates the strategy/taxonomy/libraries/tuning params;
    this enforces the cross-field invariants the per-param check can't (a source for
    custom/spike-in, libraries for standard/spike-in, spike-in => real taxonomy,
    minimizer<=kmer) and locks Bracken's ``-k`` to the kraken2 ``--kmer-len``. Read lengths
    default to ``auto`` => derived from the sample sheet's platforms, so every length the run
    classifies at has a matching ``databaseLmers.kmer_distrib``.
    """
    if not build:
        return None
    build = dict(build)
    # orchestration keys consumed by the dbbuild script, not kraken2-build flags
    source = build.pop("source", None)
    read_lengths = build.pop("read_lengths", "auto")
    auto = bool(build.pop("auto", True))
    # build the aligned BLAST validation DB together with kraken2? None = undecided (build_config
    # may default it on when modules.validate is enabled). Not a kraken2-build flag.
    blast = build.pop("blast", None)
    # NCBI deprecated rsync, so kraken2-build's default rsync downloads now fail; --use-ftp
    # (wget) is the working path, hence the default. Override to False only on a host where
    # rsync to NCBI still works and you want its speed.
    use_ftp = bool(build.pop("use_ftp", True))
    download_on = str(build.pop("download_on", "rule"))
    if download_on not in {"rule", "login"}:
        raise registry.ValidationError("db.build.download_on must be 'rule' or 'login'")

    params = registry.validate("kraken2-build", build)  # strategy/taxonomy/libraries + tuning
    strategy = params.get("strategy", "standard")
    taxonomy = params.get("taxonomy", "real")

    if strategy in {"custom-fasta", "custom-folder", "spike-in"} and not source:
        raise registry.ValidationError(
            f"db.build.strategy={strategy} needs db.build.source (a FASTA, or a folder of FASTAs)")
    if strategy in {"standard", "spike-in"} and not params.get("libraries"):
        raise registry.ValidationError(
            f"db.build.strategy={strategy} needs db.build.libraries (e.g. 'viral' or 'bacteria,viral')")
    if strategy == "spike-in" and taxonomy != "real":
        raise registry.ValidationError(
            "db.build.strategy=spike-in requires taxonomy: real — synthetic taxids cannot merge "
            "into a standard library's taxonomy tree")
    if params.get("libraries"):
        bad = [l for l in (x.strip() for x in str(params["libraries"]).split(","))
               if l and l not in _KRAKEN2_LIBRARIES]
        if bad:
            raise registry.ValidationError(
                f"db.build.libraries has unknown NCBI libraries {bad}; choose from "
                f"{sorted(_KRAKEN2_LIBRARIES)}")
    kmer = int(params.get("kmer_len", 35))
    if int(params.get("minimizer_len", 31)) > kmer:
        raise registry.ValidationError("db.build.minimizer_len must be <= kmer_len")

    if read_lengths == "auto":
        read_lengths = _platform_read_lengths(samples)
    elif isinstance(read_lengths, (list, tuple)) and all(isinstance(x, int) for x in read_lengths):
        read_lengths = sorted(set(read_lengths))
    else:
        raise registry.ValidationError(
            "db.build.read_lengths must be 'auto' or a list of integer read lengths")

    out = dict(params)
    # Materialize the resolved orchestration choices so the stored config is self-describing
    # (registry.validate only returns keys the user passed; defaults live in the registry).
    out.update(strategy=strategy, taxonomy=taxonomy,
               source=source, read_lengths=read_lengths, auto=auto, download_on=download_on,
               use_ftp=use_ftp,
               bracken_kmer_len=kmer)  # Bracken -k is locked to the kraken2 --kmer-len
    if blast is not None:
        out["blast"] = bool(blast)
    return out


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
    probe: Any = None,
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
    phylogenetics: Dict[str, Any] | None = None,
    mafft: Dict[str, Any] | None = None,
    iqtree: Dict[str, Any] | None = None,
    fasttree: Dict[str, Any] | None = None,
    validate: Dict[str, Any] | None = None,
    blastn: Dict[str, Any] | None = None,
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

    # Probe-conditioned routing: backfill any MISSING per-sample platform from measured
    # inference (a gap-fill — a declared platform is never overridden; mismatches surface as
    # warnings in cfg["probe"]). Parameter values are not auto-set here; the probe context
    # promotes interview questions instead (no silent tuning).
    probe_report = _load_probe(probe)
    if isinstance(samples, list):
        samples = [dict(r) for r in samples]  # don't mutate the caller's records
    backfilled = _backfill_platforms(samples, probe_report)

    # Contig-consuming modules (binning/reconcile/bgc/strain/domain/functional-AMR) are
    # satisfied either by running the assembler OR by a sample supplying pre-assembled
    # contigs (an isolate genome, a prior assembly, or references) — see the `contigs` column.
    have_contigs = _any_provided_contigs(samples)
    assembly_ok = bool(mods.get("assembly") or have_contigs)

    db = dict(db or {})
    db_build = _validate_db_build(db.get("build"), samples)
    # A configured db.build writes the DB to db.kraken2; default that path so the user need
    # not repeat it. The DB need not exist yet — the build step produces it.
    if db_build and not db.get("kraken2"):
        db["kraken2"] = os.path.join(outdir, "dbs", str(db_build.get("strategy", "standard")))

    # Build the ALIGNED BLAST validation DB together with kraken2 when validating — while the
    # library genomes are guaranteed present (before any --clean). Default it on for validate;
    # the user can set db.build.blast: false explicitly (doctor then warns strongly). A
    # co-located in-scope db.blast means validate needs no build_from / external nt DB.
    if db_build:
        if mods.get("validate") and db_build.get("blast") is None:
            db_build["blast"] = True
        if db_build.get("blast") and not db.get("blast"):
            db["blast"] = os.path.join(db["kraken2"], "blast", "insync")

    if mods.get("classify") and not db.get("kraken2"):
        raise registry.ValidationError(
            "classify is enabled but db.kraken2 is missing (set a path, or add a db.build block)")
    if mods.get("abundance") and not db.get("bracken", db.get("kraken2")):
        raise registry.ValidationError("abundance is enabled but no Bracken db is set")
    if mods.get("binning") and not assembly_ok:
        raise registry.ValidationError(
            "binning requires assembly to be enabled (or a sample with pre-assembled contigs)")
    if mods.get("bin_refinement") and not mods.get("binning"):
        raise registry.ValidationError(
            "bin_refinement (MaxBin2+CONCOCT→DAS_Tool→dRep) requires binning to be enabled"
        )
    if mods.get("reconcile") and not (assembly_ok and mods.get("classify")):
        raise registry.ValidationError(
            "reconcile requires contigs (assembly or provided) and classify (for read calls)"
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
    if mods.get("bgc") and not assembly_ok:
        raise registry.ValidationError(
            "bgc (antiSMASH biosynthetic gene clusters) requires contigs "
            "(enable assembly or provide pre-assembled contigs) — it mines contigs")
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
    if mods.get("strain") and not assembly_ok:
        raise registry.ValidationError(
            "strain (inStrain) requires contigs (assembly or provided) — it profiles SNVs "
            "over the contigs+mapping"
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
    if mods.get("phylogenetics"):
        _validate_phylogenetics(phylogenetics or {})
    doms = [d.lower() for d in (domains or [])]
    bad = [d for d in doms if d not in KNOWN_DOMAINS]
    if bad:
        raise registry.ValidationError(
            f"unknown domain(s) {bad}; choose from {sorted(KNOWN_DOMAINS)}"
        )
    if mods.get("domain_taxonomy"):
        if not doms:
            raise registry.ValidationError("domain_taxonomy needs a non-empty `domains` list")
        if ("viral" in doms or "eukaryote" in doms) and not assembly_ok:
            raise registry.ValidationError(
                "viral/eukaryote domain taxonomy requires contigs (assembly or provided)")
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
        "mafft": mafft, "iqtree": iqtree, "fasttree": fasttree,
        "blastn": blastn,
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
    for extra in DB_EXTRA_KEYS:
        if db.get(extra):
            cfg["db"][extra] = db[extra]
    if db_build:
        cfg["db"]["build"] = db_build
    # db.provision: module DBs to auto-fetch at run time (idempotent; skip if present).
    provision = db.get("provision")
    if provision:
        from . import dbprovision
        if not isinstance(provision, list):
            raise registry.ValidationError("db.provision must be a list of tool names")
        bad = [t for t in provision if t not in dbprovision.SPECS]
        if bad:
            raise registry.ValidationError(
                f"db.provision: no provisioner for {bad}; known: {sorted(dbprovision.SPECS)}")
        manual = [t for t in provision if dbprovision.SPECS[t].get("manual")]
        if manual:
            raise registry.ValidationError(
                f"db.provision can't auto-fetch {manual} (no CLI downloader); download them "
                f"manually and set db.<tool> to the path instead")
        missing = [t for t in provision if not cfg["db"].get(t)]
        if missing:
            raise registry.ValidationError(
                f"db.provision lists {missing} but their db.<tool> output path(s) are unset")
        cfg["db"]["provision"] = list(provision)
    if probe_report:  # provenance: what the probe measured + any backfills/warnings it drove
        proj = probe_report.get("project", {})
        cfg["probe"] = {
            "measured": probe_report.get("measured", False),
            "platform_consensus": proj.get("platform_consensus"),
            "warnings": proj.get("warnings", []),
            "backfilled_platforms": dict(backfilled),
        }
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
    brpp: Dict[str, int] = {}
    # A db.build builds databaseLmers for exactly the lengths each platform implies, so the
    # classify-time Bracken -r MUST match — backfill the per-platform map from the same source.
    if db_build:
        for rec in (cfg["samples"] if isinstance(cfg["samples"], list) else []):
            p = str(rec.get("platform", "")).lower()
            if p in _PLATFORM_READ_LENGTH:
                brpp[p] = _PLATFORM_READ_LENGTH[p]
        if isinstance(cfg["samples"], str):  # TSV: derive from the sheet's platforms
            for p in _sheet_platforms(cfg["samples"]):
                if p in _PLATFORM_READ_LENGTH:
                    brpp[p] = _PLATFORM_READ_LENGTH[p]
    if bracken_read_length_by_platform:
        bad = [k for k in bracken_read_length_by_platform if k not in KNOWN_PLATFORMS]
        if bad:
            raise registry.ValidationError(
                f"bracken_read_length_by_platform: unknown platform(s) {bad}")
        brpp.update({k: int(v) for k, v in bracken_read_length_by_platform.items()})  # user wins
    if brpp:
        cfg["bracken_read_length_by_platform"] = brpp
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
    if mods.get("phylogenetics"):
        cfg["phylogenetics"] = _validate_phylogenetics(phylogenetics or {})
    if mods.get("validate"):
        if not mods.get("classify"):
            raise registry.ValidationError(
                "validate cross-checks classifier calls, so it needs classify enabled")
        cfg["validate"] = _validate_validate(validate, db)
    cfg.update(cleaned_sections)
    return cfg


def write_config(cfg: Dict[str, Any], path: str = "config.yaml") -> str:
    """Write the config dict to YAML and return the absolute path."""
    with open(path, "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False, default_flow_style=False)
    return os.path.abspath(path)
