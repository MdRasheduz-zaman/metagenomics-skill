"""Metagenomics skill — MCP server (primary surface) + HTTP wrapper for web agents.

The tools are thin: all real logic lives in the ``metagx`` package, whose parameter
registries are the single source of truth. The intended LLM flow is:

  1. list_tools()                  -> which steps exist
  2. get_interview(tool)           -> questions to ask the user, with types/ranges
  3. ...conduct the interview in natural language...
  4. build_config(...)             -> validate answers, write config.yaml
  5. run_pipeline(config)          -> execute the Snakemake workflow
  6. get_results(config)           -> read back the comparison matrices

Run as MCP (stdio):     python mcp_server.py
Run as HTTP for web:    uvicorn mcp_server:app --host 0.0.0.0 --port 8000

NOTE: do NOT add ``from __future__ import annotations`` here. FastMCP introspects each
tool's signature and skips generic annotations via ``get_origin(param.annotation)``; with
stringized (PEP 563) annotations that guard misses and ``issubclass(str, Context)`` raises
"arg 1 must be a class". Real annotation objects keep tool registration working.
"""

import glob
import inspect
import json
import os
from typing import Any, Dict, List, Optional

from mcp.server.fastmcp import FastMCP

from metagx import (
    advise,
    catalog,
    config_builder,
    dbbuild,
    evidence_pack,
    history,
    paper,
    presets,
    probe,
    registry,
    report,
    runner,
    schedulers,
    sync_help,
    tool_advisor,
)

mcp = FastMCP("metagx", json_response=True)


# --------------------------------------------------------------------------- #
# Discovery / interview                                                        #
# --------------------------------------------------------------------------- #
@mcp.tool()
def list_pipeline_tools() -> str:
    """List every pipeline tool that has a parameter registry, with a description.

    Use this first to decide which modules (qc/classify/abundance/assembly/binning)
    the user wants.
    """
    return json.dumps(
        {t: registry.load_registry(t)["description"].strip() for t in registry.list_tools()},
        indent=2,
    )


@mcp.tool()
def get_parameters(tool: str) -> str:
    """Return the COMPLETE parameter registry for a tool (every flag, type, range).

    Use when the user wants fine control or asks 'what can I configure for X'.
    """
    try:
        return json.dumps(registry.load_registry(tool), indent=2)
    except registry.ValidationError as e:
        return f"error: {e}"


@mcp.tool()
def list_presets() -> str:
    """List ready-to-run workflow presets (name, description, when to use).

    Offer these at the START of the interview: a preset pre-fills sensible modules and
    parameters, then the user only adjusts what they care about. Pass the chosen name as
    `preset` to build_config.
    """
    return json.dumps(presets.describe_presets(), indent=2)


@mcp.tool()
def get_preset(name: str) -> str:
    """Return the template config a preset contributes (so you can show/explain it)."""
    try:
        return json.dumps(presets.get_preset_config(name), indent=2)
    except KeyError as e:
        return f"error: {e}"


@mcp.tool()
def list_schedulers() -> str:
    """List HPC scheduler backends for run_pipeline's `executor` (local/slurm/lsf/sge/pbs/
    generic): the executor each uses, the plugin it needs, and what to edit before first use.
    """
    return json.dumps(schedulers.describe(), indent=2)


@mcp.tool()
def get_interview(tool: str, max_tier: int = 2, context: dict | None = None) -> str:
    """Return the questions to ask the user for a tool, with answer constraints.

    tier 1 = essential, 2 = common (default), 3 = advanced. Each item includes the
    type, default, min/max/choices, whether it is sweepable, and a ready-to-use
    natural-language question. Drive your interview from this — do not invent flags.

    Pass ``context`` with what you already know about the run (e.g.
    {"goal": "strain_resolved", "estimated_bases": 6e10}) to surface normally-quiet
    params that the goal/data make relevant; promoted items carry a ``promoted`` note.
    """
    return json.dumps(registry.interview_spec(tool, max_tier=max_tier, context=context), indent=2)


@mcp.tool()
def run_probe(samples: str, consent: bool = False, max_reads: int = 100_000,
              max_samples: int | None = None, out: str = "", host_index: str = "") -> str:
    """Measure read stats from the user's samples to drive data-conditioned promotion.

    LOCAL and consent-gated: profiles a bounded head subsample of EVERY sample (read length,
    quality, GC, duplication, inferred platform), entirely on this machine, emitting only
    non-reconstructive aggregates — never read sequences or IDs. Returns per-sample profiles,
    a reconciled project view with sample-sheet mismatch warnings, and a ``context`` dict to
    pass straight into get_interview(context=).

    You MUST obtain the user's explicit permission before setting consent=True (it reads their
    read files locally). Without consent this returns {"measured": false} and you should fall
    back to a-priori suggestions. The probe never sends anything off the machine.
    """
    res = probe.run(samples, max_reads=max_reads, max_samples=max_samples,
                    out=(out or None), assume_yes=consent, host_index=(host_index or None))
    return json.dumps(res, indent=2)


# --------------------------------------------------------------------------- #
# Config building                                                             #
# --------------------------------------------------------------------------- #
@mcp.tool()
def build_config(
    samples: Any,
    db: Dict[str, str],
    project: str = "run",
    outdir: str = "results",
    threads: int = 8,
    preset: Optional[str] = None,
    modules: Optional[Dict[str, bool]] = None,
    sweep: Optional[Dict[str, Any]] = None,
    subsample: Optional[Dict[str, Any]] = None,
    read_filter: Optional[Dict[str, Any]] = None,
    fastp: Optional[Dict[str, Any]] = None,
    kraken2: Optional[Dict[str, Any]] = None,
    bracken: Optional[Dict[str, Any]] = None,
    megahit: Optional[Dict[str, Any]] = None,
    metabat2: Optional[Dict[str, Any]] = None,
    porechop_abi: Optional[Dict[str, Any]] = None,
    chopper: Optional[Dict[str, Any]] = None,
    flye: Optional[Dict[str, Any]] = None,
    domains: Optional[list] = None,
    genomad: Optional[Dict[str, Any]] = None,
    checkv: Optional[Dict[str, Any]] = None,
    checkm2: Optional[Dict[str, Any]] = None,
    gtdbtk: Optional[Dict[str, Any]] = None,
    eukrep: Optional[Dict[str, Any]] = None,
    eukcc: Optional[Dict[str, Any]] = None,
    amplicon: Optional[Dict[str, Any]] = None,
    cutadapt: Optional[Dict[str, Any]] = None,
    vsearch: Optional[Dict[str, Any]] = None,
    emu: Optional[Dict[str, Any]] = None,
    host_removal: Optional[Dict[str, Any]] = None,
    # functional layer (modules.functional)
    humann: Optional[Dict[str, Any]] = None,
    amrfinderplus: Optional[Dict[str, Any]] = None,
    abricate: Optional[Dict[str, Any]] = None,
    bakta: Optional[Dict[str, Any]] = None,
    eggnog: Optional[Dict[str, Any]] = None,
    # bin refinement (modules.bin_refinement)
    maxbin2: Optional[Dict[str, Any]] = None,
    concoct: Optional[Dict[str, Any]] = None,
    das_tool: Optional[Dict[str, Any]] = None,
    drep: Optional[Dict[str, Any]] = None,
    # assembler choice + metaSPAdes/hybrid
    assembly: Optional[Dict[str, Any]] = None,
    metaspades: Optional[Dict[str, Any]] = None,
    # second-classifier consensus (modules.classify_consensus)
    consensus: Optional[Dict[str, Any]] = None,
    metaphlan: Optional[Dict[str, Any]] = None,
    kaiju: Optional[Dict[str, Any]] = None,
    # aggregate report (modules.aggregate)
    multiqc: Optional[Dict[str, Any]] = None,
    krona: Optional[Dict[str, Any]] = None,
    # Tier 3: aDNA damage, strain, per-sample Bracken length
    mapdamage: Optional[Dict[str, Any]] = None,
    instrain: Optional[Dict[str, Any]] = None,
    bracken_read_length_by_platform: Optional[Dict[str, int]] = None,
    # gap-closing modules: BGC mining, ASV amplicon, differential abundance
    antismash: Optional[Dict[str, Any]] = None,
    dada2: Optional[Dict[str, Any]] = None,
    differential: Optional[Dict[str, Any]] = None,
    probe: Any = None,
    config_path: str = "config.yaml",
) -> str:
    """Validate interview answers and write config.yaml for the workflow.

    samples: a path to a TSV (columns: sample, r1, r2, platform, layout) or a list of
             {"sample","r1","r2"?,"platform"?,"layout"?} records. platform =
             illumina|mgi (short) | ont | pacbio_hifi | pacbio_clr (long). layout =
             se|pe|interleaved (default pe if r2 else se). Long-read platforms are
             single-end and select porechop_abi+chopper / chopper QC and Flye assembly.
    db:      {"kraken2": "<path>", "bracken": "<path>"} (bracken defaults to kraken2).
    preset:  optional preset name (see list_presets); user sections override it.
    sweep:   optional {"param":"confidence","values":[0.0,0.1,0.5]} to run the matrix.
    subsample: optional {"fraction":0.2,"seed":42} to classify a random read subset
             (faster/cheaper; single-end only). Works on FASTA or FASTQ.
    The tool-named sections (kraken2, fastp, ...) take ONLY user params from
    get_interview/get_parameters; managed flags (db/threads/io) are added automatically.

    Returns the written config as JSON, or a validation error to relay and re-ask.
    """
    try:
        cfg = config_builder.build_config(
            project=project, outdir=outdir, threads=threads, samples=samples, db=db,
            preset=preset, modules=modules, sweep=sweep, subsample=subsample,
            read_filter=read_filter, probe=probe,
            fastp=fastp, kraken2=kraken2, bracken=bracken, megahit=megahit, metabat2=metabat2,
            porechop_abi=porechop_abi, chopper=chopper, flye=flye,
            domains=domains, genomad=genomad, checkv=checkv, checkm2=checkm2,
            gtdbtk=gtdbtk, eukrep=eukrep, eukcc=eukcc,
            amplicon=amplicon, cutadapt=cutadapt, vsearch=vsearch, emu=emu,
            host_removal=host_removal,
            humann=humann, amrfinderplus=amrfinderplus, abricate=abricate, bakta=bakta,
            eggnog=eggnog, maxbin2=maxbin2, concoct=concoct, das_tool=das_tool, drep=drep,
            assembly=assembly, metaspades=metaspades, consensus=consensus,
            metaphlan=metaphlan, kaiju=kaiju, multiqc=multiqc, krona=krona,
            mapdamage=mapdamage, instrain=instrain,
            antismash=antismash, dada2=dada2, differential=differential,
            bracken_read_length_by_platform=bracken_read_length_by_platform,
        )
    except registry.ValidationError as e:
        return json.dumps({"ok": False, "error": str(e)}, indent=2)
    path = config_builder.write_config(cfg, config_path)
    return json.dumps({"ok": True, "path": path, "config": cfg}, indent=2)


# --------------------------------------------------------------------------- #
# Execution / results                                                         #
# --------------------------------------------------------------------------- #
@mcp.tool()
def build_database(genomes: str, db_dir: str, read_length: int = 150, threads: int = 4,
                   dry_run: bool = False) -> str:
    """Build a custom kraken2 + Bracken database from a FASTA of reference genomes.

    Assigns each genome a synthetic taxid and a minimal taxonomy (no NCBI download), then
    runs kraken2-build and bracken-build. Use when the user has their own reference genomes
    rather than a prebuilt index. dry_run writes the taxonomy/library and returns the
    commands without building. Returns the taxid map, commands, and per-step logs.
    """
    result = dbbuild.build_db(genomes=genomes, db_dir=db_dir, read_length=read_length,
                              threads=threads, run=not dry_run)
    return json.dumps(result, indent=2)


@mcp.tool()
def build_kaiju_database(genomes: str, db_dir: str, taxonomy_dir: str, threads: int = 4,
                         dry_run: bool = False) -> str:
    """Build a custom Kaiju (protein) database from reference genomes — no NCBI download.

    Predicts proteins (prodigal), labels each with its genome's synthetic taxid, and builds
    the Kaiju FM-index. The output dir is a drop-in `db.kaiju` for the consensus module
    (`modules.classify_consensus` with `consensus.classifier: kaiju`). `taxonomy_dir` is a
    dir with names.dmp + nodes.dmp (e.g. the kraken2 db's taxonomy/ from build_database, so
    the taxids line up for the kraken2-vs-kaiju cross-check). Returns commands + per-step status.
    """
    result = dbbuild.build_kaiju_db(genomes=genomes, db_dir=db_dir, taxonomy_dir=taxonomy_dir,
                                    threads=threads, run=not dry_run)
    return json.dumps(result, indent=2)


@mcp.tool()
def run_pipeline(config_path: str = "config.yaml", cores: str = "all", dry_run: bool = False,
                 use_conda: bool = False, executor: str | None = None) -> str:
    """Run the Snakemake workflow against a config. Set dry_run to preview the plan.

    use_conda lets Snakemake auto-provision per-rule tools (workflow/envs/) — needed for the
    domain-taxonomy tools (geNomad/CheckV, GTDB-Tk/CheckM2, EukRep/EukCC) if not installed.

    executor submits to an HPC scheduler via a bundled profile (local/slurm/lsf/sge/pbs/
    generic — see list_schedulers). The bundled profile must be edited for your site first.
    """
    profile = None
    if executor:
        try:
            profile = schedulers.profile_path(executor)
        except (KeyError, FileNotFoundError) as e:
            return json.dumps({"ok": False, "error": str(e)}, indent=2)
    proc = runner.run(config=config_path, cores=cores, dry_run=dry_run,
                      use_conda=use_conda, profile=profile)
    ok = proc.returncode == 0
    tail = (proc.stdout or "")[-3000:] + "\n" + (proc.stderr or "")[-3000:]
    return json.dumps({"ok": ok, "returncode": proc.returncode, "log": tail}, indent=2)


@mcp.tool()
def get_results(config_path: str = "config.yaml") -> str:
    """Read back the JSON summary matrices produced by a completed run."""
    try:
        import yaml
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
    except FileNotFoundError:
        return f"error: {config_path} not found"
    out = os.path.join(cfg.get("outdir", "results"), cfg.get("project", "run"))
    summaries = {}
    for path in sorted(glob.glob(os.path.join(out, "summary", "*.json"))):
        with open(path) as fh:
            summaries[os.path.basename(path)] = json.load(fh)
    if not summaries:
        return f"No summaries under {out}/summary/. Run the pipeline first."
    return json.dumps(summaries, indent=2)


@mcp.tool()
def compare_platforms(manifest: str = "", outdir: str = "", paper: bool = False) -> str:
    """Compare the SAME sample sequenced on multiple platforms (ONT/Illumina/PacBio/...).

    Reads a TSV manifest (columns: label, platform, platform_class, kreport, reads, contigs,
    reference; blank cells skip that block) and writes an integrated comparison_table.tsv +
    comparison.json + fig_*.png comparing classification (classified %, species recovered,
    diversity, genome-length bias) and assembly (contigs, N50, reference breadth via minimap2,
    read concordance). Defaults to config/cross_platform.manifest.tsv ->
    results/experiments/cross_platform_comparison/. Needs minimap2/samtools on PATH for the
    assembly-recovery metrics. Set paper=True to also write an IMRaD comparison manuscript
    (comparison_paper.{tex,pdf}; needs pdflatex). Returns the output directory.
    """
    from metagx import compare
    try:
        out = compare.run(manifest=manifest or None, outdir=outdir or None, paper=paper)
    except FileNotFoundError as e:
        return f"error: {e}"
    return json.dumps({"outdir": out}, indent=2)


@mcp.tool()
def generate_report(config_path: str = "config.yaml", fmt: str = "md") -> str:
    """Generate provenance + a Methods paragraph + a full report from a finished run.

    fmt in {md, latex, pdf} (latex/pdf need pandoc). Captures tool versions, exact
    commands, database identity, and QC/classification metrics into manifest.json, and
    writes methods.md (ready to paste into a paper) + report.md under
    results/<project>/report/. Returns the written paths and the Methods text.
    """
    import yaml
    try:
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
    except FileNotFoundError:
        return f"error: {config_path} not found"
    result = report.generate(cfg, fmt=fmt)
    return json.dumps(result, indent=2)


@mcp.tool()
def generate_paper(config_path: str = "config.yaml", compile_pdf: bool = True) -> str:
    """Write a full IMRaD manuscript (Introduction/Methods/Results/Discussion) from a run.

    Elaborates the interview-captured design + the run's result files into a structured
    LaTeX paper and compiles it to PDF with pdflatex (if installed), under
    results/<project>/report/paper.{tex,pdf}. Every number is read back from the results;
    interpretation is framed as caveat-aware discussion to refine, not fabricated claims.
    Returns the written paths and whether the PDF compiled. Run the pipeline (and ideally
    generate_report) first.
    """
    import yaml
    try:
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
    except FileNotFoundError:
        return f"error: {config_path} not found"
    result = paper.generate(cfg, compile_pdf=compile_pdf)
    return json.dumps(result, indent=2)


@mcp.tool()
def get_recommendations(
    config_path: str = "",
    tool: str = "kraken2",
    platform: str = "illumina",
    param: str = "confidence",
) -> str:
    """Evidence-based parameter suggestions for enabled tools or one tool/platform/param.

    Pass config_path for full multi-tool routing (QC, Bracken read length, kraken2, optional
    modules). Otherwise returns a single param recommendation from registries + evidence/.
    """
    import yaml

    if config_path:
        try:
            with open(config_path) as fh:
                cfg = yaml.safe_load(fh)
        except FileNotFoundError:
            return f"error: {config_path} not found"
        return json.dumps(tool_advisor.recommend_config(cfg), indent=2)
    return json.dumps(evidence_pack.recommend(tool, platform, param=param), indent=2)


@mcp.tool()
def advise_run(
    config_path: str = "config.yaml",
    write_outputs: bool = True,
    record_history: bool = False,
) -> str:
    """Post-run advisor: classification metrics, warnings, and next-config hints.

    Reads finished results under results/<project>/ and returns rules-first suggestions.
    Optionally writes results/<project>/advisor/advisor.json and appends to history.jsonl.
    """
    import yaml

    try:
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
    except FileNotFoundError:
        return f"error: {config_path} not found"
    analysis = advise.analyze(cfg)
    paths = advise.write_advisor_outputs(cfg, analysis) if write_outputs else {}
    if record_history:
        history.record_from_run(cfg, config_path, analysis, success=True, returncode=0)
    return json.dumps({"analysis": analysis, "paths": paths}, indent=2)


@mcp.tool()
def get_run_history(limit: int = 20, best_metric: str = "") -> str:
    """List recent pipeline trials from .metagx/history.jsonl.

    Set best_metric (e.g. mean_percent_classified) to return only the top trial by that metric.
    """
    if best_metric:
        best = history.best_trial(metric=best_metric)
        return json.dumps(best or {}, indent=2)
    return json.dumps(history.read_entries(limit=limit), indent=2)


@mcp.tool()
def sync_tool_help(tool: str = "") -> str:
    """Diff live tool --help output against parameter registries (maintainer drift check)."""
    if tool:
        return json.dumps(sync_help.diff_registry(tool), indent=2)
    return json.dumps(sync_help.sync_all(), indent=2)


@mcp.tool()
def get_catalog() -> str:
    """Index of tools, evidence files, and workflow scripts."""
    return json.dumps(catalog.build_catalog(), indent=2)


# --------------------------------------------------------------------------- #
# HTTP wrapper for web agents (ChatGPT Actions, Gemini webhooks, Perplexity).  #
# Optional: only the MCP stdio surface (above) is needed for Claude Desktop /  #
# Cursor. FastAPI is imported lazily so missing the `serve` HTTP extras does    #
# not break the primary MCP surface (`pip install metagx[serve]` adds it).      #
# --------------------------------------------------------------------------- #
try:
    from fastapi import FastAPI
    from pydantic import BaseModel
    _HAVE_FASTAPI = True
except ModuleNotFoundError:          # MCP-only install — HTTP surface unavailable.
    _HAVE_FASTAPI = False

if _HAVE_FASTAPI:
    app = FastAPI(title="metagx universal bioinformatics agent server")

    _BUILD_PARAMS = set(inspect.signature(config_builder.build_config).parameters)

    class BuildRequest(BaseModel):
        # Accept ANY build_config section so web agents reach the full feature set
        # (Tier 2/3 tools included); the handler filters to build_config's params.
        # Keeps this surface in lockstep with the builder — no per-tool edits here.
        model_config = {"extra": "allow"}
        samples: Any
        db: Dict[str, str]
        project: str = "run"
        outdir: str = "results"
        threads: int = 8
        preset: Optional[str] = None
        modules: Optional[Dict[str, bool]] = None
        sweep: Optional[Dict[str, Any]] = None
        run: bool = False
        executor: Optional[str] = None   # HPC backend when run=True (see list_schedulers)

    @app.get("/api/v1/tools")
    def http_tools():
        return {t: registry.load_registry(t) for t in registry.list_tools()}

    @app.get("/api/v1/presets")
    def http_presets():
        return presets.describe_presets()

    @app.post("/api/v1/report")
    def http_report(config_path: str = "config.yaml", fmt: str = "md"):
        import yaml
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh)
        return report.generate(cfg, fmt=fmt)

    @app.get("/api/v1/interview/{tool}")
    def http_interview(tool: str, max_tier: int = 2):
        return registry.interview_spec(tool, max_tier=max_tier)

    @app.post("/api/v1/build-and-run")
    def http_build_and_run(req: BuildRequest):
        payload = {k: v for k, v in req.model_dump(exclude={"run", "executor"}).items()
                   if k in _BUILD_PARAMS}
        try:
            cfg = config_builder.build_config(**payload)
        except registry.ValidationError as e:
            return {"status": "error", "error": str(e)}
        path = config_builder.write_config(cfg)
        if not req.run:
            return {"status": "config_written", "path": path, "config": cfg}
        profile = None
        if req.executor:
            try:
                profile = schedulers.profile_path(req.executor)
            except (KeyError, FileNotFoundError) as e:
                return {"status": "error", "error": str(e)}
        proc = runner.run(config=path, profile=profile)
        return {
            "status": "success" if proc.returncode == 0 else "error",
            "returncode": proc.returncode,
            "log": (proc.stdout or "")[-3000:] + (proc.stderr or "")[-3000:],
        }

    app.mount("/mcp", mcp.streamable_http_app())


if __name__ == "__main__":
    mcp.run()
