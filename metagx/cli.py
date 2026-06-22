"""metagx command-line interface.

A universal surface for any agent that can run a shell command (Codex, Ollama,
plain terminals) and a debugging aid for the MCP/skill paths. Every subcommand
maps onto the same core functions used by the MCP server.

    metagx tools                       # list tools that have a registry
    metagx params kraken2              # full parameter registry (JSON)
    metagx interview kraken2 --tier 2  # questions an LLM should ask (JSON)
    metagx build-config answers.json   # validate answers -> config.yaml
    metagx validate config.yaml        # validate an existing config
    metagx run --config config.yaml    # run the Snakemake workflow
    metagx results --config config.yaml  # print result summaries (JSON)
"""

from __future__ import annotations

import argparse
import glob
import inspect
import json
import os
import sys

import yaml

from . import (
    advise,
    catalog,
    config_builder,
    dbbuild,
    dbfetch,
    doctor,
    evidence_pack,
    formats,
    history,
    paper,
    presets,
    probe,
    registry,
    report,
    runner,
    scaffold,
    schedulers,
    sync_help,
    tool_advisor,
)


def _print_json(obj) -> None:
    print(json.dumps(obj, indent=2))


def cmd_tools(_args) -> int:
    _print_json(
        {t: registry.load_registry(t)["description"].strip() for t in registry.list_tools()}
    )
    return 0


def cmd_params(args) -> int:
    _print_json(registry.load_registry(args.tool))
    return 0


def cmd_interview(args) -> int:
    context = json.loads(args.context) if args.context else {}
    if args.probe:
        with open(args.probe) as fh:
            context.update((json.load(fh).get("context") or {}))
    if args.goal:
        context["goal"] = args.goal
    _print_json(registry.interview_spec(args.tool, max_tier=args.tier, context=context or None))
    return 0


def cmd_probe(args) -> int:
    res = probe.run(args.samples, max_reads=args.max_reads, max_samples=args.max_samples,
                    out=args.out, assume_yes=args.yes, host_index=args.host_index)
    if not res.get("measured"):
        print(res.get("reason", "probe did not run"), file=sys.stderr)
    _print_json(res.get("context") if args.context_only else res)
    return 0


def cmd_presets(args) -> int:
    _print_json(presets.describe_presets())
    return 0


def cmd_preset(args) -> int:
    try:
        _print_json(presets.get_preset_config(args.name))
    except KeyError as e:
        print(str(e), file=sys.stderr)
        return 2
    return 0


def cmd_build_db(args) -> int:
    lengths = [int(x) for x in str(args.read_length).split(",") if x.strip()]
    if args.strategy:  # multi-strategy builder (standard / custom-* / spike-in)
        source = args.source or args.genomes
        if args.strategy in ("custom-fasta", "custom-folder", "spike-in") and not source:
            sys.stderr.write(f"build-db --strategy {args.strategy} needs --source (or --genomes)\n")
            return 2
        if args.strategy in ("standard", "spike-in") and not args.libraries:
            sys.stderr.write(f"build-db --strategy {args.strategy} needs --libraries\n")
            return 2
        result = dbbuild.build_database(
            db_dir=args.db, strategy=args.strategy, taxonomy=args.taxonomy,
            libraries=args.libraries, source=source,
            read_lengths=lengths, threads=args.threads,
            use_ftp=args.use_ftp, run=not args.dry_run,
        )
    else:  # legacy custom-fasta synthetic path
        if not args.genomes:
            sys.stderr.write("build-db needs --genomes (or --strategy with --source/--libraries)\n")
            return 2
        result = dbbuild.build_db(
            genomes=args.genomes, db_dir=args.db,
            read_length=lengths if len(lengths) > 1 else lengths[0],
            threads=args.threads, run=not args.dry_run,
        )
    _print_json(result)
    return 0 if result.get("ok", not result.get("ran", False)) else 1


def cmd_fetch_db(args) -> int:
    if args.list:
        _print_json(dbfetch.describe())
        return 0
    if args.name not in dbfetch.INDICES and not args.url:
        print(f"unknown index '{args.name}'. Run `metagx fetch-db --list` to see options, "
              f"or pass --url for a custom prebuilt index.", file=sys.stderr)
        return 2
    result = dbfetch.fetch(name=args.name, db_dir=args.dir, url=args.url,
                           run=not args.dry_run, force=args.force)
    _print_json(result)
    if not args.dry_run and result.get("ran") and not result.get("ok"):
        return 1
    return 0


def cmd_readlen(args) -> int:
    _print_json({f: formats.estimate_read_length(f) for f in args.files})
    return 0


def cmd_build_cat_db(args) -> int:
    result = dbbuild.build_cat_db(genomes=args.genomes, db_dir=args.db,
                                  taxonomy_dir=args.taxonomy, run=not args.dry_run)
    _print_json(result)
    return 0 if result.get("ok", not result.get("ran", False)) else 1


def cmd_build_kaiju_db(args) -> int:
    result = dbbuild.build_kaiju_db(genomes=args.genomes, db_dir=args.db,
                                    taxonomy_dir=args.taxonomy, threads=args.threads,
                                    run=not args.dry_run)
    _print_json(result)
    return 0 if result.get("ok", not result.get("ran", False)) else 1


def cmd_report(args) -> int:
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    result = report.generate(cfg, fmt=args.format)
    _print_json(result["paths"])
    return 0


def cmd_paper(args) -> int:
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    result = paper.generate(cfg, compile_pdf=not args.no_pdf)
    _print_json(result["paths"])
    if not args.no_pdf and not result["compiled"]:
        print("note: pdflatex not found or compile failed — wrote paper.tex only "
              "(install a LaTeX distribution to get the PDF).", file=sys.stderr)
    return 0


def cmd_build_config(args) -> int:
    with open(args.answers) as fh:
        answers = json.load(fh) if args.answers.endswith(".json") else yaml.safe_load(fh)
    try:
        cfg = config_builder.build_config(**answers)
    except registry.ValidationError as e:
        print(f"validation error: {e}", file=sys.stderr)
        return 2
    path = config_builder.write_config(cfg, args.out)
    print(f"wrote {path}")
    return 0


def cmd_validate(args) -> int:
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    # Round-trip through build_config to reuse ALL validation. Forward every
    # build_config parameter that appears in the config by name so this stays in
    # lockstep with the builder — adding a new tool section (registry-as-truth)
    # needs no edit here and cannot silently skip validation.
    params = inspect.signature(config_builder.build_config).parameters
    kwargs = {name: cfg[name] for name in params if name in cfg}
    kwargs.setdefault("samples", cfg.get("samples"))
    kwargs.setdefault("db", cfg.get("db", {}))
    try:
        config_builder.build_config(**kwargs)
    except (registry.ValidationError, KeyError, TypeError) as e:
        print(f"invalid: {e}", file=sys.stderr)
        return 2
    print("config is valid")
    return 0


def cmd_doctor(args) -> int:
    """Environment preflight: detect arch/env hazards and missing tools/DB, print remedies."""
    db_paths = None
    if args.config and os.path.isfile(args.config):
        with open(args.config) as fh:
            cfg = yaml.safe_load(fh) or {}
        db_paths = cfg.get("db") or None
    checks = doctor.run(db_paths=db_paths)
    if args.json:
        _print_json([c.as_dict() for c in checks])
    else:
        print(doctor.format_report(checks))
    n_fail = sum(1 for c in checks if c.status == "fail")
    n_warn = sum(1 for c in checks if c.status == "warn")
    if n_fail or (args.strict and n_warn):
        return 1
    return 0


def cmd_schedulers(_args) -> int:
    """List the bundled HPC scheduler backends (for `metagx run --executor`)."""
    _print_json(schedulers.describe())
    return 0


def _resolve_profile(args) -> str | None:
    """Pick the Snakemake profile dir from --profile / --executor / --slurm.

    Precedence: an explicit external --profile wins; then --executor <name>;
    then the legacy --slurm alias. Returns None for a plain local run.
    """
    if args.profile:
        return args.profile
    name = args.executor
    if not name and args.slurm:
        name = "slurm"          # legacy alias, kept for back-compat
    if not name or name == "local":
        # `local` still routes through its bundled profile (thread/mem caps);
        # truly plain local runs (no --executor) get None.
        if name == "local":
            return schedulers.profile_path("local")
        return None
    try:
        return schedulers.profile_path(name)
    except (KeyError, FileNotFoundError) as e:
        raise SystemExit(f"metagx run: {e}")


def cmd_run(args) -> int:
    profile = _resolve_profile(args)
    try:
        proc = runner.run(config=args.config, cores=args.cores, dry_run=args.dry_run,
                          use_conda=args.use_conda, profile=profile)
    except runner.CondaFrontendError as e:
        sys.stderr.write(f"error: {e}\n")
        return 1
    sys.stdout.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    if proc.returncode == 0 and not args.dry_run and not args.no_history:
        try:
            with open(args.config) as fh:
                cfg = yaml.safe_load(fh)
            analysis = advise.analyze(cfg)
            if not args.no_advisor:
                advise.write_advisor_outputs(cfg, analysis)
            history.record_from_run(
                cfg, args.config, analysis,
                success=True, returncode=0, path=args.history_file,
            )
        except OSError:
            pass
    return proc.returncode


def cmd_recommend(args) -> int:
    cfg = None
    if args.config:
        with open(args.config) as fh:
            cfg = yaml.safe_load(fh)
        if args.all or args.tool is None:
            _print_json(tool_advisor.recommend_config(cfg))
            return 0
    plat = args.platform
    if not plat and cfg:
        plats = advise.platforms_from_config(cfg)
        plat = plats[0] if plats else "illumina"
    plat = plat or "illumina"
    _print_json(evidence_pack.recommend(args.tool or "kraken2", plat, param=args.param))
    return 0


def cmd_advise(args) -> int:
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    analysis = advise.analyze(cfg)
    paths = {}
    if args.write:
        paths = advise.write_advisor_outputs(cfg, analysis)
    if args.record_history:
        history.record_from_run(
            cfg, args.config, analysis,
            success=True, returncode=0, path=args.history_file,
        )
    out = {"analysis": analysis, "paths": paths}
    _print_json(out)
    return 0


def cmd_history(args) -> int:
    if args.best:
        best = history.best_trial(path=args.history_file, metric=args.best)
        if not best:
            print("no history entries with metrics", file=sys.stderr)
            return 1
        _print_json(best)
        return 0
    _print_json(history.read_entries(path=args.history_file, limit=args.limit))
    return 0


def cmd_sync_help(args) -> int:
    if args.tool:
        _print_json(sync_help.diff_registry(args.tool))
    else:
        _print_json(sync_help.sync_all())
    return 0


def cmd_scaffold(args) -> int:
    res = scaffold.scaffold(args.command, name=args.name)
    if not res.get("ok"):
        print(f"scaffold failed: {res.get('error')}", file=sys.stderr)
        return 1
    if res.get("registry_exists"):
        print(f"note: a curated registry for '{res['tool']}' already exists — "
              "this stub is for diffing/extending, do not blindly overwrite it.", file=sys.stderr)
    if args.out:
        if os.path.exists(args.out):
            print(f"refusing to overwrite existing {args.out}", file=sys.stderr)
            return 1
        with open(args.out, "w") as fh:
            fh.write(res["yaml"])
        print(f"wrote {args.out} ({res['tool']}, version: {res.get('version')})", file=sys.stderr)
    else:
        sys.stdout.write(res["yaml"])
    return 0


def cmd_catalog(_args) -> int:
    _print_json(catalog.build_catalog())
    return 0


def cmd_results(args) -> int:
    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    outdir = cfg.get("outdir", "results")
    summaries = {}
    for path in sorted(glob.glob(os.path.join(outdir, "**", "*summary*.json"), recursive=True)):
        with open(path) as fh:
            summaries[path] = json.load(fh)
    if not summaries:
        print(f"no result summaries found under {outdir}/", file=sys.stderr)
        return 1
    _print_json(summaries)
    return 0


def cmd_compare(args) -> int:
    from . import compare  # lazy: pulls matplotlib/pandas, only when actually comparing
    compare.run(manifest=args.manifest, outdir=args.outdir, paper=args.paper)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="metagx", description="Schema-driven metagenomics pipeline.")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("tools", help="list tools with a parameter registry").set_defaults(func=cmd_tools)

    sp = sub.add_parser("params", help="dump a tool's full registry as JSON")
    sp.add_argument("tool")
    sp.set_defaults(func=cmd_params)

    sp = sub.add_parser("interview", help="questions an LLM should ask for a tool")
    sp.add_argument("tool")
    sp.add_argument("--tier", type=int, default=2, help="max tier: 1 core, 2 common, 3 advanced")
    sp.add_argument("--goal", default=None,
                    help="experimental goal (e.g. strain_resolved); may promote quiet params")
    sp.add_argument("--context", default=None,
                    help='JSON of goal/data facts for promote_when, e.g. \'{"estimated_bases": 6e10}\'')
    sp.add_argument("--probe", default=None,
                    help="probe.json from `metagx probe`; loads MEASURED context for promotion")
    sp.set_defaults(func=cmd_interview)

    sp = sub.add_parser("probe",
                        help="measure read stats from your samples (local, consent-gated) -> context")
    sp.add_argument("--samples", required=True, help="sample sheet TSV (or inline-validated path)")
    sp.add_argument("--max-reads", type=int, default=100_000, help="reads sampled per file (cap)")
    sp.add_argument("--max-samples", type=int, default=None, help="cap number of samples scanned")
    sp.add_argument("--host-index", default=None,
                    help="host reference (FASTA/minimap2 index); measures host fraction if minimap2 on PATH")
    sp.add_argument("--out", default=None, help="dir to write probe.json + probe.md")
    sp.add_argument("--yes", action="store_true",
                    help="grant + remember LOCAL probe consent (data never leaves the machine)")
    sp.add_argument("--context-only", action="store_true",
                    help="print just the context dict (for piping into `interview --probe`)")
    sp.set_defaults(func=cmd_probe)

    sub.add_parser("presets", help="list workflow presets with descriptions").set_defaults(func=cmd_presets)

    sp = sub.add_parser("preset", help="dump a preset's template config as JSON")
    sp.add_argument("name")
    sp.set_defaults(func=cmd_preset)

    sp = sub.add_parser("build-config", help="validate interview answers into config.yaml")
    sp.add_argument("answers", help="JSON/YAML file of build_config kwargs")
    sp.add_argument("--out", default="config.yaml")
    sp.set_defaults(func=cmd_build_config)

    sp = sub.add_parser("validate", help="validate an existing config.yaml")
    sp.add_argument("config")
    sp.set_defaults(func=cmd_validate)

    sp = sub.add_parser("run", help="run the Snakemake workflow")
    sp.add_argument("--config", default="config.yaml")
    sp.add_argument("--cores", default="all")
    sp.add_argument("--dry-run", action="store_true")
    sp.add_argument("--use-conda", action="store_true",
                    help="auto-provision per-rule tools via conda (workflow/envs/)")
    sp.add_argument("--executor", default=None, metavar="NAME",
                    help="HPC scheduler backend: "
                         f"{'|'.join(schedulers.list_schedulers())} "
                         "(see `metagx schedulers`; edit the bundled profile first)")
    sp.add_argument("--slurm", action="store_true",
                    help="alias for --executor slurm (back-compat)")
    sp.add_argument("--profile", default=None,
                    help="path to a custom Snakemake profile dir (overrides --executor)")
    sp.add_argument("--no-history", action="store_true",
                    help="do not append to .metagx/history.jsonl after a successful run")
    sp.add_argument("--no-advisor", action="store_true",
                    help="skip writing results/<project>/advisor/ after a successful run")
    sp.add_argument("--history-file", default=None,
                    help="history log path (default: .metagx/history.jsonl)")
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("recommend", help="evidence-based parameter suggestions for a platform")
    sp.add_argument("--tool", default=None, help="single-tool mode (default: all tools when --config)")
    sp.add_argument("--param", default="confidence")
    sp.add_argument("--platform", default=None,
                    help="illumina|ont|pacbio_hifi|pacbio_clr (default: from --config or illumina)")
    sp.add_argument("--config", default=None, help="full multi-tool recommendations from a run config")
    sp.add_argument("--all", action="store_true",
                    help="with --config: recommend all enabled modules/tools (default when --config alone)")
    sp.set_defaults(func=cmd_recommend)

    sp = sub.add_parser("advise", help="post-run advisor: metrics, warnings, next-config hints")
    sp.add_argument("--config", default="config.yaml")
    sp.add_argument("--write", action="store_true",
                    help="write advisor.json, trial_log.md, next_config.suggested.yaml")
    sp.add_argument("--record-history", action="store_true",
                    help="append this analysis to history.jsonl")
    sp.add_argument("--history-file", default=None)
    sp.set_defaults(func=cmd_advise)

    sp = sub.add_parser("history", help="list prior run trials (.metagx/history.jsonl)")
    sp.add_argument("--limit", type=int, default=20)
    sp.add_argument("--history-file", default=None)
    sp.add_argument("--best", default=None, metavar="METRIC",
                    help="show best trial by metric (e.g. mean_percent_classified)")
    sp.set_defaults(func=cmd_history)

    sp = sub.add_parser("sync-help", help="diff live tool --help against parameter registries")
    sp.add_argument("--tool", default=None, help="single tool; default: all registries")
    sp.set_defaults(func=cmd_sync_help)

    sp = sub.add_parser("scaffold",
                        help="generate a capability-complete registry stub from a tool's --help")
    sp.add_argument("command", help="the executable to probe (e.g. flye)")
    sp.add_argument("--name", default=None,
                    help="registry tool name (default: the command's first token)")
    sp.add_argument("--out", default=None,
                    help="write YAML to this path (refuses to overwrite); default: stdout")
    sp.set_defaults(func=cmd_scaffold)

    sub.add_parser("catalog", help="index of tools, evidence, and workflow scripts").set_defaults(
        func=cmd_catalog
    )

    sub.add_parser("schedulers", help="list HPC scheduler backends for `run --executor`"
                   ).set_defaults(func=cmd_schedulers)

    sp = sub.add_parser("doctor",
                        help="preflight: detect arch/env/tool/DB hazards and print remedies")
    sp.add_argument("--config", default=None,
                    help="also check the kraken2 DB referenced by this config.yaml")
    sp.add_argument("--json", action="store_true", help="emit checks as JSON")
    sp.add_argument("--strict", action="store_true", help="exit non-zero on warnings too")
    sp.set_defaults(func=cmd_doctor)

    sp = sub.add_parser("results", help="print result summary JSON files")
    sp.add_argument("--config", default="config.yaml")
    sp.set_defaults(func=cmd_results)

    sp = sub.add_parser("build-db", help="build a kraken2+Bracken db (custom genomes, NCBI libraries, or spike-in)")
    sp.add_argument("--genomes", help="FASTA of reference genomes (custom-fasta; legacy synthetic path)")
    sp.add_argument("--db", required=True, help="output database directory")
    sp.add_argument("--read-length", default="150",
                    help="Bracken build length(s); comma-separated for multiple (e.g. 150,1000)")
    sp.add_argument("--threads", type=int, default=4)
    sp.add_argument("--strategy", choices=["standard", "custom-fasta", "custom-folder", "spike-in"],
                    help="db.build strategy; given => use the multi-strategy builder")
    sp.add_argument("--taxonomy", choices=["real", "synthetic"], default="real",
                    help="taxonomy for custom/spike-in (real=NCBI, synthetic=flat); default real")
    sp.add_argument("--libraries", help="NCBI libraries for standard/spike-in (e.g. 'viral' or 'bacteria,viral')")
    sp.add_argument("--source", help="FASTA file or folder of FASTAs for custom-*/spike-in")
    sp.add_argument("--no-use-ftp", dest="use_ftp", action="store_false",
                    help="use rsync instead of FTP for NCBI downloads (rsync is deprecated; default FTP)")
    sp.add_argument("--dry-run", action="store_true", help="write taxonomy/library, print commands, don't build")
    sp.set_defaults(func=cmd_build_db, use_ftp=True)

    sp = sub.add_parser("fetch-db",
                        help="download a prebuilt standard kraken2+Bracken index (onboarding)")
    sp.add_argument("name", nargs="?", default=dbfetch.DEFAULT,
                    help=f"index name (default: {dbfetch.DEFAULT}; see --list)")
    sp.add_argument("--list", action="store_true", help="list curated indices with sizes + URLs")
    sp.add_argument("--dir", default="local_databases/kraken2", help="output database directory")
    sp.add_argument("--url", default=None, help="custom prebuilt-index tarball URL (overrides name)")
    sp.add_argument("--force", action="store_true", help="re-download even if a built index exists")
    sp.add_argument("--dry-run", action="store_true", help="print the plan/command, don't download")
    sp.set_defaults(func=cmd_fetch_db)

    sp = sub.add_parser("readlen", help="estimate read-length stats (to pick Bracken length)")
    sp.add_argument("files", nargs="+", help="FASTA/FASTQ (±gz) read files")
    sp.set_defaults(func=cmd_readlen)

    sp = sub.add_parser("build-cat-db", help="build a custom CAT contig-annotation db from genomes")
    sp.add_argument("--genomes", required=True, help="FASTA of reference genomes")
    sp.add_argument("--db", required=True, help="output CAT build directory")
    sp.add_argument("--taxonomy", required=True, help="dir with names.dmp + nodes.dmp (e.g. the kraken2 db's taxonomy/)")
    sp.add_argument("--dry-run", action="store_true")
    sp.set_defaults(func=cmd_build_cat_db)

    sp = sub.add_parser("build-kaiju-db",
                        help="build a custom Kaiju (protein) db from genomes for the consensus module")
    sp.add_argument("--genomes", required=True, help="FASTA of reference genomes")
    sp.add_argument("--db", required=True, help="output Kaiju db directory (becomes db.kaiju)")
    sp.add_argument("--taxonomy", required=True,
                    help="dir with names.dmp + nodes.dmp (e.g. the kraken2 db's taxonomy/)")
    sp.add_argument("--threads", type=int, default=4)
    sp.add_argument("--dry-run", action="store_true")
    sp.set_defaults(func=cmd_build_kaiju_db)

    sp = sub.add_parser("compare", help="compare the same sample across sequencing platforms (manifest TSV)")
    sp.add_argument("--manifest", default=None,
                    help="TSV: label,platform,platform_class,kreport,reads,contigs,reference "
                         "(default: config/cross_platform.manifest.tsv)")
    sp.add_argument("--outdir", default=None,
                    help="output dir (default: results/experiments/cross_platform_comparison)")
    sp.add_argument("--paper", action="store_true",
                    help="also write an IMRaD comparison manuscript (comparison_paper.{tex,pdf}; needs pdflatex)")
    sp.set_defaults(func=cmd_compare)

    sp = sub.add_parser("report", help="generate provenance manifest + Methods + report")
    sp.add_argument("--config", default="config.yaml")
    sp.add_argument("--format", default="md", choices=["md", "latex", "pdf"])
    sp.set_defaults(func=cmd_report)

    sp = sub.add_parser("paper", help="write a full IMRaD manuscript (LaTeX -> PDF via pdflatex)")
    sp.add_argument("--config", default="config.yaml")
    sp.add_argument("--no-pdf", action="store_true", help="write paper.tex only, skip pdflatex")
    sp.set_defaults(func=cmd_paper)

    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
