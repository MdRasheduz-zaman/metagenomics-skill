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

from . import config_builder, dbbuild, formats, paper, presets, registry, report, runner


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
    _print_json(registry.interview_spec(args.tool, max_tier=args.tier))
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
    result = dbbuild.build_db(
        genomes=args.genomes, db_dir=args.db,
        read_length=lengths if len(lengths) > 1 else lengths[0],
        threads=args.threads, run=not args.dry_run,
    )
    _print_json(result)
    return 0 if result.get("ok", not result.get("ran", False)) else 1


def cmd_readlen(args) -> int:
    _print_json({f: formats.estimate_read_length(f) for f in args.files})
    return 0


def cmd_build_cat_db(args) -> int:
    result = dbbuild.build_cat_db(genomes=args.genomes, db_dir=args.db,
                                  taxonomy_dir=args.taxonomy, run=not args.dry_run)
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


def cmd_run(args) -> int:
    profile = args.profile
    if args.slurm and not profile:
        profile = os.path.join(os.path.dirname(runner.workflow_path()), "profiles", "slurm")
    proc = runner.run(config=args.config, cores=args.cores, dry_run=args.dry_run,
                      use_conda=args.use_conda or args.slurm, profile=profile)
    sys.stdout.write(proc.stdout)
    sys.stderr.write(proc.stderr)
    return proc.returncode


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
    compare.run(manifest=args.manifest, outdir=args.outdir)
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
    sp.set_defaults(func=cmd_interview)

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
    sp.add_argument("--slurm", action="store_true",
                    help="submit jobs via SLURM using the bundled profile (edit partition/account)")
    sp.add_argument("--profile", default=None,
                    help="path to a Snakemake profile dir (overrides --slurm)")
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("results", help="print result summary JSON files")
    sp.add_argument("--config", default="config.yaml")
    sp.set_defaults(func=cmd_results)

    sp = sub.add_parser("build-db", help="build a custom kraken2+Bracken db from a genomes FASTA")
    sp.add_argument("--genomes", required=True, help="FASTA of reference genomes")
    sp.add_argument("--db", required=True, help="output database directory")
    sp.add_argument("--read-length", default="150",
                    help="Bracken build length(s); comma-separated for multiple (e.g. 150,1000)")
    sp.add_argument("--threads", type=int, default=4)
    sp.add_argument("--dry-run", action="store_true", help="write taxonomy/library, print commands, don't build")
    sp.set_defaults(func=cmd_build_db)

    sp = sub.add_parser("readlen", help="estimate read-length stats (to pick Bracken length)")
    sp.add_argument("files", nargs="+", help="FASTA/FASTQ (±gz) read files")
    sp.set_defaults(func=cmd_readlen)

    sp = sub.add_parser("build-cat-db", help="build a custom CAT contig-annotation db from genomes")
    sp.add_argument("--genomes", required=True, help="FASTA of reference genomes")
    sp.add_argument("--db", required=True, help="output CAT build directory")
    sp.add_argument("--taxonomy", required=True, help="dir with names.dmp + nodes.dmp (e.g. the kraken2 db's taxonomy/)")
    sp.add_argument("--dry-run", action="store_true")
    sp.set_defaults(func=cmd_build_cat_db)

    sp = sub.add_parser("compare", help="compare the same sample across sequencing platforms (manifest TSV)")
    sp.add_argument("--manifest", default=None,
                    help="TSV: label,platform,platform_class,kreport,reads,contigs,reference "
                         "(default: config/cross_platform.manifest.tsv)")
    sp.add_argument("--outdir", default=None,
                    help="output dir (default: results/experiments/cross_platform_comparison)")
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
