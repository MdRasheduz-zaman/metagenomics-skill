"""metagx doctor — environment preflight.

Turns the bioconda landmines that used to live as tribal knowledge (the abricate
samtools-0.1.x downgrade, missing tools, a missing database) into machine-checked
diagnostics. Each check yields a ``Check`` with a status and, when something is wrong, the
*exact* remedy — so a stranger on their own machine gets steered, not stranded.

Statuses:
  ok    — verified good.
  info  — context, no action needed.
  warn  — works, but a known footgun is armed (e.g. a tool present but below its version floor).
  fail  — will break a real run; remedy provided.

`run()` returns a list of Checks; `cli()` prints them and exits non-zero if any failed
(or, with --strict, if any warned).
"""

from __future__ import annotations

import os
import platform
import re
import shutil
import subprocess
from dataclasses import dataclass
from typing import Dict, List, Optional

from . import registry, report, runner

# Minimum acceptable (major, minor) per load-bearing tool — the single source of truth,
# mirrored from environment.yml. tests/test_tool_versions.py imports this so the floors
# can't drift between the doctor and the test. samtools is the regression guard: the
# abricate/mapdamage2 dependency chain can silently drag in 0.1.x, whose `sort -o` breaks
# the whole pipeline.
VERSION_FLOORS: Dict[str, tuple] = {
    "kraken2": (2, 1),
    "bracken": (2, 9),
    "fastp": (0, 23),
    "megahit": (1, 2),
    "flye": (2, 9),
    "minimap2": (2, 26),
    "samtools": (1, 18),
    "metabat2": (2, 15),
    "diamond": (2, 0),
    "mafft": (7, 4),
    "checkv": (1, 0),
    "genomad": (1, 7),
    "kaiju": (1, 9),
    "multiqc": (1, 0),
    "mapDamage": (2, 2),
}

# The tools a *default* run (QC → classify → abundance → assembly → bin → reconcile) needs.
# Their absence is a fail, not an info; the rest of VERSION_FLOORS are optional modules.
CORE_TOOLS = ["kraken2", "bracken", "fastp", "megahit", "minimap2", "samtools", "metabat2"]

_OK, _INFO, _WARN, _FAIL = "ok", "info", "warn", "fail"


@dataclass
class Check:
    name: str
    status: str
    message: str
    remedy: Optional[str] = None

    def as_dict(self) -> Dict[str, Optional[str]]:
        return {"name": self.name, "status": self.status,
                "message": self.message, "remedy": self.remedy}


def _parse_xy(version_str: Optional[str]):
    m = re.search(r"(\d+)\.(\d+)", version_str or "")
    return (int(m.group(1)), int(m.group(2))) if m else None


def check_platform() -> List[Check]:
    sys_name, machine = platform.system(), platform.machine()
    out = [Check("platform", _INFO, f"{sys_name} / {machine} / Python {platform.python_version()}")]
    return out


def check_workflow() -> Check:
    try:
        p = runner.workflow_path()
    except FileNotFoundError as e:
        return Check("workflow", _FAIL, "Snakemake workflow not found.", remedy=str(e))
    wf = os.path.dirname(p)
    missing = [s for s in ("rules", "scripts", "envs") if not os.path.isdir(os.path.join(wf, s))]
    if missing:
        return Check("workflow", _FAIL,
                     f"Snakefile found but incomplete (missing {', '.join(missing)}).",
                     remedy="Reinstall metagx so the full workflow/ tree ships "
                            "(`pip install .` from the repo).")
    return Check("workflow", _OK, f"workflow resolves: {p}")


def check_tools(floors: Optional[Dict[str, tuple]] = None,
                core: Optional[List[str]] = None) -> List[Check]:
    floors = floors or VERSION_FLOORS
    core = core or CORE_TOOLS
    captured = report.tool_versions(sorted(floors))
    checks: List[Check] = []
    for tool in sorted(floors):
        raw = captured.get(tool, "not found on PATH")
        is_core = tool in core
        if raw == "not found on PATH":
            if is_core:
                checks.append(Check(
                    f"tool:{tool}", _FAIL, f"{tool} (core) not on PATH.",
                    remedy="Install the core stack: get the spec with `metagx env-file --write` "
                           "(a wheel install has no environment.yml in your cwd), then "
                           "`conda env create -f environment.yml` (or `mamba`), and activate it."))
            else:
                checks.append(Check(
                    f"tool:{tool}", _INFO,
                    f"{tool} not on PATH (optional module; provisioned via `metagx run --use-conda`)."))
            continue
        got = _parse_xy(raw)
        floor = floors[tool]
        if got is None:
            checks.append(Check(f"tool:{tool}", _OK, f"{tool} present ({raw}; version not machine-readable)."))
        elif got < floor:
            checks.append(Check(
                f"tool:{tool}", _FAIL,
                f"{tool} {got[0]}.{got[1]} is below the required floor {floor[0]}.{floor[1]} (got {raw!r}).",
                remedy=(("samtools 0.1.x is the abricate/mapdamage2 downgrade trap — keep abricate "
                         "out of the core env (it lives in workflow/envs/amr.yaml, use --use-conda) "
                         "and re-pin `samtools >=1.18` in environment.yml, then reinstall.")
                        if tool == "samtools"
                        else f"Upgrade {tool} to >= {floor[0]}.{floor[1]} (re-pin in environment.yml and reinstall).")))
        else:
            checks.append(Check(f"tool:{tool}", _OK, f"{tool} {got[0]}.{got[1]} OK ({raw})."))
    return checks


def check_active_tools(cfg: Optional[Dict] = None) -> List[Check]:
    """Presence-check the tools THIS config actually activates that aren't in VERSION_FLOORS
    (blastn/spades/vsearch/genomad/...). ``check_tools`` only covers the fixed floor list, so a
    config that enables e.g. validation but lacks blastn used to sail through preflight and fail
    mid-run. A tool that lives in the core env is a fail when missing; one that is provisioned via
    ``--use-conda`` (not in the core env) is an info, not a failure."""
    out: List[Check] = []
    if not cfg:
        return out
    active = set(report.active_tools(cfg)) - set(VERSION_FLOORS) - {"snakemake"}
    if not active:
        return out
    core_pkgs = set(report._env_yaml_packages(runner.environment_file_path()))  # noqa: PLC2701
    for tool in sorted(active):
        exe = registry.resolve_command(tool)
        if shutil.which(exe):
            checks_msg = f"{tool} present ({exe})." if exe == tool else f"{tool} present (as {exe})."
            out.append(Check(f"tool:{tool}", _OK, checks_msg))
        elif report.conda_package(tool) in core_pkgs:
            out.append(Check(
                f"tool:{tool}", _FAIL, f"{tool} (enabled by this config, core env) not on PATH.",
                remedy="install the core stack: `metagx env-file --write` then "
                       "`conda env create -f environment.yml` (or mamba), and activate it."))
        else:
            out.append(Check(
                f"tool:{tool}", _INFO,
                f"{tool} not on PATH (enabled by this config; provisioned on first use via "
                f"`metagx run --use-conda`)."))
    return out


def check_bracken_runs() -> Optional[Check]:
    """A present-but-unrunnable bracken binary is worse than a missing one because it passes a
    naive `which` check. Actually invoke it."""
    exe = shutil.which("bracken")
    if not exe:
        return None  # absence is already reported by check_tools
    try:
        p = subprocess.run(["bracken", "-h"], capture_output=True, text=True, timeout=20)
    except (subprocess.SubprocessError, OSError) as e:
        return Check("bracken-runs", _FAIL, f"bracken is on PATH but failed to execute: {e}",
                     remedy="Reinstall bracken (e.g. from environment.yml / bioconda), or use Docker.")
    blob = (p.stdout or "") + (p.stderr or "")
    if p.returncode != 0 and "usage" not in blob.lower() and "bracken" not in blob.lower():
        return Check("bracken-runs", _FAIL,
                     f"bracken is on PATH but `bracken -h` exited {p.returncode} with no usage text.",
                     remedy="Reinstall bracken (e.g. from environment.yml / bioconda), or use Docker.")
    return Check("bracken-runs", _OK, "bracken executes.")


def check_conda_frontend() -> Check:
    """Whether `metagx run --use-conda` (the route to the heavy optional tools) can work."""
    frontend = runner.pick_conda_frontend()
    if shutil.which(frontend) is None:
        return Check("conda-frontend", _WARN,
                     "Neither mamba nor conda is on PATH.",
                     remedy="--use-conda (GTDB-Tk/CheckM2/antiSMASH/AMR/…) needs a conda frontend. "
                            "Install mamba (`conda install -n base -c conda-forge mamba`). "
                            "Not needed for the core kraken2/Bracken pipeline.")
    problem = runner.conda_preflight(frontend)
    if problem:
        return Check("conda-frontend", _WARN, f"{frontend} present but can't drive --use-conda.",
                     remedy=problem)
    return Check("conda-frontend", _OK, f"{frontend} can drive `--use-conda`.")


def check_database(db_paths: Optional[Dict[str, str]] = None) -> Check:
    """A configured kraken2 DB is the #1 real-user blocker. If a config supplies db paths,
    verify they exist; otherwise point at the onboarding command."""
    if not db_paths:
        return Check("database", _INFO,
                     "No database configured (pass --config to check one).",
                     remedy="Get a usable kraken2/Bracken index with `metagx fetch-db --list` "
                            "then `metagx fetch-db <name> --dir <path>`.")
    kdb = db_paths.get("kraken2")
    has_build = bool(db_paths.get("build"))
    if not kdb:
        return Check("database", _WARN, "Config has no db.kraken2 entry.",
                     remedy="Set db.kraken2 to a kraken2 index directory, or run `metagx fetch-db`.")
    info = report.db_info(kdb)
    if not info.get("present"):
        if has_build:  # the DB doesn't exist yet, but db.build will produce it — not a failure
            return Check("database", _INFO,
                         f"kraken2 db {kdb} not present yet — will be built by db.build.",
                         remedy="Run `metagx build-db` (or `metagx run`, which auto-builds it).")
        return Check("database", _FAIL, f"Configured kraken2 db not found: {kdb}",
                     remedy="Download one with `metagx fetch-db standard-8 --dir <path>` "
                            "or build a custom one with `metagx build-db`.")
    if not os.path.isfile(os.path.join(kdb, "hash.k2d")):
        return Check("database", _FAIL, f"{kdb} exists but has no hash.k2d (not a built kraken2 db).",
                     remedy="Finish the build (`metagx build-db`) or re-download (`metagx fetch-db`).")
    size_gb = info.get("size_bytes", 0) / 1e9
    return Check("database", _OK, f"kraken2 db OK: {kdb} ({size_gb:.1f} GB).")


def _source_has_taxid_headers(source: str):
    """Scan a db.build source (FASTA file or folder of FASTAs) and report whether its headers
    carry `kraken:taxid|` tags. Returns True/False, or None if the source can't be read yet
    (e.g. a path that doesn't exist at preflight time) so we don't false-alarm."""
    import glob
    import gzip

    files = []
    if os.path.isdir(source):
        for ext in ("*.fa", "*.fna", "*.fasta", "*.fa.gz", "*.fna.gz", "*.fasta.gz"):
            files += glob.glob(os.path.join(source, ext))
    elif os.path.isfile(source):
        files = [source]
    if not files:
        return None
    for fp in files[:5]:                      # sample a few files; headers are uniform per build
        opener = gzip.open if fp.endswith(".gz") else open
        try:
            with opener(fp, "rt") as fh:
                for line in fh:
                    if line.startswith(">"):
                        if "kraken:taxid|" in line or "taxid|" in line:
                            return True
                        break                 # first header per file is representative enough
        except OSError:
            return None
    return False


def check_db_build(db_paths: Optional[Dict[str, str]] = None) -> List[Check]:
    """When db.build is configured, surface the build tooling, the masking dependency, and
    the air-gapped-HPC download caveat — so a user catches a no-internet compute node before
    a multi-hour job, or tells us up front that their cluster is air-gapped."""
    out: List[Check] = []
    build = (db_paths or {}).get("build")
    if not build:
        return out
    strategy = build.get("strategy", "standard")
    taxonomy = build.get("taxonomy", "real")
    needs_download = strategy in {"standard", "spike-in"} or taxonomy == "real"

    for t in ("kraken2-build", "bracken-build"):
        if shutil.which(t) is None:
            out.append(Check(f"db-build:{t}", _FAIL, f"{t} not on PATH (needed to build the DB).",
                             remedy="Install the core stack (kraken2 + bracken ship these): "
                                    "`conda env create -f environment.yml`."))
        else:
            out.append(Check(f"db-build:{t}", _OK, f"{t} present."))

    # Real-taxonomy custom/spike-in builds need each sequence to carry a real NCBI taxid in its
    # header (kraken:taxid|<id>); otherwise the build "succeeds" but those sequences map nowhere.
    # Catch it here, before a long build, rather than discovering a useless DB afterward.
    src = build.get("source")
    if strategy in {"custom-fasta", "custom-folder", "spike-in"} and taxonomy == "real" and src:
        tagged = _source_has_taxid_headers(src)
        if tagged is False:
            out.append(Check(
                "db-build:taxids", _WARN,
                f"taxonomy: real but the headers in {src} carry no `kraken:taxid|<id>` tag — "
                "those sequences won't map into the NCBI taxonomy.",
                remedy="Tag headers as `>acc|kraken:taxid|<ncbi_taxid> ...`, or use "
                       "taxonomy: synthetic if you only need to detect/quantify these genomes."))
    if needs_download and not build.get("no_masking", False) and shutil.which("dustmasker") is None:
        out.append(Check("db-build:masking", _WARN,
                         "dustmasker (BLAST+) not on PATH — low-complexity masking will fail.",
                         remedy="Install `blast`, or set db.build.no_masking: true (a few more "
                                "false positives, no BLAST+ dependency)."))
    if needs_download:
        kind = "taxonomy + libraries" if strategy in {"standard", "spike-in"} else "taxonomy"
        on = build.get("download_on", "rule")
        out.append(Check(
            "db-build:network", _INFO,
            f"db.build downloads NCBI {kind} (download_on={on}, use_ftp={build.get('use_ftp', True)}). "
            "HPC check: make sure the node that runs the build can reach the internet — many "
            "clusters allow it from compute nodes, some need http_proxy/https_proxy set, a few "
            "air-gap compute nodes (then set db.build.download_on: login, or pre-stage the DB "
            "from a login/data-transfer node)."))
        # NCBI deprecated rsync, so a from-scratch standard build now fetches every genome
        # individually over FTP/wget — fine for small libraries (viral, UniVec_Core), but slow
        # for large ones (bacteria alone is ~15k+ genomes / hours). Prefer the prebuilt index.
        libs = {l.strip() for l in str(build.get("libraries") or "").split(",") if l.strip()}
        big = libs - {"viral", "UniVec_Core", "plasmid"}
        if strategy in {"standard", "spike-in"} and big:
            out.append(Check(
                "db-build:slow-download", _WARN,
                f"building {sorted(big)} from NCBI is slow now — rsync is deprecated, so "
                "kraken2-build fetches genomes one-by-one over FTP (bacteria/nt = hours).",
                remedy="Prefer a prebuilt index: `metagx fetch-db --list` "
                       "(standard-8 ~6GB, standard ~76GB) and set db.kraken2 to it; reserve "
                       "db.build for custom/spike-in or small libraries (viral)."))
    # The aligned BLAST DB (db.build.blast) runs makeblastdb on the SAME genomes; for a big
    # standard library that's a large, slow build + a large second on-disk DB. Flag it so the
    # user sizes for it (the validation reference can be huge — e.g. a full standard library).
    if build.get("blast"):
        libs = {l.strip() for l in str(build.get("libraries") or "").split(",") if l.strip()}
        big = libs - {"viral", "UniVec_Core", "plasmid"}
        if strategy in {"standard", "spike-in"} and big:
            out.append(Check(
                "db-build:blast-cost", _WARN,
                f"db.build.blast will makeblastdb the {sorted(big)} library too — a large, slow "
                "build and a second multi-GB on-disk DB alongside the kraken2 index.",
                remedy="Fine for custom/small/viral DBs. For a big standard library, either accept "
                       "the cost (size your disk), or validate a representative subset, or use "
                       "validate.remote for a handful of sequences."))
    return out


def check_module_dbs(cfg: Optional[Dict] = None) -> List[Check]:
    """Fail-fast when an enabled module needs a reference DB that's unset or empty — with the
    exact `metagx fetch-db --tool` command — so a run doesn't crash mid-pipeline on a missing
    domain/functional DB (genomad/checkv/checkm2/gtdbtk/bakta/...)."""
    out: List[Check] = []
    if not cfg:
        return out
    from . import dbprovision
    db = cfg.get("db", {}) or {}
    for tool, key in dbprovision.needed_dbs(cfg).items():
        spec = dbprovision.SPECS.get(tool, {})
        self_gates = bool(spec.get("self_gates"))
        manual = bool(spec.get("manual"))
        remedy_fetch = (f"download manually — {spec.get('docs')}" if manual
                        else f"`metagx fetch-db --tool {tool} --dir <dir>` ({spec.get('size')})")
        path = db.get(key)
        if not path:
            status = _INFO if self_gates else _FAIL
            msg = (f"{tool} DB not set (db.{key}); its step will be skipped." if self_gates
                   else f"module needs {tool} ({spec.get('needed_by')}) but db.{key} is unset.")
            out.append(Check(f"moduledb:{tool}", status, msg,
                             remedy=f"{remedy_fetch}, then set db.{key}."))
        elif not dbprovision.is_provisioned(tool, path):
            out.append(Check(f"moduledb:{tool}", _FAIL,
                             f"db.{key}={path} has no recognizable {tool} DB files.",
                             remedy=f"{remedy_fetch} (or --use-conda to provision the tool first)."))
        else:
            out.append(Check(f"moduledb:{tool}", _OK, f"{tool} DB present: {path}."))
    return out


def check_config_flags(cfg: Optional[Dict] = None) -> List[Check]:
    """Validate that flags a config sets actually exist in the *installed* tool versions —
    catching a renamed/removed flag (or wrong tool version) before the run, not mid-pipeline.
    Only fires for tools on PATH whose --help parses; otherwise stays silent."""
    out: List[Check] = []
    if not cfg:
        return out
    from . import toollock
    for f in toollock.config_flag_check(cfg):
        out.append(Check(f"toolflag:{f['tool']}", _FAIL, f["message"],
                         remedy=f"check the installed {f['tool']} version, pin it (see "
                                f"`metagx lock`), or move the flag to {f['tool']}.extra_args."))
    return out


def check_param_conflicts(cfg: Optional[Dict] = None) -> List[Check]:
    """Catch flag combinations that each validate fine in isolation but are jointly broken — e.g.
    kraken2 ``use_mpa_style`` set while the Bracken/abundance module is on (the mpa report format
    breaks the kreport parser). Registry-driven: a param declares ``conflicts: [{when:, message:}]``
    and this evaluates them against the run context (module toggles + the tool's set params)."""
    out: List[Check] = []
    if not cfg:
        return out
    # Evaluate against the *effective* module map, not the literal `modules:` block. Several
    # modules (abundance, classify, qc) default ON and run whether or not the config lists them,
    # so a hand-written config that omits `modules:` must still trip a conflict against a
    # default-on module — otherwise the guard misses the very case it exists for.
    from .config_builder import DEFAULT_MODULES
    mods = {**DEFAULT_MODULES, **(cfg.get("modules", {}) or {})}
    ctx = {f"module_{name}": bool(on) for name, on in mods.items()}
    for tool in registry.list_tools():
        section = cfg.get(tool)
        if not isinstance(section, dict) or not section:
            continue
        try:
            conflicts = registry.param_conflicts(tool, section, context=ctx)
        except registry.ValidationError:
            continue  # a malformed section is config-validation's job, not ours
        for c in conflicts:
            out.append(Check(f"conflict:{tool}.{c['param']}", _FAIL, c["message"],
                             remedy="resolve the conflict: unset the flag, or turn off the "
                                    "module it is incompatible with."))
    return out


def check_registry_version_drift(probe=None) -> List[Check]:
    """INFO when an installed tool's version differs from the ``tested_version`` its registry was
    curated against — a nudge to run ``metagx refresh <tool>``. Only fires for tools that both
    record a ``tested_version`` and are installed; a missing/uninstalled tool is silent (other
    checks own PATH). ``probe`` is injectable for tests (default: ``toollock.probe_tool``)."""
    from . import refresh, toollock
    probe = probe or toollock.probe_tool
    out: List[Check] = []
    for tool in registry.list_tools():
        tested = registry.version_info(tool).get("tested_version")
        if not tested:
            continue
        info = probe(tool)
        if not info.get("installed") or not info.get("version"):
            continue
        inst_tok = refresh.version_token(info["version"])
        test_tok = refresh.version_token(tested)
        if inst_tok and test_tok and inst_tok != test_tok:
            out.append(Check(
                f"version-drift:{tool}", _INFO,
                f"{tool} installed {inst_tok}, registry curated for {test_tok}.",
                remedy=f"run `metagx refresh {tool}` to review flag changes, then bump tested_version."))
    return out


def check_validate_alignment(cfg: Optional[Dict] = None) -> List[Check]:
    """STRONG guard on kraken2↔blastn alignment for the `validate` module. The validation is
    only meaningful if the BLAST reference covers the SAME organisms as the classifier — best
    guaranteed by building both together (db.build.blast) while the kraken2 library is on disk."""
    out: List[Check] = []
    if not cfg:
        return out
    mods = cfg.get("modules", {}) or {}
    if not mods.get("validate"):
        return out
    db = cfg.get("db", {}) or {}
    build = db.get("build") or {}
    val = cfg.get("validate", {}) or {}
    if build:
        if build.get("blast"):
            out.append(Check("validate:alignment", _OK,
                             "kraken2 + aligned BLAST DB are built together (db.build.blast) — "
                             "validation is in scope and taxid-tagged."))
        else:
            out.append(Check(
                "validate:alignment", _WARN,
                "STRONG WARNING: db.build will build the kraken2 DB but NOT the aligned BLAST DB "
                "(db.build.blast: false). Validation will be out-of-scope — and IMPOSSIBLE if the "
                "kraken2 library is later cleaned (kraken2-build --clean deletes the genomes blastn "
                "needs). You are validating against a different reference than you classified with.",
                remedy="Set db.build.blast: true to build both together (recommended), or keep the "
                       "kraken2 library (never run kraken2-build --clean) and set "
                       "validate.build_from: classifier. Use db.blast/validate.remote only to "
                       "compare against a deliberately broader reference."))
    elif not (db.get("blast") or val.get("build_from") or val.get("remote")):
        out.append(Check("validate:alignment", _WARN,
                         "validate is on but no BLAST reference is configured.",
                         remedy="db.blast, or validate.build_from (same genomes as the classifier), "
                                "or validate.remote."))
    elif db.get("blast") and not val.get("build_from"):
        out.append(Check("validate:alignment", _INFO,
                         "validating against db.blast — ensure it covers the SAME organisms as the "
                         "kraken2 DB, or the agreement rate reflects DB-scope differences, not "
                         "classifier error. validate.build_from: classifier guarantees alignment."))
    return out


def run(db_paths: Optional[Dict[str, str]] = None, cfg: Optional[Dict] = None) -> List[Check]:
    """Run every preflight check and return the results in display order."""
    checks: List[Check] = []
    checks += check_platform()
    checks.append(check_workflow())
    checks += check_tools()
    checks += check_active_tools(cfg)
    b = check_bracken_runs()
    if b:
        checks.append(b)
    checks.append(check_conda_frontend())
    checks.append(check_database(db_paths))
    checks += check_db_build(db_paths)
    checks += check_module_dbs(cfg)
    checks += check_config_flags(cfg)
    checks += check_param_conflicts(cfg)
    checks += check_registry_version_drift()
    checks += check_validate_alignment(cfg)
    return checks


_GLYPH = {_OK: "✓", _INFO: "·", _WARN: "!", _FAIL: "✗"}


def format_report(checks: List[Check]) -> str:
    lines = []
    for c in checks:
        lines.append(f"  {_GLYPH.get(c.status, '?')} [{c.status:>4}] {c.name}: {c.message}")
        if c.remedy and c.status in (_WARN, _FAIL):
            lines.append(f"        ↳ {c.remedy}")
    n_fail = sum(1 for c in checks if c.status == _FAIL)
    n_warn = sum(1 for c in checks if c.status == _WARN)
    lines.append("")
    if n_fail:
        lines.append(f"  {n_fail} failure(s), {n_warn} warning(s) — fix the failures before running.")
    elif n_warn:
        lines.append(f"  No failures, {n_warn} warning(s) — review the footguns above.")
    else:
        lines.append("  All checks passed.")
    return "\n".join(lines)
