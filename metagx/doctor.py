"""metagx doctor — environment preflight.

Turns the macOS/arm64 + bioconda landmines that used to live as tribal knowledge in
CLAUDE.md (Rosetta / CONDA_SUBDIR scoping, the broken Bracken osx-64 build, the abricate
samtools-0.1.x downgrade, missing tools, a missing database) into machine-checked
diagnostics. Each check yields a ``Check`` with a status and, when something is wrong, the
*exact* remedy — so a stranger on their own machine gets steered, not stranded.

Statuses:
  ok    — verified good.
  info  — context, no action needed.
  warn  — works, but a known footgun is armed (e.g. CONDA_SUBDIR leaking to base).
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

from . import report, runner

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
    if sys_name == "Darwin" and machine == "arm64":
        out.append(Check(
            "apple-silicon", _INFO,
            "Apple Silicon (arm64): bioconda lacks native builds for several tools, so the "
            "supported install runs them x86_64 under Rosetta (CONDA_SUBDIR=osx-64).",
            remedy="The most reliable path on this machine is Docker or a Linux box. If you "
                   "install natively, follow scripts/install_bio_macos_arm64.sh and scope "
                   "CONDA_SUBDIR to that single command — never export it into your base shell.",
        ))
    return out


def check_conda_subdir_leak() -> Check:
    """A globally-exported CONDA_SUBDIR=osx-64 will pull x86_64 packages into the arm64 base
    env and corrupt it — the single nastiest documented footgun."""
    val = os.environ.get("CONDA_SUBDIR")
    if not val:
        return Check("conda-subdir", _OK, "CONDA_SUBDIR is not exported globally.")
    if val == "osx-64" and platform.machine() == "arm64":
        return Check(
            "conda-subdir", _WARN,
            f"CONDA_SUBDIR={val} is exported in this shell.",
            remedy="On Apple Silicon this can drag x86_64 packages into your arm64 base env and "
                   "break conda. Unset it (`unset CONDA_SUBDIR`) and scope it inline only to the "
                   "specific bioconda install command that needs it.",
        )
    return Check("conda-subdir", _INFO, f"CONDA_SUBDIR={val} is exported.")


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
                    remedy="Install the core stack: `conda env create -f environment.yml` "
                           "(or `mamba env create -f environment.yml`), then activate it."))
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


def check_bracken_runs() -> Optional[Check]:
    """Bracken's osx-64 conda build is known-broken; a present-but-unrunnable binary is worse
    than a missing one because it passes a naive `which` check. Actually invoke it."""
    exe = shutil.which("bracken")
    if not exe:
        return None  # absence is already reported by check_tools
    try:
        p = subprocess.run(["bracken", "-h"], capture_output=True, text=True, timeout=20)
    except (subprocess.SubprocessError, OSError) as e:
        return Check("bracken-runs", _FAIL, f"bracken is on PATH but failed to execute: {e}",
                     remedy="On Apple Silicon the osx-64 Bracken conda build is broken — "
                            "see scripts/install_bio_macos_arm64.sh, or use Docker/Linux.")
    blob = (p.stdout or "") + (p.stderr or "")
    if p.returncode != 0 and "usage" not in blob.lower() and "bracken" not in blob.lower():
        return Check("bracken-runs", _FAIL,
                     f"bracken is on PATH but `bracken -h` exited {p.returncode} with no usage text.",
                     remedy="Likely the broken osx-64 build. Reinstall per "
                            "scripts/install_bio_macos_arm64.sh, or use Docker/Linux.")
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
    if not kdb:
        return Check("database", _WARN, "Config has no db.kraken2 entry.",
                     remedy="Set db.kraken2 to a kraken2 index directory, or run `metagx fetch-db`.")
    info = report.db_info(kdb)
    if not info.get("present"):
        return Check("database", _FAIL, f"Configured kraken2 db not found: {kdb}",
                     remedy="Download one with `metagx fetch-db standard-8 --dir <path>` "
                            "or build a custom one with `metagx build-db`.")
    if not os.path.isfile(os.path.join(kdb, "hash.k2d")):
        return Check("database", _FAIL, f"{kdb} exists but has no hash.k2d (not a built kraken2 db).",
                     remedy="Finish the build (`metagx build-db`) or re-download (`metagx fetch-db`).")
    size_gb = info.get("size_bytes", 0) / 1e9
    return Check("database", _OK, f"kraken2 db OK: {kdb} ({size_gb:.1f} GB).")


def run(db_paths: Optional[Dict[str, str]] = None) -> List[Check]:
    """Run every preflight check and return the results in display order."""
    checks: List[Check] = []
    checks += check_platform()
    checks.append(check_conda_subdir_leak())
    checks.append(check_workflow())
    checks += check_tools()
    b = check_bracken_runs()
    if b:
        checks.append(b)
    checks.append(check_conda_frontend())
    checks.append(check_database(db_paths))
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
