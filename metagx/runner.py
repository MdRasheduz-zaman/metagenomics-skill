"""Thin wrapper around invoking the Snakemake workflow."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from importlib import resources
from typing import List, Optional, Tuple

# Snakemake 8+ refuses to drive --use-conda with an older conda; mamba has no such gate.
MIN_CONDA = (24, 7, 1)


class CondaFrontendError(RuntimeError):
    """The conda/mamba frontend can't drive ``--use-conda``; raised before launching Snakemake."""


def _parse_version(text: Optional[str]) -> Optional[Tuple[int, int, int]]:
    """Extract a leading x.y.z from e.g. 'conda 23.10.0' / 'mamba 1.5.8'."""
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text or "")
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2)), int(m.group(3) or 0))


def _frontend_version(frontend: str) -> Optional[str]:
    try:
        out = subprocess.run([frontend, "--version"], capture_output=True, text=True)
    except (FileNotFoundError, OSError):
        return None
    return (out.stdout or out.stderr).strip() or None


def pick_conda_frontend() -> str:
    """Prefer mamba (faster, and sidesteps Snakemake's conda-version gate), else conda."""
    return "mamba" if shutil.which("mamba") else "conda"


def conda_preflight(frontend: str, version_str: Optional[str] = None) -> Optional[str]:
    """Return an actionable error string if the frontend can't run ``--use-conda``, else None."""
    if version_str is None and shutil.which(frontend) is None:
        return (f"--use-conda needs '{frontend}' on PATH, but it was not found. "
                f"Install mamba (`conda install -n base -c conda-forge mamba`) or drop --use-conda.")
    if version_str is None:
        version_str = _frontend_version(frontend)
    if frontend == "mamba":
        return None  # mamba satisfies Snakemake's solver requirement regardless of version
    ver = _parse_version(version_str)
    if ver is None:
        return None  # unrecognized format — let Snakemake make the call
    if ver < MIN_CONDA:
        need = ".".join(map(str, MIN_CONDA))
        got = ".".join(map(str, ver))
        return (f"--use-conda needs conda >= {need} (a Snakemake 8+ requirement), but found {got}. "
                f"Update it (`conda update -n base conda`) or install mamba "
                f"(`conda install -n base -c conda-forge mamba`), then retry.")
    return None


def workflow_path() -> str:
    """Absolute path to the bundled Snakefile.

    Resolves in priority order so the tool works from a git clone (editable
    install), from a real wheel install (``workflow/`` shipped as package data
    under ``metagx/workflow/``), or when invoked inside the repo:

      1. sibling of the package — ``<repo>/workflow/Snakefile`` (clone / editable)
      2. inside the installed package — ``<site-packages>/metagx/workflow/Snakefile``
      3. the current working directory — ``./workflow/Snakefile``
    """
    pkg_dir = os.path.dirname(os.path.abspath(__file__))           # .../metagx
    for candidate in (
        os.path.join(os.path.dirname(pkg_dir), "workflow", "Snakefile"),  # (1) repo sibling
        os.path.join(pkg_dir, "workflow", "Snakefile"),                   # (2) packaged data
        os.path.join(os.getcwd(), "workflow", "Snakefile"),              # (3) cwd
    ):
        if os.path.isfile(candidate):
            return candidate
    raise FileNotFoundError(
        "Could not locate workflow/Snakefile. metagx ships the Snakemake workflow "
        "as package data, so a normal install should find it. If you installed from "
        "a source checkout with `pip install -e .`, run metagx from inside the repo "
        "(the workflow/ directory must sit next to the metagx/ package), or reinstall "
        "with a build that includes package data (`pip install .`). "
        "See the README 'Installation' section."
    )


def environment_file_path() -> Optional[str]:
    """Absolute path to the bundled ``environment.yml`` (the core conda env spec), or None.

    Same resolution as ``workflow_path``: repo-sibling (editable/clone) then packaged
    (``metagx/environment.yml`` in a wheel). A wheel-only end user has no ``environment.yml``
    in their cwd, so doctor's "conda env create -f environment.yml" remedy would otherwise point
    at a file they don't have — ``metagx env-file`` surfaces/copies the packaged one.
    """
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    for candidate in (
        os.path.join(os.path.dirname(pkg_dir), "environment.yml"),  # repo sibling (editable)
        os.path.join(pkg_dir, "environment.yml"),                   # packaged in the wheel
    ):
        if os.path.isfile(candidate):
            return candidate
    return None


def _build_command(config: str, cores: str | int, dry_run: bool, use_conda: bool,
                   profile: str | None, extra: List[str] | None) -> List[str]:
    """The Snakemake argv, shared by ``run`` (foreground) and ``start_run`` (background).

    Invoked under *this* interpreter (``sys.executable -m snakemake``), not a bare ``snakemake``
    off PATH: the workflow's common.smk does ``from metagx import ...`` in the Snakemake process,
    so PATH could pick an env without metagx → ModuleNotFoundError at load. Raises
    ``CondaFrontendError`` when ``use_conda`` is requested but the frontend can't drive it.
    """
    cmd = [
        sys.executable, "-m", "snakemake",
        "--snakefile", workflow_path(),
        "--configfile", os.path.abspath(config),
        "--cores", str(cores),
        "--printshellcmds",
    ]
    if dry_run:
        cmd.append("--dry-run")
    if use_conda:
        frontend = pick_conda_frontend()
        if not dry_run:                       # a plain dry-run never provisions envs
            problem = conda_preflight(frontend)
            if problem:
                raise CondaFrontendError(problem)
        cmd += ["--use-conda", "--conda-frontend", frontend]
    if profile:
        cmd += ["--workflow-profile", os.path.abspath(profile)]
    if extra:
        cmd.extend(extra)
    return cmd


def _jobs_root(root: str | None) -> str:
    return root or os.path.join(".metagx", "jobs")


def _pid_alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except (OSError, ProcessLookupError):
        return False
    return True


def start_run(
    config: str = "config.yaml",
    cores: str | int = "all",
    use_conda: bool = False,
    profile: str | None = None,
    extra: List[str] | None = None,
    jobs_root: str | None = None,
) -> dict:
    """Launch the workflow **detached** and return a job handle immediately.

    For the MCP/HTTP surface: a real run is minutes-to-hours, so blocking the request until it
    finishes times out the transport (and only the tail of the log survives). This spawns Snakemake
    in its own session, tees all output to ``<jobs_root>/<job_id>/run.log``, and records the exit
    code to ``returncode`` on completion. Poll with ``run_status(job_id)``. Raises
    ``CondaFrontendError`` before launching if ``use_conda`` can't be driven.
    """
    import json
    import shlex
    import time

    cmd = _build_command(config, cores, False, use_conda, profile, extra)
    root = _jobs_root(jobs_root)
    job_id = time.strftime("%Y%m%dT%H%M%S") + f"-{os.getpid()}"
    d = os.path.join(root, job_id)
    os.makedirs(d, exist_ok=True)
    log = os.path.join(d, "run.log")
    rc = os.path.join(d, "returncode")
    inner = " ".join(shlex.quote(a) for a in cmd)
    # tee to the log and record the exit code so run_status can report done+returncode after the
    # detached process is gone (we can't wait() on a reparented child).
    wrapper = f"{inner} > {shlex.quote(log)} 2>&1; printf '%s' \"$?\" > {shlex.quote(rc)}"
    proc = subprocess.Popen(["sh", "-c", wrapper], start_new_session=True)
    meta = {"job_id": job_id, "pid": proc.pid, "config": os.path.abspath(config),
            "dir": d, "log": log, "started": time.time()}
    with open(os.path.join(d, "job.json"), "w") as fh:
        json.dump(meta, fh)
    return meta


def run_status(job_id: str, jobs_root: str | None = None, tail_lines: int = 40) -> dict:
    """Status of a ``start_run`` job: ``running`` | ``done`` (with ``returncode``) | ``stopped`` |
    ``unknown``, plus the last ``tail_lines`` of the log."""
    import json

    d = os.path.join(_jobs_root(jobs_root), job_id)
    meta_path = os.path.join(d, "job.json")
    if not os.path.isfile(meta_path):
        return {"job_id": job_id, "status": "unknown", "error": "no such job"}
    with open(meta_path) as fh:
        meta = json.load(fh)
    log = meta.get("log", os.path.join(d, "run.log"))
    tail = ""
    if os.path.isfile(log):
        with open(log) as fh:
            tail = "".join(fh.readlines()[-tail_lines:])
    rc_path = os.path.join(d, "returncode")
    if os.path.isfile(rc_path):
        with open(rc_path) as fh:
            raw = fh.read().strip()
        rc = int(raw) if raw.isdigit() else 1
        return {"job_id": job_id, "status": "done", "returncode": rc, "log": log, "log_tail": tail}
    status = "running" if _pid_alive(meta.get("pid", -1)) else "stopped"
    return {"job_id": job_id, "status": status, "log": log, "log_tail": tail}


def run(
    config: str = "config.yaml",
    cores: str | int = "all",
    dry_run: bool = False,
    use_conda: bool = False,
    profile: str | None = None,
    extra: List[str] | None = None,
    stream: bool = False,
) -> subprocess.CompletedProcess:
    """Run the workflow. Returns the completed process (check returncode/stdout).

    ``use_conda`` makes Snakemake create/use the per-rule conda envs under workflow/envs/,
    so the heavy domain tools are provisioned automatically on first run.

    ``stream``: inherit this process's stdout/stderr so Snakemake's progress is shown **live**
    instead of being captured and dumped only at the end. A metagenomics run is minutes-to-hours;
    a captured run looks hung and loses everything on Ctrl-C. The CLI streams; callers that need the
    log text back (the MCP/HTTP surface) leave ``stream`` off and read ``stdout``/``stderr`` (which
    are ``None`` when streaming). The advisor/history read result *files*, not this output, so
    streaming costs nothing there.
    """
    cmd = _build_command(config, cores, dry_run, use_conda, profile, extra)
    if stream:
        # Inherit stdio → live progress. stdout/stderr on the result are None by design.
        return subprocess.run(cmd, text=True)
    return subprocess.run(cmd, capture_output=True, text=True)
