"""Thin wrapper around invoking the Snakemake workflow."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
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

    Resolves whether running from the source tree or an installed package.
    """
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidate = os.path.join(here, "workflow", "Snakefile")
    if os.path.isfile(candidate):
        return candidate
    # fallback: cwd
    candidate = os.path.join(os.getcwd(), "workflow", "Snakefile")
    if os.path.isfile(candidate):
        return candidate
    raise FileNotFoundError("Could not locate workflow/Snakefile")


def run(
    config: str = "config.yaml",
    cores: str | int = "all",
    dry_run: bool = False,
    use_conda: bool = False,
    profile: str | None = None,
    extra: List[str] | None = None,
) -> subprocess.CompletedProcess:
    """Run the workflow. Returns the completed process (check returncode/stdout).

    ``use_conda`` makes Snakemake create/use the per-rule conda envs under workflow/envs/,
    so the heavy domain tools are provisioned automatically on first run.
    """
    cmd = [
        "snakemake",
        "--snakefile", workflow_path(),
        "--configfile", os.path.abspath(config),
        "--cores", str(cores),
        "--printshellcmds",
    ]
    if dry_run:
        cmd.append("--dry-run")
    if use_conda:
        frontend = pick_conda_frontend()
        # A plain dry-run never provisions envs, so don't block it on the frontend version.
        if not dry_run:
            problem = conda_preflight(frontend)
            if problem:
                raise CondaFrontendError(problem)
        cmd += ["--use-conda", "--conda-frontend", frontend]
    if profile:
        # a directory containing config.yaml (e.g. the bundled SLURM profile)
        cmd += ["--workflow-profile", os.path.abspath(profile)]
    if extra:
        cmd.extend(extra)
    return subprocess.run(cmd, capture_output=True, text=True)
