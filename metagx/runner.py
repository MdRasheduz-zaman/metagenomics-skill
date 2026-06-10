"""Thin wrapper around invoking the Snakemake workflow."""

from __future__ import annotations

import os
import subprocess
from importlib import resources
from typing import List


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
        cmd += ["--use-conda", "--conda-frontend", "conda"]
    if profile:
        # a directory containing config.yaml (e.g. the bundled SLURM profile)
        cmd += ["--workflow-profile", os.path.abspath(profile)]
    if extra:
        cmd.extend(extra)
    return subprocess.run(cmd, capture_output=True, text=True)
