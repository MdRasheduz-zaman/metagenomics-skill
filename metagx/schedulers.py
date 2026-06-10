"""Cluster/scheduler execution backends — the single source of truth for how
metagx submits the Snakemake workflow to an HPC scheduler.

Each entry maps a short name (what the user passes to ``metagx run --executor``) to:
  * the Snakemake v8 executor it uses,
  * the bundled profile directory under ``workflow/profiles/<name>/``,
  * the plugin package that must be installed for it to work, and
  * a one-line note on what to edit before first use.

The CLI resolver, ``metagx schedulers`` discovery, and the docs all read this dict,
so adding a backend is a one-place change (mirrors the parameter-registry philosophy).
"""

from __future__ import annotations

import os
from typing import Dict, List

# name -> metadata. `profile` is the directory name under workflow/profiles/.
SCHEDULERS: Dict[str, Dict[str, str]] = {
    "local": {
        "profile": "local",
        "executor": "(none — local cores)",
        "plugin": "(built in)",
        "summary": "A single fat node / workstation with no scheduler. Per-rule "
                   "thread + memory caps, no job submission.",
        "edit": "Optionally lower mem/threads to fit your machine.",
    },
    "slurm": {
        "profile": "slurm",
        "executor": "slurm",
        "plugin": "snakemake-executor-plugin-slurm",
        "summary": "SLURM (sbatch). Native Snakemake v8 plugin — preferred on SLURM clusters.",
        "edit": "Set slurm_partition and slurm_account.",
    },
    "lsf": {
        "profile": "lsf",
        "executor": "lsf",
        "plugin": "snakemake-executor-plugin-lsf",
        "summary": "IBM LSF / OpenLAVA (bsub). Native Snakemake v8 plugin.",
        "edit": "Set lsf_queue (and lsf_project if your site requires it).",
    },
    "sge": {
        "profile": "sge",
        "executor": "cluster-generic",
        "plugin": "snakemake-executor-plugin-cluster-generic",
        "summary": "Sun Grid Engine / SGE / UGE / OGS (qsub) via the generic-cluster executor.",
        "edit": "Set the parallel-environment name (-pe) and memory resource "
                "(h_vmem vs mem_free) for your site in cluster-generic-submit-cmd.",
    },
    "pbs": {
        "profile": "pbs",
        "executor": "cluster-generic",
        "plugin": "snakemake-executor-plugin-cluster-generic",
        "summary": "PBS Pro / TORQUE / OpenPBS (qsub) via the generic-cluster executor.",
        "edit": "Set the queue (-q) and confirm the -l resource syntax (nodes=1:ppn vs "
                "select=1:ncpus) for your PBS flavour.",
    },
    "generic": {
        "profile": "generic",
        "executor": "cluster-generic",
        "plugin": "snakemake-executor-plugin-cluster-generic",
        "summary": "Any other scheduler (HTCondor, Moab, OAR, Flux, …) — fill in your "
                   "own submit command.",
        "edit": "Replace cluster-generic-submit-cmd with your scheduler's submit command.",
    },
}


def list_schedulers() -> List[str]:
    return list(SCHEDULERS)


def profiles_dir() -> str:
    """Absolute path to the bundled ``workflow/profiles`` directory."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(here, "workflow", "profiles")


def profile_path(name: str) -> str:
    """Resolve a scheduler name to its bundled profile directory.

    Raises ``KeyError`` (with the valid names) for an unknown scheduler, and
    ``FileNotFoundError`` if the bundled profile is missing.
    """
    key = (name or "").strip().lower()
    if key not in SCHEDULERS:
        raise KeyError(
            f"Unknown scheduler '{name}'. Choose from: {', '.join(SCHEDULERS)}."
        )
    path = os.path.join(profiles_dir(), SCHEDULERS[key]["profile"])
    if not os.path.isdir(path):
        raise FileNotFoundError(f"Bundled profile for '{key}' not found at {path}")
    return path


def describe() -> List[Dict[str, str]]:
    """Rows for ``metagx schedulers`` (name + metadata, in declared order)."""
    rows = []
    for name, meta in SCHEDULERS.items():
        rows.append({"name": name, **meta})
    return rows
