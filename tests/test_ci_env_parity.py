"""CI must exercise what we ship (EVALUATION-2026-06-22, §2).

The CI e2e job used to install a hand-picked subset of tools that silently drifted from
environment.yml — so QC/CAT/kaiju/genomad/amplicon/consensus modules were advertised but
never run. These tests pin the contract: the e2e job builds its env *from* environment.yml
(parity by construction), and the version floors the doctor enforces don't reference a tool
that isn't declared somewhere we install.
"""
import os
import re

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _load(path):
    with open(os.path.join(REPO, path)) as fh:
        return yaml.safe_load(fh)


def test_ci_e2e_builds_from_environment_yml():
    ci = _load(".github/workflows/ci.yml")
    e2e = ci["jobs"]["e2e"]["steps"]
    setup = [s for s in e2e if "setup-micromamba" in (s.get("uses") or "")]
    assert setup, "e2e job no longer uses setup-micromamba"
    envfiles = [s.get("with", {}).get("environment-file") for s in setup]
    assert "environment.yml" in envfiles, (
        "CI e2e must install from environment.yml so its tool stack can't drift from what we "
        "ship/document. If you switch to create-args, you reintroduce the coverage-overselling "
        "bug — update this test deliberately if that's truly intended."
    )


def test_environment_yml_declares_core_classify_stack():
    env = _load("environment.yml")
    deps = [d for d in env["dependencies"] if isinstance(d, str)]
    names = {re.split(r"[<>=\s]", d, 1)[0] for d in deps}
    for tool in ("kraken2", "bracken", "fastp", "megahit", "minimap2", "samtools", "metabat2"):
        assert tool in names, f"{tool} missing from environment.yml core stack"


def test_doctor_floors_are_installed_somewhere():
    """Every tool the doctor enforces a floor on must be declared in environment.yml or a
    workflow/envs/*.yaml (so a user following our docs actually gets it). Guards against the
    doctor warning about a tool we never tell anyone to install."""
    from metagx.doctor import VERSION_FLOORS

    declared = set()
    env = _load("environment.yml")
    for d in env["dependencies"]:
        if isinstance(d, str):
            declared.add(re.split(r"[<>=\s]", d, 1)[0])
    envs_dir = os.path.join(REPO, "workflow", "envs")
    for fn in os.listdir(envs_dir):
        if fn.endswith(".yaml"):
            spec = _load(os.path.join("workflow", "envs", fn))
            for d in (spec or {}).get("dependencies", []) or []:
                if isinstance(d, str):
                    declared.add(re.split(r"[<>=\s:]", d, 1)[0])

    # Map the doctor's display names to their conda package names where they differ.
    alias = {"mapDamage": "mapdamage2", "metabat2": "metabat2"}
    for tool in VERSION_FLOORS:
        pkg = alias.get(tool, tool)
        assert pkg in declared, (
            f"doctor enforces a version floor for {tool!r} ({pkg}) but it is declared in neither "
            "environment.yml nor any workflow/envs/*.yaml — users can't get it from our docs."
        )
