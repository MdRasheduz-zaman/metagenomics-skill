"""Packaging guards: the Snakemake workflow + MCP server must ship with the package.

The "works only locally" trap (EVALUATION-2026-06-22, P1): a real `pip install metagx`
outside the repo used to miss workflow/ and mcp_server.py, so runner.run() died with
FileNotFoundError. These tests assert the data is declared for the wheel and that
workflow_path() resolves a *complete* workflow tree (Snakefile + every include dir),
not just the Snakefile alone.
"""
import os
import tomllib

from metagx import runner

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_workflow_path_resolves_to_a_real_file():
    p = runner.workflow_path()
    assert os.path.isfile(p)
    assert os.path.basename(p) == "Snakefile"


def test_workflow_tree_is_complete_next_to_snakefile():
    # Shipping the Snakefile but not its rules/scripts/envs would still FileNotFound at run
    # time. Assert the whole tree travels together.
    wf = os.path.dirname(runner.workflow_path())
    for sub in ("rules", "scripts", "envs", "profiles"):
        d = os.path.join(wf, sub)
        assert os.path.isdir(d), f"workflow/{sub} missing next to Snakefile"
        assert os.listdir(d), f"workflow/{sub} is empty"
    assert os.path.isfile(os.path.join(wf, "rules", "common.smk"))


def test_pyproject_force_includes_workflow_and_mcp():
    with open(os.path.join(REPO, "pyproject.toml"), "rb") as fh:
        cfg = tomllib.load(fh)
    fi = cfg["tool"]["hatch"]["build"]["targets"]["wheel"]["force-include"]
    # workflow/ and mcp_server.py must be force-included so a non-editable wheel carries them.
    assert fi.get("workflow") == "metagx/workflow"
    assert fi.get("mcp_server.py") == "metagx/mcp_server.py"


def test_workflow_path_error_is_actionable(monkeypatch):
    # When nothing resolves, the error must name the fix, not just "not found".
    monkeypatch.setattr(os.path, "isfile", lambda _p: False)
    try:
        runner.workflow_path()
        assert False, "expected FileNotFoundError"
    except FileNotFoundError as e:
        msg = str(e)
        assert "package data" in msg
        assert "pip install" in msg
