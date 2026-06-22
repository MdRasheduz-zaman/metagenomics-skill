"""Conda-frontend preflight for `metagx run --use-conda` (`metagx/runner.py`).

Snakemake 8+ refuses an older conda; the runner now prefers mamba and, on conda, fails fast
with an actionable message instead of a cryptic Snakemake traceback.
"""
import sys

from metagx import runner


def test_parse_version_extracts_xyz():
    assert runner._parse_version("conda 23.10.0") == (23, 10, 0)
    assert runner._parse_version("mamba 1.5.8") == (1, 5, 8)
    assert runner._parse_version("24.7") == (24, 7, 0)
    assert runner._parse_version("no version here") is None


def test_conda_too_old_returns_actionable_message():
    msg = runner.conda_preflight("conda", version_str="conda 23.10.0")
    assert msg is not None
    assert "24.7.1" in msg and "23.10.0" in msg
    assert "mamba" in msg  # offers the faster way out


def test_conda_new_enough_passes():
    assert runner.conda_preflight("conda", version_str="conda 24.7.1") is None
    assert runner.conda_preflight("conda", version_str="conda 26.5.2") is None


def test_mamba_always_ok_regardless_of_version():
    # mamba has no version gate; any parseable version (or none) passes
    assert runner.conda_preflight("mamba", version_str="mamba 1.5.8") is None


def test_unrecognized_version_defers_to_snakemake():
    assert runner.conda_preflight("conda", version_str="conda weird-build") is None


def test_run_raises_before_launching_when_conda_too_old(monkeypatch, tmp_path):
    # use_conda + old conda + no mamba -> CondaFrontendError, and snakemake is never spawned.
    monkeypatch.setattr(runner, "pick_conda_frontend", lambda: "conda")
    monkeypatch.setattr(runner, "conda_preflight", lambda f, version_str=None: "conda too old")
    monkeypatch.setattr(runner.subprocess, "run",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("snakemake spawned!")))
    cfg = tmp_path / "config.yaml"
    cfg.write_text("project: x\n")
    try:
        runner.run(config=str(cfg), use_conda=True)
        assert False, "expected CondaFrontendError"
    except runner.CondaFrontendError as e:
        assert "too old" in str(e)


def test_dry_run_skips_preflight(monkeypatch, tmp_path):
    # a dry-run never provisions envs, so it must not be blocked by the frontend check
    monkeypatch.setattr(runner, "pick_conda_frontend", lambda: "conda")
    monkeypatch.setattr(runner, "conda_preflight",
                        lambda f, version_str=None: (_ for _ in ()).throw(
                            AssertionError("preflight ran on dry-run")))
    captured = {}

    def fake_run(cmd, **k):
        captured["cmd"] = cmd
        class P:
            returncode, stdout, stderr = 0, "", ""
        return P()

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    cfg = tmp_path / "config.yaml"
    cfg.write_text("project: x\n")
    runner.run(config=str(cfg), use_conda=True, dry_run=True)
    assert "--use-conda" in captured["cmd"] and "--dry-run" in captured["cmd"]


def test_snakemake_runs_under_this_interpreter(monkeypatch, tmp_path):
    """The Snakemake subprocess must run as `sys.executable -m snakemake`, not a bare
    `snakemake` off PATH. The workflow's common.smk does `from metagx import ...` in the
    Snakemake process's Python; resolving snakemake via PATH can pick an interpreter that
    lacks metagx (e.g. running `metagx` by absolute path without activating its env), which
    failed with `ModuleNotFoundError: No module named 'metagx'` at DAG-load time. Pinning the
    invocation to this interpreter guarantees the subprocess shares metagx's environment.
    """
    captured = {}

    def fake_run(cmd, **k):
        captured["cmd"] = cmd
        class P:
            returncode, stdout, stderr = 0, "", ""
        return P()

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    cfg = tmp_path / "config.yaml"
    cfg.write_text("project: x\n")
    runner.run(config=str(cfg), dry_run=True)

    cmd = captured["cmd"]
    assert cmd[:3] == [sys.executable, "-m", "snakemake"], cmd
    # never invoke a bare PATH-resolved snakemake
    assert cmd[0] != "snakemake"
