"""C1b: `runner.start_run` launches Snakemake detached and `run_status` reports progress +
final returncode — so the MCP/HTTP run tools return a job handle instead of blocking (and timing
out the transport) for the whole run.
"""
import time

from metagx import runner


def _poll(job_id, jobs_root, timeout=10.0):
    deadline = time.time() + timeout
    st = runner.run_status(job_id, jobs_root=jobs_root)
    while st["status"] not in {"done", "stopped", "unknown"} and time.time() < deadline:
        time.sleep(0.05)
        st = runner.run_status(job_id, jobs_root=jobs_root)
    return st


def test_start_run_completes_zero(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "_build_command", lambda *a, **k: ["sh", "-c", "echo hi; exit 0"])
    root = str(tmp_path / "jobs")
    job = runner.start_run(config="c.yaml", jobs_root=root)
    assert "job_id" in job and "log" in job
    st = _poll(job["job_id"], root)
    assert st["status"] == "done" and st["returncode"] == 0
    assert "hi" in st["log_tail"]


def test_start_run_records_nonzero(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "_build_command",
                        lambda *a, **k: ["sh", "-c", "echo boom 1>&2; exit 7"])
    root = str(tmp_path / "jobs")
    job = runner.start_run(config="c.yaml", jobs_root=root)
    st = _poll(job["job_id"], root)
    assert st["status"] == "done" and st["returncode"] == 7
    assert "boom" in st["log_tail"]


def test_run_status_unknown_job(tmp_path):
    st = runner.run_status("nope", jobs_root=str(tmp_path / "jobs"))
    assert st["status"] == "unknown"
