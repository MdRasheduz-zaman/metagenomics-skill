"""Custom database builders (kraken2 / CAT / Kaiju) — no NCBI download.

The plan-only tests need no external tools. The real-build test runs prodigal +
kaiju-mkbwt/mkfmi when present (skips in CI), verifying the Kaiju index that the
consensus module's db.kaiju consumes.
"""
import os
import shutil
import subprocess

import pytest

from metagx import dbbuild

_GENOMES = os.path.join(os.path.dirname(__file__), "..", "data", "genomes.fasta")


def _write_genomes(path):
    # two tiny "genomes" — enough to exercise the taxid mapping / planning logic
    path.write_text(">genomeA description A\nACGTACGTACGTACGT\n"
                    ">genomeB description B\nTTTTGGGGCCCCAAAA\n")
    return str(path)


def test_build_kaiju_db_plan_only(tmp_path):
    g = _write_genomes(tmp_path / "g.fasta")
    res = dbbuild.build_kaiju_db(g, str(tmp_path / "kdb"),
                                 taxonomy_dir=str(tmp_path / "tax"), run=False)
    assert res["n_genomes"] == 2
    assert res["ran"] is False
    # the planned commands name the right tools + the consensus-expected .fmi output
    assert "prodigal" in res["commands"]["prodigal"]
    assert "kaiju-mkbwt" in res["commands"]["mkbwt"]
    assert res["fmi"].endswith("kaiju_db.fmi")


def test_build_kaiju_db_reports_missing_tools(tmp_path, monkeypatch):
    g = _write_genomes(tmp_path / "g.fasta")
    monkeypatch.setattr(dbbuild, "_have", lambda t: False)
    res = dbbuild.build_kaiju_db(g, str(tmp_path / "kdb"),
                                 taxonomy_dir=str(tmp_path / "tax"), run=True)
    assert res["ran"] is False and "not on PATH" in res["note"]


def test_build_db_recovers_from_kraken2_build_sigpipe(tmp_path, monkeypatch):
    """kraken2-build exits 64 via a SIGPIPE in its internal `cat | build_db` pipe on small
    DBs, yet writes a valid database. build_db must trust the artifacts over the exit code.

    This is the exact failure that broke the CI e2e job (the build step returned 64 with
    "xargs: cat: terminated by signal 13") while the same DB built fine on macOS where the
    prebuilt DB meant the build never ran. Regression guard for that false negative.
    """
    g = _write_genomes(tmp_path / "g.fasta")
    db_dir = tmp_path / "db"
    monkeypatch.setattr(dbbuild, "_have", lambda t: True)

    def fake_run(cmd, *a, **k):
        db = cmd[cmd.index("--db") + 1] if "--db" in cmd else cmd[cmd.index("-d") + 1]
        if "--build" in cmd:
            # emulate the real DB artifacts being written, then a non-zero SIGPIPE exit
            for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
                (tmp_path / "db" / f).write_text("x")
            return subprocess.CompletedProcess(cmd, 64, "Building database files (step 3)...",
                                               "xargs: cat: terminated by signal 13")
        if "bracken-build" in cmd[0]:
            length = cmd[cmd.index("-l") + 1]
            (tmp_path / "db" / f"database{length}mers.kmer_distrib").write_text("x")
            return subprocess.CompletedProcess(cmd, 0, "ok", "")
        return subprocess.CompletedProcess(cmd, 0, "ok", "")  # add-to-library

    monkeypatch.setattr(dbbuild.subprocess, "run", fake_run)
    res = dbbuild.build_db(g, str(db_dir), read_length=[150, 1000], run=True)
    assert res["ok"] is True, res
    assert "build" in res.get("recovered", [])
    assert "SIGPIPE" in res.get("note", "")


def test_build_db_recovers_via_single_threaded_retry(tmp_path, monkeypatch):
    """On a low-core runner, multithreaded `kraken2-build --build` aborts in step 3 *before*
    writing any `*.k2d` (the `cat | build_db` pipe races, `cat` dies with SIGPIPE). The
    artifact-trust path can't help — there are no artifacts. build_db must retry once with
    `--threads 1`, which removes the race.

    This is the second-round CI e2e failure (identical "xargs: cat: terminated by signal 13"
    but with NO database written), which the artifact-only recovery did not catch.
    """
    g = _write_genomes(tmp_path / "g.fasta")
    db_dir = tmp_path / "db"
    monkeypatch.setattr(dbbuild, "_have", lambda t: True)

    def fake_run(cmd, *a, **k):
        if "--build" in cmd:
            threads = cmd[cmd.index("--threads") + 1]
            if threads != "1":
                # multithreaded build aborts mid-step-3 writing nothing
                return subprocess.CompletedProcess(
                    cmd, 64, "Building database files (step 3)...",
                    "build_db: OMP only wants you to use 2 threads\nxargs: cat: terminated by signal 13")
            for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):  # single-threaded retry succeeds
                (db_dir / f).write_text("x")
            return subprocess.CompletedProcess(cmd, 0, "Database construction complete.", "")
        if "bracken-build" in cmd[0]:
            length = cmd[cmd.index("-l") + 1]
            (db_dir / f"database{length}mers.kmer_distrib").write_text("x")
            return subprocess.CompletedProcess(cmd, 0, "ok", "")
        return subprocess.CompletedProcess(cmd, 0, "ok", "")  # add-to-library

    monkeypatch.setattr(dbbuild.subprocess, "run", fake_run)
    res = dbbuild.build_db(g, str(db_dir), read_length=[150], threads=4, run=True)
    assert res["ok"] is True, res
    assert "build" in res.get("recovered", [])
    assert res["logs"]["build"].get("retry_threads1", {}).get("returncode") == 0


def test_build_db_clamps_threads_to_cpu_count(tmp_path, monkeypatch):
    """threads must never exceed online CPUs — bracken-build's kmer2read_distr aborts (rc 1,
    "thread count exceeds number of processors") instead of reducing. This was the real CI
    e2e failure on the 2-core runner once the kraken2 build step started passing.
    """
    g = _write_genomes(tmp_path / "g.fasta")
    monkeypatch.setattr(dbbuild, "_have", lambda t: True)
    monkeypatch.setattr(dbbuild, "_usable_cpus", lambda: 2)  # emulate a 2-core runner
    seen = {}

    def fake_run(cmd, *a, **k):
        if "bracken-build" in cmd[0]:
            seen["bracken_t"] = cmd[cmd.index("-t") + 1]
            length = cmd[cmd.index("-l") + 1]
            (tmp_path / "db" / f"database{length}mers.kmer_distrib").write_text("x")
        elif "--build" in cmd:
            seen["kraken_threads"] = cmd[cmd.index("--threads") + 1]
            for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
                (tmp_path / "db" / f).write_text("x")
        return subprocess.CompletedProcess(cmd, 0, "ok", "")

    monkeypatch.setattr(dbbuild.subprocess, "run", fake_run)
    res = dbbuild.build_db(g, str(tmp_path / "db"), read_length=[150], threads=4, run=True)
    assert res["ok"] is True
    assert res["threads"] == 2
    assert seen["bracken_t"] == "2" and seen["kraken_threads"] == "2"
    assert "clamped" in res.get("note_threads", "")


def test_build_db_real_failure_still_fails(tmp_path, monkeypatch):
    """A non-zero exit with *no* artifacts is a genuine failure — must not be swallowed."""
    g = _write_genomes(tmp_path / "g.fasta")
    monkeypatch.setattr(dbbuild, "_have", lambda t: True)

    def fake_run(cmd, *a, **k):
        if "--build" in cmd:  # exits non-zero and writes nothing
            return subprocess.CompletedProcess(cmd, 1, "", "boom: out of memory")
        return subprocess.CompletedProcess(cmd, 0, "ok", "")

    monkeypatch.setattr(dbbuild.subprocess, "run", fake_run)
    res = dbbuild.build_db(g, str(tmp_path / "db"), run=True)
    assert res["ok"] is False and res["failed_step"] == "build"


@pytest.mark.skipif(
    not (shutil.which("prodigal") and shutil.which("kaiju-mkbwt")
         and shutil.which("kaiju-mkfmi") and os.path.isfile(_GENOMES)),
    reason="prodigal/kaiju build tools or bundled genomes absent (skips in CI)")
def test_build_kaiju_db_real(tmp_path):
    """Build a real Kaiju index from the bundled genomes; assert the db.kaiju layout."""
    db = tmp_path / "kdb"
    base = tmp_path / "ktax"
    dbbuild.write_library_and_taxonomy(_GENOMES, str(base))  # writes base/taxonomy/{nodes,names}.dmp
    tax = base / "taxonomy"
    res = dbbuild.build_kaiju_db(_GENOMES, str(db), taxonomy_dir=str(tax), threads=2)
    assert res["ran"] and res["ok"], res.get("tail")
    assert res["n_proteins"] > 0
    # the directory is a drop-in db.kaiju for rules/consensus.smk
    for f in ("kaiju_db.fmi", "nodes.dmp", "names.dmp"):
        assert (db / f).is_file(), f"missing {f}"
