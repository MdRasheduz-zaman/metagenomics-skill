"""PR 1 regression tests — the CLI surface + hostile-input seams the suite under-guarded.

These close the exact test gaps that let the findings ship (see the codebase assessment, Part 5):
the CLI entry layer was tested only by proxy (via the underlying functions), and no test fed a
BOM/Excel sheet. Each test maps to one fix:

  C1  runner streams instead of capturing (live progress on a long run)
  C2  `doctor --config <missing>` fails loud instead of silently running config-less
  H1  a bad path prints a one-line error, not a raw traceback
  H2  sample-sheet readers tolerate a UTF-8 BOM (Excel) + header whitespace
"""
import subprocess
import types

import pytest

from metagx import cli, config_builder, formats, probe, runner


# --------------------------------------------------------------------------- #
# H2 — BOM / whitespace tolerance on every sample-sheet reader                 #
# --------------------------------------------------------------------------- #
_BOM = "﻿"


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return str(path)


def test_read_tsv_dicts_strips_utf8_bom(tmp_path):
    # Excel "UTF-8" export prepends a BOM; without utf-8-sig the first header becomes "﻿sample"
    # and row["sample"] raises KeyError.
    p = _write(tmp_path / "s.tsv", f"{_BOM}sample\tr1\tplatform\ns1\treads.fq.gz\tont\n")
    rows = formats.read_tsv_dicts(p)
    assert rows[0]["sample"] == "s1"           # not "﻿sample"
    assert "sample" in rows[0] and f"{_BOM}sample" not in rows[0]


def test_read_tsv_dicts_normalizes_header_whitespace(tmp_path):
    p = _write(tmp_path / "s.tsv", " sample \t r1 \nX\treads.fq.gz\n")
    rows = formats.read_tsv_dicts(p)
    assert rows[0]["sample"] == "X"


def test_probe_load_sheet_handles_bom(tmp_path):
    p = _write(tmp_path / "s.tsv", f"{_BOM}sample\tr1\tplatform\ns1\treads.fq.gz\tont\n")
    rows = probe.load_sheet(p)
    assert rows[0]["sample"] == "s1"           # would KeyError-free; values stripped


def test_config_builder_readers_handle_bom(tmp_path):
    # _sheet_platforms and _any_provided_contigs both went through the shared reader.
    p = _write(tmp_path / "s.tsv",
               f"{_BOM}sample\tr1\tplatform\tcontigs\ns1\t\tillumina\tgenome.fa\n")
    assert config_builder._sheet_platforms(p) == {"illumina"}
    assert config_builder._any_provided_contigs(p) is True


# --------------------------------------------------------------------------- #
# C1 — `metagx run` streams (does not capture) so progress is live            #
# --------------------------------------------------------------------------- #
def _fake_completed(**kw):
    return subprocess.CompletedProcess(args=["x"], returncode=0, stdout=None, stderr=None)


def test_runner_stream_true_does_not_capture(monkeypatch):
    seen = {}

    def fake_run(cmd, **kwargs):
        seen.update(kwargs)
        return _fake_completed()

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    runner.run(config="c.yaml", dry_run=True, stream=True)
    # streaming must inherit stdio, i.e. NOT pass capture_output=True
    assert seen.get("capture_output") is not True


def test_runner_stream_false_still_captures(monkeypatch):
    seen = {}

    def fake_run(cmd, **kwargs):
        seen.update(kwargs)
        return _fake_completed()

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    runner.run(config="c.yaml", dry_run=True, stream=False)
    assert seen.get("capture_output") is True   # the MCP/HTTP surface still needs the text back


def test_cmd_run_streams_and_tolerates_none_output(monkeypatch):
    passed = {}

    def fake_run(**kwargs):
        passed.update(kwargs)
        return _fake_completed()                 # stdout/stderr None, as when streaming

    monkeypatch.setattr(cli.runner, "run", fake_run)
    args = types.SimpleNamespace(config="c.yaml", cores="all", dry_run=True, use_conda=False,
                                 profile=None, executor=None, slurm=False,
                                 no_history=True, no_advisor=True, history_file=None)
    rc = cli.cmd_run(args)                        # must not crash writing a None stdout
    assert rc == 0
    assert passed.get("stream") is True


# --------------------------------------------------------------------------- #
# C2 — doctor --config <missing> FAILS instead of silently running config-less #
# --------------------------------------------------------------------------- #
def test_doctor_missing_config_fails_loud(monkeypatch, capsys):
    monkeypatch.setattr(cli.doctor, "run", lambda **kw: [])   # skip real tool probing (fast)
    rc = cli.main(["doctor", "--config", "/no/such/config.yaml"])
    out = capsys.readouterr().out
    assert rc == 1
    assert "not found" in out and "/no/such/config.yaml" in out


def test_doctor_no_config_is_not_a_failure(monkeypatch):
    monkeypatch.setattr(cli.doctor, "run", lambda **kw: [])
    assert cli.main(["doctor"]) == 0             # omitting --config is fine (env-only check)


# --------------------------------------------------------------------------- #
# H1 — a bad path is a one-line error, not a traceback                         #
# --------------------------------------------------------------------------- #
def test_main_missing_file_is_clean_error(capsys):
    rc = cli.main(["validate", "/no/such/config.yaml"])
    err = capsys.readouterr().err
    assert rc == 1
    assert err.startswith("error:")
    assert "Traceback" not in err


def test_main_debug_reraises(monkeypatch):
    with pytest.raises(FileNotFoundError):
        cli.main(["--debug", "validate", "/no/such/config.yaml"])
