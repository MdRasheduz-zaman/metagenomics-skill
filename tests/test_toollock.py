"""Tool version + valid-arg locking: drift diff + config-aware flag validation.

The diff/check logic is pure (operates on dicts + a stubbed --help capture), so it runs with
no bio tools installed. A renamed/removed flag or a vanished tool must be caught.
"""
import json

from metagx import registry, sync_help, toollock


def test_conda_version_reads_owning_package(tmp_path):
    """conda_version resolves a binary's version from the conda-meta record that owns it — the
    fallback for tools with no usable --version flag (FastTree/Kaiju/Krona/vsearch)."""
    (tmp_path / "bin").mkdir()
    exe = tmp_path / "bin" / "FastTree"
    exe.write_text("")
    (tmp_path / "conda-meta").mkdir()
    (tmp_path / "conda-meta" / "fasttree-2.2.0-h7b50bb2_1.json").write_text(
        json.dumps({"version": "2.2.0", "files": ["bin/FastTree", "share/x"]}))
    assert sync_help.conda_version(str(exe)) == "2.2.0"


def test_conda_version_none_when_not_a_conda_install(tmp_path):
    (tmp_path / "bin").mkdir()
    exe = tmp_path / "bin" / "tool"
    exe.write_text("")  # no conda-meta dir in the prefix
    assert sync_help.conda_version(str(exe)) is None
    assert sync_help.conda_version(None) is None


def test_conda_version_matches_by_basename_for_symlinked_binaries(tmp_path):
    """A wrapper binary whose conda record lists the file under a different path still matches."""
    (tmp_path / "bin").mkdir()
    exe = tmp_path / "bin" / "ktImportText"
    exe.write_text("")
    (tmp_path / "conda-meta").mkdir()
    (tmp_path / "conda-meta" / "krona-2.8.1-pl5321hdfd78af_1.json").write_text(
        json.dumps({"version": "2.8.1", "files": ["opt/krona/bin/ktImportText"]}))
    assert sync_help.conda_version(str(exe)) == "2.8.1"


def _cap(version, flags):
    """A fake sync_help.capture_help result."""
    def _f(_command):
        return {"ok": True, "version": version, "flags": [{"flag": x} for x in flags]}
    return _f


# --- drift between a lock and the current install ----------------------------------------
def test_diff_lock_flags_removed_is_error():
    locked = {"tools": {"kraken2": {"command": "kraken2", "installed": True,
                                    "version": "2.1.3", "flags": ["--confidence", "--threads"]}}}
    current = {"tools": {"kraken2": {"command": "kraken2", "installed": True, "help_ok": True,
                                     "version": "2.1.3", "flags": ["--threads"]}}}  # lost --confidence
    drift = toollock.diff_lock(locked, current)
    assert any(d["kind"] == "flags_removed" and d["severity"] == "error" for d in drift)


def test_diff_lock_missing_tool_is_error():
    locked = {"tools": {"flye": {"command": "flye", "installed": True, "version": "2.9", "flags": []}}}
    current = {"tools": {"flye": {"command": "flye", "installed": False}}}
    drift = toollock.diff_lock(locked, current)
    assert drift and drift[0]["kind"] == "missing" and drift[0]["severity"] == "error"


def test_diff_lock_version_change_is_warn_not_error():
    locked = {"tools": {"fastp": {"command": "fastp", "installed": True,
                                  "version": "0.23.4", "flags": ["--thread"]}}}
    current = {"tools": {"fastp": {"command": "fastp", "installed": True, "help_ok": True,
                                   "version": "1.0.0", "flags": ["--thread"]}}}
    drift = toollock.diff_lock(locked, current)
    assert drift and drift[0]["kind"] == "version" and drift[0]["severity"] == "warn"


def test_diff_lock_clean_when_identical():
    snap = {"tools": {"kraken2": {"command": "kraken2", "installed": True, "help_ok": True,
                                  "version": "2.1.3", "flags": ["--confidence"]}}}
    assert toollock.diff_lock(snap, snap) == []


# --- config-aware flag validity against the installed --help -----------------------------
def test_config_flag_check_flags_unknown_flag(monkeypatch):
    # kraken2 on PATH, but its --help no longer lists --confidence -> the config using it must fail.
    monkeypatch.setattr(toollock.shutil, "which", lambda _exe: "/usr/bin/kraken2")
    cfg = {"kraken2": {"confidence": 0.1}}
    findings = toollock.config_flag_check(cfg, capture=_cap("2.9", ["--threads", "--db"]))
    assert findings and findings[0]["tool"] == "kraken2"
    assert "--confidence" in findings[0]["message"]


def test_config_flag_check_passes_when_flag_present(monkeypatch):
    monkeypatch.setattr(toollock.shutil, "which", lambda _exe: "/usr/bin/kraken2")
    cfg = {"kraken2": {"confidence": 0.1}}
    findings = toollock.config_flag_check(cfg, capture=_cap("2.9", ["--confidence", "--threads"]))
    assert findings == []


def test_config_flag_check_skips_tools_not_on_path(monkeypatch):
    monkeypatch.setattr(toollock.shutil, "which", lambda _exe: None)
    cfg = {"kraken2": {"confidence": 0.1}}
    assert toollock.config_flag_check(cfg, capture=_cap("x", [])) == []


def test_probe_tool_shape(monkeypatch):
    monkeypatch.setattr(toollock.shutil, "which", lambda _exe: "/usr/bin/kraken2")
    p = toollock.probe_tool("kraken2", capture=_cap("2.1.3", ["--confidence"]))
    assert p["tool"] == "kraken2" and p["installed"] and p["version"] == "2.1.3"
    assert "--confidence" in p["flags"]
