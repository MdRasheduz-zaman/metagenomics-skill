"""Database onboarding (`metagx fetch-db`, metagx/dbfetch.py).

`plan()` is pure so it's fully tested offline; `fetch()` is exercised with a fake downloader
that writes the expected artifacts. A network-gated test (METAGX_NET_TESTS=1) HEAD-checks the
real URLs so link-rot — like the dead Snakefile default this work uncovered — is caught.
"""
import os
import urllib.request

import pytest

from metagx import dbfetch


def test_describe_lists_all_indices():
    names = {d["name"] for d in dbfetch.describe()}
    assert names == set(dbfetch.INDICES)
    for d in dbfetch.describe():
        assert d["url"].startswith("https://") and d["url"].endswith(".tar.gz")
        assert d["size"] and d["description"]


def test_default_index_is_known():
    assert dbfetch.DEFAULT in dbfetch.INDICES


def test_index_url_unknown_raises():
    with pytest.raises(KeyError):
        dbfetch.index_url("does-not-exist")


def test_plan_is_pure_and_curls_to_dir(tmp_path):
    p = dbfetch.plan("standard-8", str(tmp_path / "kdb"))
    assert p["url"].endswith("k2_standard_08gb_20241228.tar.gz")
    assert p["db"] == os.path.abspath(str(tmp_path / "kdb"))
    assert "curl" in p["command"] and "tar -xzf" in p["command"]
    assert p["config_hint"]["db"]["kraken2"] == p["db"]


def test_plan_custom_url_overrides_name(tmp_path):
    p = dbfetch.plan("standard-8", str(tmp_path / "k"), url="https://example.org/custom.tar.gz")
    assert p["url"] == "https://example.org/custom.tar.gz"


def test_fetch_dry_run_does_not_download(tmp_path):
    r = dbfetch.fetch("viral", str(tmp_path / "k"), run=False)
    assert r["ran"] is False
    assert not os.path.exists(os.path.join(str(tmp_path / "k"), "hash.k2d"))


def test_fetch_reuses_existing_built_index(tmp_path):
    db = tmp_path / "k"
    db.mkdir()
    for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (db / f).write_bytes(b"x")
    r = dbfetch.fetch("viral", str(db))
    assert r["ran"] is False and r["ok"] is True
    assert "reusing" in r["note"]


def test_fetch_success_when_artifacts_appear(tmp_path, monkeypatch):
    db = tmp_path / "k"

    def fake_run(cmd, capture_output, text):
        os.makedirs(str(db), exist_ok=True)
        for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
            open(os.path.join(str(db), f), "wb").write(b"x")

        class P:
            returncode, stdout, stderr = 0, "", ""
        return P()

    monkeypatch.setattr(dbfetch.shutil, "which", lambda _: "/usr/bin/x")
    monkeypatch.setattr(dbfetch.subprocess, "run", fake_run)
    r = dbfetch.fetch("viral", str(db))
    assert r["ran"] is True and r["ok"] is True


def test_fetch_fails_when_no_artifacts(tmp_path, monkeypatch):
    db = tmp_path / "k"

    def fake_run(cmd, capture_output, text):
        class P:
            returncode, stdout, stderr = 1, "", "curl: (6) could not resolve host"
        return P()

    monkeypatch.setattr(dbfetch.shutil, "which", lambda _: "/usr/bin/x")
    monkeypatch.setattr(dbfetch.subprocess, "run", fake_run)
    r = dbfetch.fetch("viral", str(db))
    assert r["ran"] is True and r["ok"] is False
    assert "missing" in r["note"]


@pytest.mark.skipif(os.environ.get("METAGX_NET_TESTS") != "1",
                    reason="network test; set METAGX_NET_TESTS=1 to check live DB URLs")
@pytest.mark.parametrize("name", sorted(dbfetch.INDICES))
def test_index_url_is_live(name):
    req = urllib.request.Request(dbfetch.index_url(name), method="HEAD")
    with urllib.request.urlopen(req, timeout=30) as resp:
        assert resp.status == 200
        assert int(resp.headers.get("Content-Length", "0")) > 0
