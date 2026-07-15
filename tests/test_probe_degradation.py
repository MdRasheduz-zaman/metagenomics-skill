"""Probe spec §3.3: a consented probe that can't read the data must DEGRADE to advisory
(measured: false) with a warning — not report a successful, empty measurement, and never infer
a platform from zero reads.
"""
import gzip

from metagx import probe


def _sheet(tmp_path, rows):
    p = tmp_path / "s.tsv"
    p.write_text("sample\tr1\tplatform\n" + "".join(f"{n}\t{r1}\t{plat}\n" for n, r1, plat in rows))
    return str(p)


def test_unreadable_paths_degrade_to_advisory(tmp_path, monkeypatch):
    monkeypatch.setattr(probe.consent, "set", lambda *a: "local")
    sheet = _sheet(tmp_path, [("s1", str(tmp_path / "nope.fastq.gz"), "ont")])
    r = probe.run(sheet, assume_yes=True)
    assert r["measured"] is False
    assert r["context"] == {"measured": False}
    assert any("s1" in u for u in r["unreadable"])


def test_empty_file_not_inferred_illumina(tmp_path, monkeypatch):
    monkeypatch.setattr(probe.consent, "set", lambda *a: "local")
    f = tmp_path / "empty.fastq.gz"
    with gzip.open(f, "wt"):
        pass                                          # 0 reads
    sheet = _sheet(tmp_path, [("s1", str(f), "ont")])
    r = probe.run(sheet, assume_yes=True)
    assert r["measured"] is False                     # was silently -> illumina before
    assert "s1" not in r.get("samples", {})


def test_partial_readable_measures_and_warns(tmp_path, monkeypatch):
    monkeypatch.setattr(probe.consent, "set", lambda *a: "local")
    good = tmp_path / "g.fastq"
    good.write_text("@r\nACGTACGTACGT\n+\nIIIIIIIIIIII\n" * 30)   # 30 reads
    sheet = _sheet(tmp_path, [("good", str(good), "illumina"),
                              ("bad", str(tmp_path / "missing.fq"), "ont")])
    r = probe.run(sheet, assume_yes=True)
    assert r["measured"] is True
    assert "good" in r["samples"] and "bad" not in r["samples"]
    assert any("skipped" in w for w in r["project"]["warnings"])
