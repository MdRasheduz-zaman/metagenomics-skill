"""Custom database builders (kraken2 / CAT / Kaiju) — no NCBI download.

The plan-only tests need no external tools. The real-build test runs prodigal +
kaiju-mkbwt/mkfmi when present (skips in CI), verifying the Kaiju index that the
consensus module's db.kaiju consumes.
"""
import os
import shutil

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
