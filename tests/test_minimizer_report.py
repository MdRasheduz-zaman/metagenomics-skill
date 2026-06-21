"""kraken2 ``--report-minimizer-data`` changes the report's column layout: it inserts two
columns (distinct minimizers, distinct k-mers) between ``taxon_reads`` and ``rank``, so the
report has 8 columns instead of 6 and ``name`` moves from index 5 to 7.

Every kreport consumer must handle both layouts, or it silently misreads the rank/name when
minimizer reporting is on (e.g. Krona would treat the rank code ``S`` as a taxon name). These
tests pin that: the standard and minimizer-augmented reports must yield identical results from
``formats.kreport_row`` and from each downstream parser (Krona text, classifier consensus,
reconcile read report).
"""
import importlib.util
from pathlib import Path

from metagx import formats

REPO = Path(__file__).resolve().parents[1]


def _load(script: str):
    """Import a workflow script module by path (they're not a package)."""
    path = REPO / "workflow" / "scripts" / script
    spec = importlib.util.spec_from_file_location(script[:-3], path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Standard 6-column report (pct, clade_reads, taxon_reads, rank, taxid, name).
STD = (
    " 50.00\t1000\t1000\tU\t0\tunclassified\n"
    " 50.00\t1000\t0\tR\t1\troot\n"
    " 30.00\t600\t600\tS\t1001\t  Dengue virus 1\n"
    " 20.00\t400\t400\tS\t1002\t  Zika virus\n"
)
# Same data with --report-minimizer-data: two extra columns after taxon_reads (8 total).
MIN = (
    " 50.00\t1000\t1000\t0\t0\tU\t0\tunclassified\n"
    " 50.00\t1000\t0\t1200\t900\tR\t1\troot\n"
    " 30.00\t600\t600\t800\t640\tS\t1001\t  Dengue virus 1\n"
    " 20.00\t400\t400\t500\t410\tS\t1002\t  Zika virus\n"
)


def test_kreport_row_handles_both_layouts():
    for text in (STD, MIN):
        rows = [formats.kreport_row(l) for l in text.splitlines()]
        species = [r for r in rows if r and r["rank"] == "S"]
        names = {r["name"].strip() for r in species}
        assert names == {"Dengue virus 1", "Zika virus"}
        # taxid/reads read correctly regardless of the inserted minimizer columns
        by_taxid = {r["taxid"]: r for r in species}
        assert by_taxid["1001"]["taxon_reads"] == "600"
        assert by_taxid["1002"]["taxid"] == "1002"


def test_krona_identical_for_both_layouts(tmp_path):
    krona = _load("kreport2krona.py")
    out_std, out_min = tmp_path / "std.txt", tmp_path / "min.txt"
    (tmp_path / "std.kreport").write_text(STD)
    (tmp_path / "min.kreport").write_text(MIN)
    n_std = krona.convert(str(tmp_path / "std.kreport"), str(out_std))
    n_min = krona.convert(str(tmp_path / "min.kreport"), str(out_min))
    assert n_std == n_min > 0
    assert out_std.read_text() == out_min.read_text()
    # the rank code must never leak in as a taxon name
    assert "\tS\n" not in out_min.read_text() and "\tS\t" not in out_min.read_text()
    assert "Dengue virus 1" in out_min.read_text()


def test_consensus_and_reconcile_identical_for_both_layouts(tmp_path):
    cc = _load("classifier_consensus.py")
    rc = _load("reconcile.py")
    (tmp_path / "std.kreport").write_text(STD)
    (tmp_path / "min.kreport").write_text(MIN)

    assert cc.parse_kraken_species(str(tmp_path / "std.kreport")) == \
        cc.parse_kraken_species(str(tmp_path / "min.kreport"))

    rr_std = rc.parse_read_report(str(tmp_path / "std.kreport"))
    rr_min = rc.parse_read_report(str(tmp_path / "min.kreport"))
    assert rr_std == rr_min
    assert set(rr_min) == {1001, 1002}
    assert rr_min[1001]["name"] == "Dengue virus 1"
