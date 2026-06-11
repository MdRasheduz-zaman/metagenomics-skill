"""Unit tests for the kraken2-vs-second-classifier consensus parser/scoring helpers.

The logic lives in workflow/scripts/classifier_consensus.py (a Snakemake `script:`), so it
is loaded by path rather than imported as a package.
"""
import importlib.util
from pathlib import Path

import pytest

_SCRIPT = (Path(__file__).resolve().parents[1]
           / "workflow" / "scripts" / "classifier_consensus.py")
_spec = importlib.util.spec_from_file_location("classifier_consensus", _SCRIPT)
cc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cc)


def test_norm_unifies_label_styles():
    assert cc._norm("s__Escherichia_coli") == "escherichia coli"
    assert cc._norm("Escherichia coli K-12") == "escherichia coli"
    assert cc._norm("g__Bacteroides") == "bacteroides"


def test_parse_kraken_species(tmp_path):
    kr = tmp_path / "s.kreport"
    kr.write_text(
        "50.0\t100\t100\tU\t0\tunclassified\n"
        "10.0\t20\t20\tS\t562\tEscherichia coli\n"
        " 5.0\t10\t10\tS\t1280\tStaphylococcus aureus\n"
    )
    d = cc.parse_kraken_species(str(kr))
    assert d["escherichia coli"] == 10.0
    assert "staphylococcus aureus" in d
    assert "unclassified" not in d  # rank U is skipped


def test_parse_metaphlan(tmp_path):
    mp = tmp_path / "s.metaphlan.tsv"
    mp.write_text(
        "#mpa_v\n"
        "k__Bacteria\t2\t100.0\n"
        "k__Bacteria|p__Proteobacteria|...|s__Escherichia_coli\t562\t63.5\n"
        "k__Bacteria|...|s__Escherichia_coli|t__SGB\t0\t63.5\n"  # strain row ignored
    )
    d = cc.parse_metaphlan(str(mp))
    assert d == {"escherichia coli": 63.5}


def test_parse_kaiju_table_drops_nonspecies_rows(tmp_path):
    # Exact rows kaiju2table emits on a viral DB (verified by real run). Three are NOT species
    # and must be dropped: "unclassified", the DB-dependent "cannot be assigned to a (non-viral)
    # species", and the 0-read higher-rank catch-all "Viruses" (taxid 10239).
    kj = tmp_path / "s.kaiju.species.tsv"
    kj.write_text(
        "file\tpercent\treads\ttaxon_id\ttaxon_name\n"
        "s.kaiju\t4.235000\t847\t1008\tPowassan virus, complete genome\n"
        "s.kaiju\t3.680000\t736\t1001\tDengue virus 1, complete genome\n"
        "s.kaiju\t0.000000\t0\t10239\tViruses\n"
        "s.kaiju\t0.445000\t89\tNA\tcannot be assigned to a (non-viral) species\n"
        "s.kaiju\t51.670003\t10334\tNA\tunclassified\n"
    )
    d = cc.parse_kaiju_table(str(kj))
    assert len(d) == 2                    # only the two real species survive
    assert any("powassan" in k for k in d)
    assert any("dengue" in k for k in d)
    assert "cannot be" not in d           # prefix-matched and dropped
    assert "viruses" not in d             # 0-read catch-all dropped
    assert "unclassified" not in d


def test_concordance_metrics():
    kraken = {"escherichia coli": 60.0, "klebsiella pneumoniae": 30.0, "ghost taxon": 1.0}
    other = {"escherichia coli": 55.0, "klebsiella pneumoniae": 33.0}
    r = cc.concordance(kraken, other, "metaphlan", "s1")
    assert r["n_shared"] == 2
    assert r["jaccard"] == pytest.approx(2 / 3, abs=1e-3)
    assert "ghost taxon" in r["kraken2_only"]
    assert "escherichia coli" in r["top_overlap"]
