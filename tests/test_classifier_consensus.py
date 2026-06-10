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


def test_concordance_metrics():
    kraken = {"escherichia coli": 60.0, "klebsiella pneumoniae": 30.0, "ghost taxon": 1.0}
    other = {"escherichia coli": 55.0, "klebsiella pneumoniae": 33.0}
    r = cc.concordance(kraken, other, "metaphlan", "s1")
    assert r["n_shared"] == 2
    assert r["jaccard"] == pytest.approx(2 / 3, abs=1e-3)
    assert "ghost taxon" in r["kraken2_only"]
    assert "escherichia coli" in r["top_overlap"]
