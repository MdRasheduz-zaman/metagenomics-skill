"""Prevalence-based decontamination (controls) — `workflow/scripts/decontam.py`.

Pure-Python, so directly testable. Covers the prevalence flag semantics and that the cleaned
table drops control samples + flagged taxa while retaining real biology.
"""
import csv
import importlib.util
import pathlib

import pytest

# decontam.py lives under workflow/scripts (not an installed package) — load it by path.
_SPEC = importlib.util.spec_from_file_location(
    "decontam",
    pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "decontam.py",
)
dc = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(dc)


def _rows(triples):
    return [{"sample": s, "name": n, "fraction_total_reads": f} for s, n, f in triples]


def test_flag_semantics_contaminant_vs_real():
    rows = _rows([
        ("blank", "Contaminant sp", "0.40"),
        ("blank", "Other reagent", "0.10"),     # control-only
        ("s1", "Contaminant sp", "0.05"), ("s1", "Real bug A", "0.50"),
        ("s2", "Contaminant sp", "0.04"), ("s2", "Real bug A", "0.60"),
        ("s3", "Real bug A", "0.70"),
    ])
    flags = dc.flag_contaminants(rows, controls={"blank"})
    # in control AND samples, control-prevalence >= sample-prevalence -> contaminant
    assert flags["Contaminant sp"]["flagged"] is True
    # control-only taxon: more prevalent in control than samples -> flagged (spec), but it is
    # absent from real samples so removing it is a no-op there.
    assert flags["Other reagent"]["flagged"] is True
    assert flags["Other reagent"]["sample_prevalence"] == 0.0
    # sample-only taxon never appears in a control -> never flagged
    assert flags["Real bug A"]["flagged"] is False


def test_run_cleaned_table_drops_controls_and_flagged(tmp_path):
    src = tmp_path / "bracken.tsv"
    with open(src, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "name", "fraction_total_reads"])
        for r in [("blank", "Contaminant sp", "0.40"),
                  ("s1", "Contaminant sp", "0.05"), ("s1", "Real bug A", "0.50"),
                  ("s2", "Contaminant sp", "0.04"), ("s2", "Real bug A", "0.60")]:
            w.writerow(r)
    res = dc.run(str(src), controls=["blank"],
                 out_flagged=str(tmp_path / "f.tsv"), out_cleaned=str(tmp_path / "c.tsv"))
    assert "Contaminant sp" in res["flagged"]
    with open(tmp_path / "c.tsv") as fh:
        cleaned = list(csv.DictReader(fh, delimiter="\t"))
    assert {r["sample"] for r in cleaned} == {"s1", "s2"}      # control removed
    assert {r["name"] for r in cleaned} == {"Real bug A"}      # contaminant removed, biology kept


def test_no_controls_flags_nothing():
    rows = _rows([("s1", "Bug", "0.5"), ("s2", "Bug", "0.5")])
    flags = dc.flag_contaminants(rows, controls=set())
    assert all(not d["flagged"] for d in flags.values())
