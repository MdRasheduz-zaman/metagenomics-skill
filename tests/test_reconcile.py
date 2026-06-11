"""Unit tests for contig↔read reconciliation (`workflow/scripts/reconcile.py`).

Pure-Python parsing + join logic, loaded by path (it is a Snakemake `script:`). Covers the
fiddly real-format parsers (kraken2 --use-names contig calls, CAT add_names lineages) and the
concordance/flag/CAT-cross-check join in main(). Fixtures mirror real tool output.
"""
import csv
import importlib.util
import pathlib

_SPEC = importlib.util.spec_from_file_location(
    "reconcile",
    pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "reconcile.py",
)
rc = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rc)


def test_parse_contig_calls_names_taxids_and_paired_length(tmp_path):
    p = tmp_path / "c.kraken"
    p.write_text(
        "C\tcontig_1\tYellow fever virus, complete genome (taxid 1004)\t5977\t1004:50\n"
        "U\tcontig_2\tunclassified (taxid 0)\t1200\t0:30\n"
        "C\tcontig_3\tEscherichia coli (taxid 562)\t300|310\t562:5\n"   # paired length -> first
    )
    rows = {r["contig"]: r for r in rc.parse_contig_calls(str(p))}
    assert rows["contig_1"]["taxid"] == 1004
    assert rows["contig_1"]["taxon"] == "Yellow fever virus, complete genome"
    assert rows["contig_1"]["length"] == 5977 and rows["contig_1"]["status"] == "C"
    assert rows["contig_2"]["taxid"] == 0 and rows["contig_2"]["status"] == "U"
    assert rows["contig_3"]["length"] == 300


def test_parse_cat_extracts_taxid_and_species(tmp_path):
    p = tmp_path / "cat.named.txt"
    p.write_text(
        "# contig\tclassification\treason\tlineage\tscores\tsuperkingdom\t...\tspecies\n"
        "contig_1\ttaxid assigned\tbased on 1/1 ORFs\t1;1004\t1.00;1.00\tNA\t"
        "Yellow fever virus, complete genome: 1.00\n"
    )
    cat = rc.parse_cat(str(p))
    assert cat["contig_1"]["taxid"] == 1004
    assert cat["contig_1"]["taxon"] == "Yellow fever virus, complete genome"


def test_parse_read_report_species_only(tmp_path):
    p = tmp_path / "s.kreport"
    p.write_text(
        "50.0\t100\t100\tU\t0\tunclassified\n"
        " 4.10\t820\t10\tS\t1008\tPowassan virus\n"
        "10.0\t999\t0\tG\t11\tGenusRow\n"          # not rank S -> skipped
    )
    reads = rc.parse_read_report(str(p))
    assert set(reads) == {1008}
    assert reads[1008]["reads"] == 820 and reads[1008]["pct"] == 4.10


def _write(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def test_main_concordance_flags_and_cat_cross_check(tmp_path):
    import json
    contigs = tmp_path / "c.kraken"
    contigs.write_text(
        "C\tcontig_1\tPowassan virus (taxid 1008)\t9000\t1008:50\n"   # in reads too -> both
        "C\tcontig_2\tNovel bug (taxid 777)\t4000\t777:40\n"          # not in reads -> contigs_only
        "U\tcontig_3\tunclassified (taxid 0)\t2000\t0:20\n"           # high-cov unclassified -> flag
    )
    _write(tmp_path / "depth.txt", ["contigName", "contigLen", "totalAvgDepth"],
           [["contig_1", 9000, 5.0], ["contig_2", 4000, 3.0], ["contig_3", 2000, 8.0]])
    _write(tmp_path / "breadth.txt", ["#rname", "coverage"],
           [["contig_1", 100.0], ["contig_2", 95.0], ["contig_3", 100.0]])
    (tmp_path / "s.kreport").write_text(
        " 4.10\t820\t10\tS\t1008\tPowassan virus\n"
        " 2.00\t400\t5\tS\t999\tReads-only species\n"   # not assembled -> reads_only (pct>=1 -> flag)
    )
    (tmp_path / "cat.named.txt").write_text(
        "# header\n"
        "contig_1\ttaxid assigned\tx\t1;1008\t1.0\tNA\tPowassan virus: 1.00\n"  # agrees with kraken2
        "contig_2\ttaxid assigned\tx\t1;888\t1.0\tNA\tConflicting call: 1.00\n"  # conflicts (777 vs 888)
    )
    out = {k: str(tmp_path / f"{k}") for k in ("ct", "rt", "fl", "js")}
    rc.main(str(contigs), str(tmp_path / "depth.txt"), str(tmp_path / "s.kreport"),
            "samp", "confidence_0.1", out["ct"], out["rt"], out["fl"], out["js"],
            cat_named=str(tmp_path / "cat.named.txt"), breadth_path=str(tmp_path / "breadth.txt"))
    s = json.load(open(out["js"]))
    assert s["n_contigs"] == 3 and s["n_contigs_classified"] == 2
    assert s["taxa_concordance"] == {"both": 1, "reads_only": 1, "contigs_only": 1}
    # CAT cross-check: 2 CAT-classified, 1 agrees (contig_1), 1 conflicts (contig_2: 777 vs 888)
    assert s["cat_cross_check"] == {"n_cat_classified": 2, "kraken2_cat_agree": 1,
                                    "kraken2_cat_conflict": 1}
    flag_types = {r["type"] for r in csv.DictReader(open(out["fl"]), delimiter="\t",
                  fieldnames=["type", "k", "v", "e", "note"]) if r["type"] != "type"}
    assert {"contig_only_taxon", "reads_only_taxon", "unclassified_contig"} <= flag_types
