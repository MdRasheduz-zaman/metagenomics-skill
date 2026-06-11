"""Unit tests for read-vs-contig accuracy scoring (`workflow/scripts/read_contig_accuracy.py`).

The samtools-driven main() is verified by real execution (ROADMAP); here we pin the pure
parsers, the nodes.dmp lineage walk, and the read/contig bucketing — the scientific core.
"""
import importlib.util
import pathlib

_SPEC = importlib.util.spec_from_file_location(
    "read_contig_accuracy",
    pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "read_contig_accuracy.py",
)
rca = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(rca)

# Tiny taxonomy: species 1008 -> genus 11 -> root 1 ; 999 is an unrelated species.
PARENT = {1008: 11, 11: 1, 1: 1, 999: 1}


def test_parse_contig_and_read_taxa(tmp_path):
    ck = tmp_path / "c.kraken"
    ck.write_text(
        "C\tcontig_1\tPowassan virus (taxid 1008)\t9000\t1008:5\n"
        "U\tcontig_2\tunclassified (taxid 0)\t500\t0:5\n"
    )
    assert rca.parse_contig_taxa(str(ck)) == {"contig_1": 1008, "contig_2": 0}
    rk = tmp_path / "r.kraken"
    rk.write_text("C\tread_1\t1008\t250\t1008:50\nU\tread_2\t0\t250\t0:50\n")
    assert rca.parse_read_taxa(str(rk)) == {"read_1": 1008, "read_2": 0}


def test_lineage_walk():
    assert rca._is_ancestor(11, 1008, PARENT) is True       # genus is ancestor of species
    assert rca._is_ancestor(1008, 11, PARENT) is False
    assert rca._related(1008, 11, PARENT) is True
    assert rca._related(1008, 999, PARENT) is False


def test_bucket_all_five_categories():
    assert rca.bucket(1008, 1008, PARENT) == "concordant_exact"
    assert rca.bucket(1008, 11, PARENT) == "concordant_lineage"   # read=species, contig=genus
    assert rca.bucket(1008, 999, PARENT) == "discordant"
    assert rca.bucket(0, 1008, PARENT) == "read_unclassified"
    assert rca.bucket(1008, 0, PARENT) == "contig_unclassified"   # contig unclassified -> can't judge


def test_load_nodes(tmp_path):
    nd = tmp_path / "nodes.dmp"
    nd.write_text("1008\t|\t11\t|\tspecies\t|\n11\t|\t1\t|\tgenus\t|\n1\t|\t1\t|\tno rank\t|\n")
    parent = rca.load_nodes(str(nd))
    assert parent[1008] == 11 and parent[11] == 1
