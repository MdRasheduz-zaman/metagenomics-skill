"""BLAST validation module: the blastn registry, validate.py core logic, and config wiring.

A kraken2/Bracken call is not proof; the `validate` module BLASTs a subsample of reads for the
top taxa and checks the best alignment's organism against the classifier. These tests pin the
parsing + agreement logic (no BLAST needed) and the config/db gating.
"""
import os
import shutil
import subprocess

import pytest

from metagx import config_builder as cb
from metagx import registry
from metagx import validation as validate

_FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "viral")

SAMPLES = [{"sample": "s1", "r1": "x.fastq", "platform": "illumina", "layout": "se"}]


# --- registry ----------------------------------------------------------------------------
def test_blastn_registry_renders_with_managed_outfmt_intact():
    assert "blastn" in registry.list_tools()
    args = registry.render_args(
        "blastn", {"evalue": 1e-10, "perc_identity": 90, "max_target_seqs": 5},
        managed={"query": "q.fa", "db": "nt", "out": "o.tsv", "num_threads": 4,
                 "outfmt": "6 qseqid sseqid sscinames"})
    assert "-evalue" in args and "-query" in args
    # the multi-token outfmt stays one argv element (subprocess-safe, no shell quoting)
    assert "6 qseqid sseqid sscinames" in args


# --- kraken / kreport parsing ------------------------------------------------------------
def test_parse_kraken_assignments_both_id_styles():
    text = (
        "C\tread1\t562\t150\t562:100\n"               # bare taxid
        "C\tread2\tEscherichia coli (taxid 562)\t150\t...\n"  # --use-names style
        "U\tread3\t0\t150\t0:150\n"                    # unclassified -> skipped
    )
    a = validate.parse_kraken_assignments(text)
    assert a == {"read1": "562", "read2": "562"}


def test_top_taxa_reads_kreport(tmp_path):
    kr = tmp_path / "s.kreport"
    kr.write_text(
        "50.0\t500\t500\tS\t562\tEscherichia coli\n"
        "30.0\t300\t300\tS\t1280\tStaphylococcus aureus\n"
        "10.0\t100\t100\tG\t561\tEscherichia\n"
    )
    taxa = validate.top_taxa(str(kr), level="S", top_n=10)
    assert [t["taxid"] for t in taxa] == ["562", "1280"]  # species only, by reads desc
    assert taxa[0]["name"] == "Escherichia coli"


# --- agreement logic ---------------------------------------------------------------------
def test_names_agree_genus_and_species():
    assert validate.names_agree("Escherichia coli", "Escherichia coli strain K-12", "genus")
    assert validate.names_agree("Escherichia coli", "Escherichia coli O157:H7", "species")
    # genus matches but species differs -> species-level disagreement
    assert validate.names_agree("Escherichia coli", "Escherichia albertii", "genus")
    assert not validate.names_agree("Escherichia coli", "Escherichia albertii", "species")
    # different genus -> disagree at any level
    assert not validate.names_agree("Escherichia coli", "Salmonella enterica", "genus")


def test_best_hit_and_assess():
    rows = validate.parse_blast6(
        "r1\tNC_1\t99.0\t150\t1e-50\t280\t562\tEscherichia coli\tE. coli genome\n"
        "r1\tNC_2\t90.0\t150\t1e-20\t180\t28901\tSalmonella enterica\tSalmonella\n"
        "r2\tNC_3\t99.0\t150\t1e-50\t290\t1280\tStaphylococcus aureus\tS. aureus\n"
    )
    best = validate.best_hit_per_query(rows)
    assert best["r1"]["sseqid"] == "NC_1"  # higher bitscore wins
    a = validate.assess({"r1": "Escherichia coli", "r2": "Escherichia coli", "r3": "Escherichia coli"},
                        best, level="genus")
    assert a["n_queries"] == 3 and a["n_with_hits"] == 2
    assert a["n_agree"] == 1            # r1 agrees (E. coli), r2 is S. aureus, r3 no hit
    assert a["agreement_rate"] == 0.5   # over queries WITH hits
    assert a["hit_rate"] == round(2 / 3, 4)


def test_assess_falls_back_to_stitle_when_sscinames_na():
    # No NCBI taxdb installed -> blastn writes "N/A" in sscinames; the organism is in stitle.
    rows = validate.parse_blast6(
        "r1\tNC_1\t99.0\t150\t1e-50\t280\tN/A\tN/A\tYellow fever virus, complete genome\n")
    best = validate.best_hit_per_query(rows)
    a = validate.assess({"r1": "Yellow fever virus"}, best, level="genus")
    assert a["n_agree"] == 1  # matched via stitle, not the "N/A" sscinames
    assert a["per_query"][0]["blast"].startswith("Yellow fever")


def test_verdict_bands():
    assert validate.verdict(0.9, 0.8) == "corroborated"
    assert validate.verdict(0.6, 0.8) == "partial"
    assert validate.verdict(0.1, 0.8) == "discordant"
    assert validate.verdict(0.9, 0.1) == "inconclusive"  # too few hits to judge


# --- read extraction ---------------------------------------------------------------------
def test_extract_sequences_fastq_and_paired_suffix(tmp_path):
    fq = tmp_path / "r.fastq"
    fq.write_text("@read1/1 desc\nACGT\n+\nIIII\n@read2/1\nTTTT\n+\nIIII\n")
    got = validate.extract_sequences([str(fq)], {"read1", "read2"})
    assert got == {"read1": "ACGT", "read2": "TTTT"}  # /1 suffix + description stripped


# --- config gating -----------------------------------------------------------------------
def test_validate_requires_blast_db_or_remote():
    with pytest.raises(registry.ValidationError) as e:
        cb.build_config(samples=SAMPLES, db={"kraken2": "DB"},
                        modules={"validate": True})
    assert "db.blast" in str(e.value)


def test_validate_accepts_local_db():
    cfg = cb.build_config(samples=SAMPLES, db={"kraken2": "DB", "blast": "/blast/nt"},
                          modules={"validate": True})
    assert cfg["validate"]["target"] == "reads"
    assert cfg["db"]["blast"] == "/blast/nt"
    assert cfg["modules"]["validate"] is True


def test_validate_remote_needs_no_local_db():
    cfg = cb.build_config(samples=SAMPLES, db={"kraken2": "DB"},
                          modules={"validate": True}, validate={"remote": True})
    assert cfg["validate"]["remote"] is True


def test_validate_needs_classify():
    with pytest.raises(registry.ValidationError) as e:
        cb.build_config(samples=SAMPLES, db={"kraken2": "DB", "blast": "/b"},
                        modules={"validate": True, "classify": False})
    assert "classify" in str(e.value)


# --- gated real-BLAST e2e (skips when BLAST+ / fixture absent) ----------------------------
@pytest.mark.skipif(
    not (shutil.which("blastn") and shutil.which("makeblastdb")
         and os.path.isdir(_FIXTURE)),
    reason="needs BLAST+ on PATH and the viral fixture")
def test_blast_validate_end_to_end_distinguishes_right_from_wrong(tmp_path):
    """Build a custom BLAST DB from the viral genomes, BLAST a real read, and confirm the
    agreement check corroborates the TRUE organism and flags a WRONG label as discordant."""
    db = str(tmp_path / "viraldb")
    subprocess.run(["makeblastdb", "-in", os.path.join(_FIXTURE, "genomes.fasta"),
                    "-dbtype", "nucl", "-out", db, "-parse_seqids"],
                   check=True, capture_output=True, text=True)
    rid = "read_1_forward_ref4"  # ref4 == 4th genome == Yellow fever virus
    seqs = validate.extract_sequences([os.path.join(_FIXTURE, "ont_reads.fasta")], {rid})
    qfa = str(tmp_path / "q.fasta")
    assert validate.write_fasta({rid: seqs[rid]}, qfa) == 1
    out = str(tmp_path / "blast.tsv")
    outfmt = "6 " + " ".join(validate.BLAST6_FIELDS)
    argv = ["blastn"] + registry.render_args(
        "blastn", {"evalue": 1e-5, "task": "blastn"},
        managed={"query": qfa, "db": db, "out": out, "num_threads": 2, "outfmt": outfmt})
    subprocess.run(argv, check=True, capture_output=True, text=True)
    hits = validate.best_hit_per_query(validate.parse_blast6(open(out).read()))
    assert rid in hits, "the viral read should hit its source genome"

    right = validate.assess({rid: "Yellow fever virus"}, hits, level="genus")
    wrong = validate.assess({rid: "Dengue virus"}, hits, level="genus")
    assert right["agreement_rate"] == 1.0
    assert validate.verdict(right["agreement_rate"], right["hit_rate"]) == "corroborated"
    assert wrong["n_agree"] == 0  # mislabel caught (the whole point of validation)
