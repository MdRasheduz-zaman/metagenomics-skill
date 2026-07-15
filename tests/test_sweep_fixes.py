"""Config-builder fixes found by the all-modules real-data sweep.

  * `db` is optional — DB-free module runs (phylogenetics-only, amplicon-only) must build.
  * filtered_assembly is single-end only (filter.smk); a paired sample must be rejected up front
    with a clear message, not die later with a cryptic Snakemake MissingInputException.
"""
import pytest

from metagx import config_builder as c


def test_build_config_db_optional_for_phylogenetics():
    cfg = c.build_config(
        project="p",
        samples=[{"sample": "x", "r1": "x.fasta", "platform": "ont", "layout": "se"}],
        modules={"classify": False, "abundance": False, "phylogenetics": True},
        phylogenetics={"input": "genomes.fasta", "method": "auto"})
    assert "db" in cfg and not cfg["db"].get("kraken2")   # no classifier DB needed


def test_build_config_db_optional_for_amplicon():
    cfg = c.build_config(
        project="p",
        samples=[{"sample": "a", "r1": "a_1.fq", "r2": "a_2.fq", "platform": "illumina",
                  "layout": "pe", "library": "amplicon"}],
        modules={"classify": False, "abundance": False})
    assert cfg["samples"][0]["library"] == "amplicon"


def test_filtered_assembly_rejects_paired():
    with pytest.raises(c.registry.ValidationError, match="single-end only"):
        c.build_config(
            project="p",
            samples=[{"sample": "a", "r1": "a_1.fq", "r2": "a_2.fq",
                      "platform": "illumina", "layout": "pe"}],
            db={"kraken2": "k"},
            modules={"assembly": True, "classify": True, "filtered_assembly": True})


def test_filtered_assembly_allows_single_end():
    cfg = c.build_config(
        project="p",
        samples=[{"sample": "a", "r1": "a.fq", "platform": "illumina", "layout": "se"}],
        db={"kraken2": "k"},
        modules={"assembly": True, "classify": True, "filtered_assembly": True})
    assert cfg["modules"]["filtered_assembly"] is True
