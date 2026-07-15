"""`metagx tag-taxids` helper + the synthetic-taxonomy validate-default fix (RT-1).

The real-taxonomy path needs `>seqid|kraken:taxid|<taxid>` headers; this locks in the tagging
(offline map + online resolution) and the graceful-degradation defaults for a synthetic DB.
"""
import json

import pytest

from metagx import config_builder, taxid_tag, validation


# --------------------------------------------------------------------------- #
# header parsing                                                              #
# --------------------------------------------------------------------------- #
def test_parse_accession_plain():
    assert taxid_tag.parse_accession(">NC_001477.1 Dengue virus 1, complete genome") == "NC_001477.1"


def test_parse_accession_already_tagged():
    assert taxid_tag.parse_accession(">NC_001477.1|kraken:taxid|11053 Dengue") == "NC_001477.1"


def test_already_tagged():
    assert taxid_tag.already_tagged(">x|kraken:taxid|9") is True
    assert taxid_tag.already_tagged(">NC_1.1 name") is False


# --------------------------------------------------------------------------- #
# offline map                                                                 #
# --------------------------------------------------------------------------- #
def test_load_map_tolerates_header_comments_ws(tmp_path):
    p = tmp_path / "m.tsv"
    p.write_text("accession\ttaxid\n# a comment\nNC_001477.1\t11053\nNC_012532.1  64320\n\n")
    m = taxid_tag.load_map(str(p))
    assert m == {"NC_001477.1": "11053", "NC_012532.1": "64320"}


def test_tag_fasta_offline_and_idempotent(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">NC_001477.1 Dengue\nACGT\n>NC_012532.1 Zika\nTTTT\n")
    out = tmp_path / "out.fa"
    mp = {"NC_001477.1": "11053", "NC_012532.1": "64320"}
    s = taxid_tag.tag_fasta(str(fa), str(out), mp)
    assert s["n_tagged"] == 2 and s["n_missing"] == 0
    text = out.read_text()
    assert ">NC_001477.1|kraken:taxid|11053 Dengue" in text
    assert "ACGT" in text  # sequence untouched
    # re-tagging the output is a no-op (idempotent)
    out2 = tmp_path / "out2.fa"
    s2 = taxid_tag.tag_fasta(str(out), str(out2), mp)
    assert s2["n_already"] == 2 and s2["n_tagged"] == 0


def test_tag_fasta_version_stripped_lookup(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">NC_001477.1 Dengue\nACGT\n")
    out = tmp_path / "out.fa"
    s = taxid_tag.tag_fasta(str(fa), str(out), {"NC_001477": "11053"})  # map lacks the .1 version
    assert s["n_tagged"] == 1
    assert "kraken:taxid|11053" in out.read_text()


def test_tag_fasta_missing_fails_loud(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">NC_999999.9 unknown\nACGT\n")
    with pytest.raises(ValueError, match="no taxid"):
        taxid_tag.tag_fasta(str(fa), str(tmp_path / "o.fa"), {"NC_001477.1": "11053"})


def test_tag_fasta_allow_missing(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">NC_999999.9 unknown\nACGT\n")
    out = tmp_path / "o.fa"
    s = taxid_tag.tag_fasta(str(fa), str(out), {}, allow_missing=True)
    assert s["n_missing"] == 1 and s["n_tagged"] == 0
    assert ">NC_999999.9 unknown" in out.read_text()  # left untagged, sequence preserved


# --------------------------------------------------------------------------- #
# online resolution (mocked transport — no real network in the suite)         #
# --------------------------------------------------------------------------- #
def test_resolve_online_parses_esummary(monkeypatch):
    canned = {"result": {"uids": ["9626685"],
                         "9626685": {"accessionversion": "NC_001477.1",
                                     "caption": "NC_001477", "taxid": 11053}}}
    monkeypatch.setattr(taxid_tag, "_http_get", lambda url, timeout=30: json.dumps(canned))
    monkeypatch.setattr(taxid_tag.time, "sleep", lambda *_: None)
    mapping, failures = taxid_tag.resolve_online(["NC_001477.1"])
    assert failures == []
    assert mapping["NC_001477.1"] == "11053" and mapping["NC_001477"] == "11053"


def test_resolve_online_collects_failures(monkeypatch):
    def boom(url, timeout=30):
        raise OSError("network down")
    monkeypatch.setattr(taxid_tag, "_http_get", boom)
    monkeypatch.setattr(taxid_tag.time, "sleep", lambda *_: None)
    mapping, failures = taxid_tag.resolve_online(["NC_1.1", "NC_2.1"])
    assert mapping == {} and set(failures) == {"NC_1.1", "NC_2.1"}


# --------------------------------------------------------------------------- #
# ranks_present + RT-1: validate default rank follows the DB taxonomy         #
# --------------------------------------------------------------------------- #
def test_ranks_present(tmp_path):
    kr = tmp_path / "k.kreport"
    kr.write_text(" 50\t10\t10\tU\t0\tunclassified\n 50\t10\t1\tR\t1\troot\n"
                  "  4\t9\t9\tS\t1008\t  Powassan\n")
    assert validation.ranks_present(str(kr)) == ["R", "S", "U"]


def _cfg(db_build, validate_block=None):
    return config_builder.build_config(
        project="t", outdir="results", threads=4,
        samples=[{"sample": "s", "r1": "s.fastq.gz", "platform": "ont", "layout": "se"}],
        db={"build": db_build}, modules={"classify": True, "abundance": True, "validate": True},
        validate=validate_block)


def test_synthetic_db_defaults_validate_to_species():
    cfg = _cfg({"strategy": "custom-fasta", "source": "g.fasta", "taxonomy": "synthetic",
                "blast": True})
    assert cfg["validate"]["level"] == "species" and cfg["validate"]["rank"] == "S"


def test_custom_real_db_also_defaults_to_species():
    # a custom DB (user's own genomes) validates the exact calls at species — reads are assigned
    # at the leaf and names are often non-binomial (viruses), so genus name-matching is unreliable.
    cfg = _cfg({"strategy": "custom-fasta", "source": "g.fasta", "taxonomy": "real", "blast": True})
    assert cfg["validate"]["level"] == "species" and cfg["validate"]["rank"] == "S"


def test_prebuilt_db_keeps_genus_default():
    # no db.build (a broad prebuilt/standard DB) keeps genus — robust to species-level BLAST ambiguity.
    cfg = config_builder.build_config(
        project="t", outdir="results", threads=4,
        samples=[{"sample": "s", "r1": "s.fastq.gz", "platform": "ont", "layout": "se"}],
        db={"kraken2": "k", "blast": "b"},
        modules={"classify": True, "abundance": True, "validate": True})
    assert cfg["validate"]["level"] == "genus" and cfg["validate"]["rank"] == "G"


def test_explicit_level_overrides_default():
    cfg = _cfg({"strategy": "custom-fasta", "source": "g.fasta", "taxonomy": "synthetic",
                "blast": True}, validate_block={"level": "genus"})
    assert cfg["validate"]["level"] == "genus"
