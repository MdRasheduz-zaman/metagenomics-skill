"""Sample-sheet validation now covers TSV paths (H3), rejects dotted ids (M1), and allows
contigs-only samples — so a bad sheet fails `metagx validate` clearly instead of as a raw
Snakemake traceback at run time.
"""
import pytest

from metagx import config_builder as c


def test_tsv_bad_platform_rejected(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("sample\tr1\tplatform\ns1\tr.fq.gz\tbadplat\n")
    with pytest.raises(c.registry.ValidationError, match="unknown platform"):
        c._validate_samples(str(s))


def test_tsv_long_read_with_r2_rejected(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("sample\tr1\tr2\tplatform\ns1\ta.fq\tb.fq\tont\n")
    with pytest.raises(c.registry.ValidationError, match="single-end"):
        c._validate_samples(str(s))


def test_tsv_good_passes(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("sample\tr1\tplatform\tlayout\ns1\tr.fq.gz\tont\tse\n")
    assert c._validate_samples(str(s)) == str(s)


def test_tsv_bom_tolerated(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("﻿sample\tr1\tplatform\ns1\tr.fq.gz\tont\n", encoding="utf-8")
    assert c._validate_samples(str(s)) == str(s)   # Excel BOM, no KeyError


def test_dotted_sample_id_rejected_inline():
    with pytest.raises(c.registry.ValidationError, match="contains"):
        c._validate_samples([{"sample": "pat.01", "r1": "a.fq", "platform": "ont"}])


def test_dotted_sample_id_rejected_tsv(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("sample\tr1\tplatform\npat.01\tr.fq.gz\tont\n")
    with pytest.raises(c.registry.ValidationError, match="contains"):
        c._validate_samples(str(s))


def test_contigs_only_sample_allowed():
    assert c._validate_samples([{"sample": "g", "contigs": "genome.fa", "platform": "illumina"}])


def test_grouped_tsv_without_group_on_every_row_ok(tmp_path):
    # empty 'group' cells must NOT be rejected (a common differential sheet shape)
    s = tmp_path / "s.tsv"
    s.write_text("sample\tr1\tplatform\tgroup\na\ta.fq\tont\tcase\nb\tb.fq\tont\t\n")
    assert c._validate_samples(str(s)) == str(s)
