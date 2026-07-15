"""Cross-sample checks (differential groups, decontam control, damage ancient) now enforce for a
TSV sample sheet too — matching the structure the tool generates — not just inline lists. An
absent (not-yet-created) sheet still defers to run time.
"""
import pytest

from metagx import config_builder as c

_DB = {"kraken2": "k"}


def _sheet(tmp_path, text):
    p = tmp_path / "s.tsv"
    p.write_text(text)
    return str(p)


def _build(tmp_path, text, modules, **extra):
    return c.build_config(project="t", samples=_sheet(tmp_path, text), db=_DB,
                          modules=modules, **extra)


# --- differential ---------------------------------------------------------- #
_DIFF_OK = ("sample\tr1\tplatform\tgroup\n"
            "a\ta.fq.gz\tont\tcase\nb\tb.fq.gz\tont\tcase\n"
            "c\tc.fq.gz\tont\tcontrol\nd\td.fq.gz\tont\tcontrol\n")


def test_differential_tsv_valid_passes(tmp_path):
    assert _build(tmp_path, _DIFF_OK, {"abundance": True, "differential": True})


def test_differential_tsv_one_group_rejected(tmp_path):
    one = ("sample\tr1\tplatform\tgroup\n"
           "a\ta.fq.gz\tont\tcase\nb\tb.fq.gz\tont\tcase\n")
    with pytest.raises(c.registry.ValidationError, match=">=2 sample groups"):
        _build(tmp_path, one, {"abundance": True, "differential": True})


def test_differential_tsv_under_replicated_rejected(tmp_path):
    ur = ("sample\tr1\tplatform\tgroup\n"
          "a\ta.fq.gz\tont\tcase\nb\tb.fq.gz\tont\tcase\nc\tc.fq.gz\tont\tcontrol\n")
    with pytest.raises(c.registry.ValidationError, match="under-replicated"):
        _build(tmp_path, ur, {"abundance": True, "differential": True})


def test_differential_tsv_custom_group_column(tmp_path):
    ok = ("sample\tr1\tplatform\ttreatment\n"
          "a\ta.fq.gz\tont\tx\nb\tb.fq.gz\tont\tx\nc\tc.fq.gz\tont\ty\nd\td.fq.gz\tont\ty\n")
    assert _build(tmp_path, ok, {"abundance": True, "differential": True},
                  differential={"group_column": "treatment"})


# --- decontam -------------------------------------------------------------- #
def test_decontam_tsv_without_control_rejected(tmp_path):
    no_ctrl = "sample\tr1\tplatform\na\ta.fq.gz\tont\nb\tb.fq.gz\tont\n"
    with pytest.raises(c.registry.ValidationError, match="control: true"):
        _build(tmp_path, no_ctrl, {"abundance": True, "decontam": True})


def test_decontam_tsv_with_control_passes(tmp_path):
    ctrl = "sample\tr1\tplatform\tcontrol\na\ta.fq.gz\tont\t\nb\tb.fq.gz\tont\ttrue\n"
    assert _build(tmp_path, ctrl, {"abundance": True, "decontam": True})


# --- damage (aDNA) --------------------------------------------------------- #
def test_damage_tsv_without_ancient_rejected(tmp_path):
    no_anc = "sample\tr1\tplatform\tlibrary\na\ta.fq.gz\tillumina\twgs\n"
    with pytest.raises(c.registry.ValidationError, match="library=ancient"):
        _build(tmp_path, no_anc, {"classify": True, "abundance": True,
                                  "assembly": True, "damage": True})


def test_damage_tsv_with_ancient_passes(tmp_path):
    anc = "sample\tr1\tplatform\tlibrary\na\ta.fq.gz\tillumina\tancient\n"
    assert _build(tmp_path, anc, {"classify": True, "abundance": True,
                                  "assembly": True, "damage": True})


# --- absent sheet still defers --------------------------------------------- #
def test_absent_tsv_defers_group_check(tmp_path):
    missing = str(tmp_path / "not_created_yet.tsv")
    # no group error even though differential is on — the sheet doesn't exist yet
    assert c.build_config(project="t", samples=missing, db=_DB,
                          modules={"abundance": True, "differential": True})
