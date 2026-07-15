"""The generated sample sheet's STRUCTURE is derived from the analysis (deterministic, no data):
the tool provides the right columns/shape from the interview answers, it doesn't just validate a
sheet the user brings. differential -> group column + 2 groups x 2; decontam -> control row;
aDNA -> library=ancient.
"""
from metagx import project


def _header(t):
    return t.splitlines()[0].split("\t")


def _body(t):
    return t.splitlines()[1:]


def test_plain_run_is_generic_two_rows():
    t = project._samples_template("ont", {"modules": {"classify": True}})
    assert _header(t) == ["sample", "r1", "r2", "platform", "layout", "library"]
    assert len(_body(t)) == 2


def test_differential_has_group_column_and_two_by_two():
    t = project._samples_template("illumina", {"modules": {"differential": True}})
    assert "group" in _header(t)
    body = _body(t)
    assert len(body) == 4
    groups = [r.split("\t")[-1] for r in body]
    assert groups.count("case") == 2 and groups.count("control") == 2


def test_differential_respects_custom_group_column():
    t = project._samples_template(
        "illumina", {"modules": {"differential": True}, "differential": {"group_column": "treatment"}})
    assert "treatment" in _header(t) and "group" not in _header(t)


def test_decontam_adds_control_row():
    t = project._samples_template("illumina", {"modules": {"decontam": True}})
    assert "control" in _header(t)
    assert any(r.split("\t")[-1] == "true" for r in _body(t))   # a BLANK negative control


def test_adna_marks_library_ancient():
    t = project._samples_template("illumina", {"modules": {"damage": True, "assembly": True}})
    assert all(r.split("\t")[5] == "ancient" for r in _body(t))


def test_short_read_shows_paired_example():
    t = project._samples_template("illumina", {"modules": {"classify": True}})
    assert any("_R2.fastq.gz" in r and "\tpe\t" in r for r in _body(t))


def test_long_read_is_single_end():
    t = project._samples_template("ont", {"modules": {"classify": True}})
    assert all("\tse\t" in r and "_R2" not in r for r in _body(t))
