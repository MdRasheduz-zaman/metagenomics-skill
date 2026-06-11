"""Unit tests for ancient-DNA damage authentication (`workflow/scripts/damage_authenticate.py`).

Pure-Python: parses mapDamage2 terminal-deamination frequency tables and emits an
authenticity verdict. Fixtures mirror real 5pCtoT_freq.txt / 3pGtoA_freq.txt format.
"""
import importlib.util
import json
import pathlib

_SPEC = importlib.util.spec_from_file_location(
    "damage_authenticate",
    pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "damage_authenticate.py",
)
da = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(da)


def test_first_position_freq_reads_terminal_row(tmp_path):
    f = tmp_path / "5pCtoT_freq.txt"
    f.write_text("pos\t5pC>T\n1\t0.2800\n2\t0.1600\n3\t0.1000\n")
    assert da.first_position_freq(str(f)) == 0.28


def test_authenticate_requires_both_ends():
    # both ends elevated -> authentic
    assert da.authenticate(0.28, 0.24, 0.05)["damage_present"] is True
    # only one end elevated -> not authentic (conservative, double-stranded signature needs both)
    assert da.authenticate(0.28, 0.01, 0.05)["damage_present"] is False
    assert da.authenticate(0.01, 0.24, 0.05)["damage_present"] is False
    # neither -> modern/low-coverage
    res = da.authenticate(0.015, 0.012, 0.05)
    assert res["damage_present"] is False
    assert "modern contamination" in res["verdict"]


def test_run_writes_json(tmp_path):
    ct = tmp_path / "ct.txt"; ct.write_text("pos\t5pC>T\n1\t0.30\n2\t0.10\n")
    ga = tmp_path / "ga.txt"; ga.write_text("pos\t3pG>A\n1\t0.26\n2\t0.09\n")
    out = tmp_path / "auth.json"
    res = da.run(str(ct), str(ga), "anc1", threshold=0.05, out_json=str(out))
    assert res["sample"] == "anc1" and res["damage_present"] is True
    on_disk = json.load(open(out))
    assert on_disk["ct_5prime_pos1"] == 0.30 and on_disk["ga_3prime_pos1"] == 0.26
