import json
import os

import pytest

from metagx import consent, probe, registry


def _fastq(path, reads, qual_char):
    with open(path, "w") as fh:
        for i, seq in enumerate(reads):
            fh.write(f"@r{i}\n{seq}\n+\n{qual_char * len(seq)}\n")


def _chdir_tmp(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # isolate .metagx/consent.json per test


@pytest.fixture
def short_clean(tmp_path):
    p = tmp_path / "ill.fastq"
    _fastq(p, ["ACGT" * 38] * 200, "I")  # 152 bp, Q40 -> illumina, low error
    return str(p)


@pytest.fixture
def long_noisy(tmp_path):
    p = tmp_path / "ont.fastq"
    _fastq(p, ["ACGTACGTAC" * 100] * 50, "+")  # 1000 bp, Q10 -> ont, high error
    return str(p)


def test_profile_short_clean(short_clean):
    pr = probe.profile_file(short_clean)
    assert pr["read_length"]["median"] == 152
    assert pr["inferred_platform_class"] == "illumina"
    assert pr["mean_q"] == 40.0 and pr["q20_frac"] == 1.0
    assert pr["est_error"] < 0.01


def test_profile_long_noisy(long_noisy):
    pr = probe.profile_file(long_noisy)
    assert pr["read_length"]["median"] == 1000
    assert pr["inferred_platform_class"] == "ont"
    assert pr["q20_frac"] == 0.0  # Q10 -> nothing passes Q20


def test_run_profiles_all_samples(tmp_path, monkeypatch, short_clean, long_noisy):
    _chdir_tmp(tmp_path, monkeypatch)
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tr1\tplatform\nA\t%s\tillumina\nB\t%s\tont\n" % (short_clean, long_noisy))
    rep = probe.run(str(sheet), assume_yes=True)
    assert rep["measured"] and set(rep["samples"]) == {"A", "B"}
    assert rep["project"]["platform_consensus"] == "mixed"


def test_run_flags_sheet_mismatch(tmp_path, monkeypatch, short_clean):
    _chdir_tmp(tmp_path, monkeypatch)
    # declare ONT but the reads are short+clean -> short/long mismatch warning
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tr1\tplatform\nX\t%s\tont\n" % short_clean)
    rep = probe.run(str(sheet), assume_yes=True)
    assert any("mismatch" in w for w in rep["project"]["warnings"])
    assert rep["context"]["platform_mismatch"] is True


def test_context_drives_promotion(tmp_path, monkeypatch):
    _chdir_tmp(tmp_path, monkeypatch)
    deep = tmp_path / "deep.fastq"
    _fastq(deep, ["ACGTACGTAC" * 100] * 200, "I")  # long + lots of reads
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tr1\tplatform\nD\t%s\tont\n" % deep)
    ctx = probe.run(str(sheet), assume_yes=True, max_reads=50)["context"]
    # force the deep threshold regardless of fixture size, to prove the wiring end-to-end
    ctx["estimated_bases"] = 6e10
    names = {p["name"] for p in registry.interview_spec("flye", max_tier=2, context=ctx)}
    assert "asm_coverage" in names


def test_consent_off_degrades(tmp_path, monkeypatch):
    _chdir_tmp(tmp_path, monkeypatch)
    consent.set("probe", "off")
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tr1\tplatform\nA\tnope.fastq\tillumina\n")
    rep = probe.run(str(sheet))  # no assume_yes, stored 'off'
    assert rep["measured"] is False and rep["context"] == {"measured": False}


def test_output_is_non_reconstructive(tmp_path, monkeypatch):
    _chdir_tmp(tmp_path, monkeypatch)
    secret = "GATTACAGATTACAGATTACAGATTACAGATTACA"  # a distinctive read sequence
    fq = tmp_path / "x.fastq"
    _fastq(fq, [secret] * 20, "I")
    sheet = tmp_path / "s.tsv"
    sheet.write_text("sample\tr1\tplatform\nA\t%s\tillumina\n" % fq)
    out = tmp_path / "probe_out"
    probe.run(str(sheet), assume_yes=True, out=str(out))
    blob = (out / "probe.json").read_text() + (out / "probe.md").read_text()
    assert secret not in blob and "@r0" not in blob  # no sequence, no read ID


def test_thresholds_come_from_evidence():
    from metagx import evidence_pack
    ev = evidence_pack.load_evidence("platform_inference")
    th = probe._load_thresholds()
    assert th["long_min_len"] == ev["long_min_median_len"]
    assert th["hifi_max_err"] == ev["accurate_max_est_error"]
    assert th["low_q20"] == ev["low_q20_fraction"]


def test_bounded_by_max_reads(tmp_path):
    fq = tmp_path / "big.fastq"
    _fastq(fq, ["ACGT" * 38] * 1000, "I")
    pr = probe.profile_file(str(fq), max_reads=100)
    assert pr["n_sampled"] == 100  # stopped at the cap, did not read all 1000
