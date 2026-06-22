"""The db.build pipeline-step foundation: the kraken2-build/bracken-build registries, the
config_builder db.build validation (strategy/taxonomy/library/read-length invariants), and
the doctor air-gap + build-tool advisories. Execution (actually building a DB) is exercised
by the e2e tests; this layer is tool-free.
"""
import os

import pytest

from metagx import config_builder as cb
from metagx import dbbuild, doctor, registry

SAMPLES = [{"sample": "s1", "r1": "x.fasta", "platform": "ont", "layout": "se"},
           {"sample": "s2", "r1": "y.fastq", "platform": "illumina", "layout": "se"}]


# --- registries --------------------------------------------------------------------------
def test_build_registries_exist_and_render():
    for t in ("kraken2-build", "bracken-build"):
        assert t in registry.list_tools()
    args = registry.render_args(
        "kraken2-build",
        {"strategy": "standard", "taxonomy": "real", "libraries": "viral", "no_masking": True},
        managed={"db": "DB", "threads": 4})
    assert "--db" in args and "--no-masking" in args
    # interpreted orchestration params never leak as CLI tokens
    assert "standard" not in args and "viral" not in args and "real" not in args


def test_bracken_build_read_length_is_managed():
    b = registry.render_args("bracken-build", {"kmer_len": 35},
                             managed={"db": "DB", "threads": 4, "read_length": 1000})
    assert b.count("-l") == 1 and "1000" in b
    with pytest.raises(registry.ValidationError):     # -l is managed, user may not set it
        registry.validate("bracken-build", {"read_length": 1000})


# --- config_builder db.build -------------------------------------------------------------
def test_standard_build_derives_read_lengths_from_platforms():
    cfg = cb.build_config(project="p", samples=SAMPLES,
                          db={"build": {"strategy": "standard", "libraries": "viral,bacteria"}},
                          modules={"classify": True, "abundance": True})
    b = cfg["db"]["build"]
    assert b["taxonomy"] == "real"                    # default when unsure
    assert b["read_lengths"] == [150, 1000]           # illumina + ont
    assert b["bracken_kmer_len"] == b.get("kmer_len", 35)
    assert cfg["db"]["kraken2"].endswith("dbs/standard")  # default output path


def test_explicit_read_lengths_respected():
    cfg = cb.build_config(project="p", samples=SAMPLES,
                          db={"kraken2": "d", "build": {"strategy": "standard",
                              "libraries": "viral", "read_lengths": [100, 250, 100]}},
                          modules={"classify": True})
    assert cfg["db"]["build"]["read_lengths"] == [100, 250]   # deduped + sorted


@pytest.mark.parametrize("bad,msg", [
    ({"strategy": "spike-in", "taxonomy": "synthetic", "libraries": "bacteria", "source": "g.fa"},
     "requires taxonomy: real"),
    ({"strategy": "custom-folder"}, "needs db.build.source"),
    ({"strategy": "standard", "libraries": ""}, "needs db.build.libraries"),
    ({"strategy": "standard", "libraries": "viral,bogus"}, "unknown NCBI libraries"),
    ({"strategy": "standard", "libraries": "viral", "kmer_len": 20, "minimizer_len": 31},
     "minimizer_len must be <= kmer_len"),
])
def test_db_build_invariants_reject(bad, msg):
    with pytest.raises(registry.ValidationError) as e:
        cb.build_config(project="p", samples=SAMPLES, modules={"classify": True}, db={"build": bad})
    assert msg in str(e.value)


def test_synthetic_custom_build_ok():
    cfg = cb.build_config(project="p", samples=SAMPLES,
                          db={"kraken2": "d", "build": {"strategy": "custom-fasta",
                              "taxonomy": "synthetic", "source": "genomes.fasta"}},
                          modules={"classify": True})
    assert cfg["db"]["build"]["taxonomy"] == "synthetic"
    assert cfg["db"]["build"]["source"] == "genomes.fasta"


# --- doctor advisories -------------------------------------------------------------------
def test_doctor_surfaces_airgap_for_download_builds():
    checks = doctor.check_db_build({"kraken2": "x", "build": {
        "strategy": "standard", "taxonomy": "real", "libraries": "viral"}})
    net = [c for c in checks if c.name == "db-build:network"]
    assert net and "air-gap" in net[0].message
    assert any(c.name.startswith("db-build:kraken2-build") for c in checks)


def test_doctor_no_airgap_for_synthetic_build():
    checks = doctor.check_db_build({"build": {
        "strategy": "custom-fasta", "taxonomy": "synthetic", "source": "g.fa"}})
    assert not any(c.name == "db-build:network" for c in checks)


def test_doctor_missing_db_is_info_when_build_configured():
    c = doctor.check_database({"kraken2": "/no/such/db", "build": {"strategy": "standard"}})
    assert c.status == "info" and "will be built" in c.message


def test_bracken_length_backfilled_from_db_build():
    """db.build derives databaseLmers from platforms; the classify-time Bracken -r map must
    match so abundance asks for a length the DB actually built."""
    cfg = cb.build_config(project="p", samples=SAMPLES,
                          db={"kraken2": "d", "build": {"strategy": "standard", "libraries": "viral"}},
                          modules={"classify": True, "abundance": True})
    assert cfg["bracken_read_length_by_platform"] == {"ont": 1000, "illumina": 150}


# --- build_database execution planning (run=False: no tools needed) ----------------------
def test_build_database_standard_plan(tmp_path):
    r = dbbuild.build_database(db_dir=str(tmp_path), strategy="standard",
                               libraries="viral,archaea", read_lengths=[150, 1000], run=False)
    cmds = r["commands"]
    assert any("--download-taxonomy" in v for v in cmds.values())
    assert "download-library-viral" in cmds and "download-library-archaea" in cmds
    assert sum(k.startswith("bracken-build") for k in cmds) == 2


def test_build_database_custom_folder_synthetic_writes_taxonomy(tmp_path):
    folder = tmp_path / "genomes"
    folder.mkdir()
    (folder / "a.fasta").write_text(">a\nACGTACGTACGT\n")
    (folder / "b.fna").write_text(">b\nTTTTGGGGCCCC\n")
    db = tmp_path / "db"
    r = dbbuild.build_database(db_dir=str(db), strategy="custom-folder", taxonomy="synthetic",
                               source=str(folder), read_lengths=[150], run=False)
    assert r["n_sequences"] == 2
    assert os.path.isfile(db / "taxonomy" / "names.dmp")
    assert "--no-masking" in r["commands"]["add-to-library"]   # synthetic avoids dustmasker


def test_build_database_spike_in_orders_taxonomy_then_library_then_custom(tmp_path):
    src = tmp_path / "mine.fasta"
    src.write_text(">x|kraken:taxid|9606\nACGTACGTACGT\n")
    r = dbbuild.build_database(db_dir=str(tmp_path / "db"), strategy="spike-in", taxonomy="real",
                               libraries="viral", source=str(src), read_lengths=[150], run=False)
    keys = list(r["commands"])
    assert keys[0] == "download-taxonomy"
    assert any(k.startswith("download-library") for k in keys)
    assert any(k.startswith("add-to-library") for k in keys)
    # real build masks by default -> the add step carries no --no-masking
    assert "--no-masking" not in r["commands"][next(k for k in keys if k.startswith("add-to-library"))]


def test_db_is_built_requires_core_and_bracken(tmp_path):
    assert dbbuild.db_is_built(str(tmp_path), [150]) is False
    for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (tmp_path / f).write_text("x")
    assert dbbuild.db_is_built(str(tmp_path), []) is True       # core present, no bracken needed
    assert dbbuild.db_is_built(str(tmp_path), [150]) is False   # bracken distrib missing
    (tmp_path / "database150mers.kmer_distrib").write_text("x")
    assert dbbuild.db_is_built(str(tmp_path), [150]) is True


def test_use_ftp_flows_through_config():
    cfg = cb.build_config(project="p", samples=SAMPLES,
                          db={"kraken2": "d", "build": {"strategy": "standard", "libraries": "viral"}},
                          modules={"classify": True})
    assert cfg["db"]["build"]["use_ftp"] is True
    cfg2 = cb.build_config(project="p", samples=SAMPLES,
                           db={"kraken2": "d", "build": {"strategy": "standard", "libraries": "viral",
                                                         "use_ftp": False}},
                           modules={"classify": True})
    assert cfg2["db"]["build"]["use_ftp"] is False


def test_doctor_warns_slow_download_for_large_libraries():
    big = doctor.check_db_build({"build": {"strategy": "standard", "libraries": "bacteria,viral"}})
    assert any(c.name == "db-build:slow-download" for c in big)
    small = doctor.check_db_build({"build": {"strategy": "standard", "libraries": "viral"}})
    assert not any(c.name == "db-build:slow-download" for c in small)


def test_doctor_taxid_precheck_real_custom(tmp_path):
    untagged = tmp_path / "u.fasta"
    untagged.write_text(">seq1 plain header\nACGTACGT\n")
    checks = doctor.check_db_build({"build": {"strategy": "custom-fasta", "taxonomy": "real",
                                              "source": str(untagged)}})
    assert any(c.name == "db-build:taxids" for c in checks)
    tagged = tmp_path / "t.fasta"
    tagged.write_text(">acc|kraken:taxid|9606 human\nACGTACGT\n")
    ok = doctor.check_db_build({"build": {"strategy": "custom-fasta", "taxonomy": "real",
                                          "source": str(tagged)}})
    assert not any(c.name == "db-build:taxids" for c in ok)
    # synthetic taxonomy never needs taxid headers
    syn = doctor.check_db_build({"build": {"strategy": "custom-fasta", "taxonomy": "synthetic",
                                           "source": str(untagged)}})
    assert not any(c.name == "db-build:taxids" for c in syn)


def test_report_db_info_carries_build_manifest(tmp_path):
    import json as _json
    from metagx import report
    (tmp_path / "hash.k2d").write_bytes(b"x" * 16)
    (tmp_path / ".metagx_db.json").write_text(_json.dumps(
        {"strategy": "standard", "taxonomy": "real", "kraken2_version": "2.17.1"}))
    info = report.db_info(str(tmp_path))
    assert info["build"]["strategy"] == "standard"
    assert info["build"]["kraken2_version"] == "2.17.1"


def test_db_build_auto_false_skips_autobuild_dependency():
    """auto: false must drop the build manifest from classify's inputs (no surprise build)."""
    cfg_auto = cb.build_config(project="p", samples=SAMPLES,
                               db={"kraken2": "d", "build": {"strategy": "standard", "libraries": "viral"}},
                               modules={"classify": True})
    cfg_manual = cb.build_config(project="p", samples=SAMPLES,
                                 db={"kraken2": "d", "build": {"strategy": "standard", "libraries": "viral",
                                                               "auto": False}},
                                 modules={"classify": True})
    assert cfg_auto["db"]["build"]["auto"] is True
    assert cfg_manual["db"]["build"]["auto"] is False


def test_doctor_warns_blast_cost_for_big_standard_library():
    # db.build.blast on a big standard library => makeblastdb is heavy; doctor should warn.
    checks = doctor.check_db_build({"kraken2": "x", "build": {
        "strategy": "standard", "taxonomy": "real", "libraries": "bacteria", "blast": True}})
    assert any(c.name == "db-build:blast-cost" and c.status == "warn" for c in checks)
    # ...but NOT for a small viral library
    small = doctor.check_db_build({"kraken2": "x", "build": {
        "strategy": "standard", "taxonomy": "real", "libraries": "viral", "blast": True}})
    assert not any(c.name == "db-build:blast-cost" for c in small)
