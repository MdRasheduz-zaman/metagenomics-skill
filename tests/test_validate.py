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


def test_normalize_seqid2taxid_strips_kraken_cruft(tmp_path):
    m = tmp_path / "seqid2taxid.map"
    m.write_text("NC_001477.1|kraken:taxid|1001\t1001\nNC_012532.1|kraken:taxid|1002\t1002\n")
    got = validate.normalize_seqid2taxid(str(m))
    assert got == {"NC_001477.1": "1001", "NC_012532.1": "1002"}  # bare accession keys


def test_parse_names_dmp_uses_name_txt_not_unique_name(tmp_path):
    # regression: name is column 2 (name_txt), not column 3 (unique_name, usually empty)
    n = tmp_path / "names.dmp"
    n.write_text("1004\t|\tYellow fever virus\t|\t\t|\tscientific name\t|\n"
                 "1004\t|\tYFV\t|\t\t|\tsynonym\t|\n")
    names = validate.parse_names_dmp(str(n))
    assert names["1004"] == "Yellow fever virus"  # scientific name only, real text


def test_kraken2_db_sources_finds_custom_library(tmp_path):
    d = tmp_path / "k2"
    (d / "taxonomy").mkdir(parents=True)
    (d / "custom_library.fasta").write_text(">NC_1|kraken:taxid|1001 x\nACGT\n")
    (d / "seqid2taxid.map").write_text("NC_1|kraken:taxid|1001\t1001\n")
    (d / "taxonomy" / "names.dmp").write_text("1001\t|\tX\t|\t\t|\tscientific name\t|\n")
    src = validate.kraken2_db_sources(str(d))
    assert src["fastas"] and src["fastas"][0].endswith("custom_library.fasta")
    assert src["seqid2taxid"] and src["names_dmp"]


def test_kraken2_db_sources_raises_for_prebuilt_or_cleaned(tmp_path):
    # a prebuilt/fetched (or --clean'd) index has only the opaque *.k2d hash, no genomes
    d = tmp_path / "prebuilt"
    d.mkdir()
    for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (d / f).write_bytes(b"\x00\x01")
    with pytest.raises(ValueError) as e:
        validate.kraken2_db_sources(str(d))
    assert "no source genomes" in str(e.value)


def test_assess_uses_name_resolver_for_in_sync_taxonomy():
    # BLAST hit carries kraken2's taxid in staxids; resolver (names.dmp) gives the in-sync name
    rows = validate.parse_blast6("r1\tNC_1\t99.0\t150\t1e-50\t280\t1004\tN/A\tsome title\n")
    best = validate.best_hit_per_query(rows)
    resolver = {"1004": "Yellow fever virus"}.get
    a = validate.assess({"r1": "Yellow fever virus"}, best, level="genus", name_resolver=resolver)
    assert a["n_agree"] == 1 and a["per_query"][0]["blast"] == "Yellow fever virus"


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


def test_validate_build_from_fasta_satisfies_db_requirement():
    # building the BLAST DB from a FASTA (the in-scope design) means no db.blast is needed up front
    cfg = cb.build_config(samples=SAMPLES, db={"kraken2": "DB"},
                          modules={"validate": True}, validate={"build_from": "refs.fasta"})
    assert cfg["validate"]["build_from"] == "refs.fasta"


def test_validate_build_from_classifier_requires_kraken2_db_path():
    # "classifier" resolves the BLAST DB from the kraken2 DB dir, so a config with no kraken2
    # DB at all is rejected (the classify gate fires first; either way it points at db.kraken2).
    with pytest.raises(registry.ValidationError) as e:
        cb.build_config(samples=SAMPLES, db={},
                        modules={"validate": True, "classify": True},
                        validate={"build_from": "classifier"})
    assert "db.kraken2" in str(e.value)


def test_validate_build_from_classifier_ok_with_kraken2_db():
    cfg = cb.build_config(samples=SAMPLES, db={"kraken2": "/my/k2db"},
                          modules={"validate": True}, validate={"build_from": "classifier"})
    assert cfg["validate"]["build_from"] == "classifier"


def test_validate_build_from_classifier_ok_via_db_build():
    # a db.build defaults db.kraken2 to the build output dir, satisfying "classifier"
    cfg = cb.build_config(
        samples=SAMPLES,
        db={"build": {"strategy": "custom-fasta", "source": "genomes.fasta", "taxonomy": "synthetic"}},
        modules={"validate": True}, validate={"build_from": "classifier"})
    assert cfg["validate"]["build_from"] == "classifier"
    assert cfg["db"]["kraken2"]  # defaulted from db.build


def test_db_build_blast_defaults_on_with_validate_and_colocates():
    # validate + db.build => build the aligned BLAST DB together, db.blast co-located, no build_from
    cfg = cb.build_config(
        samples=SAMPLES,
        db={"build": {"strategy": "custom-fasta", "source": "g.fasta", "taxonomy": "synthetic"}},
        modules={"validate": True})
    assert cfg["db"]["build"]["blast"] is True
    assert cfg["db"]["blast"].endswith("blast/insync")
    assert "build_from" not in cfg["validate"]  # uses the co-located db.blast directly


def test_db_build_blast_explicit_false_does_not_colocate():
    # opt-out: user supplies their own (out-of-scope) db.blast; we don't co-locate
    cfg = cb.build_config(
        samples=SAMPLES,
        db={"blast": "/ext/nt",
            "build": {"strategy": "custom-fasta", "source": "g.fasta", "taxonomy": "synthetic",
                      "blast": False}},
        modules={"validate": True})
    assert cfg["db"]["build"]["blast"] is False
    assert cfg["db"]["blast"] == "/ext/nt"


def test_db_build_blast_opt_in_without_validate():
    # "I want both DBs" even though validate isn't enabled in this config
    cfg = cb.build_config(
        samples=SAMPLES,
        db={"build": {"strategy": "custom-fasta", "source": "g.fasta", "taxonomy": "synthetic",
                      "blast": True}},
        modules={"classify": True})
    assert cfg["db"]["build"]["blast"] is True
    assert cfg["db"]["blast"].endswith("blast/insync")


def test_doctor_strong_warns_validate_without_aligned_blast():
    from metagx import doctor
    cfg = {"modules": {"validate": True, "classify": True},
           "db": {"blast": "/ext/nt", "build": {"strategy": "custom-fasta", "blast": False}},
           "validate": {}}
    checks = doctor.check_validate_alignment(cfg)
    assert checks and checks[0].status == "warn" and "STRONG" in checks[0].message


def test_doctor_ok_when_kraken2_and_blast_built_together():
    from metagx import doctor
    cfg = {"modules": {"validate": True},
           "db": {"build": {"strategy": "custom-fasta", "blast": True}}, "validate": {}}
    checks = doctor.check_validate_alignment(cfg)
    assert checks and checks[0].status == "ok"


@pytest.mark.skipif(not (shutil.which("kraken2-build") and shutil.which("makeblastdb")
                         and shutil.which("blastdbcmd") and os.path.isdir(_FIXTURE)),
                    reason="needs kraken2-build + BLAST+ and the viral fixture")
def test_joint_build_produces_aligned_taxid_tagged_blast_db(tmp_path):
    """db.build with build_blast builds kraken2 + an aligned, taxid-tagged BLAST DB together."""
    from metagx import dbbuild
    db_dir = str(tmp_path / "db")
    res = dbbuild.build_database(db_dir=db_dir, strategy="custom-fasta", taxonomy="synthetic",
                                 source=os.path.join(_FIXTURE, "genomes.fasta"),
                                 read_lengths=[150], threads=2, build_blast=True, run=True)
    assert res.get("ok") and res["blast"]["ok"] and res["blast"].get("taxid_mapped")
    prefix = os.path.join(db_dir, "blast", "insync")
    assert validate.blast_db_present(prefix)
    out = subprocess.run(["blastdbcmd", "-db", prefix, "-entry", "all", "-outfmt", "%T"],
                         capture_output=True, text=True).stdout.split()
    assert "1001" in out  # subjects carry kraken2's exact (synthetic) taxids


def test_build_from_drops_blast_from_needed_dbs():
    from metagx import dbprovision
    cfg = {"modules": {"classify": True, "validate": True}, "validate": {"build_from": "x.fasta"}}
    assert "blast" not in dbprovision.needed_dbs(cfg)
    cfg_need = {"modules": {"classify": True, "validate": True}, "validate": {}}
    assert "blast" in dbprovision.needed_dbs(cfg_need)


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


@pytest.mark.skipif(not (shutil.which("makeblastdb") and os.path.isdir(_FIXTURE)),
                    reason="needs BLAST+ makeblastdb and the viral fixture")
def test_build_blast_db_from_fixture_is_in_scope_and_idempotent(tmp_path):
    """build_blast_db builds an in-scope DB from the classifier's own genomes and is idempotent."""
    prefix = str(tmp_path / "insync")
    res = validate.build_blast_db(os.path.join(_FIXTURE, "genomes.fasta"), prefix)
    assert res["ok"] and res.get("ran") is True
    assert validate.blast_db_present(prefix)
    again = validate.build_blast_db(os.path.join(_FIXTURE, "genomes.fasta"), prefix)
    assert again["ok"] and again.get("ran") is False  # skipped — already present


@pytest.mark.skipif(not (shutil.which("makeblastdb") and shutil.which("blastdbcmd")
                         and os.path.isdir(_FIXTURE)),
                    reason="needs BLAST+ (makeblastdb/blastdbcmd) and the viral fixture")
def test_taxid_map_attaches_taxids_to_blast_subjects(tmp_path):
    """A normalized acc->taxid map makes the BLAST subjects carry those exact taxids — the
    bridge that lets validation compare taxid-to-taxid instead of fuzzy strings."""
    fasta = tmp_path / "g.fasta"
    fasta.write_text(">NC_001477.1 Dengue virus 1\nACGTACGTACGTACGTACGTACGTACGT\n"
                     ">NC_002031.1 Yellow fever virus\nTTTTGGGGCCCCAAAATTTTGGGGCCCC\n")
    tmap = tmp_path / "acc2taxid.tsv"
    tmap.write_text("NC_001477.1\t1001\nNC_002031.1\t1004\n")
    prefix = str(tmp_path / "db")
    res = validate.build_blast_db(str(fasta), prefix, taxid_map=str(tmap))
    assert res["ok"] and res.get("taxid_mapped")
    out = subprocess.run(["blastdbcmd", "-db", prefix, "-entry", "all", "-outfmt", "%a %T"],
                         capture_output=True, text=True).stdout
    pairs = dict(line.split()[:2] for line in out.splitlines() if line.split())
    # accession -> taxid carried through (note: -parse_seqids may keep version in %a)
    assert any(t == "1001" for t in pairs.values()) and any(t == "1004" for t in pairs.values())
