"""Real end-to-end execution of the workflow against the bundled viral database.

This is the test the project lacked: it runs the *actual* Snakemake pipeline —
kraken2 classification + Bracken abundance — on real reads for each sequencing
platform (ONT / Illumina / PacBio) against the 30-genome custom viral database in
``local_databases/viral_custom``, and asserts the answers are *correct*, not merely
that files appeared. It is the difference between "the code imports" and "the
pipeline produces the right taxa."

It skips cleanly when the bio tools, the database, or the data are absent — so CI
(which has none of them) stays green — and runs fully on a machine where the
``metagx-bio`` conda env is on PATH and ``data/`` + ``local_databases/`` are present:

    export PATH="$HOME/miniconda3/envs/metagx-bio/bin:$PATH"
    pytest tests/test_pipeline_e2e.py -q

The database was built from ``data/genomes.fasta`` (Dengue, Zika, Chikungunya,
West Nile, Yellow fever, ... 30 viral genomes), so classified reads must land on
those species and nothing outside them.
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest
import yaml

from metagx import dbbuild, runner

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
# Committed fixtures so the database can be (re)built from scratch in CI, where the large
# prebuilt DB and data/ are gitignored and absent.
FIXTURE_GENOMES = REPO / "tests" / "fixtures" / "viral" / "genomes.fasta"
FIXTURE_ONT = REPO / "tests" / "fixtures" / "viral" / "ont_reads.fasta"
PREBUILT_DB = REPO / "local_databases" / "viral_custom"

# Tools the classify+abundance path needs.
_HAVE_TOOLS = bool(shutil.which("kraken2") and shutil.which("bracken"))
# The DB is available if it's already built, OR we can build it from the committed genomes
# fixture (kraken2-build + bracken-build present). This is what lets CI run the real pipeline.
# Set METAGX_FORCE_DB_BUILD=1 to ignore any prebuilt DB and always build from the committed
# fixture — gives the from-scratch build path coverage even on a dev box that *has* a prebuilt
# DB (the local state that hid the SIGPIPE/thread bugs for four CI rounds). CI builds from
# scratch anyway; this lets a developer reproduce the CI build path locally without deleting
# their DB.
_FORCE_BUILD = bool(os.environ.get("METAGX_FORCE_DB_BUILD"))
_PREBUILT_OK = (not _FORCE_BUILD) and (PREBUILT_DB / "hash.k2d").is_file() and \
    (PREBUILT_DB / "database150mers.kmer_distrib").is_file()
_CAN_BUILD = bool(shutil.which("kraken2-build") and shutil.which("bracken-build")
                  and FIXTURE_GENOMES.is_file())

requires_stack = pytest.mark.skipif(
    not (_HAVE_TOOLS and (_PREBUILT_OK or _CAN_BUILD)),
    reason="kraken2/bracken not on PATH, and no prebuilt viral DB nor kraken2-build to build "
    "one from the genomes fixture (skips in CI without bio tools)",
)


@pytest.fixture(scope="session")
def viral_db(tmp_path_factory):
    """Path to the viral kraken2+Bracken DB: the prebuilt one if present, else built fresh
    from the 30-genome fixture (with 150- and 1000-mer Bracken distributions). Built once
    per session. This is the same `metagx build-db` path users run, so building it here also
    exercises the custom-DB build recipe end to end."""
    if _PREBUILT_OK:
        return PREBUILT_DB
    db_dir = tmp_path_factory.mktemp("viral_db")
    res = dbbuild.build_db(str(FIXTURE_GENOMES), str(db_dir),
                           read_length=[150, 1000], threads=4, run=True)
    if not res.get("ok"):
        # Surface the actual tool output untruncated — a bare `{res}` gets abbreviated by
        # pytest's repr (the `...`), which has hidden the real cause across CI rounds.
        step = res.get("failed_step")
        log = (res.get("logs") or {}).get(step, {})
        retry = log.get("retry_threads1")
        raise AssertionError(
            "DB build failed.\n"
            f"  failed_step : {step}\n"
            f"  note        : {res.get('note')}\n"
            f"  returncode  : {log.get('returncode')}\n"
            f"  tail        : {log.get('tail')}\n"
            f"  retry(-t1)  : {retry}"
        )
    return db_dir

# The assembly→map→bin→reconcile path needs the long-read assembler + mapping stack.
_ASM_TOOLS = ("flye", "minimap2", "samtools", "metabat2",
              "jgi_summarize_bam_contig_depths")
_HAVE_ASM = all(shutil.which(t) for t in _ASM_TOOLS)
requires_assembly = pytest.mark.skipif(
    not (_HAVE_TOOLS and (_PREBUILT_OK or _CAN_BUILD) and _HAVE_ASM),
    reason="assembly/binning stack (flye/minimap2/samtools/metabat2/jgi) not all present "
    "(skips in CI; run with the metagx-bio env active)",
)

# A handful of species that are definitely in the 30-genome DB; classified reads
# from any platform should recover at least some of them.
KNOWN_SPECIES = ("Dengue", "Zika", "Chikungunya", "West Nile", "Yellow fever")


# --------------------------------------------------------------------------- #
# Parsers (the report formats are stable kraken2/Bracken text)                #
# --------------------------------------------------------------------------- #
def _parse_kreport(path: Path):
    """Return (classified_fraction, {species_name: clade_reads})."""
    unclassified = 0
    classified = 0
    species = {}
    for line in path.read_text().splitlines():
        cols = line.split("\t")
        if len(cols) < 6:
            continue
        clade_reads, rank, name = int(cols[1]), cols[3].strip(), cols[5].strip()
        if rank == "U":
            unclassified += clade_reads
        if rank == "R":  # root clade = everything classified
            classified += clade_reads
        if rank == "S":
            species[name] = clade_reads
    total = unclassified + classified
    frac = classified / total if total else 0.0
    return frac, species


def _parse_bracken(path: Path):
    """Return list of (name, fraction) rows from a Bracken output table."""
    rows = []
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    fcol = header.index("fraction_total_reads")
    for line in lines[1:]:
        c = line.split("\t")
        if len(c) > fcol:
            rows.append((c[0], float(c[fcol])))
    return rows


# --------------------------------------------------------------------------- #
# Fixtures: tiny per-platform inputs derived from the bundled data            #
# --------------------------------------------------------------------------- #
def _subsample_fasta(src: Path, dst: Path, n_reads: int) -> Path:
    """Copy the first ``n_reads`` FASTA records (keeps the e2e run fast)."""
    out, count = [], 0
    for line in src.read_text().splitlines():
        if line.startswith(">"):
            count += 1
            if count > n_reads:
                break
        out.append(line)
    dst.write_text("\n".join(out) + "\n")
    return dst


# Per-platform spec: where the reads live, how to declare the sample, the Bracken
# read length, and the correctness band expected for *this* data.
#   - Illumina/PacBio reads were simulated straight from the DB genomes -> high recall.
#   - ONT is a noisier real metagenomic FASTA -> partial classification.
PLATFORM_SPECS = {
    "ont": dict(
        reads=[DATA / "simulated_metagenomic_reads.fasta"],
        subsample=4000, platform="ont", layout="se",
        bracken_len=1000, frac_min=0.20, frac_max=0.85, min_species=10,
    ),
    "illumina": dict(
        reads=[DATA / "illumina_sim" / "illumina_sim_R1.fastq.gz",
               DATA / "illumina_sim" / "illumina_sim_R2.fastq.gz"],
        subsample=None, platform="illumina", layout="pe",
        bracken_len=150, frac_min=0.80, frac_max=1.01, min_species=20,
    ),
    "pacbio_hifi": dict(
        reads=[DATA / "pacbio_sim" / "pacbio_hifi.fastq.gz"],
        subsample=None, platform="pacbio_hifi", layout="se",
        bracken_len=1000, frac_min=0.30, frac_max=1.01, min_species=10,
    ),
}


def _build_config(tmp_path: Path, project: str, row: dict, bracken_len: int, db: Path) -> Path:
    sheet = tmp_path / f"{project}.samples.tsv"
    cols = ["sample", "r1", "r2", "platform", "layout"]
    row = {**{"sample": "s", "r2": ""}, **row}
    sheet.write_text(
        "\t".join(cols) + "\n" + "\t".join(str(row.get(c, "")) for c in cols) + "\n"
    )
    cfg = {
        "project": project,
        "outdir": str(tmp_path / "out"),
        "threads": 4,
        "samples": str(sheet),
        "db": {"kraken2": str(db), "bracken": str(db)},
        "modules": {"qc": row["platform"] == "illumina", "classify": True,
                    "abundance": True, "assembly": False},
        "kraken2": {"minimum_hit_groups": 2},
        "bracken": {"read_length": bracken_len, "level": "S", "threshold": 0},
    }
    cfg_path = tmp_path / f"{project}.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    return cfg_path


# --------------------------------------------------------------------------- #
# The test                                                                    #
# --------------------------------------------------------------------------- #
@requires_stack
@pytest.mark.parametrize("platform", sorted(PLATFORM_SPECS))
def test_platform_classifies_correctly(platform, tmp_path, viral_db):
    spec = PLATFORM_SPECS[platform]
    reads = spec["reads"]
    # ONT can fall back to the committed reads fixture so it runs in CI (data/ is gitignored).
    if platform == "ont" and not reads[0].is_file() and FIXTURE_ONT.is_file():
        reads = [FIXTURE_ONT]
    if not all(p.is_file() for p in reads):
        pytest.skip(f"{platform}: reads not present (no data/ and no fixture)")

    # Stage inputs (subsampling the large ONT FASTA to keep the run fast).
    if spec["subsample"]:
        r1 = _subsample_fasta(reads[0], tmp_path / "in.fasta", spec["subsample"])
        row = {"r1": str(r1), "platform": spec["platform"], "layout": spec["layout"]}
    else:
        row = {"r1": str(reads[0]), "platform": spec["platform"],
               "layout": spec["layout"]}
        if len(reads) > 1:
            row["r2"] = str(reads[1])

    project = f"e2e_{platform}"
    cfg = _build_config(tmp_path, project, row, spec["bracken_len"], viral_db)

    proc = runner.run(config=str(cfg), cores=4, dry_run=False)
    assert proc.returncode == 0, (
        f"[{platform}] pipeline failed:\n{(proc.stderr or proc.stdout)[-3000:]}"
    )

    out = tmp_path / "out" / project
    kreport = out / "kraken2" / "s.confidence_0.0.kreport"
    bracken = out / "bracken" / "s.confidence_0.0.S.bracken"
    assert kreport.is_file(), f"[{platform}] no kreport produced"
    assert bracken.is_file(), f"[{platform}] no bracken output produced"

    frac, species = _parse_kreport(kreport)
    # 1) classified fraction sits in the platform's expected band
    assert spec["frac_min"] <= frac <= spec["frac_max"], (
        f"[{platform}] classified fraction {frac:.3f} outside "
        f"[{spec['frac_min']}, {spec['frac_max']}]"
    )
    # 2) enough distinct species recovered (DB has 30)
    assert len(species) >= spec["min_species"], (
        f"[{platform}] only {len(species)} species recovered "
        f"(expected >= {spec['min_species']})"
    )
    # 3) the hits are real DB members, not noise: at least one known virus
    assert any(any(k in name for k in KNOWN_SPECIES) for name in species), (
        f"[{platform}] no known DB virus among {list(species)[:5]}"
    )
    # 4) Bracken redistributed into a valid distribution summing to ~1
    rows = _parse_bracken(bracken)
    assert rows, f"[{platform}] empty Bracken table"
    total = sum(f for _, f in rows)
    assert 0.98 <= total <= 1.02, f"[{platform}] Bracken fractions sum to {total:.3f}"


# --------------------------------------------------------------------------- #
# Deep path: assembly (Flye) + mapping + binning + read/contig reconciliation #
# --------------------------------------------------------------------------- #
import json
import re


@requires_assembly
def test_ont_assembly_binning_reconcile(tmp_path, viral_db):
    """Run the heavy modules end-to-end on ONT reads and check they cohere.

    Exercises Flye assembly, minimap2/samtools mapping, jgi depth, MetaBAT2 binning,
    kraken2-on-contigs, and the reconcile script — the render_args call sites that the
    classify-only test never touches. Correctness checks:
      * contigs are actually assembled,
      * read-level and contig-level kraken2 calls agree on most reads (concordance),
      * binning completes and reports a bin count (zero MAGs is fine on this toy data —
        what matters is the pipeline does not crash on it).
    """
    db = viral_db
    ont_fa = DATA / "simulated_metagenomic_reads.fasta"
    if not ont_fa.is_file():
        ont_fa = FIXTURE_ONT     # CI fallback (committed reads)
    if not ont_fa.is_file():
        pytest.skip("no ONT reads (data/ and fixture both absent)")
    reads = _subsample_fasta(ont_fa, tmp_path / "ont.fasta", 8000)

    sheet = tmp_path / "deep.samples.tsv"
    sheet.write_text(
        "sample\tr1\tr2\tplatform\tlayout\n"
        f"s\t{reads}\t\tont\tse\n"
    )
    cfg = {
        "project": "deep",
        "outdir": str(tmp_path / "out"),
        "threads": 4,
        "samples": str(sheet),
        "db": {"kraken2": str(db), "bracken": str(db)},
        "modules": {"qc": False, "classify": True, "abundance": True,
                    "assembly": True, "binning": True, "reconcile": True},
        "kraken2": {"minimum_hit_groups": 2},
        "bracken": {"read_length": 1000, "level": "S", "threshold": 0},
        "flye": {"meta": True, "min_overlap": 1000},
    }
    cfg_path = tmp_path / "deep.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=4, dry_run=False)
    assert proc.returncode == 0, f"deep pipeline failed:\n{(proc.stderr or proc.stdout)[-3000:]}"

    out = tmp_path / "out" / "deep"
    contigs = out / "assembly" / "s" / "final.contigs.fa"
    recon = out / "reconcile" / "s.read_accuracy.json"
    per_taxon = out / "reconcile" / "s.reconciliation.tsv"
    bins_marker = out / "binning" / "s" / "bins.done"

    # 1) assembly produced contigs
    assert contigs.is_file(), "no contigs assembled"
    n_contigs = contigs.read_text().count(">")
    assert n_contigs >= 5, f"only {n_contigs} contigs assembled"

    # 2) read calls agree with contig calls on the vast majority of aligned reads
    acc = json.loads(recon.read_text())
    rates = acc["rates_on_classified_contigs"]
    assert acc["reads_aligned"] > 0
    assert rates["concordant_pct"] >= 80.0, (
        f"read/contig concordance only {rates['concordant_pct']}%"
    )

    # 3) the per-taxon reconciliation recovers known DB viruses with assembly support
    taxa_text = per_taxon.read_text()
    assert any(k in taxa_text for k in KNOWN_SPECIES), "no known DB virus in reconciliation"

    # 4) binning completed without crashing and recorded a (possibly zero) bin count
    assert bins_marker.is_file(), "binning marker missing"
    m = re.search(r"bins=(\d+)", bins_marker.read_text())
    assert m, f"binning marker has no count: {bins_marker.read_text()!r}"
    assert int(m.group(1)) >= 0


# --------------------------------------------------------------------------- #
# Differential abundance + diversity on the case/control demo data            #
# --------------------------------------------------------------------------- #
@requires_stack
def test_differential_and_diversity_module(tmp_path, viral_db):
    """Stats path: per-sample classify+abundance -> α/β-diversity + differential abundance.

    The four demo samples are random subsamples of one population, so the statistically
    correct answer is *zero* significant taxa — this verifies the FDR control does not
    manufacture false positives, plus that the diversity metrics are well-formed.
    """
    db = viral_db
    diff_dir = DATA / "diff_demo"
    samples = {"case1": "case", "case2": "case", "ctrl1": "control", "ctrl2": "control"}
    if not all((diff_dir / f"{s}.fasta").is_file() for s in samples):
        pytest.skip("diff_demo data not present under data/")

    sheet = tmp_path / "diff.samples.tsv"
    lines = ["sample\tr1\tr2\tplatform\tlayout\tgroup"]
    for s, grp in samples.items():
        lines.append(f"{s}\t{diff_dir / f'{s}.fasta'}\t\tillumina\tse\t{grp}")
    sheet.write_text("\n".join(lines) + "\n")

    cfg = {
        "project": "diff",
        "outdir": str(tmp_path / "out"),
        "threads": 4,
        "samples": str(sheet),
        "db": {"kraken2": str(db), "bracken": str(db)},
        "modules": {"qc": False, "classify": True, "abundance": True,
                    "stats": True, "differential": True},
        "differential": {"n_permutations": 999, "fdr": 0.1, "reference_group": "control"},
        "kraken2": {"minimum_hit_groups": 2},
        "bracken": {"read_length": 150, "level": "S", "threshold": 0},
    }
    cfg_path = tmp_path / "diff.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=4, dry_run=False)
    assert proc.returncode == 0, f"differential pipeline failed:\n{(proc.stderr or proc.stdout)[-3000:]}"

    out = tmp_path / "out" / "diff"
    da = json.loads((out / "stats" / "differential_abundance.json").read_text())
    # 1) ran over all taxa across both groups, with the configured method
    assert da["summary"]["n_a"] == 2 and da["summary"]["n_b"] == 2
    assert da["summary"]["n_taxa"] > 0
    assert "CLR" in da["summary"]["method"]
    # 2) FDR control yields no false positives on these null (single-population) samples
    assert da["summary"]["n_significant"] == 0, (
        f"expected 0 significant taxa on null data, got {da['summary']['n_significant']}: "
        f"{da['significant_taxa']}"
    )

    # 3) α-diversity table is well-formed with sane Shannon entropies
    alpha = (out / "stats" / "alpha_diversity.tsv").read_text().splitlines()
    header = alpha[0].split("\t")
    assert {"richness", "shannon", "simpson"} <= set(header)
    sh = header.index("shannon")
    for row in alpha[1:]:
        cols = row.split("\t")
        assert float(cols[sh]) > 0, f"non-positive Shannon entropy in {cols}"


# --------------------------------------------------------------------------- #
# Additional modules through the real workflow: phylogenetics, second-classifier #
# consensus, and report aggregation (Krona + MultiQC).                          #
# --------------------------------------------------------------------------- #
PHYLO_FIXTURE = REPO / "tests" / "fixtures" / "phylo_demo.fasta"
KAIJU_DB = REPO / "local_databases" / "kaiju_custom"
# Kaiju's index can be built from the same committed genomes fixture, so the consensus
# cross-check runs in CI too (no gitignored DB needed) — mirrors the viral_db approach.
_CAN_BUILD_KAIJU = bool(shutil.which("prodigal") and shutil.which("kaiju-mkbwt")
                        and shutil.which("kaiju-mkfmi") and FIXTURE_GENOMES.is_file())


@pytest.fixture(scope="session")
def kaiju_db(tmp_path_factory, viral_db):
    """Path to a Kaiju protein DB: the prebuilt one if present, else built fresh from the
    30-genome fixture, reusing viral_db's taxonomy/ so taxids line up for the cross-check."""
    if (KAIJU_DB / "kaiju_db.fmi").is_file():
        return KAIJU_DB
    db_dir = tmp_path_factory.mktemp("kaiju_db")
    res = dbbuild.build_kaiju_db(str(FIXTURE_GENOMES), str(db_dir),
                                 taxonomy_dir=str(Path(viral_db) / "taxonomy"),
                                 threads=4, run=True)
    assert res.get("ok"), f"Kaiju DB build failed: {res}"
    return db_dir

_HAVE_PHYLO = bool(shutil.which("mafft") and (shutil.which("iqtree2")
                   or shutil.which("iqtree3") or shutil.which("iqtree")
                   or shutil.which("fasttree") or shutil.which("FastTree")))
_HAVE_KAIJU = bool(shutil.which("kaiju") and shutil.which("kaiju2table")
                   and ((KAIJU_DB / "kaiju_db.fmi").is_file() or _CAN_BUILD_KAIJU))
_HAVE_AGG = bool(shutil.which("ktImportText") and shutil.which("multiqc"))


@pytest.mark.skipif(not (_HAVE_PHYLO and PHYLO_FIXTURE.is_file()),
                    reason="MAFFT + IQ-TREE/FastTree not available (skips in CI without them)")
def test_phylogenetics_module(tmp_path):
    """phylogenetics module through `metagx run`: MAFFT align -> IQ-TREE/FastTree -> tree.

    Also guards against tool-version drift: the script resolves the IQ-TREE binary across
    iqtree2/iqtree3/iqtree rather than hardcoding `iqtree2` (which broke on IQ-TREE 3), and
    avoids `from __future__` (which a Snakemake `script:` preamble turns into a SyntaxError).
    """
    sheet = tmp_path / "p.samples.tsv"
    (tmp_path / "dummy.fasta").write_text(">d\nACGTACGTACGTACGT\n")
    sheet.write_text("sample\tr1\tr2\tplatform\tlayout\n"
                     f"d\t{tmp_path/'dummy.fasta'}\t\tillumina\tse\n")
    cfg = {
        "project": "phylo", "outdir": str(tmp_path / "out"), "threads": 2,
        "samples": str(sheet), "db": {"kraken2": "x", "bracken": "x"},
        "modules": {"qc": False, "classify": False, "abundance": False,
                    "assembly": False, "phylogenetics": True},
        "phylogenetics": {"input": str(PHYLO_FIXTURE), "method": "iqtree",
                          "sequence_type": "nt", "trim": False},
    }
    cfg_path = tmp_path / "p.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=2, dry_run=False)
    assert proc.returncode == 0, f"phylogenetics failed:\n{(proc.stderr or proc.stdout)[-3000:]}"

    out = tmp_path / "out" / "phylo" / "phylogenetics"
    tree = (out / "tree.nwk").read_text()
    js = json.loads((out / "phylogenetics.json").read_text())
    assert (out / "aligned.fasta").is_file()
    assert js["n_sequences"] == 6
    # all six fixture taxa made it into the tree
    for taxon in ("taxonA1", "taxonA2", "taxonA3", "taxonB1", "taxonB2", "taxonB3"):
        assert taxon in tree, f"{taxon} missing from tree"
    assert "iqtree" in js["tree_method"] or js["tree_method"] == "fasttree"


@requires_stack
@pytest.mark.skipif(not _HAVE_KAIJU,
                    reason="kaiju/kaiju2table or kaiju_custom DB absent (skips in CI)")
def test_consensus_module_kaiju(tmp_path, viral_db, kaiju_db):
    """classify_consensus module: kraken2 cross-checked against an independent protein
    classifier (Kaiju). On the viral set the two orthogonal methods should agree strongly."""
    ont = DATA / "simulated_metagenomic_reads.fasta"
    if not ont.is_file():
        ont = FIXTURE_ONT
    reads = _subsample_fasta(ont, tmp_path / "ont.fasta", 1500)
    sheet = tmp_path / "c.samples.tsv"
    sheet.write_text("sample\tr1\tr2\tplatform\tlayout\n"
                     f"c\t{reads}\t\tont\tse\n")
    cfg = {
        "project": "cons", "outdir": str(tmp_path / "out"), "threads": 4,
        "samples": str(sheet),
        "db": {"kraken2": str(viral_db), "bracken": str(viral_db), "kaiju": str(kaiju_db)},
        "modules": {"qc": False, "classify": True, "abundance": False,
                    "assembly": False, "classify_consensus": True},
        "consensus": {"classifier": "kaiju"}, "kraken2": {"minimum_hit_groups": 2},
    }
    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=4, dry_run=False)
    assert proc.returncode == 0, f"consensus failed:\n{(proc.stderr or proc.stdout)[-3000:]}"

    js = json.loads((tmp_path / "out" / "cons" / "consensus" / "c.consensus.json").read_text())
    assert js["second_classifier"] == "kaiju"
    assert js["n_species_kraken2"] > 10 and js["n_species_kaiju"] > 10
    assert js["n_shared"] >= 10
    assert js["jaccard"] >= 0.5, f"two classifiers disagree too much: jaccard={js['jaccard']}"


@requires_stack
@pytest.mark.skipif(not _HAVE_AGG,
                    reason="Krona (ktImportText)/MultiQC absent (skips in CI)")
def test_aggregate_module_krona_with_minimizer_report(tmp_path, viral_db):
    """aggregate module (Krona + MultiQC) — run with kraken2 --report-minimizer-data ON, so it
    also proves the minimizer-format Krona fix in the *real* rule: the chart must show taxon
    names, not the rank codes a fixed-index parser would emit from the shifted columns."""
    ont = DATA / "simulated_metagenomic_reads.fasta"
    if not ont.is_file():
        ont = FIXTURE_ONT
    reads = _subsample_fasta(ont, tmp_path / "ont.fasta", 1500)
    sheet = tmp_path / "a.samples.tsv"
    sheet.write_text("sample\tr1\tr2\tplatform\tlayout\n"
                     f"a\t{reads}\t\tont\tse\n")
    cfg = {
        "project": "agg", "outdir": str(tmp_path / "out"), "threads": 4,
        "samples": str(sheet),
        "db": {"kraken2": str(viral_db), "bracken": str(viral_db)},
        "modules": {"qc": False, "classify": True, "abundance": True,
                    "assembly": False, "aggregate": True},
        "kraken2": {"minimum_hit_groups": 2, "report_minimizer_data": True},
        "bracken": {"read_length": 1000, "level": "S", "threshold": 0},
    }
    cfg_path = tmp_path / "a.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=4, dry_run=False)
    assert proc.returncode == 0, f"aggregate failed:\n{(proc.stderr or proc.stdout)[-3000:]}"

    report = tmp_path / "out" / "agg" / "report"
    krona = (report / "krona.html").read_text()
    assert (report / "multiqc" / "multiqc_report.html").is_file()
    # real taxon names present (would be absent if minimizer columns mis-parsed the name)
    assert any(v in krona for v in KNOWN_SPECIES), "Krona chart has no known virus name"


# --------------------------------------------------------------------------- #
# Pre-assembled contigs input (users with isolate genomes / MAGs / references) #
# --------------------------------------------------------------------------- #
@requires_stack
def test_provided_contigs_skip_assembly(tmp_path, viral_db):
    """A sample can supply pre-assembled contigs; the assembler is skipped and downstream
    modules consume them. Here the 30 viral genomes are fed as `contigs` and reconciled
    against ONT reads — proving the staging rule wires contigs into the contig-consuming path
    without running MEGAHIT/Flye (which also dodges the MEGAHIT-on-osx-64 segfault)."""
    ont = DATA / "simulated_metagenomic_reads.fasta"
    if not ont.is_file():
        ont = FIXTURE_ONT
    reads = _subsample_fasta(ont, tmp_path / "ont.fasta", 1500)
    sheet = tmp_path / "pc.samples.tsv"
    sheet.write_text(
        "sample\tr1\tr2\tplatform\tlayout\tlibrary\tcontigs\n"
        f"s\t{reads}\t\tont\tse\twgs\t{FIXTURE_GENOMES}\n")
    cfg = {
        "project": "pc", "outdir": str(tmp_path / "out"), "threads": 4,
        "samples": str(sheet),
        "db": {"kraken2": str(viral_db), "bracken": str(viral_db)},
        "modules": {"qc": False, "classify": True, "abundance": False,
                    "assembly": True, "binning": True, "reconcile": True},
        "kraken2": {"minimum_hit_groups": 2}, "bracken": {"read_length": 1000},
    }
    cfg_path = tmp_path / "pc.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=4, dry_run=False)
    assert proc.returncode == 0, f"provided-contigs run failed:\n{(proc.stderr or proc.stdout)[-3000:]}"

    out = tmp_path / "out" / "pc"
    staged = out / "assembly" / "s" / "final.contigs.fa"
    assert staged.is_file(), "provided contigs were not staged as the assembly output"
    # staged contigs are exactly the user's genomes (30 records), not a re-assembly
    assert staged.read_text().count(">") == 30
    # the assembler must NOT have run (no Flye/MEGAHIT log dirs)
    assert not (out / "assembly" / "s" / "final.contigs.fa.flye_tmp").exists()
    # reconcile classified the provided contigs and recovered known DB viruses
    per_taxon = (out / "reconcile" / "s.reconciliation.tsv").read_text()
    assert any(k in per_taxon for k in KNOWN_SPECIES)


# --------------------------------------------------------------------------- #
# Bacterial AMR screening (functional module) via the contigs input          #
# --------------------------------------------------------------------------- #
AMR_PLASMID = DATA / "bacteria" / "amr_plasmid.fasta"
_HAVE_ABRICATE = bool(shutil.which("abricate") and shutil.which("blastn"))


@pytest.mark.skipif(not (_HAVE_ABRICATE and AMR_PLASMID.is_file()),
                    reason="abricate+blastn or the AMR plasmid fixture absent "
                    "(provisioned via `metagx run --use-conda`; skips otherwise)")
def test_amr_screening_on_provided_genome(tmp_path):
    """functional/AMR (ABRicate) on a pre-assembled bacterial genome — the realistic
    'I have a genome, screen it for resistance' end-user flow. The pNDM-HK plasmid carries
    blaNDM-1, so a working ABRicate must report at least one resistance gene."""
    sheet = tmp_path / "amr.samples.tsv"
    sheet.write_text(
        "sample\tr1\tr2\tplatform\tlayout\tlibrary\tcontigs\n"
        f"plasmid\t\t\tillumina\tse\twgs\t{AMR_PLASMID}\n")
    cfg = {
        "project": "amr", "outdir": str(tmp_path / "out"), "threads": 4,
        "samples": str(sheet), "db": {"kraken2": "x", "bracken": "x"},
        "modules": {"qc": False, "classify": False, "abundance": False,
                    "assembly": True, "functional": True},
        "functional": {"amr": True, "pathways": False, "annotation": False},
    }
    cfg_path = tmp_path / "amr.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=4, dry_run=False)
    assert proc.returncode == 0, f"AMR run failed:\n{(proc.stderr or proc.stdout)[-3000:]}"

    tsv = (tmp_path / "out" / "amr" / "functional" / "plasmid" / "amr" / "abricate.tsv").read_text()
    rows = [l for l in tsv.splitlines() if l and not l.startswith("#")]
    assert len(rows) >= 1, f"ABRicate found no AMR genes on a blaNDM-1 plasmid:\n{tsv[:500]}"
    assert "NDM" in tsv.upper() or "BLA" in tsv.upper(), f"expected a bla/NDM hit:\n{tsv[:800]}"
