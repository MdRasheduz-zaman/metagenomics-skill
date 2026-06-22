"""CI dry-run gate over the actual Snakemake workflow.

This is the test that was missing: the unit suite verifies the registries and the
Python helpers in isolation, but nothing executed the *workflow* — so a renamed
registry param left dangling in a ``.smk`` rule, a broken ``include:`` graph, a
mistyped managed key, or a wildcard regression all passed CI green.

A Snakemake dry-run with ``--printshellcmds`` (which ``runner.run`` always passes)
builds the full DAG **and resolves every ``params:`` lambda**, which is where the
rules call ``registry.render_args``. So this exercises the registry→rule→command
seam for kraken2/Bracken/fastp/MEGAHIT/Flye without needing the bio tools or any
database on disk — only tiny synthetic read files, which we create in ``tmp_path``.

It needs only ``snakemake`` (a core dependency), so it runs in CI. For a real
execution test against the viral DB, see ``test_pipeline_e2e.py``.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pytest
import yaml

from metagx import runner

# A short but legal read; sequence content is irrelevant to DAG construction.
_SEQ = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"


def _write_fastq(path: Path, n: int = 4, gz: bool = False) -> None:
    rec = "".join(
        f"@r{i}\n{_SEQ}\n+\n{'I' * len(_SEQ)}\n" for i in range(n)
    )
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(rec)
    else:
        path.write_text(rec)


def _write_fasta(path: Path, n: int = 4) -> None:
    path.write_text("".join(f">r{i}\n{_SEQ}\n" for i in range(n)))


# Each scenario exercises a different platform/layout path through common.smk's
# routing and a different set of rules' render_args calls.
SCENARIOS = {
    "illumina_pe": {
        "reads": {"r1": ("ill_R1.fastq.gz", "fastq_gz"),
                  "r2": ("ill_R2.fastq.gz", "fastq_gz")},
        "platform": "illumina", "layout": "pe",
        "modules": {"qc": True, "classify": True, "abundance": True, "assembly": True},
        "extra": {"kraken2": {"minimum_base_quality": 20}, "bracken": {"read_length": 150}},
        "expect_rules": {"fastp_pe", "kraken2", "bracken", "megahit"},
    },
    "ont_fasta": {
        "reads": {"r1": ("ont.fasta", "fasta")},
        "platform": "ont", "layout": "se",
        "modules": {"qc": False, "classify": True, "abundance": True,
                    "assembly": True, "binning": True, "reconcile": True},
        "extra": {"bracken": {"read_length": 1000}, "flye": {"meta": True},
                  "db_cat": True},
        "expect_rules": {"kraken2", "bracken", "flye", "metabat2", "reconcile"},
    },
    "pacbio_hifi_se": {
        "reads": {"r1": ("pb.fastq.gz", "fastq_gz")},
        "platform": "pacbio_hifi", "layout": "se",
        "modules": {"qc": True, "classify": True, "abundance": True},
        "extra": {"bracken": {"read_length": 1000}},
        "expect_rules": {"pacbio_qc", "kraken2", "bracken"},
    },
    # Domain taxonomy/quality: the rules that were rewired from hardcoded shell flags to
    # registry.render_args (gtdbtk/checkm2/checkv/eukcc). This guards that wiring.
    "ont_domains": {
        "reads": {"r1": ("ont.fasta", "fasta")},
        "platform": "ont", "layout": "se",
        "modules": {"qc": False, "classify": True, "abundance": True,
                    "assembly": True, "binning": True, "domain_taxonomy": True},
        "extra": {"bracken": {"read_length": 1000}, "flye": {"meta": True},
                  "domains": ["prokaryote", "viral", "eukaryote"],
                  "gtdbtk": {"pplacer_cpus": 1}, "checkm2": {"lowmem": True}},
        "expect_rules": {"flye", "gtdbtk", "checkm2", "checkv", "eukcc"},
    },
    # Pre-assembled contigs + functional/AMR: the staging rule replaces the assembler, and the
    # functional sub-flag selection (amr only) is exercised. No reads needed for this sample.
    "provided_contigs_amr": {
        "reads": {"r1": ("reads.fasta", "fasta")},   # reads present but assembly is staged
        "platform": "illumina", "layout": "se",
        "modules": {"qc": False, "classify": False, "abundance": False,
                    "assembly": True, "functional": True},
        "extra": {"contigs": True, "functional": {"amr": True, "pathways": False,
                                                   "annotation": False}},
        "expect_rules": {"stage_provided_contigs", "abricate"},
    },
    # Amplicon (marker-gene) branch — was registry/config-only, never DAG-tested. A sample
    # with library=amplicon routes through cutadapt primer trim -> VSEARCH OTUs (the default).
    "amplicon_otu": {
        "reads": {"r1": ("amp.fastq.gz", "fastq_gz")},
        "platform": "illumina", "layout": "se",
        "modules": {"qc": False, "classify": False, "abundance": False, "assembly": False},
        "extra": {"library": "amplicon"},
        "expect_rules": {"vsearch_otus"},
    },
    # Strain (inStrain) + BGC (antiSMASH) on an assembled WGS sample — both were config-only.
    # Exercises their render_args wiring through the real DAG.
    "strain_bgc": {
        "reads": {"r1": ("wgs.fastq.gz", "fastq_gz")},
        "platform": "illumina", "layout": "se",
        "modules": {"qc": False, "classify": False, "abundance": False,
                    "assembly": True, "binning": True, "strain": True, "bgc": True},
        "extra": {},
        "expect_rules": {"instrain", "antismash"},
    },
    # Cross-sample statistics (α/β diversity) over the combined Bracken table — was script-unit
    # tested but its rule was never resolved through the workflow.
    "stats_diversity": {
        "reads": {"r1": ("s.fastq.gz", "fastq_gz")},
        "platform": "illumina", "layout": "se",
        "modules": {"qc": False, "classify": True, "abundance": True,
                    "assembly": False, "stats": True},
        "extra": {"bracken": {"read_length": 150}},
        "expect_rules": {"diversity"},
    },
}


def _materialise(sc: dict, tmp_path: Path) -> Path:
    """Write synthetic reads + sample sheet + config; return the config path."""
    row = {"sample": "s", "platform": sc["platform"], "layout": sc["layout"]}
    for key, (fname, kind) in sc["reads"].items():
        fp = tmp_path / fname
        if kind == "fasta":
            _write_fasta(fp)
        else:
            _write_fastq(fp, gz=kind.endswith("gz"))
        row[key] = str(fp)
    row.setdefault("r2", "")
    if sc["extra"].get("library"):
        row["library"] = sc["extra"]["library"]

    # optional pre-assembled contigs input
    if sc["extra"].get("contigs"):
        fp = tmp_path / "contigs.fasta"
        _write_fasta(fp, n=8)
        row["contigs"] = str(fp)

    sheet = tmp_path / "samples.tsv"
    cols = ["sample", "r1", "r2", "platform", "layout", "library", "contigs"]
    sheet.write_text(
        "\t".join(cols) + "\n" + "\t".join(str(row.get(c, "")) for c in cols) + "\n"
    )

    db = {"kraken2": str(tmp_path / "db"), "bracken": str(tmp_path / "db")}
    if sc["extra"].get("db_cat"):
        db["cat"] = str(tmp_path / "catdb")
    if sc["modules"].get("domain_taxonomy"):
        for tool in ("gtdbtk", "checkm2", "checkv", "eukcc", "genomad"):
            db[tool] = str(tmp_path / f"{tool}_db")
    cfg = {
        "project": "dryrun",
        "outdir": str(tmp_path / "out"),
        "threads": 2,
        "samples": str(sheet),
        "db": db,
        "modules": sc["modules"],
    }
    for k, v in sc["extra"].items():
        if k not in ("db_cat", "contigs", "library"):  # sample-sheet/marker keys, not config
            cfg[k] = v
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))
    return cfg_path


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_workflow_dry_run_builds_and_renders(name, tmp_path):
    """The DAG builds and every params lambda (render_args) resolves cleanly."""
    sc = SCENARIOS[name]
    cfg = _materialise(sc, tmp_path)

    proc = runner.run(config=str(cfg), cores=2, dry_run=True)
    combined = (proc.stdout or "") + (proc.stderr or "")

    assert proc.returncode == 0, (
        f"snakemake dry-run failed for scenario '{name}':\n{combined[-3000:]}"
    )
    # The shell commands are resolved in a dry-run (--printshellcmds), so the
    # expected rules — and thus their render_args output — must appear.
    for rule in sc["expect_rules"]:
        assert rule in combined, (
            f"scenario '{name}': expected rule '{rule}' missing from the plan.\n"
            f"{combined[-2000:]}"
        )


def test_db_build_step_wires_into_dag(tmp_path):
    """A configured db.build adds build_kraken2_db as an upstream dependency of kraken2, and
    the classify-time Bracken -r matches the per-platform length the DB is built for (the
    coupling that makes a detached CLI build error-prone)."""
    from metagx import config_builder as cb

    genomes = tmp_path / "genomes.fasta"
    _write_fasta(genomes, n=3)
    reads = tmp_path / "ont.fasta"
    _write_fasta(reads, n=4)
    cfg = cb.build_config(
        project="dbb", outdir=str(tmp_path / "out"),
        samples=[{"sample": "s", "r1": str(reads), "platform": "ont", "layout": "se"}],
        db={"kraken2": str(tmp_path / "viraldb"),
            "build": {"strategy": "custom-fasta", "taxonomy": "synthetic", "source": str(genomes)}},
        modules={"classify": True, "abundance": True})
    cfg_path = tmp_path / "config.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg))

    proc = runner.run(config=str(cfg_path), cores=2, dry_run=True)
    combined = (proc.stdout or "") + (proc.stderr or "")
    assert proc.returncode == 0, combined[-3000:]
    assert "build_kraken2_db" in combined            # the build runs as a pipeline step
    assert ".metagx_db.json" in combined             # kraken2 depends on the build manifest
    assert "bracken -r 1000" in combined             # ont -> 1000, matching the built DB


def test_managed_key_typo_raises_loudly():
    """The safety net behind the gate: an orphaned/mistyped managed key must raise.

    The parametrized dry-runs above prove the *real* rules pass only valid managed
    keys (returncode 0). This proves the failure mode itself is loud — a registry
    rename that leaves a ``.smk`` referencing a dead managed key cannot silently drop
    the flag and emit a corrupt command line.
    """
    from metagx import registry

    with pytest.raises(registry.ValidationError) as exc:
        registry.render_args("kraken2", {}, managed={"reprot": "x"})  # typo of `report`
    assert "reprot" in str(exc.value)

    # the valid spelling renders without error
    args = registry.render_args("kraken2", {}, managed={"report": "x"})
    assert "--report" in args and "x" in args
