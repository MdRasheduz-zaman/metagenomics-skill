"""`metagx env-file --config` — a MINIMAL conda env tailored to a config's active tools.

The whole point: a gut-profiling run must not drag in the full assembly/binning/domain stack.
active_core_packages intersects report.active_tools(cfg) with the bundled environment.yml, so
heavy/domain tools (which self-provision per-rule via --use-conda) are excluded by construction.
"""
import re
import textwrap

import yaml

from metagx import report, runner


def _fake_env_yml(tmp_path):
    p = tmp_path / "environment.yml"
    p.write_text(textwrap.dedent("""\
        name: metagx
        channels: [conda-forge, bioconda]
        dependencies:
          - python >=3.10
          - kraken2 >=2.1
          - bracken >=2.9
          - fastp >=0.23
          - chopper >=0.7
          - porechop_abi >=0.5
          - megahit >=1.2
          - flye >=2.9
          - minimap2 >=2.26
          - samtools >=1.18
          - metabat2 >=2.15
          - snakemake-minimal >=8.10
    """))
    return str(p)


def _gut_ont_cfg():
    return {"project": "gut", "modules": {"qc": True, "classify": True, "abundance": True},
            "samples": [{"sample": "s", "r1": "x.fq", "platform": "ont"}]}


def test_gut_profiling_env_is_minimal(tmp_path):
    env = _fake_env_yml(tmp_path)
    pkgs = report.active_core_packages(_gut_ont_cfg(), env)
    bases = {p.split()[0] for p in pkgs}
    assert bases == {"kraken2", "bracken", "chopper", "porechop_abi"}   # ONT QC + classify + abundance
    assert "megahit" not in bases and "flye" not in bases and "fastp" not in bases  # no junk
    # version floors from environment.yml are preserved
    assert "kraken2 >=2.1" in pkgs


def test_short_read_assembly_adds_only_relevant_tools(tmp_path):
    env = _fake_env_yml(tmp_path)
    cfg = {"project": "p", "modules": {"qc": True, "classify": True, "assembly": True, "binning": True},
           "samples": [{"sample": "s", "r1": "x.fq", "platform": "illumina"}]}
    bases = {p.split()[0] for p in report.active_core_packages(cfg, env)}
    assert "fastp" in bases and "megahit" in bases              # short-read QC + assembler
    assert {"minimap2", "samtools", "metabat2"} <= bases        # binning
    assert "flye" not in bases                                  # long-read assembler not selected


def test_config_env_yaml_is_valid_and_named(tmp_path):
    env = _fake_env_yml(tmp_path)
    text = report.config_env_yaml(_gut_ont_cfg(), env)
    doc = yaml.safe_load(text)
    assert doc["name"] == "metagx-gut"
    assert doc["channels"] == ["conda-forge", "bioconda"]
    deps = {d.split()[0] for d in doc["dependencies"]}
    assert "python" in deps and "kraken2" in deps and "bracken" in deps
    assert "megahit" not in deps
    # snakemake is NOT in the tool env (it lives in metagx's own venv)
    assert not any(d.startswith("snakemake") for d in doc["dependencies"])


def _core_env_bases():
    """Package bases in the REAL bundled environment.yml (the curated collision-free core set)."""
    path = runner.environment_file_path()
    assert path, "bundled environment.yml must resolve in a repo checkout"
    deps = (yaml.safe_load(open(path)) or {}).get("dependencies", []) or []
    return {re.split(r"[ <>=]", d, 1)[0] for d in deps if isinstance(d, str)}, deps


def test_validate_config_pulls_blast_but_gut_does_not():
    """Per-config selectivity + the fixed validate gap: BLAST appears only when modules.validate
    is on. (blast/blastn used to be provisioned by NO env — see wiring invariant I.)"""
    env = runner.environment_file_path()
    gut = {"project": "g", "modules": {"qc": True, "classify": True, "abundance": True},
           "samples": [{"sample": "s", "r1": "x", "platform": "ont"}]}
    val = {"project": "v", "modules": {"qc": True, "classify": True, "validate": True},
           "samples": [{"sample": "s", "r1": "x", "platform": "illumina"}]}
    gut_bases = {p.split()[0] for p in report.active_core_packages(gut, env)}
    val_bases = {p.split()[0] for p in report.active_core_packages(val, env)}
    assert "blast" in val_bases
    assert "blast" not in gut_bases


def test_known_collider_stays_out_of_core_env():
    """Collision safety: abricate hard-pins samtools 0.1.x (breaks samtools sort -o for the whole
    pipeline), so it MUST be isolated in its own per-rule env, never the shared core env — and the
    core env must keep the samtools >=1.18 pin."""
    bases, deps = _core_env_bases()
    assert "abricate" not in bases, "abricate must NOT be in the core env (it collides on samtools)"
    assert any(str(d).startswith("samtools") and ">=1.18" in str(d) for d in deps), \
        "core env must pin samtools >=1.18"
