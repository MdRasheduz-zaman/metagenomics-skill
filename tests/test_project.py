"""`metagx project` — the tool generating collision-safe, runnable analysis folders for any config.

The core guarantee: a config's tools are SPLIT into a tailored core env (tools that coexist in
environment.yml) and `--use-conda` tools (isolated per-rule envs), and a known collider (abricate,
which hard-pins samtools 0.1.x) is NEVER placed in the shared tailored env.
"""
import shutil
import subprocess

import pytest

from metagx import project, runner

ENV = runner.environment_file_path()


def _cfg(modules, platform="illumina"):
    return {"project": "p", "modules": modules,
            "samples": [{"sample": "s", "r1": "x", "platform": platform}]}


def test_pure_core_config_needs_no_per_rule_envs():
    plan = project.plan_env(_cfg({"qc": True, "classify": True, "abundance": True}, "ont"), ENV)
    bases = {p.split()[0] for p in plan["core_packages"]}
    assert bases == {"kraken2", "bracken", "chopper", "porechop_abi"}
    assert plan["per_rule_tools"] == []
    assert plan["use_conda"] is False       # tailored env covers everything


def test_collider_goes_to_use_conda_never_the_tailored_env():
    plan = project.plan_env(_cfg({"qc": True, "classify": True, "assembly": True, "functional": True}), ENV)
    core = {p.split()[0] for p in plan["core_packages"]}
    assert "abricate" not in core           # the collider must NOT be in the shared env
    assert "abricate" in plan["per_rule_tools"]
    assert plan["use_conda"] is True
    # the coexisting core tools still make the tailored env
    assert {"kraken2", "fastp", "megahit"} <= core


def test_scaffold_writes_full_folder_and_matches_plan():
    written = {}
    cfg = _cfg({"qc": True, "classify": True, "abundance": True}, "ont")
    cfg["project"] = "gut"
    res = project.scaffold(cfg, "OUT", executor="slurm", platform="ont", env_yml_path=ENV,
                           writer=lambda p, t: written.__setitem__(p, t))
    assert {"config.yaml", "samples.tsv", "env.yaml", "00_setup.sh", "run.sh", "README.md"} \
        <= set(res["files"])
    assert any(f.endswith("profile/config.yaml") for f in res["files"])
    # pure-core config -> profile use-conda false, and no --use-conda in the run command
    prof = next(t for p, t in written.items() if p.endswith("profile/config.yaml"))
    assert "use-conda: false" in prof
    run = next(t for p, t in written.items() if p.endswith("run.sh"))
    assert "--use-conda" not in run


def test_scaffold_amr_sets_use_conda_true_in_profile_and_run():
    written = {}
    cfg = _cfg({"qc": True, "classify": True, "assembly": True, "functional": True})
    project.scaffold(cfg, "OUT", executor="slurm", platform="illumina", env_yml_path=ENV,
                     writer=lambda p, t: written.__setitem__(p, t))
    prof = next(t for p, t in written.items() if p.endswith("profile/config.yaml"))
    run = next(t for p, t in written.items() if p.endswith("run.sh"))
    env = next(t for p, t in written.items() if p.endswith("env.yaml"))
    assert "use-conda: true" in prof
    assert "--use-conda" in run
    assert "abricate" not in env             # collider isolated, not in the tailored env


def _gen_scripts(modules, platform="ont", executor="slurm"):
    written = {}
    cfg = _cfg(modules, platform); cfg["project"] = "t"
    project.scaffold(cfg, "OUT", executor=executor, platform=platform, env_yml_path=ENV,
                     writer=lambda p, t: written.__setitem__(p, t))
    return written


def test_setup_notes_are_echo_lines_not_bare_text():
    """Regression: the env note was once inserted as bare shell text (bash ran 'Every' as a
    command, and a --use-conda branch had backticks that would command-substitute)."""
    for mods in ({"qc": True, "classify": True, "abundance": True},                       # use_conda False
                 {"qc": True, "classify": True, "assembly": True, "functional": True}):    # use_conda True
        setup = next(t for p, t in _gen_scripts(mods).items() if p.endswith("00_setup.sh"))
        for line in setup.splitlines():
            if "per-rule env" in line:
                assert line.lstrip().startswith("echo"), f"bare-text note: {line!r}"
        assert "`" not in setup, "no backticks (command substitution) in the generated script"


@pytest.mark.skipif(not shutil.which("bash"), reason="bash not available")
def test_generated_scripts_pass_bash_syntax_check():
    for mods in ({"qc": True, "classify": True, "abundance": True},
                 {"qc": True, "classify": True, "assembly": True, "functional": True}):
        for name in ("00_setup.sh", "run.sh"):
            text = next(t for p, t in _gen_scripts(mods).items() if p.endswith(name))
            r = subprocess.run(["bash", "-n"], input=text, capture_output=True, text=True)
            assert r.returncode == 0, f"{name} failed bash -n: {r.stderr}"
