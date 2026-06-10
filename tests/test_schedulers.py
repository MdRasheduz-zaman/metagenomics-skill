"""Scheduler-backend registry + bundled profile integrity."""

import os

import pytest
import yaml

from metagx import schedulers


def test_known_schedulers_present():
    names = schedulers.list_schedulers()
    for n in ("local", "slurm", "lsf", "sge", "pbs", "generic"):
        assert n in names


def test_every_scheduler_has_a_valid_bundled_profile():
    """Each declared backend resolves to a real profile dir with a parseable config."""
    for name in schedulers.list_schedulers():
        path = schedulers.profile_path(name)
        cfg_file = os.path.join(path, "config.yaml")
        assert os.path.isfile(cfg_file), f"{name}: missing {cfg_file}"
        cfg = yaml.safe_load(open(cfg_file))
        assert isinstance(cfg, dict) and cfg, f"{name}: empty/invalid profile"


def test_unknown_scheduler_raises_with_choices():
    with pytest.raises(KeyError) as e:
        schedulers.profile_path("does-not-exist")
    # the error lists the valid names so the user can recover
    assert "slurm" in str(e.value)


def test_generic_cluster_submit_cmds_have_no_embedded_comments():
    """A `#` inside the quoted submit string would truncate the shell command."""
    for name in ("sge", "pbs", "generic"):
        cfg = yaml.safe_load(open(os.path.join(schedulers.profile_path(name), "config.yaml")))
        cmd = cfg.get("cluster-generic-submit-cmd", "")
        assert cmd, f"{name}: cluster-generic profile needs a submit cmd"
        assert "#" not in cmd, f"{name}: submit cmd contains '#' (would be a shell comment)"


def test_cluster_profiles_declare_executor_and_threads():
    """Submit-based backends must name an executor and pass {threads} through."""
    for name in ("slurm", "lsf", "sge", "pbs", "generic"):
        cfg = yaml.safe_load(open(os.path.join(schedulers.profile_path(name), "config.yaml")))
        assert cfg.get("executor"), f"{name}: profile must set an executor"
        # cluster-generic backends must thread {threads} into the submit command
        if cfg["executor"] == "cluster-generic":
            assert "{threads}" in cfg["cluster-generic-submit-cmd"], \
                f"{name}: submit cmd should request {{threads}} cores"


def test_describe_rows_carry_metadata():
    for row in schedulers.describe():
        for key in ("name", "executor", "plugin", "summary", "edit"):
            assert row.get(key), f"{row.get('name')}: missing {key}"
