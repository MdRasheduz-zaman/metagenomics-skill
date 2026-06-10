"""`metagx validate` must forward EVERY build_config parameter present in the config.

Regression guard for a drift bug: cmd_validate round-trips the YAML through
build_config to reuse validation, but used to forward only a hardcoded subset of
kwargs. Tier 2/3 sections (assembly/metaspades/consensus/instrain/...) were
dropped, so valid hybrid configs were falsely rejected and bad params in those
sections passed validation silently. cmd_validate now forwards generically via
inspect.signature; these tests lock that in.
"""
import inspect
import types

import yaml

from metagx import cli, config_builder


def _validate(tmp_path, cfg) -> int:
    p = tmp_path / "config.yaml"
    p.write_text(yaml.safe_dump(cfg))
    return cli.cmd_validate(types.SimpleNamespace(config=str(p)))


BASE = {"project": "t", "db": {"kraken2": "k"}}
SPE = [{"sample": "a", "r1": "a_1.fq.gz", "r2": "a_2.fq.gz"}]


def test_validate_hybrid_metaspades_is_accepted(tmp_path):
    # assembly.assembler must survive the round-trip (was dropped -> false reject).
    hyb = [{"sample": "a", "r1": "a_1.fq.gz", "r2": "a_2.fq.gz",
            "long_reads": "a_ont.fq.gz", "long_platform": "ont"}]
    cfg = {**BASE, "samples": hyb, "modules": {"assembly": True},
           "assembly": {"assembler": "metaspades"}}
    assert _validate(tmp_path, cfg) == 0


def test_validate_catches_bad_tier3_tool_param(tmp_path):
    # instrain is a Tier-3 section the old validate never forwarded -> never checked.
    cfg = {**BASE, "samples": SPE, "modules": {"assembly": True, "strain": True},
           "instrain": {"not_a_real_flag": 1}}
    assert _validate(tmp_path, cfg) == 2


def test_cmd_validate_forwards_all_build_config_params():
    """Completeness guard: don't regress to a hardcoded forward list.

    Every keyword of build_config (minus the ones cmd_validate fills itself)
    must be a config key it can forward. With the inspect-based implementation
    this is true by construction; this test fails loudly if someone reverts to
    enumerating kwargs by hand and forgets one.
    """
    src = inspect.getsource(cli.cmd_validate)
    assert "inspect.signature(config_builder.build_config)" in src, (
        "cmd_validate must forward build_config params generically, not via a "
        "hardcoded list (that list silently drifts as tools are added)."
    )
