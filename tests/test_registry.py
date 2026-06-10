import pytest
from metagx import registry


def test_list_tools_has_core():
    tools = registry.list_tools()
    for t in ("kraken2", "bracken", "fastp", "flye", "genomad", "cutadapt"):
        assert t in tools


def test_render_args_bool_scalar_managed():
    args = registry.render_args(
        "kraken2", {"confidence": 0.1, "use_names": True, "quick": False},
        managed={"db": "DB", "threads": 8, "paired": True})
    assert "--confidence" in args and "0.1" in args
    assert "--use-names" in args            # bool True -> flag present
    assert "--quick" not in args            # bool False -> omitted
    assert "--db" in args and "DB" in args
    assert "--paired" in args


def test_validate_range_and_managed_and_sweep():
    with pytest.raises(registry.ValidationError):
        registry.validate("kraken2", {"confidence": 2.0})       # > max
    with pytest.raises(registry.ValidationError):
        registry.validate("kraken2", {"db": "x"})               # managed
    with pytest.raises(registry.ValidationError):
        registry.validate("kraken2", {"minimum_hit_groups": [1, 2]})  # not sweepable
    assert registry.validate("kraken2", {"confidence": [0.0, 0.1]}) == {"confidence": [0.0, 0.1]}


def test_interview_excludes_managed():
    spec = registry.interview_spec("kraken2", max_tier=3)
    names = {p["name"] for p in spec}
    assert "confidence" in names
    assert "db" not in names and "threads" not in names


VALID_TYPES = {"int", "float", "bool", "str", "path", "enum"}


@pytest.mark.parametrize("tool", registry.list_tools())
def test_every_registry_is_well_formed(tool):
    """Schema integrity for ALL registries (auto-covers newly added tools)."""
    reg = registry.load_registry(tool)
    assert reg.get("command"), f"{tool}: missing command"
    for name, spec in reg["params"].items():
        assert spec.get("type") in VALID_TYPES, f"{tool}.{name}: bad type {spec.get('type')}"
        if spec.get("type") == "enum":
            assert spec.get("choices"), f"{tool}.{name}: enum needs choices"
        # A user param needs a `flag` unless it is `interpreted` (consumed by a
        # workflow script, not rendered to a CLI flag). `managed` params are
        # injected by the workflow and also flagless. Requiring one of the three
        # keeps a genuinely forgotten flag a failure.
        if not spec.get("managed") and not spec.get("interpreted"):
            assert spec.get("flag"), f"{tool}.{name}: user param needs a flag (or interpreted: true)"
        # `interpreted` and `managed` are mutually exclusive intents.
        assert not (spec.get("managed") and spec.get("interpreted")), \
            f"{tool}.{name}: cannot be both managed and interpreted"
    # interview + render must not raise for any registry
    registry.interview_spec(tool, max_tier=3)
    registry.render_args(tool, {})


def test_interpreted_params_are_never_rendered():
    """`interpreted` params are consumed by workflow scripts, not emitted as flags."""
    for tool in registry.list_tools():
        reg = registry.load_registry(tool)
        interpreted = [n for n, s in reg["params"].items() if s.get("interpreted")]
        for name in interpreted:
            spec = reg["params"][name]
            # use the default (or a choice) as a plausible value
            val = spec.get("default")
            if val is None and spec.get("choices"):
                val = spec["choices"][0]
            args = registry.render_args(tool, {name: val})
            assert name.replace("_", "-") not in " ".join(args), \
                f"{tool}.{name} is interpreted but was rendered into args"


def test_new_gap_closing_registries_present():
    tools = registry.list_tools()
    assert "antismash" in tools and "dada2" in tools
