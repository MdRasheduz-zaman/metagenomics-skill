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
        if not spec.get("managed"):
            assert spec.get("flag"), f"{tool}.{name}: user param needs a flag"
    # interview + render must not raise for any registry
    registry.interview_spec(tool, max_tier=3)
    registry.render_args(tool, {})


def test_new_gap_closing_registries_present():
    tools = registry.list_tools()
    assert "antismash" in tools and "dada2" in tools
