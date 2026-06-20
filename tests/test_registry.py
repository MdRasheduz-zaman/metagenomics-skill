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


def test_when_matches_equality_and_numeric():
    assert registry.when_matches({"goal": "strain_resolved"}, {"goal": "Strain_Resolved"})
    assert not registry.when_matches({"goal": "strain_resolved"}, {"goal": "diversity"})
    assert registry.when_matches({"estimated_bases_gte": 5e10}, {"estimated_bases": 6e10})
    assert not registry.when_matches({"estimated_bases_gte": 5e10}, {"estimated_bases": 1e9})
    # unknown facts never match — we don't promote on what we don't know
    assert not registry.when_matches({"goal": "strain_resolved"}, {})
    assert not registry.when_matches({}, {"goal": "strain_resolved"})


def test_promote_when_surfaces_quiet_param_with_reason():
    base = {p["name"] for p in registry.interview_spec("flye", max_tier=2)}
    assert "keep_haplotypes" not in base  # quiet by default (ask:false/tier:3)

    promoted = registry.interview_spec("flye", max_tier=2, context={"goal": "strain_resolved"})
    kh = next(p for p in promoted if p["name"] == "keep_haplotypes")
    assert kh["tier"] == 1                      # pulled to the promoted tier
    assert "promoted" in kh and kh["promoted"]["reason"]
    assert kh == promoted[0] or promoted[0]["tier"] == 1  # promoted sorts to the front


def test_promote_when_numeric_threshold():
    none_ctx = {p["name"] for p in registry.interview_spec("flye", max_tier=2)}
    assert "asm_coverage" not in none_ctx
    deep = {p["name"] for p in
            registry.interview_spec("flye", max_tier=2, context={"estimated_bases": 6e10})}
    assert "asm_coverage" in deep
    shallow = {p["name"] for p in
               registry.interview_spec("flye", max_tier=2, context={"estimated_bases": 1e9})}
    assert "asm_coverage" not in shallow


@pytest.mark.parametrize("tool", registry.list_tools())
def test_every_tool_has_extra_args_valve(tool):
    """Capability-completeness floor: every registry exposes a working passthrough."""
    spec = registry.load_registry(tool)["params"].get("extra_args")
    assert spec and spec.get("passthrough") and spec["type"] == "str"
    assert registry.render_args(tool, {"extra_args": "--x 1"})[-2:] == ["--x", "1"]


def test_kraken2_promote_and_passthrough():
    base = {p["name"] for p in registry.interview_spec("kraken2", max_tier=2)}
    assert "report_minimizer_data" not in base and "quick" not in base
    audit = {p["name"] for p in
             registry.interview_spec("kraken2", max_tier=2, context={"goal": "false_positive_audit"})}
    assert "report_minimizer_data" in audit
    fast = {p["name"] for p in
            registry.interview_spec("kraken2", max_tier=2, context={"goal": "fast_screen"})}
    assert "quick" in fast
    # passthrough valve renders raw, bzip2 is managed-rejected
    assert registry.render_args("kraken2", {"extra_args": "--minimizer-spaces 7"})[-2:] == \
        ["--minimizer-spaces", "7"]
    with pytest.raises(registry.ValidationError):
        registry.validate("kraken2", {"bzip2_compressed": True})


def test_promotion_respects_max_tier_ceiling():
    """A param promoted to tier 1 still must not appear if to_tier exceeds max_tier."""
    # keep_haplotypes promotes to tier 1, so it shows even at max_tier=1
    got = {p["name"] for p in
           registry.interview_spec("flye", max_tier=1, context={"goal": "strain_resolved"})}
    assert "keep_haplotypes" in got


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
        # A flagless param is only legitimate as `managed`, `interpreted`, or `passthrough`
        # (the raw escape hatch). Anything else with no flag is a forgotten flag.
        if not any(spec.get(k) for k in ("managed", "interpreted", "passthrough")):
            assert spec.get("flag"), f"{tool}.{name}: user param needs a flag (or interpreted/passthrough)"
        # `interpreted` and `managed` are mutually exclusive intents.
        assert not (spec.get("managed") and spec.get("interpreted")), \
            f"{tool}.{name}: cannot be both managed and interpreted"
        # passthrough is a raw string of CLI tokens and must not carry a flag of its own.
        if spec.get("passthrough"):
            assert spec["type"] == "str", f"{tool}.{name}: passthrough must be type str"
            assert not spec.get("flag"), f"{tool}.{name}: passthrough must not declare a flag"
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


def test_passthrough_renders_raw_tokens_after_flags():
    """A passthrough value is split shell-style and appended verbatim, after real flags."""
    args = registry.render_args("flye", {"meta": True, "extra_args": "--resume --stop-after consensus"})
    assert "--meta" in args
    assert args[-3:] == ["--resume", "--stop-after", "consensus"]


def test_passthrough_empty_is_skipped():
    assert registry.render_args("flye", {"extra_args": ""}) == []


def test_passthrough_value_validates_as_str():
    """The valve goes through normal validation as a plain string (not managed-rejected)."""
    cleaned = registry.validate("flye", {"extra_args": "--resume"})
    assert cleaned == {"extra_args": "--resume"}


def test_flye_read_type_flags_are_managed():
    """User cannot set a read-type input flag — the workflow picks it by platform."""
    with pytest.raises(registry.ValidationError):
        registry.validate("flye", {"nano_raw": "reads.fq"})


def test_new_gap_closing_registries_present():
    tools = registry.list_tools()
    assert "antismash" in tools and "dada2" in tools
