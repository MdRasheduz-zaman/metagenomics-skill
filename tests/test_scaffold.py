import yaml

from metagx import scaffold

# A realistic slice of `flye --help` covering the argparse formats the parser must handle:
# short+long with metavar, long-only with metavar, bare bool, wrapped description, enum braces.
FLYE_HELP = """\
usage: flye [options]

optional arguments:
  -h, --help            show this help message and exit
  --nano-raw path [path ...]
                        ONT regular reads, pre-Guppy5 (<20% error)
  -g SIZE, --genome-size SIZE
                        estimated genome size (for example, 5m or 2.6g)
  -t INT, --threads INT
                        number of parallel threads [1]
  -i INT, --iterations INT
                        number of polishing iterations [1]
  --asm-coverage INT    reduced coverage for initial assembly [not set]
  --read-error FLOAT    adjust parameters for given read error rate
  --meta                metagenome / uneven coverage mode
  --keep-haplotypes     do not collapse alternative haplotypes
  --mode {fast,normal}  assembly mode
  -v, --version         show version and exit
"""


def _params():
    text = scaffold.from_help_text(FLYE_HELP, "flye")
    return yaml.safe_load(text)["params"]


def test_help_and_version_are_skipped():
    p = _params()
    assert "help" not in p and "version" not in p


def test_bool_flag_inferred():
    p = _params()
    assert p["meta"]["type"] == "bool"
    assert p["meta"]["flag"] == "--meta"
    assert p["keep_haplotypes"]["type"] == "bool"


def test_metavar_types_inferred():
    p = _params()
    assert p["threads"]["type"] == "int"
    assert p["iterations"]["type"] == "int"
    assert p["asm_coverage"]["type"] == "int"
    assert p["read_error"]["type"] == "float"
    assert p["genome_size"]["type"] == "str"
    assert p["nano_raw"]["type"] == "path"


def test_long_flag_is_canonical_over_short():
    p = _params()
    # "-t INT, --threads INT" -> param keyed by the long flag
    assert p["threads"]["flag"] == "--threads"
    assert p["genome_size"]["flag"] == "--genome-size"


def test_numeric_default_extracted():
    p = _params()
    assert p["iterations"].get("default") == 1
    # "[not set]" is not a number -> no default emitted
    assert "default" not in p["asm_coverage"]


def test_enum_choices_parsed():
    p = _params()
    assert p["mode"]["type"] == "enum"
    assert p["mode"]["choices"] == ["fast", "normal"]


def test_wrapped_description_captured():
    p = _params()
    assert "ONT regular reads" in p["nano_raw"]["question"]


def test_everything_is_quiet_by_default():
    p = _params()
    modeled = {k: v for k, v in p.items() if k != "extra_args"}
    assert all(v["ask"] is False and v["tier"] == 3 for v in modeled.values())


def test_extra_args_passthrough_valve_present():
    p = _params()
    assert p["extra_args"]["passthrough"] is True
    assert p["extra_args"]["type"] == "str"
    assert "flag" not in p["extra_args"]


def test_stub_obeys_registry_wellformedness():
    """The scaffold output must satisfy the same invariants the registry test enforces."""
    reg = yaml.safe_load(scaffold.from_help_text(FLYE_HELP, "flye"))
    assert reg["command"] == "flye"
    for name, spec in reg["params"].items():
        assert spec["type"] in {"bool", "int", "float", "str", "path", "enum"}
        if spec["type"] == "enum":
            assert spec.get("choices")
        flagless_ok = spec.get("managed") or spec.get("interpreted") or spec.get("passthrough")
        if not flagless_ok:
            assert spec.get("flag"), f"{name} needs a flag"
