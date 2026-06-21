"""Guards for Snakemake `script:`-directive files — ours and any a user adds.

Snakemake's `script:` directive is a first-class way to attach Python (or R) logic to a rule:
the file is executed with a `snakemake` object injected (input/output/params/wildcards/threads).
It is great for users extending metagx, but it has one sharp edge that is invisible to ordinary
unit tests (which import the function directly): Snakemake **prepends a preamble** to the file at
run time, so:

  * `from __future__ import ...` is no longer the first statement -> SyntaxError, and
  * a top-level `import` of a name the script later rebinds locally can shadow oddly.

We hit exactly the first one (it silently broke the phylogenetics and consensus modules until
they were run end to end). This test scans every file referenced by a `script:` directive and
fails on the `from __future__` antipattern, so neither we nor a contributor can reintroduce it.
It also sanity-checks that each script actually uses the injected `snakemake` object.
"""
import ast
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RULES = REPO / "workflow" / "rules"
SCRIPTS = REPO / "workflow" / "scripts"


def _script_files():
    """Every distinct script referenced by a `script:` directive in the rule files."""
    refs = set()
    for smk in RULES.glob("*.smk"):
        for m in re.finditer(r'script:\s*\n\s*"(\.\./scripts/[^"]+)"', smk.read_text()):
            refs.add((SCRIPTS / Path(m.group(1)).name))
    return sorted(refs)


PY_SCRIPTS = [p for p in _script_files() if p.suffix == ".py"]


def test_some_script_files_are_discovered():
    assert PY_SCRIPTS, "no `script:` python files found — did the rule layout change?"


@pytest.mark.parametrize("script", PY_SCRIPTS, ids=lambda p: p.name)
def test_no_future_import_in_script_directive_files(script):
    """`from __future__` becomes a SyntaxError under Snakemake's `script:` preamble injection."""
    tree = ast.parse(script.read_text(), filename=str(script))
    offenders = [n for n in tree.body
                 if isinstance(n, ast.ImportFrom) and n.module == "__future__"]
    assert not offenders, (
        f"{script.name} has `from __future__ import ...`, which Snakemake's `script:` preamble "
        f"turns into a SyntaxError at run time. Remove it (plain annotations work on py>=3.10)."
    )


@pytest.mark.parametrize("script", PY_SCRIPTS, ids=lambda p: p.name)
def test_script_uses_injected_snakemake_object(script):
    """Each `script:` file should consume the injected `snakemake` object (its whole purpose)."""
    text = script.read_text()
    assert "snakemake" in text, f"{script.name} never references the injected `snakemake` object"
