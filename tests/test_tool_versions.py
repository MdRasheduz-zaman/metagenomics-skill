"""Tool-version floors — because in bioinformatics, tool versions matter.

Bioconda packages drift, and an unrelated install can silently *downgrade* a tool (we hit
exactly this: `mapdamage2`/`abricate` dragged in samtools 0.1.19, whose `sort -o` API breaks
the whole pipeline). A pipeline that "has the tool" but at the wrong version produces wrong or
broken output with no error. This test pins a minimum version for every load-bearing tool and
fails loudly if the environment regresses below it — the floors mirror `environment.yml`.

Each tool is checked only if it is on PATH (so the test is informative locally and in the CI
`e2e` job, and simply skips tools a given machine doesn't have). It also prints the full
version table, which feeds the validation report.
"""
import re

import pytest

from metagx import report
from metagx.doctor import VERSION_FLOORS

# Single source of truth for the floors lives in metagx.doctor (so `metagx doctor` and this
# test can't drift). A tool below its floor is a hard failure, not a warning: it would have
# caught the samtools 0.1.19 regression.
MIN_VERSIONS = VERSION_FLOORS


def _parse_xy(version_str: str):
    """First (major, minor) in a captured version line, or None."""
    m = re.search(r"(\d+)\.(\d+)", version_str or "")
    return (int(m.group(1)), int(m.group(2))) if m else None


# Capture once for the whole module (subprocess calls are the slow part).
_CAPTURED = report.tool_versions(sorted(MIN_VERSIONS))


def test_print_version_table():
    """Always-on: emit the validated tool versions (visible with `pytest -s`)."""
    print("\n=== validated tool versions ===")
    for tool in sorted(_CAPTURED):
        print(f"  {tool:12s} {_CAPTURED[tool]}")


@pytest.mark.parametrize("tool", sorted(MIN_VERSIONS))
def test_tool_meets_version_floor(tool):
    raw = _CAPTURED[tool]
    if raw == "not found on PATH":
        pytest.skip(f"{tool} not on PATH in this shell — activate the bio env "
                    f"(e.g. `conda activate metagx-bio` or append its bin) before running")
    got = _parse_xy(raw)
    if got is None:
        # present but exposes no parseable version (e.g. bracken has no --version) — can't
        # assert a floor, but it's installed. Don't fail; downgrades that DO parse still fail.
        pytest.skip(f"{tool} installed but version not machine-readable: {raw!r}")
    floor = MIN_VERSIONS[tool]
    assert got >= floor, (
        f"{tool} {got[0]}.{got[1]} is below the required floor {floor[0]}.{floor[1]} "
        f"(captured: {raw!r}). An install likely downgraded it — re-pin in environment.yml."
    )
