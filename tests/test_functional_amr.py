"""Functional/AMR module — ABRicate command wiring + a tool-gated real run.

The real-run test executes ABRicate (with its bundled CARD/NCBI/ResFinder/... DBs) on a
known resistance-gene fixture and asserts it is detected — so the functional/AMR path is
verified by execution, not just dry-run. It skips cleanly where ABRicate isn't installed
(CI without --use-conda); run it with the bioconda env (e.g. metagx-amr) on PATH.
"""
import os
import shutil
import subprocess

import pytest

from metagx import registry

_FIX = os.path.join(os.path.dirname(__file__), "fixtures", "amr_ermA.fasta")


def test_abricate_command_renders():
    """The exact tokens the abricate rule builds from the registry."""
    args = registry.render_args("abricate", {"db": "card", "minid": 80.0},
                                managed={"threads": 4})
    assert args == ["--db", "card", "--minid", "80.0", "--threads", "4"]


@pytest.mark.skipif(not (shutil.which("abricate") and shutil.which("blastn")),
                    reason="ABRicate/blastn not installed (skips in CI; run with the AMR conda env)")
def test_abricate_detects_known_resistance_gene(tmp_path):
    """Real execution: ABRicate must detect the ErmA resistance gene in the fixture."""
    args = registry.render_args("abricate", {"db": "card", "minid": 90.0},
                                managed={"threads": 2})
    out = tmp_path / "abricate.tsv"
    # mirror the rule's shell: abricate <args> <contigs> > out
    with open(out, "w") as fh:
        proc = subprocess.run(["abricate"] + args + [_FIX], stdout=fh,
                              stderr=subprocess.PIPE, text=True)
    assert proc.returncode == 0, proc.stderr
    lines = out.read_text().splitlines()
    assert lines and lines[0].startswith("#FILE")          # valid ABRicate TSV header
    hits = [ln for ln in lines[1:] if ln.strip()]
    assert hits, "ABRicate found no AMR gene in the ErmA positive control"
    assert any("rm" in ln for ln in hits)                  # erm(A)/ErmA gene name in a hit row
