"""Aggregate module (MultiQC + Krona) — command construction + kreport2krona conversion.

MultiQC is verified by real execution (ROADMAP); here we pin the command tokens both rules
build, and unit-test the pure-Python kreport2krona converter (which replaces ktImportTaxonomy
so Krona works with custom kraken2 DBs — no taxonomy database needed).
"""
import importlib.util
import pathlib

from metagx import registry

_SPEC = importlib.util.spec_from_file_location(
    "kreport2krona",
    pathlib.Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "kreport2krona.py",
)
k2k = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(k2k)


def test_multiqc_command_tokens():
    args = registry.render_args(
        "multiqc", {},
        managed={"outdir": "out/report/multiqc", "filename": "multiqc_report", "force": True})
    assert args == ["-o", "out/report/multiqc", "-n", "multiqc_report", "-f"]


def test_krona_uses_ktimporttext_output_only():
    # ktImportText needs only -o; the lineage comes from the kreport2krona text, not a tax DB.
    assert registry.load_registry("krona")["command"] == "ktImportText"
    args = registry.render_args("krona", {}, managed={"output": "out/report/krona.html"})
    assert args == ["-o", "out/report/krona.html"]


def test_kreport2krona_builds_lineage_from_indentation(tmp_path):
    # 2 spaces per rank in the name column => indentation is the lineage. Custom-DB taxids fine.
    kr = tmp_path / "s.kreport"
    kr.write_text(
        " 50.0\t5000\t5000\tU\t0\tunclassified\n"
        " 49.0\t4900\t0\tR\t1\troot\n"
        "  4.5\t450\t450\tS\t1008\t  Powassan virus, complete genome\n"
        "  3.7\t370\t370\tS\t1018\t  Acute bee paralysis virus, complete genome\n"
    )
    out = tmp_path / "s.krona.txt"
    n = k2k.convert(str(kr), str(out))
    lines = [ln.rstrip("\n").split("\t") for ln in open(out)]
    assert n == 2                                  # unclassified + 0-read root dropped
    assert ["450", "root", "Powassan virus, complete genome"] in lines
    assert ["370", "root", "Acute bee paralysis virus, complete genome"] in lines
    assert not any("unclassified" in row for row in lines)
