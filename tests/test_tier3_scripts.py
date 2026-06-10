"""Unit tests for the Tier 3 Snakemake script helpers (damage authentication + decontam).

Both live under workflow/scripts/ (Snakemake `script:` files), so they are loaded by path.
"""
import importlib.util
from pathlib import Path

_SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


da = _load("damage_authenticate")
dc = _load("decontam")


# ---------------- damage authentication ----------------
def test_first_position_freq(tmp_path):
    f = tmp_path / "5pCtoT_freq.txt"
    f.write_text("pos\t5pC>T\n1\t0.18\n2\t0.09\n3\t0.04\n")
    assert da.first_position_freq(str(f)) == 0.18


def test_authenticate_verdict():
    pos = da.authenticate(0.18, 0.16, 0.05)
    assert pos["damage_present"] is True
    assert "authentic" in pos["verdict"]
    neg = da.authenticate(0.01, 0.00, 0.05)
    assert neg["damage_present"] is False
    assert "no clear" in neg["verdict"]


# ---------------- decontam (prevalence) ----------------
def _rows():
    # taxon "Contaminant sp" appears in the blank; "Real sp" only in real samples
    return [
        {"sample": "blank", "name": "Contaminant sp", "fraction_total_reads": "0.4"},
        {"sample": "s1", "name": "Contaminant sp", "fraction_total_reads": "0.1"},
        {"sample": "s1", "name": "Real sp", "fraction_total_reads": "0.5"},
        {"sample": "s2", "name": "Real sp", "fraction_total_reads": "0.6"},
    ]


def test_flag_contaminants_prevalence():
    flags = dc.flag_contaminants(_rows(), controls={"blank"})
    assert flags["Contaminant sp"]["flagged"] is True    # in the blank, not more prevalent in reals
    assert flags["Real sp"]["flagged"] is False          # never in a control


def test_decontam_run_writes_cleaned(tmp_path):
    import csv
    src = tmp_path / "bracken_combined.tsv"
    with open(src, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["sample", "name", "fraction_total_reads"],
                           delimiter="\t")
        w.writeheader()
        for r in _rows():
            w.writerow(r)
    flagged = tmp_path / "flagged.tsv"
    cleaned = tmp_path / "cleaned.tsv"
    res = dc.run(str(src), ["blank"], str(flagged), str(cleaned))
    assert res["n_flagged"] == 1 and res["flagged"] == ["Contaminant sp"]
    cleaned_rows = list(csv.DictReader(open(cleaned), delimiter="\t"))
    # control rows and the flagged taxon are gone; "Real sp" remains for the real samples
    assert all(r["sample"] != "blank" for r in cleaned_rows)
    assert all(r["name"] != "Contaminant sp" for r in cleaned_rows)
    assert {r["name"] for r in cleaned_rows} == {"Real sp"}
