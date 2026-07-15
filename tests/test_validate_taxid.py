"""RT-3: taxonomy roll-up + taxid-based rank agreement in the validate module.

Genus-level validation was empty/unreliable on non-binomial taxa (a genus name and its species
name share no tokens). These lock in rolling a leaf assignment / BLAST subject up to the requested
rank via nodes.dmp and comparing taxids.
"""
from metagx import validation

_NODES = (
    "11053\t|\t12637\t|\tspecies\t|\t\n"        # Dengue virus 1 (species)
    "12637\t|\t3044782\t|\tno rank\t|\t\n"      # Denguevirus grouping (no rank)
    "3044782\t|\t11050\t|\tgenus\t|\t\n"        # Orthoflavivirus (genus)
    "11050\t|\t1\t|\tfamily\t|\t\n"             # Flaviviridae (family)
    "1\t|\t1\t|\tno rank\t|\t\n"
)


def _tree(tmp_path):
    p = tmp_path / "nodes.dmp"
    p.write_text(_NODES)
    return validation.parse_nodes_dmp(str(p))


def test_ancestor_at_rank_rolls_species_to_genus(tmp_path):
    tree = _tree(tmp_path)
    assert validation.ancestor_at_rank("11053", "genus", tree) == "3044782"
    assert validation.ancestor_at_rank("11053", "species", tree) == "11053"
    assert validation.ancestor_at_rank("11053", "family", tree) == "11050"
    assert validation.ancestor_at_rank("11053", "kingdom", tree) is None  # not in this lineage


def test_assess_agrees_by_taxid_when_names_differ(tmp_path):
    # classifier says genus Orthoflavivirus (3044782); BLAST subject is a species (Dengue, 11053)
    # whose name shares no token with the genus — name matching would FAIL, taxid roll-up AGREES.
    tree = _tree(tmp_path)
    q = {"r1": "Orthoflavivirus"}
    hits = {"r1": {"qseqid": "r1", "staxids": "11053", "sscinames": "Dengue virus 1",
                   "bitscore": 100.0}}
    a = validation.assess(q, hits, level="genus",
                          query_taxids={"r1": "3044782"}, tree=tree)
    assert a["n_agree"] == 1 and a["agreement_rate"] == 1.0
    # without the tree, the name-token comparison disagrees (the pre-RT-3 behavior)
    a2 = validation.assess(q, hits, level="genus")
    assert a2["n_agree"] == 0


def test_assess_taxid_disagrees_on_wrong_genus(tmp_path):
    tree = _tree(tmp_path)
    q = {"r1": "Alphavirus"}
    hits = {"r1": {"qseqid": "r1", "staxids": "11053", "bitscore": 100.0}}  # rolls to 3044782
    a = validation.assess(q, hits, level="genus", query_taxids={"r1": "999999"}, tree=tree)
    assert a["n_agree"] == 0
