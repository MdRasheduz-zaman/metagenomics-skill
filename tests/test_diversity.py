import math
import numpy as np
from metagx import diversity as dv


def test_build_matrix():
    rows = [{"sample": "s1", "name": "A", "new_est_reads": "10"},
            {"sample": "s1", "name": "B", "new_est_reads": "30"},
            {"sample": "s2", "name": "A", "new_est_reads": "20"}]
    samples, taxa, mat = dv.build_matrix(rows)
    assert samples == ["s1", "s2"] and set(taxa) == {"A", "B"}
    assert mat.sum() == 60


def test_alpha_indices_known_values():
    even = np.array([1.0, 1.0, 1.0, 1.0])
    assert dv.richness(even) == 4
    assert math.isclose(dv.shannon(even), math.log(4), rel_tol=1e-9)
    assert math.isclose(dv.simpson(even), 0.75, rel_tol=1e-9)       # 1 - 4*(1/4)^2
    assert math.isclose(dv.pielou(even), 1.0, rel_tol=1e-9)         # perfectly even
    assert dv.shannon(np.array([10.0, 0, 0])) == 0.0                # single taxon


def test_relative_abundance_sums_to_one():
    mat = np.array([[1.0, 3.0], [2.0, 2.0]])
    rel = dv.relative_abundance(mat)
    assert np.allclose(rel.sum(axis=1), [1.0, 1.0])


def test_braycurtis_identity_and_disjoint():
    mat = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 0.0]])
    d = dv.braycurtis(mat)
    assert d[0, 0] == 0.0
    assert math.isclose(d[0, 1], 1.0)      # disjoint communities -> 1
    assert d[0, 2] == 0.0                   # identical


def test_pcoa_shapes():
    d = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]], dtype=float)
    coords, expl = dv.pcoa(d, n_axes=2)
    assert coords.shape == (3, 2)
    assert len(expl) == 2 and expl[0] >= expl[1]


def test_chao1_bias_corrected_known_values():
    # counts [1,1,2,5]: S_obs=4, f1=2, f2=1 -> 4 + 2*1/(2*(1+1)) = 4.5
    assert math.isclose(dv.chao1(np.array([1.0, 1, 2, 5])), 4.5, rel_tol=1e-9)
    # no doubletons (f2=0) stays defined: 4 + 3*2/(2*1) = 7
    assert math.isclose(dv.chao1(np.array([1.0, 1, 1, 5])), 7.0, rel_tol=1e-9)
    # no singletons -> equals observed richness
    assert math.isclose(dv.chao1(np.array([5.0, 5, 5])), 3.0, rel_tol=1e-9)


def test_ace_known_values():
    # all abundant (>10): ACE == observed richness
    assert math.isclose(dv.ace(np.array([20.0, 30, 40])), 3.0, rel_tol=1e-9)
    # mixed case computed by hand -> ~7.6122
    assert math.isclose(dv.ace(np.array([1.0, 1, 2, 5, 20])), 7.612245, rel_tol=1e-5)
    # degenerate coverage (all rare are singletons) falls back to Chao1
    assert math.isclose(dv.ace(np.array([1.0, 1, 1])), dv.chao1(np.array([1.0, 1, 1])))


def test_goods_coverage():
    # N=9, f1=2 -> 1 - 2/9
    assert math.isclose(dv.goods_coverage(np.array([1.0, 1, 2, 5])), 1 - 2 / 9, rel_tol=1e-9)
    # no singletons -> full coverage
    assert math.isclose(dv.goods_coverage(np.array([3.0, 4])), 1.0, rel_tol=1e-9)


def test_rarefaction_curve_invariants():
    curve = dict(dv.rarefaction_curve(np.array([2.0, 3.0]), depths=[1, 2, 5]))
    assert math.isclose(curve[1], 1.0, rel_tol=1e-9)   # one read -> exactly 1 taxon
    assert math.isclose(curve[2], 1.6, rel_tol=1e-9)   # hand-computed Hurlbert
    assert math.isclose(curve[5], 2.0, rel_tol=1e-9)   # full depth -> S_obs
    # monotonic non-decreasing
    ys = [y for _, y in dv.rarefaction_curve(np.array([5.0, 3, 2, 1]))]
    assert all(b >= a - 1e-9 for a, b in zip(ys, ys[1:]))


def test_jaccard_and_sorensen_presence_absence():
    mat = np.array([[1.0, 1, 0], [1, 0, 1]])
    assert math.isclose(dv.jaccard(mat)[0, 1], 1 - 1 / 3, rel_tol=1e-9)
    assert math.isclose(dv.sorensen(mat)[0, 1], 0.5, rel_tol=1e-9)
    # identical communities -> 0 dissimilarity
    same = np.array([[1.0, 2, 3], [4, 5, 6]])
    assert dv.jaccard(same)[0, 1] == 0.0


def test_core_taxa_prevalence():
    mat = np.array([[1.0, 1, 0], [1, 0, 1], [1, 1, 1]])
    taxa = ["A", "B", "C"]
    samples = ["s1", "s2", "s3"]
    core = dv.core_taxa(samples, taxa, mat, prevalence=0.8)
    assert [c["taxon"] for c in core] == ["A"]           # only A is in all 3
    assert math.isclose(core[0]["prevalence"], 1.0)
    loose = dv.core_taxa(samples, taxa, mat, prevalence=0.6)
    assert {c["taxon"] for c in loose} == {"A", "B", "C"}
    assert loose[0]["taxon"] == "A"                       # sorted by prevalence desc
