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
