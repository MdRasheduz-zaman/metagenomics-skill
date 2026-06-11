"""Unit tests for differential abundance (CLR + permutation + BH FDR)."""
import numpy as np
import pytest

from metagx import differential as da


def test_benjamini_hochberg_monotone_and_bounded():
    q = da.benjamini_hochberg(np.array([0.01, 0.02, 0.03, 0.5]))
    assert np.all((q >= 0) & (q <= 1))
    # q-values are monotone in the same order as p-values
    assert q[0] <= q[1] <= q[2] <= q[3]
    # an all-ones input stays at one
    assert np.allclose(da.benjamini_hochberg(np.ones(5)), 1.0)


def test_detects_a_planted_signal_with_replication():
    # taxon0 strongly enriched in group A; 6 vs 6 gives the permutation test power.
    rng = np.random.default_rng(0)
    samples = [f"a{i}" for i in range(6)] + [f"b{i}" for i in range(6)]
    groups = {s: ("A" if s.startswith("a") else "B") for s in samples}
    base = rng.integers(40, 60, size=(12, 4)).astype(float)
    base[:6, 0] += 400          # taxon0 way up in group A
    rows, summary = da.differential_abundance(
        samples, ["t0", "t1", "t2", "t3"], base, groups, n_perm=999, fdr=0.1, seed=1)
    by_taxon = {r["taxon"]: r for r in rows}
    assert by_taxon["t0"]["significant"] is True
    assert summary["n_significant"] >= 1
    assert by_taxon["t0"]["clr_diff"] > 0          # enriched in A (group_a)


def test_random_groups_control_false_positives():
    # no real difference -> expect ~no significant taxa (FDR control)
    rng = np.random.default_rng(2)
    samples = [f"s{i}" for i in range(10)]
    groups = {s: ("A" if i < 5 else "B") for i, s in enumerate(samples)}
    mat = rng.integers(10, 100, size=(10, 8)).astype(float)
    _, summary = da.differential_abundance(samples, [f"t{i}" for i in range(8)],
                                           mat, groups, n_perm=999, fdr=0.05, seed=3)
    assert summary["n_significant"] == 0


def test_requires_exactly_two_groups():
    samples = ["a", "b", "c"]
    mat = np.ones((3, 2))
    with pytest.raises(ValueError):
        da.differential_abundance(samples, ["t0", "t1"], mat,
                                  {"a": "A", "b": "B", "c": "C"})


def test_reference_group_sets_direction():
    samples = [f"a{i}" for i in range(4)] + [f"b{i}" for i in range(4)]
    groups = {s: ("hi" if s.startswith("a") else "lo") for s in samples}
    mat = np.ones((8, 2)) * 10
    mat[:4, 0] = 100  # taxon0 high in 'hi'
    rows, summary = da.differential_abundance(
        samples, ["t0", "t1"], mat, groups, n_perm=199, reference_group="lo")
    assert summary["group_b"] == "lo"          # reference is the subtrahend
    assert next(r for r in rows if r["taxon"] == "t0")["clr_diff"] > 0


def test_dirichlet_clr_instances_shape_and_centering():
    rng = np.random.default_rng(0)
    counts = np.array([[100.0, 50.0, 1.0], [10.0, 80.0, 60.0]])
    inst = da.dirichlet_clr_instances(counts, mc_samples=16, rng=rng)
    assert inst.shape == (16, 2, 3)
    # CLR rows are mean-centred (sum to ~0) by construction
    assert np.allclose(inst.sum(axis=2), 0.0, atol=1e-9)


def test_mc_detects_planted_signal_and_records_method():
    # same planted signal as the point-estimate test, now through Dirichlet Monte-Carlo.
    rng = np.random.default_rng(0)
    samples = [f"a{i}" for i in range(6)] + [f"b{i}" for i in range(6)]
    groups = {s: ("A" if s.startswith("a") else "B") for s in samples}
    base = rng.integers(40, 60, size=(12, 4)).astype(float)
    base[:6, 0] += 400
    rows, summary = da.differential_abundance(
        samples, ["t0", "t1", "t2", "t3"], base, groups,
        n_perm=499, fdr=0.1, seed=1, mc_samples=32)
    by_taxon = {r["taxon"]: r for r in rows}
    assert by_taxon["t0"]["significant"] is True
    assert by_taxon["t0"]["clr_diff"] > 0
    assert summary["mc_samples"] == 32
    assert "Dirichlet" in summary["method"]


def test_mc_controls_false_positives_on_null():
    rng = np.random.default_rng(2)
    samples = [f"s{i}" for i in range(10)]
    groups = {s: ("A" if i < 5 else "B") for i, s in enumerate(samples)}
    mat = rng.integers(10, 100, size=(10, 8)).astype(float)
    _, summary = da.differential_abundance(
        samples, [f"t{i}" for i in range(8)], mat, groups,
        n_perm=499, fdr=0.05, seed=3, mc_samples=32)
    assert summary["n_significant"] == 0


def test_mc_one_falls_back_to_point_estimate():
    # mc_samples<=1 must reproduce the legacy single-CLR permutation path exactly.
    rng = np.random.default_rng(5)
    samples = [f"a{i}" for i in range(4)] + [f"b{i}" for i in range(4)]
    groups = {s: ("A" if s.startswith("a") else "B") for s in samples}
    mat = rng.integers(20, 80, size=(8, 5)).astype(float)
    rows1, summary1 = da.differential_abundance(
        samples, [f"t{i}" for i in range(5)], mat, groups, n_perm=299, seed=7, mc_samples=1)
    # deterministic point-estimate numbers match a hand-run of the legacy formula
    clr_mat = da.clr(mat)
    is_a = np.array([True] * 4 + [False] * 4)
    obs, pvals = da.permutation_pvalues(clr_mat, is_a, 299, seed=7)
    assert summary1["mc_samples"] == 1
    assert "point estimate" in summary1["method"]
    by_taxon = {r["taxon"]: r for r in rows1}
    assert by_taxon["t0"]["clr_diff"] == round(float(obs[0]), 4)
    assert by_taxon["t0"]["p_value"] == round(float(pvals[0]), 5)
