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
