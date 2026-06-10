"""Differential abundance testing between two sample groups (CLR + permutation).

Diversity (``metagx.diversity``) describes communities; this answers the comparative
question — *which taxa differ between conditions?* — that a study usually exists to answer.

Approach (an ALDEx2-lite, dependency-light): centred-log-ratio transform the count matrix
(compositional-data correct), then for each taxon compare the two groups' mean CLR with a
nonparametric **permutation test** (no normality assumption) and control the false-discovery
rate with Benjamini-Hochberg. Pure numpy — no scipy/R — so it is unit-testable and runs
anywhere the core install runs.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from .diversity import clr


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """BH step-up FDR correction. Returns q-values aligned to the input order."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n, dtype=float)
    q[order] = np.clip(ranked, 0, 1)
    return q


def _group_diff(clr_mat: np.ndarray, is_a: np.ndarray) -> np.ndarray:
    """Mean-CLR difference per taxon: group A minus group B."""
    return clr_mat[is_a].mean(axis=0) - clr_mat[~is_a].mean(axis=0)


def permutation_pvalues(clr_mat: np.ndarray, is_a: np.ndarray, n_perm: int,
                        seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
    """Two-sided permutation p-values for the per-taxon group difference.

    Returns ``(observed_diff, pvalues)``. Labels are shuffled ``n_perm`` times; the p-value
    is the fraction of permutations whose |difference| meets or exceeds the observed (with the
    standard +1 so a p-value is never exactly 0).
    """
    rng = np.random.default_rng(seed)
    obs = _group_diff(clr_mat, is_a)
    n = clr_mat.shape[0]
    n_a = int(is_a.sum())
    ge = np.zeros(clr_mat.shape[1], dtype=int)
    abs_obs = np.abs(obs)
    for _ in range(n_perm):
        perm = rng.permutation(n)
        mask = np.zeros(n, dtype=bool)
        mask[perm[:n_a]] = True
        ge += (np.abs(_group_diff(clr_mat, mask)) >= abs_obs).astype(int)
    pvals = (ge + 1) / (n_perm + 1)
    return obs, pvals


def differential_abundance(
    samples: List[str], taxa: List[str], mat: np.ndarray,
    sample_group: Dict[str, str], n_perm: int = 999, fdr: float = 0.05,
    seed: int = 42, reference_group: Optional[str] = None,
) -> Tuple[List[Dict], Dict]:
    """Per-taxon differential abundance between exactly two groups.

    ``sample_group`` maps sample name -> group label; samples without a label are dropped.
    Returns ``(rows, summary)`` where rows carry mean CLR per group, the difference, a
    standardized effect size, raw p and BH q, and a significance flag.
    """
    keep = [i for i, s in enumerate(samples) if sample_group.get(s)]
    if not keep:
        raise ValueError("no samples carry a group label")
    grp = [sample_group[samples[i]] for i in keep]
    groups = sorted(set(grp))
    if len(groups) != 2:
        raise ValueError(f"differential abundance needs exactly 2 groups, got {groups}")
    if reference_group and reference_group in groups:
        b = reference_group
        a = [g for g in groups if g != b][0]
    else:
        a, b = groups  # sorted: difference is (first - second)
    sub = mat[keep, :]
    clr_mat = clr(sub)
    is_a = np.array([g == a for g in grp], dtype=bool)
    if is_a.sum() < 2 or (~is_a).sum() < 2:
        raise ValueError("each group needs >=2 samples for the permutation test")

    obs, pvals = permutation_pvalues(clr_mat, is_a, n_perm, seed)
    qvals = benjamini_hochberg(pvals)
    mean_a = clr_mat[is_a].mean(axis=0)
    mean_b = clr_mat[~is_a].mean(axis=0)
    var_a = clr_mat[is_a].var(axis=0, ddof=1)
    var_b = clr_mat[~is_a].var(axis=0, ddof=1)
    pooled = np.sqrt((var_a + var_b) / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        effect = np.where(pooled > 0, obs / pooled, 0.0)

    rows: List[Dict] = []
    for t in range(len(taxa)):
        rows.append({
            "taxon": taxa[t],
            f"mean_clr_{a}": round(float(mean_a[t]), 4),
            f"mean_clr_{b}": round(float(mean_b[t]), 4),
            "clr_diff": round(float(obs[t]), 4),
            "effect_size": round(float(effect[t]), 4),
            "p_value": round(float(pvals[t]), 5),
            "q_value": round(float(qvals[t]), 5),
            "significant": bool(qvals[t] < fdr),
        })
    rows.sort(key=lambda r: (r["q_value"], r["p_value"]))
    summary = {
        "group_a": a, "group_b": b,
        "n_a": int(is_a.sum()), "n_b": int((~is_a).sum()),
        "n_taxa": len(taxa), "n_permutations": n_perm, "fdr": fdr,
        "n_significant": int(sum(r["significant"] for r in rows)),
        "diff_definition": f"clr_diff = mean_clr_{a} - mean_clr_{b}",
    }
    return rows, summary
