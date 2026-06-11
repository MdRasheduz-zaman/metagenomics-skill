"""Differential abundance testing between two sample groups (ALDEx2-style).

Diversity (``metagx.diversity``) describes communities; this answers the comparative
question — *which taxa differ between conditions?* — that a study usually exists to answer.

Approach (an ALDEx2-style test, dependency-light): a Bracken count vector is one multinomial
draw from an unknown composition, so a single CLR point estimate hides the sampling-depth
uncertainty. Following ALDEx2, we model each sample's composition with a **Dirichlet posterior**
``Dir(counts + 0.5)`` (Jeffreys prior), draw ``mc_samples`` Monte-Carlo instances of the
proportions, centred-log-ratio transform each, run a nonparametric **permutation test** per
instance, and report the **expected** p-value / CLR difference across instances (the median for
the effect size). The expected p-values are then FDR-controlled with Benjamini-Hochberg. Setting
``mc_samples <= 1`` recovers the old single point-estimate behaviour (fixed-pseudocount CLR).
Pure numpy — no scipy/R — so it is unit-testable and runs anywhere the core install runs.
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


def _standardized_effect(clr_mat: np.ndarray, is_a: np.ndarray, obs: np.ndarray) -> np.ndarray:
    """Per-taxon effect size: group-difference / pooled within-group SD."""
    var_a = clr_mat[is_a].var(axis=0, ddof=1)
    var_b = clr_mat[~is_a].var(axis=0, ddof=1)
    pooled = np.sqrt((var_a + var_b) / 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(pooled > 0, obs / pooled, 0.0)


def dirichlet_clr_instances(counts: np.ndarray, mc_samples: int, rng: np.random.Generator,
                            pseudocount: float = 0.5) -> np.ndarray:
    """Monte-Carlo CLR instances drawn from each sample's Dirichlet posterior.

    Each row of ``counts`` is treated as a multinomial draw; its composition posterior is
    ``Dir(counts + pseudocount)``. Returns an array of shape ``(mc_samples, n_samples, n_taxa)``
    where every ``[m]`` slice is one fully-resampled, CLR-transformed count matrix.
    """
    n_s, n_t = counts.shape
    inst = np.empty((mc_samples, n_s, n_t), dtype=float)
    alpha = counts + pseudocount
    for i in range(n_s):
        props = rng.dirichlet(alpha[i], size=mc_samples)  # (mc_samples, n_taxa), rows sum to 1
        logp = np.log(props)
        inst[:, i, :] = logp - logp.mean(axis=1, keepdims=True)
    return inst


def aldex_monte_carlo(counts: np.ndarray, is_a: np.ndarray, n_perm: int, mc_samples: int,
                      seed: int = 42):
    """ALDEx2-style test: permutation test over Dirichlet Monte-Carlo CLR instances.

    Returns ``(exp_diff, exp_p, med_effect, mean_a, mean_b)`` aggregated across instances —
    the expected (mean) CLR difference and p-value, the median standardized effect, and the
    expected per-group mean CLR. Each instance uses its own permutation seed for independence.
    """
    rng = np.random.default_rng(seed)
    inst = dirichlet_clr_instances(counts, mc_samples, rng)
    n_t = counts.shape[1]
    P = np.empty((mc_samples, n_t)); D = np.empty((mc_samples, n_t))
    E = np.empty((mc_samples, n_t)); MA = np.empty((mc_samples, n_t)); MB = np.empty((mc_samples, n_t))
    for m in range(mc_samples):
        clr_mat = inst[m]
        obs, pvals = permutation_pvalues(clr_mat, is_a, n_perm, seed=seed + m + 1)
        D[m] = obs
        P[m] = pvals
        E[m] = _standardized_effect(clr_mat, is_a, obs)
        MA[m] = clr_mat[is_a].mean(axis=0)
        MB[m] = clr_mat[~is_a].mean(axis=0)
    return D.mean(0), P.mean(0), np.median(E, axis=0), MA.mean(0), MB.mean(0)


def differential_abundance(
    samples: List[str], taxa: List[str], mat: np.ndarray,
    sample_group: Dict[str, str], n_perm: int = 999, fdr: float = 0.05,
    seed: int = 42, reference_group: Optional[str] = None, mc_samples: int = 128,
) -> Tuple[List[Dict], Dict]:
    """Per-taxon differential abundance between exactly two groups.

    ``sample_group`` maps sample name -> group label; samples without a label are dropped.
    With ``mc_samples > 1`` (default) the test propagates sampling-depth uncertainty via
    ALDEx2-style Dirichlet Monte-Carlo (see module docstring); ``mc_samples <= 1`` falls back
    to a single fixed-pseudocount CLR point estimate. Returns ``(rows, summary)`` where rows
    carry mean CLR per group, the difference, a standardized effect size, p and BH q, and a
    significance flag.
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
    is_a = np.array([g == a for g in grp], dtype=bool)
    if is_a.sum() < 2 or (~is_a).sum() < 2:
        raise ValueError("each group needs >=2 samples for the permutation test")

    if mc_samples and mc_samples > 1:
        obs, pvals, effect, mean_a, mean_b = aldex_monte_carlo(
            sub, is_a, n_perm, int(mc_samples), seed)
        method = f"ALDEx2-style Dirichlet Monte-Carlo (mc={int(mc_samples)}) + CLR + permutation + BH-FDR"
    else:
        clr_mat = clr(sub)
        obs, pvals = permutation_pvalues(clr_mat, is_a, n_perm, seed)
        effect = _standardized_effect(clr_mat, is_a, obs)
        mean_a = clr_mat[is_a].mean(axis=0)
        mean_b = clr_mat[~is_a].mean(axis=0)
        method = "CLR point estimate (fixed pseudocount) + permutation + BH-FDR"
    qvals = benjamini_hochberg(pvals)

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
        "mc_samples": int(mc_samples) if (mc_samples and mc_samples > 1) else 1,
        "method": method,
        "n_significant": int(sum(r["significant"] for r in rows)),
        "diff_definition": f"clr_diff = mean_clr_{a} - mean_clr_{b}",
    }
    return rows, summary
