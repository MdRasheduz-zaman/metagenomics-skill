"""Cross-sample community statistics — normalization, alpha/beta diversity, ordination,
richness estimation, rarefaction, and core-microbiome analysis.

Pure-Python (numpy only) so it runs anywhere and is unit-testable. Operates on a
sample x taxon count matrix built from the combined Bracken table. Covers:
  * normalization: relative abundance (TSS), CLR
  * alpha: Shannon, Gini-Simpson, observed richness, Pielou evenness, Chao1 & ACE
    (asymptotic richness estimators), Good's coverage (sampling completeness)
  * rarefaction: analytic Hurlbert expected richness (no random subsampling)
  * beta: Bray-Curtis (abundance), Jaccard & Sorensen (presence/absence)
  * ordination: classical PCoA
  * core microbiome: taxa shared across a prevalence threshold of samples
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np


def build_matrix(rows: List[Dict], value_key: str = "new_est_reads",
                 taxon_key: str = "name", sample_key: str = "sample"):
    """rows -> (samples, taxa, matrix[n_samples, n_taxa]) of counts."""
    samples, taxa = [], []
    s_index, t_index = {}, {}
    for r in rows:
        s, t = r[sample_key], r[taxon_key]
        if s not in s_index:
            s_index[s] = len(samples); samples.append(s)
        if t not in t_index:
            t_index[t] = len(taxa); taxa.append(t)
    mat = np.zeros((len(samples), len(taxa)), dtype=float)
    for r in rows:
        try:
            mat[s_index[r[sample_key]], t_index[r[taxon_key]]] += float(r[value_key])
        except (ValueError, TypeError):
            pass
    return samples, taxa, mat


def relative_abundance(mat: np.ndarray) -> np.ndarray:
    """Total-sum scaling (each sample sums to 1)."""
    tot = mat.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    return mat / tot


def clr(mat: np.ndarray, pseudocount: float = 0.5) -> np.ndarray:
    """Centered log-ratio transform (compositional)."""
    m = mat + pseudocount
    logm = np.log(m)
    return logm - logm.mean(axis=1, keepdims=True)


def shannon(vec: np.ndarray) -> float:
    p = vec[vec > 0]
    s = p.sum()
    if s <= 0:
        return 0.0
    p = p / s
    return float(-(p * np.log(p)).sum())


def simpson(vec: np.ndarray) -> float:
    """Gini-Simpson index (1 - sum p^2)."""
    p = vec[vec > 0]
    s = p.sum()
    if s <= 0:
        return 0.0
    p = p / s
    return float(1.0 - (p * p).sum())


def richness(vec: np.ndarray) -> int:
    return int((vec > 0).sum())


def pielou(vec: np.ndarray) -> float:
    r = richness(vec)
    return float(shannon(vec) / math.log(r)) if r > 1 else 0.0


# --------------------------------------------------------------------------- #
# Asymptotic richness estimators + sampling completeness                       #
# These answer "how many taxa are really there (incl. unseen)?" and "did I     #
# sequence deeply enough?" — count-based, so they round est. reads to integers.#
# --------------------------------------------------------------------------- #
def _int_counts(vec: np.ndarray) -> np.ndarray:
    """Non-negative integer counts for the present taxa (rounds est. reads)."""
    c = np.rint(np.asarray(vec, dtype=float)).astype(int)
    return c[c > 0]


def _freq_of_freqs(counts: np.ndarray) -> Dict[int, int]:
    """f_i = number of taxa observed exactly i times."""
    fof: Dict[int, int] = {}
    for c in counts:
        fof[int(c)] = fof.get(int(c), 0) + 1
    return fof


def chao1(vec: np.ndarray) -> float:
    """Bias-corrected Chao1 estimate of true richness (uses singletons/doubletons).

    S_obs + f1(f1-1) / (2(f2+1)). The bias-corrected form is defined even when
    there are no doubletons (f2=0), unlike the classic f1^2/(2 f2).
    """
    counts = _int_counts(vec)
    s_obs = counts.size
    if s_obs == 0:
        return 0.0
    fof = _freq_of_freqs(counts)
    f1, f2 = fof.get(1, 0), fof.get(2, 0)
    return float(s_obs + (f1 * (f1 - 1)) / (2 * (f2 + 1)))


def ace(vec: np.ndarray, rare_threshold: int = 10) -> float:
    """Abundance-based Coverage Estimator of true richness.

    Splits taxa into 'abundant' (>threshold) and 'rare' (<=threshold), then corrects
    the rare count by the estimated sample coverage. Falls back to Chao1 when the
    coverage estimate is degenerate (all rare taxa are singletons).
    """
    counts = _int_counts(vec)
    if counts.size == 0:
        return 0.0
    fof = _freq_of_freqs(counts)
    s_abund = int((counts > rare_threshold).sum())
    rare = counts[counts <= rare_threshold]
    s_rare = int(rare.size)
    if s_rare == 0:
        return float(s_abund)
    n_rare = int(rare.sum())
    f1 = fof.get(1, 0)
    c_ace = 1.0 - (f1 / n_rare) if n_rare > 0 else 0.0
    if c_ace <= 0 or n_rare <= 1:
        return chao1(vec)            # coverage estimate degenerate
    summ = sum(i * (i - 1) * fof.get(i, 0) for i in range(1, rare_threshold + 1))
    gamma2 = max((s_rare / c_ace) * summ / (n_rare * (n_rare - 1)) - 1.0, 0.0)
    return float(s_abund + s_rare / c_ace + (f1 / c_ace) * gamma2)


def goods_coverage(vec: np.ndarray) -> float:
    """Good's coverage: estimated fraction of the community that was observed (1 - f1/N)."""
    counts = _int_counts(vec)
    n = int(counts.sum())
    if n == 0:
        return 0.0
    f1 = _freq_of_freqs(counts).get(1, 0)
    return float(1.0 - f1 / n)


def rarefaction_curve(vec: np.ndarray, depths: List[int] | None = None,
                      n_points: int = 10) -> List[Tuple[int, float]]:
    """Analytic rarefaction (Hurlbert's expected richness) — no random subsampling.

    Returns [(depth, expected_richness), ...]. E[S_m] = sum_i (1 - C(N-n_i, m)/C(N, m)),
    computed with log-gamma for numerical stability. The curve flattening toward S_obs
    means the community was sampled deeply enough.
    """
    counts = _int_counts(vec)
    n = int(counts.sum())
    s_obs = counts.size
    if n == 0 or s_obs == 0:
        return []
    if depths is None:
        step = max(1, n // n_points)
        depths = list(range(step, n, step)) + [n]
    depths = [int(m) for m in depths if 0 < int(m) <= n]
    lg = math.lgamma
    log_choose_N = lg(n + 1)
    out: List[Tuple[int, float]] = []
    for m in depths:
        log_denom = log_choose_N - lg(m + 1) - lg(n - m + 1)
        exp_s = 0.0
        for ni in counts:
            rem = n - int(ni)
            if rem >= m:
                log_num = lg(rem + 1) - lg(m + 1) - lg(rem - m + 1)
                exp_s += 1.0 - math.exp(log_num - log_denom)
            else:
                exp_s += 1.0      # taxon cannot be missed at this depth
        out.append((m, round(exp_s, 4)))
    return out


def alpha_table(samples: List[str], mat: np.ndarray) -> List[Dict]:
    out = []
    for i, s in enumerate(samples):
        v = mat[i]
        out.append({"sample": s, "richness": richness(v),
                    "shannon": round(shannon(v), 4), "simpson": round(simpson(v), 4),
                    "pielou_evenness": round(pielou(v), 4),
                    "chao1": round(chao1(v), 2), "ace": round(ace(v), 2),
                    "goods_coverage": round(goods_coverage(v), 4)})
    return out


# --------------------------------------------------------------------------- #
# Presence/absence beta diversity + core microbiome                            #
# --------------------------------------------------------------------------- #
def _presence(mat: np.ndarray, detection: float = 0.0) -> np.ndarray:
    return mat > detection


def jaccard(mat: np.ndarray, detection: float = 0.0) -> np.ndarray:
    """Pairwise Jaccard dissimilarity (presence/absence) — complements Bray-Curtis."""
    pres = _presence(mat, detection)
    n = pres.shape[0]
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a = int((pres[i] & pres[j]).sum())
            union = int((pres[i] | pres[j]).sum())
            val = 1.0 - a / union if union > 0 else 0.0
            d[i, j] = d[j, i] = val
    return d


def sorensen(mat: np.ndarray, detection: float = 0.0) -> np.ndarray:
    """Pairwise Sørensen–Dice dissimilarity (presence/absence)."""
    pres = _presence(mat, detection)
    n = pres.shape[0]
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            a = int((pres[i] & pres[j]).sum())
            denom = int(pres[i].sum() + pres[j].sum())
            val = 1.0 - (2 * a) / denom if denom > 0 else 0.0
            d[i, j] = d[j, i] = val
    return d


def core_taxa(samples: List[str], taxa: List[str], mat: np.ndarray,
              prevalence: float = 0.8, detection: float = 0.0) -> List[Dict]:
    """Taxa present (count > detection) in >= `prevalence` fraction of samples.

    Returns [{taxon, prevalence, mean_rel_abundance}, ...] sorted by prevalence then
    mean abundance — the stable 'core' community shared across samples.
    """
    if not samples or not taxa:
        return []
    pres = _presence(mat, detection)
    rel = relative_abundance(mat)
    n = len(samples)
    out = []
    for t in range(len(taxa)):
        prev = float(pres[:, t].sum()) / n
        if prev >= prevalence:
            out.append({"taxon": taxa[t], "prevalence": round(prev, 4),
                        "mean_rel_abundance": round(float(rel[:, t].mean()), 6)})
    out.sort(key=lambda d: (-d["prevalence"], -d["mean_rel_abundance"]))
    return out


def braycurtis(mat: np.ndarray) -> np.ndarray:
    """Pairwise Bray-Curtis dissimilarity on the (count) matrix."""
    n = mat.shape[0]
    d = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(mat[i] - mat[j]).sum()
            den = (mat[i].sum() + mat[j].sum())
            val = num / den if den > 0 else 0.0
            d[i, j] = d[j, i] = val
    return d


def pcoa(dist: np.ndarray, n_axes: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """Classical (Torgerson) PCoA. Returns (coords[n, n_axes], explained_fraction)."""
    n = dist.shape[0]
    if n < 2:
        return np.zeros((n, n_axes)), np.zeros(n_axes)
    d2 = dist ** 2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j.dot(d2).dot(j)
    vals, vecs = np.linalg.eigh(b)
    idx = np.argsort(vals)[::-1]
    vals, vecs = vals[idx], vecs[:, idx]
    pos = np.clip(vals, 0, None)
    k = min(n_axes, n)
    coords = vecs[:, :k] * np.sqrt(pos[:k])
    total = pos.sum() or 1.0
    explained = (pos[:k] / total)
    if coords.shape[1] < n_axes:  # pad
        coords = np.pad(coords, ((0, 0), (0, n_axes - coords.shape[1])))
        explained = np.pad(explained, (0, n_axes - explained.shape[0]))
    return coords, explained
