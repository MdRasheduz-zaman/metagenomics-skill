"""Cross-sample community statistics — normalization, alpha/beta diversity, ordination.

Pure-Python (numpy only) so it runs anywhere and is unit-testable. Operates on a
sample x taxon count matrix built from the combined Bracken table.
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


def alpha_table(samples: List[str], mat: np.ndarray) -> List[Dict]:
    out = []
    for i, s in enumerate(samples):
        v = mat[i]
        out.append({"sample": s, "richness": richness(v),
                    "shannon": round(shannon(v), 4), "simpson": round(simpson(v), 4),
                    "pielou_evenness": round(pielou(v), 4)})
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
