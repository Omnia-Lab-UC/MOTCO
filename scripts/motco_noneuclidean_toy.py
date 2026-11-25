"""
Toy probe: what happens to MOTCO's magnitude / angle / shape when the latent space is non-linear?

Two omics blocks, 2 groups x 4 stages, 40 samples per cell. Both groups move along straight lines
in FEATURE space, so the ground truth is known exactly. The same stage-mean statistics (MOTCO's own
estimators) are then computed in several "spaces" built from the SNF fused network.

Run from the repository root:
    uv run python scripts/motco_noneuclidean_toy.py > results/snf-distance-statistics-toy-2026-11-25/summary.txt

Report: docs/reports/snf-distance-statistics-toy-2026-11-25.md

This is a sanity toy (1-D linear gradients, n = 320, 20 replicates, no RRPP), not a harness study.
It calls the private estimators in motco.stats.trajectory directly on cell means.
"""
import warnings
import numpy as np
from scipy.sparse.csgraph import shortest_path
from motco.stats.snf import get_affinity_matrix, SNF, get_spectral
from motco.stats import trajectory as T

warnings.filterwarnings("ignore")
K, NCELL, P = 4, 40, 50
CONTRAST = [[0, 1, 2, 3], [4, 5, 6, 7]]          # rows of the stage-mean matrix: A t0..t3, B t0..t3


# ----------------------------------------------------------------------------- simulation
def simulate(rng, sB=(0, 1, 2, 3), theta_deg=0.0, step=3.0, sd=0.2, jitter=0.35):
    """A sits at s = 0,1,2,3 along direction w. B sits at sB along w rotated by theta (feature space)."""
    s_stage = np.concatenate([np.repeat([0, 1, 2, 3.0], NCELL), np.repeat(np.asarray(sB, float), NCELL)])
    s = s_stage + rng.normal(0, jitter, s_stage.size)       # individual spread along the gradient
    grp = np.repeat([0, 1], K * NCELL)
    cell = np.repeat(np.arange(2 * K), NCELL)
    th, blocks = np.radians(theta_deg), []
    for _ in range(2):
        Q, _ = np.linalg.qr(rng.normal(size=(P, 2)))
        w, wp = Q[:, 0], Q[:, 1]
        dirs = np.where(grp[:, None] == 0, w[None, :], (np.cos(th) * w + np.sin(th) * wp)[None, :])
        blocks.append(step * s[:, None] * dirs + rng.normal(0, sd, (s.size, P)))
    return blocks, cell


# ----------------------------------------------------------------------------- estimators
def cell_means(Y, cell):
    return np.vstack([Y[cell == c].mean(0) for c in range(2 * K)])


def stats_from_points(M):
    """MOTCO's estimators on the 8 stage means -> (length ratio B/A, angle in degrees, shape distance)."""
    la, lb = T._estimate_size(M, CONTRAST[0]), T._estimate_size(M, CONTRAST[1])
    ua, ub = T._estimate_net_displacement(M, CONTRAST[0]), T._estimate_net_displacement(M, CONTRAST[1])
    return lb / la, np.degrees(np.arccos(np.clip(ua @ ub, -1, 1))), T._estimate_shape(M, CONTRAST)[0, 1]


def centroid_sqdist(D2, cell):
    """Gower: squared centroid-to-centroid distances from pairwise sample distances only."""
    C = 2 * K
    within = [D2[np.ix_(cell == a, cell == a)].mean() for a in range(C)]
    out = np.array([[D2[np.ix_(cell == a, cell == b)].mean() - 0.5 * within[a] - 0.5 * within[b]
                     for b in range(C)] for a in range(C)])
    np.fill_diagonal(out, 0)
    return np.clip(out, 0, None)


def stats_from_distances(D, cell):
    """Distance-only MOTCO: centroid distances -> PCoA of the 8 centroids -> the same estimators.
    Also returns the Berg-Nikolaev quadrilateral-cosine angle, which needs no coordinates at all."""
    cD2 = centroid_sqdist(D ** 2, cell)
    n = cD2.shape[0]
    J = np.eye(n) - 1 / n
    w, V = np.linalg.eigh(-0.5 * J @ cD2 @ J)
    M = V[:, w > 1e-10] * np.sqrt(w[w > 1e-10])
    a1, aK, b1, bK = 0, 3, 4, 7
    cosq = (cD2[a1, bK] + cD2[aK, b1] - cD2[a1, b1] - cD2[aK, bK]) / (2 * np.sqrt(cD2[a1, aK] * cD2[b1, bK]))
    return (*stats_from_points(M), np.degrees(np.arccos(np.clip(cosq, -1, 1))))


# ----------------------------------------------------------------------------- spaces
def _walk(W):
    W = (W + W.T) / 2
    np.fill_diagonal(W, 0)                     # SNF puts 0.5 on the diagonal; drop it before diffusing
    return W


def diffusion_map(W, t):
    """Eigenvalue-scaled spectral coordinates: Euclidean distance here = diffusion distance at time t."""
    W = _walk(W.copy()); d = W.sum(1)
    lam, v = np.linalg.eigh(W / np.sqrt(np.outer(d, d)))
    order = np.argsort(-lam)[1:]               # drop the trivial component
    lam, psi = lam[order], (v / np.sqrt(d)[:, None])[:, order]
    return (psi * lam ** t)[:, lam > 0]


def graph_geodesic(W, edge_len, k=10):
    """Shortest paths on the fused kNN graph (Isomap-style), edges weighted by feature-space length."""
    W = _walk(W.copy()); n = W.shape[0]; L = np.zeros((n, n))
    for i in range(n):
        nb = np.argsort(-W[i])[:k]
        L[i, nb] = edge_len[i, nb]
    return shortest_path(np.maximum(L, L.T), method="D", directed=False)


def graph_smooth(W, X, t=1):
    """Hybrid: use the fused graph only as a denoiser, keep measuring in feature space."""
    W = _walk(W.copy())
    return np.linalg.matrix_power(W / W.sum(1, keepdims=True), t) @ X


def one_replicate(seed, **sim):
    blocks, cell = simulate(np.random.default_rng(seed), **sim)
    X = np.hstack(blocks)
    Wf = SNF(get_affinity_matrix(blocks, K=20, eps=0.5), k=20, t=20)
    De = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    geo = stats_from_distances(graph_geodesic(Wf, De), cell)
    return {
        "feature space (ground truth)":             stats_from_points(cell_means(X, cell)),
        "SNF + get_spectral(10)  [MOTCO today]":    stats_from_points(cell_means(get_spectral(Wf.copy(), 10), cell)),
        "SNF + diffusion map, t=8":                 stats_from_points(cell_means(diffusion_map(Wf, 8), cell)),
        "SNF + diffusion map, t=32":                stats_from_points(cell_means(diffusion_map(Wf, 32), cell)),
        "SNF graph geodesic, distance-only":        geo[:3],
        "   same, angle by quadrilateral cosine":   (np.nan, geo[3], np.nan),
        "feature space after SNF-graph smoothing":  stats_from_points(cell_means(graph_smooth(Wf, X), cell)),
        "sanity: Euclidean D, distance-only":       stats_from_distances(De, cell)[:3],
    }


if __name__ == "__main__":
    R = 20
    scenarios = [
        ("NULL: B identical to A                      (truth: 1.0 /  0 / 0)", dict()),
        ("PURE MAGNITUDE: B goes 1.5x further         (truth: 1.5 /  0 / 0)", dict(sB=(0, 1.5, 3, 4.5))),
        ("PURE ORIENTATION: B rotated 30 deg          (truth: 1.0 / 30 / 0)", dict(theta_deg=30)),
        ("PURE ORIENTATION: B rotated 60 deg          (truth: 1.0 / 60 / 0)", dict(theta_deg=60)),
    ]
    for label, sim in scenarios:
        res = [one_replicate(100 + r, **sim) for r in range(R)]
        print(f"\n### {label}   mean of {R} replicates")
        print(f"{'space':42s} {'len B/A':>8s} {'angle':>8s} {'shape':>8s}")
        for key in res[0]:
            m = np.nanmean(np.array([r[key] for r in res], float), 0)
            print(f"{key:42s} " + " ".join(f"{v:8.2f}" if np.isfinite(v) else f"{'':8s}" for v in m))
