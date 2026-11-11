# tests/test_snf.py
from __future__ import annotations

import numpy as np
import pytest

from motco.stats.snf import SNF, get_affinity_matrix, get_spectral


def _datasets(n: int = 20, p: int = 5, n_datasets: int = 2, seed: int = 0) -> list[np.ndarray]:
    rng = np.random.default_rng(seed)
    return [rng.standard_normal((n, p)) for _ in range(n_datasets)]


def test_affinity_matrix_returns_correct_count():
    dats = _datasets()
    Ws = get_affinity_matrix(dats, K=5)
    assert len(Ws) == len(dats)


def test_affinity_matrix_shapes_are_square():
    dats = _datasets(n=20)
    Ws = get_affinity_matrix(dats, K=5)
    for W in Ws:
        assert W.shape == (20, 20)


def test_affinity_matrix_non_negative():
    dats = _datasets()
    Ws = get_affinity_matrix(dats, K=5)
    for W in Ws:
        assert np.all(W >= 0)


def test_affinity_matrix_symmetric():
    dats = _datasets()
    Ws = get_affinity_matrix(dats, K=5)
    for W in Ws:
        assert np.allclose(W, W.T)


def test_snf_output_shape():
    n = 20
    dats = _datasets(n=n)
    Ws = get_affinity_matrix(dats, K=5)
    fused = SNF(Ws, k=5, t=5)
    assert fused.shape == (n, n)


def test_snf_output_approximate_symmetric():
    dats = _datasets()
    Ws = get_affinity_matrix(dats, K=5)
    fused = SNF(Ws, k=5, t=5)
    # SNF output is approximately symmetric due to averaging steps
    assert np.allclose(fused, fused.T, atol=0.05)


def test_get_spectral_default_shape():
    dats = _datasets(n=20)
    Ws = get_affinity_matrix(dats, K=5)
    fused = SNF(Ws, k=5, t=5)
    emb = get_spectral(fused)
    assert emb.shape == (20, 10)  # default 10 components


def test_get_spectral_custom_components():
    dats = _datasets(n=20)
    Ws = get_affinity_matrix(dats, K=5)
    fused = SNF(Ws, k=5, t=5)
    emb = get_spectral(fused, n_components=5)
    assert emb.shape == (20, 5)


def _square_affinity(n=10):
    rng = np.random.default_rng(0)
    W = rng.uniform(0, 1, (n, n))
    return (W + W.T) / 2


# --- get_affinity_matrix: misaligned datasets ---
def test_affinity_matrix_row_mismatch():
    rng = np.random.default_rng(0)
    d1 = rng.standard_normal((10, 5))
    d2 = rng.standard_normal((8, 5))
    with pytest.raises(ValueError, match="8 rows"):
        get_affinity_matrix([d1, d2], K=3)


# --- get_affinity_matrix: K >= n_samples ---
def test_affinity_matrix_K_too_large():
    rng = np.random.default_rng(0)
    d1 = rng.standard_normal((5, 3))
    d2 = rng.standard_normal((5, 3))
    with pytest.raises(ValueError, match="K=10"):
        get_affinity_matrix([d1, d2], K=10)


# --- get_affinity_matrix: eps <= 0 ---
def test_affinity_matrix_invalid_eps():
    rng = np.random.default_rng(0)
    d1 = rng.standard_normal((10, 3))
    d2 = rng.standard_normal((10, 3))
    with pytest.raises(ValueError, match="eps"):
        get_affinity_matrix([d1, d2], K=3, eps=-0.1)


# --- get_affinity_matrix: NaN in data ---
def test_affinity_matrix_nan():
    rng = np.random.default_rng(0)
    d1 = rng.standard_normal((10, 3))
    d2 = rng.standard_normal((10, 3))
    d2[0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN"):
        get_affinity_matrix([d1, d2], K=3)


# --- SNF: non-square matrix ---
def test_snf_non_square():
    W1 = _square_affinity(5)
    W2 = np.ones((5, 6))
    with pytest.raises(ValueError, match="not square"):
        SNF([W1, W2])


# --- SNF: shape mismatch ---
def test_snf_shape_mismatch():
    W1 = _square_affinity(5)
    W2 = _square_affinity(6)
    with pytest.raises(ValueError, match="shape"):
        SNF([W1, W2])


# --- SNF: k >= n_samples ---
def test_snf_k_too_large():
    W1 = _square_affinity(5)
    W2 = _square_affinity(5)
    with pytest.raises(ValueError, match="k=10"):
        SNF([W1, W2], k=10)
