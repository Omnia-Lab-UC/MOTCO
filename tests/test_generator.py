"""Tests for the numpy generative core and its realism vs InterSIM.

The realism gate compares the numpy generator against a committed InterSIM
fixture (``tests/data/intersim_realism_fixture.npz``) generated once in R with
``delta=0`` -- which reduces InterSIM to pure baseline MVN + ``rev.logit``
sampling, giving an exact model-equivalence check with no R dependency in CI.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from motco.simulations.generator import (
    GeneratorError,
    bernoulli_indicators,
    derive_coupled_indicators,
    generate_omics,
    rev_logit,
)
from motco.simulations.reference import (
    ReferenceCacheMissingError,
    load_reference,
)

FIXTURE = Path(__file__).parent / "data" / "intersim_realism_fixture.npz"

# Documented tolerances for the delta=0 numpy-vs-InterSIM comparison at
# n_sample=4000. Per-feature means agree tightly; per-feature variances carry
# finite-sample Monte-Carlo noise (relative SE ~ sqrt(2/N)), so a single
# feature may drift up to ~30% while the average stays well under 10%.
MEAN_ABS_MAX = 0.15
MEAN_ABS_AVG = 0.05
VAR_REL_MAX = 0.30
VAR_REL_AVG = 0.10


@pytest.fixture(scope="module")
def reference():
    return load_reference()


def test_rev_logit_maps_into_unit_interval():
    x = np.array([-10.0, 0.0, 10.0])
    out = rev_logit(x)
    assert np.all((out > 0) & (out < 1))
    assert out[1] == pytest.approx(0.5)


def test_load_reference_missing_cache_raises(tmp_path):
    with pytest.raises(ReferenceCacheMissingError, match="reference cache not found"):
        load_reference(tmp_path / "absent.npz")


def test_reference_shapes_and_provenance(reference):
    assert reference.n_cpg == 367
    assert reference.n_gene == 131
    assert reference.n_protein == 160
    assert reference.cov_M.shape == (367, 367)
    assert reference.incidence_cpg_gene.shape == (367, 131)
    assert reference.incidence_gene_protein.shape == (131, 160)
    assert reference.provenance["intersim_version"]


def test_generate_omics_shapes_and_alignment(reference):
    rng = np.random.default_rng(0)
    ind_m = bernoulli_indicators(rng, reference.n_cpg, 3, 0.2)
    ind_e, ind_p = derive_coupled_indicators(ind_m, reference)
    g = generate_omics(
        cell_sizes=[100, 100, 150],
        indicators_methyl=ind_m,
        indicators_expr=ind_e,
        indicators_protein=ind_p,
        rng=rng,
        reference=reference,
    )
    assert g.methylation.shape == (350, 367)
    assert g.expression.shape == (350, 131)
    assert g.proteomics.shape == (350, 160)
    assert g.cell_ids.shape == (350,)
    assert np.bincount(g.cell_ids).tolist() == [100, 100, 150]


def test_methylation_stays_in_unit_interval(reference):
    rng = np.random.default_rng(1)
    ind_m = bernoulli_indicators(rng, reference.n_cpg, 2, 0.5)
    ind_e, ind_p = derive_coupled_indicators(ind_m, reference)
    g = generate_omics(
        cell_sizes=[50, 50],
        indicators_methyl=ind_m,
        indicators_expr=ind_e,
        indicators_protein=ind_p,
        delta_methyl=5.0,
        rng=rng,
        reference=reference,
    )
    assert np.all((g.methylation > 0) & (g.methylation < 1))


def test_coupled_indicators_follow_incidence(reference):
    rng = np.random.default_rng(2)
    ind_m = bernoulli_indicators(rng, reference.n_cpg, 2, 0.3)
    ind_e, ind_p = derive_coupled_indicators(ind_m, reference)
    # A gene is differential iff some mapped CpG is differential.
    expected_e = (reference.incidence_cpg_gene.T @ ind_m > 0).astype(float)
    assert np.array_equal(ind_e, expected_e)
    expected_p = (reference.incidence_gene_protein.T @ ind_e > 0).astype(float)
    assert np.array_equal(ind_p, expected_p)


def test_generation_is_reproducible(reference):
    def run():
        rng = np.random.default_rng(7)
        ind_m = bernoulli_indicators(rng, reference.n_cpg, 3, 0.2)
        ind_e, ind_p = derive_coupled_indicators(ind_m, reference)
        return generate_omics(
            cell_sizes=[40, 40, 40],
            indicators_methyl=ind_m,
            indicators_expr=ind_e,
            indicators_protein=ind_p,
            rng=rng,
            reference=reference,
        )

    a, b = run(), run()
    assert np.array_equal(a.methylation, b.methylation)
    assert np.array_equal(a.expression, b.expression)
    assert np.array_equal(a.proteomics, b.proteomics)


def test_bad_indicator_shape_raises(reference):
    rng = np.random.default_rng(0)
    with pytest.raises(GeneratorError, match="methylation indicators must have shape"):
        generate_omics(
            cell_sizes=[10, 10],
            indicators_methyl=np.zeros((reference.n_cpg, 3)),
            indicators_expr=np.zeros((reference.n_gene, 2)),
            indicators_protein=np.zeros((reference.n_protein, 2)),
            rng=rng,
            reference=reference,
        )


def _shifted_run(reference, **shifts):
    rng = np.random.default_rng(11)
    ind_m = bernoulli_indicators(rng, reference.n_cpg, 3, 0.2)
    ind_e, ind_p = derive_coupled_indicators(ind_m, reference)
    return generate_omics(
        cell_sizes=[30, 0, 25],
        indicators_methyl=ind_m,
        indicators_expr=ind_e,
        indicators_protein=ind_p,
        rng=rng,
        reference=reference,
        **shifts,
    )


def test_zero_cell_shift_is_byte_identical(reference):
    plain = _shifted_run(reference)
    zero = _shifted_run(
        reference,
        cell_shift_methyl=np.zeros((reference.n_cpg, 3)),
        cell_shift_expr=np.zeros((reference.n_gene, 3)),
        cell_shift_protein=np.zeros((reference.n_protein, 3)),
    )
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(plain, layer).tobytes() == getattr(zero, layer).tobytes()


def test_cell_shift_moves_samples_by_exactly_the_shift(reference):
    shift_rng = np.random.default_rng(3)
    shift_m = shift_rng.normal(size=(reference.n_cpg, 3))
    shift_e = shift_rng.normal(size=(reference.n_gene, 3))
    plain = _shifted_run(reference)
    shifted = _shifted_run(reference, cell_shift_methyl=shift_m, cell_shift_expr=shift_e)
    cells = plain.cell_ids
    assert np.array_equal(cells, shifted.cell_ids)
    # Methylation is shifted in M-values, before the inverse logit.
    to_m = lambda beta: np.log(beta / (1.0 - beta))  # noqa: E731
    np.testing.assert_allclose(to_m(shifted.methylation) - to_m(plain.methylation), shift_m.T[cells], atol=1e-8)
    np.testing.assert_allclose(shifted.expression - plain.expression, shift_e.T[cells], atol=1e-10)
    assert np.array_equal(shifted.proteomics, plain.proteomics)


def test_bad_cell_shift_shape_raises(reference):
    with pytest.raises(GeneratorError, match="cell shift must have shape"):
        _shifted_run(reference, cell_shift_expr=np.zeros((reference.n_gene, 2)))


@pytest.mark.slow
def test_realism_matches_intersim_fixture(reference):
    """numpy delta=0 output matches the committed InterSIM fixture within tolerance."""

    fixture = dict(np.load(FIXTURE))
    sizes = [1200, 1200, 1600]  # cluster.sample.prop=(0.3,0.3,0.4) of 4000
    rng = np.random.default_rng(int(fixture["seed"]))
    g = generate_omics(
        cell_sizes=sizes,
        indicators_methyl=np.zeros((reference.n_cpg, 3)),
        indicators_expr=np.zeros((reference.n_gene, 3)),
        indicators_protein=np.zeros((reference.n_protein, 3)),
        delta_methyl=0.0,
        delta_expr=0.0,
        delta_protein=0.0,
        rng=rng,
        reference=reference,
    )
    for name, mat in (
        ("methyl", g.methylation),
        ("expr", g.expression),
        ("protein", g.proteomics),
    ):
        mean_abs = np.abs(mat.mean(0) - fixture[f"{name}_mean"])
        var_rel = np.abs(mat.var(0) - fixture[f"{name}_var"]) / np.maximum(
            np.abs(fixture[f"{name}_var"]), 1e-9
        )
        assert mean_abs.max() < MEAN_ABS_MAX, f"{name} mean drift {mean_abs.max()}"
        assert mean_abs.mean() < MEAN_ABS_AVG, f"{name} mean avg drift {mean_abs.mean()}"
        assert var_rel.max() < VAR_REL_MAX, f"{name} var drift {var_rel.max()}"
        assert var_rel.mean() < VAR_REL_AVG, f"{name} var avg drift {var_rel.mean()}"
