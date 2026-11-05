"""The stage-span orientation's alignment axis (``span_alignment``, ``span_support``)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from motco.simulations.diagnostics import geometry_from_means
from motco.simulations.evaluation import SimulationEvaluationParams
from motco.simulations.generator import logit
from motco.simulations.grid import _to_jsonable, make_simulation_cell, parameter_signature
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import (
    OUT_OF_SPAN_TOLERANCE,
    SemiSyntheticTrajectoryError,
    SemiSyntheticTrajectoryParams,
    generate_semisynthetic_trajectory,
    null_pooled_scale,
    out_of_span_direction,
)

TWO_BLOCKS = ("methylation", "expression")


@pytest.fixture(scope="module")
def reference():
    return load_reference()


def _span(**overrides) -> SemiSyntheticTrajectoryParams:
    base = {
        "seed": 23,
        "n_samples": 120,
        "p_dmp": 0.1,
        "construction": "stage_span",
        "span_layers": TWO_BLOCKS,
        "trajectory_mode": "orientation",
        "group_effect_size": 0.5,
    }
    return SemiSyntheticTrajectoryParams(**{**base, **overrides})


def _null_configuration(dataset, reference) -> tuple[np.ndarray, np.ndarray]:
    """Group A's (null) standardized stage configuration over the span, and the scale."""

    layers = tuple(dataset.truth["span_layers"])
    population = dataset.population_trajectories.layers
    null_means = {layer: population[layer].loc["A"].to_numpy() for layer in layers}
    weights = dataset.metadata["stage"].value_counts().sort_index().to_numpy(dtype=float)
    scale = null_pooled_scale(null_means, weights, layers=layers, reference=reference)
    configuration = np.hstack([null_means[layer] for layer in layers])
    return (configuration - weights @ configuration / weights.sum()) / scale, scale


def _standardized_population(dataset, reference) -> pd.DataFrame:
    layers = tuple(dataset.truth["span_layers"])
    _, scale = _null_configuration(dataset, reference)
    population = dataset.population_trajectories.layers
    return pd.concat([population[layer] for layer in layers], axis=1) / scale


# --------------------------------------------------------------------------- #
# 1.1 Parameters, validation and identity
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"span_alignment": -1.0}, r"span_alignment must be an angle in degrees in \[0, 90\]"),
        ({"span_alignment": 90.5}, r"span_alignment must be an angle in degrees in \[0, 90\]"),
        ({"span_alignment": float("nan")}, r"span_alignment must be an angle in degrees in \[0, 90\]"),
        ({"span_alignment": "45"}, r"span_alignment must be an angle in degrees in \[0, 90\]"),
        ({"span_alignment": True}, r"span_alignment must be an angle in degrees in \[0, 90\]"),
        ({"span_alignment": 30.0, "construction": None, "span_layers": None}, "construction is not 'stage_span'"),
        ({"span_alignment": 30.0, "trajectory_mode": "magnitude"}, "trajectory_mode='magnitude'"),
        ({"span_alignment": 30.0, "trajectory_mode": "shape"}, "trajectory_mode='shape'"),
        ({"span_alignment": 30.0, "trajectory_mode": "translation"}, "trajectory_mode='translation'"),
        ({"span_support": 20}, "span_support=20 is set but span_alignment is absent or 0"),
        ({"span_support": 20, "span_alignment": 0}, "span_support=20 is set but span_alignment is absent or 0"),
        ({"span_support": 0, "span_alignment": 30.0}, "span_support must be a positive integer"),
        ({"span_support": -3, "span_alignment": 30.0}, "span_support must be a positive integer"),
        ({"span_support": 2.5, "span_alignment": 30.0}, "span_support must be a positive integer"),
        ({"span_support": 10_000, "span_alignment": 30.0}, r"span_support=10000 exceeds the \d+ span-layer"),
    ],
)
def test_invalid_alignment_is_rejected(reference, overrides, match) -> None:
    with pytest.raises(SemiSyntheticTrajectoryError, match=match):
        generate_semisynthetic_trajectory(_span(**overrides), reference=reference)


def test_zero_alignment_normalizes_to_absent() -> None:
    assert _span(span_alignment=0).span_alignment is None
    assert _span(span_alignment=0.0) == _span()


def _cell(generator: SemiSyntheticTrajectoryParams):
    return make_simulation_cell(
        phase="power_primary",
        generator_params=generator,
        evaluation_params=SimulationEvaluationParams(integration_method="pls", permutations=9),
    )


def test_absent_and_zero_alignment_keep_the_cell_identity(reference) -> None:
    plain, zero = _span(), _span(span_alignment=0.0)
    assert _cell(zero).cell_id == _cell(plain).cell_id
    assert parameter_signature(_cell(zero)) == parameter_signature(_cell(plain))
    payload = _to_jsonable(plain)
    assert "span_alignment" not in payload and "span_support" not in payload
    tilted = _span(span_alignment=30.0)
    assert _cell(tilted).cell_id != _cell(plain).cell_id
    assert _cell(replace(tilted, span_support=20)).cell_id != _cell(tilted).cell_id

    expected = generate_semisynthetic_trajectory(plain, reference=reference)
    dataset = generate_semisynthetic_trajectory(zero, reference=reference)
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(dataset, layer).to_numpy().tobytes() == getattr(expected, layer).to_numpy().tobytes()
    assert dataset.truth.keys() == expected.truth.keys()
    assert "span_alignment" not in dataset.truth
    assert "out_of_span" not in dataset.truth["stage_span"]
    assert dataset.truth["transform"] == expected.truth["transform"]


# --------------------------------------------------------------------------- #
# 1.2 The out-of-span direction
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("support", [None, 20])
def test_out_of_span_direction_lives_on_unused_features(reference, support) -> None:
    dataset = generate_semisynthetic_trajectory(_span(trajectory_mode="none", group_effect_size=0.0),
                                                reference=reference)
    z, _ = _null_configuration(dataset, reference)
    direction = out_of_span_direction(z, seed=23, support=support)
    o = direction.vector
    centered = z - z.mean(axis=0)
    unused = np.max(np.abs(centered), axis=0) <= OUT_OF_SPAN_TOLERANCE

    assert np.linalg.norm(o) == pytest.approx(1.0, abs=1e-12)
    np.testing.assert_allclose(centered @ o, 0.0, atol=1e-10)
    assert not np.any(o[~unused])
    assert direction.available_features == int(unused.sum())
    if support is None:
        assert np.count_nonzero(o) == direction.support_size == int(unused.sum())
    else:
        assert np.count_nonzero(o) == direction.support_size == support
    # It is a separate stream: the same call repeats exactly.
    again = out_of_span_direction(z, seed=23, support=support)
    assert again.vector.tobytes() == o.tobytes() and again.direction_id == direction.direction_id
    assert out_of_span_direction(z, seed=24, support=support).direction_id != direction.direction_id


def test_too_few_unused_features_is_named() -> None:
    z = np.array([[0.0, 1.0, 5.0, 5.0], [1.0, 0.0, 5.0, 5.0], [2.0, 2.0, 5.0, 5.0]])
    with pytest.raises(SemiSyntheticTrajectoryError, match="span_support=3 exceeds the 2 span-layer features"):
        out_of_span_direction(z, seed=1, support=3)
    with pytest.raises(SemiSyntheticTrajectoryError, match="none of the 2 features qualifies"):
        out_of_span_direction(z[:, :2], seed=1)


@pytest.mark.parametrize("support", [None, 20])
def test_the_direction_is_shared_across_effects_and_alignments(reference, support) -> None:
    identifiers = set()
    for effect, alignment in ((0.25, 30.0), (0.5, 30.0), (0.5, 90.0), (0.0, 60.0)):
        dataset = generate_semisynthetic_trajectory(
            _span(group_effect_size=effect, span_alignment=alignment, span_support=support), reference=reference
        )
        block = dataset.truth["stage_span"]["out_of_span"]
        identifiers.add((block["direction_id"], block["support_size"]))
    assert len(identifiers) == 1


# --------------------------------------------------------------------------- #
# 1.3 The tilted rotation: fixed geometry, out-of-span difference, null + plant
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("n_stages", [3, 4])
@pytest.mark.parametrize("alignment", [30.0, 60.0, 90.0])
@pytest.mark.parametrize("effect", [0.25, 0.5])
def test_intended_geometry_does_not_depend_on_alignment(reference, n_stages, alignment, effect) -> None:
    dataset = generate_semisynthetic_trajectory(
        _span(n_stages=n_stages, group_effect_size=effect, span_alignment=alignment), reference=reference
    )
    intended = dataset.truth["stage_span"]["intended_geometry"]
    assert intended["angle"] == pytest.approx(90.0 * effect)
    assert intended["delta"] == pytest.approx(0.0, abs=1e-9)
    assert intended["shape"] == pytest.approx(0.0, abs=1e-10)

    stages = [str(stage) for stage in range(n_stages)]
    measured = geometry_from_means(
        _standardized_population(dataset, reference), ["A", "B"], stages, orientation_estimator="net_displacement"
    )
    assert measured.angle == pytest.approx(90.0 * effect)
    assert measured.delta == pytest.approx(0.0, abs=1e-9)
    assert measured.shape == pytest.approx(0.0, abs=1e-10)
    assert dataset.truth["transform"] == {
        "construction": "stage_span",
        "rotation_degrees": 90.0 * effect,
        "span_alignment": alignment,
    }


@pytest.mark.parametrize("n_stages", [3, 4])
@pytest.mark.parametrize("support", [None, 20])
def test_full_alignment_leaves_the_pooled_stage_span(reference, n_stages, support) -> None:
    dataset = generate_semisynthetic_trajectory(
        _span(n_stages=n_stages, n_samples=40 * n_stages, span_alignment=90.0, span_support=support),
        reference=reference,
    )
    sizes = dataset.metadata.groupby(["group", "stage"]).size()
    assert sizes.nunique() == 1  # balanced cells
    joint = _standardized_population(dataset, reference)
    group_a, group_b = joint.loc["A"].to_numpy(), joint.loc["B"].to_numpy()
    difference = group_b - group_a
    varying = difference - difference.mean(axis=0)
    pooled = (group_a + group_b) / 2.0
    pooled = pooled - pooled.mean(axis=0)
    assert np.linalg.norm(varying) > 0.1
    np.testing.assert_allclose(varying @ pooled.T, 0.0, atol=1e-10)

    # At φ = 0 the same rotation lies inside the span instead.
    in_span = generate_semisynthetic_trajectory(
        _span(n_stages=n_stages, n_samples=40 * n_stages), reference=reference
    )
    joint = _standardized_population(in_span, reference)
    difference = joint.loc["B"].to_numpy() - joint.loc["A"].to_numpy()
    varying = difference - difference.mean(axis=0)
    pooled = (joint.loc["A"].to_numpy() + joint.loc["B"].to_numpy()) / 2.0
    basis, _ = np.linalg.qr((pooled - pooled.mean(axis=0)).T)
    assert np.sum((varying @ basis) ** 2) / np.sum(varying**2) == pytest.approx(1.0, abs=1e-8)


@pytest.mark.parametrize(("alignment", "support"), [(30.0, None), (90.0, None), (90.0, 20)])
def test_aligned_dataset_is_the_null_dataset_plus_the_plant(reference, alignment, support) -> None:
    params = _span(span_alignment=alignment, span_support=support)
    null = generate_semisynthetic_trajectory(
        replace(params, trajectory_mode="none", group_effect_size=0.0, construction=None, span_layers=None,
                span_alignment=None, span_support=None),
        reference=reference,
    )
    dataset = generate_semisynthetic_trajectory(params, reference=reference)

    group_a = dataset.metadata["group"] == "A"
    group_b = ~group_a
    stages_b = dataset.metadata.loc[group_b, "stage"].to_numpy()
    assert dataset.metadata.equals(null.metadata)
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(dataset, layer).loc[group_a].equals(getattr(null, layer).loc[group_a])
    assert dataset.proteomics.equals(null.proteomics)

    population = dataset.population_trajectories.layers
    null_population = null.population_trajectories.layers
    for layer in TWO_BLOCKS:
        plant = population[layer].loc["B"].to_numpy() - null_population[layer].loc["B"].to_numpy()
        assert np.abs(plant).max() > 0
        if layer == "methylation":
            difference = logit(dataset.methylation.loc[group_b].to_numpy()) - logit(
                null.methylation.loc[group_b].to_numpy()
            )
            atol = 1e-8
        else:
            difference = dataset.expression.loc[group_b].to_numpy() - null.expression.loc[group_b].to_numpy()
            atol = 1e-10
        np.testing.assert_allclose(difference, plant[stages_b], atol=atol)


def test_truth_names_the_alignment(reference) -> None:
    dataset = generate_semisynthetic_trajectory(_span(span_alignment=60, span_support=20), reference=reference)
    truth = dataset.truth
    assert truth["span_alignment"] == 60.0
    assert truth["span_support"] == 20
    block = truth["stage_span"]["out_of_span"]
    assert block["support_size"] == 20
    assert block["span_layers"] == list(TWO_BLOCKS)
    assert block["available_features"] >= 20
    assert isinstance(block["direction_id"], str) and len(block["direction_id"]) == 12
    dense = generate_semisynthetic_trajectory(_span(span_alignment=60), reference=reference).truth
    assert dense["span_support"] is None
    assert dense["stage_span"]["out_of_span"]["support_size"] == block["available_features"]
