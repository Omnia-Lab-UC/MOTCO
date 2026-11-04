"""The opt-in stage-span construction of the semi-synthetic generator."""

from __future__ import annotations

import hashlib
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from motco.simulations.diagnostics import geometry_from_means
from motco.simulations.evaluation import SimulationEvaluationParams
from motco.simulations.generator import logit, omic_population_means
from motco.simulations.grid import make_simulation_cell, parameter_signature
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import (
    STAGE_SPAN_RANK_TOLERANCE,
    SemiSyntheticTrajectoryError,
    SemiSyntheticTrajectoryParams,
    edit_stage_configuration,
    expected_surgery_headroom,
    generate_semisynthetic_trajectory,
    null_pooled_scale,
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
    }
    return SemiSyntheticTrajectoryParams(**{**base, **overrides})


# --------------------------------------------------------------------------- #
# 2.1 Parameters, validation and identity
# --------------------------------------------------------------------------- #


def _dataset_digest(dataset) -> str:
    """Digest of the BLAS-free parts: indicators (the RNG stream), metadata and population means.

    The sampled matrices are left out because ``multivariate_normal``'s SVD makes
    their last bits depend on the BLAS thread count; they are compared in-process
    below instead.
    """

    digest = hashlib.sha256()
    for group in ("A", "B"):
        for layer in ("methylation", "expression", "proteomics"):
            digest.update(np.ascontiguousarray(dataset.truth["indicators"][group][layer], dtype=float).tobytes())
    digest.update(dataset.metadata.to_csv().encode())
    for frame in dataset.population_trajectories.layers.values():
        digest.update(frame.to_numpy().tobytes())
    return digest.hexdigest()[:16]


# Captured at 85b74bb, before the construction existed (identical at one BLAS thread
# and unpinned): BLAS-free dataset digest, cell id and the signature prefix of a PLS
# cell at n = 90, p_dmp 0.1, seed 17. The full sampled matrices were also checked
# byte-identical against 85b74bb at a fixed thread count.
_PRE_CHANGE = {
    ("none", 0.0): ("49ba4b7a4eec8993", "power_primary-8e219c6c2efb", "4fffdbf1c49a5343"),
    ("magnitude", 0.1): ("85efc1e99532b355", "power_primary-2742d5055026", "39c80e3d6d1519aa"),
    ("orientation", 0.5): ("adc74d6b517a86f9", "power_primary-f222e49a3500", "345a4d6aa3d17684"),
    ("shape", 0.5): ("517921913ba0523b", "power_primary-985b701ea292", "1812c159801168ec"),
    ("translation", 0.5): ("611b4dd900dd33bb", "power_primary-4d4ff49f7694", "049822c065b15f79"),
}


@pytest.mark.parametrize(("mode", "effect"), sorted(_PRE_CHANGE))
def test_default_construction_is_byte_identical_to_before(reference, mode, effect) -> None:
    params = SemiSyntheticTrajectoryParams(
        seed=17, trajectory_mode=mode, group_effect_size=effect, n_samples=90, p_dmp=0.1
    )
    surgery = replace(params, construction="surgery")
    assert surgery.construction is None and surgery == params
    cell = make_simulation_cell(
        phase="power_primary",
        generator_params=params,
        evaluation_params=SimulationEvaluationParams(integration_method="pls", permutations=9),
    )
    digest, cell_id, signature = _PRE_CHANGE[(mode, effect)]
    dataset = generate_semisynthetic_trajectory(params, reference=reference)
    assert _dataset_digest(dataset) == digest
    assert "construction" not in dataset.truth and "stage_span" not in dataset.truth
    assert cell.cell_id == cell_id
    assert parameter_signature(cell)[:16] == signature
    spelled = generate_semisynthetic_trajectory(surgery, reference=reference)
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(spelled, layer).to_numpy().tobytes() == getattr(dataset, layer).to_numpy().tobytes()


def test_span_layers_normalize_to_canonical_order() -> None:
    params = _span(span_layers=["expression", "methylation"])
    assert params.span_layers == TWO_BLOCKS
    assert params == _span()


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"construction": "lattice"}, "Unknown construction: 'lattice'"),
        ({"span_layers": ("methylation", "lipids")}, r"unknown layer\(s\) \['lipids'\]"),
        ({"span_layers": ("expression", "expression")}, "repeats a layer"),
        ({"span_layers": ()}, "must not be empty"),
        ({"construction": None}, "span_layers .* is set but construction is not 'stage_span'"),
        ({"construction": "surgery"}, "span_layers .* is set but construction is not 'stage_span'"),
        ({"trajectory_mode": "orientation", "group_effect_size": 2.5}, "must be at most 2"),
        ({"trajectory_mode": "orientation", "group_effect_size": 0.5, "n_stages": 2}, "no second direction"),
        ({"trajectory_mode": "shape", "group_effect_size": 0.5, "n_stages": 2}, "no second direction"),
    ],
)
def test_invalid_stage_span_params_are_rejected(reference, overrides, match) -> None:
    with pytest.raises(SemiSyntheticTrajectoryError, match=match):
        generate_semisynthetic_trajectory(_span(**overrides), reference=reference)


def test_absent_span_layers_span_all_three_blocks(reference) -> None:
    dataset = generate_semisynthetic_trajectory(
        _span(span_layers=None, trajectory_mode="magnitude", group_effect_size=0.1), reference=reference
    )
    assert dataset.truth["span_layers"] == ["methylation", "expression", "proteomics"]


def test_stage_span_cells_have_no_surgery_headroom(reference) -> None:
    for mode in ("orientation", "shape", "translation"):
        params = _span(trajectory_mode=mode, group_effect_size=2.0)
        assert expected_surgery_headroom(params, reference=reference) is None


# --------------------------------------------------------------------------- #
# 2.2 Null population pooled scale
# --------------------------------------------------------------------------- #


def test_null_pooled_scale_matches_a_large_sample(reference) -> None:
    params = SemiSyntheticTrajectoryParams(seed=4, n_samples=9000, p_dmp=0.1, stage_sample_prop=(0.5, 0.3, 0.2))
    dataset = generate_semisynthetic_trajectory(params, reference=reference)
    indicators = dataset.truth["indicators"]["A"]
    null_means = omic_population_means(
        indicators_methyl=indicators["methylation"],
        indicators_expr=indicators["expression"],
        indicators_protein=indicators["proteomics"],
        reference=reference,
    )
    weights = dataset.metadata["stage"].value_counts().sort_index().to_numpy(dtype=float)
    layers = ("methylation", "expression", "proteomics")
    expected = null_pooled_scale(null_means, weights, layers=layers, reference=reference)
    observed = np.concatenate(
        [
            logit(dataset.methylation.to_numpy()).std(axis=0),
            dataset.expression.to_numpy().std(axis=0),
            dataset.proteomics.to_numpy().std(axis=0),
        ]
    )
    relative = np.abs(observed / expected - 1.0)
    # The SD of a sample SD is about SD / sqrt(2n): 0.75 % at n = 9,000.
    assert np.median(relative) < 0.01
    assert relative.max() < 0.04


# --------------------------------------------------------------------------- #
# 2.3 The four edits and their degeneracy errors
# --------------------------------------------------------------------------- #


def _standardized_population(dataset, reference) -> pd.DataFrame:
    """B's and A's population means divided by the null pooled scale over the span."""

    layers = tuple(dataset.truth["span_layers"])
    population = dataset.population_trajectories.layers
    null_means = {layer: population[layer].loc["A"].to_numpy() for layer in layers}
    weights = dataset.metadata["stage"].value_counts().sort_index().to_numpy(dtype=float)
    scale = null_pooled_scale(null_means, weights, layers=layers, reference=reference)
    joint = pd.concat([population[layer] for layer in layers], axis=1)
    return joint / scale


def _population_geometry(dataset, reference):
    stages = [str(stage) for stage in range(dataset.truth["n_stages"])]
    return geometry_from_means(
        _standardized_population(dataset, reference), ["A", "B"], stages, orientation_estimator="net_displacement"
    )


@pytest.mark.parametrize("n_stages", [3, 4])
@pytest.mark.parametrize("effect", [0.25, 1.0])
def test_stage_span_edits_realize_their_intended_geometry(reference, n_stages, effect) -> None:
    base = {"n_stages": n_stages, "group_effect_size": effect}
    magnitude = generate_semisynthetic_trajectory(_span(trajectory_mode="magnitude", **base), reference=reference)
    geometry = _population_geometry(magnitude, reference)
    assert geometry.path_lengths["B"] / geometry.path_lengths["A"] == pytest.approx(1.0 + effect)
    assert geometry.angle == pytest.approx(0.0, abs=1e-5)
    assert geometry.shape == pytest.approx(0.0, abs=1e-10)

    orientation = generate_semisynthetic_trajectory(_span(trajectory_mode="orientation", **base), reference=reference)
    geometry = _population_geometry(orientation, reference)
    assert geometry.angle == pytest.approx(90.0 * effect)
    assert geometry.delta == pytest.approx(0.0, abs=1e-9)
    assert geometry.shape == pytest.approx(0.0, abs=1e-10)

    shape = generate_semisynthetic_trajectory(_span(trajectory_mode="shape", **base), reference=reference)
    geometry = _population_geometry(shape, reference)
    assert geometry.shape > 1e-3
    assert geometry.delta == pytest.approx(0.0, abs=1e-9)
    assert geometry.angle == pytest.approx(0.0, abs=1e-5)
    assert shape.truth["transform"]["shape_stage"] == list(range(1, n_stages - 1))[(n_stages - 2) // 2]

    translation = generate_semisynthetic_trajectory(_span(trajectory_mode="translation", **base), reference=reference)
    geometry = _population_geometry(translation, reference)
    assert geometry.delta == pytest.approx(0.0, abs=1e-9)
    assert geometry.angle == pytest.approx(0.0, abs=1e-5)
    assert geometry.shape == pytest.approx(0.0, abs=1e-10)

    for dataset, measured in (
        (magnitude, _population_geometry(magnitude, reference)),
        (orientation, _population_geometry(orientation, reference)),
        (shape, _population_geometry(shape, reference)),
    ):
        intended = dataset.truth["stage_span"]["intended_geometry"]
        assert intended["delta"] == pytest.approx(measured.delta, abs=1e-8)
        assert intended["angle"] == pytest.approx(measured.angle, abs=1e-5)
        assert intended["shape"] == pytest.approx(measured.shape, abs=1e-8)


@pytest.mark.parametrize("mode", ["orientation", "shape"])
def test_rank_one_configuration_has_no_second_direction(mode) -> None:
    collinear = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [3.0, 6.0, 0.0]])
    with pytest.raises(SemiSyntheticTrajectoryError, match="no second direction"):
        edit_stage_configuration(collinear, mode, 0.5)
    nearly = collinear + np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1e-3 * STAGE_SPAN_RANK_TOLERANCE], [0.0, 0.0, 0.0]])
    with pytest.raises(SemiSyntheticTrajectoryError, match="no second direction"):
        edit_stage_configuration(nearly, mode, 0.5)


def test_rank_one_configuration_still_scales_and_translates() -> None:
    collinear = np.array([[0.0, 0.0], [1.0, 2.0], [3.0, 6.0]])
    scaled, _ = edit_stage_configuration(collinear, "magnitude", 0.5)
    shifted, effect = edit_stage_configuration(collinear, "translation", 1.0)
    assert np.allclose(scaled - scaled.mean(0), 1.5 * (collinear - collinear.mean(0)))
    assert np.allclose(shifted - collinear, effect["offset_length"] * np.array([1.0, 2.0]) / np.sqrt(5.0))


def test_blocks_outside_the_span_are_untouched(reference) -> None:
    null = generate_semisynthetic_trajectory(_span(construction=None, span_layers=None), reference=reference)
    for mode in ("magnitude", "orientation", "shape", "translation"):
        dataset = generate_semisynthetic_trajectory(
            _span(trajectory_mode=mode, group_effect_size=0.5), reference=reference
        )
        assert dataset.proteomics.equals(null.proteomics)
        pd.testing.assert_frame_equal(
            dataset.population_trajectories.layers["proteomics"], null.population_trajectories.layers["proteomics"]
        )


# --------------------------------------------------------------------------- #
# 2.4 The dataset is the null dataset plus the plant
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("mode", ["magnitude", "orientation", "shape", "translation"])
def test_stage_span_dataset_is_the_null_dataset_plus_the_plant(reference, mode) -> None:
    params = _span(trajectory_mode=mode, group_effect_size=0.5, group_stage_sizes=((11, 10, 28), (9, 10, 12)),
                   n_samples=600)
    null = generate_semisynthetic_trajectory(
        replace(params, trajectory_mode="none", group_effect_size=0.0, construction=None, span_layers=None),
        reference=reference,
    )
    dataset = generate_semisynthetic_trajectory(params, reference=reference)

    group_a = dataset.metadata["group"] == "A"
    group_b = ~group_a
    stages_b = dataset.metadata.loc[group_b, "stage"].to_numpy()
    assert dataset.metadata.equals(null.metadata)
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(dataset, layer).loc[group_a].equals(getattr(null, layer).loc[group_a])

    population = dataset.population_trajectories.layers
    null_population = null.population_trajectories.layers
    for layer in ("methylation", "expression"):
        plant = population[layer].loc["B"].to_numpy() - null_population[layer].loc["B"].to_numpy()
        assert np.abs(plant).max() > 0
        assert population[layer].loc["A"].equals(null_population[layer].loc["A"])
        if layer == "methylation":
            difference = logit(dataset.methylation.loc[group_b].to_numpy()) - logit(
                null.methylation.loc[group_b].to_numpy()
            )
            atol = 1e-8
        else:
            difference = dataset.expression.loc[group_b].to_numpy() - null.expression.loc[group_b].to_numpy()
            atol = 1e-10
        np.testing.assert_allclose(difference, plant[stages_b], atol=atol)


def test_stage_span_truth_records_the_construction(reference) -> None:
    dataset = generate_semisynthetic_trajectory(
        _span(trajectory_mode="orientation", group_effect_size=0.5), reference=reference
    )
    truth = dataset.truth
    assert truth["construction"] == "stage_span"
    assert truth["span_layers"] == list(TWO_BLOCKS)
    assert truth["transform"] == {"construction": "stage_span", "rotation_degrees": 45.0}
    block = truth["stage_span"]
    assert block["units"] == "null_standardized"
    assert block["mean_step"] == pytest.approx(block["path_length"] / 2)
    assert block["intended_geometry"]["angle"] == pytest.approx(45.0)
    assert block["intended_geometry"]["orientation_estimator"] == "net_displacement"
    assert len(block["plant_norm"]) == 3 and min(block["plant_norm"]) > 0
    # No surgery is recorded for a stage-span cell.
    assert not {"orientation_relocated", "censored"} & set(truth["transform"])


@pytest.mark.parametrize(("mode", "effect"), [("none", 0.5), ("orientation", 0.0)])
def test_null_stage_span_is_the_null_dataset(reference, mode, effect) -> None:
    null = generate_semisynthetic_trajectory(_span(construction=None, span_layers=None), reference=reference)
    dataset = generate_semisynthetic_trajectory(
        _span(trajectory_mode=mode, group_effect_size=effect), reference=reference
    )
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(dataset, layer).equals(getattr(null, layer))
    assert dataset.truth["transform"] == {}
    assert dataset.truth["stage_span"]["plant_norm"] == [0.0, 0.0, 0.0]
