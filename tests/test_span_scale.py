"""The planted-scale stage-span plant (``span_scale``)."""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace

import numpy as np
import pytest

import motco.simulations.semisynthetic as semisynthetic
from motco.simulations.evaluation import SimulationEvaluationParams
from motco.simulations.generator import logit
from motco.simulations.grid import _to_jsonable, make_simulation_cell, parameter_signature
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import (
    SemiSyntheticTrajectoryError,
    SemiSyntheticTrajectoryParams,
    edit_stage_configuration,
    generate_semisynthetic_trajectory,
    null_pooled_scale,
    planted_pooled_scale,
)

TWO_BLOCKS = ("methylation", "expression")
SEAAD_TABLE = ((11, 10, 28), (9, 10, 12))


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
        "trajectory_mode": "shape",
        "group_effect_size": 0.5,
        "span_scale": "planted",
    }
    return SemiSyntheticTrajectoryParams(**{**base, **overrides})


def _cell(generator: SemiSyntheticTrajectoryParams):
    return make_simulation_cell(
        phase="power_primary",
        generator_params=generator,
        evaluation_params=SimulationEvaluationParams(integration_method="pls", permutations=9),
    )


def _group_means(dataset) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray, np.ndarray]:
    layers = tuple(dataset.truth["span_layers"])
    population = dataset.population_trajectories.layers
    means_a = {layer: population[layer].loc["A"].to_numpy() for layer in layers}
    means_b = {layer: population[layer].loc["B"].to_numpy() for layer in layers}
    sizes = dataset.metadata.groupby(["group", "stage"]).size()
    return means_a, means_b, sizes.loc["A"].to_numpy(dtype=float), sizes.loc["B"].to_numpy(dtype=float)


def _planted_configurations(dataset, reference) -> tuple[np.ndarray, np.ndarray]:
    """A's and B's population stage means over the span, divided by the planted pooled SD.

    The SD is recomputed here from the population trajectories and the realized
    cell sizes, independently of the generator's fixed point.
    """

    layers = tuple(dataset.truth["span_layers"])
    means_a, means_b, sizes_a, sizes_b = _group_means(dataset)
    scale = planted_pooled_scale(means_a, means_b, sizes_a, sizes_b, layers=layers, reference=reference)
    a = np.hstack([means_a[layer] for layer in layers]) / scale
    b = np.hstack([means_b[layer] for layer in layers]) / scale
    return a, b


def _path_length(configuration: np.ndarray) -> float:
    return float(np.linalg.norm(np.diff(configuration, axis=0), axis=1).sum())


# --------------------------------------------------------------------------- #
# 2.1 Parameters and validation
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    ("overrides", "match"),
    [
        ({"span_scale": "sample"}, r"Unknown span_scale: 'sample'"),
        ({"span_scale": "null"}, r"Unknown span_scale: 'null'"),
        ({"span_scale": True}, r"Unknown span_scale: True"),
        ({"construction": None, "span_layers": None}, r"span_scale='planted' is set but construction is not"),
        ({"construction": "surgery", "span_layers": None}, r"span_scale='planted' is set but construction is not"),
        (
            {"trajectory_mode": "orientation", "span_alignment": 30.0},
            r"span_scale='planted' and span_alignment are not supported together",
        ),
    ],
)
def test_invalid_span_scale_is_rejected(reference, overrides, match) -> None:
    with pytest.raises(SemiSyntheticTrajectoryError, match=match):
        generate_semisynthetic_trajectory(_span(**overrides), reference=reference)


def test_absent_span_scale_is_the_default() -> None:
    assert SemiSyntheticTrajectoryParams(seed=1).span_scale is None
    assert _span(span_scale=None) == replace(_span(), span_scale=None)


# --------------------------------------------------------------------------- #
# 2.4 Identity without the setting
# --------------------------------------------------------------------------- #


def _dataset_digest(dataset) -> str:
    """Digest of the BLAS-free parts: indicators, metadata, population means and the truth without indicators."""

    digest = hashlib.sha256()
    for group in ("A", "B"):
        for layer in ("methylation", "expression", "proteomics"):
            digest.update(np.ascontiguousarray(dataset.truth["indicators"][group][layer], dtype=float).tobytes())
    digest.update(dataset.metadata.to_csv().encode())
    for frame in dataset.population_trajectories.layers.values():
        digest.update(frame.to_numpy().tobytes())
    truth = {key: value for key, value in dataset.truth.items() if key != "indicators"}
    digest.update(json.dumps(truth, sort_keys=True, default=str).encode())
    return digest.hexdigest()[:16]


# Captured at dc59d2f, before span_scale existed: the BLAS-free dataset digest
# (truth included), the cell id and the signature prefix of a null-scale
# stage-span PLS cell at n = 90, p_dmp 0.1, seed 17, methylation + expression.
_PRE_CHANGE = {
    ("magnitude", 0.1): ("c6afc1fa37b76f98", "power_primary-d6f52a79a876", "ba9329c8c7a08767"),
    ("orientation", 0.5): ("7eed0b84df3cb06b", "power_primary-698f92231645", "7ffc0586f9d02c8e"),
    ("shape", 0.5): ("d31c9e24e3dada0f", "power_primary-d9b75022e33e", "d3d8e96c90553b5d"),
    ("translation", 0.5): ("b3f715ea3db3d013", "power_primary-75ec0c961942", "0a22ef6d81e64217"),
}


@pytest.mark.parametrize(("mode", "effect"), sorted(_PRE_CHANGE))
def test_null_scale_stage_span_is_byte_identical_to_before(reference, mode, effect) -> None:
    params = SemiSyntheticTrajectoryParams(
        seed=17,
        n_samples=90,
        p_dmp=0.1,
        construction="stage_span",
        span_layers=TWO_BLOCKS,
        trajectory_mode=mode,
        group_effect_size=effect,
    )
    digest, cell_id, signature = _PRE_CHANGE[(mode, effect)]
    dataset = generate_semisynthetic_trajectory(params, reference=reference)
    assert _dataset_digest(dataset) == digest
    assert _cell(params).cell_id == cell_id
    assert parameter_signature(_cell(params))[:16] == signature

    explicit = replace(params, span_scale=None)
    assert _cell(explicit).cell_id == cell_id
    assert "span_scale" not in _to_jsonable(params)
    assert "span_scale" not in dataset.truth
    assert dataset.truth["stage_span"]["units"] == "null_standardized"
    assert not {"scale_iterations", "scale_residual"} & set(dataset.truth["stage_span"])

    planted = replace(params, span_scale="planted")
    assert _cell(planted).cell_id != cell_id
    assert parameter_signature(_cell(planted)) != parameter_signature(_cell(params))
    assert _to_jsonable(planted)["span_scale"] == "planted"


# --------------------------------------------------------------------------- #
# 2.2 The fixed point makes the edit exact in the planted coordinates
# --------------------------------------------------------------------------- #


_COLUMNS = {
    "n200": {"n_samples": 200},
    "n80": {"n_samples": 80},
    "seaad": {"n_samples": 600, "group_stage_sizes": SEAAD_TABLE},
    "three_blocks": {"n_samples": 200, "span_layers": None},
    "four_stages": {"n_samples": 200, "n_stages": 4},
}

_MODE_EFFECTS = (
    [("magnitude", effect) for effect in (0.1, 0.5, 1.0)]
    + [("orientation", effect) for effect in (0.1, 0.5, 1.0, 2.0)]
    + [("shape", effect) for effect in (0.1, 0.5, 1.0)]
    + [("translation", effect) for effect in (0.1, 0.5, 1.0)]
)


@pytest.mark.parametrize("column", sorted(_COLUMNS))
def test_planted_edit_is_exact_in_the_planted_coordinates(reference, column) -> None:
    for mode, effect in _MODE_EFFECTS:
        params = _span(seed=11, trajectory_mode=mode, group_effect_size=effect, **_COLUMNS[column])
        dataset = generate_semisynthetic_trajectory(params, reference=reference)
        a, b = _planted_configurations(dataset, reference)
        path = _path_length(a)
        tolerance = 1e-8 * path
        label = f"{column} {mode} e={effect}"

        # B's planted-standardized configuration is the mode's edit of A's.
        expected, _ = edit_stage_configuration(a, mode, effect)
        assert np.max(np.abs(b - expected)) <= tolerance, label

        block = dataset.truth["stage_span"]
        intended = block["intended_geometry"]
        assert block["path_length"] == pytest.approx(path, rel=1e-10), label
        assert block["scale_residual"] <= semisynthetic.PLANTED_SCALE_TOLERANCE, label
        assert 1 <= block["scale_iterations"] <= semisynthetic.PLANTED_SCALE_MAX_ITERATIONS, label
        target_delta = effect * path if mode == "magnitude" else 0.0
        assert abs(intended["delta"] - target_delta) <= tolerance, label
        assert abs(_path_length(b) - path - target_delta) <= tolerance, label
        target_angle = 90.0 * effect if mode == "orientation" else 0.0
        assert intended["angle"] == pytest.approx(target_angle, abs=1e-5), label
        if mode == "shape":
            assert intended["shape"] > 1e-3, label
        else:
            assert intended["shape"] <= 1e-8, label


def test_the_null_scale_edit_is_not_exact_in_the_planted_coordinates(reference) -> None:
    # The leak the setting removes: the null-scale bend shortens B under the planted SD.
    dataset = generate_semisynthetic_trajectory(_span(span_scale=None, n_samples=200), reference=reference)
    a, b = _planted_configurations(dataset, reference)
    assert (_path_length(b) - _path_length(a)) / _path_length(a) < -1e-3


def test_a_non_converging_fixed_point_fails_loudly(reference, monkeypatch) -> None:
    monkeypatch.setattr(semisynthetic, "PLANTED_SCALE_MAX_ITERATIONS", 3)
    with pytest.raises(
        SemiSyntheticTrajectoryError,
        match=r"did not converge for trajectory_mode='shape' at group_effect_size=0\.5; the last update still "
        r"moved it by [0-9.e-]+ \(relative\) after 3 iterations",
    ):
        generate_semisynthetic_trajectory(_span(), reference=reference)


# --------------------------------------------------------------------------- #
# 2.3 Truth
# --------------------------------------------------------------------------- #


def test_planted_truth_names_its_scale_and_convergence(reference) -> None:
    dataset = generate_semisynthetic_trajectory(_span(trajectory_mode="orientation"), reference=reference)
    truth = dataset.truth
    assert truth["span_scale"] == "planted"
    assert truth["construction"] == "stage_span"
    assert truth["transform"] == {"construction": "stage_span", "rotation_degrees": 45.0}
    block = truth["stage_span"]
    assert block["units"] == "planted_standardized"
    assert isinstance(block["scale_iterations"], int) and block["scale_iterations"] > 1
    assert 0.0 <= block["scale_residual"] <= semisynthetic.PLANTED_SCALE_TOLERANCE
    assert block["mean_step"] == pytest.approx(block["path_length"] / 2)

    a, b = _planted_configurations(dataset, reference)
    assert block["path_length"] == pytest.approx(_path_length(a), rel=1e-10)
    np.testing.assert_allclose(block["plant_norm"], np.linalg.norm(b - a, axis=1), rtol=1e-9)
    assert block["intended_geometry"]["path_lengths"]["B"] == pytest.approx(_path_length(b), rel=1e-10)


# --------------------------------------------------------------------------- #
# 2.5 The planted dataset is the null dataset plus the plant
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("mode", ["magnitude", "orientation", "shape", "translation"])
def test_planted_dataset_is_the_null_dataset_plus_the_plant(reference, mode) -> None:
    params = _span(trajectory_mode=mode, group_stage_sizes=SEAAD_TABLE, n_samples=600)
    null = generate_semisynthetic_trajectory(
        replace(
            params, trajectory_mode="none", group_effect_size=0.0, construction=None, span_layers=None, span_scale=None
        ),
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


@pytest.mark.parametrize(("mode", "effect"), [("none", 0.5), ("shape", 0.0), ("orientation", 0.0)])
def test_zero_effect_under_the_planted_scale_is_the_null_dataset(reference, mode, effect) -> None:
    null = generate_semisynthetic_trajectory(
        _span(construction=None, span_layers=None, span_scale=None, trajectory_mode="none", group_effect_size=0.0),
        reference=reference,
    )
    dataset = generate_semisynthetic_trajectory(
        _span(trajectory_mode=mode, group_effect_size=effect), reference=reference
    )
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(dataset, layer).equals(getattr(null, layer))
    assert dataset.metadata.equals(null.metadata)
    assert dataset.truth["transform"] == {}
    assert dataset.truth["stage_span"]["plant_norm"] == [0.0, 0.0, 0.0]


@pytest.mark.parametrize("column", ["n200", "seaad", "three_blocks"])
def test_planted_magnitude_is_the_null_scale_magnitude(reference, column) -> None:
    # Scaling about the centroid commutes with a per-feature scale, so the plant,
    # and every sampled value, is the null-scale construction's, bit for bit.
    planted = _span(trajectory_mode="magnitude", group_effect_size=0.05, **_COLUMNS[column])
    null_scale = replace(planted, span_scale=None)
    dataset = generate_semisynthetic_trajectory(planted, reference=reference)
    expected = generate_semisynthetic_trajectory(null_scale, reference=reference)
    for layer in ("methylation", "expression", "proteomics"):
        assert getattr(dataset, layer).to_numpy().tobytes() == getattr(expected, layer).to_numpy().tobytes()
    for layer, frame in dataset.population_trajectories.layers.items():
        assert frame.equals(expected.population_trajectories.layers[layer])
    assert dataset.truth["stage_span"]["scale_iterations"] == 2
    assert dataset.truth["stage_span"]["scale_residual"] == 0.0


def test_planted_pooled_scale_of_equal_groups_is_the_null_scale(reference) -> None:
    dataset = generate_semisynthetic_trajectory(
        _span(trajectory_mode="none", group_effect_size=0.0, group_stage_sizes=SEAAD_TABLE, n_samples=600),
        reference=reference,
    )
    means_a, _, sizes_a, sizes_b = _group_means(dataset)
    planted = planted_pooled_scale(means_a, means_a, sizes_a, sizes_b, layers=TWO_BLOCKS, reference=reference)
    null = null_pooled_scale(means_a, sizes_a + sizes_b, layers=TWO_BLOCKS, reference=reference)
    np.testing.assert_allclose(planted, null, rtol=1e-12)
