"""The opt-in group dispersion scale (``group_dispersion_scale``, ``generate_omics(noise_scale=)``)."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from motco.simulations.evaluation import SimulationEvaluationParams
from motco.simulations.generator import (
    GeneratorError,
    derive_coupled_indicators,
    generate_omics,
    logit,
    markov_indicators,
    omic_population_means,
    rev_logit,
)
from motco.simulations.grid import (
    _to_jsonable,
    derive_replicate_seed,
    enumerate_type_i_grid,
    make_simulation_cell,
    parameter_signature,
)
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import (
    SemiSyntheticTrajectoryError,
    SemiSyntheticTrajectoryParams,
    generate_semisynthetic_trajectory,
    group_dispersion_within_multiplier,
    null_pooled_scale,
    planted_pooled_scale,
)
from motco.simulations.study import (
    MatchedSeedPolicy,
    StudyConfig,
    dump_study_config,
    enumerate_study,
    load_study_config,
)
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

LAYERS = ("methylation", "expression", "proteomics")
TABLE = ((5, 6, 7), (8, 4, 9))
SCALE_AXIS = "generator.group_dispersion_scale"
RATIO_AXIS = "generator.group_ratio"


@pytest.fixture(scope="module")
def reference():
    return load_reference()


def _params(**overrides: Any) -> SemiSyntheticTrajectoryParams:
    base: dict[str, Any] = {"seed": 31, "n_samples": 60, "p_dmp": 0.1}
    return SemiSyntheticTrajectoryParams(**{**base, **overrides})


def _analysed(dataset, layer: str) -> np.ndarray:
    values = getattr(dataset, layer).to_numpy()
    return logit(values) if layer == "methylation" else values


def _row_means(dataset, layer: str, group: str) -> np.ndarray:
    """Each sample's exact population cell mean, in analysed units."""

    frame = dataset.population_trajectories.layers[layer].loc[group]
    stages = dataset.metadata.loc[dataset.metadata["group"] == group, "stage"].astype(str)
    return frame.loc[stages.to_list()].to_numpy()


def _assert_truth_equal(left: Any, right: Any) -> None:
    if isinstance(left, dict):
        assert isinstance(right, dict) and set(left) == set(right)
        for key in left:
            _assert_truth_equal(left[key], right[key])
    elif isinstance(left, np.ndarray):
        assert isinstance(right, np.ndarray) and left.tobytes() == right.tobytes()
    else:
        assert left == right


# --------------------------------------------------------------------------- #
# 1.1 generate_omics(noise_scale=)
# --------------------------------------------------------------------------- #


def _indicators(reference, n_cells: int, seed: int = 5) -> dict[str, Any]:
    methyl = markov_indicators(np.random.default_rng(seed), reference.n_cpg, n_cells, 0.1)
    expr, protein = derive_coupled_indicators(methyl, reference)
    return {"indicators_methyl": methyl, "indicators_expr": expr, "indicators_protein": protein}


@pytest.mark.parametrize("scale", [0.5, 1.7, 3.0])
def test_noise_scale_stretches_the_draw_about_the_cell_means(reference, scale) -> None:
    indicators = _indicators(reference, 3)
    cell_sizes = [5, 0, 7]  # an empty cell keeps the stacking order honest
    rng_plain, rng_scaled = np.random.default_rng(11), np.random.default_rng(11)
    plain = generate_omics(cell_sizes=cell_sizes, rng=rng_plain, reference=reference, **indicators)
    scaled = generate_omics(
        cell_sizes=cell_sizes, rng=rng_scaled, reference=reference, noise_scale=scale, **indicators
    )
    means = omic_population_means(reference=reference, **indicators)
    rows = {layer: np.repeat(means[layer], cell_sizes, axis=0) for layer in LAYERS}

    plain_m = logit(plain.methylation)
    expected_m = rows["methylation"] + scale * (plain_m - rows["methylation"])
    np.testing.assert_allclose(scaled.methylation, rev_logit(expected_m), rtol=0, atol=1e-10)
    for layer in ("expression", "proteomics"):
        expected = rows[layer] + scale * (getattr(plain, layer) - rows[layer])
        np.testing.assert_allclose(getattr(scaled, layer), expected, rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(scaled.cell_ids, plain.cell_ids)
    # No random numbers drawn: both generators stand at the same stream position.
    assert rng_scaled.random() == rng_plain.random()


def test_noise_scale_keeps_methylation_inside_the_unit_interval(reference) -> None:
    indicators = _indicators(reference, 2)
    generated = generate_omics(
        cell_sizes=[40, 40], rng=np.random.default_rng(2), reference=reference, noise_scale=4.0, **indicators
    )
    assert np.all(generated.methylation > 0.0)
    assert np.all(generated.methylation < 1.0)


def test_noise_scale_multiplies_the_dispersion_in_m_values(reference) -> None:
    indicators = _indicators(reference, 1)
    plain = generate_omics(cell_sizes=[3000], rng=np.random.default_rng(9), reference=reference, **indicators)
    scaled = generate_omics(
        cell_sizes=[3000], rng=np.random.default_rng(9), reference=reference, noise_scale=2.0, **indicators
    )
    # Features the clipped logit cannot invert (values within 1e-6 of 0 or 1) are left out.
    inside = np.all((scaled.methylation > 1e-6) & (scaled.methylation < 1 - 1e-6), axis=0)
    assert inside.mean() > 0.9
    ratio = logit(scaled.methylation).std(axis=0) / logit(plain.methylation).std(axis=0)
    np.testing.assert_allclose(ratio[inside], 2.0, rtol=1e-6)
    for layer in ("expression", "proteomics"):
        np.testing.assert_allclose(getattr(scaled, layer).std(axis=0) / getattr(plain, layer).std(axis=0), 2.0)


def test_absent_noise_scale_is_the_unscaled_call(reference) -> None:
    indicators = _indicators(reference, 2)
    default = generate_omics(cell_sizes=[4, 6], rng=np.random.default_rng(3), reference=reference, **indicators)
    explicit = generate_omics(
        cell_sizes=[4, 6], rng=np.random.default_rng(3), reference=reference, noise_scale=None, **indicators
    )
    for layer in LAYERS:
        assert getattr(default, layer).tobytes() == getattr(explicit, layer).tobytes()


@pytest.mark.parametrize("scale", [0, 0.0, -1.0, float("nan"), float("inf"), True, "2"])
def test_invalid_noise_scale_is_rejected(reference, scale) -> None:
    with pytest.raises(GeneratorError, match="noise_scale"):
        generate_omics(
            cell_sizes=[2, 2], rng=np.random.default_rng(1), reference=reference, noise_scale=scale,
            **_indicators(reference, 2),
        )


# --------------------------------------------------------------------------- #
# 1.2 SemiSyntheticTrajectoryParams.group_dispersion_scale
# --------------------------------------------------------------------------- #


CASES = {
    "none": {},
    "orientation_surgery": {"trajectory_mode": "orientation", "group_effect_size": 0.3},
    "magnitude_surgery": {"trajectory_mode": "magnitude", "group_effect_size": 0.5},
    "size_table": {"n_samples": 600, "group_stage_sizes": TABLE},
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_second_group_is_stretched_about_its_population_cell_means(reference, case) -> None:
    overrides = CASES[case]
    plain = generate_semisynthetic_trajectory(_params(**overrides), reference=reference)
    scaled = generate_semisynthetic_trajectory(_params(**overrides, group_dispersion_scale=2.0), reference=reference)
    groups = plain.metadata["group"].to_numpy()
    assert scaled.metadata.equals(plain.metadata)
    for layer in LAYERS:
        plain_values = getattr(plain, layer).to_numpy()
        scaled_values = getattr(scaled, layer).to_numpy()
        # Group A is untouched, byte for byte.
        assert scaled_values[groups == "A"].tobytes() == plain_values[groups == "A"].tobytes()
        rows = _row_means(plain, layer, "B")
        expected = rows + 2.0 * (_analysed(plain, layer)[groups == "B"] - rows)
        if layer == "methylation":
            np.testing.assert_allclose(scaled_values[groups == "B"], rev_logit(expected), rtol=0, atol=1e-10)
        else:
            np.testing.assert_allclose(scaled_values[groups == "B"], expected, rtol=1e-12, atol=1e-12)
    # The means and the indicators do not move.
    for layer in LAYERS:
        assert scaled.population_trajectories.layers[layer].equals(plain.population_trajectories.layers[layer])
    _assert_truth_equal(scaled.truth["indicators"], plain.truth["indicators"])
    assert scaled.truth["deltas"] == plain.truth["deltas"]


def test_unit_scale_is_absent(reference) -> None:
    assert SemiSyntheticTrajectoryParams(seed=1).group_dispersion_scale is None
    for one in (1, 1.0, np.float64(1.0)):
        params = _params(group_dispersion_scale=one)
        assert params.group_dispersion_scale is None
        assert params == _params()
    plain = generate_semisynthetic_trajectory(_params(), reference=reference)
    unit = generate_semisynthetic_trajectory(_params(group_dispersion_scale=1.0), reference=reference)
    for layer in LAYERS:
        assert getattr(unit, layer).to_numpy().tobytes() == getattr(plain, layer).to_numpy().tobytes()
    _assert_truth_equal(unit.truth, plain.truth)
    assert "group_dispersion" not in plain.truth


def test_truth_names_the_scale_its_covariance_multiplier_and_group(reference) -> None:
    dataset = generate_semisynthetic_trajectory(_params(group_dispersion_scale=1.5), reference=reference)
    assert dataset.truth["group_dispersion"] == {"scale": 1.5, "covariance_multiplier": 2.25, "group": "B"}
    relabelled = generate_semisynthetic_trajectory(
        _params(group_dispersion_scale=1.5, group_labels=("F", "M")), reference=reference
    )
    assert relabelled.truth["group_dispersion"]["group"] == "M"


@pytest.mark.parametrize("scale", [0, -1.5, float("nan"), float("inf"), True, "2"])
def test_invalid_group_dispersion_scale_is_rejected(reference, scale) -> None:
    with pytest.raises(SemiSyntheticTrajectoryError, match="group_dispersion_scale"):
        generate_semisynthetic_trajectory(_params(group_dispersion_scale=scale), reference=reference)


# --------------------------------------------------------------------------- #
# 1.3 Stage-span scales with one group wider
# --------------------------------------------------------------------------- #


def _span(**overrides: Any) -> SemiSyntheticTrajectoryParams:
    base: dict[str, Any] = {
        "seed": 23,
        "n_samples": 120,
        "p_dmp": 0.1,
        "group_ratio": 0.65,
        "construction": "stage_span",
        "trajectory_mode": "translation",
        "group_effect_size": 1.0,
        "span_scale": "planted",
        "group_dispersion_scale": 2.0,
    }
    return SemiSyntheticTrajectoryParams(**{**base, **overrides})


def _group_means(dataset):
    population = dataset.population_trajectories.layers
    means_a = {layer: population[layer].loc["A"].to_numpy() for layer in LAYERS}
    means_b = {layer: population[layer].loc["B"].to_numpy() for layer in LAYERS}
    sizes = dataset.metadata.groupby(["group", "stage"]).size()
    return means_a, means_b, sizes.loc["A"].to_numpy(dtype=float), sizes.loc["B"].to_numpy(dtype=float)


def test_within_multiplier_weights_each_group_by_its_realized_size() -> None:
    assert group_dispersion_within_multiplier(None, [1, 2], [3, 4]) is None
    assert group_dispersion_within_multiplier(2.0, [10, 10, 10], [5, 5, 5]) == pytest.approx((30 + 4 * 15) / 45)
    assert group_dispersion_within_multiplier(1.5, [7, 8, 9], [3, 2, 1]) == pytest.approx((24 + 2.25 * 6) / 30)


def test_pooled_scales_carry_the_within_cell_factor(reference) -> None:
    dataset = generate_semisynthetic_trajectory(
        _params(trajectory_mode="none", group_ratio=0.8, group_dispersion_scale=2.0), reference=reference
    )
    means_a, _, sizes_a, sizes_b = _group_means(dataset)
    multiplier = group_dispersion_within_multiplier(2.0, sizes_a, sizes_b)
    assert multiplier is not None
    n_a, n_b = sizes_a.sum(), sizes_b.sum()
    assert n_a != n_b
    assert multiplier == pytest.approx((n_a + 4.0 * n_b) / (n_a + n_b))
    weights = sizes_a + sizes_b
    shares = weights / weights.sum()
    diagonals = {
        "methylation": np.diag(reference.cov_M),
        "expression": np.diag(reference.cov_expr),
        "proteomics": np.diag(reference.cov_protein),
    }
    direct = []
    for layer in LAYERS:
        centred = means_a[layer] - shares @ means_a[layer]
        direct.append(np.sqrt(multiplier * diagonals[layer] + shares @ centred**2))
    expected = np.concatenate(direct)
    null = null_pooled_scale(means_a, weights, layers=LAYERS, reference=reference, within_multiplier=multiplier)
    planted = planted_pooled_scale(
        means_a, means_a, sizes_a, sizes_b, layers=LAYERS, reference=reference, within_multiplier=multiplier
    )
    np.testing.assert_allclose(null, expected, rtol=1e-12)
    np.testing.assert_allclose(planted, expected, rtol=1e-12)
    # Absent, both scales are the reference-variance scales.
    np.testing.assert_array_equal(
        null_pooled_scale(means_a, weights, layers=LAYERS, reference=reference, within_multiplier=None),
        null_pooled_scale(means_a, weights, layers=LAYERS, reference=reference),
    )


@pytest.mark.parametrize("span_scale", [None, "planted"])
@pytest.mark.parametrize("null", [{"trajectory_mode": "none"}, {"group_effect_size": 0.0}])
def test_null_plant_with_a_scale_is_the_null_dataset(reference, span_scale, null) -> None:
    plant = generate_semisynthetic_trajectory(_span(span_scale=span_scale, **null), reference=reference)
    plain = generate_semisynthetic_trajectory(
        _span(construction=None, span_scale=None, trajectory_mode="none", group_effect_size=0.0),
        reference=reference,
    )
    for layer in LAYERS:
        assert getattr(plant, layer).to_numpy().tobytes() == getattr(plain, layer).to_numpy().tobytes()


def test_planted_translation_with_a_scale_moves_every_stage_by_one_vector(reference) -> None:
    dataset = generate_semisynthetic_trajectory(_span(), reference=reference)
    assert dataset.truth["stage_span"]["units"] == "planted_standardized"
    assert dataset.truth["stage_span"]["scale_residual"] <= 1e-12
    means_a, means_b, sizes_a, sizes_b = _group_means(dataset)
    shift = np.hstack([means_b[layer] - means_a[layer] for layer in LAYERS])
    np.testing.assert_allclose(shift, np.broadcast_to(shift[0], shift.shape), rtol=0, atol=1e-12)
    # Exact in the planted scale with each group's own noise: |B − A| is the offset there.
    multiplier = group_dispersion_within_multiplier(2.0, sizes_a, sizes_b)
    scale = planted_pooled_scale(
        means_a, means_b, sizes_a, sizes_b, layers=LAYERS, reference=reference, within_multiplier=multiplier
    )
    offset = dataset.truth["transform"]["offset_length"]
    np.testing.assert_allclose(np.linalg.norm(shift / scale, axis=1), offset, rtol=1e-9)
    # The reference-variance scale would state the offset in a scale the data do not have.
    unaware = planted_pooled_scale(means_a, means_b, sizes_a, sizes_b, layers=LAYERS, reference=reference)
    assert abs(float(np.linalg.norm(shift[0] / unaware)) / offset - 1.0) > 0.05


def test_planted_scale_matches_a_large_sample(reference) -> None:
    params = _span(n_samples=600, group_stage_sizes=((2000, 2000, 2000), (1000, 1000, 1000)), group_ratio=0.5)
    dataset = generate_semisynthetic_trajectory(params, reference=reference)
    means_a, means_b, sizes_a, sizes_b = _group_means(dataset)
    multiplier = group_dispersion_within_multiplier(2.0, sizes_a, sizes_b)
    expected = planted_pooled_scale(
        means_a, means_b, sizes_a, sizes_b, layers=LAYERS, reference=reference, within_multiplier=multiplier
    )
    observed = np.concatenate([_analysed(dataset, layer).std(axis=0) for layer in LAYERS])
    relative = np.abs(observed / expected - 1.0)
    # The SD of a sample SD is about SD / sqrt(2n): 0.5 % at n = 9,000.
    assert np.median(relative) < 0.01
    assert relative.max() < 0.04
    unaware = planted_pooled_scale(means_a, means_b, sizes_a, sizes_b, layers=LAYERS, reference=reference)
    assert np.median(np.abs(observed / unaware - 1.0)) > 0.1


# --------------------------------------------------------------------------- #
# 1.4 Identifiers
# --------------------------------------------------------------------------- #


def _cell(generator: SemiSyntheticTrajectoryParams):
    return make_simulation_cell(
        phase="power_primary",
        generator_params=generator,
        evaluation_params=SimulationEvaluationParams(integration_method="concat", permutations=0),
    )


def test_a_set_scale_changes_the_cell_identity() -> None:
    plain, scaled = _cell(_params()), _cell(_params(group_dispersion_scale=1.5))
    assert plain.cell_id != scaled.cell_id
    assert parameter_signature(plain) != parameter_signature(scaled)
    assert "group_dispersion_scale" not in _to_jsonable(plain.generator_params)
    assert _to_jsonable(scaled.generator_params)["group_dispersion_scale"] == 1.5
    unit = _cell(_params(group_dispersion_scale=1))
    assert (unit.cell_id, parameter_signature(unit)) == (plain.cell_id, parameter_signature(plain))


def test_type_i_grid_keeps_the_scale() -> None:
    grid = enumerate_type_i_grid(
        baseline_generator_params=_params(
            trajectory_mode="orientation", group_effect_size=0.5, group_dispersion_scale=2.0
        ),
        evaluation_params=SimulationEvaluationParams(integration_method="concat", permutations=0),
        axes={"generator.group_ratio": (0.5, 0.8)},
    )
    assert len(grid.cells) == 2
    for cell in grid.cells:
        assert cell.generator_params.group_effect_size == 0.0
        assert cell.generator_params.group_dispersion_scale == 2.0


def _config(**overrides: Any) -> StudyConfig:
    defaults: dict[str, Any] = {
        "generator": SemiSyntheticTrajectoryParams(
            seed=2, trajectory_mode="magnitude", n_samples=60, p_dmp=0.1, group_dispersion_scale=1.5
        ),
        "evaluation": SimulationEvaluationParams(integration_method="concat", permutations=0, seed=3),
        "trajectory_modes": ("magnitude", "orientation"),
        "effect_sizes": (0.0, 0.5),
        "n_replicates": 2,
        "base_seed": 100,
        "matched_seeds": MatchedSeedPolicy(enabled=True, primary_family="fam"),
    }
    defaults.update(overrides)
    return StudyConfig(**defaults)


def _identity(config: StudyConfig) -> dict[str, tuple[str, list[int]]]:
    return {
        cell.cell_id: (parameter_signature(cell), [derive_replicate_seed(cell, i) for i in range(cell.n_replicates)])
        for cell in enumerate_study(config).cells
    }


def test_study_cells_carry_the_scale_anchors_and_negative_controls_included() -> None:
    config = _config()
    cells = enumerate_study(config).cells
    assert any(cell.phase.startswith("type_i") for cell in cells)
    assert any(cell.metadata.get("negative_control") for cell in cells)
    assert any(cell.metadata.get("zero_effect_anchor") for cell in cells)
    for cell in cells:
        assert cell.generator_params.group_dispersion_scale == 1.5
    without = _config(generator=replace(config.generator, group_dispersion_scale=None))
    assert not set(_identity(config)) & set(_identity(without))


def _write_config(tmp_path: Path, generator_extra: dict[str, Any], design_axes=None) -> Path:
    payload: dict[str, Any] = {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "n_samples": 60, "p_dmp": 0.1, **generator_extra},
        "evaluation": {"integration_method": "concat", "permutations": 0, "seed": 3},
        "trajectory_modes": ["magnitude", "orientation"],
        "effect_sizes": [0.0, 0.5],
        "n_replicates": 2,
        "base_seed": 100,
        "matched_seeds": {"enabled": True, "primary_family": "fam"},
    }
    if design_axes is not None:
        payload["design_grid"] = {"axes": design_axes}
    path = tmp_path / f"config-{len(list(tmp_path.iterdir()))}.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


@pytest.mark.parametrize("spelling", [None, 1, 1.0])
def test_null_or_unit_scale_equals_an_omitted_key(tmp_path: Path, spelling) -> None:
    omitted = load_study_config(_write_config(tmp_path, {}))
    spelled = load_study_config(_write_config(tmp_path, {"group_dispersion_scale": spelling}))
    assert spelled == omitted
    assert _identity(spelled) == _identity(omitted)
    dumped = tmp_path / "dumped.json"
    dump_study_config(spelled, dumped)
    assert "group_dispersion_scale" not in json.loads(dumped.read_text(encoding="utf-8"))["generator"]


def test_a_set_scale_loads_enumerates_and_dumps(tmp_path: Path) -> None:
    config = load_study_config(_write_config(tmp_path, {"group_dispersion_scale": 1.5}))
    assert config.generator.group_dispersion_scale == 1.5
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert json.loads(dumped.read_text(encoding="utf-8"))["generator"]["group_dispersion_scale"] == 1.5
    assert load_study_config(dumped) == config
    assert _identity(load_study_config(dumped)) == _identity(config)


def test_scale_as_a_design_axis_crossed_with_the_group_ratio(tmp_path: Path) -> None:
    path = _write_config(tmp_path, {}, design_axes={SCALE_AXIS: [None, 2.0], RATIO_AXIS: [0.5, 0.8]})
    config = load_study_config(path)
    points = config.design_points()
    assert len(points) == 4
    design = [cell for cell in enumerate_study(config).cells if cell.phase == DESIGN_PHASE]
    by_point: dict[tuple[Any, Any], set[Any]] = {}
    for cell in design:
        point = cell.metadata[DESIGN_POINT_KEY]
        params = cell.generator_params
        assert params.group_dispersion_scale == point[SCALE_AXIS]
        assert params.group_ratio == point[RATIO_AXIS]
        by_point.setdefault((point[SCALE_AXIS], point[RATIO_AXIS]), set()).add(cell.cell_id)
    assert set(by_point) == {(None, 0.8), (2.0, 0.5), (2.0, 0.8)}
    anchors = {}
    for cell in design:
        if cell.metadata.get("zero_effect_anchor"):
            point = cell.metadata[DESIGN_POINT_KEY]
            anchors[(point[SCALE_AXIS], point[RATIO_AXIS])] = cell.generator_params
    assert anchors[(None, 0.8)] != anchors[(2.0, 0.8)]
    assert replace(anchors[(None, 0.8)], group_dispersion_scale=2.0) == anchors[(2.0, 0.8)]
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert load_study_config(dumped) == config
