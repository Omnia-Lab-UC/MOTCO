"""The opt-in stage dispersion scale (``stage_dispersion_scale``, ``generate_omics(noise_scale=[...])``)."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import motco.simulations.semisynthetic as semisynthetic
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
    cell_dispersion_within_multiplier,
    dispersion_within_multiplier,
    edit_stage_configuration,
    generate_semisynthetic_trajectory,
    group_dispersion_within_multiplier,
    null_pooled_scale,
    planted_pooled_scale,
)
from motco.simulations.study import (
    MatchedSeedPolicy,
    StudyConfig,
    StudyConfigError,
    dump_study_config,
    enumerate_study,
    load_study_config,
)
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

LAYERS = ("methylation", "expression", "proteomics")
SEAAD_TABLE = ((11, 10, 28), (9, 10, 12))
STAGES = (1.12, 1.278, 0.802)
SCALE_AXIS = "generator.stage_dispersion_scale"


@pytest.fixture(scope="module")
def reference():
    return load_reference()


def _params(**overrides: Any) -> SemiSyntheticTrajectoryParams:
    base: dict[str, Any] = {"seed": 31, "n_samples": 60, "p_dmp": 0.1}
    return SemiSyntheticTrajectoryParams(**{**base, **overrides})


def _analysed(dataset, layer: str) -> np.ndarray:
    values = getattr(dataset, layer).to_numpy()
    return logit(values) if layer == "methylation" else values


def _row_means(dataset, layer: str) -> np.ndarray:
    """Each sample's exact population cell mean, in analysed units, in row order."""

    frame = dataset.population_trajectories.layers[layer]
    keys = list(zip(dataset.metadata["group"].astype(str), dataset.metadata["stage"].astype(str), strict=True))
    return frame.loc[keys].to_numpy()


def _row_multipliers(dataset, stage_scales, group_scale: float | None = None) -> np.ndarray:
    stages = dataset.metadata["stage"].astype(int).to_numpy()
    factors = np.asarray(stage_scales, dtype=float)[stages]
    if group_scale is not None:
        factors = factors * np.where(dataset.metadata["group"].to_numpy() == "B", group_scale, 1.0)
    return factors[:, None]


def _assert_truth_equal(left: Any, right: Any) -> None:
    if isinstance(left, dict):
        assert isinstance(right, dict) and set(left) == set(right)
        for key in left:
            _assert_truth_equal(left[key], right[key])
    elif isinstance(left, np.ndarray):
        assert isinstance(right, np.ndarray) and left.tobytes() == right.tobytes()
    else:
        assert left == right


def _assert_stretched(scaled, plain, factors: np.ndarray) -> None:
    assert scaled.metadata.equals(plain.metadata)
    for layer in LAYERS:
        rows = _row_means(plain, layer)
        expected = rows + factors * (_analysed(plain, layer) - rows)
        observed = getattr(scaled, layer).to_numpy()
        if layer == "methylation":
            np.testing.assert_allclose(observed, rev_logit(expected), rtol=0, atol=1e-10)
        else:
            np.testing.assert_allclose(observed, expected, rtol=1e-12, atol=1e-12)


# --------------------------------------------------------------------------- #
# 1.2 generate_omics(noise_scale=) with one scale per cell
# --------------------------------------------------------------------------- #


def _indicators(reference, n_cells: int, seed: int = 5) -> dict[str, Any]:
    methyl = markov_indicators(np.random.default_rng(seed), reference.n_cpg, n_cells, 0.1)
    expr, protein = derive_coupled_indicators(methyl, reference)
    return {"indicators_methyl": methyl, "indicators_expr": expr, "indicators_protein": protein}


@pytest.mark.parametrize("scales", [(0.5, 2.0, 1.3), [1.7, 0.9, 3.0], np.array([2.0, 1.0, 0.25])])
def test_per_cell_scales_stretch_each_cell_about_its_mean(reference, scales) -> None:
    indicators = _indicators(reference, 3)
    cell_sizes = [5, 0, 7]  # an empty cell keeps the stacking order honest
    rng_plain, rng_scaled = np.random.default_rng(11), np.random.default_rng(11)
    plain = generate_omics(cell_sizes=cell_sizes, rng=rng_plain, reference=reference, **indicators)
    scaled = generate_omics(
        cell_sizes=cell_sizes, rng=rng_scaled, reference=reference, noise_scale=scales, **indicators
    )
    means = omic_population_means(reference=reference, **indicators)
    rows = {layer: np.repeat(means[layer], cell_sizes, axis=0) for layer in LAYERS}
    factors = np.repeat(np.asarray(scales, dtype=float), cell_sizes)[:, None]

    plain_m = logit(plain.methylation)
    expected_m = rows["methylation"] + factors * (plain_m - rows["methylation"])
    np.testing.assert_allclose(scaled.methylation, rev_logit(expected_m), rtol=0, atol=1e-10)
    for layer in ("expression", "proteomics"):
        expected = rows[layer] + factors * (getattr(plain, layer) - rows[layer])
        np.testing.assert_allclose(getattr(scaled, layer), expected, rtol=1e-12, atol=1e-12)
    np.testing.assert_array_equal(scaled.cell_ids, plain.cell_ids)
    # No random numbers drawn: both generators stand at the same stream position.
    assert rng_scaled.random() == rng_plain.random()


def test_per_cell_scales_multiply_each_cells_dispersion(reference) -> None:
    indicators = _indicators(reference, 2)
    plain = generate_omics(cell_sizes=[2000, 2000], rng=np.random.default_rng(9), reference=reference, **indicators)
    scaled = generate_omics(
        cell_sizes=[2000, 2000], rng=np.random.default_rng(9), reference=reference, noise_scale=[0.5, 2.0],
        **indicators,
    )
    for cell, scale in ((0, 0.5), (1, 2.0)):
        rows = scaled.cell_ids == cell
        for layer in ("expression", "proteomics"):
            ratio = getattr(scaled, layer)[rows].std(axis=0) / getattr(plain, layer)[rows].std(axis=0)
            np.testing.assert_allclose(ratio, scale, rtol=1e-9)


@pytest.mark.parametrize("scale", [0.5, 1.0, 1.7])
def test_equal_per_cell_scales_are_the_single_scale(reference, scale) -> None:
    indicators = _indicators(reference, 3)
    single = generate_omics(
        cell_sizes=[4, 0, 6], rng=np.random.default_rng(3), reference=reference, noise_scale=scale, **indicators
    )
    per_cell = generate_omics(
        cell_sizes=[4, 0, 6], rng=np.random.default_rng(3), reference=reference, noise_scale=(scale,) * 3,
        **indicators,
    )
    for layer in LAYERS:
        assert getattr(per_cell, layer).tobytes() == getattr(single, layer).tobytes()


@pytest.mark.parametrize(
    ("scales", "match"),
    [
        ((1.0, 2.0), r"one scale per cell \(3\); got 2"),
        ((1.0, 2.0, 3.0, 4.0), r"one scale per cell \(3\); got 4"),
        ((1.0, 0.0, 2.0), "positive finite numbers"),
        ((1.0, -1.0, 2.0), "positive finite numbers"),
        ((1.0, float("nan"), 2.0), "positive finite numbers"),
        ((1.0, float("inf"), 2.0), "positive finite numbers"),
        ((1.0, True, 2.0), "positive finite numbers"),
        ((1.0, "2", 2.0), "positive finite numbers"),
    ],
)
def test_invalid_per_cell_scales_are_rejected_before_any_draw(reference, scales, match) -> None:
    rng = np.random.default_rng(1)
    with pytest.raises(GeneratorError, match=rf"noise_scale must .*{match}"):
        generate_omics(
            cell_sizes=[2, 2, 2], rng=rng, reference=reference, noise_scale=scales, **_indicators(reference, 3)
        )
    assert rng.random() == np.random.default_rng(1).random()


# --------------------------------------------------------------------------- #
# 1.3 SemiSyntheticTrajectoryParams.stage_dispersion_scale
# --------------------------------------------------------------------------- #


CASES = {
    "none": {},
    "orientation_surgery": {"trajectory_mode": "orientation", "group_effect_size": 0.3},
    "magnitude_surgery": {"trajectory_mode": "magnitude", "group_effect_size": 0.5},
    "size_table": {"n_samples": 600, "group_stage_sizes": SEAAD_TABLE},
    "four_stages": {"n_stages": 4, "stage_dispersion_scale": (1.3, 0.8, 1.0, 2.0)},
    "stage_span": {
        "construction": "stage_span",
        "trajectory_mode": "translation",
        "group_effect_size": 0.7,
        "group_stage_sizes": SEAAD_TABLE,
        "n_samples": 600,
    },
}


@pytest.mark.parametrize("case", sorted(CASES))
def test_dataset_is_the_unscaled_one_stretched_by_stage(reference, case) -> None:
    overrides = dict(CASES[case])
    stages = overrides.pop("stage_dispersion_scale", STAGES)
    if case == "stage_span":
        # Under stage_span the plant's scale reads the noise; the null dataset does not.
        overrides.update(group_effect_size=0.0)
    plain = generate_semisynthetic_trajectory(_params(**overrides), reference=reference)
    scaled = generate_semisynthetic_trajectory(_params(**overrides, stage_dispersion_scale=stages), reference=reference)
    _assert_stretched(scaled, plain, _row_multipliers(plain, stages))
    # The means and the indicators do not move.
    for layer in LAYERS:
        assert scaled.population_trajectories.layers[layer].equals(plain.population_trajectories.layers[layer])
    _assert_truth_equal(scaled.truth["indicators"], plain.truth["indicators"])
    assert scaled.truth["deltas"] == plain.truth["deltas"]


@pytest.mark.parametrize("case", ["none", "orientation_surgery", "size_table"])
def test_it_composes_with_the_group_scale_as_a_product(reference, case) -> None:
    overrides = CASES[case]
    plain = generate_semisynthetic_trajectory(_params(**overrides), reference=reference)
    both = generate_semisynthetic_trajectory(
        _params(**overrides, stage_dispersion_scale=STAGES, group_dispersion_scale=1.5), reference=reference
    )
    _assert_stretched(both, plain, _row_multipliers(plain, STAGES, group_scale=1.5))
    assert both.truth["group_dispersion"]["scale"] == 1.5
    assert both.truth["stage_dispersion"]["scale"] == list(STAGES)


def test_equal_stage_multipliers_are_the_uniform_stretch(reference) -> None:
    plain = generate_semisynthetic_trajectory(_params(), reference=reference)
    scaled = generate_semisynthetic_trajectory(_params(stage_dispersion_scale=(2.0, 2.0, 2.0)), reference=reference)
    _assert_stretched(scaled, plain, np.full((len(plain.metadata), 1), 2.0))


def test_truth_names_each_stage_multiplier_and_its_square(reference) -> None:
    dataset = generate_semisynthetic_trajectory(_params(stage_dispersion_scale=[1.5, 1, 0.5]), reference=reference)
    assert dataset.truth["stage_dispersion"] == {
        "scale": [1.5, 1.0, 0.5],
        "covariance_multiplier": [2.25, 1.0, 0.25],
    }


def test_all_ones_is_absent(reference) -> None:
    assert SemiSyntheticTrajectoryParams(seed=1).stage_dispersion_scale is None
    for ones in ((1, 1, 1), [1.0, 1.0, 1.0], np.ones(3)):
        params = _params(stage_dispersion_scale=ones)
        assert params.stage_dispersion_scale is None
        assert params == _params()
    plain = generate_semisynthetic_trajectory(_params(), reference=reference)
    unit = generate_semisynthetic_trajectory(_params(stage_dispersion_scale=(1, 1, 1)), reference=reference)
    for layer in LAYERS:
        assert getattr(unit, layer).to_numpy().tobytes() == getattr(plain, layer).to_numpy().tobytes()
    _assert_truth_equal(unit.truth, plain.truth)
    assert "stage_dispersion" not in plain.truth


def test_a_sequence_is_held_as_a_tuple_of_floats() -> None:
    params = _params(stage_dispersion_scale=[2, 1, np.float32(0.5)])
    assert params.stage_dispersion_scale == (2.0, 1.0, 0.5)
    assert all(type(value) is float for value in params.stage_dispersion_scale)
    assert params == _params(stage_dispersion_scale=(2.0, 1.0, 0.5))


@pytest.mark.parametrize(
    ("scales", "match"),
    [
        ((1.0, 2.0), r"one multiplier per stage \(n_stages=3\); got 2"),
        ((1.0, 2.0, 3.0, 0.5), r"one multiplier per stage \(n_stages=3\); got 4"),
        ((), r"one multiplier per stage \(n_stages=3\); got 0"),
        ((1.0, 0.0, 2.0), r"positive finite numbers; got \(1.0, 0.0, 2.0\)"),
        ((1.0, -2.0, 2.0), "positive finite numbers"),
        ((1.0, float("nan"), 2.0), "positive finite numbers"),
        ((1.0, float("inf"), 2.0), "positive finite numbers"),
        ((1.0, True, 2.0), "positive finite numbers"),
        ((1.0, "2", 2.0), "positive finite numbers"),
        (2.0, "be a sequence"),
        ("2", "be a sequence"),
    ],
)
def test_invalid_stage_dispersion_scale_is_rejected(reference, scales, match) -> None:
    with pytest.raises(SemiSyntheticTrajectoryError, match=rf"stage_dispersion_scale must .*{match}"):
        generate_semisynthetic_trajectory(_params(stage_dispersion_scale=scales), reference=reference)


@pytest.mark.parametrize(
    "null", [{"trajectory_mode": "none"}, {"trajectory_mode": "orientation", "group_effect_size": 0.0}]
)
def test_a_zero_effect_cell_keeps_the_scale(reference, null) -> None:
    dataset = generate_semisynthetic_trajectory(_params(**null, stage_dispersion_scale=STAGES), reference=reference)
    assert dataset.truth["stage_dispersion"]["scale"] == list(STAGES)
    plain = generate_semisynthetic_trajectory(_params(**null), reference=reference)
    _assert_stretched(dataset, plain, _row_multipliers(plain, STAGES))


# --------------------------------------------------------------------------- #
# 1.4 Pooled scales with one multiplier per cell
# --------------------------------------------------------------------------- #


def _span(**overrides: Any) -> SemiSyntheticTrajectoryParams:
    base: dict[str, Any] = {
        "seed": 23,
        "n_samples": 600,
        "p_dmp": 0.1,
        "group_stage_sizes": SEAAD_TABLE,
        "construction": "stage_span",
        "trajectory_mode": "translation",
        "group_effect_size": 1.0,
        "span_scale": "planted",
        "stage_dispersion_scale": STAGES,
    }
    return SemiSyntheticTrajectoryParams(**{**base, **overrides})


def _group_means(dataset, layers=LAYERS):
    population = dataset.population_trajectories.layers
    means_a = {layer: population[layer].loc["A"].to_numpy() for layer in layers}
    means_b = {layer: population[layer].loc["B"].to_numpy() for layer in layers}
    sizes = dataset.metadata.groupby(["group", "stage"]).size()
    return means_a, means_b, sizes.loc["A"].to_numpy(dtype=float), sizes.loc["B"].to_numpy(dtype=float)


def test_the_cell_multiplier_is_the_size_weighted_mean_square() -> None:
    value = cell_dispersion_within_multiplier((1.0, 2.0), (0.5, 1.5), (3, 1), (2, 4))
    assert value == pytest.approx((3 * 1 + 1 * 4 + 2 * 0.25 + 4 * 2.25) / 10)
    with pytest.raises(SemiSyntheticTrajectoryError, match="one dispersion multiplier per cell"):
        cell_dispersion_within_multiplier((1.0, 2.0), (0.5,), (3, 1), (2, 4))


@pytest.mark.parametrize(("sizes_a", "sizes_b"), [((10, 10, 10), (5, 5, 5)), ((7, 8, 9), (3, 2, 1))])
@pytest.mark.parametrize("group", [1.5, 2.0])
def test_a_group_scale_alone_is_the_0_26_multiplier(sizes_a, sizes_b, group) -> None:
    expected = group_dispersion_within_multiplier(group, sizes_a, sizes_b)
    # The general formula agrees with the group formula...
    general = cell_dispersion_within_multiplier((1.0,) * 3, (group,) * 3, sizes_a, sizes_b)
    assert general == pytest.approx(expected, rel=1e-14)
    # ...and without a stage scale the dispatch is the 0.26.0 computation, bit for bit.
    params = _params(group_dispersion_scale=group)
    assert dispersion_within_multiplier(params, sizes_a, sizes_b) == expected
    assert dispersion_within_multiplier(_params(), sizes_a, sizes_b) is None


def test_equal_stage_multipliers_give_their_square() -> None:
    params = _params(stage_dispersion_scale=(1.5, 1.5, 1.5))
    assert dispersion_within_multiplier(params, (11, 10, 28), (9, 10, 12)) == pytest.approx(2.25, rel=1e-14)


def test_the_multiplier_combines_stage_and_group_per_cell() -> None:
    params = _params(stage_dispersion_scale=(1.0, 2.0, 0.5), group_dispersion_scale=3.0)
    sizes_a, sizes_b = (11, 10, 28), (9, 10, 12)
    expected = (11 * 1 + 10 * 4 + 28 * 0.25 + 9 * 9 + 10 * 36 + 12 * 2.25) / 80
    assert dispersion_within_multiplier(params, sizes_a, sizes_b) == pytest.approx(expected, rel=1e-14)


def test_pooled_scales_carry_the_stage_weighted_factor(reference) -> None:
    dataset = generate_semisynthetic_trajectory(
        _params(n_samples=600, group_stage_sizes=SEAAD_TABLE, stage_dispersion_scale=STAGES), reference=reference
    )
    means_a, _, sizes_a, sizes_b = _group_means(dataset)
    weights = sizes_a + sizes_b
    multiplier = float(weights @ np.asarray(STAGES) ** 2 / weights.sum())
    computed = dispersion_within_multiplier(_params(stage_dispersion_scale=STAGES), sizes_a, sizes_b)
    assert computed == pytest.approx(multiplier, rel=1e-14)
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


def test_planted_scale_matches_a_large_sample_with_unequal_stage_sizes(reference) -> None:
    table = ((3000, 1500, 600), (1200, 2400, 900))
    params = _span(group_stage_sizes=table, stage_dispersion_scale=(1.0, 1.6, 0.6))
    dataset = generate_semisynthetic_trajectory(params, reference=reference)
    means_a, means_b, sizes_a, sizes_b = _group_means(dataset)
    np.testing.assert_array_equal(sizes_a, table[0])
    multiplier = dispersion_within_multiplier(params, sizes_a, sizes_b)
    expected = planted_pooled_scale(
        means_a, means_b, sizes_a, sizes_b, layers=LAYERS, reference=reference, within_multiplier=multiplier
    )
    observed = np.concatenate([_analysed(dataset, layer).std(axis=0) for layer in LAYERS])
    relative = np.abs(observed / expected - 1.0)
    # The SD of a sample SD is about SD / sqrt(2n): 0.5 % at n = 9,600.
    assert np.median(relative) < 0.01
    assert relative.max() < 0.04
    # The reference-variance scale is not the scale these data have.
    unaware = planted_pooled_scale(means_a, means_b, sizes_a, sizes_b, layers=LAYERS, reference=reference)
    assert np.median(np.abs(observed / unaware - 1.0)) > 0.05


_MODE_EFFECTS = (("magnitude", 0.5), ("orientation", 0.5), ("orientation", 1.0), ("shape", 0.5), ("translation", 1.0))


@pytest.mark.parametrize("group", [None, 1.5])
@pytest.mark.parametrize(("mode", "effect"), _MODE_EFFECTS)
def test_planted_edit_is_exact_with_a_stage_scale(reference, mode, effect, group) -> None:
    params = _span(seed=11, trajectory_mode=mode, group_effect_size=effect, group_dispersion_scale=group)
    dataset = generate_semisynthetic_trajectory(params, reference=reference)
    block = dataset.truth["stage_span"]
    assert block["units"] == "planted_standardized"
    assert block["scale_residual"] <= semisynthetic.PLANTED_SCALE_TOLERANCE
    means_a, means_b, sizes_a, sizes_b = _group_means(dataset)
    multiplier = dispersion_within_multiplier(params, sizes_a, sizes_b)
    scale = planted_pooled_scale(
        means_a, means_b, sizes_a, sizes_b, layers=LAYERS, reference=reference, within_multiplier=multiplier
    )
    a = np.hstack([means_a[layer] for layer in LAYERS]) / scale
    b = np.hstack([means_b[layer] for layer in LAYERS]) / scale
    expected, _ = edit_stage_configuration(a, mode, effect)
    path = float(np.linalg.norm(np.diff(a, axis=0), axis=1).sum())
    assert np.max(np.abs(b - expected)) <= 1e-8 * path
    assert block["path_length"] == pytest.approx(path, rel=1e-10)
    if mode == "orientation":
        assert block["intended_geometry"]["angle"] == pytest.approx(90.0 * effect, abs=1e-5)
    # The scale without the stage multipliers states the edit in a scale the data do not have.
    unaware = planted_pooled_scale(means_a, means_b, sizes_a, sizes_b, layers=LAYERS, reference=reference)
    if mode != "magnitude":
        a_u = np.hstack([means_a[layer] for layer in LAYERS]) / unaware
        b_u = np.hstack([means_b[layer] for layer in LAYERS]) / unaware
        expected_u, _ = edit_stage_configuration(a_u, mode, effect)
        assert np.max(np.abs(b_u - expected_u)) > 1e-6 * path


# --------------------------------------------------------------------------- #
# 1.5 Identifiers and study configurations
# --------------------------------------------------------------------------- #


def _cell(generator: SemiSyntheticTrajectoryParams):
    return make_simulation_cell(
        phase="power_primary",
        generator_params=generator,
        evaluation_params=SimulationEvaluationParams(integration_method="concat", permutations=0),
    )


def test_a_set_scale_changes_the_cell_identity() -> None:
    plain, scaled = _cell(_params()), _cell(_params(stage_dispersion_scale=STAGES))
    assert plain.cell_id != scaled.cell_id
    assert parameter_signature(plain) != parameter_signature(scaled)
    assert "stage_dispersion_scale" not in _to_jsonable(plain.generator_params)
    assert _to_jsonable(scaled.generator_params)["stage_dispersion_scale"] == list(STAGES)
    unit = _cell(_params(stage_dispersion_scale=(1, 1, 1)))
    assert (unit.cell_id, parameter_signature(unit)) == (plain.cell_id, parameter_signature(plain))
    spelled = _cell(_params(stage_dispersion_scale=list(STAGES)))
    assert (spelled.cell_id, parameter_signature(spelled)) == (scaled.cell_id, parameter_signature(scaled))


def test_type_i_grid_keeps_the_scale() -> None:
    grid = enumerate_type_i_grid(
        baseline_generator_params=_params(
            trajectory_mode="orientation", group_effect_size=0.5, stage_dispersion_scale=STAGES
        ),
        evaluation_params=SimulationEvaluationParams(integration_method="concat", permutations=0),
        axes={"generator.group_ratio": (0.5, 0.8)},
    )
    assert len(grid.cells) == 2
    for cell in grid.cells:
        assert cell.generator_params.group_effect_size == 0.0
        assert cell.generator_params.stage_dispersion_scale == STAGES


def _config(**overrides: Any) -> StudyConfig:
    defaults: dict[str, Any] = {
        "generator": SemiSyntheticTrajectoryParams(
            seed=2, trajectory_mode="magnitude", n_samples=60, p_dmp=0.1, stage_dispersion_scale=STAGES
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
        assert cell.generator_params.stage_dispersion_scale == STAGES
    without = _config(generator=replace(config.generator, stage_dispersion_scale=None))
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


@pytest.mark.parametrize("spelling", [None, [1, 1, 1], [1.0, 1.0, 1.0]])
def test_null_or_unit_scale_equals_an_omitted_key(tmp_path: Path, spelling) -> None:
    omitted = load_study_config(_write_config(tmp_path, {}))
    spelled = load_study_config(_write_config(tmp_path, {"stage_dispersion_scale": spelling}))
    assert spelled == omitted
    assert _identity(spelled) == _identity(omitted)
    dumped = tmp_path / "dumped.json"
    dump_study_config(spelled, dumped)
    assert "stage_dispersion_scale" not in json.loads(dumped.read_text(encoding="utf-8"))["generator"]


def test_a_set_scale_loads_enumerates_and_dumps(tmp_path: Path) -> None:
    config = load_study_config(_write_config(tmp_path, {"stage_dispersion_scale": list(STAGES)}))
    assert config.generator.stage_dispersion_scale == STAGES
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert json.loads(dumped.read_text(encoding="utf-8"))["generator"]["stage_dispersion_scale"] == list(STAGES)
    assert load_study_config(dumped) == config
    assert _identity(load_study_config(dumped)) == _identity(config)


def test_a_wrong_length_fails_when_the_config_loads(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match=r"stage_dispersion_scale must have one multiplier per stage"):
        load_study_config(_write_config(tmp_path, {"stage_dispersion_scale": [1.5, 0.5]}))


def test_scale_as_a_design_axis(tmp_path: Path) -> None:
    path = _write_config(tmp_path, {}, design_axes={SCALE_AXIS: [None, list(STAGES), [0.8, 0.9, 1.2]]})
    config = load_study_config(path)
    assert len(config.design_points()) == 3
    design = [cell for cell in enumerate_study(config).cells if cell.phase == DESIGN_PHASE]
    seen = set()
    for cell in design:
        point = cell.metadata[DESIGN_POINT_KEY]
        value = point[SCALE_AXIS]
        expected = None if value is None else tuple(float(v) for v in value)
        assert cell.generator_params.stage_dispersion_scale == expected
        seen.add(expected)
    assert seen == {STAGES, (0.8, 0.9, 1.2)}
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert load_study_config(dumped) == config


def test_a_wrong_length_on_an_axis_fails_at_enumeration(tmp_path: Path) -> None:
    path = _write_config(tmp_path, {}, design_axes={SCALE_AXIS: [None, [0.8, 1.2]]})
    with pytest.raises(StudyConfigError, match=r"generator\.stage_dispersion_scale must have one multiplier per stage"):
        enumerate_study(load_study_config(path))


def test_crossing_the_stage_count_fails_at_enumeration(tmp_path: Path) -> None:
    path = _write_config(
        tmp_path, {"stage_dispersion_scale": list(STAGES)}, design_axes={"generator.n_stages": [3, 4]}
    )
    with pytest.raises(StudyConfigError, match=r"generator\.stage_dispersion_scale .*n_stages=4"):
        enumerate_study(load_study_config(path))
