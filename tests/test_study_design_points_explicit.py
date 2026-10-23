"""Explicit design points: loading, validation, enumeration, crossed-form agreement.

The explicit form lists design points directly instead of crossing axis values,
so non-rectangular measurement columns (a forced-rank PLS column, a PCA column
with its own rank, a concat column) can share one study. These tests pin that
contract on small configs; nothing is evaluated.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from motco.simulations import SemiSyntheticTrajectoryParams, SimulationEvaluationError, SimulationEvaluationParams
from motco.simulations.grid import derive_replicate_seed, parameter_signature
from motco.simulations.study import (
    DesignGrid,
    MatchedSeedPolicy,
    StudyConfig,
    StudyConfigError,
    dump_study_config,
    enumerate_study,
    load_study_config,
)
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

METHOD = "evaluation.integration_method"
FORCED = "evaluation.integration_params.forced_components"
NCOMP = "evaluation.integration_params.n_components"
EST = "evaluation.orientation_estimator"
FLOOR = "evaluation.min_components"
RHO = "generator.baseline_continuity"
N = "generator.n_samples"

# Loaded (normalized) points: the PC1 estimator is held as ``None``.
BASELINE = {METHOD: "pls", FORCED: None, NCOMP: None, EST: None}
FIVE_POINTS = [
    BASELINE,
    {METHOD: "pls", FORCED: None, NCOMP: None, EST: "net_displacement"},
    {METHOD: "pls", FORCED: 4, NCOMP: None, EST: "net_displacement"},
    {METHOD: "pca", FORCED: None, NCOMP: 4, EST: "net_displacement"},
    {METHOD: "concat", FORCED: None, NCOMP: None, EST: "net_displacement"},
]
# Config spelling: PC1 is ``"pc1"``; a ``null`` estimator names the evaluator's default.
RAW_BASELINE = {**BASELINE, EST: "pc1"}
RAW_FIVE_POINTS = [RAW_BASELINE, *FIVE_POINTS[1:]]


def _raw(design_grid: dict, **overrides) -> dict:
    raw: dict = {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "n_samples": 60, "p_dmp": 0.1},
        "evaluation": {"integration_method": "pls", "permutations": 0, "seed": 3, "orientation_estimator": "pc1"},
        "trajectory_modes": ["translation"],
        "effect_sizes": [0.0, 1.0],
        "n_replicates": 2,
        "base_seed": 100,
        "matched_seeds": {"enabled": True, "primary_family": "fam"},
        "design_grid": design_grid,
    }
    raw.update(overrides)
    return raw


def _load(tmp_path: Path, raw: dict) -> StudyConfig:
    path = tmp_path / "config.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    return load_study_config(path)


# ── Loading and round trip ────────────────────────────────────────────────────


def test_explicit_points_load_in_listed_order(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"points": RAW_FIVE_POINTS}))
    assert config.design_grid.explicit
    assert config.design_points() == FIVE_POINTS
    assert config.baseline_design_point() == BASELINE
    assert list(config.design_grid.axes) == [METHOD, FORCED, NCOMP, EST]
    assert config.design_grid.axes[METHOD] == ("pls", "pca", "concat")
    assert config.design_grid.axes[FORCED] == (None, 4)
    assert config.design_grid.axes[NCOMP] == (None, 4)
    assert config.design_grid.axes[EST] == (None, "net_displacement")


def test_explicit_points_round_trip_through_dump(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"points": RAW_FIVE_POINTS}))
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    payload = json.loads(dumped.read_text(encoding="utf-8"))
    assert payload["design_grid"] == {"points": RAW_FIVE_POINTS}
    assert payload["evaluation"]["orientation_estimator"] == "pc1"
    reloaded = load_study_config(dumped)
    assert reloaded.design_points() == FIVE_POINTS
    assert [cell.cell_id for cell in enumerate_study(reloaded).cells] == [
        cell.cell_id for cell in enumerate_study(config).cells
    ]


def test_crossed_form_still_dumps_axes(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"axes": {EST: ["pc1", "net_displacement"]}}))
    assert not config.design_grid.explicit
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert json.loads(dumped.read_text(encoding="utf-8"))["design_grid"] == {
        "axes": {EST: ["pc1", "net_displacement"]}
    }


# ── Rejections ────────────────────────────────────────────────────────────────


def test_both_forms_are_rejected(tmp_path: Path) -> None:
    raw = _raw({"axes": {EST: ["pc1", "net_displacement"]}, "points": RAW_FIVE_POINTS})
    with pytest.raises(StudyConfigError, match="both 'axes' and 'points'"):
        _load(tmp_path, raw)


def test_point_with_a_differing_axis_set_is_rejected(tmp_path: Path) -> None:
    short = {METHOD: "concat", FORCED: None, EST: "net_displacement"}
    with pytest.raises(StudyConfigError, match=r"points\[2\].*missing axis.*n_components"):
        _load(tmp_path, _raw({"points": [RAW_BASELINE, FIVE_POINTS[1], short]}))
    extra = {**FIVE_POINTS[1], RHO: 0.5}
    with pytest.raises(StudyConfigError, match=r"points\[1\].*unexpected axis.*baseline_continuity"):
        _load(tmp_path, _raw({"points": [RAW_BASELINE, extra]}))


def test_duplicate_points_are_rejected(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match=r"points\[2\].*duplicates design_grid.points\[1\]"):
        _load(tmp_path, _raw({"points": [RAW_BASELINE, FIVE_POINTS[1], FIVE_POINTS[1]]}))


def test_null_and_net_displacement_estimator_are_the_same_point(tmp_path: Path) -> None:
    default = {**BASELINE, EST: None}
    with pytest.raises(StudyConfigError, match=r"points\[2\].*duplicates design_grid.points\[1\].*null names"):
        _load(tmp_path, _raw({"points": [RAW_BASELINE, FIVE_POINTS[1], default]}))
    # Alone, null normalizes to the evaluator's default, net displacement.
    config = _load(tmp_path, _raw({"points": [RAW_BASELINE, default]}))
    assert config.design_points()[1][EST] == "net_displacement"


def test_null_and_pc1_estimator_are_distinct_points(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"points": [RAW_BASELINE, {**BASELINE, EST: None}]}))
    assert [point[EST] for point in config.design_points()] == [None, "net_displacement"]
    assert config.baseline_design_point()[EST] is None


def test_null_point_is_the_baseline_when_the_baseline_evaluation_omits_the_estimator(tmp_path: Path) -> None:
    raw = _raw({"points": [{EST: None}, {EST: "pc1"}]})
    del raw["evaluation"]["orientation_estimator"]
    config = _load(tmp_path, raw)
    assert config.evaluation.orientation_estimator == "net_displacement"
    assert config.baseline_design_point() == {EST: "net_displacement"}
    cells = enumerate_study(config).cells
    primary = [cell for cell in cells if cell.phase != DESIGN_PHASE]
    design = [cell for cell in cells if cell.phase == DESIGN_PHASE]
    assert primary and design
    assert all(cell.evaluation_params.orientation_estimator == "net_displacement" for cell in primary)
    assert {json.dumps(cell.metadata[DESIGN_POINT_KEY]) for cell in design} == {json.dumps({EST: None})}
    assert all(cell.evaluation_params.orientation_estimator is None for cell in design)


def test_pc1_baseline_with_a_net_displacement_point_enumerates_as_the_committed_profiles(
    tmp_path: Path,
) -> None:
    config = _load(tmp_path, _raw({"points": [{EST: "pc1"}, {EST: "net_displacement"}]}))
    assert config.evaluation.orientation_estimator is None
    assert config.baseline_design_point() == {EST: None}
    cells = enumerate_study(config).cells
    primary = [cell for cell in cells if cell.phase != DESIGN_PHASE]
    design = [cell for cell in cells if cell.phase == DESIGN_PHASE]
    assert primary and design
    assert all(cell.evaluation_params.orientation_estimator is None for cell in primary)
    assert all(cell.evaluation_params.orientation_estimator == "net_displacement" for cell in design)
    # The PC1 primary column evaluates as a config without the axis; a config
    # that omits the estimator evaluates (and hashes) with the default instead.
    plain = _raw({})
    del plain["design_grid"]
    plain_cells = enumerate_study(_load(tmp_path, plain)).cells
    assert [c.evaluation_params for c in primary] == [c.evaluation_params for c in plain_cells]
    del plain["evaluation"]["orientation_estimator"]
    unnamed_cells = enumerate_study(_load(tmp_path, plain)).cells
    assert all(c.evaluation_params.orientation_estimator == "net_displacement" for c in unnamed_cells)
    assert {c.cell_id for c in unnamed_cells}.isdisjoint({c.cell_id for c in plain_cells})


# ── Rank floor: field and axis ────────────────────────────────────────────────


def test_floor_spelling_in_a_config(tmp_path: Path) -> None:
    plain = _raw({})
    del plain["design_grid"]
    assert "min_components" not in plain["evaluation"]
    assert _load(tmp_path, plain).evaluation.min_components == 2  # omitted: the evaluator's default
    plain["evaluation"]["min_components"] = 1
    assert _load(tmp_path, plain).evaluation.min_components is None  # the unfloored rule, hash-omitted
    plain["evaluation"]["min_components"] = None
    assert _load(tmp_path, plain).evaluation.min_components == 2  # null names the default
    plain["evaluation"]["min_components"] = 3
    assert _load(tmp_path, plain).evaluation.min_components == 3
    plain["evaluation"]["min_components"] = 0
    # Anything but null and 1 passes through for the evaluator to validate.
    with pytest.raises(SimulationEvaluationError, match="min_components must be an integer of at least 1"):
        _load(tmp_path, plain)


def test_floor_one_config_enumerates_as_a_config_from_before_the_floor(tmp_path: Path) -> None:
    plain = _raw({})
    del plain["design_grid"]
    before = enumerate_study(_load(tmp_path, plain)).cells
    plain["evaluation"]["min_components"] = 1
    pinned = enumerate_study(_load(tmp_path, plain)).cells
    # ``before`` is the default (floor 2): the pinned config must hash differently from it ...
    assert {c.cell_id for c in pinned}.isdisjoint({c.cell_id for c in before})
    # ... and identically to the same config evaluated with the field held as None.
    assert all(c.evaluation_params.min_components is None for c in pinned)
    assert all("min_components" not in parameter_signature(c) for c in pinned)


def test_null_and_one_floor_are_distinct_points(tmp_path: Path) -> None:
    raw = _raw({"points": [{FLOOR: 1}, {FLOOR: None}]}, evaluation={**_raw({})["evaluation"], "min_components": 1})
    config = _load(tmp_path, raw)
    assert config.evaluation.min_components is None
    assert [point[FLOOR] for point in config.design_points()] == [None, 2]
    assert config.baseline_design_point()[FLOOR] is None
    with pytest.raises(StudyConfigError, match=r"duplicates design_grid.points\[1\].*min_components.*null names"):
        _load(tmp_path, _raw({"points": [{FLOOR: 1}, {FLOOR: 2}, {FLOOR: None}]},
                             evaluation={**_raw({})["evaluation"], "min_components": 1}))


def test_null_floor_point_is_the_baseline_when_the_baseline_omits_the_floor(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"points": [{FLOOR: None}, {FLOOR: 1}]}))
    assert config.evaluation.min_components == 2
    assert config.baseline_design_point() == {FLOOR: 2}
    cells = enumerate_study(config).cells
    primary = [cell for cell in cells if cell.phase != DESIGN_PHASE]
    design = [cell for cell in cells if cell.phase == DESIGN_PHASE]
    assert primary and design
    assert all(cell.evaluation_params.min_components == 2 for cell in primary)
    assert {json.dumps(cell.metadata[DESIGN_POINT_KEY]) for cell in design} == {json.dumps({FLOOR: None})}
    assert all(cell.evaluation_params.min_components is None for cell in design)


def test_floor_axis_round_trips_through_dump(tmp_path: Path) -> None:
    raw = _raw({"points": [{FLOOR: 1}, {FLOOR: 2}]}, evaluation={**_raw({})["evaluation"], "min_components": 1})
    config = _load(tmp_path, raw)
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    payload = json.loads(dumped.read_text(encoding="utf-8"))
    assert payload["evaluation"]["min_components"] == 1  # never written as null, which would mean the default
    assert payload["design_grid"] == {"points": [{FLOOR: 1}, {FLOOR: 2}]}
    reloaded = load_study_config(dumped)
    assert reloaded.design_points() == [{FLOOR: None}, {FLOOR: 2}]
    assert [cell.cell_id for cell in enumerate_study(reloaded).cells] == [
        cell.cell_id for cell in enumerate_study(config).cells
    ]


def test_missing_baseline_point_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match="must include the baseline point"):
        _load(tmp_path, _raw({"points": FIVE_POINTS[1:]}))


def test_explicit_axis_also_declared_ofat_is_rejected(tmp_path: Path) -> None:
    raw = _raw({"points": RAW_FIVE_POINTS}, axes={EST: ["pc1", "net_displacement"]})
    with pytest.raises(StudyConfigError, match="also declared under axes"):
        _load(tmp_path, raw)


def test_unknown_namespace_in_a_point_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match="unsupported namespace"):
        _load(tmp_path, _raw({"points": [{"study.foo": 1}]}))


# ── Enumeration ───────────────────────────────────────────────────────────────


def test_explicit_points_enumerate_their_own_anchored_grids(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"points": RAW_FIVE_POINTS}))
    design = [cell for cell in enumerate_study(config).cells if cell.phase == DESIGN_PHASE]
    by_point: dict[str, list] = {}
    for cell in design:
        by_point.setdefault(json.dumps(cell.metadata[DESIGN_POINT_KEY], sort_keys=True), []).append(cell)
    assert len(by_point) == 4
    for point in FIVE_POINTS[1:]:
        cells = by_point[json.dumps(point, sort_keys=True)]
        assert sorted((c.metadata["trajectory_mode"], c.metadata["effect_size"]) for c in cells) == [
            ("none", 0.0),
            ("translation", 1.0),
        ]
        for cell in cells:
            evaluation = cell.evaluation_params
            assert evaluation.integration_method == point[METHOD]
            assert evaluation.orientation_estimator == point[EST]
            params = dict(evaluation.integration_params)
            assert params.get("forced_components") == point[FORCED]
            assert params.get("n_components") == point[NCOMP]
    assert not any(cell.metadata.get(DESIGN_POINT_KEY) == BASELINE for cell in design)


def _signature_view(config: StudyConfig) -> list[tuple]:
    return [
        (
            cell.cell_id,
            parameter_signature(cell),
            tuple(derive_replicate_seed(cell, i) for i in range(cell.n_replicates)),
        )
        for cell in enumerate_study(config).cells
    ]


def test_explicit_list_of_a_rectangle_equals_the_crossed_form() -> None:
    generator = SemiSyntheticTrajectoryParams(seed=2, trajectory_mode="magnitude", n_samples=60, p_dmp=0.1)
    common: dict = {
        "generator": generator,
        "evaluation": SimulationEvaluationParams(integration_method="concat", permutations=0, seed=3),
        "trajectory_modes": ("magnitude", "orientation"),
        "effect_sizes": (0.0, 0.5, 1.0),
        "n_replicates": 2,
        "base_seed": 100,
        "matched_seeds": MatchedSeedPolicy(enabled=True, primary_family="fam"),
    }
    crossed = StudyConfig(**common, design_grid=DesignGrid(axes={RHO: (0.0, 0.5), N: (60, 120)}))
    explicit = StudyConfig(
        **common,
        design_grid=DesignGrid(
            points=tuple({RHO: rho, N: n} for rho in (0.0, 0.5) for n in (60, 120))
        ),
    )
    assert explicit.design_grid.axes == crossed.design_grid.axes
    assert _signature_view(explicit) == _signature_view(crossed)
