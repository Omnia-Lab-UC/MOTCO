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

from motco.simulations import SemiSyntheticTrajectoryParams, SimulationEvaluationParams
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
RHO = "generator.baseline_continuity"
N = "generator.n_samples"

BASELINE = {METHOD: "pls", FORCED: None, NCOMP: None, EST: None}
FIVE_POINTS = [
    BASELINE,
    {METHOD: "pls", FORCED: None, NCOMP: None, EST: "net_displacement"},
    {METHOD: "pls", FORCED: 4, NCOMP: None, EST: "net_displacement"},
    {METHOD: "pca", FORCED: None, NCOMP: 4, EST: "net_displacement"},
    {METHOD: "concat", FORCED: None, NCOMP: None, EST: "net_displacement"},
]


def _raw(design_grid: dict, **overrides) -> dict:
    raw: dict = {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "n_samples": 60, "p_dmp": 0.1},
        "evaluation": {"integration_method": "pls", "permutations": 0, "seed": 3},
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
    config = _load(tmp_path, _raw({"points": FIVE_POINTS}))
    assert config.design_grid.explicit
    assert config.design_points() == FIVE_POINTS
    assert config.baseline_design_point() == BASELINE
    assert list(config.design_grid.axes) == [METHOD, FORCED, NCOMP, EST]
    assert config.design_grid.axes[METHOD] == ("pls", "pca", "concat")
    assert config.design_grid.axes[FORCED] == (None, 4)
    assert config.design_grid.axes[NCOMP] == (None, 4)
    assert config.design_grid.axes[EST] == (None, "net_displacement")


def test_explicit_points_round_trip_through_dump(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"points": FIVE_POINTS}))
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    payload = json.loads(dumped.read_text(encoding="utf-8"))
    assert payload["design_grid"] == {"points": FIVE_POINTS}
    reloaded = load_study_config(dumped)
    assert reloaded.design_points() == FIVE_POINTS
    assert [cell.cell_id for cell in enumerate_study(reloaded).cells] == [
        cell.cell_id for cell in enumerate_study(config).cells
    ]


def test_crossed_form_still_dumps_axes(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"axes": {EST: [None, "net_displacement"]}}))
    assert not config.design_grid.explicit
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert json.loads(dumped.read_text(encoding="utf-8"))["design_grid"] == {
        "axes": {EST: [None, "net_displacement"]}
    }


# ── Rejections ────────────────────────────────────────────────────────────────


def test_both_forms_are_rejected(tmp_path: Path) -> None:
    raw = _raw({"axes": {EST: [None, "net_displacement"]}, "points": FIVE_POINTS})
    with pytest.raises(StudyConfigError, match="both 'axes' and 'points'"):
        _load(tmp_path, raw)


def test_point_with_a_differing_axis_set_is_rejected(tmp_path: Path) -> None:
    short = {METHOD: "concat", FORCED: None, EST: "net_displacement"}
    with pytest.raises(StudyConfigError, match=r"points\[2\].*missing axis.*n_components"):
        _load(tmp_path, _raw({"points": [BASELINE, FIVE_POINTS[1], short]}))
    extra = {**FIVE_POINTS[1], RHO: 0.5}
    with pytest.raises(StudyConfigError, match=r"points\[1\].*unexpected axis.*baseline_continuity"):
        _load(tmp_path, _raw({"points": [BASELINE, extra]}))


def test_duplicate_points_are_rejected(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match=r"points\[2\].*duplicates design_grid.points\[1\]"):
        _load(tmp_path, _raw({"points": [BASELINE, FIVE_POINTS[1], FIVE_POINTS[1]]}))


def test_pc1_and_null_estimator_are_the_same_point(tmp_path: Path) -> None:
    pc1 = {**BASELINE, EST: "pc1"}
    with pytest.raises(StudyConfigError, match=r"points\[1\].*duplicates design_grid.points\[0\]"):
        _load(tmp_path, _raw({"points": [BASELINE, pc1]}))
    # Alone, "pc1" normalizes to null and is the baseline point.
    config = _load(tmp_path, _raw({"points": [pc1, FIVE_POINTS[1]]}))
    assert config.design_points()[0][EST] is None


def test_missing_baseline_point_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match="must include the baseline point"):
        _load(tmp_path, _raw({"points": FIVE_POINTS[1:]}))


def test_explicit_axis_also_declared_ofat_is_rejected(tmp_path: Path) -> None:
    raw = _raw({"points": FIVE_POINTS}, axes={EST: [None, "net_displacement"]})
    with pytest.raises(StudyConfigError, match="also declared under axes"):
        _load(tmp_path, raw)


def test_unknown_namespace_in_a_point_is_rejected(tmp_path: Path) -> None:
    with pytest.raises(StudyConfigError, match="unsupported namespace"):
        _load(tmp_path, _raw({"points": [{"study.foo": 1}]}))


# ── Enumeration ───────────────────────────────────────────────────────────────


def test_explicit_points_enumerate_their_own_anchored_grids(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw({"points": FIVE_POINTS}))
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
