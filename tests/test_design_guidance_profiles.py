"""Structure of the four committed design-guidance pilot profiles.

Checks the spec scenarios of ``design-guidance-grid``: the core profiles cross
sample size x stages x blocks into sixteen anchored columns whose block columns
share data; the arm profiles move one generator setting at a time off the core
baseline (the three per-omic stage effects counting as one, set jointly); every
profile names the production measurement explicitly; and every cell fits the
surgery headroom. Identity (cell ids, signatures, seeds) is pinned separately
by ``tests/test_post_floor_config_identity.py``.
"""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

from motco.simulations.grid import SimulationCell, derive_replicate_seed
from motco.simulations.semisynthetic import POOL_LIMITED_MODES, expected_surgery_headroom
from motco.simulations.study import StudyConfig, enumerate_study, load_study_config
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"

CORE = ("design_guidance_core_magnitude_pilot.json", "design_guidance_core_geometry_pilot.json")
ARMS = ("design_guidance_arms_magnitude_pilot.json", "design_guidance_arms_geometry_pilot.json")
MAGNITUDE = ("design_guidance_core_magnitude_pilot.json", "design_guidance_arms_magnitude_pilot.json")
GEOMETRY = ("design_guidance_core_geometry_pilot.json", "design_guidance_arms_geometry_pilot.json")
PROFILES = CORE + ARMS

N = "generator.n_samples"
STAGES = "generator.n_stages"
LAYERS = "evaluation.integration_params.layers"
DELTAS = ("generator.delta_methyl", "generator.delta_expr", "generator.delta_protein")

#: Every arm column besides the baseline, as the setting it moves and the value.
EXPECTED_ARMS: tuple[tuple[str, Any], ...] = (
    ("generator.group_ratio", 0.65),
    ("generator.group_ratio", 0.8),
    ("generator.stage_sample_prop", (0.5, 0.3, 0.2)),
    ("generator.p_dmp", 0.05),
    ("generator.p_dmp", 0.15),
    ("stage_effect", 0.25),
    ("stage_effect", 0.5),
    ("stage_effect", 1.0),
    ("generator.baseline_continuity", 0.5),
    ("generator.baseline_continuity", 0.8),
)


def _load(name: str) -> StudyConfig:
    return load_study_config(CONFIG_DIR / name)


def _raw(name: str) -> dict[str, Any]:
    return json.loads((CONFIG_DIR / name).read_text(encoding="utf-8"))


def _frozen(value: Any) -> Any:
    return tuple(_frozen(item) for item in value) if isinstance(value, list | tuple) else value


def _columns(config: StudyConfig) -> dict[tuple[tuple[str, Any], ...], list[SimulationCell]]:
    """Design cells grouped by design point, the baseline column included."""

    columns: dict[tuple[tuple[str, Any], ...], list[SimulationCell]] = {}
    for cell in enumerate_study(config).cells:
        if cell.phase not in ("power_primary", DESIGN_PHASE):
            continue
        point = cell.metadata[DESIGN_POINT_KEY]
        key = tuple((axis, _frozen(value)) for axis, value in point.items())
        columns.setdefault(key, []).append(cell)
    return columns


def _mode_effect(cell: SimulationCell) -> tuple[str, float]:
    return str(cell.metadata["trajectory_mode"]), float(cell.metadata["effect_size"])


# ── Profiles name the production measurement ─────────────────────────────────


@pytest.mark.parametrize("name", PROFILES)
def test_profile_names_the_production_measurement_explicitly(name: str) -> None:
    raw = _raw(name)
    assert raw["evaluation"]["integration_method"] == "pls"
    assert raw["evaluation"]["orientation_estimator"] == "net_displacement"
    assert raw["evaluation"]["min_components"] == 2
    assert raw["generator"]["magnitude_kind"] == "joint"
    assert raw["attribution"] == {"enabled": False}
    assert raw["acceptance"]["anchor_decision"]["statistics"] == ["delta", "angle", "shape"]
    assert raw["acceptance"]["gate"] == {"enabled": False}
    assert raw["acceptance"]["type_i"]

    config = _load(name)
    assert config.evaluation.orientation_estimator == "net_displacement"
    assert config.evaluation.min_components == 2
    assert config.generator.magnitude_kind == "joint"
    assert config.n_replicates == 100 and config.evaluation.permutations == 199
    assert config.base_seed == 1028 and config.matched_seeds.primary_family == "design-guidance"
    assert config.acceptance.anchor_decision is not None and config.acceptance.anchor_decision.alpha == 0.05


@pytest.mark.parametrize("name", PROFILES)
def test_profile_baseline_is_the_shared_design_guidance_baseline(name: str) -> None:
    generator = _load(name).generator
    assert (generator.n_samples, generator.n_stages, generator.group_ratio, generator.p_dmp) == (200, 3, 0.5, 0.1)
    assert generator.baseline_continuity == 0.0 and generator.stage_sample_prop is None
    assert (generator.delta_methyl, generator.delta_expr, generator.delta_protein) == (2.0, 2.0, 2.0)
    assert generator.surgery_censoring == "error"
    assert _load(name).evaluation.integration_params["layers"] == ["methylation", "expression"]


# ── Mode and effect lists ────────────────────────────────────────────────────


@pytest.mark.parametrize("name", MAGNITUDE)
def test_magnitude_profiles_run_magnitude_on_the_small_effect_axis(name: str) -> None:
    config = _load(name)
    assert config.trajectory_modes == ("magnitude",)
    assert config.effect_sizes == (0.0, 0.01, 0.02, 0.05, 0.1, 0.25)


@pytest.mark.parametrize("name", GEOMETRY)
def test_geometry_profiles_run_the_geometry_modes_on_the_large_effect_axis(name: str) -> None:
    config = _load(name)
    assert config.trajectory_modes == ("orientation", "shape", "translation")
    assert config.effect_sizes == (0.0, 0.25, 0.5, 1.0)


# ── Column structure ─────────────────────────────────────────────────────────


@pytest.mark.parametrize(("name", "n_columns"), [(n, 16) for n in CORE] + [(n, 11) for n in ARMS])
def test_every_column_has_one_anchor_and_one_cell_per_mode_and_effect(name: str, n_columns: int) -> None:
    config = _load(name)
    columns = _columns(config)
    assert len(columns) == n_columns
    nonzero = [float(e) for e in config.effect_sizes if e != 0.0]
    expected = Counter({("none", 0.0): 1, **{(m, e): 1 for m in config.trajectory_modes for e in nonzero}})
    for cells in columns.values():
        assert Counter(_mode_effect(c) for c in cells) == expected
    # The baseline column is the primary grid; every other column is a design-point grid.
    phases = {key: {c.phase for c in cells} for key, cells in columns.items()}
    assert sum(p == {"power_primary"} for p in phases.values()) == 1
    assert sum(p == {DESIGN_PHASE} for p in phases.values()) == n_columns - 1


@pytest.mark.parametrize("name", CORE)
def test_core_profile_crosses_sample_size_stages_and_blocks(name: str) -> None:
    config = _load(name)
    assert set(config.design_grid.axes) == {N, STAGES, LAYERS}
    assert config.design_grid.axes[N] == (80, 200, 500, 1000)
    assert config.design_grid.axes[STAGES] == (3, 4)
    assert config.design_grid.axes[LAYERS] == (["methylation", "expression"], None)
    assert config.baseline_design_point() == {N: 200, STAGES: 3, LAYERS: ["methylation", "expression"]}
    points = {(p[N], p[STAGES], _frozen(p[LAYERS])) for p in config.design_points()}
    assert len(points) == 16


@pytest.mark.parametrize("name", CORE)
def test_core_block_columns_draw_the_same_generator_seeds(name: str) -> None:
    """Columns differing only in the measured blocks evaluate the same dataset at every replicate index."""

    config = _load(name)
    by_data_point: dict[tuple[Any, ...], dict[Any, SimulationCell]] = {}
    for key, cells in _columns(config).items():
        point = dict(key)
        for cell in cells:
            data_key = (point[N], point[STAGES], *_mode_effect(cell))
            by_data_point.setdefault(data_key, {})[point[LAYERS]] = cell
    assert len(by_data_point) == 8 * (1 + len(config.trajectory_modes) * (len(config.effect_sizes) - 1))
    for data_key, by_layers in by_data_point.items():
        assert set(by_layers) == {("methylation", "expression"), None}, data_key
        two, three = by_layers[("methylation", "expression")], by_layers[None]
        assert two.cell_id != three.cell_id
        assert two.generator_params == three.generator_params
        assert [derive_replicate_seed(two, i) for i in range(two.n_replicates)] == [
            derive_replicate_seed(three, i) for i in range(three.n_replicates)
        ]


@pytest.mark.parametrize("name", ARMS)
def test_arm_columns_move_exactly_one_setting_off_the_baseline(name: str) -> None:
    config = _load(name)
    assert config.design_grid.explicit
    baseline = config.baseline_design_point()
    assert baseline == {
        "generator.group_ratio": 0.5,
        "generator.stage_sample_prop": None,
        "generator.p_dmp": 0.1,
        **{axis: 2.0 for axis in DELTAS},
        "generator.baseline_continuity": 0.0,
    }
    moved: list[tuple[str, Any]] = []
    for point in config.design_points():
        if config.is_baseline_design_point(point):
            continue
        changed = {axis for axis, value in point.items() if _frozen(value) != _frozen(baseline[axis])}
        if changed & set(DELTAS):
            assert changed == set(DELTAS), f"a stage-effect column must move all three deltas: {point}"
            values = {point[axis] for axis in DELTAS}
            assert len(values) == 1, f"the stage effect moves jointly: {point}"
            moved.append(("stage_effect", values.pop()))
        else:
            assert len(changed) == 1, f"an arm column moves exactly one setting: {point}"
            (axis,) = changed
            moved.append((axis, _frozen(point[axis])))
    assert sorted(moved, key=repr) == sorted(EXPECTED_ARMS, key=repr)


# ── Headroom ─────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name", PROFILES)
def test_every_cell_fits_the_surgery_headroom(name: str) -> None:
    cells = enumerate_study(_load(name)).cells  # enumeration itself refuses an over-headroom cell
    pool_limited = [c for c in cells if c.generator_params.trajectory_mode in POOL_LIMITED_MODES]
    # The Type I baseline always carries the translation control; only the
    # geometry profiles put pool-limited surgeries in their design columns.
    in_columns = [c for c in pool_limited if c.phase in ("power_primary", DESIGN_PHASE)]
    assert bool(in_columns) == (name in GEOMETRY)
    for cell in pool_limited:
        headroom = expected_surgery_headroom(cell.generator_params)
        assert headroom is not None and headroom.fits, cell.cell_id
