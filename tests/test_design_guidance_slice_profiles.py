"""Structure of the four committed design-guidance paper-grade profiles.

Checks the spec scenarios of ``design-guidance-paper-grade-slice``: the slice
profiles cross sample size x one stage count x blocks into six anchored columns
whose block columns share data and whose design cells reproduce the pilot core
cells' seeds at replicate indices 0-99; the stage-effect profiles move the
joint stage effect through three explicit points off a delta = 0.25 baseline
with every other setting equal to the slice baseline's; every profile names the
production measurement, runs 500 x 999, predeclares the size rules only, and
fits the surgery headroom. Identity (cell ids, signatures, seeds) is pinned
separately by ``tests/test_post_floor_config_identity.py``.
"""

from __future__ import annotations

import dataclasses
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

import pytest

from motco.simulations.grid import SimulationCell, derive_replicate_seed
from motco.simulations.semisynthetic import POOL_LIMITED_MODES, expected_surgery_headroom
from motco.simulations.study import StudyConfig, enumerate_study, load_study_config
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"

SLICE = ("design_guidance_core_magnitude_study.json", "design_guidance_core_geometry_study.json")
ARM = ("design_guidance_stage_effect_magnitude_study.json", "design_guidance_stage_effect_geometry_study.json")
MAGNITUDE = ("design_guidance_core_magnitude_study.json", "design_guidance_stage_effect_magnitude_study.json")
GEOMETRY = ("design_guidance_core_geometry_study.json", "design_guidance_stage_effect_geometry_study.json")
PROFILES = SLICE + ARM
#: The pilot profile each slice column reproduces, and the pilot arms profile each arm column reproduces.
PILOT_OF = {
    "design_guidance_core_magnitude_study.json": "design_guidance_core_magnitude_pilot.json",
    "design_guidance_core_geometry_study.json": "design_guidance_core_geometry_pilot.json",
    "design_guidance_stage_effect_magnitude_study.json": "design_guidance_arms_magnitude_pilot.json",
    "design_guidance_stage_effect_geometry_study.json": "design_guidance_arms_geometry_pilot.json",
}
SLICE_OF = {ARM[0]: SLICE[0], ARM[1]: SLICE[1]}

N = "generator.n_samples"
STAGES = "generator.n_stages"
LAYERS = "evaluation.integration_params.layers"
TWO = ["methylation", "expression"]
DELTAS = ("generator.delta_methyl", "generator.delta_expr", "generator.delta_protein")
STAGE_EFFECTS = (0.25, 0.5, 1.0)

#: The predeclared rules, exactly as written in every profile (task 4.1 / design D4).
EXPECTED_ACCEPTANCE: dict[str, Any] = {
    "type_i": [
        {"alpha": 0.05, "se_tolerance": 2.0},
        {"alpha": 0.05, "se_tolerance": 2.0, "bound": "upper", "trajectory_mode": "translation", "statistic": "angle"},
    ],
    "power": [],
    "specificity": [],
    "gate": {"enabled": False},
    "anchor_decision": {"alpha": 0.05, "se_multiplier": 2.0, "statistics": ["delta", "angle", "shape"]},
}


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


def _seeds(cell: SimulationCell, n: int) -> list[int]:
    return [derive_replicate_seed(cell, i) for i in range(n)]


# ── Profiles name the production measurement, at paper grade ─────────────────


@pytest.mark.parametrize("name", PROFILES)
def test_profile_names_the_production_measurement_at_paper_grade(name: str) -> None:
    raw = _raw(name)
    assert raw["evaluation"]["integration_method"] == "pls"
    assert raw["evaluation"]["orientation_estimator"] == "net_displacement"
    assert raw["evaluation"]["min_components"] == 2
    assert raw["generator"]["magnitude_kind"] == "joint"
    assert raw["attribution"] == {"enabled": False}

    config = _load(name)
    assert config.evaluation.orientation_estimator == "net_displacement"
    assert config.evaluation.min_components == 2
    assert config.generator.magnitude_kind == "joint"
    assert config.n_replicates == 500 and config.evaluation.permutations == 999
    assert config.base_seed == 1028 and config.matched_seeds.primary_family == "design-guidance"
    assert config.matched_seeds.enabled and config.matched_seeds.shared_zero_effect_anchor
    assert config.evaluation.seed == 3
    assert config.evaluation.integration_params["layers"] == TWO
    assert config.generator.surgery_censoring == "error"


@pytest.mark.parametrize("name", SLICE)
def test_slice_baseline_is_the_shared_design_guidance_baseline(name: str) -> None:
    generator = _load(name).generator
    assert (generator.n_samples, generator.n_stages, generator.group_ratio, generator.p_dmp) == (200, 3, 0.5, 0.1)
    assert generator.baseline_continuity == 0.0 and generator.stage_sample_prop is None
    assert (generator.delta_methyl, generator.delta_expr, generator.delta_protein) == (2.0, 2.0, 2.0)
    pilot = _load(PILOT_OF[name])
    # The slice is the pilot profile at paper grade: same generator, same evaluation but the permutation count.
    assert generator == pilot.generator
    assert dataclasses.replace(pilot.evaluation, permutations=999) == _load(name).evaluation


# ── Mode and effect lists ────────────────────────────────────────────────────


@pytest.mark.parametrize("name", MAGNITUDE)
def test_magnitude_profiles_run_magnitude_between_the_resolvable_and_the_saturated_effects(name: str) -> None:
    config = _load(name)
    assert config.trajectory_modes == ("magnitude",)
    assert config.effect_sizes == (0.0, 0.02, 0.05, 0.1)


@pytest.mark.parametrize("name", GEOMETRY)
def test_geometry_profiles_run_the_geometry_modes_and_translation_at_three_effects(name: str) -> None:
    config = _load(name)
    assert config.trajectory_modes == ("orientation", "shape", "translation")
    assert config.effect_sizes == (0.0, 0.25, 0.5, 1.0)


# ── Column structure ─────────────────────────────────────────────────────────


@pytest.mark.parametrize(("name", "n_columns"), [(n, 6) for n in SLICE] + [(n, 3) for n in ARM])
def test_every_column_has_one_anchor_and_one_cell_per_mode_and_effect(name: str, n_columns: int) -> None:
    config = _load(name)
    columns = _columns(config)
    assert len(columns) == n_columns
    nonzero = [float(e) for e in config.effect_sizes if e != 0.0]
    expected = Counter({("none", 0.0): 1, **{(m, e): 1 for m in config.trajectory_modes for e in nonzero}})
    for cells in columns.values():
        assert Counter(_mode_effect(c) for c in cells) == expected
        assert all(c.n_replicates == 500 and c.evaluation_params.permutations == 999 for c in cells)
    phases = {key: {c.phase for c in cells} for key, cells in columns.items()}
    assert sum(p == {"power_primary"} for p in phases.values()) == 1
    assert sum(p == {DESIGN_PHASE} for p in phases.values()) == n_columns - 1
    type_i = [c for c in enumerate_study(config).cells if c.phase == "type_i_baseline"]
    assert sorted(c.generator_params.trajectory_mode for c in type_i) == ["none", "translation"]


@pytest.mark.parametrize(
    ("name", "n_cells", "n_units"),
    [(SLICE[0], 26, 13_000), (SLICE[1], 62, 31_000), (ARM[0], 14, 7_000), (ARM[1], 32, 16_000)],
)
def test_cell_and_unit_counts(name: str, n_cells: int, n_units: int) -> None:
    cells = enumerate_study(_load(name)).cells
    assert len(cells) == n_cells
    assert sum(c.n_replicates for c in cells) == n_units


@pytest.mark.parametrize("name", SLICE)
def test_slice_crosses_sample_size_one_stage_count_and_blocks(name: str) -> None:
    config = _load(name)
    assert not config.design_grid.explicit
    assert set(config.design_grid.axes) == {N, STAGES, LAYERS}
    assert config.design_grid.axes[N] == (80, 200, 500)
    assert config.design_grid.axes[STAGES] == (3,)
    assert config.design_grid.axes[LAYERS] == (TWO, None)
    assert config.baseline_design_point() == {N: 200, STAGES: 3, LAYERS: TWO}
    assert len({(p[N], p[STAGES], _frozen(p[LAYERS])) for p in config.design_points()}) == 6


@pytest.mark.parametrize("name", SLICE)
def test_slice_block_columns_draw_the_same_generator_seeds_at_every_replicate(name: str) -> None:
    config = _load(name)
    by_data_point: dict[tuple[Any, ...], dict[Any, SimulationCell]] = {}
    for key, cells in _columns(config).items():
        point = dict(key)
        for cell in cells:
            by_data_point.setdefault((point[N], *_mode_effect(cell)), {})[point[LAYERS]] = cell
    assert len(by_data_point) == 3 * (1 + len(config.trajectory_modes) * (len(config.effect_sizes) - 1))
    for data_key, by_layers in by_data_point.items():
        assert set(by_layers) == {_frozen(TWO), None}, data_key
        two, three = by_layers[_frozen(TWO)], by_layers[None]
        assert two.cell_id != three.cell_id
        assert two.generator_params == three.generator_params
        assert _seeds(two, 500) == _seeds(three, 500)


@pytest.mark.parametrize("name", SLICE + ARM)
def test_design_cells_reproduce_the_pilot_seeds_at_the_first_hundred_replicates(name: str) -> None:
    """Every primary / design-point cell draws the pilot cell's seeds at 0-99; the Type I baseline cells do not."""

    config, pilot = _load(name), _load(PILOT_OF[name])
    axes = (N, STAGES, LAYERS) if name in SLICE else DELTAS

    #: The pilot arms profiles' columns that are not stage-effect columns (another factor moved).
    other_arms = {
        "generator.group_ratio": 0.5,
        "generator.stage_sample_prop": None,
        "generator.p_dmp": 0.1,
        "generator.baseline_continuity": 0.0,
    }

    def keyed(cfg: StudyConfig) -> dict[tuple[Any, ...], SimulationCell]:
        out: dict[tuple[Any, ...], SimulationCell] = {}
        for key, cells in _columns(cfg).items():
            point = dict(key)
            if any(point.get(axis, value) != value for axis, value in other_arms.items()):
                continue
            for cell in cells:
                out[tuple(point[a] for a in axes) + _mode_effect(cell)] = cell
        return out

    study_cells, pilot_cells = keyed(config), keyed(pilot)
    assert len(study_cells) == {SLICE[0]: 24, SLICE[1]: 60, ARM[0]: 12, ARM[1]: 30}[name]
    for key, cell in study_cells.items():
        assert key in pilot_cells, key
        assert _seeds(cell, 100) == _seeds(pilot_cells[key], 100), key
        assert cell.n_replicates == 500 and pilot_cells[key].n_replicates == 100

    def type_i(cfg: StudyConfig) -> dict[str, SimulationCell]:
        cells = enumerate_study(cfg).cells
        return {c.generator_params.trajectory_mode: c for c in cells if c.phase == "type_i_baseline"}

    for mode, cell in type_i(config).items():
        assert _seeds(cell, 3) != _seeds(type_i(pilot)[mode], 3), mode


# ── Predeclared rules ────────────────────────────────────────────────────────


@pytest.mark.parametrize("name", PROFILES)
def test_acceptance_block_predeclares_the_size_rules_only(name: str) -> None:
    assert _raw(name)["acceptance"] == EXPECTED_ACCEPTANCE
    assert "report_contract" not in _raw(name)
    acceptance = _load(name).acceptance
    phase_scoped, column_scoped = acceptance.type_i
    assert not phase_scoped.column_scoped and phase_scoped.bound == "two_sided"
    assert column_scoped.column_scoped
    assert column_scoped.trajectory_mode == "translation"
    assert (column_scoped.statistic, column_scoped.bound) == ("angle", "upper")
    assert column_scoped.alpha == 0.05 and column_scoped.se_tolerance == 2.0
    assert acceptance.power == () and acceptance.specificity == ()
    assert not acceptance.gate.enabled
    assert acceptance.anchor_decision is not None
    assert (acceptance.anchor_decision.alpha, acceptance.anchor_decision.se_multiplier) == (0.05, 2.0)
    assert acceptance.anchor_decision.statistics == ("delta", "angle", "shape")
    assert acceptance.design_point is None and acceptance.rank_decision is None and acceptance.column_decision is None
    assert _load(name).report_contract is None
    assert not _load(name).attribution.enabled


# ── Stage-effect arm ─────────────────────────────────────────────────────────


@pytest.mark.parametrize("name", ARM)
def test_arm_columns_move_the_joint_stage_effect_off_a_quarter_baseline(name: str) -> None:
    config = _load(name)
    assert config.design_grid.explicit
    assert set(config.design_grid.axes) == set(DELTAS)
    assert config.baseline_design_point() == {axis: 0.25 for axis in DELTAS}
    values = []
    for point in config.design_points():
        assert set(point) == set(DELTAS)
        assert len({point[axis] for axis in DELTAS}) == 1, f"the stage effect moves jointly: {point}"
        values.append(point[DELTAS[0]])
    assert sorted(values) == list(STAGE_EFFECTS)

    # Every other setting equals the slice baseline's: same generator but the three deltas, same evaluation.
    slice_config = _load(SLICE_OF[name])
    assert config.generator == dataclasses.replace(
        slice_config.generator, delta_methyl=0.25, delta_expr=0.25, delta_protein=0.25
    )
    assert config.evaluation == slice_config.evaluation
    assert (config.trajectory_modes, config.effect_sizes) == (slice_config.trajectory_modes, slice_config.effect_sizes)
    assert (config.n_replicates, config.base_seed, config.matched_seeds) == (
        slice_config.n_replicates, slice_config.base_seed, slice_config.matched_seeds
    )
    assert config.acceptance == slice_config.acceptance
    for point, cells in _columns(config).items():
        for cell in cells:
            assert (cell.generator_params.n_samples, cell.generator_params.n_stages) == (200, 3)
            assert cell.evaluation_params.integration_params["layers"] == TWO
            assert cell.generator_params.delta_methyl == dict(point)[DELTAS[0]]
    assert "second priority" in _raw(name)["metadata"]["intent"]
    assert any("Baseline deviation" in note for note in _raw(name)["metadata"]["notes"])


# ── Headroom and pinning ─────────────────────────────────────────────────────


@pytest.mark.parametrize("name", PROFILES)
def test_every_cell_fits_the_surgery_headroom(name: str) -> None:
    cells = enumerate_study(_load(name)).cells  # enumeration itself refuses an over-headroom cell
    pool_limited = [c for c in cells if c.generator_params.trajectory_mode in POOL_LIMITED_MODES]
    in_columns = [c for c in pool_limited if c.phase in ("power_primary", DESIGN_PHASE)]
    assert bool(in_columns) == (name in GEOMETRY)
    for cell in pool_limited:
        headroom = expected_surgery_headroom(cell.generator_params)
        assert headroom is not None and headroom.fits, cell.cell_id


def test_profiles_are_post_floor_configs() -> None:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from post_floor_fixture_capture import POST_FLOOR_CONFIGS
    finally:
        sys.path.pop(0)
    assert set(PROFILES) <= set(POST_FLOOR_CONFIGS)
    for name in PROFILES:
        notes = " ".join(_raw(name)["metadata"]["notes"])
        assert "post-floor config" in notes and "net_displacement" in notes
