"""The three stage-span paper-grade profiles: structure, shared anchors and reproduction of earlier runs.

``stage_span_paper_grade_{span_geometry,span_magnitude,surgery_geometry}.json``
run the planted-scale stage-span construction at five columns and the surgeries
at the two columns no earlier paper-grade run covers, 500 x 999 in the
``design-guidance`` family. Identity (cell ids, signatures, seeds) is pinned by
``tests/data/post_floor_config_identity.json``; this module checks the
structure the decision and the reproduction checks rely on.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest

from motco.simulations.grid import SimulationCell, derive_replicate_seed
from motco.simulations.study import StudyConfig, enumerate_study, load_study_config

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
TWO = ("methylation", "expression")
SEA_AD = ((11, 10, 28), (9, 10, 12))
SPAN_GEOMETRY = "stage_span_paper_grade_span_geometry.json"
SPAN_MAGNITUDE = "stage_span_paper_grade_span_magnitude.json"
SURGERY_GEOMETRY = "stage_span_paper_grade_surgery_geometry.json"
PROFILES = (SPAN_GEOMETRY, SPAN_MAGNITUDE, SURGERY_GEOMETRY)
COLUMNS = ("n200", "n80", "seaad", "n200_three_blocks", "n200_four_stages")
EFFECTS = {
    SPAN_GEOMETRY: {(m, e) for m in ("orientation", "shape", "translation") for e in (0.1, 0.25, 0.5)},
    SPAN_MAGNITUDE: {("magnitude", e) for e in (0.02, 0.05, 0.1)},
    SURGERY_GEOMETRY: {(m, e) for m in ("orientation", "shape", "translation") for e in (0.25, 0.5, 1.0)},
}
REPLICATES = 100


def _load(name: str) -> StudyConfig:
    return load_study_config(CONFIG_DIR / name)


def _column(cell: SimulationCell) -> str:
    generator = cell.generator_params
    name = "seaad" if generator.group_stage_sizes is not None else f"n{generator.n_samples}"
    if cell.evaluation_params.integration_params.get("layers") is None:
        name += "_three_blocks"
    if generator.n_stages == 4:
        name += "_four_stages"
    return name


def _cells(name: str) -> dict[tuple[str, str, float], SimulationCell]:
    cells = {}
    for cell in enumerate_study(_load(name)).cells:
        if cell.phase not in ("power_primary", "power_design"):
            continue
        key = (_column(cell), str(cell.metadata["trajectory_mode"]), float(cell.metadata["effect_size"]))
        assert key not in cells, key
        cells[key] = cell
    return cells


def _seeds(cell: SimulationCell) -> list[int]:
    return [derive_replicate_seed(cell, index) for index in range(REPLICATES)]


def _same_data(cell: SimulationCell, other: SimulationCell, *, ignore_scale: bool = False) -> None:
    generator = replace(cell.generator_params, span_scale=None) if ignore_scale else cell.generator_params
    assert generator == other.generator_params
    assert replace(cell.evaluation_params, permutations=0, span_capture=None) == replace(
        other.evaluation_params, permutations=0, span_capture=None
    )
    assert _seeds(cell) == _seeds(other)


# ── Structure ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("name", PROFILES)
def test_profile_runs_the_production_measurement_at_paper_grade(name: str) -> None:
    config = _load(name)
    assert config.n_replicates == 500 and config.evaluation.permutations == 999
    assert config.base_seed == 1028 and config.matched_seeds.primary_family == "design-guidance"
    assert config.matched_seeds.enabled and config.matched_seeds.shared_zero_effect_anchor
    assert config.evaluation.integration_method == "pls"
    assert config.evaluation.orientation_estimator == "net_displacement"
    assert config.evaluation.min_components == 2
    assert config.evaluation.span_capture is True
    assert config.evaluation.n_jobs == 1 and config.evaluation.seed == 3
    assert config.generator.magnitude_kind == "joint" and config.generator.p_dmp == 0.1
    assert config.acceptance.anchor_decision is not None
    assert config.metadata["submission"]["template"] == "scripts/motco_study_array.sbatch"
    if name == SURGERY_GEOMETRY:
        assert config.generator.construction is None
        assert config.generator.span_layers is None and config.generator.span_scale is None
    else:
        assert config.generator.construction == "stage_span"
        assert config.generator.span_scale == "planted"


@pytest.mark.parametrize(
    ("name", "columns", "n_cells"),
    [
        (SPAN_GEOMETRY, COLUMNS, 52),
        (SPAN_MAGNITUDE, COLUMNS, 22),
        (SURGERY_GEOMETRY, ("seaad", "n200_four_stages"), 22),
    ],
)
def test_profile_enumerates_its_columns_effects_and_counts(name: str, columns: tuple[str, ...], n_cells: int) -> None:
    cells = _cells(name)
    assert {key[0] for key in cells} == set(columns)
    for column in columns:
        assert {key[1:] for key in cells if key[0] == column} == {("none", 0.0)} | EFFECTS[name]
    everything = enumerate_study(_load(name)).cells
    assert len(everything) == n_cells
    assert sum(cell.n_replicates for cell in everything) == 500 * n_cells


@pytest.mark.parametrize("name", [SPAN_GEOMETRY, SPAN_MAGNITUDE])
def test_stage_span_columns_measure_their_span_layers(name: str) -> None:
    for (column, mode, _), cell in _cells(name).items():
        generator, layers = cell.generator_params, cell.evaluation_params.integration_params.get("layers")
        if mode == "none":
            assert generator.construction is None and generator.span_scale is None
            continue
        assert generator.construction == "stage_span" and generator.span_scale == "planted"
        if column == "n200_three_blocks":
            assert layers is None and generator.span_layers is None
        else:
            assert tuple(layers) == TWO and generator.span_layers == TWO
        assert generator.n_stages == (4 if column == "n200_four_stages" else 3)
        if column == "seaad":
            assert generator.group_stage_sizes == SEA_AD and generator.n_samples == 600


def test_the_magnitude_profile_differs_from_the_geometry_profile_only_in_modes_and_effects() -> None:
    geometry = json.loads((CONFIG_DIR / SPAN_GEOMETRY).read_text())
    magnitude = json.loads((CONFIG_DIR / SPAN_MAGNITUDE).read_text())
    for key in set(geometry) | set(magnitude):
        if key in ("trajectory_modes", "effect_sizes", "metadata"):
            continue
        assert geometry[key] == magnitude[key], key


# ── Shared anchors ──────────────────────────────────────────────────────────


def test_the_three_profiles_anchor_on_the_same_datasets() -> None:
    anchors: dict[str, list[SimulationCell]] = {}
    for name in PROFILES:
        for (column, mode, _), cell in _cells(name).items():
            if mode == "none":
                anchors.setdefault(column, []).append(cell)
    assert {column: len(cells) for column, cells in anchors.items()} == {
        "n200": 2,
        "n80": 2,
        "seaad": 3,
        "n200_three_blocks": 2,
        "n200_four_stages": 3,
    }
    for cells in anchors.values():
        first = cells[0]
        # Same generator, evaluation and seeds; the cell ids may differ (an anchor records the modes it resolves).
        for other in cells[1:]:
            _same_data(other, first)


def test_effect_cells_share_the_anchor_seeds() -> None:
    for name in PROFILES:
        cells = _cells(name)
        for (column, mode, _), cell in cells.items():
            assert _seeds(cell) == _seeds(cells[(column, "none", 0.0)])


# ── Reproduction of earlier runs ────────────────────────────────────────────


def test_shared_column_anchors_are_the_design_guidance_slice_anchors() -> None:
    slice_cells = {
        key: cell for key, cell in _cells("design_guidance_core_geometry_study.json").items() if key[1] == "none"
    }
    for column in ("n200", "n80", "n200_three_blocks"):
        _same_data(_cells(SPAN_GEOMETRY)[(column, "none", 0.0)], slice_cells[(column, "none", 0.0)])
        assert _cells(SPAN_GEOMETRY)[(column, "none", 0.0)].evaluation_params.permutations == 999
        assert slice_cells[(column, "none", 0.0)].evaluation_params.permutations == 999


def test_seaad_surgery_cells_are_the_stage_span_pilot_surgery_cells() -> None:
    pilot = _cells("stage_span_pilot_surgery_geometry.json")
    for (column, mode, effect), cell in _cells(SURGERY_GEOMETRY).items():
        if column == "seaad":
            _same_data(cell, pilot[(column, mode, effect)])


def test_four_stage_surgery_cells_are_the_design_guidance_pilot_cells() -> None:
    pilot = _cells("design_guidance_core_geometry_pilot.json")
    for (column, mode, effect), cell in _cells(SURGERY_GEOMETRY).items():
        if column == "n200_four_stages":
            _same_data(cell, pilot[(column, mode, effect)])


def test_planted_magnitude_cells_are_the_stage_span_pilot_magnitude_cells() -> None:
    pilot = _cells("stage_span_pilot_span_magnitude.json")
    matched = 0
    for (column, mode, effect), cell in _cells(SPAN_MAGNITUDE).items():
        if column in ("n200", "n80", "seaad"):
            _same_data(cell, pilot[(column, mode, effect)], ignore_scale=True)
            matched += 1
    assert matched == 12
