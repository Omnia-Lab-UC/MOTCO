"""The committed Phase 6 small-n pilot profile reproduces the SEA-AD case-study design."""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import pytest

from motco.simulations.grid import derive_replicate_seed, run_simulation_replicate
from motco.simulations.semisynthetic import generate_semisynthetic_trajectory
from motco.simulations.study import enumerate_study, load_study_config

PILOT = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study" / "phase6_small_n_pilot.json"
PILOT_GRID = (0.0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0)
SEA_AD = ((11, 10, 28), (9, 10, 12))
MEASURED = ["methylation", "expression"]


def test_pilot_declares_the_case_study_design() -> None:
    raw = json.loads(PILOT.read_text(encoding="utf-8"))
    for proportional in ("n_samples", "stage_sample_prop", "group_ratio", "surgery_censoring"):
        assert proportional not in raw["generator"]
    assert "design_grid" not in raw and raw["axes"] == {}

    config = load_study_config(PILOT)
    assert config.generator.n_stages == 3
    assert config.generator.group_stage_sizes == SEA_AD
    assert config.generator.magnitude_kind == "joint"
    assert config.generator.surgery_censoring == "error"
    assert config.evaluation.integration_method == "pls"
    assert config.evaluation.integration_params["layers"] == MEASURED
    assert "forced_components" not in config.evaluation.integration_params
    assert config.trajectory_modes == ("magnitude", "orientation", "shape", "translation")
    assert config.effect_sizes == PILOT_GRID
    assert (config.n_replicates, config.evaluation.permutations) == (100, 199)
    assert config.matched_seeds.enabled and config.matched_seeds.shared_zero_effect_anchor
    assert config.matched_seeds.primary_family == "phase6-small-n"
    assert not config.acceptance.gate.enabled


def test_pilot_enumerates_under_headroom_with_the_declared_cell_sizes() -> None:
    config = load_study_config(PILOT)
    grid = enumerate_study(config)  # rejects any over-headroom cell
    assert len(grid.cells) == 2 + 1 + 4 * 8
    anchors = [cell for cell in grid.cells if cell.metadata.get("zero_effect_anchor")]
    assert len(anchors) == 1
    cell = grid.cells[0]
    params = replace(cell.generator_params, seed=derive_replicate_seed(cell, 0))
    counts = generate_semisynthetic_trajectory(params).metadata.groupby(["group", "stage"]).size()
    assert [[int(counts[(g, s)]) for s in range(3)] for g in ("A", "B")] == [list(row) for row in SEA_AD]


PILOT_CELLS = enumerate_study(load_study_config(PILOT)).cells


@pytest.mark.slow
@pytest.mark.parametrize("cell", PILOT_CELLS, ids=[cell.cell_id for cell in PILOT_CELLS])
def test_every_pilot_cell_runs_end_to_end(cell) -> None:
    """One replicate of every cell at a few permutations, including CV at the 9-sample cells.

    One test per cell, so pytest-xdist spreads the cells over its workers.
    """

    fast = replace(cell, evaluation_params=replace(cell.evaluation_params, permutations=9))
    record = run_simulation_replicate(fast, 0)
    assert record.status == "completed", (cell.cell_id, record.error_message)
    integration = record.integration_metadata
    assert integration["component_selection"] == "cv"
    assert integration["integration_params"]["layers"] == MEASURED
    assert integration["layer_feature_counts"] == {"methylation": 367, "expression": 131}
    assert record.truth_metadata["group_stage_sizes"] == {"A": [11, 10, 28], "B": [9, 10, 12]}
    assert set(record.p_values) == {"delta", "angle", "shape"}


# ── Paper-grade split profiles ───────────────────────────────────────────────

CONFIG_DIR = PILOT.parent
MAGNITUDE = CONFIG_DIR / "phase6_small_n_magnitude.json"
STUDY = CONFIG_DIR / "phase6_small_n_study.json"
SEA_AD_50 = ((10, 9, 25), (9, 9, 12))
SIZES = "generator.group_stage_sizes"
LAYERS = "evaluation.integration_params.layers"
PILOT_DIR = "results/phase6-small-n-pilot-2026-10-09"


@pytest.mark.parametrize(
    ("path", "modes", "effects", "n_cells"),
    [
        (MAGNITUDE, ("magnitude",), (0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0), 2 + 4 * (1 + 6)),
        (STUDY, ("orientation", "shape", "translation"), (0.0, 0.25, 0.5, 0.75, 1.0), 2 + 4 * (1 + 3 * 4)),
    ],
    ids=["magnitude", "orientation-shape-translation"],
)
def test_paper_grade_profiles_follow_the_pilot(path, modes, effects, n_cells) -> None:
    raw = json.loads(path.read_text(encoding="utf-8"))
    pilot = json.loads(PILOT.read_text(encoding="utf-8"))
    assert raw["generator"] == pilot["generator"]
    evaluation = dict(raw["evaluation"])
    assert evaluation.pop("permutations") == 999
    assert evaluation == {key: value for key, value in pilot["evaluation"].items() if key != "permutations"}
    assert (raw["base_seed"], raw["matched_seeds"]) == (pilot["base_seed"], pilot["matched_seeds"])
    assert raw["report_contract"] == {
        "driver_component": "observed",
        "cross_replicate_driver_agreement": "descriptive",
        "n_jobs_override": "forbid",
    }
    assert raw["metadata"]["pilot"]["directory"] == PILOT_DIR
    assert set(raw["metadata"]["effect_axis"]) >= {str(effect) for effect in effects}

    config = load_study_config(path)
    assert config.trajectory_modes == modes
    assert config.effect_sizes == effects
    assert (config.n_replicates, config.evaluation.permutations) == (500, 999)
    assert config.acceptance.gate.enabled
    assert config.design_grid.axes == {
        SIZES: (SEA_AD, SEA_AD_50),
        LAYERS: (MEASURED, None),
    }
    grid = enumerate_study(config)  # rejects any over-headroom cell
    assert len(grid.cells) == n_cells
    assert grid.metadata["design_grid"]["n_points"] == 4
    assert sum(1 for cell in grid.cells if cell.metadata.get("zero_effect_anchor")) == 4


def test_split_profiles_and_pilot_share_their_anchor_datasets() -> None:
    def anchor(path):
        grid = enumerate_study(load_study_config(path))
        (cell,) = [c for c in grid.cells if c.metadata.get("zero_effect_anchor") and c.phase == "power_primary"]
        return cell

    anchors = [anchor(path) for path in (PILOT, MAGNITUDE, STUDY)]
    reference = anchors[0]
    for cell in anchors[1:]:
        assert cell.generator_params == reference.generator_params
        for index in range(3):
            assert derive_replicate_seed(cell, index) == derive_replicate_seed(reference, index)


def test_design_columns_realize_their_cell_sizes() -> None:
    grid = enumerate_study(load_study_config(STUDY))
    for cell in grid.cells:
        if not cell.metadata.get("zero_effect_anchor"):
            continue
        params = replace(cell.generator_params, seed=derive_replicate_seed(cell, 0))
        counts = generate_semisynthetic_trajectory(params).metadata.groupby(["group", "stage"]).size()
        expected = cell.generator_params.group_stage_sizes
        assert [[int(counts[(g, s)]) for s in range(3)] for g in ("A", "B")] == [list(row) for row in expected]
        assert int(counts.sum()) in (80, 74)
