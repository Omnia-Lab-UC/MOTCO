"""The committed net-displacement Type I profiles: structure, rule, shared datasets.

Both profiles are null-only (the shared zero-effect anchor and the translation
control at e = 1.0) across five explicit measurement columns. Their baseline
column sits in the parent study's matched-seed family, so it re-evaluates the
parent's datasets at every replicate index.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from motco.simulations.grid import derive_replicate_seed
from motco.simulations.semisynthetic import generate_semisynthetic_trajectory
from motco.simulations.study import AnchorDecisionRule, enumerate_study, load_study_config
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
PHASE6 = CONFIG_DIR / "net_displacement_type_i_phase6.json"
PHASE5 = CONFIG_DIR / "net_displacement_type_i_phase5.json"
PHASE6_PARENT = CONFIG_DIR / "phase6_small_n_study.json"
PHASE5_PARENT = CONFIG_DIR / "phase5_power_study.json"

METHOD = "evaluation.integration_method"
FORCED = "evaluation.integration_params.forced_components"
NCOMP = "evaluation.integration_params.n_components"
EST = "evaluation.orientation_estimator"


def _points(rank: int) -> list[dict]:
    return [
        {METHOD: "pls", FORCED: None, NCOMP: None, EST: None},
        {METHOD: "pls", FORCED: None, NCOMP: None, EST: "net_displacement"},
        {METHOD: "pls", FORCED: rank, NCOMP: None, EST: "net_displacement"},
        {METHOD: "pca", FORCED: None, NCOMP: rank, EST: "net_displacement"},
        {METHOD: "concat", FORCED: None, NCOMP: None, EST: "net_displacement"},
    ]


def _primary(path: Path, mode: str, effect: float):
    (cell,) = [
        cell
        for cell in enumerate_study(load_study_config(path)).cells
        if cell.phase == "power_primary"
        and cell.metadata.get("trajectory_mode") == mode
        and cell.metadata.get("effect_size") == effect
    ]
    return cell


@pytest.mark.parametrize(
    ("path", "rank", "stages", "family", "base_seed"),
    [(PHASE6, 4, 3, "phase6-small-n", 800), (PHASE5, 6, 4, "phase5-primary", 600)],
)
def test_profile_enumerates_twelve_null_cells(
    path: Path, rank: int, stages: int, family: str, base_seed: int
) -> None:
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert "surgery_censoring" not in raw["generator"]
    config = load_study_config(path)
    assert config.generator.surgery_censoring == "error"
    assert config.generator.n_stages == stages
    assert config.generator.p_dmp == 0.1 and config.generator.baseline_continuity == 0.0
    assert config.generator.magnitude_kind == "joint"
    assert config.trajectory_modes == ("translation",) and config.effect_sizes == (0.0, 1.0)
    assert (config.n_replicates, config.evaluation.permutations) == (500, 999)
    assert config.base_seed == base_seed and config.matched_seeds.primary_family == family
    assert config.matched_seeds.shared_zero_effect_anchor
    assert not config.attribution.enabled and not config.acceptance.gate.enabled
    assert config.acceptance.power == ()
    assert config.acceptance.anchor_decision == AnchorDecisionRule(alpha=0.05, se_multiplier=2.0)
    assert config.design_points() == _points(rank)

    grid = enumerate_study(config)  # rejects any over-headroom cell
    assert len(grid.cells) == 12
    phases = [cell.phase for cell in grid.cells]
    assert phases.count("type_i_baseline") == 2 and phases.count("power_primary") == 2
    design = [cell for cell in grid.cells if cell.phase == DESIGN_PHASE]
    columns: dict[str, list] = {}
    for cell in design:
        columns.setdefault(json.dumps(cell.metadata[DESIGN_POINT_KEY], sort_keys=True), []).append(cell)
    assert len(columns) == 4
    for cells in columns.values():
        assert sorted((c.metadata["trajectory_mode"], c.metadata["effect_size"]) for c in cells) == [
            ("none", 0.0),
            ("translation", 1.0),
        ]
    assert all(cell.n_replicates == 500 and cell.evaluation_params.permutations == 999 for cell in grid.cells)


def test_phase6_profile_layout() -> None:
    config = load_study_config(PHASE6)
    assert config.generator.group_stage_sizes == ((11, 10, 28), (9, 10, 12))
    assert config.evaluation.integration_params["layers"] == ["methylation", "expression"]


def test_phase5_profile_layout() -> None:
    config = load_study_config(PHASE5)
    assert config.generator.n_samples == 1200 and config.generator.group_ratio == 0.5
    assert "layers" not in config.evaluation.integration_params


@pytest.mark.parametrize(("path", "parent"), [(PHASE6, PHASE6_PARENT), (PHASE5, PHASE5_PARENT)])
@pytest.mark.parametrize(("mode", "effect"), [("none", 0.0), ("translation", 1.0)])
def test_baseline_cells_draw_the_parent_study_seeds(path: Path, parent: Path, mode: str, effect: float) -> None:
    cell = _primary(path, mode, effect)
    parent_cell = _primary(parent, mode, effect)
    assert cell.evaluation_params == parent_cell.evaluation_params
    assert [derive_replicate_seed(cell, i) for i in range(cell.n_replicates)] == [
        derive_replicate_seed(parent_cell, i) for i in range(parent_cell.n_replicates)
    ]


def _dataset_digest(cell, index: int) -> str:
    dataset = generate_semisynthetic_trajectory(
        replace(cell.generator_params, seed=derive_replicate_seed(cell, index))
    )
    digest = hashlib.sha256()
    for frame in (dataset.methylation, dataset.expression, dataset.proteomics, dataset.metadata):
        digest.update(pd.util.hash_pandas_object(frame, index=True).to_numpy().tobytes())
    return digest.hexdigest()


@pytest.mark.slow
def test_phase5_anchor_datasets_equal_the_paper_grade_anchor() -> None:
    """Six n = 1200 datasets (about two minutes): seeds equal is not enough at Phase 5,
    because the two configs differ in ``magnitude_kind``, which must be inert at e = 0."""

    cell = _primary(PHASE5, "none", 0.0)
    parent_cell = _primary(PHASE5_PARENT, "none", 0.0)
    assert cell.generator_params.magnitude_kind != parent_cell.generator_params.magnitude_kind
    for index in range(3):
        assert _dataset_digest(cell, index) == _dataset_digest(parent_cell, index)
