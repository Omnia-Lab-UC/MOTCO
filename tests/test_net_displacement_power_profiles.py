"""The committed net-displacement power comparison profiles: structure, rule, shared datasets.

Both profiles measure the four trajectory modes over two explicit columns on
matched datasets: ``pls_cv × pc1`` (baseline) and ``pls_cv × net_displacement``
(the candidate). Their baseline column sits in the parent studies' matched-seed
families, so it re-evaluates the parents' datasets at every replicate index.
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
from motco.simulations.study import (
    ColumnDecisionRule,
    StatisticPair,
    TypeIBound,
    enumerate_study,
    load_study_config,
)
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
PHASE6 = CONFIG_DIR / "net_displacement_power_phase6.json"
PHASE5 = CONFIG_DIR / "net_displacement_power_phase5.json"
PHASE6_PARENT = CONFIG_DIR / "phase6_small_n_study.json"
PHASE6_MAGNITUDE_PARENT = CONFIG_DIR / "phase6_small_n_magnitude.json"
PHASE5_PARENT = CONFIG_DIR / "phase5_power_study.json"
PHASE5_MAGNITUDE_PARENT = CONFIG_DIR / "phase5_magnitude_remeasurement.json"

METHOD = "evaluation.integration_method"
EST = "evaluation.orientation_estimator"
MODES = ("magnitude", "orientation", "shape", "translation")
EFFECTS = (0.25, 0.5, 0.75, 1.0)


def _primary(path: Path, mode: str, effect: float):
    (cell,) = [
        cell
        for cell in enumerate_study(load_study_config(path)).cells
        if cell.phase == "power_primary"
        and cell.metadata.get("trajectory_mode") == mode
        and cell.metadata.get("effect_size") == effect
    ]
    return cell


def _rule(criterion: str) -> ColumnDecisionRule:
    return ColumnDecisionRule(
        candidate={METHOD: "pls", EST: "net_displacement"},
        target=StatisticPair("orientation", "angle"),
        criterion=criterion,
        controls=(
            StatisticPair("translation", "angle"),
            StatisticPair("magnitude", "angle"),
            StatisticPair("shape", "angle"),
        ),
        protected=(StatisticPair("magnitude", "delta"), StatisticPair("shape", "shape")),
        type_i_bound=TypeIBound(alpha=0.05, se_tolerance=2.0),
        gain_se_multiplier=2.0,
        loss_se_multiplier=2.0,
    )


@pytest.mark.parametrize(
    ("path", "stages", "family", "base_seed", "criterion"),
    [(PHASE6, 3, "phase6-small-n", 800, "gain"), (PHASE5, 4, "phase5-primary", 600, "no_loss")],
)
def test_profile_enumerates_thirty_six_cells(
    path: Path, stages: int, family: str, base_seed: int, criterion: str
) -> None:
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert "surgery_censoring" not in raw["generator"]
    config = load_study_config(path)
    assert config.generator.surgery_censoring == "error"
    assert config.generator.n_stages == stages
    assert config.generator.p_dmp == 0.1 and config.generator.baseline_continuity == 0.0
    assert config.generator.magnitude_kind == "joint"
    assert config.trajectory_modes == MODES and config.effect_sizes == (0.0, *EFFECTS)
    assert (config.n_replicates, config.evaluation.permutations) == (500, 999)
    assert config.base_seed == base_seed and config.matched_seeds.primary_family == family
    assert config.matched_seeds.shared_zero_effect_anchor
    assert not config.attribution.enabled and not config.acceptance.gate.enabled
    assert config.acceptance.power == ()
    assert config.report_contract is None
    assert config.design_points() == [{METHOD: "pls", EST: None}, {METHOD: "pls", EST: "net_displacement"}]

    targets = (*config.acceptance.type_i, *config.acceptance.specificity)
    assert targets and all(target.bound == "upper" for target in targets)
    assert [(t.trajectory_mode, t.statistic) for t in config.acceptance.specificity] == [
        ("translation", "delta"),
        ("translation", "angle"),
        ("translation", "shape"),
        ("magnitude", "angle"),
        ("magnitude", "shape"),
    ]
    assert config.acceptance.column_decision == _rule(criterion)

    grid = enumerate_study(config)  # rejects any over-headroom cell
    assert len(grid.cells) == 36
    phases = [cell.phase for cell in grid.cells]
    assert phases.count("type_i_baseline") == 2 and phases.count("power_primary") == 17
    expected = sorted([("none", 0.0)] + [(mode, effect) for mode in MODES for effect in EFFECTS])
    primary = [cell for cell in grid.cells if cell.phase == "power_primary"]
    assert sorted((c.metadata["trajectory_mode"], c.metadata["effect_size"]) for c in primary) == expected
    design = [cell for cell in grid.cells if cell.phase == DESIGN_PHASE]
    assert {json.dumps(cell.metadata[DESIGN_POINT_KEY], sort_keys=True) for cell in design} == {
        json.dumps({METHOD: "pls", EST: "net_displacement"}, sort_keys=True)
    }
    assert sorted((c.metadata["trajectory_mode"], c.metadata["effect_size"]) for c in design) == expected
    assert all(cell.n_replicates == 500 and cell.evaluation_params.permutations == 999 for cell in grid.cells)
    assert sum(cell.n_replicates for cell in grid.cells) == 18_000


def test_phase6_profile_layout() -> None:
    config = load_study_config(PHASE6)
    assert config.generator.group_stage_sizes == ((11, 10, 28), (9, 10, 12))
    assert config.evaluation.integration_params["layers"] == ["methylation", "expression"]


def test_phase5_profile_layout() -> None:
    config = load_study_config(PHASE5)
    assert config.generator.n_samples == 1200 and config.generator.group_ratio == 0.5
    assert "layers" not in config.evaluation.integration_params


def _parent_cases() -> list[tuple[Path, Path, str, float]]:
    cases = [(PHASE6, PHASE6_PARENT, "none", 0.0), (PHASE5, PHASE5_PARENT, "none", 0.0)]
    for mode in ("orientation", "shape", "translation"):
        for effect in EFFECTS:
            cases.append((PHASE6, PHASE6_PARENT, mode, effect))
            cases.append((PHASE5, PHASE5_PARENT, mode, effect))
    for effect in (0.25, 0.5, 1.0):
        cases.append((PHASE6, PHASE6_MAGNITUDE_PARENT, "magnitude", effect))
    for effect in (0.25, 1.0):
        cases.append((PHASE5, PHASE5_MAGNITUDE_PARENT, "magnitude", effect))
    return cases


@pytest.mark.parametrize(("path", "parent", "mode", "effect"), _parent_cases())
def test_baseline_cells_draw_the_parent_study_seeds(path: Path, parent: Path, mode: str, effect: float) -> None:
    cell = _primary(path, mode, effect)
    parent_cell = _primary(parent, mode, effect)
    # The parents attach driver attribution to their orientation cells; it is a
    # post-hoc diagnostic that leaves the p-values alone, and the only difference.
    assert replace(parent_cell.evaluation_params, attribution=cell.evaluation_params.attribution) == (
        cell.evaluation_params
    )
    assert [derive_replicate_seed(cell, i) for i in range(cell.n_replicates)] == [
        derive_replicate_seed(parent_cell, i) for i in range(parent_cell.n_replicates)
    ]


def test_candidate_column_draws_the_baseline_seeds() -> None:
    grid = enumerate_study(load_study_config(PHASE6))
    for cell in (c for c in grid.cells if c.phase == DESIGN_PHASE):
        baseline = _primary(PHASE6, cell.metadata["trajectory_mode"], cell.metadata["effect_size"])
        assert [derive_replicate_seed(cell, i) for i in range(3)] == [
            derive_replicate_seed(baseline, i) for i in range(3)
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
def test_phase5_orientation_datasets_equal_the_paper_grade_datasets() -> None:
    """Six n = 1200 datasets: seeds equal is not enough at Phase 5, because the
    two configs differ in ``magnitude_kind``, which must be inert under orientation."""

    cell = _primary(PHASE5, "orientation", 1.0)
    parent_cell = _primary(PHASE5_PARENT, "orientation", 1.0)
    assert cell.generator_params.magnitude_kind != parent_cell.generator_params.magnitude_kind
    for index in range(3):
        assert _dataset_digest(cell, index) == _dataset_digest(parent_cell, index)
