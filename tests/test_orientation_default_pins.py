"""Committed study configs keep their identity across the orientation-estimator default flip.

The fixture was captured on the pre-flip code with the configs unedited, so it
records what was actually committed and run, all under the PC1 estimator. Each
config now pins ``"orientation_estimator": "pc1"`` (and the net-displacement
profiles name their baseline column ``"pc1"`` rather than ``null``); these tests
prove the pins are behavior-neutral.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from motco.simulations.grid import derive_replicate_seed, parameter_signature
from motco.simulations.study import enumerate_study, load_study_config

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
FIXTURE = Path(__file__).resolve().parent / "data" / "orientation_default_identity.json"
CAPTURED = json.loads(FIXTURE.read_text(encoding="utf-8"))["configs"]
ESTIMATOR_AXIS = "evaluation.orientation_estimator"


def test_fixture_covers_every_committed_config() -> None:
    assert sorted(CAPTURED) == sorted(path.name for path in CONFIG_DIR.glob("*.json"))
    assert len(CAPTURED) == 19
    assert all(CAPTURED.values())


@pytest.mark.parametrize("name", sorted(CAPTURED))
def test_committed_config_identity_survives_the_default_flip(name: str) -> None:
    config = load_study_config(CONFIG_DIR / name)
    assert config.evaluation.orientation_estimator is None

    grid = enumerate_study(config)
    actual = {
        cell.cell_id: {
            "signature": parameter_signature(cell),
            "seeds": [derive_replicate_seed(cell, i) for i in range(min(3, cell.n_replicates))],
        }
        for cell in grid.cells
    }
    expected = CAPTURED[name]
    assert set(actual) == set(expected)
    assert actual == expected


@pytest.mark.parametrize("path", sorted(CONFIG_DIR.glob("*.json")), ids=lambda p: p.name)
def test_every_committed_config_names_its_orientation_estimator(path: Path) -> None:
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert "orientation_estimator" in raw["evaluation"], (
        f"{path.name} must choose evaluation.orientation_estimator deliberately; "
        "do not copy the historical 'pc1' pin without reason."
    )

    design_grid = raw.get("design_grid") or {}
    points = list(design_grid.get("points") or [])
    axes = design_grid.get("axes") or {}
    for point in points:
        if ESTIMATOR_AXIS in point:
            assert point[ESTIMATOR_AXIS] is not None, (
                f"{path.name}: a design point writes null for {ESTIMATOR_AXIS}; "
                "null now names the evaluator's default, so name the PC1 column 'pc1'."
            )
    if ESTIMATOR_AXIS in axes:
        assert None not in axes[ESTIMATOR_AXIS], f"{path.name}: {ESTIMATOR_AXIS} axis lists null"
