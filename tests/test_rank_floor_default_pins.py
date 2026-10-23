"""Committed study configs keep their identity across the rank-floor default flip.

``tests/data/orientation_default_identity.json`` was captured on 2026-10-15
with the configs unedited, before ``evaluation.min_components`` existed, so it
records what was actually committed and run: every cell of all 19 configs under
the unfloored rank rule (``tests/data/committed_config_identity.json`` predates
the net-displacement profiles and covers 16; ``tests/test_committed_config_identity.py``
keeps checking those). Each config now pins ``"min_components": 1`` beside
``"orientation_estimator": "pc1"``; these tests prove the pins are
behavior-neutral (the loaded field is the hash-omitted ``None``, and every cell
id, signature and matched seed is unchanged) and that every committed config
names the floor deliberately.
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
FLOOR_AXIS = "evaluation.min_components"


def test_fixture_covers_every_committed_config() -> None:
    assert sorted(CAPTURED) == sorted(path.name for path in CONFIG_DIR.glob("*.json"))
    assert len(CAPTURED) == 19
    assert all(CAPTURED.values())


@pytest.mark.parametrize("name", sorted(CAPTURED))
def test_committed_config_identity_survives_the_floor_default_flip(name: str) -> None:
    config = load_study_config(CONFIG_DIR / name)
    assert config.evaluation.min_components is None  # the hash-omitted spelling of the unfloored rule
    assert config.evaluation.resolved_min_components == 1

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
    assert all(cell.evaluation_params.min_components is None for cell in grid.cells)


@pytest.mark.parametrize("path", sorted(CONFIG_DIR.glob("*.json")), ids=lambda p: p.name)
def test_every_committed_config_names_its_rank_floor(path: Path) -> None:
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert "min_components" in raw["evaluation"], (
        f"{path.name} must choose evaluation.min_components deliberately; "
        "do not copy the historical floor-1 pin without reason."
    )
    assert raw["evaluation"]["min_components"] == 1, (
        f"{path.name}: every committed config predates the rank floor and pins 1 (historical)."
    )

    design_grid = raw.get("design_grid") or {}
    points = list(design_grid.get("points") or [])
    axes = design_grid.get("axes") or {}
    for point in points:
        if FLOOR_AXIS in point:
            assert point[FLOOR_AXIS] is not None, (
                f"{path.name}: a design point writes null for {FLOOR_AXIS}; "
                "null names the evaluator's default, so name the unfloored column 1."
            )
    if FLOOR_AXIS in axes:
        assert None not in axes[FLOOR_AXIS], f"{path.name}: {FLOOR_AXIS} axis lists null"
