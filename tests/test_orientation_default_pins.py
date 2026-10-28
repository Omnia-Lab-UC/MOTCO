"""Committed study configs keep their identity across the orientation-estimator default flip.

The fixture was captured on the pre-flip code with the configs unedited, so it
records what was actually committed and run, all under the PC1 estimator. Each
historical config pins ``"orientation_estimator": "pc1"`` (and the
net-displacement profiles name their baseline column ``"pc1"`` rather than
``null``); these tests prove the pins are behavior-neutral.

Configs committed after the flip are **not** historical: they name their
estimator explicitly (never ``null``, which would follow a later default) and
are pinned by ``tests/data/post_floor_config_identity.json`` instead
(``tests/test_post_floor_config_identity.py``). The two fixtures are disjoint
and together cover every committed config, so a config can be in neither.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from motco.simulations.grid import derive_replicate_seed, parameter_signature
from motco.simulations.study import enumerate_study, load_study_config

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
FIXTURE = Path(__file__).resolve().parent / "data" / "orientation_default_identity.json"
POST_FLOOR_FIXTURE = Path(__file__).resolve().parent / "data" / "post_floor_config_identity.json"
HISTORICAL = json.loads(FIXTURE.read_text(encoding="utf-8"))["configs"]
POST_FLOOR = json.loads(POST_FLOOR_FIXTURE.read_text(encoding="utf-8"))["configs"]
ESTIMATOR_AXIS = "evaluation.orientation_estimator"


def test_fixture_covers_every_committed_config() -> None:
    committed = sorted(path.name for path in CONFIG_DIR.glob("*.json"))
    assert not set(HISTORICAL) & set(POST_FLOOR), "a config is pinned as both historical and post-flip"
    assert sorted(set(HISTORICAL) | set(POST_FLOOR)) == committed, (
        "every committed config must be in exactly one identity fixture: the historical "
        f"{FIXTURE.name} or the post-floor {POST_FLOOR_FIXTURE.name} "
        "(append it to tests/post_floor_fixture_capture.py and regenerate)"
    )
    assert len(HISTORICAL) == 19  # the pre-flip set is closed
    assert all(HISTORICAL.values())


@pytest.mark.parametrize("name", sorted(HISTORICAL))
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
    expected = HISTORICAL[name]
    assert set(actual) == set(expected)
    assert actual == expected


@pytest.mark.parametrize("path", sorted(CONFIG_DIR.glob("*.json")), ids=lambda p: p.name)
def test_every_committed_config_names_its_orientation_estimator(path: Path) -> None:
    raw = json.loads(path.read_text(encoding="utf-8"))
    assert "orientation_estimator" in raw["evaluation"], (
        f"{path.name} must choose evaluation.orientation_estimator deliberately; "
        "do not copy the historical 'pc1' pin without reason."
    )
    if path.name in HISTORICAL:
        assert raw["evaluation"]["orientation_estimator"] == "pc1", (
            f"{path.name}: a historical config was run under PC1 and pins 'pc1' (historical)."
        )
    else:
        assert raw["evaluation"]["orientation_estimator"] is not None, (
            f"{path.name}: a post-flip config names its estimator; "
            "null would name the evaluator's default and silently follow a later flip."
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
