"""Every post-floor study config keeps its cell ids, signatures and matched seeds.

``tests/data/post_floor_config_identity.json`` pins the configs committed after
the rank floor (the design-guidance pilot profiles first); they name their
estimator and floor explicitly and are not in the historical pre-flip fixture.
A diff here means a committed study would no longer resume, or that a config
was edited without regenerating the fixture (``tests/post_floor_fixture_capture.py``).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from motco.simulations.grid import derive_replicate_seed, parameter_signature
from motco.simulations.study import enumerate_study, load_study_config

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
FIXTURE = Path(__file__).resolve().parent / "data" / "post_floor_config_identity.json"
POST_FLOOR = json.loads(FIXTURE.read_text(encoding="utf-8"))["configs"]


def test_fixture_lists_the_capture_script_configs() -> None:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    try:
        from post_floor_fixture_capture import POST_FLOOR_CONFIGS
    finally:
        sys.path.pop(0)

    assert sorted(POST_FLOOR) == sorted(POST_FLOOR_CONFIGS)
    assert all(POST_FLOOR.values())


@pytest.mark.parametrize("name", sorted(POST_FLOOR))
def test_post_floor_config_identity_is_unchanged(name: str) -> None:
    grid = enumerate_study(load_study_config(CONFIG_DIR / name))
    actual = {
        cell.cell_id: {
            "signature": parameter_signature(cell),
            "seeds": [derive_replicate_seed(cell, i) for i in range(min(3, cell.n_replicates))],
        }
        for cell in grid.cells
    }
    expected = POST_FLOOR[name]
    assert set(actual) == set(expected)
    assert actual == expected
