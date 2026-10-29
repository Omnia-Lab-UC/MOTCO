"""Builds the post-floor config identity fixture compared by ``tests/test_post_floor_config_identity.py``.

Run from the repo root to regenerate ``tests/data/post_floor_config_identity.json``:

    uv run python tests/post_floor_fixture_capture.py > tests/data/post_floor_config_identity.json

The fixture pins every study config committed **after** the rank floor
(``minimum-latent-rank``, 0.10.0): the cell ids, parameter signatures and the
first three matched seeds of every enumerated cell, exactly as
``tests/data/orientation_default_identity.json`` does for the 19 historical
configs. Those historical configs record what was run under PC1 and the
unfloored rank rule; the configs listed here name their estimator and floor
explicitly and are not held to the historical pins. A new config committed under
``examples/trajectory_power_study/`` is appended to ``POST_FLOOR_CONFIGS`` and the
fixture is regenerated; the pin tests require every committed config to be in
exactly one of the two fixtures.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any

from motco.simulations.grid import derive_replicate_seed, parameter_signature
from motco.simulations.study import enumerate_study, load_study_config

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"

#: Configs committed after the rank floor, in the order they were added.
POST_FLOOR_CONFIGS: tuple[str, ...] = (
    "design_guidance_core_magnitude_pilot.json",
    "design_guidance_core_geometry_pilot.json",
    "design_guidance_arms_magnitude_pilot.json",
    "design_guidance_arms_geometry_pilot.json",
    # design-guidance-paper-grade-slice (0.12.0): the paper-grade slice and the stage-effect arm.
    "design_guidance_core_magnitude_study.json",
    "design_guidance_core_geometry_study.json",
    "design_guidance_stage_effect_magnitude_study.json",
    "design_guidance_stage_effect_geometry_study.json",
)


def config_identity(name: str) -> dict[str, dict[str, Any]]:
    """Cell id -> ``{signature, seeds (replicates 0-2)}`` for one committed config."""

    grid = enumerate_study(load_study_config(CONFIG_DIR / name))
    return {
        cell.cell_id: {
            "signature": parameter_signature(cell),
            "seeds": [derive_replicate_seed(cell, i) for i in range(min(3, cell.n_replicates))],
        }
        for cell in grid.cells
    }


def main() -> None:
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    payload = {
        "note": "Identity of every study config committed after the rank floor (minimum-latent-rank, 0.10.0), "
        "starting with the four design-guidance pilot profiles (design-guidance-grid), then the paper-grade "
        "slice and stage-effect profiles (design-guidance-paper-grade-slice). Maps each config "
        "under examples/trajectory_power_study/ to {cell_id: {signature, seeds (replicates 0-2)}}. These "
        "configs name evaluation.orientation_estimator and evaluation.min_components explicitly and are not "
        "held to the historical PC1 / floor-1 pins. Checked by tests/test_post_floor_config_identity.py; "
        "regenerate with tests/post_floor_fixture_capture.py after appending a new config.",
        "capture_revision": revision,
        "configs": {name: config_identity(name) for name in POST_FLOOR_CONFIGS},
    }
    print(json.dumps(payload, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
