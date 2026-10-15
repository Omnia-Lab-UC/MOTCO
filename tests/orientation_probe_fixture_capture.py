"""Builds the pre-change fixture compared by ``tests/test_orientation_probe_byte_identity.py``.

Run from the repo root on the revision named in the fixture's
``capture_revision`` to regenerate ``tests/data/orientation_probe_prechange_fixture.json``:

    uv run python tests/orientation_probe_fixture_capture.py > tests/data/orientation_probe_prechange_fixture.json

The fixture records PC1-orientation evaluations: every entry in ``EVALUATIONS``
pins ``orientation_estimator="pc1"``, so the net-displacement harness default
does not change what the byte-identity test compares.
"""

from __future__ import annotations

import json
import subprocess
from typing import Any

import numpy as np

from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationParams,
    generate_semisynthetic_trajectory,
)
from motco.simulations.evaluation import evaluate_semisynthetic_trajectory

DATASET_PARAMS = SemiSyntheticTrajectoryParams(
    seed=21, trajectory_mode="orientation", n_samples=80, n_stages=3, group_effect_size=0.5, p_dmp=0.1
)

EVALUATIONS = {
    "concat": SimulationEvaluationParams(
        integration_method="concat",
        permutations=19,
        seed=7,
        include_null_distributions=True,
        orientation_estimator="pc1",
    ),
    "snf": SimulationEvaluationParams(
        integration_method="snf",
        permutations=19,
        seed=7,
        include_null_distributions=True,
        orientation_estimator="pc1",
    ),
    "pls": SimulationEvaluationParams(
        integration_method="pls",
        integration_params={"n_repeats": 2},
        permutations=19,
        seed=7,
        include_null_distributions=True,
        orientation_estimator="pc1",
    ),
}


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [jsonable(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def evaluation_summary(params: SimulationEvaluationParams) -> dict[str, Any]:
    """Every non-runtime field of one evaluation result, JSON-safe."""

    dataset = generate_semisynthetic_trajectory(DATASET_PARAMS)
    result = evaluate_semisynthetic_trajectory(dataset, params)
    runtime = {k: v for k, v in result.runtime_metadata.items() if k != "runtime_seconds"}
    return jsonable(
        {
            "observed_angles": result.observed_angles,
            "observed_deltas": result.observed_deltas,
            "observed_shapes": result.observed_shapes,
            "pair_statistics": result.pair_statistics,
            "p_values": result.p_values,
            "null_distributions": result.null_distributions,
            "null_summary": result.null_summary,
            "config_spectrum": result.config_spectrum,
            "latent_matrix_metadata": result.latent_matrix_metadata,
            "realized_geometry": result.realized_geometry,
            "runtime_metadata": runtime,
        }
    )


def main() -> None:
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    payload = {
        "note": "Evaluation outputs captured before the net-displacement orientation estimator and "
        "the diagnostic pca / pls_within_group integration methods existed "
        "(orientation-measurement-probe).",
        "capture_revision": revision,
        "evaluations": {name: evaluation_summary(params) for name, params in EVALUATIONS.items()},
    }
    print(json.dumps(payload, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
