"""Capture the harness evaluations the orthonormal-coordinates change must leave unchanged.

Run before ``latent_coordinates`` existed (orthonormal-span-coordinates), from the
repository root::

    uv run python -m tests.orthonormal_coordinates_fixture_capture \
        > tests/data/orthonormal_coordinates_prechange_fixture.json

Every entry is evaluated under the current harness defaults (net displacement,
rank floor 2) on one stage-span dataset, with RRPP and span capture where the
space supports it. ``tests/test_orthonormal_coordinates.py`` compares an
evaluation without the setting against it, at the tight tolerance of
``tests/test_phase6_byte_identity.py``.
"""

from __future__ import annotations

import json
import subprocess
from typing import Any

from motco.simulations.evaluation import SimulationEvaluationParams, evaluate_semisynthetic_trajectory
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory
from tests.phase6_fixture_capture import jsonable

TWO_BLOCKS = ("methylation", "expression")
FAST_CV = {"layers": list(TWO_BLOCKS), "n_repeats": 1, "max_components": 6}

DATASET = SemiSyntheticTrajectoryParams(
    seed=31,
    n_samples=90,
    p_dmp=0.1,
    trajectory_mode="shape",
    group_effect_size=0.5,
    construction="stage_span",
    span_layers=TWO_BLOCKS,
    span_scale="planted",
)

EVALUATIONS: dict[str, SimulationEvaluationParams] = {
    "pls": SimulationEvaluationParams(
        integration_method="pls", integration_params=FAST_CV, permutations=19, seed=5, span_capture=True
    ),
    "pls_forced": SimulationEvaluationParams(
        integration_method="pls",
        integration_params={**FAST_CV, "forced_components": 3},
        permutations=19,
        seed=5,
        span_capture=True,
    ),
    "pls_residual_pca": SimulationEvaluationParams(
        integration_method="pls_residual_pca", integration_params=FAST_CV, permutations=19, seed=5, span_capture=True
    ),
    "pca": SimulationEvaluationParams(
        integration_method="pca",
        integration_params={"layers": list(TWO_BLOCKS), "n_components": 3},
        permutations=19,
        seed=5,
        span_capture=True,
    ),
    "pls_within_group": SimulationEvaluationParams(
        integration_method="pls_within_group",
        integration_params={"layers": list(TWO_BLOCKS)},
        permutations=19,
        seed=5,
    ),
    "concat": SimulationEvaluationParams(
        integration_method="concat",
        integration_params={"layers": list(TWO_BLOCKS)},
        permutations=19,
        seed=5,
        span_capture=True,
    ),
    "snf": SimulationEvaluationParams(
        integration_method="snf",
        integration_params={"layers": list(TWO_BLOCKS), "spectral_components": 4},
        permutations=19,
        seed=5,
    ),
}


def evaluation_summary(name: str) -> dict[str, Any]:
    dataset = generate_semisynthetic_trajectory(DATASET, reference=load_reference())
    result = evaluate_semisynthetic_trajectory(dataset, EVALUATIONS[name])
    runtime = dict(result.runtime_metadata)
    runtime.pop("runtime_seconds")
    return jsonable(
        {
            "pair_statistics": result.pair_statistics,
            "p_values": result.p_values,
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
        "note": "Harness evaluations captured before latent_coordinates existed (orthonormal-span-coordinates).",
        "capture_revision": revision,
        "evaluations": {name: evaluation_summary(name) for name in EVALUATIONS},
    }
    print(json.dumps(payload, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
