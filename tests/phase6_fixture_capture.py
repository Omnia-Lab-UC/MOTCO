"""Builds the pre-change Phase 6 fixture compared by ``tests/test_phase6_byte_identity.py``.

Run from the repo root on the revision named in the fixture's
``capture_revision`` to regenerate ``tests/data/phase6_prechange_fixture.json``:

    uv run python tests/phase6_fixture_capture.py > tests/data/phase6_prechange_fixture.json

The fixture records PC1-orientation evaluations: every entry in ``EVALUATIONS``
pins ``orientation_estimator="pc1"``, so the net-displacement harness default
does not change what the byte-identity test compares.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from typing import Any

import numpy as np

from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationParams,
    generate_semisynthetic_trajectory,
)
from motco.simulations.evaluation import AttributionDiagnosticSettings, evaluate_semisynthetic_trajectory

DATASET_PARAMS = {
    "magnitude": SemiSyntheticTrajectoryParams(
        seed=11, trajectory_mode="magnitude", n_samples=90, n_stages=3, group_effect_size=0.5, group_ratio=0.4
    ),
    "orientation": SemiSyntheticTrajectoryParams(
        seed=12, trajectory_mode="orientation", n_samples=80, n_stages=3, group_effect_size=0.3, p_dmp=0.1
    ),
}

EVALUATIONS = {
    "concat": SimulationEvaluationParams(
        integration_method="concat", permutations=9, seed=5, orientation_estimator="pc1"
    ),
    "snf": SimulationEvaluationParams(integration_method="snf", permutations=9, seed=5, orientation_estimator="pc1"),
    "pls": SimulationEvaluationParams(
        integration_method="pls",
        integration_params={"n_repeats": 2},
        permutations=9,
        seed=5,
        attribution=AttributionDiagnosticSettings(enabled=True, bootstrap_replicates=3, top_k=5),
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


def indicator_digests(indicators: dict[str, dict[str, np.ndarray]]) -> dict[str, dict[str, str]]:
    """Exact integer-valued indicators, digested to keep the fixture small."""

    return {
        group: {
            layer: hashlib.sha256(np.ascontiguousarray(values, dtype=np.int64).tobytes()).hexdigest()
            for layer, values in layers.items()
        }
        for group, layers in indicators.items()
    }


def dataset_summary(params: SemiSyntheticTrajectoryParams) -> dict[str, Any]:
    dataset = generate_semisynthetic_trajectory(params)
    truth = dict(dataset.truth)
    truth["indicators"] = indicator_digests(truth["indicators"])
    return {
        "metadata": jsonable(dataset.metadata.to_dict(orient="list")),
        "layers": {
            layer: {
                "columns": list(getattr(dataset, layer).columns.astype(str)),
                "column_sums": jsonable(getattr(dataset, layer).to_numpy().sum(axis=0)),
                "first_row": jsonable(getattr(dataset, layer).to_numpy()[0]),
            }
            for layer in ("methylation", "expression", "proteomics")
        },
        "truth": jsonable(truth),
    }


def evaluation_summary(name: str) -> dict[str, Any]:
    dataset = generate_semisynthetic_trajectory(DATASET_PARAMS["orientation"])
    result = evaluate_semisynthetic_trajectory(dataset, EVALUATIONS[name])
    attribution = dict(result.attribution_diagnostics)
    attribution.pop("runtime", None)
    return jsonable(
        {
            "pair_statistics": result.pair_statistics,
            "p_values": result.p_values,
            "null_summary": result.null_summary,
            "config_spectrum": result.config_spectrum,
            "latent_matrix_metadata": result.latent_matrix_metadata,
            "realized_geometry": result.realized_geometry,
            "attribution_diagnostics": attribution,
        }
    )


def main() -> None:
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    payload = {
        "note": "Generator and evaluation outputs captured before group_stage_sizes and "
        "integration_params.layers existed (phase6-small-n-operating-study).",
        "capture_revision": revision,
        "datasets": {name: dataset_summary(params) for name, params in DATASET_PARAMS.items()},
        "evaluations": {name: evaluation_summary(name) for name in EVALUATIONS},
    }
    print(json.dumps(payload, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
