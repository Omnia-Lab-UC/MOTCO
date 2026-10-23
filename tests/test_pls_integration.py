from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from motco.simulations import (
    SemiSyntheticTrajectoryDataset,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
)

# Small CV knobs keep ``plsda_doubleCV`` cheap for the fast test suite while
# still exercising the full double-nested selection path.
SMALL_CV = {"cv2_splits": 3, "cv1_splits": 2, "n_repeats": 2, "max_components": 5}

_LAYERS = ("methylation", "expression", "proteomics")


def make_pls_dataset(
    *, n_stages: int = 2, n_per_cell: int = 8, n_features: int = 12, seed: int = 0
) -> SemiSyntheticTrajectoryDataset:
    """A dataset with a clear stage signal so PLS-DA is well-posed and stable."""

    rng = np.random.default_rng(seed)
    rows: list[tuple[str, str, int]] = []
    blocks: dict[str, list[np.ndarray]] = {layer: [] for layer in _LAYERS}
    for stage in range(n_stages):
        for group in ("A", "B"):
            for replicate in range(n_per_cell):
                rows.append((f"s{stage}_{group}_{replicate}", group, stage))
                for layer_index, layer in enumerate(_LAYERS):
                    signal = float(stage) + 0.3 * layer_index
                    blocks[layer].append(signal + rng.normal(scale=0.5, size=n_features))
    sample_ids = [row[0] for row in rows]
    metadata = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "group": [row[1] for row in rows],
            "stage": [row[2] for row in rows],
            "cluster": [row[2] for row in rows],
        },
        index=sample_ids,
    )
    frames = {
        layer: pd.DataFrame(
            np.asarray(mats),
            index=sample_ids,
            columns=[f"{layer[0]}{j}" for j in range(n_features)],
        )
        for layer, mats in blocks.items()
    }
    return SemiSyntheticTrajectoryDataset(
        methylation=frames["methylation"],
        expression=frames["expression"],
        proteomics=frames["proteomics"],
        metadata=metadata,
        truth={"trajectory_mode": "magnitude", "seed": seed},
    )


def test_pls_integration_builds_latent_space_and_records_metadata() -> None:
    result = integrate_semisynthetic_dataset(
        make_pls_dataset(),
        SimulationEvaluationParams(integration_method="pls", integration_params=SMALL_CV),
    )

    assert result.metadata["integration_method"] == "pls"
    assert result.metadata["integration_role"] == "latent_space"
    selected_lv = result.metadata["integration_params"]["selected_lv"]
    assert selected_lv >= 2  # the harness's default rank floor
    assert result.metadata["integration_params"]["min_components"] == 2
    assert result.matrix.shape == (32, selected_lv)
    assert list(result.matrix.columns) == [f"pls_{i}" for i in range(selected_lv)]
    assert result.metadata["integration_params"]["stage_col"] == "stage"
    assert "cv_mean_auroc" in result.metadata


def test_pls_integration_is_deterministic() -> None:
    params = SimulationEvaluationParams(integration_method="pls", integration_params=SMALL_CV)
    first = integrate_semisynthetic_dataset(make_pls_dataset(), params)
    second = integrate_semisynthetic_dataset(make_pls_dataset(), params)

    assert (
        first.metadata["integration_params"]["selected_lv"]
        == second.metadata["integration_params"]["selected_lv"]
    )
    np.testing.assert_array_equal(first.matrix.to_numpy(), second.matrix.to_numpy())


def test_pls_integration_handles_multiple_stages() -> None:
    result = integrate_semisynthetic_dataset(
        make_pls_dataset(n_stages=3, n_per_cell=8),
        SimulationEvaluationParams(integration_method="pls", integration_params=SMALL_CV),
    )

    assert result.matrix.shape[0] == 48
    assert result.metadata["integration_method"] == "pls"


def test_pls_integration_rejects_too_few_samples_per_stage() -> None:
    # 1 sample per (stage, group) -> 2 per stage cannot support the CV folds.
    dataset = make_pls_dataset(n_stages=2, n_per_cell=1)
    params = SimulationEvaluationParams(integration_method="pls", integration_params=SMALL_CV)

    with pytest.raises(SimulationEvaluationError, match="PLS integration"):
        integrate_semisynthetic_dataset(dataset, params)


def test_pls_integration_end_to_end_trajectory_statistics() -> None:
    result = evaluate_semisynthetic_trajectory(
        make_pls_dataset(n_stages=3, n_per_cell=8),
        SimulationEvaluationParams(
            integration_method="pls",
            integration_params=SMALL_CV,
            permutations=2,
            seed=7,
        ),
    )

    assert "pls_latent" in result.realized_geometry["checkpoints"]
    latent_diagnostic = result.realized_geometry["checkpoints"]["pls_latent"]["joint"]
    assert latent_diagnostic["delta"] == pytest.approx(result.pair_statistics["delta"])
    assert latent_diagnostic["angle"] == pytest.approx(result.pair_statistics["angle"])

    assert set(result.pair_statistics) == {"delta", "angle", "shape"}
    assert result.latent_matrix_metadata["integration_method"] == "pls"
    assert result.latent_matrix_metadata["integration_role"] == "latent_space"
