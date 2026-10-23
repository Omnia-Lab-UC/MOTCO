"""The PLS rank floor in the simulation harness: spelling, hashing, records.

``SimulationEvaluationParams.min_components`` defaults to 2 and is hashed into
cell ids, signatures and ``runtime_metadata``; ``1`` normalizes to ``None``, the
hash-omitted spelling of the unfloored rule, so floor-1 cells and records equal
those produced before the field existed (the orientation-estimator mechanism).
The core-level floor tests live in ``tests/test_latent_adjustment.py``.
"""

from __future__ import annotations

from dataclasses import asdict

import pytest

from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
)
from motco.simulations.grid import SimulationCell, parameter_signature, run_simulation_replicate
from motco.simulations.study import MatchedSeedPolicy, StudyConfig, enumerate_study
from tests.test_pls_integration import SMALL_CV, make_pls_dataset

_IDENTITY_FIELDS = {"cell_id", "cell_metadata", "parameter_signature", "runtime_metadata"}


def _cell(evaluation: SimulationEvaluationParams) -> SimulationCell:
    config = StudyConfig(
        generator=SemiSyntheticTrajectoryParams(
            seed=2, trajectory_mode="magnitude", n_samples=60, n_stages=3, p_dmp=0.1
        ),
        evaluation=evaluation,
        trajectory_modes=("translation",),
        effect_sizes=(0.0, 1.0),
        n_replicates=2,
        base_seed=11,
        matched_seeds=MatchedSeedPolicy(enabled=True, primary_family="floor"),
    )
    for cell in enumerate_study(config).cells:
        if cell.metadata.get("trajectory_mode") == "none" and cell.phase == "power_primary":
            return cell
    raise AssertionError("no anchor cell")


def _comparable(record) -> dict:
    payload = asdict(record)
    runtime = {k: v for k, v in payload["runtime_metadata"].items() if k != "runtime_seconds"}
    return {**{k: v for k, v in payload.items() if k not in _IDENTITY_FIELDS}, "runtime": runtime}


# ── Spelling ────────────────────────────────────────────────────────────────


def test_default_floor_is_two_and_one_normalizes_to_none() -> None:
    assert SimulationEvaluationParams().min_components == 2
    assert SimulationEvaluationParams().resolved_min_components == 2
    one = SimulationEvaluationParams(min_components=1)
    assert one.min_components is None
    assert one.resolved_min_components == 1
    assert SimulationEvaluationParams(min_components=None).resolved_min_components == 1
    assert SimulationEvaluationParams(min_components=3).min_components == 3


@pytest.mark.parametrize("floor", [0, 2.5, -1, "2", True])
def test_invalid_floor_is_rejected(floor) -> None:
    with pytest.raises(SimulationEvaluationError, match="min_components must be an integer of at least 1"):
        SimulationEvaluationParams(min_components=floor)


# ── Hashing ─────────────────────────────────────────────────────────────────


def test_floor_one_is_hash_omitted_and_the_default_is_hashed() -> None:
    common = {"integration_method": "pls", "integration_params": dict(SMALL_CV), "permutations": 0, "seed": 3}
    plain = SimulationEvaluationParams(**common)
    one = SimulationEvaluationParams(**common, min_components=1)
    two = SimulationEvaluationParams(**common, min_components=2)
    none_cell, one_cell, two_cell = _cell(plain), _cell(one), _cell(two)
    # The default is 2, so ``plain`` and ``two`` are the same cell.
    assert none_cell.cell_id == two_cell.cell_id
    assert parameter_signature(none_cell) == parameter_signature(two_cell)
    # Floor 1 hashes as the field never existed: neither the id nor the signature mentions it.
    assert one_cell.cell_id != two_cell.cell_id
    assert "min_components" not in parameter_signature(one_cell)
    assert one_cell.evaluation_params.min_components is None


def test_floor_one_record_equals_a_pre_floor_record() -> None:
    """A floor-1 PLS evaluation records no floor entry; the default records 2 and selects rank >= 2."""

    dataset = make_pls_dataset(n_stages=3, n_per_cell=6, n_features=8, seed=4)
    one = evaluate_semisynthetic_trajectory(
        dataset,
        SimulationEvaluationParams(integration_method="pls", integration_params=SMALL_CV, min_components=1),
    )
    two = evaluate_semisynthetic_trajectory(
        dataset, SimulationEvaluationParams(integration_method="pls", integration_params=SMALL_CV)
    )
    assert "min_components" not in one.runtime_metadata
    assert "min_components" not in one.latent_matrix_metadata["integration_params"]
    assert two.runtime_metadata["min_components"] == 2
    assert two.latent_matrix_metadata["integration_params"]["min_components"] == 2
    assert two.latent_matrix_metadata["selected_lv"] >= 2


def test_concat_hashes_the_floor_but_measures_nothing_differently() -> None:
    common = {"integration_method": "concat", "permutations": 9, "seed": 2026}
    one_cell = _cell(SimulationEvaluationParams(**common, min_components=1))
    two_cell = _cell(SimulationEvaluationParams(**common))
    assert one_cell.cell_id != two_cell.cell_id
    one = run_simulation_replicate(one_cell, 0)
    two = run_simulation_replicate(two_cell, 0)
    assert "min_components" not in one.runtime_metadata
    assert two.runtime_metadata["min_components"] == 2
    assert one.generator_seed == two.generator_seed
    assert one.pair_statistics == two.pair_statistics
    assert one.p_values == two.p_values
    one_payload, two_payload = _comparable(one), _comparable(two)
    one_payload["runtime"].pop("min_components", None)
    two_payload["runtime"].pop("min_components", None)
    assert one_payload == two_payload
