"""Within-stage permutation draws in the simulation harness (``permutation_strata``)."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

import motco.simulations.evaluation as evaluation_module
from motco.simulations.evaluation import (
    PERMUTATION_STRATA,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    _empirical_p_value,
    _extract_null_distributions,
    _summarize_null_distributions,
    build_simulation_trajectory_design,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
    pls_coordinates_from_pls,
    stage_strata,
)
from motco.simulations.grid import _to_jsonable, make_simulation_cell, parameter_signature
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import generate_semisynthetic_trajectory
from motco.stats.permutation import RRPP
from tests.orthonormal_coordinates_fixture_capture import DATASET, FAST_CV

PERMUTATIONS = 29
SEED = 17


@pytest.fixture(scope="module")
def dataset():
    return generate_semisynthetic_trajectory(DATASET, reference=load_reference())


def _params(**overrides) -> SimulationEvaluationParams:
    base = {
        "integration_method": "pls",
        "integration_params": FAST_CV,
        "permutations": PERMUTATIONS,
        "seed": SEED,
        "orientation_estimator": "net_displacement",
    }
    return SimulationEvaluationParams(**{**base, **overrides})


@pytest.fixture(scope="module")
def scores(dataset):
    return integrate_semisynthetic_dataset(dataset, _params())


def _primitive(dataset, latent, params: SimulationEvaluationParams, *, strata: bool):
    """The p-values and null summary of RRPP called directly on the evaluation's latent matrix and models."""

    design = build_simulation_trajectory_design(dataset.metadata, group_col="group", stage_col="stage")
    deltas, angles, shapes, _ = RRPP(
        latent.matrix, design.model_full, design.model_reduced, design.ls_means, design.contrast,
        permutations=PERMUTATIONS, progress=False, seed=SEED, return_eigengaps=True,
        orientation_estimator=params.resolved_orientation_estimator,
        strata=stage_strata(dataset.metadata, "stage", design.stage_levels) if strata else None,
    )
    return _extract_null_distributions(deltas, angles, shapes, shape_available=len(design.stage_levels) >= 3)


def _check_against_primitive(result, null) -> None:
    for statistic, observed in result.pair_statistics.items():
        if statistic in null and np.isfinite(observed):
            assert result.p_values[statistic] == _empirical_p_value(null[statistic], observed)
    assert result.null_summary == _summarize_null_distributions(null)


def test_strata_are_named() -> None:
    assert PERMUTATION_STRATA == ("stage",)


def test_unknown_strata_are_rejected_before_any_fit() -> None:
    with pytest.raises(SimulationEvaluationError, match=r"Unknown permutation_strata 'group'; allowed: \['stage'\]"):
        SimulationEvaluationParams(permutation_strata="group")


def test_stage_strata_follow_the_design_order() -> None:
    import pandas as pd

    metadata = pd.DataFrame({"stage": ["10", "2", "1", "10", "2"]})
    np.testing.assert_array_equal(stage_strata(metadata, "stage", ["1", "2", "10"]), [2, 1, 0, 2, 1])
    with pytest.raises(SimulationEvaluationError, match="not among the design's stage levels"):
        stage_strata(metadata, "stage", ["1", "2"])


def test_only_the_null_moves(dataset, scores) -> None:
    plain = evaluate_semisynthetic_trajectory(dataset, _params(), latent=scores)
    within = evaluate_semisynthetic_trajectory(dataset, _params(permutation_strata="stage"), latent=scores)
    assert within.pair_statistics == plain.pair_statistics
    for ours, theirs in (
        (within.observed_deltas, plain.observed_deltas),
        (within.observed_angles, plain.observed_angles),
        (within.observed_shapes, plain.observed_shapes),
    ):
        assert np.asarray(ours).tobytes() == np.asarray(theirs).tobytes()
    assert within.latent_matrix_metadata == plain.latent_matrix_metadata
    assert within.null_summary != plain.null_summary
    _check_against_primitive(within, _primitive(dataset, scores, _params(), strata=True))
    _check_against_primitive(plain, _primitive(dataset, scores, _params(), strata=False))


def test_the_setting_does_not_change_the_fit(dataset) -> None:
    plain = evaluate_semisynthetic_trajectory(dataset, _params(permutations=0))
    within = evaluate_semisynthetic_trajectory(dataset, _params(permutations=0, permutation_strata="stage"))
    assert within.latent_matrix_metadata == plain.latent_matrix_metadata
    assert within.pair_statistics == plain.pair_statistics
    assert within.p_values == plain.p_values == {}


def test_a_prebuilt_integration_is_not_refitted(dataset, scores, monkeypatch) -> None:
    calls: list[str] = []
    real = evaluation_module.integrate_semisynthetic_dataset

    def counting(*args, **kwargs):
        calls.append("fit")
        return real(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "integrate_semisynthetic_dataset", counting)
    evaluate_semisynthetic_trajectory(dataset, _params(permutation_strata="stage"), latent=scores)
    assert calls == []


def test_predictive_coordinates_under_within_stage_draws(dataset, scores) -> None:
    params = _params(latent_coordinates="predictive_orthonormal", permutation_strata="stage")
    latent = pls_coordinates_from_pls(dataset, scores, params)
    result = evaluate_semisynthetic_trajectory(dataset, params, latent=latent)
    assert latent.matrix.shape[1] == 2
    assert result.latent_matrix_metadata["selected_lv"] == scores.metadata["selected_lv"]
    _check_against_primitive(result, _primitive(dataset, latent, params, strata=True))


def test_concat_under_within_stage_draws(dataset) -> None:
    params = _params(integration_method="concat", integration_params={"layers": FAST_CV["layers"]},
                     permutation_strata="stage")
    latent = integrate_semisynthetic_dataset(dataset, params)
    result = evaluate_semisynthetic_trajectory(dataset, params, latent=latent)
    _check_against_primitive(result, _primitive(dataset, latent, params, strata=True))


def test_runtime_metadata_names_the_strata_only_when_set(dataset, scores) -> None:
    within = evaluate_semisynthetic_trajectory(dataset, _params(permutation_strata="stage"), latent=scores)
    plain = evaluate_semisynthetic_trajectory(dataset, _params(), latent=scores)
    assert within.runtime_metadata["permutation_strata"] == "stage"
    assert "permutation_strata" not in plain.runtime_metadata


def test_no_permutations_with_strata_completes(dataset, scores) -> None:
    result = evaluate_semisynthetic_trajectory(
        dataset, _params(permutations=0, permutation_strata="stage"), latent=scores
    )
    assert result.p_values == {}
    assert result.null_summary == {}


def test_cell_identity_moves_only_when_set() -> None:
    plain = make_simulation_cell(phase="power_primary", generator_params=DATASET)
    within = make_simulation_cell(
        phase="power_primary",
        generator_params=DATASET,
        evaluation_params=SimulationEvaluationParams(permutation_strata="stage"),
    )
    assert within.cell_id != plain.cell_id
    assert parameter_signature(within) != parameter_signature(plain)
    assert _to_jsonable(within.evaluation_params)["permutation_strata"] == "stage"
    assert "permutation_strata" not in _to_jsonable(SimulationEvaluationParams())
    explicit_none = make_simulation_cell(
        phase="power_primary",
        generator_params=DATASET,
        evaluation_params=replace(SimulationEvaluationParams(), permutation_strata=None),
    )
    assert explicit_none.cell_id == plain.cell_id
    assert parameter_signature(explicit_none) == parameter_signature(plain)
