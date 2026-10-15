from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from motco.simulations import (
    SemiSyntheticTrajectoryDataset,
    SimulationEvaluationError,
    SimulationEvaluationParams,
    build_simulation_trajectory_design,
    evaluate_semisynthetic_trajectory,
    integrate_semisynthetic_dataset,
)


def make_dataset(*, n_stages: int = 3) -> SemiSyntheticTrajectoryDataset:
    rows = []
    for stage in range(n_stages):
        for group in ("A", "B"):
            for replicate in range(2):
                rows.append((f"s{stage}_{group}_{replicate}", group, stage, replicate))
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

    values = []
    for _, group, stage, replicate in rows:
        group_offset = 1.0 if group == "B" else 0.0
        values.append(
            [
                float(stage + group_offset + replicate * 0.1),
                float((stage + 1) * (group_offset + 1.0)),
                float(replicate - stage * 0.2),
            ]
        )
    base = pd.DataFrame(values, index=sample_ids, columns=["f0", "f1", "f2"])
    return SemiSyntheticTrajectoryDataset(
        methylation=base.add_prefix("m"),
        expression=(base + 0.5).add_prefix("g"),
        proteomics=(base * 0.5).add_prefix("p"),
        metadata=metadata,
        truth={"trajectory_mode": "magnitude", "seed": 17},
    )


# Captured from ``concat`` / 5 permutations / seed 2026 at commit fab94d6 — i.e.
# before the permutation-null summary existed. Pinned so the summary can be shown
# to leave the test it describes untouched.
PRE_CHANGE_P_VALUES = {
    "angle": 0.3333333333333333,
    "delta": 0.6666666666666666,
    "shape": 0.3333333333333333,
}
PRE_CHANGE_PAIR_STATISTICS = {
    "angle": 42.194748480162865,
    "delta": 0.23766293337564282,
    "shape": 0.2704919212869626,
}


def test_concat_integration_standardizes_and_records_metadata() -> None:
    result = integrate_semisynthetic_dataset(make_dataset(), SimulationEvaluationParams(integration_method="concat"))

    assert result.matrix.shape == (12, 9)
    assert result.metadata["integration_method"] == "concat"
    assert result.metadata["integration_params"] == {"standardize": True}
    np.testing.assert_allclose(result.matrix.mean(axis=0).to_numpy(), 0.0, atol=1e-12)


def test_snf_integration_uses_existing_helpers_and_records_resolved_params() -> None:
    result = integrate_semisynthetic_dataset(
        make_dataset(),
        SimulationEvaluationParams(
            integration_method="snf",
            integration_params={"K": 3, "k": 3, "t": 2, "spectral_components": 2},
        ),
    )

    assert result.matrix.shape == (12, 2)
    assert result.metadata["integration_method"] == "snf"
    assert result.metadata["integration_params"]["spectral_components"] == 2
    assert result.metadata["fused_shape"] == (12, 12)


def _pls_params(**integration_params: Any) -> SimulationEvaluationParams:
    # ``max_components`` is held below the fixture's inner-fold rank so the CV
    # path is exercisable on this deliberately small dataset.
    params: dict[str, Any] = {
        "stage_col": "stage",
        "n_repeats": 1,
        "max_components": 4,
        **integration_params,
    }
    return SimulationEvaluationParams(integration_method="pls", integration_params=params)


def test_pls_integration_selects_components_by_cross_validation_by_default() -> None:
    result = integrate_semisynthetic_dataset(make_dataset(), _pls_params())

    assert result.metadata["component_selection"] == "cv"
    assert result.metadata["cv_repeats_completed"] > 0
    assert result.metadata["cv_mean_auroc"] is not None
    assert "forced_components" not in result.metadata["integration_params"]
    assert result.matrix.shape[1] == result.metadata["selected_lv"]


def test_pls_integration_honors_and_records_a_forced_component_count() -> None:
    result = integrate_semisynthetic_dataset(make_dataset(), _pls_params(forced_components=5))

    assert result.matrix.shape[1] == 5
    assert result.metadata["selected_lv"] == 5
    assert result.metadata["component_selection"] == "forced"
    assert result.metadata["integration_params"]["forced_components"] == 5
    # No cross-validation ran, so none of its outputs are reported.
    assert result.metadata["cv_mean_auroc"] is None
    assert result.metadata["cv_repeats_completed"] == 0


@pytest.mark.parametrize("forced", [1, 10])
def test_pls_integration_rejects_an_infeasible_forced_component_count(forced: int) -> None:
    """Rejected, never clamped: a diagnostic that moved its own rank is useless."""

    with pytest.raises(SimulationEvaluationError, match="outside the feasible range"):
        integrate_semisynthetic_dataset(make_dataset(), _pls_params(forced_components=forced))


def test_pls_integration_rejects_a_non_integer_forced_component_count() -> None:
    with pytest.raises(SimulationEvaluationError, match="must be an integer"):
        integrate_semisynthetic_dataset(make_dataset(), _pls_params(forced_components=3.5))


def test_unsupported_integration_method_is_rejected() -> None:
    params = SimulationEvaluationParams(integration_method="bad")  # type: ignore[arg-type]

    with pytest.raises(SimulationEvaluationError, match="Unsupported integration_method"):
        evaluate_semisynthetic_trajectory(make_dataset(), params)


def test_invalid_snf_parameters_are_rejected() -> None:
    params = SimulationEvaluationParams(integration_method="snf", integration_params={"K": 12})

    with pytest.raises(SimulationEvaluationError, match="'K' must be at most 11"):
        integrate_semisynthetic_dataset(make_dataset(), params)


def test_missing_metadata_columns_are_rejected() -> None:
    dataset = make_dataset()
    bad = SemiSyntheticTrajectoryDataset(
        methylation=dataset.methylation,
        expression=dataset.expression,
        proteomics=dataset.proteomics,
        metadata=dataset.metadata.drop(columns=["group"]),
        truth=dataset.truth,
    )

    with pytest.raises(SimulationEvaluationError, match="missing required column"):
        evaluate_semisynthetic_trajectory(bad)


def test_omics_metadata_row_mismatch_is_rejected() -> None:
    dataset = make_dataset()
    bad = SemiSyntheticTrajectoryDataset(
        methylation=dataset.methylation.iloc[::-1],
        expression=dataset.expression,
        proteomics=dataset.proteomics,
        metadata=dataset.metadata,
        truth=dataset.truth,
    )

    with pytest.raises(SimulationEvaluationError, match="not aligned"):
        evaluate_semisynthetic_trajectory(bad)


def test_design_objects_and_contrast_are_built_from_sorted_levels() -> None:
    design = build_simulation_trajectory_design(make_dataset().metadata)

    assert design.group_levels == ["A", "B"]
    assert design.stage_levels == ["0", "1", "2"]
    assert design.model_full.shape == (12, 6)
    assert design.model_reduced.shape == (12, 4)
    assert design.ls_means.shape == (6, 6)
    assert design.contrast == [[0, 1, 2], [3, 4, 5]]


def test_observed_only_evaluation_returns_statistics_without_p_values() -> None:
    result = evaluate_semisynthetic_trajectory(
        make_dataset(),
        SimulationEvaluationParams(integration_method="concat", permutations=0),
    )

    assert result.observed_deltas.shape == (2, 2)
    assert result.observed_angles.shape == (2, 2)
    assert result.observed_shapes.shape == (2, 2)
    assert set(result.pair_statistics) == {"delta", "angle", "shape"}
    assert result.p_values == {}
    assert result.truth_metadata["trajectory_mode"] == "magnitude"
    assert result.runtime_metadata["permutations"] == 0


def test_rrpp_p_values_use_plus_one_correction_and_seed() -> None:
    params = SimulationEvaluationParams(
        integration_method="concat",
        permutations=3,
        seed=123,
        include_null_distributions=True,
    )
    first = evaluate_semisynthetic_trajectory(make_dataset(), params)
    second = evaluate_semisynthetic_trajectory(make_dataset(), params)

    assert first.p_values == second.p_values
    assert first.null_distributions == second.null_distributions
    assert set(first.p_values) == {"delta", "angle", "shape"}
    for statistic, observed in first.pair_statistics.items():
        assert first.null_distributions is not None
        null_values = np.asarray(first.null_distributions[statistic], dtype=float)
        expected = float((1.0 + np.sum(null_values >= observed)) / 4.0)
        assert first.p_values[statistic] == expected
        assert 0.0 < first.p_values[statistic] <= 1.0


def test_shape_pair_statistic_is_unavailable_for_two_stages() -> None:
    result = evaluate_semisynthetic_trajectory(
        make_dataset(n_stages=2),
        SimulationEvaluationParams(integration_method="concat", permutations=2, seed=321),
    )

    assert np.isnan(result.pair_statistics["shape"])
    assert "shape" not in result.p_values
    assert not result.runtime_metadata["shape_available"]


def test_truth_metadata_is_json_serializable_without_indicator_arrays() -> None:
    """The persisted truth_metadata must be JSON-safe (no raw numpy indicator arrays)."""

    import json

    dataset = make_dataset()
    # mimic the real generator's truth: JSON-able summary plus raw ndarray indicators
    dataset = SemiSyntheticTrajectoryDataset(
        methylation=dataset.methylation,
        expression=dataset.expression,
        proteomics=dataset.proteomics,
        metadata=dataset.metadata,
        truth={
            "trajectory_mode": "orientation",
            "indicator_counts": {"A": {"methylation": [1, 2, 1]}},
            "indicators": {"A": {"methylation": np.zeros((4, 3))}},
        },
    )
    result = evaluate_semisynthetic_trajectory(dataset, SimulationEvaluationParams(permutations=0))
    assert "indicators" not in result.truth_metadata
    assert result.truth_metadata["indicator_counts"] == {"A": {"methylation": [1, 2, 1]}}
    # round-trips through JSON (this is what the study's JSONL persistence does)
    json.dumps(result.truth_metadata)


def test_null_summary_reports_moments_and_quantiles_matching_the_draws() -> None:
    params = SimulationEvaluationParams(
        integration_method="concat",
        permutations=8,
        seed=2026,
        include_null_distributions=True,
    )
    result = evaluate_semisynthetic_trajectory(make_dataset(), params)

    assert result.null_distributions is not None
    assert set(result.null_summary) == set(result.null_distributions)
    for statistic, draws in result.null_distributions.items():
        values = np.asarray(draws, dtype=float)
        entry = result.null_summary[statistic]
        assert entry["count"] == float(values.size)
        assert entry["mean"] == pytest.approx(float(values.mean()))
        assert entry["sd"] == pytest.approx(float(values.std(ddof=1)))
        for name, q in (("q50", 0.50), ("q90", 0.90), ("q95", 0.95), ("q99", 0.99)):
            assert entry[name] == pytest.approx(float(np.quantile(values, q)))
        # JSON-safe scalars only — this is what the study's JSONL persistence writes.
        assert all(isinstance(value, float) and np.isfinite(value) for value in entry.values())


def test_null_summary_excludes_non_finite_draws_and_reduces_the_count() -> None:
    from motco.simulations.evaluation import _summarize_null_distributions

    clean = [1.0, 2.0, 3.0, 4.0]
    dirty = [1.0, np.nan, 2.0, np.inf, 3.0, -np.inf, 4.0]
    summary = _summarize_null_distributions({"angle": dirty, "delta": clean})

    assert len(dirty) == 7
    assert summary["angle"]["count"] == 4.0  # three non-finite draws dropped
    assert summary["delta"]["count"] == 4.0
    # The retained draws are exactly the finite ones, so both summaries agree.
    assert summary["angle"] == summary["delta"]
    assert summary["angle"]["mean"] == pytest.approx(2.5)


def test_null_summary_omits_undefined_moments_rather_than_emitting_nan() -> None:
    import json

    from motco.simulations.evaluation import _summarize_null_distributions

    summary = _summarize_null_distributions({"empty": [np.nan], "single": [7.0]})

    assert summary["empty"] == {"count": 0.0}
    assert "sd" not in summary["single"]
    assert summary["single"]["mean"] == pytest.approx(7.0)
    # NaN would round-trip through ``json`` as a non-conforming literal.
    assert "NaN" not in json.dumps(summary)


def test_null_summary_is_independent_of_full_distribution_retention() -> None:
    base = dict(integration_method="concat", permutations=4, seed=99)
    without = evaluate_semisynthetic_trajectory(
        make_dataset(), SimulationEvaluationParams(**base, include_null_distributions=False)
    )
    with_full = evaluate_semisynthetic_trajectory(
        make_dataset(), SimulationEvaluationParams(**base, include_null_distributions=True)
    )

    assert without.null_distributions is None
    assert without.null_summary
    assert with_full.null_distributions is not None
    assert with_full.null_summary == without.null_summary


def test_null_summary_is_inert_and_absent_without_permutations() -> None:
    """The summary must not perturb the test it describes.

    The expected values are the ones this configuration produced before the
    summary existed; they are pinned as literals so a future change to the
    summary cannot silently move the statistics or p-values it summarizes.
    """

    # The literals are PC1 values: pin the estimator they were produced with.
    params = SimulationEvaluationParams(
        integration_method="concat", permutations=5, seed=2026, orientation_estimator="pc1"
    )
    first = evaluate_semisynthetic_trajectory(make_dataset(), params)
    second = evaluate_semisynthetic_trajectory(make_dataset(), params)

    assert first.p_values == second.p_values
    assert first.pair_statistics == second.pair_statistics
    assert first.null_summary == second.null_summary
    assert first.p_values == pytest.approx(PRE_CHANGE_P_VALUES)
    assert first.pair_statistics == pytest.approx(PRE_CHANGE_PAIR_STATISTICS)
    assert first.runtime_metadata["permutations"] == 5
    assert first.group_levels == ["A", "B"]
    assert first.contrast == [[0, 1, 2], [3, 4, 5]]

    no_perms = evaluate_semisynthetic_trajectory(
        make_dataset(), SimulationEvaluationParams(integration_method="concat", permutations=0)
    )
    assert no_perms.null_summary == {}
    assert no_perms.p_values == {}


# ── Latent configuration spectrum ─────────────────────────────────────────────


def test_config_spectrum_accompanies_every_evaluation() -> None:
    import json

    result = evaluate_semisynthetic_trajectory(
        make_dataset(), SimulationEvaluationParams(integration_method="concat", permutations=0)
    )

    block = result.config_spectrum
    assert set(block) == {"version", "pooled", "groups"}
    assert set(block["groups"]) == {"A", "B"}
    for entry in [block["pooled"], *block["groups"].values()]:
        assert 0.0 <= entry["relative_eigengap"] <= 1.0
        assert sum(entry["spectrum"]) == pytest.approx(1.0)
    # JSON-safe: no NaN/Inf literals reach the persisted record.
    assert "NaN" not in json.dumps(block)


def test_config_spectrum_matches_an_independent_recomputation() -> None:
    from motco.stats.trajectory import configuration_spectrum, estimate_betas

    dataset = make_dataset()
    params = SimulationEvaluationParams(integration_method="concat", permutations=0)
    result = evaluate_semisynthetic_trajectory(dataset, params)

    latent = integrate_semisynthetic_dataset(dataset, params).matrix
    design = build_simulation_trajectory_design(dataset.metadata)
    obs_vect = np.asarray(design.ls_means, dtype=float) @ np.asarray(
        estimate_betas(design.model_full, latent), dtype=float
    )
    pooled = (obs_vect[design.contrast[0], :] + obs_vect[design.contrast[1], :]) / 2.0

    assert result.config_spectrum["pooled"] == configuration_spectrum(pooled)
    for level, levels in zip(result.group_levels, design.contrast):
        assert result.config_spectrum["groups"][level] == configuration_spectrum(obs_vect[levels, :])


def test_permutation_eigengap_summary_accompanies_rrpp_runs() -> None:
    result = evaluate_semisynthetic_trajectory(
        make_dataset(),
        SimulationEvaluationParams(integration_method="concat", permutations=6, seed=5),
    )

    summary = result.config_spectrum["permutation_pooled_eigengap"]
    assert summary["count"] == 6.0
    assert set(summary) == {"count", "mean", "sd", "q05", "q50", "q95"}
    assert 0.0 <= summary["mean"] <= 1.0
    # Only the summary survives — never the per-permutation vectors.
    assert not any(
        isinstance(value, list) for value in result.config_spectrum.values()
    )


def test_permutation_eigengap_summary_is_absent_without_permutations() -> None:
    result = evaluate_semisynthetic_trajectory(
        make_dataset(), SimulationEvaluationParams(integration_method="concat", permutations=0)
    )

    assert "permutation_pooled_eigengap" not in result.config_spectrum
    assert result.config_spectrum["pooled"] is not None


def test_permutation_eigengap_summary_drops_undefined_draws() -> None:
    from motco.simulations.evaluation import _summarize_permutation_eigengaps

    summary = _summarize_permutation_eigengaps([0.2, None, 0.4, float("nan")])

    assert summary["count"] == 2.0
    assert summary["mean"] == pytest.approx(0.3)
    assert _summarize_permutation_eigengaps([None]) == {"count": 0.0}


def test_config_spectrum_recording_leaves_the_test_untouched() -> None:
    """Recording the covariate moves no statistic, draw, or p-value.

    The literals are the values this configuration produced before the spectrum
    existed (the same ones the null-summary change pinned), so the inertness
    claim is checked against the pre-change pipeline and not merely against a
    second run of the current one.
    """

    params = SimulationEvaluationParams(
        integration_method="concat",
        permutations=5,
        seed=2026,
        include_null_distributions=True,
        orientation_estimator="pc1",
    )
    result = evaluate_semisynthetic_trajectory(make_dataset(), params)

    assert result.p_values == pytest.approx(PRE_CHANGE_P_VALUES)
    assert result.pair_statistics == pytest.approx(PRE_CHANGE_PAIR_STATISTICS)
    assert result.null_summary["angle"]["count"] == 5.0
    assert result.config_spectrum["permutation_pooled_eigengap"]["count"] == 5.0
