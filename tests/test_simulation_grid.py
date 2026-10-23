from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationParams,
    SimulationEvaluationResult,
    SimulationGrid,
    SimulationGridError,
    SimulationReplicateResult,
    SimulationRunConfig,
    append_replicate_results,
    derive_replicate_seed,
    enumerate_power_grid,
    enumerate_type_i_grid,
    make_simulation_cell,
    parameter_signature,
    read_replicate_results,
    rejection_indicator,
    run_simulation_grid,
    run_simulation_replicate,
    summarize_rejection_rates,
)


def baseline_generator() -> SemiSyntheticTrajectoryParams:
    return SemiSyntheticTrajectoryParams(
        seed=2,
        trajectory_mode="magnitude",
        n_samples=60,
        group_effect_size=0.2,
        group_ratio=0.5,
    )


def baseline_evaluation() -> SimulationEvaluationParams:
    return SimulationEvaluationParams(integration_method="concat", permutations=0, seed=3)


def fake_result(
    p_values: dict[str, float] | None = None,
    pair_statistics: dict[str, float] | None = None,
    truth_seed: int = 0,
    null_summary: dict[str, dict[str, float]] | None = None,
    config_spectrum: dict | None = None,
) -> SimulationEvaluationResult:
    matrix = np.zeros((2, 2), dtype=float)
    return SimulationEvaluationResult(
        null_summary=null_summary or {},
        config_spectrum=config_spectrum or {},
        observed_deltas=matrix,
        observed_angles=matrix,
        observed_shapes=matrix,
        pair_statistics=pair_statistics or {"delta": 1.0, "angle": 2.0, "shape": float("nan")},
        p_values=p_values or {"delta": 0.01, "angle": 0.20},
        latent_matrix_metadata={"integration_method": "concat"},
        truth_metadata={"seed": truth_seed},
        runtime_metadata={"runtime_seconds": 0.1},
        evaluation_params=baseline_evaluation(),
        group_levels=["A", "B"],
        stage_levels=["0", "1"],
        contrast=[[0, 1], [2, 3]],
        realized_geometry={
            "schema_version": 1,
            "requested": {"trajectory_mode": "orientation", "group_effect_size": 0.5},
            "checkpoints": {
                "observed_standardized": {
                    "joint": {
                        "path_lengths": {"A": 1.0, "B": 1.0},
                        "delta": 0.0,
                        "angle": 20.0,
                        "shape": None,
                        "availability": {"delta": True, "angle": True, "shape": False},
                    }
                }
            },
        },
    )


def test_type_i_grid_enumeration_is_stable_and_null() -> None:
    axes = {
        "generator.n_samples": [60, 120],
        "generator.group_ratio": [0.5, 0.7],
        "evaluation.permutations": [0, 2],
    }
    first = enumerate_type_i_grid(
        baseline_generator_params=baseline_generator(),
        evaluation_params=baseline_evaluation(),
        axes=axes,
        n_replicates=2,
        base_seed=101,
    )
    second = enumerate_type_i_grid(
        baseline_generator_params=baseline_generator(),
        evaluation_params=baseline_evaluation(),
        axes=axes,
        n_replicates=2,
        base_seed=101,
    )

    assert [cell.cell_id for cell in first.cells] == [cell.cell_id for cell in second.cells]
    assert len(first.cells) == 4
    assert first.cells[0].phase == "type_i_baseline"
    assert all(cell.generator_params.trajectory_mode == "none" for cell in first.cells)
    assert all(cell.generator_params.group_effect_size == 0.0 for cell in first.cells)
    assert {cell.metadata["varied_axis"] for cell in first.cells[1:]} == {
        "generator.n_samples",
        "generator.group_ratio",
        "evaluation.permutations",
    }


def test_power_grid_enumeration_includes_modes_effects_and_axes() -> None:
    grid = enumerate_power_grid(
        baseline_generator_params=baseline_generator(),
        evaluation_params=baseline_evaluation(),
        trajectory_modes=["magnitude", "orientation"],
        effect_sizes=[0.1, 0.2],
        axes={"generator.n_samples": [60, 120]},
    )

    assert len(grid.cells) == 8
    primary = [cell for cell in grid.cells if cell.phase == "power_primary"]
    ofat = [cell for cell in grid.cells if cell.phase == "power_ofat"]
    assert len(primary) == 4
    assert len(ofat) == 4
    assert {(cell.generator_params.trajectory_mode, cell.generator_params.group_effect_size) for cell in primary} == {
        ("magnitude", 0.1),
        ("magnitude", 0.2),
        ("orientation", 0.1),
        ("orientation", 0.2),
    }


def test_invalid_grid_inputs_are_rejected() -> None:
    with pytest.raises(SimulationGridError, match="n_replicates"):
        make_simulation_cell(
            phase="type_i_baseline",
            generator_params=baseline_generator(),
            n_replicates=0,
        )

    with pytest.raises(SimulationGridError, match="namespace prefix"):
        enumerate_type_i_grid(
            baseline_generator_params=baseline_generator(),
            axes={"n_samples": [60, 120]},
        )


def test_replicate_seed_derivation_is_deterministic() -> None:
    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        n_replicates=2,
        base_seed=99,
    )

    assert derive_replicate_seed(cell, 0) == derive_replicate_seed(cell, 0)
    assert derive_replicate_seed(cell, 0) != derive_replicate_seed(cell, 1)
    with pytest.raises(SimulationGridError, match="replicate_index"):
        derive_replicate_seed(cell, 2)


def test_replicate_seed_always_fits_unsigned_int31() -> None:
    SEED_MAX = 2**31 - 1

    seeds: list[int] = []
    for base_seed in (0, 1, 42, 100, 2_797_983_684, -1):
        for cell_index in range(8):
            cell = make_simulation_cell(
                cell_id=f"probe-{cell_index}",
                phase="power_primary",
                generator_params=baseline_generator(),
                n_replicates=16,
                base_seed=base_seed,
            )
            for replicate_index in range(cell.n_replicates):
                seeds.append(derive_replicate_seed(cell, replicate_index))

    assert all(0 <= s <= SEED_MAX for s in seeds)
    # Sanity: at least some draws land in the high half, exercising the mask.
    assert max(seeds) > 2**30


def test_replicate_seed_masks_known_high_bit_value() -> None:
    # An unmasked 32-bit derivation produced 2_797_983_684 (high bit set,
    # = 0xa6b8f084). Clearing the high bit gives 0x26b8f084 = 650_500_036.
    unmasked = 2_797_983_684
    assert unmasked & 0x7FFFFFFF == 650_500_036
    assert 650_500_036 <= 2**31 - 1


def test_parameter_signature_includes_seed_derivation_version(monkeypatch: pytest.MonkeyPatch) -> None:
    import motco.simulations.grid as grid_module

    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        n_replicates=2,
        base_seed=77,
    )
    current_signature = parameter_signature(cell)

    assert parameter_signature(cell) == current_signature

    original_stable_digest = grid_module._stable_digest

    def stable_digest_without_version(payload: object, *, length: int | None = None) -> str:
        if isinstance(payload, dict) and "seed_derivation_version" in payload:
            payload = {k: v for k, v in payload.items() if k != "seed_derivation_version"}
        return original_stable_digest(payload, length=length)

    monkeypatch.setattr(grid_module, "_stable_digest", stable_digest_without_version)
    assert parameter_signature(cell) != current_signature


def test_parameter_signature_includes_geometry_schema_version(monkeypatch: pytest.MonkeyPatch) -> None:
    import motco.simulations.grid as grid_module

    cell = make_simulation_cell(phase="type_i_baseline", generator_params=baseline_generator())
    current_signature = parameter_signature(cell)

    monkeypatch.setattr(grid_module, "DIAGNOSTIC_SCHEMA_VERSION", 2)

    assert parameter_signature(cell) != current_signature


def test_run_replicate_uses_injectable_evaluator_and_records_seeds() -> None:
    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        evaluation_params=baseline_evaluation(),
    )

    def evaluator(
        generator_params: SemiSyntheticTrajectoryParams,
        evaluation_params: SimulationEvaluationParams,
    ) -> SimulationEvaluationResult:
        assert evaluation_params.seed == 3
        return fake_result(truth_seed=generator_params.seed)

    record = run_simulation_replicate(cell, 0, evaluator=evaluator)

    assert record.status == "completed"
    assert record.generator_seed == record.replicate_seed
    assert record.evaluation_seed == 3
    assert record.truth_metadata["seed"] == record.replicate_seed
    assert record.parameter_signature == parameter_signature(cell)


def test_run_replicate_error_policy_records_failures() -> None:
    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
    )

    def evaluator(
        generator_params: SemiSyntheticTrajectoryParams,
        evaluation_params: SimulationEvaluationParams,
    ) -> SimulationEvaluationResult:
        raise RuntimeError("boom")

    record = run_simulation_replicate(cell, 0, evaluator=evaluator, error_policy="record")

    assert record.status == "failed"
    assert record.error_type == "RuntimeError"
    assert record.error_message == "boom"


def test_jsonl_persistence_read_and_resume(tmp_path) -> None:
    path = tmp_path / "results.jsonl"
    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        n_replicates=2,
    )
    grid = SimulationGrid(cells=(cell,))
    calls = 0

    def evaluator(
        generator_params: SemiSyntheticTrajectoryParams,
        evaluation_params: SimulationEvaluationParams,
    ) -> SimulationEvaluationResult:
        nonlocal calls
        calls += 1
        return fake_result(truth_seed=generator_params.seed)

    first = run_simulation_grid(grid, config=SimulationRunConfig(output_path=path), evaluator=evaluator)
    second = run_simulation_grid(grid, config=SimulationRunConfig(output_path=path), evaluator=evaluator)
    loaded = read_replicate_results(path)

    assert len(first) == 2
    assert second == []
    assert calls == 2
    assert len(loaded) == 2
    assert {record.replicate_index for record in loaded} == {0, 1}
    assert loaded[0].realized_geometry["schema_version"] == 1
    assert loaded[0].realized_geometry["checkpoints"]["observed_standardized"]["joint"]["shape"] is None


def test_legacy_jsonl_record_loads_without_geometry(tmp_path) -> None:
    import json

    path = tmp_path / "legacy.jsonl"
    record = run_simulation_replicate(
        make_simulation_cell(phase="type_i_baseline", generator_params=baseline_generator()),
        0,
        evaluator=lambda generator, evaluation: fake_result(),
    )
    payload = record.__dict__.copy()
    payload.pop("realized_geometry")
    path.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    loaded = read_replicate_results(path)

    assert len(loaded) == 1
    assert loaded[0].realized_geometry == {}


def test_resume_detects_parameter_mismatch(tmp_path) -> None:
    path = tmp_path / "results.jsonl"
    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        cell_id="shared",
    )
    append_replicate_results(
        path,
        [
            run_simulation_replicate(
                cell,
                0,
                evaluator=lambda generator, evaluation: fake_result(),
            )
        ],
    )
    changed = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=SemiSyntheticTrajectoryParams(seed=2, n_samples=60, group_ratio=0.7),
        cell_id="shared",
    )

    with pytest.raises(SimulationGridError, match="different parameter signature"):
        run_simulation_grid(SimulationGrid(cells=(changed,)), config=SimulationRunConfig(output_path=path))


def test_rejection_summaries_handle_available_and_missing_statistics() -> None:
    records = [
        SimulationReplicateResult(
            cell_id="cell-a",
            phase="type_i_baseline",
            replicate_index=0,
            replicate_seed=1,
            generator_seed=1,
            evaluation_seed=1,
            parameter_signature="sig",
            status="completed",
            p_values={"delta": 0.01, "angle": 0.20},
        ),
        SimulationReplicateResult(
            cell_id="cell-a",
            phase="type_i_baseline",
            replicate_index=1,
            replicate_seed=2,
            generator_seed=2,
            evaluation_seed=2,
            parameter_signature="sig",
            status="completed",
            p_values={"delta": 0.10, "angle": 0.03},
        ),
    ]

    summaries = summarize_rejection_rates(records, alpha=0.05, statistics=("delta", "angle", "shape"))
    by_stat = {summary.statistic: summary for summary in summaries}

    assert rejection_indicator(0.01, alpha=0.05)
    assert rejection_indicator(None, alpha=0.05) is None
    assert by_stat["delta"].rejection_rate == 0.5
    assert by_stat["delta"].monte_carlo_se == pytest.approx((0.5 * 0.5 / 2) ** 0.5)
    assert by_stat["angle"].rejected_replicates == 1
    assert by_stat["shape"].available_replicates == 0
    assert by_stat["shape"].rejection_rate is None
    assert by_stat["shape"].unavailable_replicates == 2


def sample_null_summary() -> dict[str, dict[str, float]]:
    return {
        "delta": {"count": 199.0, "mean": 0.8, "sd": 0.2, "q50": 0.75, "q90": 1.0, "q95": 1.1, "q99": 1.3},
        "angle": {"count": 197.0, "mean": 41.0, "sd": 9.0, "q50": 40.0, "q90": 53.0, "q95": 57.0, "q99": 64.0},
    }


def test_null_summary_reaches_the_persisted_record(tmp_path) -> None:
    path = tmp_path / "results.jsonl"
    summary = sample_null_summary()
    cell = make_simulation_cell(phase="power", generator_params=baseline_generator(), n_replicates=1)
    run_simulation_grid(
        SimulationGrid(cells=(cell,)),
        config=SimulationRunConfig(output_path=path),
        evaluator=lambda generator, evaluation: fake_result(null_summary=summary),
    )

    loaded = read_replicate_results(path)

    assert len(loaded) == 1
    assert loaded[0].status == "completed"
    assert loaded[0].null_summary == summary
    # The pair (observed statistic, its own critical value) — the quantity the
    # pivotality question is about — is computable from this single record.
    assert loaded[0].pair_statistics["angle"] == 2.0
    assert loaded[0].null_summary["angle"]["q95"] == 57.0


def test_legacy_record_loads_with_empty_null_summary_distinguishable_from_no_permutations(tmp_path) -> None:
    import json

    cell = make_simulation_cell(phase="power", generator_params=baseline_generator())
    record = run_simulation_replicate(
        cell, 0, evaluator=lambda generator, evaluation: fake_result(null_summary=sample_null_summary())
    )

    # A record written before the field existed: the key is simply absent.
    legacy_payload = record.__dict__.copy()
    legacy_payload.pop("null_summary")
    legacy_payload["runtime_metadata"] = {"runtime_seconds": 0.1, "permutations": 199}
    # A record from an evaluation that ran no permutations: the key is present
    # and empty, and the run itself is recorded as having permuted nothing.
    no_perms_payload = record.__dict__.copy()
    no_perms_payload["null_summary"] = {}
    no_perms_payload["runtime_metadata"] = {"runtime_seconds": 0.1, "permutations": 0}

    path = tmp_path / "mixed.jsonl"
    path.write_text(
        json.dumps(legacy_payload) + "\n" + json.dumps(no_perms_payload) + "\n",
        encoding="utf-8",
    )
    legacy, no_perms = read_replicate_results(path)

    assert legacy.null_summary == {}
    assert no_perms.null_summary == {}
    assert legacy.runtime_metadata["permutations"] == 199
    assert no_perms.runtime_metadata["permutations"] == 0


def test_null_summary_does_not_enter_the_parameter_signature(tmp_path) -> None:
    """The record field is derived from draws that already happened.

    ``null_summary`` itself never enters the signature — the resume below, which
    skips a record written without the field, is what pins that. The digest is
    pinned separately so *any* signature change has to be a deliberate edit
    here; it last moved when the attribution diagnostic record gained the
    principal-orientation block (``align-orientation-attribution``, schema
    version 2), whose deliberate consequence is that pre-change shards refuse to
    resume — the params hash *is* the version. It moved again when
    ``magnitude_kind`` defaulted to ``joint`` (``default-joint-magnitude-construction``):
    ``baseline_generator()`` names no kind, so it resolves to the new default.
    Committed configs pin ``"all"`` and keep their signatures
    (``tests/test_magnitude_default_pins.py``). It moved again when
    ``orientation_estimator`` defaulted to ``net_displacement``
    (``default-net-displacement-orientation``): ``baseline_evaluation()`` names no
    estimator, so the default is hashed. Committed configs pin ``"pc1"``
    (``tests/test_orientation_default_pins.py``). It moved again when the PLS
    rank floor ``min_components`` defaulted to 2 (``minimum-latent-rank``):
    ``baseline_evaluation()`` names no floor, so the default is hashed.
    Committed configs pin ``1`` (``tests/test_rank_floor_default_pins.py``).
    """

    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        evaluation_params=baseline_evaluation(),
        n_replicates=2,
        cell_id="pinned",
    )
    assert parameter_signature(cell) == "7cd15434da7da8430461ec10959099126f4e751e915b21aca989b89c1ee120fb"
    # Pinning the former rank-floor default restores the pre-floor digest exactly.
    unfloored = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        evaluation_params=replace(baseline_evaluation(), min_components=1),
        n_replicates=2,
        cell_id="pinned",
    )
    assert parameter_signature(unfloored) == "55e695c3321689aac97cd6de84c63c43f5e9afbb01d6cb96a0a413269287bae2"
    # Pinning the former orientation default too restores the pre-orientation-flip digest exactly.
    pc1_evaluation = replace(baseline_evaluation(), orientation_estimator="pc1", min_components=1)
    pinned_pc1 = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        evaluation_params=pc1_evaluation,
        n_replicates=2,
        cell_id="pinned",
    )
    assert parameter_signature(pinned_pc1) == "655d20b127fc03f6682a90d6db5820973b219e5254d144df6d10db86ca6b282c"
    # Pinning both former defaults restores the pre-magnitude-flip digest exactly.
    pinned_all = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=replace(baseline_generator(), magnitude_kind="all"),
        evaluation_params=pc1_evaluation,
        n_replicates=2,
        cell_id="pinned",
    )
    assert parameter_signature(pinned_all) == "50c2c5416fb127b814dd93b781d686b1c132727b32db6474e11910205223e6b9"

    # Resume against a record written without the field must skip, not overwrite.
    import json

    path = tmp_path / "results.jsonl"
    prior = run_simulation_replicate(
        cell, 0, evaluator=lambda generator, evaluation: fake_result(truth_seed=99)
    )
    payload = prior.__dict__.copy()
    payload.pop("null_summary")
    path.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    calls = 0

    def evaluator(
        generator_params: SemiSyntheticTrajectoryParams,
        evaluation_params: SimulationEvaluationParams,
    ) -> SimulationEvaluationResult:
        nonlocal calls
        calls += 1
        return fake_result(truth_seed=generator_params.seed, null_summary=sample_null_summary())

    run_simulation_grid(
        SimulationGrid(cells=(cell,)), config=SimulationRunConfig(output_path=path), evaluator=evaluator
    )

    loaded = read_replicate_results(path)
    assert calls == 1  # replicate 0 was skipped, only replicate 1 ran
    assert len(loaded) == 2
    by_index = {record.replicate_index: record for record in loaded}
    assert by_index[0].truth_metadata["seed"] == 99  # the pre-change record survived untouched
    assert by_index[0].null_summary == {}
    assert by_index[1].null_summary == sample_null_summary()


# ── Latent configuration spectrum ─────────────────────────────────────────────


def sample_config_spectrum() -> dict:
    return {
        "version": 1,
        "pooled": {
            "n_points": 4,
            "n_dimensions": 6,
            "total_variance": 12.5,
            "spectrum": [0.6, 0.25, 0.1, 0.05],
            "relative_eigengap": 0.35,
        },
        "groups": {
            "A": {
                "n_points": 4,
                "n_dimensions": 6,
                "total_variance": 9.0,
                "spectrum": [0.5, 0.3, 0.15, 0.05],
                "relative_eigengap": 0.2,
            },
            "B": {
                "n_points": 4,
                "n_dimensions": 6,
                "total_variance": 8.0,
                "spectrum": [0.7, 0.2, 0.1],
                "relative_eigengap": 0.5,
            },
        },
        "permutation_pooled_eigengap": {
            "count": 199.0,
            "mean": 0.31,
            "sd": 0.08,
            "q05": 0.19,
            "q50": 0.30,
            "q95": 0.45,
        },
    }


def test_config_spectrum_reaches_the_persisted_record(tmp_path) -> None:
    path = tmp_path / "results.jsonl"
    spectrum = sample_config_spectrum()
    cell = make_simulation_cell(phase="power", generator_params=baseline_generator(), n_replicates=1)

    run_simulation_grid(
        SimulationGrid(cells=(cell,)),
        config=SimulationRunConfig(output_path=path),
        evaluator=lambda generator, evaluation: fake_result(config_spectrum=spectrum),
    )

    (loaded,) = read_replicate_results(path)
    assert loaded.config_spectrum == spectrum


def test_config_spectrum_version_enters_the_parameter_signature(monkeypatch: pytest.MonkeyPatch) -> None:
    """Unlike the null summary, the spectrum schema is signature-bearing.

    A pre-change shard must refuse resume rather than produce a merged set in
    which only some records carry the covariate.
    """

    import motco.simulations.grid as grid_module

    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        n_replicates=2,
        base_seed=77,
    )
    current = parameter_signature(cell)
    assert parameter_signature(cell) == current  # stable across enumerations

    monkeypatch.setattr(grid_module, "CONFIG_SPECTRUM_VERSION", 2)
    assert parameter_signature(cell) != current


def test_resume_against_a_pre_spectrum_signature_is_refused(tmp_path) -> None:
    import json

    cell = make_simulation_cell(
        phase="power_primary",
        generator_params=baseline_generator(),
        evaluation_params=baseline_evaluation(),
        n_replicates=1,
        cell_id="resume-probe",
    )
    prior = run_simulation_replicate(
        cell, 0, evaluator=lambda generator, evaluation: fake_result(truth_seed=7)
    )
    payload = prior.__dict__.copy()
    payload.pop("config_spectrum")
    # A record written under the pre-change contract: no spectrum, and a
    # signature computed without the spectrum schema version.
    payload["parameter_signature"] = "pre-change-signature"
    path = tmp_path / "results.jsonl"
    path.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    with pytest.raises(SimulationGridError, match="different parameter signature"):
        run_simulation_grid(
            SimulationGrid(cells=(cell,)),
            config=SimulationRunConfig(output_path=path),
            evaluator=lambda generator, evaluation: fake_result(),
        )


def test_shape_statistic_version_enters_the_parameter_signature(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The contract of the statistic itself is signature-bearing, not just the record schema.

    ``shape`` under proper-rotation alignment and under full orthogonal
    alignment are different statistics; a shard must not mix them.
    """

    import motco.simulations.grid as grid_module

    cell = make_simulation_cell(
        phase="type_i_baseline",
        generator_params=baseline_generator(),
        n_replicates=2,
        base_seed=91,
    )
    current = parameter_signature(cell)
    assert parameter_signature(cell) == current  # stable across enumerations

    monkeypatch.setattr(grid_module, "SHAPE_STATISTIC_VERSION", 1)
    assert parameter_signature(cell) != current


def test_resume_against_a_pre_reflection_policy_signature_is_refused(tmp_path) -> None:
    """A shard written under the proper-rotation ``shape`` contract must refuse resume."""

    import json

    import motco.simulations.grid as grid_module

    cell = make_simulation_cell(
        phase="power_primary",
        generator_params=baseline_generator(),
        evaluation_params=baseline_evaluation(),
        n_replicates=1,
        cell_id="reflection-resume-probe",
    )
    prior = run_simulation_replicate(
        cell, 0, evaluator=lambda generator, evaluation: fake_result(truth_seed=7)
    )
    payload = prior.__dict__.copy()
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(grid_module, "SHAPE_STATISTIC_VERSION", 1)
        payload["parameter_signature"] = parameter_signature(cell)
    path = tmp_path / "results.jsonl"
    path.write_text(json.dumps(payload) + "\n", encoding="utf-8")

    with pytest.raises(SimulationGridError, match="different parameter signature"):
        run_simulation_grid(
            SimulationGrid(cells=(cell,)),
            config=SimulationRunConfig(output_path=path),
            evaluator=lambda generator, evaluation: fake_result(),
        )


def test_resume_under_the_new_schema_succeeds(tmp_path) -> None:
    path = tmp_path / "results.jsonl"
    cell = make_simulation_cell(
        phase="power_primary",
        generator_params=baseline_generator(),
        evaluation_params=baseline_evaluation(),
        n_replicates=2,
        cell_id="resume-ok",
    )
    calls = 0

    def evaluator(generator_params, evaluation_params) -> SimulationEvaluationResult:
        nonlocal calls
        calls += 1
        return fake_result(config_spectrum=sample_config_spectrum())

    run_simulation_grid(
        SimulationGrid(cells=(cell,)), config=SimulationRunConfig(output_path=path), evaluator=evaluator
    )
    run_simulation_grid(
        SimulationGrid(cells=(cell,)), config=SimulationRunConfig(output_path=path), evaluator=evaluator
    )

    assert calls == 2  # the second run skipped both completed replicates
    assert all(record.config_spectrum for record in read_replicate_results(path))


def test_legacy_record_loads_with_an_empty_spectrum_block(tmp_path) -> None:
    """An absent field and a recorded degenerate spectrum stay distinguishable."""

    import json

    cell = make_simulation_cell(phase="power", generator_params=baseline_generator(), n_replicates=1)
    degenerate = {
        "version": 1,
        "pooled": {
            "n_points": 4,
            "n_dimensions": 6,
            "total_variance": 0.0,
            "spectrum": [],
            "relative_eigengap": None,
        },
        "groups": {},
    }
    record = run_simulation_replicate(
        cell, 0, evaluator=lambda generator, evaluation: fake_result(config_spectrum=degenerate)
    )

    legacy_payload = record.__dict__.copy()
    legacy_payload.pop("config_spectrum")
    path = tmp_path / "mixed.jsonl"
    path.write_text(
        json.dumps(legacy_payload) + "\n" + json.dumps(record.__dict__) + "\n",
        encoding="utf-8",
    )
    legacy, recorded = read_replicate_results(path)

    assert legacy.config_spectrum == {}
    assert recorded.config_spectrum == degenerate
    assert recorded.config_spectrum["pooled"]["relative_eigengap"] is None
