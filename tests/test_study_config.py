from __future__ import annotations

import json
from pathlib import Path

import pytest

from motco.simulations import SemiSyntheticTrajectoryParams, SimulationEvaluationParams
from motco.simulations.study import (
    AcceptanceTargets,
    PowerMonotonicityTarget,
    ReportContract,
    SpecificityTarget,
    StudyConfig,
    StudyConfigError,
    TypeIControlTarget,
    dump_study_config,
    enumerate_study,
    load_study_config,
)
from motco.simulations.study.enumerate import NEGATIVE_CONTROL_MODES


def _baseline_config(**overrides) -> StudyConfig:
    defaults: dict = {
        "generator": SemiSyntheticTrajectoryParams(seed=2, trajectory_mode="magnitude", n_samples=60),
        "evaluation": SimulationEvaluationParams(integration_method="concat", permutations=0, seed=3),
        "trajectory_modes": ("magnitude", "orientation", "shape"),
        "effect_sizes": (0.1, 0.5),
        "axes": {"generator.n_samples": (60, 120)},
        "n_replicates": 2,
        "base_seed": 100,
        "alpha": 0.05,
        "acceptance": AcceptanceTargets(
            type_i=(TypeIControlTarget(alpha=0.05),),
            power=(PowerMonotonicityTarget(trajectory_mode="magnitude", statistic="delta", min_power_at_top=0.8),),
            specificity=(SpecificityTarget(trajectory_mode="translation", statistic="angle", alpha=0.05),),
        ),
    }
    defaults.update(overrides)
    return StudyConfig(**defaults)


def test_enumerate_study_is_deterministic_and_negative_controls_present() -> None:
    config = _baseline_config()
    first = enumerate_study(config)
    second = enumerate_study(config)

    assert [c.cell_id for c in first.cells] == [c.cell_id for c in second.cells]

    null_cells = [c for c in first.cells if c.phase.startswith("type_i_")]
    null_modes = {c.generator_params.trajectory_mode for c in null_cells}
    for mode in NEGATIVE_CONTROL_MODES:
        assert mode in null_modes, f"negative-control mode {mode!r} missing from Type I cells"


def test_enumerate_study_unique_cell_ids() -> None:
    grid = enumerate_study(_baseline_config())
    ids = [c.cell_id for c in grid.cells]
    assert len(ids) == len(set(ids))


def test_invalid_axes_namespace_is_rejected() -> None:
    with pytest.raises(StudyConfigError, match="namespace"):
        _baseline_config(axes={"intersim.n_sample": (20, 40)})


def test_unknown_trajectory_mode_is_rejected() -> None:
    with pytest.raises(StudyConfigError, match="Unknown trajectory mode"):
        _baseline_config(trajectory_modes=("not_a_mode",))


def test_negative_n_replicates_is_rejected() -> None:
    with pytest.raises(StudyConfigError, match="non-negative"):
        _baseline_config(n_replicates=-1)


def test_empty_trajectory_modes_or_effect_sizes_is_rejected() -> None:
    with pytest.raises(StudyConfigError, match="trajectory_modes"):
        _baseline_config(trajectory_modes=())
    with pytest.raises(StudyConfigError, match="effect_sizes"):
        _baseline_config(effect_sizes=())


def test_alpha_out_of_range_is_rejected() -> None:
    with pytest.raises(StudyConfigError, match="alpha"):
        _baseline_config(alpha=0.0)
    with pytest.raises(StudyConfigError, match="alpha"):
        _baseline_config(alpha=1.0)


def test_load_study_config_json(tmp_path: Path) -> None:
    path = tmp_path / "study.json"
    payload = {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "group_ratio": 0.5, "n_samples": 60},
        "evaluation": {"integration_method": "concat", "permutations": 0, "seed": 3},
        "trajectory_modes": ["magnitude", "translation"],
        "effect_sizes": [0.1, 0.5],
        "axes": {"generator.n_samples": [60, 120]},
        "n_replicates": 2,
        "base_seed": 1,
        "alpha": 0.05,
        "acceptance": {
            "type_i": [{"alpha": 0.05, "se_tolerance": 2.0}],
            "power": [{"trajectory_mode": "magnitude", "statistic": "delta", "min_power_at_top": 0.8}],
            "specificity": [{"trajectory_mode": "translation", "statistic": "angle", "alpha": 0.05}],
        },
        "metadata": {"study": "smoke"},
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    config = load_study_config(path)
    assert config.generator.seed == 2
    assert config.generator.trajectory_mode == "magnitude"
    assert config.trajectory_modes == ("magnitude", "translation")
    assert config.effect_sizes == (0.1, 0.5)
    assert config.axes["generator.n_samples"] == (60, 120)
    assert len(config.acceptance.type_i) == 1
    assert config.acceptance.power[0].trajectory_mode == "magnitude"
    assert config.metadata["study"] == "smoke"


def test_load_study_config_missing_required_field(tmp_path: Path) -> None:
    path = tmp_path / "study.json"
    path.write_text(json.dumps({"trajectory_modes": ["magnitude"], "effect_sizes": [0.1]}), encoding="utf-8")
    with pytest.raises(StudyConfigError, match="missing required field"):
        load_study_config(path)


def test_load_study_config_unknown_generator_field(tmp_path: Path) -> None:
    path = tmp_path / "study.json"
    payload = {
        "generator": {"seed": 2, "not_a_field": 7},
        "evaluation": {},
        "trajectory_modes": ["magnitude"],
        "effect_sizes": [0.1],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(StudyConfigError, match="unknown field"):
        load_study_config(path)


def _contract_payload(**contract) -> dict:
    return {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "group_ratio": 0.5, "n_samples": 60},
        "evaluation": {"integration_method": "concat", "permutations": 0, "seed": 3, "n_jobs": 1},
        "trajectory_modes": ["magnitude", "translation"],
        "effect_sizes": [0.1, 0.5],
        "report_contract": contract,
    }


def test_report_contract_is_loaded_and_round_trips(tmp_path: Path) -> None:
    path = tmp_path / "study.json"
    payload = _contract_payload(
        driver_component="observed",
        cross_replicate_driver_agreement="descriptive",
        n_jobs_override="forbid",
    )
    path.write_text(json.dumps(payload), encoding="utf-8")
    config = load_study_config(path)
    contract = config.report_contract
    assert contract is not None
    assert contract.driver_component == "observed"
    assert contract.cross_replicate_driver_agreement == "descriptive"
    assert contract.n_jobs_override == "forbid"
    assert contract.forbids_n_jobs_override

    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert json.loads(dumped.read_text(encoding="utf-8"))["report_contract"] == {
        "driver_component": "observed",
        "cross_replicate_driver_agreement": "descriptive",
        "n_jobs_override": "forbid",
    }
    assert load_study_config(dumped) == config


def test_report_contract_defaults_to_warn_and_descriptive(tmp_path: Path) -> None:
    path = tmp_path / "study.json"
    path.write_text(json.dumps(_contract_payload(driver_component="pls_captured")), encoding="utf-8")
    contract = load_study_config(path).report_contract
    assert contract == ReportContract(driver_component="pls_captured")
    assert contract.n_jobs_override == "warn"
    assert not contract.forbids_n_jobs_override


def test_absent_report_contract_is_none_and_round_trips(tmp_path: Path) -> None:
    payload = _contract_payload(driver_component="observed")
    del payload["report_contract"]
    path = tmp_path / "study.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    config = load_study_config(path)
    assert config.report_contract is None
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    assert json.loads(dumped.read_text(encoding="utf-8"))["report_contract"] is None
    assert load_study_config(dumped) == config


@pytest.mark.parametrize(
    ("contract", "match"),
    [
        ({"driver_component": "observed", "not_a_field": 1}, "unknown field.*not_a_field"),
        ({"driver_component": "observed", "n_jobs_override": "ignore"}, "n_jobs_override 'ignore'"),
        ({"driver_component": "observed", "cross_replicate_driver_agreement": "claim"}, "never claimable"),
        ({"driver_component": "latent"}, "driver_component 'latent'"),
        ({"n_jobs_override": "forbid"}, "driver_component is required"),
    ],
)
def test_report_contract_rejects_unknown_fields_and_values(tmp_path: Path, contract: dict, match: str) -> None:
    path = tmp_path / "study.json"
    path.write_text(json.dumps(_contract_payload(**contract)), encoding="utf-8")
    with pytest.raises(StudyConfigError, match=match):
        load_study_config(path)


def test_report_contract_dataclass_validates_directly() -> None:
    with pytest.raises(StudyConfigError, match="never claimable"):
        ReportContract(driver_component="observed", cross_replicate_driver_agreement="claim")
    with pytest.raises(StudyConfigError, match="n_jobs_override"):
        ReportContract(driver_component="observed", n_jobs_override="forbid ")


def test_unknown_root_key_is_rejected_by_name(tmp_path: Path) -> None:
    payload = _contract_payload(driver_component="observed")
    del payload["report_contract"]
    payload["report_contrat"] = {"driver_component": "observed"}  # misspelled block
    path = tmp_path / "study.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(StudyConfigError, match=r"unknown root key\(s\): \['report_contrat'\]"):
        load_study_config(path)


def test_every_committed_example_config_loads() -> None:
    config_dir = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
    names = sorted(p.name for p in config_dir.glob("*.json"))
    assert names, "no committed example configs found"
    for name in names:
        load_study_config(config_dir / name)


# --- acceptance target bound form ------------------------------------------------


def _bound_payload(type_i: dict, specificity: dict) -> dict:
    return {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "n_samples": 60},
        "evaluation": {"integration_method": "concat", "permutations": 0, "seed": 3},
        "trajectory_modes": ["magnitude", "translation"],
        "effect_sizes": [0.1, 0.5],
        "n_replicates": 2,
        "base_seed": 1,
        "alpha": 0.05,
        "acceptance": {
            "type_i": [type_i],
            "specificity": [specificity],
        },
    }


def test_target_bound_defaults_to_two_sided_and_loads_upper(tmp_path: Path) -> None:
    path = tmp_path / "study.json"
    payload = _bound_payload(
        {"alpha": 0.05},
        {"trajectory_mode": "translation", "statistic": "angle", "alpha": 0.05, "bound": "upper"},
    )
    path.write_text(json.dumps(payload), encoding="utf-8")
    config = load_study_config(path)
    assert config.acceptance.type_i[0].bound == "two_sided"
    assert config.acceptance.specificity[0].bound == "upper"


@pytest.mark.parametrize("bound", ["two_sided", "upper"])
def test_target_bound_round_trips_through_dump(tmp_path: Path, bound: str) -> None:
    config = _baseline_config(
        acceptance=AcceptanceTargets(
            type_i=(TypeIControlTarget(alpha=0.05, bound=bound),),
            specificity=(SpecificityTarget(trajectory_mode="translation", statistic="angle", alpha=0.05, bound=bound),),
        )
    )
    path = tmp_path / "dump.json"
    dump_study_config(config, path)
    dumped = json.loads(path.read_text(encoding="utf-8"))
    assert dumped["acceptance"]["type_i"][0]["bound"] == bound
    assert dumped["acceptance"]["specificity"][0]["bound"] == bound
    reloaded = load_study_config(path)
    assert reloaded.acceptance.type_i == config.acceptance.type_i
    assert reloaded.acceptance.specificity == config.acceptance.specificity


@pytest.mark.parametrize(
    ("type_i", "specificity", "match"),
    [
        (
            {"alpha": 0.05, "bound": "one_sided"},
            {"trajectory_mode": "translation", "statistic": "angle", "alpha": 0.05},
            r"acceptance\.type_i\[0\]\.bound .*'one_sided'",
        ),
        (
            {"alpha": 0.05},
            {"trajectory_mode": "translation", "statistic": "angle", "alpha": 0.05, "bound": "lower"},
            r"acceptance\.specificity\[0\]\.bound .*'lower'",
        ),
    ],
)
def test_unknown_target_bound_is_rejected(tmp_path: Path, type_i: dict, specificity: dict, match: str) -> None:
    path = tmp_path / "study.json"
    path.write_text(json.dumps(_bound_payload(type_i, specificity)), encoding="utf-8")
    with pytest.raises(StudyConfigError, match=match):
        load_study_config(path)


# --- column-scoped Type I entries ---------------------------------------------


def _type_i_payload(*entries: dict) -> dict:
    return {
        "generator": {"seed": 2, "trajectory_mode": "magnitude", "group_ratio": 0.5, "n_samples": 60},
        "evaluation": {"integration_method": "concat", "permutations": 0, "seed": 3, "n_jobs": 1},
        "trajectory_modes": ["magnitude", "translation"],
        "effect_sizes": [0.1, 0.5],
        "acceptance": {"type_i": list(entries)},
    }


def test_column_scoped_type_i_entry_loads_and_round_trips(tmp_path: Path) -> None:
    path = tmp_path / "study.json"
    entry = {
        "alpha": 0.05,
        "se_tolerance": 2.0,
        "bound": "upper",
        "trajectory_mode": "translation",
        "statistic": "angle",
    }
    path.write_text(json.dumps(_type_i_payload({"alpha": 0.05}, entry)), encoding="utf-8")
    config = load_study_config(path)
    phase_scoped, column_scoped = config.acceptance.type_i
    assert not phase_scoped.column_scoped
    assert phase_scoped.trajectory_mode is None and phase_scoped.statistic is None
    assert column_scoped.column_scoped
    assert column_scoped.trajectory_mode == "translation" and column_scoped.statistic == "angle"
    assert column_scoped.bound == "upper"

    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    echo = json.loads(dumped.read_text(encoding="utf-8"))["acceptance"]["type_i"]
    # The column-scoped fields are written only when set.
    assert echo[0] == {
        "alpha": 0.05,
        "se_tolerance": 2.0,
        "name": "type_i_control",
        "bound": "two_sided",
        "kind": "type_i_control",
    }
    assert echo[1]["trajectory_mode"] == "translation" and echo[1]["statistic"] == "angle"
    assert load_study_config(dumped) == config


def test_column_scoped_none_entry_without_statistic_loads() -> None:
    target = TypeIControlTarget(alpha=0.05, trajectory_mode="none")
    assert target.column_scoped and target.statistic is None


@pytest.mark.parametrize(
    ("entry", "message"),
    [
        ({"alpha": 0.05, "trajectory_mode": "magnitude"}, "negative-control mode"),
        ({"alpha": 0.05, "trajectory_mode": "orientation", "statistic": "angle"}, "'orientation'"),
        ({"alpha": 0.05, "statistic": "angle"}, "requires trajectory_mode"),
        ({"alpha": 0.05, "trajectory_mode": "translation", "statistic": "size"}, "'size'"),
    ],
)
def test_invalid_column_scoped_type_i_entry_is_rejected_naming_the_entry(
    tmp_path: Path, entry: dict, message: str
) -> None:
    path = tmp_path / "study.json"
    path.write_text(json.dumps(_type_i_payload({"alpha": 0.05}, entry)), encoding="utf-8")
    with pytest.raises(StudyConfigError, match=r"acceptance\.type_i\[1\]") as excinfo:
        load_study_config(path)
    assert message in str(excinfo.value)


def test_type_i_entries_without_a_mode_dump_without_the_column_scoped_keys(tmp_path: Path) -> None:
    path = tmp_path / "study.json"
    path.write_text(
        json.dumps(_type_i_payload({"alpha": 0.05}, {"alpha": 0.05, "bound": "upper", "name": "u"})), encoding="utf-8"
    )
    dumped = tmp_path / "dumped.json"
    dump_study_config(load_study_config(path), dumped)
    for entry in json.loads(dumped.read_text(encoding="utf-8"))["acceptance"]["type_i"]:
        assert "trajectory_mode" not in entry and "statistic" not in entry
        assert set(entry) == {"alpha", "se_tolerance", "name", "bound", "kind"}
