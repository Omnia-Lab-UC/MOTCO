"""Record-level consistency checks for the net-displacement power comparison."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from motco.simulations.grid import SimulationReplicateResult
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("net_displacement_power_checks")
    finally:
        sys.path.pop(0)


checks = _script()
BASE = {"evaluation.integration_method": "pls", "evaluation.orientation_estimator": None}
ND = {"evaluation.integration_method": "pls", "evaluation.orientation_estimator": "net_displacement"}


def _record(
    point, phase, mode, effect, index, *, angle=0.3, delta=0.04, shape=0.5, seed=None, stat=1.0, magnitude_kind=None
):
    meta = {"trajectory_mode": mode, "effect_size": effect}
    transform = {} if magnitude_kind is None else {"magnitude_kind": magnitude_kind}
    if point is not None:
        meta[DESIGN_POINT_KEY] = dict(point)
    return SimulationReplicateResult(
        cell_id=f"{phase}-{mode}-{effect}",
        phase=phase,
        replicate_index=index,
        replicate_seed=index,
        generator_seed=1000 + index if seed is None else seed,
        evaluation_seed=3,
        parameter_signature="sig",
        status="completed",
        p_values={"delta": delta, "angle": angle, "shape": shape},
        pair_statistics={"delta": stat, "angle": 2.0 * stat, "shape": 0.1 * stat},
        truth_metadata={"magnitude_kind": magnitude_kind or "joint", "transform": transform},
        cell_metadata=meta,
    )


def _run(n: int = 5, **candidate_overrides) -> list[SimulationReplicateResult]:
    records = []
    for mode, effect in (("none", 0.0), ("orientation", 1.0)):
        for index in range(n):
            records.append(_record(BASE, "power_primary", mode, effect, index))
            records.append(_record(ND, DESIGN_PHASE, mode, effect, index, angle=0.01, **candidate_overrides))
    return records


def test_shared_statistics_identical_with_counts() -> None:
    result = checks.check_shared_statistics(_run())
    assert result["status"] == "identical"
    (column,) = result["columns"]
    assert column["design_point"] == ND
    assert column["pairs_compared"] == 10 and column["comparisons"] == 20
    assert column["p_value_mismatches"] == 0 and column["observed_statistic_mismatches"] == 0
    assert column["unpaired_records"] == 0


def test_shared_statistics_mismatch_is_counted() -> None:
    records = _run()
    records[3].p_values["shape"] = 0.9  # a candidate record
    records[5].pair_statistics["delta"] = 1.5
    result = checks.check_shared_statistics(records)
    assert result["status"] == "mismatch"
    (column,) = result["columns"]
    assert column["p_value_mismatches"] == 1 and column["observed_statistic_mismatches"] == 1
    assert len(column["examples"]) == 2


def test_shared_observed_rounding_is_within_tolerance() -> None:
    records = _run()
    records[3].pair_statistics["delta"] *= 1 + 1e-14  # a candidate record, BLAS rounding
    result = checks.check_shared_statistics(records)
    assert result["status"] == "identical"
    (column,) = result["columns"]
    assert column["observed_statistic_mismatches"] == 0
    assert column["observed_statistics_bit_identical"] == column["comparisons"] - 1
    assert column["observed_statistics_max_relative_difference"] == pytest.approx(1e-14, rel=1e-2)
    records[3].pair_statistics["delta"] *= 1 + 1e-6
    assert checks.check_shared_statistics(records)["status"] == "mismatch"


def test_angle_differences_are_not_shared_statistics() -> None:
    # The candidate's angle p-values (0.01) differ from the baseline's (0.3) by design.
    assert checks.check_shared_statistics(_run())["status"] == "identical"


def test_parent_reproduction_counts_and_float_tolerance() -> None:
    records = _run()
    parent = [_record(None, "power_primary", "none", 0.0, index) for index in range(5)]
    parent[0].pair_statistics["delta"] *= 1 + 1e-12
    result = checks.check_parent_reproduction(records, {"parent": parent})
    assert result["status"] == "reproduced"
    assert result["baseline_records"] == 10 and result["baseline_records_with_parent"] == 5
    assert result["baseline_cells_with_uncovered_records"] == [{"trajectory_mode": "orientation", "effect_size": 1.0}]
    (entry,) = result["parents"]
    assert entry["records_compared"] == 5
    assert entry["observed_statistics_bit_identical"] == 4
    assert entry["observed_statistics_max_relative_difference"] == pytest.approx(1e-12, rel=1e-3)
    assert entry["cells_compared"] == [{"trajectory_mode": "none", "effect_size": 0.0}]


def test_parent_seed_and_p_value_mismatches_fail() -> None:
    records = _run()
    parent = [_record(None, "power_primary", "orientation", 1.0, index) for index in range(5)]
    parent[1] = _record(None, "power_primary", "orientation", 1.0, 1, seed=7)
    parent[2].p_values["angle"] = 0.31
    result = checks.check_parent_reproduction(records, {"parent": parent})
    assert result["status"] == "mismatch"
    (entry,) = result["parents"]
    assert entry["generator_seed_mismatches"] == 1 and entry["p_value_mismatches"] == 1
    assert [example["replicate_index"] for example in entry["examples"]] == [1, 2]


def test_parent_with_other_magnitude_construction_is_excluded() -> None:
    records = []
    for index in range(4):
        records.append(_record(BASE, "power_primary", "magnitude", 1.0, index, magnitude_kind="joint"))
        records.append(_record(BASE, "power_primary", "orientation", 1.0, index))
    # The historical parent: same seeds, but methylation-only magnitude datasets and so other p-values.
    historical = [
        _record(None, "power_primary", "magnitude", 1.0, i, magnitude_kind="all", delta=0.5) for i in range(4)
    ]
    historical += [_record(None, "power_primary", "orientation", 1.0, i) for i in range(4)]
    joint = [_record(None, "power_primary", "magnitude", 1.0, i, magnitude_kind="joint") for i in range(4)]
    result = checks.check_parent_reproduction(records, {"historical": historical, "joint": joint})
    assert result["status"] == "reproduced"
    assert result["baseline_records_with_parent"] == 8
    first, second = result["parents"]
    assert first["records_compared"] == 4 and first["records_excluded_magnitude_construction"] == 4
    assert first["cells_compared"] == [{"trajectory_mode": "orientation", "effect_size": 1.0}]
    assert second["records_compared"] == 4 and second["records_excluded_magnitude_construction"] == 0
    only_historical = checks.check_parent_reproduction(records, {"historical": historical})
    assert only_historical["baseline_cells_with_uncovered_records"] == [
        {"trajectory_mode": "magnitude", "effect_size": 1.0}
    ]


def test_parent_design_columns_are_not_read() -> None:
    records = _run()
    parent = [_record(ND, DESIGN_PHASE, "none", 0.0, index, seed=1) for index in range(5)]
    result = checks.check_parent_reproduction(records, {"parent": parent})
    assert result["baseline_records_with_parent"] == 0 and result["status"] == "reproduced"


def test_cli_exit_status_and_output(tmp_path: Path) -> None:
    def write(path: Path, records) -> Path:
        path.write_text("".join(json.dumps(r.__dict__) + "\n" for r in records), encoding="utf-8")
        return path

    merged = write(tmp_path / "merged.jsonl", _run())
    good = write(tmp_path / "good.jsonl", [_record(None, "power_primary", "none", 0.0, i) for i in range(5)])
    bad_records = [_record(None, "power_primary", "none", 0.0, i) for i in range(5)]
    bad_records[0].p_values["delta"] = 0.5
    bad = write(tmp_path / "bad.jsonl", bad_records)
    out = tmp_path / "checks.json"
    assert checks.main(["--merged", str(merged), "--parent", str(good), "--out", str(out)]) == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["records"] == 20 and payload["shared_statistics"]["status"] == "identical"
    assert checks.main(["--merged", str(merged), "--parent", str(bad)]) == 1
    assert checks.main(["--merged", str(merged)]) == 0


def test_cli_missing_file_is_an_error(tmp_path: Path) -> None:
    assert checks.main(["--merged", str(tmp_path / "absent.jsonl")]) == 2
