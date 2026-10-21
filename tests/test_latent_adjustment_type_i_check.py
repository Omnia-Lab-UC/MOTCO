"""``scripts/latent_adjustment_type_i_check.py``: reproducibility, arms, resume and the decision rule."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("latent_adjustment_type_i_check")
    finally:
        sys.path.pop(0)


check = _script()
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}


def _settings(out_dir: Path) -> dict:
    args = check.build_parser().parse_args(["run", "--reps", "2", "--perms", "9", "--out-dir", str(out_dir)])
    settings = check.build_settings(args)
    settings["cv"] = dict(SMALL_CV)
    return settings


def _records(out_dir: Path) -> list[dict]:
    lines = check.shard_path(out_dir, 0, 1).read_text(encoding="utf-8").splitlines()
    records = [json.loads(line) for line in lines]
    for record in records:
        record.pop("runtime_seconds")
    return records


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    outputs = []
    for name in ("first", "second"):
        out_dir = tmp_path_factory.mktemp(name)
        settings = _settings(out_dir)
        check.ensure_settings(out_dir, settings)
        assert check.run_shard(out_dir, settings, 0, 1) == 2
        outputs.append((out_dir, settings))
    return outputs


def test_replicates_are_reproducible(runs) -> None:
    (first, _), (second, _) = runs
    assert _records(first) == _records(second)


def test_only_residualized_arms_carry_the_block(runs) -> None:
    out_dir, settings = runs[0]
    for record in _records(out_dir):
        assert set(record["arms"]) == {"baseline", "stage", "group_stage"}
        assert record["arms"]["baseline"]["latent_adjustment"] is None
        for arm in ("stage", "group_stage"):
            block = record["arms"][arm]["latent_adjustment"]
            assert block["conditioning"] == arm
            assert block["covariates"] == check.declared_covariates(settings)
        for arm in check.ARMS:
            entry = record["arms"][arm]
            assert entry["error"] is None
            assert set(entry["p_values"]) == {"delta", "angle", "shape"}
            assert isinstance(entry["selected_lv"], int)
            assert -1.0 <= entry["first_component_spearman_z"] <= 1.0


def test_arms_share_the_data_and_differ_in_the_basis(runs) -> None:
    out_dir, settings = runs[0]
    record = _records(out_dir)[0]
    built = check.build_replicate(record["replicate"], settings)
    assert set(check.declared_covariates(settings)) <= set(built["metadata"].columns)
    assert len(check.declared_covariates(settings)) == 8
    stage = record["arms"]["stage"]["latent_adjustment"]
    group_stage = record["arms"]["group_stage"]["latent_adjustment"]
    n = built["metadata"].shape[0]
    q = 1 + 2 + 6  # z, two category indicators, six noise covariates
    assert stage["residual_df"] == n - 1 - 2 - q
    assert group_stage["residual_df"] == n - 1 - 5 - q
    # The D10 replicate is reused unchanged: same seeds, blocks and z / category.
    reference = check.d10.build_replicate(record["replicate"], settings)
    for name in settings["layers"]:
        assert built["blocks"][name].equals(reference["blocks"][name])
    assert built["metadata"]["z"].equals(reference["metadata"]["z"])


def test_resume_skips_recorded_replicates(runs) -> None:
    out_dir, settings = runs[0]
    before = _records(out_dir)
    assert check.run_shard(out_dir, settings, 0, 1) == 0
    assert _records(out_dir) == before


def test_merge_writes_the_decision(runs) -> None:
    out_dir, settings = runs[0]
    summary = check.merge(out_dir, settings)
    assert summary["replicates"] == 2 and summary["missing"] == []
    decision = json.loads((out_dir / "latent_adjustment_type_i_decision.json").read_text())
    assert decision["decision"] in {"adopt_stage", "adopt_group_stage", "diagnose_latent_excess"}
    assert set(decision["rates"]) == {"baseline", "stage", "group_stage"}
    assert decision["calibration"]["by_arm"]["baseline"]["replicates"] == 2
    assert decision["construction"]["declared_covariates"] == check.declared_covariates(settings)
    for name in ("pvalues", "replicates", "rates"):
        assert (out_dir / f"latent_adjustment_type_i_{name}.csv").exists()


def test_changed_settings_are_refused(runs) -> None:
    out_dir, settings = runs[0]
    with pytest.raises(SystemExit, match="different settings"):
        check.ensure_settings(out_dir, {**settings, "permutations": 19})


# ── Decision rule ────────────────────────────────────────────────────────────


def _uniform(n: int, rejections: int) -> list[float]:
    return [0.01] * rejections + [0.5] * (n - rejections)


def _rows(rates: dict[str, dict[str, int]], n: int = 500) -> list[dict]:
    rows = []
    for arm, per_statistic in rates.items():
        for statistic, rejections in per_statistic.items():
            for replicate, p_value in enumerate(_uniform(n, rejections)):
                rows.append({"replicate": replicate, "arm": arm, "statistic": statistic, "p_value": p_value})
    return rows


HOLD = {"delta": 25, "angle": 30, "shape": 34}  # 0.068 <= 0.0695 at R = 500
EXCESS_SHAPE = {"delta": 25, "angle": 30, "shape": 36}  # 0.072 > 0.0695


def test_both_conditionings_hold_adopts_stage() -> None:
    decision = check.latent_adjustment_type_i_decision(
        _rows({"baseline": HOLD, "stage": HOLD, "group_stage": HOLD}), alpha=0.05, k=2.0
    )
    assert decision["bound"]["stage"]["shape"] == pytest.approx(0.0695, abs=1e-4)
    assert decision["admissible"] == ["stage", "group_stage"]
    assert decision["decision"] == "adopt_stage"


def test_only_group_stage_holds_adopts_group_stage() -> None:
    decision = check.latent_adjustment_type_i_decision(
        _rows({"baseline": HOLD, "stage": EXCESS_SHAPE, "group_stage": HOLD}), alpha=0.05, k=2.0
    )
    assert decision["verdicts"]["stage"]["shape"] == "exceeds"
    assert decision["admissible"] == ["group_stage"]
    assert decision["decision"] == "adopt_group_stage"


def test_no_conditioning_holds_is_diagnosed() -> None:
    decision = check.latent_adjustment_type_i_decision(
        _rows({"baseline": HOLD, "stage": EXCESS_SHAPE, "group_stage": {**HOLD, "angle": 60}}), alpha=0.05, k=2.0
    )
    assert decision["admissible"] == []
    assert decision["decision"] == "diagnose_latent_excess"


def test_baseline_arm_does_not_enter_the_decision() -> None:
    held = check.latent_adjustment_type_i_decision(
        _rows({"baseline": HOLD, "stage": HOLD, "group_stage": HOLD}), alpha=0.05, k=2.0
    )
    inflated = check.latent_adjustment_type_i_decision(
        _rows({"baseline": {"delta": 200, "angle": 200, "shape": 200}, "stage": HOLD, "group_stage": HOLD}),
        alpha=0.05,
        k=2.0,
    )
    assert inflated["rates"]["baseline"]["delta"] == pytest.approx(0.4)
    assert "baseline" not in inflated["verdicts"]
    assert inflated["decision"] == held["decision"] == "adopt_stage"
    assert inflated["admissible"] == held["admissible"]


def test_missing_p_values_reduce_r_and_an_empty_arm_is_unavailable() -> None:
    rows = _rows({"baseline": HOLD, "stage": HOLD})
    rows += [{"replicate": 0, "arm": "stage", "statistic": "delta", "p_value": None}]
    decision = check.latent_adjustment_type_i_decision(rows, alpha=0.05, k=2.0)
    assert decision["replicates"]["stage"]["delta"] == 500
    assert set(decision["verdicts"]["group_stage"].values()) == {"unavailable"}
    assert decision["admissible"] == ["stage"]


def test_calibration_evidence_is_reported_from_the_baseline() -> None:
    rows = _rows({"baseline": HOLD, "stage": HOLD, "group_stage": HOLD})
    met = check.latent_adjustment_type_i_decision(
        rows, alpha=0.05, k=2.0, correlations={"baseline": [0.9, -0.85, 0.7, None], "stage": [0.1, 0.0]}
    )
    assert met["calibration"]["baseline_median_abs"] == pytest.approx(0.85)
    assert met["calibration"]["by_arm"]["baseline"]["replicates"] == 3
    assert met["calibration"]["met"] is True
    unmet = check.latent_adjustment_type_i_decision(rows, alpha=0.05, k=2.0, correlations={"baseline": [0.5, 0.6]})
    assert unmet["calibration"]["met"] is False
    assert unmet["decision"] == met["decision"]
