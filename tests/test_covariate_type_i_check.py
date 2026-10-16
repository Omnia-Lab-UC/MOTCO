"""``scripts/covariate_type_i_check.py``: reproducibility, resume, confounding and the decision rule."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest


def _script():
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__("covariate_type_i_check")
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


def test_resume_skips_recorded_replicates(runs) -> None:
    out_dir, settings = runs[0]
    before = _records(out_dir)
    assert check.run_shard(out_dir, settings, 0, 1) == 0
    assert _records(out_dir) == before


def test_arms_share_the_latent_space_and_differ_in_the_model(runs) -> None:
    for record in _records(runs[0][0]):
        arms = record["arms"]
        assert arms["adjusted"]["selected_lv"] == arms["unadjusted"]["selected_lv"]
        assert arms["adjusted"]["has_covariate_block"] is True
        assert arms["unadjusted"]["has_covariate_block"] is False
        assert arms["adjusted"]["covariate_design_columns"] == ["z", "category[c1]", "category[c2]"]
        assert set(arms["adjusted"]["p_values"]) == {"delta", "angle", "shape"}


def test_covariates_are_confounded_by_construction(runs) -> None:
    for record in _records(runs[0][0]):
        cells = record["covariate_cells"]
        means = [value for stages in cells["numeric_means"].values() for value in stages.values()]
        assert len(set(means)) == len(means)
        assert abs(cells["numeric_interaction_contrast"]) > 0.0


def test_merge_writes_the_decision(runs) -> None:
    out_dir, settings = runs[0]
    outcome = check.merge(out_dir, settings)
    assert outcome["replicates"] == 2 and outcome["missing"] == []
    decision = json.loads((out_dir / "covariate_type_i_decision.json").read_text())
    assert decision["decision"] in {"covariates_validated", "diagnose_covariate_excess"}
    assert decision["construction"]["injection"]["s_z"] == 3.0
    assert (out_dir / "covariate_type_i_rates.csv").exists()


def test_changed_settings_are_refused(runs) -> None:
    out_dir, settings = runs[0]
    with pytest.raises(SystemExit, match="differing: permutations"):
        check.ensure_settings(out_dir, {**settings, "permutations": 99})


# ── Decision rule (synthetic p-value tables) ─────────────────────────────────


def _rows(adjusted: dict[str, list[float | None]], unadjusted: dict[str, list[float | None]]) -> list[dict]:
    rows = []
    for arm, table in (("adjusted", adjusted), ("unadjusted", unadjusted)):
        for statistic, p_values in table.items():
            rows += [
                {"replicate": i, "arm": arm, "statistic": statistic, "p_value": p} for i, p in enumerate(p_values)
            ]
    return rows


def _uniform(n: int, rejections: int) -> list[float | None]:
    return [0.01] * rejections + [0.5] * (n - rejections)


def test_adjusted_arm_within_bound_validates() -> None:
    adjusted = {s: _uniform(500, 30) for s in check.STATISTICS}  # 0.060 <= 0.0695
    unadjusted = {s: _uniform(500, 400) for s in check.STATISTICS}
    decision = check.covariate_type_i_decision(_rows(adjusted, unadjusted), alpha=0.05, k=2.0)
    assert decision["decision"] == "covariates_validated"
    assert decision["bound"]["angle"] == pytest.approx(0.0695, abs=5e-5)
    assert decision["replicates"] == {s: 500 for s in check.STATISTICS}
    assert decision["unadjusted"]["delta"] == 0.8


def test_excess_in_any_statistic_blocks_validation() -> None:
    adjusted = {"delta": _uniform(500, 25), "angle": _uniform(500, 40), "shape": _uniform(500, 25)}
    unadjusted = {s: _uniform(500, 25) for s in check.STATISTICS}
    decision = check.covariate_type_i_decision(_rows(adjusted, unadjusted), alpha=0.05, k=2.0)
    assert decision["decision"] == "diagnose_covariate_excess"
    assert decision["verdicts"] == {"delta": "holds", "angle": "exceeds", "shape": "holds"}


def test_unadjusted_arm_does_not_enter_the_decision() -> None:
    adjusted = {s: _uniform(500, 25) for s in check.STATISTICS}
    clean = check.covariate_type_i_decision(_rows(adjusted, {s: _uniform(500, 25) for s in check.STATISTICS}),
                                            alpha=0.05, k=2.0)
    inflated = check.covariate_type_i_decision(_rows(adjusted, {s: _uniform(500, 500) for s in check.STATISTICS}),
                                               alpha=0.05, k=2.0)
    assert clean["decision"] == inflated["decision"] == "covariates_validated"
    assert clean["verdicts"] == inflated["verdicts"]
    assert inflated["unadjusted"]["angle"] == 1.0


def test_rejection_is_p_at_or_below_alpha_and_missing_p_values_are_unavailable() -> None:
    rows = _rows({"delta": [0.05, 0.051], "angle": [None, None], "shape": [0.5, 0.5]}, {})
    rates = {(r["arm"], r["statistic"]): r for r in check.rejection_rates(rows, alpha=0.05, k=2.0)}
    assert rates[("adjusted", "delta")]["rejections"] == 1
    assert rates[("adjusted", "angle")]["verdict"] == "unavailable"
    assert rates[("unadjusted", "delta")]["verdict"] == "descriptive"


def test_resume_skips_replicates_recorded_under_another_shard_layout(runs) -> None:
    out_dir, settings = runs[0]
    # Both replicates live in shard_000_of_001; a 2-shard layout has nothing left to do.
    assert check.run_shard(out_dir, settings, 0, 2) == 0
    assert check.run_shard(out_dir, settings, 1, 2) == 0
    assert check.merge(out_dir, settings)["replicates"] == 2
