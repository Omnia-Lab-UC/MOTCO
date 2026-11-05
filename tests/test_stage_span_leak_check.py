"""``scripts/stage_span_leak_check.py``: the six checkpoints, the merge and the pre-flight gate."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import stage_span_leak_check as leak
finally:
    sys.path.pop(0)

TWO = ["methylation", "expression"]


def _profile(tmp_path: Path) -> Path:
    raw = {
        "generator": {
            "seed": 2,
            "trajectory_mode": "magnitude",
            "n_samples": 60,
            "p_dmp": 0.1,
            "magnitude_kind": "joint",
            "construction": "stage_span",
            "span_layers": TWO,
            "span_scale": "planted",
        },
        "evaluation": {
            "integration_method": "pls",
            "orientation_estimator": "net_displacement",
            "min_components": 2,
            "integration_params": {
                "layers": TWO,
                "cv1_splits": 2,
                "cv2_splits": 2,
                "n_repeats": 1,
                "max_components": 4,
                "random_state": 1203,
            },
            "permutations": 9,
            "seed": 3,
            "span_capture": True,
        },
        "trajectory_modes": ["orientation", "shape", "translation"],
        "effect_sizes": [0.0, 0.5],
        "n_replicates": 4,
        "base_seed": 1028,
        "matched_seeds": {"enabled": True, "primary_family": "design-guidance", "shared_zero_effect_anchor": True},
    }
    path = tmp_path / "profile.json"
    path.write_text(json.dumps(raw))
    return path


@pytest.fixture(scope="module")
def arms(tmp_path_factory) -> dict[str, Path]:
    root = tmp_path_factory.mktemp("leak")
    config = _profile(root)
    out = {}
    for scale in ("planted", "null"):
        run_dir = root / scale
        for shard in (0, 1):
            leak.main(
                [
                    "run",
                    "--config",
                    str(config),
                    "--span-scale",
                    scale,
                    "--modes",
                    "shape",
                    "--effects",
                    "0.5",
                    "--replicates",
                    "1",
                    "--permutations",
                    "9",
                    "--shard",
                    str(shard),
                    "--n-shards",
                    "2",
                    "--out-dir",
                    str(run_dir),
                ]
            )
        leak.main(["merge", "--out-dir", str(run_dir)])
        out[scale] = run_dir
    out["config"] = config
    return out


def test_every_checkpoint_is_reported(arms) -> None:
    units = pd.read_csv(arms["planted"] / "leak_units.csv")
    assert len(units) == 2 and set(units["mode"]) == {"none", "shape"}
    assert units["error"].isna().all()
    for name in leak.CHECKPOINTS:
        for suffix in ("A", "B", "d"):
            assert units[f"{name}_{suffix}"].notna().all(), f"{name}_{suffix}"
    for column in ("rank", "capture", "p_delta", "p_angle", "p_shape", "seed", "seconds"):
        assert column in units
    anchor = units[units["mode"] == "none"].iloc[0]
    for name in leak.CHECKPOINTS[:-1]:
        assert anchor[f"{name}_d"] == pytest.approx(0.0, abs=1e-12)
    summary = json.loads((arms["planted"] / "summary.json").read_text())
    assert summary["n_units"] == summary["n_expected"] == 2
    assert summary["n_failed"] == 0 and summary["n_duplicates"] == 0


def test_the_bend_is_exact_in_its_own_scale(arms) -> None:
    null = pd.read_csv(arms["null"] / "leak_units.csv").set_index("mode").loc["shape"]
    planted = pd.read_csv(arms["planted"] / "leak_units.csv").set_index("mode").loc["shape"]
    assert null["span_scale"] == "null_scale" and planted["span_scale"] == "planted"
    assert abs(null["pop_null_d"]) <= 1e-8 * null["pop_null_A"]
    assert abs(planted["pop_planted_d"]) <= 1e-8 * planted["pop_planted_A"]
    # The null-scale bend leaks under the planted scale; the planted bend under the null scale.
    assert abs(null["pop_planted_d"]) > 1e-4 * null["pop_planted_A"]
    assert planted["scale_iterations"] >= 2
    assert null["seed"] == planted["seed"]


def test_a_resume_under_other_settings_is_refused(arms) -> None:
    with pytest.raises(SystemExit, match="other settings"):
        leak.main(
            [
                "run",
                "--config",
                str(arms["config"]),
                "--span-scale",
                "planted",
                "--modes",
                "shape",
                "--effects",
                "0.25",
                "--replicates",
                "1",
                "--permutations",
                "9",
                "--shard",
                "0",
                "--n-shards",
                "2",
                "--out-dir",
                str(arms["planted"]),
            ]
        )


def test_rates_use_p_strictly_below_alpha() -> None:
    units = pd.DataFrame(
        {
            "column": "n200",
            "mode": "shape",
            "effect": 0.5,
            "replicate": range(4),
            **{f"{name}_{suffix}": 1.0 for name in leak.CHECKPOINTS for suffix in ("A", "B", "d")},
            "capture": 0.9,
            "rank": 3,
            "p_delta": [0.05, 0.049, 0.2, np.nan],
            "p_angle": [0.05, 0.05, 0.05, 0.05],
            "p_shape": [0.01, 0.01, 0.01, 0.01],
            "error": [np.nan, np.nan, np.nan, np.nan],
        }
    )
    [row] = leak.summarize_cells(units).to_dict("records")
    assert (row["n_available_delta"], row["n_rejected_delta"]) == (3, 1)
    assert row["rate_angle"] == 0.0
    assert row["rate_shape"] == 1.0
    assert row["median_pop_lat_share"] == 1.0


# --------------------------------------------------------------------------- #
# The gate on fixture tables
# --------------------------------------------------------------------------- #

COLUMNS = ("n200", "n80", "seaad", "n200_three_blocks", "n200_four_stages")


def _units(errors=(), inexact=()) -> pd.DataFrame:
    rows = []
    for column in COLUMNS:
        for mode, effect in [("none", 0.0)] + [(m, e) for m in ("orientation", "shape") for e in (0.1, 0.25, 0.5)]:
            for replicate in range(3):
                key = (column, mode, effect, replicate)
                row = {"column": column, "mode": mode, "effect": effect, "replicate": replicate}
                if key in errors:
                    row["error"] = "SemiSyntheticTrajectoryError: did not converge"
                else:
                    gap = 1e-6 if key in inexact else 1e-13
                    row.update({"pop_planted_A": 10.0, "pop_planted_d": gap * 10.0, "error": np.nan})
                rows.append(row)
    return pd.DataFrame(rows)


def _cells(*, rejected=None, lat=None, scale="planted") -> pd.DataFrame:
    rejected = rejected or (lambda column, mode, effect, statistic: 5)
    lat = lat or (lambda column: -0.01 if scale == "planted" else -0.5)
    rows = []
    for column in COLUMNS:
        for mode, effect in [("none", 0.0)] + [(m, e) for m in ("orientation", "shape") for e in (0.1, 0.25, 0.5)]:
            row = {"column": column, "mode": mode, "effect": effect, "median_pop_lat_d": lat(column)}
            for statistic in leak.STATISTICS:
                row[f"n_available_{statistic}"] = 100
                row[f"n_rejected_{statistic}"] = rejected(column, mode, effect, statistic)
            rows.append(row)
    return pd.DataFrame(rows)


def test_the_gate_holds_on_clean_tables() -> None:
    result = leak.gate(_units(), _cells(), _cells(scale="null"))
    assert result["verdict"] == "gate_holds"
    assert result["conditions"]["G2"]["m"] == 60
    assert result["conditions"]["G2"]["allowed_exceedances"] == 4
    assert result["conditions"]["G2"]["ceiling_at_median"] == pytest.approx(0.13)


def test_a_failed_unit_fails_g1_alone() -> None:
    result = leak.gate(_units(errors={("n200_four_stages", "orientation", 0.5, 2)}), _cells(), _cells(scale="null"))
    assert result["failing_conditions"] == ["G1"]
    [failure] = result["failures"]
    assert (failure["column"], failure["mode"], failure["status"]) == ("n200_four_stages", "orientation", "failed")


def test_an_inexact_planted_unit_fails_g1_alone() -> None:
    result = leak.gate(_units(inexact={("seaad", "shape", 0.25, 1)}), _cells(), _cells(scale="null"))
    assert result["failing_conditions"] == ["G1"]
    assert result["failures"][0]["status"] == "inexact"


def test_an_off_target_rate_above_the_ceiling_fails_g2_alone() -> None:
    def rejected(column, mode, effect, statistic):
        return 14 if (column, mode, effect, statistic) == ("n80", "shape", 0.5, "delta") else 5

    result = leak.gate(_units(), _cells(rejected=rejected), _cells(scale="null"))
    assert result["failing_conditions"] == ["G2"]
    [failure] = result["failures"]
    assert (failure["column"], failure["trajectory_mode"], failure["statistic"], failure["status"]) == (
        "n80",
        "shape",
        "delta",
        "above_ceiling",
    )


def test_too_many_exceedances_fail_g2_and_are_all_named() -> None:
    exceeding = {(column, "orientation", 0.1, "delta") for column in COLUMNS}

    def rejected(column, mode, effect, statistic):
        return 10 if (column, mode, effect, statistic) in exceeding else 5

    result = leak.gate(_units(), _cells(rejected=rejected), _cells(scale="null"))
    assert result["failing_conditions"] == ["G2"]
    assert {f["column"] for f in result["failures"]} == set(COLUMNS)


def test_a_residual_that_does_not_shrink_fails_g3_alone() -> None:
    def planted_lat(column):
        return -0.6 if column == "n200_four_stages" else -0.01

    result = leak.gate(_units(), _cells(lat=planted_lat), _cells(scale="null"))
    assert result["failing_conditions"] == ["G3"]
    [failure] = result["failures"]
    assert failure["column"] == "n200_four_stages" and failure["status"] == "does_not_shrink"
