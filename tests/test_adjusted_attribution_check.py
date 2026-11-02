"""``scripts/adjusted_attribution_check.py``: construction, units, metrics, merge, calibration and the decision."""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


check = _script("adjusted_attribution_check")
cs = _script("case_study_type_i_check")
d10 = _script("covariate_type_i_check")
SMALL_CV = {"cv1_splits": 2, "cv2_splits": 3, "n_repeats": 2, "max_components": 5, "random_state": 1203}
SIX = ["z", "x1", "x2", "x3", "x4", "x5"]
CELLS = [("anchor", 0.0), ("orientation", 0.5), ("orientation", 1.0)]


def _args(*argv: str, command: str = "run", out_dir: str = "unused"):
    return check.build_parser().parse_args([command, *argv, "--out-dir", out_dir])


def _settings(*argv: str, reps: int = 2, perms: int = 9, bootstrap: int = 5) -> dict:
    settings = check.build_settings(
        _args(*argv, "--reps", str(reps), "--perms", str(perms), "--bootstrap", str(bootstrap))
    )
    settings["cv"] = dict(SMALL_CV)
    return settings


# ── Settings and construction ────────────────────────────────────────────────


def test_cells_and_defaults() -> None:
    settings = check.build_settings(_args())
    assert check.cells(settings) == CELLS
    assert len(check.units(settings)) == 3 * 500
    assert settings["replicates"] == 500 and settings["permutations"] == 999
    assert settings["attribution"]["bootstrap_replicates"] == 200 and settings["attribution"]["top_k"] == 20
    assert settings["base_seed"] == 2033 and settings["injection"]["s_z"] == 10.0
    assert {settings["generator"][f] for f in ("delta_methyl", "delta_expr", "delta_protein")} == {2.0}
    assert settings["generator"]["group_stage_sizes"] == [[11, 10, 28], [9, 10, 12]]
    assert settings["table"]["name"] == "astrocyte" and len(settings["table"]["sha256"]) == 64
    assert settings["covariates"]["declared"] == SIX and settings["min_components"] == 2
    assert settings["orientation_estimator"] == "net_displacement"
    assert settings["residual_correction"] == "df"
    assert settings["layers"] == ["methylation", "expression"]
    assert settings["rule"]["margins"] == check.MARGINS
    assert settings["adopted_calibration"] is None
    for parent_only in ("k", "controls", "calibration"):
        assert parent_only not in settings


def test_seeds_are_fresh_and_matched() -> None:
    settings = _settings()
    seen = []
    for replicate in (0, 1):
        built = [check.build_unit(mode, effect, replicate, settings) for mode, effect in CELLS]
        seeds = built[0]["seeds"]
        assert all(other["seeds"] == seeds for other in built[1:])
        assert set(seeds) == {"generator", "covariates", "rrpp", "noise_covariates", "bootstrap"}
        for name in SIX:
            for other in built[1:]:
                np.testing.assert_array_equal(built[0]["metadata"][name], other["metadata"][name])
        assert {key: seeds[key] for key in ("generator", "covariates", "rrpp")} == d10.replicate_seeds(2033, replicate)
        assert seeds["bootstrap"] == check.bootstrap_seed(2033, replicate)
        for parent in (2026, 2031, 2032):
            family = d10.replicate_seeds(parent, replicate)
            family["noise_covariates"] = cs._seed(parent, "noise_covariates", replicate)
            assert not set(seeds.values()) & set(family.values()), parent
        seen.append(seeds)
    assert not set(seen[0].values()) & set(seen[1].values())


@pytest.mark.parametrize(
    "argv",
    [("--s-z", "15"), ("--base-seed", "2034"), ("--perms", "19"), ("--bootstrap", "7"), ("--top-k", "10")],
)
def test_resume_under_other_settings_is_refused(tmp_path, argv) -> None:
    base = ["--reps", "2", "--perms", "9", "--bootstrap", "5", "--top-k", "20", "--s-z", "10", "--base-seed", "2033"]
    check.ensure_settings(tmp_path, check.build_settings(_args(*base)))
    changed = list(base)
    changed[changed.index(argv[0]) + 1] = argv[1]
    with pytest.raises(SystemExit, match="different settings"):
        check.ensure_settings(tmp_path, check.build_settings(_args(*changed)))
    assert not (tmp_path / "records").exists()


def test_over_headroom_cell_refuses_to_start(tmp_path) -> None:
    settings = _settings()
    settings["cells"].append(["orientation", 3.0])
    with pytest.raises(SystemExit, match=r'\["orientation", 3.0\]'):
        check.run_shard(tmp_path, settings, 0, 1)
    assert check.shard_path(tmp_path, 0, 1).read_text() == ""


# ── Units ────────────────────────────────────────────────────────────────────


def _records(out_dir: Path) -> list[dict]:
    lines = [
        line
        for path in sorted((out_dir / "records").glob("shard_*.jsonl"))
        for line in path.read_text(encoding="utf-8").splitlines()
    ]
    return sorted((json.loads(line) for line in lines), key=lambda record: record["key"])


@pytest.fixture(scope="module")
def sharded_run(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("adjusted-attribution")
    settings = _settings()
    check.ensure_settings(out_dir, settings)
    evaluated = [check.run_shard(out_dir, settings, index, 2) for index in (0, 1)]
    return out_dir, settings, evaluated


def test_two_shards_cover_every_unit_once(sharded_run) -> None:
    out_dir, settings, evaluated = sharded_run
    assert evaluated == [3, 3]
    keys = [record["key"] for record in _records(out_dir)]
    assert sorted(keys) == sorted(check.unit_key(*unit) for unit in check.units(settings))
    assert len(keys) == len(set(keys))


def test_resume_skips_recorded_units(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    assert [check.run_shard(out_dir, settings, index, 2) for index in (0, 1)] == [0, 0]
    assert check.run_shard(out_dir, settings, 0, 3) == 0


def test_records_hold_three_arms_and_the_analysis(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    for record in _records(out_dir):
        assert set(record["arms"]) == set(check.ARMS)
        assert record["arms"]["adjusted"]["mean_source"] == "covariate_adjusted"
        assert record["arms"]["oracle"]["mean_source"] == "arithmetic"
        assert record["arms"]["unadjusted"]["mean_source"] == "arithmetic"
        assert isinstance(record["arms"]["adjusted"]["rank_deficient_replicates"], int)
        assert set(record["arms"]["adjusted"]["blocks"]) == {"0->1", "1->2", "principal"}
        assert set(record["analysis"]) == {"statistics", "p_values", "selected_lv", "cv_mean_auroc"}
        assert set(record["analysis"]["p_values"]) == {"delta", "angle", "shape"}
        assert record["analysis"]["selected_lv"] >= 2
        assert record["truth"]["n_injected"] == round(0.2 * 498)
        if record["mode"] == "anchor":
            assert record["truth"]["n_drivers"] == 0
            assert record["arms"]["adjusted"]["blocks"]["principal"]["observed"]["truth_precision"] is None
        else:
            assert record["truth"]["n_drivers"] > 0
        oracle = record["arms"]["oracle"]["blocks"]["principal"]["observed"]
        assert oracle["cosine_oracle"] == pytest.approx(1.0) and oracle["overlap_oracle"] == 1.0
        assert record["seeds"]["bootstrap"] == check.bootstrap_seed(2033, record["replicate"])
        assert record["settings_hash"] == check.settings_hash(settings)


def test_oracle_removes_the_injection_exactly() -> None:
    from motco.simulations.preprocessing import apply_unit_transform

    settings = _settings()
    built = check.build_unit("orientation", 1.0, 0, settings)
    dataset, result = check.analyse(built, settings)
    oracle, injection = check.oracle_inputs(built, dataset, result, settings)
    features = result.latent_artifacts.features
    z = dataset.metadata["z"].to_numpy()
    np.testing.assert_allclose(features.to_numpy() - oracle.to_numpy(), np.outer(z, injection), atol=1e-10)

    generated = built["dataset"]
    clean = np.hstack(
        [
            apply_unit_transform("methylation", generated.methylation.to_numpy(dtype=float), "beta_to_mvalue"),
            generated.expression.to_numpy(dtype=float),
        ]
    ) / result.latent_artifacts.original_scales
    cells = list(zip(dataset.metadata["group"].astype(str), dataset.metadata["stage"].astype(str)))
    labels = sorted(set(cells))
    oracle_means = np.vstack([oracle.to_numpy()[[cell == label for cell in cells]].mean(axis=0) for label in labels])
    clean_means = np.vstack([clean[[cell == label for cell in cells]].mean(axis=0) for label in labels])
    shift = oracle_means - clean_means
    np.testing.assert_allclose(shift, np.tile(shift[0], (len(labels), 1)), atol=1e-9)


def test_stage_means_mismatch_raises_and_writes_nothing(tmp_path, monkeypatch) -> None:
    settings = _settings()
    check.ensure_settings(tmp_path, settings)
    monkeypatch.setattr(check, "model_transform", lambda model, frame: np.asarray(model.transform(frame.to_numpy())) + 1.0)
    with pytest.raises(check.ConsistencyError, match="stage_means"):
        check.run_shard(tmp_path, settings, 0, 1, limit=1)
    assert check.shard_path(tmp_path, 0, 1).read_text() == ""


# ── Metrics ──────────────────────────────────────────────────────────────────


def _bootstrap(sign, frequency, valid=10):
    return SimpleNamespace(
        valid_replicates=valid,
        sign_stability=np.asarray(sign, dtype=float),
        top_k_selection_frequency=np.asarray(frequency, dtype=float),
    )


def test_contrast_metrics_match_hand_counts() -> None:
    names = [f"f{j}" for j in range(6)]
    contrast = np.array([0.9, -0.8, 0.1, 0.7, 0.0, -0.05])
    oracle = np.array([0.9, 0.0, 0.0, 0.7, 0.6, 0.0])
    metrics = check.contrast_metrics(
        contrast,
        oracle,
        truth=np.array([1, 0, 0, 1, 0, 0], dtype=bool),
        injected=np.array([0, 1, 1, 1, 0, 0], dtype=bool),
        direction=np.array([0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
        bootstrap=_bootstrap([1.0, 0.8, 0.5, 0.9, np.nan, 0.6], [1.0, 0.5, 0.0, 0.7, 0.0, 0.1]),
        feature_names=names,
        top_k=3,
    )
    # top-3 = f0, f1, f3; truth {f0, f3}; f1 is injected and not a driver, f3 is an injected driver.
    assert metrics["truth_precision"] == 1.0  # 2 / min(3, 2)
    assert metrics["contamination"] == pytest.approx(1 / 3)
    assert metrics["projection"] == pytest.approx(-0.8)
    assert metrics["sign_stability"] == pytest.approx(0.9)
    assert metrics["top_k_frequency"] == pytest.approx(2.2 / 3)
    assert metrics["overlap_oracle"] == pytest.approx(2 / 3)  # oracle top-3 = f0, f3, f4
    expected_cosine = contrast @ oracle / (np.linalg.norm(contrast) * np.linalg.norm(oracle))
    assert metrics["cosine_oracle"] == pytest.approx(expected_cosine)
    assert metrics["valid_replicates"] == 10 and metrics["truth_size"] == 2 and metrics["available"] is True


def test_truth_precision_divides_by_the_smaller_of_k_and_the_truth_set() -> None:
    names = [f"f{j}" for j in range(5)]
    metrics = check.contrast_metrics(
        np.array([0.5, 0.4, 0.3, 0.2, 0.1]),
        None,
        truth=np.array([1, 1, 1, 1, 0], dtype=bool),
        injected=np.zeros(5, dtype=bool),
        direction=np.ones(5) / math.sqrt(5),
        bootstrap=_bootstrap(np.ones(5), np.ones(5)),
        feature_names=names,
        top_k=2,
    )
    assert metrics["truth_precision"] == 1.0  # 2 hits / min(2, 4)
    assert metrics["contamination"] == 0.0
    assert metrics["cosine_oracle"] is None and metrics["overlap_oracle"] is None


def test_unavailable_contrast_and_empty_truth() -> None:
    kwargs = dict(
        truth=np.zeros(3, dtype=bool),
        injected=np.zeros(3, dtype=bool),
        direction=np.ones(3) / math.sqrt(3),
        bootstrap=_bootstrap(np.ones(3), np.ones(3), valid=0),
        feature_names=["a", "b", "c"],
        top_k=2,
    )
    unavailable = check.contrast_metrics(None, None, **kwargs)
    assert unavailable["available"] is False and unavailable["contamination"] is None
    assert check.contrast_metrics(np.array([0.3, 0.2, 0.1]), None, **kwargs)["truth_precision"] is None


def test_top_k_ties_break_by_feature_name() -> None:
    assert check.top_k_indices(np.array([0.5, -0.5, 0.2]), ["b", "a", "c"], 1).tolist() == [1]


# ── Merge ────────────────────────────────────────────────────────────────────


def _metric_row(arm: str, replicate: int, **values) -> dict:
    row = {
        "mode": "orientation",
        "effect": 1.0,
        "replicate": replicate,
        "arm": arm,
        "block": "principal",
        "component": "observed",
        "available": True,
        "truth_size": 10,
        "rank_deficient_replicates": None,
        "angle_p_value": 0.5,
        "selected_lv": 2,
    }
    row.update({name: None for name in check.METRICS})
    row.update(values)
    return row


def test_summary_means_and_paired_ses_match_hand_computation() -> None:
    rows = [
        *(_metric_row("adjusted", r, contamination=v) for r, v in enumerate([0.2, 0.3, 0.4])),
        *(_metric_row("oracle", r, contamination=v) for r, v in enumerate([0.1, 0.1, 0.2])),
        _metric_row("unadjusted", 0, contamination=0.9),
    ]
    summary = {(row["arm"], row["quantity"]): row for row in check.summary_rows(rows)}
    adjusted = summary[("adjusted", "contamination")]
    assert adjusted["n"] == 3 and adjusted["mean"] == pytest.approx(0.3)
    assert adjusted["se"] == pytest.approx(np.std([0.2, 0.3, 0.4], ddof=1) / math.sqrt(3))
    paired = summary[("adjusted", "contamination_minus_oracle")]
    differences = [0.1, 0.2, 0.2]
    assert paired["n"] == 3 and paired["mean"] == pytest.approx(np.mean(differences))
    assert paired["se"] == pytest.approx(np.std(differences, ddof=1) / math.sqrt(3))
    single = summary[("unadjusted", "contamination_minus_oracle")]
    assert single["n"] == 1 and single["mean"] == pytest.approx(0.8) and single["se"] is None
    assert ("oracle", "contamination_minus_oracle") not in summary
    assert summary[("adjusted", "truth_precision")]["n"] == 0


def test_merge_writes_the_tables_and_is_idempotent(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    first = check.merge(out_dir, settings)
    assert first["units"] == first["expected"] == 6 and first["missing"] == []
    tables = {name: (out_dir / name).read_bytes() for name in (
        "merged.jsonl", "adjusted_attribution_metrics.csv", "adjusted_attribution_summary.csv"
    )}
    check.merge(out_dir, settings)
    assert {name: (out_dir / name).read_bytes() for name in tables} == tables
    metrics = pd.read_csv(out_dir / "adjusted_attribution_metrics.csv")
    assert len(metrics) == 6 * 3 * 3 * 3
    assert list(metrics.columns) == list(check.METRIC_COLUMNS)


# ── Calibration ──────────────────────────────────────────────────────────────


def _calibration_settings() -> dict:
    return check.calibration_settings(_args("--s-z", "10", command="calibrate"), 10.0)


def _calibration_rows(contamination: dict[str, float], precision: float) -> list[dict]:
    rows = []
    for (mode, effect), gap in contamination.items():
        for replicate in range(4):
            for arm, value in (("unadjusted", 0.1 + gap), ("oracle", 0.1)):
                row = _metric_row(arm, replicate, contamination=value, truth_precision=precision)
                rows.append({**row, "mode": mode, "effect": effect})
    return rows


def test_calibration_settings_run_the_anchor_and_the_largest_orientation() -> None:
    settings = _calibration_settings()
    assert check.cells(settings) == [("anchor", 0.0), ("orientation", 1.0)]
    assert settings["replicates"] == 20 and settings["permutations"] == 99
    assert settings["attribution"]["bootstrap_replicates"] == 50


@pytest.mark.parametrize(
    ("gaps", "precision", "contamination_met", "recoverability_met"),
    [
        ({("anchor", 0.0): 0.3, ("orientation", 1.0): 0.26}, 0.7, True, True),
        ({("anchor", 0.0): 0.3, ("orientation", 1.0): 0.1}, 0.7, False, True),
        ({("anchor", 0.0): 0.3, ("orientation", 1.0): 0.3}, 0.4, True, False),
    ],
    ids=["met", "unmet", "unrecoverable"],
)
def test_calibration_criteria(gaps, precision, contamination_met, recoverability_met) -> None:
    criteria = check.calibration_criteria(_calibration_rows(gaps, precision), _calibration_settings())
    assert criteria["contamination"]["met"] is contamination_met
    assert criteria["recoverability"]["met"] is recoverability_met
    assert criteria["recoverability"]["cell"] == "orientation@1"
    assert criteria["contamination"]["cells"]["anchor"]["mean"] == pytest.approx(0.3)


def _entry(contamination: bool, recoverability: bool) -> dict:
    return {"contamination": {"met": contamination}, "recoverability": {"met": recoverability},
            "met": contamination and recoverability}


def test_calibration_step() -> None:
    assert check.calibration_step(_entry(True, True), 10.0, 50.0) == ("adopt", None)
    assert check.calibration_step(_entry(True, False), 10.0, 50.0) == ("stop_recoverability", None)
    assert check.calibration_step(_entry(False, True), 10.0, 50.0) == ("raise", 15.0)
    assert check.calibration_step(_entry(False, True), 50.0, 50.0) == ("stop_exhausted", None)


def _calibrate(tmp_path: Path, outcomes: list[tuple[bool, bool]]) -> tuple[dict, list[float]]:
    tried: list[float] = []

    def run_once(settings, *, workers):
        tried.append(settings["injection"]["s_z"])
        contamination, recoverability = outcomes[len(tried) - 1]
        return {"s_z": settings["injection"]["s_z"], **_entry(contamination, recoverability)}

    args = _args("--s-z", "10", command="calibrate", out_dir=str(tmp_path))
    return check.calibrate(args, run_once=run_once), tried


def test_calibration_raises_s_z_until_met_and_appends_every_iteration(tmp_path) -> None:
    entry, tried = _calibrate(tmp_path, [(False, True), (False, True), (True, True)])
    assert tried == [10.0, 15.0, 20.0] and entry["outcome"] == "adopt"
    history = json.loads((tmp_path / "calibration.json").read_text())
    assert [item["s_z"] for item in history] == [10.0, 15.0, 20.0]
    assert [item["outcome"] for item in history] == ["raise", "raise", "adopt"]
    _calibrate(tmp_path, [(True, True)])
    assert len(json.loads((tmp_path / "calibration.json").read_text())) == 4


def test_calibration_stops_on_a_recoverability_failure(tmp_path) -> None:
    entry, tried = _calibrate(tmp_path, [(False, True), (True, False)])
    assert tried == [10.0, 15.0] and entry["outcome"] == "stop_recoverability"
    assert [item["outcome"] for item in json.loads((tmp_path / "calibration.json").read_text())] == [
        "raise",
        "stop_recoverability",
    ]


def test_settings_name_the_adopted_iteration(tmp_path) -> None:
    _calibrate(tmp_path, [(False, True), (True, True)])
    path = tmp_path / "calibration.json"
    settings = check.build_settings(_args("--s-z", "15", "--calibration", str(path)))
    assert settings["adopted_calibration"]["iteration"] == 1 and settings["adopted_calibration"]["s_z"] == 15.0
    with pytest.raises(SystemExit, match="no iteration at s_z = 10"):
        check.build_settings(_args("--s-z", "10", "--calibration", str(path)))


# ── Decision ─────────────────────────────────────────────────────────────────


def _summary(overrides: dict | None = None) -> list[dict]:
    """A summary where both gates and every criterion hold; ``overrides`` maps (cell, arm, quantity) → (mean, se)."""

    values: dict[tuple, tuple[float, float | None]] = {}
    for cell in CELLS:
        values[(cell, "unadjusted", "contamination_minus_oracle")] = (0.5, 0.01)
        values[(cell, "unadjusted", "projection_minus_oracle")] = (0.4, 0.01)
        values[(cell, "adjusted", "contamination_minus_oracle")] = (0.0, 0.005)
        values[(cell, "adjusted", "projection_minus_oracle")] = (0.0, 0.001)
        if cell[0] != "anchor":
            for quantity in ("truth_precision", "sign_stability", "top_k_frequency"):
                values[(cell, "adjusted", f"{quantity}_minus_oracle")] = (0.0, 0.005)
    values[(("orientation", 1.0), "oracle", "truth_precision")] = (0.8, 0.01)
    values.update(overrides or {})
    return [
        {"mode": cell[0], "effect": cell[1], "arm": arm, "block": "principal", "component": "observed",
         "quantity": quantity, "n": 500, "mean": mean, "se": se}
        for (cell, arm, quantity), (mean, se) in values.items()
    ]


def _verdict(overrides: dict | None = None) -> dict:
    return check.verdict_from_summary(_summary(overrides), check.build_settings(_args()))


def _reading(verdict: dict, criterion: str, cell: str) -> dict:
    entries = [*verdict["gates"], *verdict["criteria"]]
    return next(entry for entry in entries if entry["criterion"] == criterion and entry["cell"] == cell)


def test_validated_when_both_gates_and_every_criterion_hold() -> None:
    verdict = _verdict()
    assert verdict["verdict"] == "adjusted_attribution_validated" and verdict["failing"] == []
    assert "item 12" in verdict["consequence"] and "0.14.0" in verdict["consequence"]


def test_contamination_bound_fails_at_the_anchor() -> None:
    verdict = _verdict({(("anchor", 0.0), "adjusted", "contamination_minus_oracle"): (0.04, 0.01)})
    reading = _reading(verdict, "contamination", "anchor")
    assert reading["value"] == pytest.approx(0.06) and reading["verdict"] == "fails"
    assert verdict["verdict"] == "diagnose_adjusted_attribution"
    assert verdict["failing"] == [{"criterion": "contamination", "cell": "anchor", "verdict": "fails"}]


def test_residual_confounding_bound_holds() -> None:
    verdict = _verdict(
        {
            (("anchor", 0.0), "unadjusted", "projection_minus_oracle"): (0.40, 0.01),
            (("anchor", 0.0), "adjusted", "projection_minus_oracle"): (0.02, 0.005),
        }
    )
    reading = _reading(verdict, "residual_confounding", "anchor")
    assert reading["value"] == pytest.approx(0.03) and reading["bound"] == pytest.approx(0.04)
    assert reading["verdict"] == "holds"


def test_recovery_bound_holds() -> None:
    verdict = _verdict({(("orientation", 1.0), "adjusted", "truth_precision_minus_oracle"): (-0.03, 0.005)})
    reading = _reading(verdict, "truth_precision", "orientation@1")
    assert reading["value"] == pytest.approx(-0.04) and reading["verdict"] == "holds"
    assert verdict["verdict"] == "adjusted_attribution_validated"


def test_stability_bound_fails_and_names_top_k_frequency() -> None:
    verdict = _verdict({(("orientation", 0.5), "adjusted", "top_k_frequency_minus_oracle"): (-0.10, 0.01)})
    assert verdict["verdict"] == "diagnose_adjusted_attribution"
    assert verdict["failing"] == [{"criterion": "top_k_frequency", "cell": "orientation@0.5", "verdict": "fails"}]


def test_diagnose_names_contamination_at_the_weak_orientation() -> None:
    verdict = _verdict({(("orientation", 0.5), "adjusted", "contamination_minus_oracle"): (0.08, 0.01)})
    assert verdict["verdict"] == "diagnose_adjusted_attribution"
    assert {"criterion": "contamination", "cell": "orientation@0.5", "verdict": "fails"} in verdict["failing"]
    assert "withheld" in verdict["consequence"] and "diagnosis" in verdict["consequence"]


def test_weak_confounding_is_uninformative() -> None:
    verdict = _verdict({(("orientation", 1.0), "unadjusted", "contamination_minus_oracle"): (0.10, 0.01)})
    assert verdict["verdict"] == "construction_uninformative"
    assert verdict["failing"] == [{"criterion": "contamination_gate", "cell": "orientation@1", "verdict": "fails"}]


def test_unrecoverable_drivers_are_uninformative_whatever_the_criteria() -> None:
    verdict = _verdict(
        {
            (("orientation", 1.0), "oracle", "truth_precision"): (0.4, 0.01),
            (("anchor", 0.0), "adjusted", "contamination_minus_oracle"): (0.3, 0.01),
        }
    )
    assert verdict["verdict"] == "construction_uninformative"
    assert verdict["failing"] == [{"criterion": "recoverability_gate", "cell": "orientation@1", "verdict": "fails"}]


def test_unadjusted_readings_are_descriptive() -> None:
    settings = check.build_settings(_args())
    summary = _summary()
    unadjusted = check.criteria_readings(summary, settings, arm="unadjusted")
    assert any(entry["verdict"] == "fails" for entry in unadjusted["criteria"])
    assert check.verdict_from_summary(summary, settings)["verdict"] == "adjusted_attribution_validated"


def test_smoke_decision_is_indicative_and_reports_the_descriptive_readings(sharded_run) -> None:
    out_dir, settings, _ = sharded_run
    check.merge(out_dir, settings)
    decision = check.decide(out_dir, settings)
    assert decision["indicative"] is True
    assert any("permutations" in reason for reason in decision["indicative_reasons"])
    assert any("B = 5" in reason for reason in decision["indicative_reasons"])
    assert decision["verdict"] in check.VERDICTS
    assert decision["settings_hash"] == check.settings_hash(settings)
    assert decision["replicates"] == {"anchor": 2, "orientation@0.5": 2, "orientation@1": 2}
    descriptive = decision["descriptive"]
    assert descriptive["unadjusted_arm"]["arm"] == "unadjusted"
    assert set(descriptive["components"]) == {"pls_captured", "residual"}
    assert set(descriptive["transitions"]) == {"0->1", "1->2"}
    assert set(descriptive["angle_rejected"]["units"]) == {"anchor", "orientation@0.5", "orientation@1"}
    written = json.loads((out_dir / "adjusted_attribution_decision.json").read_text())
    assert written["verdict"] == decision["verdict"]
