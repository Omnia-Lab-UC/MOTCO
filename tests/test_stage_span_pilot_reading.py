"""``scripts/stage_span_pilot_reading.py`` on synthetic pilot reports."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import stage_span_pilot_reading as reading
finally:
    sys.path.pop(0)

N = "generator.n_samples"
SIZES = "generator.group_stage_sizes"
COLUMNS = [(200, None), (80, None), (600, "((11, 10, 28), (9, 10, 12))")]
EFFECTS = {"magnitude": (0.02, 0.05, 0.1), "geometry": (0.25, 0.5, 1.0)}
MODES = {"magnitude": ("magnitude",), "geometry": ("orientation", "shape", "translation")}
TARGET = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}


def _write_report(
    root: Path,
    construction: str,
    family: str,
    *,
    rate=lambda column, mode, effect, statistic: 0.04,
    capture=lambda column, mode, effect: 0.95,
) -> None:
    report = root / f"{construction}_{family}" / "report"
    report.mkdir(parents=True)
    operating, shares = [], []
    for n, sizes in COLUMNS:
        cells = [("none", 0.0)] + [(mode, effect) for mode in MODES[family] for effect in EFFECTS[family]]
        for mode, effect in cells:
            for statistic in reading.STATISTICS:
                value = rate((n, sizes), mode, effect, statistic)
                operating.append(
                    {
                        N: n,
                        SIZES: sizes,
                        "phase": "power_primary" if n == 200 else "power_design",
                        "trajectory_mode": mode,
                        "effect_size": effect,
                        "statistic": statistic,
                        "n_available": 100,
                        "rejection_rate": value,
                        "monte_carlo_se": (value * (1 - value) / 100) ** 0.5,
                    }
                )
            constant = mode in ("none", "translation")
            shares.append(
                {
                    N: n,
                    SIZES: sizes,
                    "phase": "power_primary" if n == 200 else "power_design",
                    "cell_id": f"{construction}-{family}-{n}-{mode}-{effect}",
                    "trajectory_mode": mode,
                    "effect_size": effect,
                    "construction": "stage_span" if construction == "span" and mode != "none" else "surgery",
                    "n_replicates": 100,
                    "n_group_difference_available": 0 if constant else 100,
                    "median_group_difference": None if constant else capture((n, sizes), mode, effect),
                    "median_pooled_stage": 0.97,
                    "median_rank": 3,
                    "median_population_latent_target": None,
                    "median_intended_target": None,
                }
            )
    # A Type I row in the capture table carries no design point and is ignored.
    shares.append({**shares[0], N: None, SIZES: None, "phase": "type_i_baseline", "cell_id": "type-i"})
    pd.DataFrame(operating).to_csv(report / "design_point_operating.csv", index=False)
    pd.DataFrame(shares).to_csv(report / "span_capture.csv", index=False)


def _pilot(tmp_path: Path, **span_overrides) -> Path:
    for construction in reading.CONSTRUCTIONS:
        for family in reading.FAMILIES:
            overrides = span_overrides if construction == "span" else {}
            if construction == "surgery":
                overrides = {
                    "rate": lambda column, mode, effect, statistic: 0.3 if mode == "magnitude" else 0.04,
                    "capture": lambda column, mode, effect: 0.02,
                }
            _write_report(tmp_path, construction, family, **overrides)
    return tmp_path


def _powerful(column, mode, effect, statistic):
    return 0.9 if TARGET.get(mode) == statistic else 0.04


def test_bound_is_the_k2_size_bound() -> None:
    assert reading.bound(100) == pytest.approx(0.05 + 2 * (0.05 * 0.95 / 100) ** 0.5)
    assert reading.bound(100) == pytest.approx(0.0936, abs=1e-4)
    assert reading.bound(0) is None


def test_all_criteria_hold_proposes_paper_grade(tmp_path: Path) -> None:
    result = reading.run(_pilot(tmp_path, rate=_powerful))
    assert result["verdict"] == "propose_paper_grade"
    assert result["failures"] == []
    assert all(result["criteria"][name]["holds"] for name in ("R1", "R2", "R3"))
    written = json.loads((tmp_path / reading.OUT_NAME).read_text())
    assert written["verdict"] == "propose_paper_grade"
    assert len(written["inputs"]) == 8
    # R1 reads 3 columns x 3 effects x 3 statistics of the stage-span translation.
    assert len(result["criteria"]["R1"]["observations"]) == 27
    # R2: 3 columns x (3 magnitude + 6 orientation/shape) cells; R3: two off-target statistics each.
    assert len(result["criteria"]["R2"]["observations"]) == 27
    assert len(result["criteria"]["R3"]["observations"]) == 54


def test_r1_failure_names_the_translation_cell_and_statistic(tmp_path: Path) -> None:
    def rate(column, mode, effect, statistic):
        if mode == "translation" and column[0] == 80 and effect == 1.0 and statistic == "angle":
            return 0.12
        return _powerful(column, mode, effect, statistic)

    result = reading.run(_pilot(tmp_path, rate=rate))
    assert result["verdict"] == "revise_construction"
    assert [f["criterion"] for f in result["failures"]] == ["R1"]
    failure = result["failures"][0]
    assert failure["column"] == f"{N}=80; {SIZES}=null"
    assert (failure["trajectory_mode"], failure["effect_size"], failure["statistic"]) == ("translation", 1.0, "angle")
    assert failure["status"] == "exceeds" and failure["rate"] == 0.12
    assert "R1" in result["rationale"]


def test_r2_failure_names_the_uncaptured_cell(tmp_path: Path) -> None:
    def capture(column, mode, effect):
        return 0.6 if column[0] == 600 and mode == "shape" else 0.95

    result = reading.run(_pilot(tmp_path, rate=_powerful, capture=capture))
    assert result["verdict"] == "revise_construction"
    failures = result["failures"]
    assert {f["criterion"] for f in failures} == {"R2"}
    assert len(failures) == 3  # the three shape effects at the SEA-AD column
    assert all(f["trajectory_mode"] == "shape" and f["status"] == "below_floor" for f in failures)
    assert all(f["column"].startswith(f"{N}=600") for f in failures)


def test_r3_failure_names_the_off_target_statistic(tmp_path: Path) -> None:
    def rate(column, mode, effect, statistic):
        if mode == "orientation" and statistic == "shape" and effect == 0.5:
            return 0.2
        return _powerful(column, mode, effect, statistic)

    result = reading.run(_pilot(tmp_path, rate=rate))
    failures = result["failures"]
    assert result["verdict"] == "revise_construction"
    assert {f["criterion"] for f in failures} == {"R3"}
    assert len(failures) == 3  # one per column
    assert all((f["trajectory_mode"], f["statistic"], f["effect_size"]) == ("orientation", "shape", 0.5)
               for f in failures)


def test_target_power_is_not_judged(tmp_path: Path) -> None:
    def weak(column, mode, effect, statistic):
        return 0.04

    result = reading.run(_pilot(tmp_path, rate=weak))
    assert result["verdict"] == "propose_paper_grade"


def test_descriptive_tables(tmp_path: Path) -> None:
    result = reading.run(_pilot(tmp_path, rate=_powerful))
    rows = {(r["family"], r["column"], r["trajectory_mode"], r["effect_size"]): r
            for r in result["target_power_and_capture"]}
    shape = rows[("geometry", f"{N}=200; {SIZES}=null", "shape", 1.0)]
    assert shape["span_target_rate"] == 0.9 and shape["surgery_target_rate"] == 0.04
    assert shape["span_median_group_difference"] == 0.95 and shape["surgery_median_group_difference"] == 0.02
    translation = rows[("geometry", f"{N}=80; {SIZES}=null", "translation", 0.25)]
    assert translation["target_statistic"] is None and translation["span_angle_rate"] == 0.04
    differences = result["magnitude_delta_difference"]
    assert len(differences) == 9
    first = differences[0]
    assert first["difference"] == pytest.approx(0.9 - 0.3)
    assert first["pooled_se"] == pytest.approx(((0.9 * 0.1 + 0.3 * 0.7) / 100) ** 0.5)
    anchors = result["anchor_size"]
    assert len(anchors) == 4 * 3
    assert all(row["delta_bound"] == pytest.approx(0.0936, abs=1e-4) for row in anchors)


def test_missing_report_or_table_is_refused(tmp_path: Path) -> None:
    _write_report(tmp_path, "span", "geometry")
    with pytest.raises(reading.ReadingError, match="design_point_operating|span_capture|Missing"):
        reading.run(tmp_path)
    (tmp_path / "span_geometry" / "report" / "span_capture.csv").unlink()
    with pytest.raises(reading.ReadingError, match="span_capture.csv"):
        reading.load_report(tmp_path / "span_geometry" / "report")


def test_missing_capture_share_is_an_unavailable_r2_failure(tmp_path: Path) -> None:
    result_dir = _pilot(tmp_path, rate=_powerful)
    path = result_dir / "span_magnitude" / "report" / "span_capture.csv"
    frame = pd.read_csv(path)
    smallest = (frame["trajectory_mode"] == "magnitude") & (frame["effect_size"] == 0.02)
    frame.loc[smallest, "median_group_difference"] = None
    frame.to_csv(path, index=False)
    result = reading.run(result_dir)
    assert result["verdict"] == "revise_construction"
    assert {f["status"] for f in result["failures"]} == {"unavailable"}
    assert {f["criterion"] for f in result["failures"]} == {"R2"}
