"""Study configs, enumeration and reporting for the stage-span construction."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pandas as pd
import pytest

from motco.simulations.grid import (
    SimulationReplicateResult,
    derive_replicate_seed,
    parameter_signature,
    summarize_realized_surgery,
)
from motco.simulations.study import (
    StudyConfig,
    StudyConfigError,
    dump_study_config,
    enumerate_study,
    load_study_config,
)
from motco.simulations.study.report import build_report_frames, write_report_csvs
from motco.simulations.study.spectrum import summarize_span_capture

TWO_BLOCKS = ["methylation", "expression"]
N = "generator.n_samples"
SIZES = "generator.group_stage_sizes"


def _raw(construction: str | None, **overrides) -> dict:
    generator = {"seed": 2, "trajectory_mode": "magnitude", "n_samples": 200, "p_dmp": 0.1, "magnitude_kind": "joint"}
    if construction is not None:
        generator.update({"construction": construction, "span_layers": TWO_BLOCKS})
    raw: dict = {
        "generator": generator,
        "evaluation": {
            "integration_method": "pls",
            "orientation_estimator": "net_displacement",
            "min_components": 2,
            "integration_params": {"layers": TWO_BLOCKS},
            "permutations": 19,
            "seed": 3,
            "span_capture": True,
        },
        "trajectory_modes": ["orientation", "shape", "translation"],
        "effect_sizes": [0.0, 0.25, 1.0],
        "n_replicates": 2,
        "base_seed": 1028,
        "matched_seeds": {"enabled": True, "primary_family": "design-guidance", "shared_zero_effect_anchor": True},
        "design_grid": {
            "points": [
                {N: 200, SIZES: None},
                {N: 80, SIZES: None},
                {N: 600, SIZES: [[11, 10, 28], [9, 10, 12]]},
            ]
        },
    }
    raw.update(overrides)
    return raw


def _load(tmp_path: Path, raw: dict, name: str = "config.json") -> StudyConfig:
    path = tmp_path / name
    path.write_text(json.dumps(raw), encoding="utf-8")
    return load_study_config(path)


# --------------------------------------------------------------------------- #
# 4.1 Loading
# --------------------------------------------------------------------------- #


def test_stage_span_config_loads_and_reloads_identically(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw("stage_span"))
    assert config.generator.construction == "stage_span"
    assert config.generator.span_layers == tuple(TWO_BLOCKS)
    assert config.evaluation.span_capture is True
    first = enumerate_study(config).cells
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    again = enumerate_study(load_study_config(dumped)).cells
    reloaded = enumerate_study(_load(tmp_path, _raw("stage_span"), "again.json")).cells
    for other in (again, reloaded):
        assert [cell.cell_id for cell in other] == [cell.cell_id for cell in first]
        assert [parameter_signature(cell) for cell in other] == [parameter_signature(cell) for cell in first]
        assert [derive_replicate_seed(cell, 1) for cell in other] == [derive_replicate_seed(cell, 1) for cell in first]


def test_dump_omits_unset_construction_settings(tmp_path: Path) -> None:
    config = _load(tmp_path, _raw(None, evaluation={**_raw(None)["evaluation"], "span_capture": False}))
    dumped = tmp_path / "dumped.json"
    dump_study_config(config, dumped)
    payload = json.loads(dumped.read_text())
    assert not {"construction", "span_layers"} & set(payload["generator"])
    assert "span_capture" not in payload["evaluation"]


def _measure_all_blocks_at_second_point(raw: dict) -> None:
    for index, point in enumerate(raw["design_grid"]["points"]):
        point["evaluation.integration_params.layers"] = None if index == 1 else TWO_BLOCKS


@pytest.mark.parametrize(
    ("mutate", "column", "span", "measured"),
    [
        (
            lambda raw: raw["evaluation"]["integration_params"].pop("layers"),
            "baseline",
            "['methylation', 'expression']",
            "['methylation', 'expression', 'proteomics']",
        ),
        (
            lambda raw: _measure_all_blocks_at_second_point(raw),
            "design point",
            "['methylation', 'expression']",
            "['methylation', 'expression', 'proteomics']",
        ),
        (
            lambda raw: raw.update({"axes": {"generator.span_layers": [["methylation"]]}}),
            "axes['generator.span_layers']",
            "['methylation']",
            "['methylation', 'expression']",
        ),
    ],
)
def test_span_layer_mismatch_is_rejected_naming_column_and_blocks(tmp_path, mutate, column, span, measured) -> None:
    raw = _raw("stage_span")
    mutate(raw)
    with pytest.raises(StudyConfigError) as excinfo:
        _load(tmp_path, raw)
    message = str(excinfo.value)
    assert column in message
    assert span in message
    assert measured in message


def test_construction_design_axis_is_rejected(tmp_path: Path) -> None:
    raw = _raw(None)
    for point in raw["design_grid"]["points"]:
        point["generator.construction"] = None
    with pytest.raises(StudyConfigError, match="sibling profiles"):
        _load(tmp_path, raw)
    raw = _raw(None, axes={"generator.construction": ["stage_span"]})
    with pytest.raises(StudyConfigError, match="sibling profiles"):
        _load(tmp_path, raw)


@pytest.mark.parametrize(
    ("generator", "match"),
    [
        ({"construction": "lattice"}, "Unknown construction"),
        ({"span_layers": TWO_BLOCKS}, "construction is not 'stage_span'"),
        ({"construction": "stage_span", "span_layers": ["lipids"]}, "unknown layer"),
    ],
)
def test_malformed_construction_settings_are_rejected_at_load(tmp_path, generator, match) -> None:
    raw = _raw(None)
    raw["generator"].update(generator)
    with pytest.raises(StudyConfigError, match=match):
        _load(tmp_path, raw)


# --------------------------------------------------------------------------- #
# 4.2 Enumeration: shared anchors and no headroom
# --------------------------------------------------------------------------- #


def _anchors(config: StudyConfig) -> dict:
    return {
        json.dumps(dict(cell.metadata.get("design_point") or {}), sort_keys=True): cell
        for cell in enumerate_study(config).cells
        if cell.metadata.get("zero_effect_anchor")
    }


def test_stage_span_and_surgery_siblings_share_their_anchors(tmp_path: Path) -> None:
    span = _anchors(_load(tmp_path, _raw("stage_span"), "span.json"))
    surgery = _anchors(_load(tmp_path, _raw(None), "surgery.json"))
    assert len(span) == 3
    assert span.keys() == surgery.keys()
    for point, cell in span.items():
        other = surgery[point]
        assert cell.cell_id == other.cell_id
        assert parameter_signature(cell) == parameter_signature(other)
        assert [derive_replicate_seed(cell, i) for i in range(2)] == [derive_replicate_seed(other, i) for i in range(2)]
        assert cell.generator_params.construction is None
        assert cell.generator_params.span_layers is None


def test_type_i_baseline_is_shared_and_effect_cells_differ(tmp_path: Path) -> None:
    span = enumerate_study(_load(tmp_path, _raw("stage_span"), "span.json")).cells
    surgery = enumerate_study(_load(tmp_path, _raw(None), "surgery.json")).cells
    null_span = [cell for cell in span if cell.phase == "type_i_baseline" and "negative_control" not in cell.metadata]
    null_surgery = [
        cell for cell in surgery if cell.phase == "type_i_baseline" and "negative_control" not in cell.metadata
    ]
    assert [cell.cell_id for cell in null_span] == [cell.cell_id for cell in null_surgery]
    effect_span = {cell.cell_id for cell in span if not cell.metadata.get("zero_effect_anchor")}
    effect_surgery = {cell.cell_id for cell in surgery if not cell.metadata.get("zero_effect_anchor")}
    assert effect_span.isdisjoint(effect_surgery - {cell.cell_id for cell in null_surgery})
    stage_span_cells = [cell for cell in span if cell.generator_params.construction == "stage_span"]
    assert stage_span_cells
    assert all(cell.generator_params.span_layers == tuple(TWO_BLOCKS) for cell in stage_span_cells)


def test_large_stage_span_effects_enumerate_without_headroom_failure(tmp_path: Path) -> None:
    raw = _raw(None, effect_sizes=[0.0, 2.0])
    raw["generator"]["p_dmp"] = 0.3
    with pytest.raises(StudyConfigError, match="surgery than the expected pool headroom"):
        enumerate_study(_load(tmp_path, raw, "surgery.json"))
    raw = _raw("stage_span", effect_sizes=[0.0, 2.0])
    raw["generator"]["p_dmp"] = 0.3
    cells = enumerate_study(_load(tmp_path, raw, "span.json")).cells
    assert any(cell.generator_params.group_effect_size == 2.0 for cell in cells)


def test_infeasible_stage_span_effect_fails_at_enumeration(tmp_path: Path) -> None:
    raw = _raw("stage_span", effect_sizes=[0.0, 2.5])
    config = _load(tmp_path, raw)
    with pytest.raises(ValueError, match="at most 2"):
        enumerate_study(config)


def _record(cell_id: str, *, construction=None, transform=None, **fields) -> SimulationReplicateResult:
    truth = {"transform": transform or {}}
    if construction:
        truth["construction"] = construction
    defaults = {
        "cell_id": cell_id,
        "phase": "power_primary",
        "replicate_index": 0,
        "replicate_seed": 1,
        "generator_seed": 1,
        "evaluation_seed": 3,
        "parameter_signature": "sig",
        "status": "completed",
        "truth_metadata": truth,
    }
    defaults.update(fields)
    return SimulationReplicateResult(**defaults)


def test_realized_surgery_reports_no_surgery_for_stage_span_cells() -> None:
    records = [
        _record("surgery", cell_metadata={"trajectory_mode": "orientation", "effect_size": 0.5},
                transform={"orientation_relocated": 9, "orientation_nominal": 9, "censored": False}),
        _record("span", construction="stage_span", cell_metadata={"trajectory_mode": "orientation", "effect_size": 0.5},
                transform={"construction": "stage_span", "rotation_degrees": 45.0, "orientation_relocated": 1}),
    ]
    summaries = {summary.cell_id: summary for summary in summarize_realized_surgery(records)}
    assert summaries["surgery"].realized_mean == 9
    assert summaries["span"].surgery_replicates == 0
    assert summaries["span"].realized_mean is None


# --------------------------------------------------------------------------- #
# 4.3 span_capture.csv
# --------------------------------------------------------------------------- #


def _capture_record(cell_id, replicate, *, mode, effect, gd, ps=0.95, rank=3, latent=None, intended=None,
                    construction=None, point=None) -> SimulationReplicateResult:
    checkpoints = {}
    if latent is not None:
        checkpoints["population_latent"] = {"joint": {"delta": latent, "angle": latent, "shape": latent}}
    truth: dict = {"transform": {}}
    if construction:
        truth["construction"] = construction
    if intended is not None:
        truth["stage_span"] = {"intended_geometry": {"delta": intended, "angle": intended, "shape": intended}}
    metadata = {"trajectory_mode": mode, "effect_size": effect}
    if point is not None:
        metadata["design_point"] = point
    return SimulationReplicateResult(
        cell_id=cell_id,
        phase="power_primary",
        replicate_index=replicate,
        replicate_seed=replicate,
        generator_seed=replicate,
        evaluation_seed=3,
        parameter_signature="sig",
        status="completed",
        p_values={"delta": 0.5, "angle": 0.5, "shape": 0.5},
        pair_statistics={"delta": 1.0, "angle": 1.0, "shape": 1.0},
        realized_geometry={
            "schema_version": 1,
            "checkpoints": checkpoints,
            "span_capture": {
                "group_difference": gd,
                "pooled_stage": ps,
                "rank": rank,
                "n_features": 498,
                "random_share": rank / 498,
            },
        },
        truth_metadata=truth,
        cell_metadata=metadata,
    )


def test_span_capture_table_summarizes_each_cell() -> None:
    point = {N: 200}
    records = [
        *[
            _capture_record("orient", i, mode="orientation", effect=0.5, gd=0.9 + 0.01 * i, latent=40.0 + i,
                            intended=45.0, construction="stage_span", point=point)
            for i in range(5)
        ],
        *[_capture_record("anchor", i, mode="none", effect=0.0, gd=None, point=point) for i in range(3)],
    ]
    frame = summarize_span_capture(records)
    assert list(frame.columns[:1]) == [N]
    assert list(frame["cell_id"]) == ["anchor", "orient"]
    orient = frame.set_index("cell_id").loc["orient"]
    assert orient["construction"] == "stage_span"
    assert orient["n_group_difference_available"] == 5
    assert orient["median_group_difference"] == pytest.approx(0.92)
    assert orient["q10_group_difference"] == pytest.approx(0.904)
    assert orient["target_statistic"] == "angle"
    assert orient["median_population_latent_target"] == pytest.approx(42.0)
    assert orient["median_intended_target"] == pytest.approx(45.0)
    assert orient["median_rank"] == 3
    anchor = frame.set_index("cell_id").loc["anchor"]
    assert anchor["construction"] == "surgery"
    assert anchor["n_group_difference_available"] == 0
    assert pd.isna(anchor["median_group_difference"])
    assert anchor["n_pooled_stage_available"] == 3
    assert pd.isna(anchor["target_statistic"])


def test_span_capture_csv_is_written_only_with_span_capture(tmp_path: Path) -> None:
    records = [_capture_record("orient", i, mode="orientation", effect=0.5, gd=0.9) for i in range(3)]
    paths = write_report_csvs(build_report_frames([], [], records), tmp_path / "with")
    assert (tmp_path / "with" / "span_capture.csv").exists()
    assert paths["span_capture"] == tmp_path / "with" / "span_capture.csv"

    plain = [copy.deepcopy(record) for record in records]
    for record in plain:
        record.realized_geometry.pop("span_capture")
    paths = write_report_csvs(build_report_frames([], [], plain), tmp_path / "without")
    assert "span_capture" not in paths
    assert not (tmp_path / "without" / "span_capture.csv").exists()
    for name in ("specificity_matrix", "power_curves", "type_i_table", "realized_surgery"):
        assert (tmp_path / "with" / f"{name}.csv").read_bytes() == (tmp_path / "without" / f"{name}.csv").read_bytes()
