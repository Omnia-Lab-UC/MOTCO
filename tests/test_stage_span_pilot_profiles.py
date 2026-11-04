"""The four stage-span pilot profiles: structure, siblings, shared anchors and design-guidance reproduction.

``stage_span_pilot_{surgery,span}_{magnitude,geometry}.json`` pair the feature-set
surgeries with the stage-span construction at three columns (n = 200, n = 80, the
SEA-AD table). The siblings differ only in ``construction``/``span_layers``, share
their zero-effect anchors cell for cell, and the surgery siblings' n = 200 and
n = 80 cells draw the design-guidance pilot cells' generator seeds and parameters.
Identity (cell ids, signatures, seeds) is pinned by
``tests/data/post_floor_config_identity.json``.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from motco.simulations.grid import SimulationCell, derive_replicate_seed, parameter_signature
from motco.simulations.study import StudyConfig, enumerate_study, load_study_config
from motco.simulations.study.enumerate import DESIGN_PHASE, DESIGN_POINT_KEY

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
TWO = ["methylation", "expression"]
N = "generator.n_samples"
SIZES = "generator.group_stage_sizes"
SEA_AD = [[11, 10, 28], [9, 10, 12]]

PROFILES = {
    (construction, family): f"stage_span_pilot_{construction}_{family}.json"
    for construction in ("surgery", "span")
    for family in ("magnitude", "geometry")
}
DESIGN_GUIDANCE = {
    "magnitude": "design_guidance_core_magnitude_pilot.json",
    "geometry": "design_guidance_core_geometry_pilot.json",
}
EXPECTED_POINTS = [{N: 200, SIZES: None}, {N: 80, SIZES: None}, {N: 600, SIZES: SEA_AD}]


def _load(name: str) -> StudyConfig:
    return load_study_config(CONFIG_DIR / name)


def _raw(name: str) -> dict[str, Any]:
    return json.loads((CONFIG_DIR / name).read_text(encoding="utf-8"))


def _frozen(value: Any) -> Any:
    return tuple(_frozen(item) for item in value) if isinstance(value, list | tuple) else value


def _columns(config: StudyConfig) -> dict[tuple[tuple[str, Any], ...], list[SimulationCell]]:
    columns: dict[tuple[tuple[str, Any], ...], list[SimulationCell]] = {}
    for cell in enumerate_study(config).cells:
        if cell.phase not in ("power_primary", DESIGN_PHASE):
            continue
        point = cell.metadata[DESIGN_POINT_KEY]
        columns.setdefault(tuple((axis, _frozen(value)) for axis, value in point.items()), []).append(cell)
    return columns


def _mode_effect(cell: SimulationCell) -> tuple[str, float]:
    return str(cell.metadata["trajectory_mode"]), float(cell.metadata["effect_size"])


def _seeds(cell: SimulationCell, n: int) -> list[int]:
    return [derive_replicate_seed(cell, i) for i in range(n)]


# ── Structure ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize(("construction", "family"), sorted(PROFILES))
def test_profile_runs_the_production_measurement_with_span_capture(construction: str, family: str) -> None:
    config = _load(PROFILES[(construction, family)])
    assert config.n_replicates == 100 and config.evaluation.permutations == 199
    assert config.base_seed == 1028 and config.matched_seeds.primary_family == "design-guidance"
    assert config.matched_seeds.enabled and config.matched_seeds.shared_zero_effect_anchor
    assert config.evaluation.integration_method == "pls"
    assert config.evaluation.orientation_estimator == "net_displacement"
    assert config.evaluation.min_components == 2
    assert config.evaluation.span_capture is True
    assert config.evaluation.integration_params["layers"] == TWO
    assert config.generator.magnitude_kind == "joint"
    assert (config.generator.n_samples, config.generator.n_stages, config.generator.p_dmp) == (200, 3, 0.1)
    if construction == "span":
        assert config.generator.construction == "stage_span"
        assert config.generator.span_layers == tuple(TWO)
    else:
        assert config.generator.construction is None and config.generator.span_layers is None


@pytest.mark.parametrize(("construction", "family"), sorted(PROFILES))
def test_profile_enumerates_three_columns_with_an_anchor_and_every_mode_effect(construction: str, family: str) -> None:
    config = _load(PROFILES[(construction, family)])
    assert [
        {axis: (_frozen(value) if value is not None else None) for axis, value in point.items()}
        for point in config.design_points()
    ] == [{axis: _frozen(value) for axis, value in point.items()} for point in EXPECTED_POINTS]
    if family == "magnitude":
        expected = {("none", 0.0)} | {("magnitude", e) for e in (0.02, 0.05, 0.1)}
    else:
        expected = {("none", 0.0)} | {(m, e) for m in ("orientation", "shape", "translation") for e in (0.25, 0.5, 1.0)}
    columns = _columns(config)
    assert len(columns) == 3
    for cells in columns.values():
        assert {_mode_effect(cell) for cell in cells} == expected
        assert len(cells) == len(expected)
        for cell in cells:
            assert cell.generator_params.n_stages == 3
    cells = enumerate_study(config).cells
    assert len(cells) == {"magnitude": 14, "geometry": 32}[family]


def test_sea_ad_column_keeps_the_default_sample_size() -> None:
    config = _load(PROFILES[("span", "geometry")])
    for key, cells in _columns(config).items():
        point = dict(key)
        if point[SIZES] is None:
            continue
        for cell in cells:
            assert cell.generator_params.group_stage_sizes == _frozen(SEA_AD)
            assert cell.generator_params.n_samples == 600


# ── Siblings ─────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("family", ["magnitude", "geometry"])
def test_siblings_differ_only_in_construction_and_span_layers(family: str) -> None:
    surgery, span = _raw(PROFILES[("surgery", family)]), _raw(PROFILES[("span", family)])
    span_generator = dict(span["generator"])
    assert span_generator.pop("construction") == "stage_span"
    assert span_generator.pop("span_layers") == TWO
    assert span_generator == surgery["generator"]
    for key in set(surgery) | set(span):
        if key in ("generator", "metadata"):
            continue
        assert surgery[key] == span[key], key


@pytest.mark.parametrize("family", ["magnitude", "geometry"])
def test_siblings_share_their_anchor_cells_and_differ_elsewhere(family: str) -> None:
    surgery = _columns(_load(PROFILES[("surgery", family)]))
    span = _columns(_load(PROFILES[("span", family)]))
    assert surgery.keys() == span.keys()
    for key in surgery:
        by_effect = {(_mode_effect(cell)): cell for cell in surgery[key]}
        for cell in span[key]:
            other = by_effect[_mode_effect(cell)]
            if cell.metadata.get("zero_effect_anchor"):
                assert cell.cell_id == other.cell_id
                assert parameter_signature(cell) == parameter_signature(other)
                assert _seeds(cell, 100) == _seeds(other, 100)
            else:
                assert cell.cell_id != other.cell_id
                assert cell.generator_params.construction == "stage_span"
                assert replace(cell.generator_params, construction=None, span_layers=None) == other.generator_params
                # Same family, same replicate index: the same baseline program and noise stream.
                assert _seeds(cell, 100) == _seeds(other, 100)


def test_all_four_profiles_anchor_on_the_same_datasets() -> None:
    by_column: dict[Any, set[tuple[int, ...]]] = {}
    for name in PROFILES.values():
        for key, cells in _columns(_load(name)).items():
            anchor = next(cell for cell in cells if cell.metadata.get("zero_effect_anchor"))
            by_column.setdefault(key, set()).add(tuple(_seeds(anchor, 100)))
            assert anchor.generator_params.construction is None
    assert len(by_column) == 3
    assert all(len(seeds) == 1 for seeds in by_column.values())


# ── Reproduction of the design-guidance pilot ───────────────────────────────


@pytest.mark.parametrize("family", ["magnitude", "geometry"])
def test_surgery_sibling_reproduces_the_design_guidance_pilot_cells(family: str) -> None:
    pilot_columns = _columns(_load(DESIGN_GUIDANCE[family]))
    pilot: dict[tuple[Any, ...], SimulationCell] = {}
    for key, cells in pilot_columns.items():
        point = dict(key)
        if point["generator.n_stages"] != 3 or point["evaluation.integration_params.layers"] != _frozen(TWO):
            continue
        for cell in cells:
            pilot[(point[N], *_mode_effect(cell))] = cell

    matched = 0
    for key, cells in _columns(_load(PROFILES[("surgery", family)])).items():
        point = dict(key)
        if point[SIZES] is not None:
            continue
        for cell in cells:
            reference = pilot[(point[N], *_mode_effect(cell))]
            assert cell.generator_params == reference.generator_params
            assert replace(cell.evaluation_params, span_capture=None) == reference.evaluation_params
            assert _seeds(cell, 100) == _seeds(reference, 100)
            assert cell.base_seed == reference.base_seed
            matched += 1
    assert matched == {"magnitude": 8, "geometry": 20}[family]
