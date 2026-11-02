"""``scripts/presentation_figures.py``: the registry, the adapters, the figures, the writer and the manifest."""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from matplotlib.colors import to_hex  # noqa: E402
from matplotlib.text import Text  # noqa: E402


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


pf = _script("presentation_figures")
ROOT = pf.REPO_ROOT
SMALL_N = 90
#: Every file stem the full set writes (spec: nine figures; mode_simulated and power_vs_n write several).
EXPECTED_STEMS = (
    "mode_geometry",
    "mode_simulated",
    "mode_simulated_anchor",
    "mode_simulated_magnitude",
    "mode_simulated_orientation",
    "mode_simulated_shape",
    "mode_simulated_translation",
    "size_atlas",
    "specificity_n1200",
    "specificity_n80",
    "power_vs_n",
    "power_vs_n_delta",
    "power_vs_n_angle",
    "power_vs_n_shape",
    "power_vs_stage_signal",
    "estimator_comparison",
    "covariates_null",
)


def _csv(path: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / path)


@pytest.fixture(scope="module")
def built() -> Any:
    """Every figure built once from the committed sources (simulated figure at a reduced n)."""

    ctx = pf.BuildContext(root=ROOT, simulated_n=SMALL_N)
    results = {}
    with plt.rc_context(pf.style_rc()):
        for figure_id, spec in pf.FIGURES.items():
            results[figure_id] = spec.build(ctx)
    yield results
    plt.close("all")


@pytest.fixture(scope="module")
def rendered(tmp_path_factory: pytest.TempPathFactory) -> Any:
    out = tmp_path_factory.mktemp("figures")
    report = pf.render(out, dpi=20, simulated_n=SMALL_N)
    return out, report


# --------------------------------------------------------------------------- style and CLI


def test_palette_constants_are_the_validated_set() -> None:
    assert pf.SURFACE == "#fcfcfb"
    assert (pf.INK, pf.INK_SECONDARY, pf.INK_MUTED) == ("#1d1d1b", "#57564f", "#8a8980")
    assert pf.STAT_COLORS == {"delta": "#5a4fcf", "angle": "#e0712c", "shape": "#1a9b7f"}
    assert pf.PALETTE_VALIDATION["colors"] == pf.STAT_COLORS
    assert pf.PALETTE_VALIDATION["surface"] == pf.SURFACE
    assert pf.PALETTE_VALIDATION["result"] == "pass"
    assert pf.CONTROL_COLOR == pf.INK_MUTED
    assert [pf.mode_color(m) for m in ("magnitude", "orientation", "shape", "translation", "none")] == [
        "#5a4fcf", "#e0712c", "#1a9b7f", pf.INK_MUTED, pf.INK_MUTED,
    ]


def test_unknown_figure_exits_nonzero_before_writing(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "out"
    assert pf.main(["--out-dir", str(out), "--figure", "size_atlas", "--figure", "no_such_figure"]) == 2
    assert not out.exists()
    error = capsys.readouterr().err
    assert "no_such_figure" in error
    assert all(name in error for name in pf.FIGURE_IDS)


def test_one_figure_selection_writes_only_that_figure(tmp_path: Path) -> None:
    report = pf.render(tmp_path, ["covariates_null"], dpi=20)
    assert report.failures == {}
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "captions.md", "covariates_null.pdf", "covariates_null.png", "manifest.json",
    ]
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert [entry["id"] for entry in manifest["figures"]] == ["covariates_null"]


def test_every_figure_renders_png_and_pdf(rendered: Any) -> None:
    out, report = rendered
    assert report.failures == {}
    for stem in EXPECTED_STEMS:
        assert (out / f"{stem}.png").is_file(), stem
        assert (out / f"{stem}.pdf").is_file(), stem
    assert [entry["id"] for entry in report.manifest["figures"]] == list(pf.FIGURE_IDS)
    assert not list(out.glob(".*.tmp"))


def test_png_is_written_at_the_requested_dpi(tmp_path: Path) -> None:
    from PIL import Image

    pf.render(tmp_path, ["specificity_n80"], dpi=40)
    with Image.open(tmp_path / "specificity_n80.png") as image:
        assert image.size == (round(pf.CANVAS[0] * 40), round(pf.CANVAS[1] * 40))
        assert round(image.info["dpi"][0]) == 40


def test_failing_builder_leaves_no_file_and_spares_the_others(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def broken(ctx: Any) -> Any:
        fig = pf.new_figure()
        fig.subplots(1, 1)
        raise pf.PresentationFigureError("builder failed on purpose")

    spec = pf.FIGURES["size_atlas"]
    monkeypatch.setitem(pf.FIGURES, "size_atlas", pf.FigureSpec(spec.id, spec.title, spec.message, (), broken))
    report = pf.render(tmp_path, ["size_atlas", "covariates_null"], dpi=20)
    assert "builder failed on purpose" in report.failures["size_atlas"]
    assert not list(tmp_path.glob("size_atlas*"))
    assert (tmp_path / "covariates_null.png").is_file()
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert [f["id"] for f in manifest["failed"]] == ["size_atlas"]


def test_failing_save_removes_every_partial_file(tmp_path: Path) -> None:
    good, bad = pf.new_figure(), pf.new_figure()

    def explode(*args: Any, **kwargs: Any) -> None:
        raise OSError("disk full")

    bad.savefig = explode  # type: ignore[method-assign]
    with pytest.raises(OSError):
        pf._save_atomic([("first", good), ("second", bad)], tmp_path, 20)
    assert list(tmp_path.iterdir()) == []
    plt.close(good)
    plt.close(bad)


def test_manifest_records_hashes_and_the_dirty_flag(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(pf, "git_state", lambda root=ROOT: {"git_revision": "abc123", "worktree_clean": False})
    pf.render(tmp_path, ["specificity_n80"], dpi=20)
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["set"]["worktree_clean"] is False
    assert manifest["set"]["git_revision"] == "abc123"
    for key in ("motco_version", "matplotlib_version", "rendered_at", "seed", "palette"):
        assert key in manifest["set"]
    assert manifest["set"]["palette"]["validation"] == pf.PALETTE_VALIDATION
    (entry,) = manifest["figures"]
    assert entry["message"] and entry["caveats"]
    assert {source["path"] for source in entry["sources"]} == {
        pf.SOURCES["core_magnitude_operating"].path, pf.SOURCES["core_geometry_operating"].path,
    }
    for source in entry["sources"]:
        assert source["sha256"] == hashlib.sha256((ROOT / source["path"]).read_bytes()).hexdigest()
        assert "net displacement" in source["measurement"] and source["run"]
    captions = (tmp_path / "captions.md").read_text(encoding="utf-8")
    assert "DIRTY worktree" in captions and entry["message"] in captions


def test_footer_carries_the_measurement_and_can_be_removed(tmp_path: Path) -> None:
    with plt.rc_context(pf.style_rc()):
        with_footer = pf.FIGURES["specificity_n80"].build(pf.BuildContext(root=ROOT))
        without = pf.FIGURES["specificity_n80"].build(pf.BuildContext(root=ROOT, footer=False))
    (_, fig), (_, bare) = with_footer.outputs[0], without.outputs[0]
    footers = [t for t in fig.findobj(Text) if t.get_gid() == "footer"]
    assert len(footers) == 1
    assert "PLS, CV rank ≥ 2 · net displacement" in footers[0].get_text()
    assert "design-guidance-slice-2026-10-29" in footers[0].get_text()
    assert not [t for t in bare.findobj(Text) if t.get_gid() == "footer"]
    plt.close(fig)
    plt.close(bare)
    pf.render(tmp_path, ["specificity_n80"], footer=False, dpi=20)
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["set"]["footer"] is False
    assert all("CV rank ≥ 2" in source["measurement"] for source in manifest["figures"][0]["sources"])


def test_text_is_at_least_12_pt_and_axes_use_plain_language(built: Any) -> None:
    internal = ("n_samples", "trajectory_mode", "generator.", "evaluation.", "effect_size")
    for figure_id, result in built.items():
        for stem, fig in result.outputs:
            for text in fig.findobj(Text):
                if text.get_visible() and text.get_text().strip():
                    assert text.get_fontsize() >= 12, (stem, text.get_text())
            for ax in fig.axes:
                labels = [ax.get_xlabel(), ax.get_ylabel(), ax.get_title()]
                labels += [t.get_text() for t in ax.get_xticklabels() + ax.get_yticklabels()]
                assert not any(token in label for label in labels for token in internal), (stem, labels)
            assert tuple(round(v, 3) for v in fig.get_size_inches()) == tuple(round(v, 3) for v in pf.CANVAS)
    power = built["power_vs_n"].outputs[0][1]
    assert {ax.get_xlabel() for ax in power.axes} == {"samples (n)"}


# --------------------------------------------------------------------------- adapters


#: Each table source's loader, the rows a figure needs from it and the column its plotted value lives in.
LOADERS: dict[str, Callable[[Path], Any]] = {
    "core_magnitude_operating": lambda root: pf.specificity_n80(root),
    "core_geometry_operating": lambda root: pf.specificity_n80(root),
    "core_magnitude_stage_signal": lambda root: pf.stage_signal_rows(root),
    "core_geometry_stage_signal": lambda root: pf.stage_signal_rows(root),
    "stage_effect_magnitude_operating": lambda root: pf.stage_signal_rows(root),
    "stage_effect_geometry_operating": lambda root: pf.stage_signal_rows(root),
    "stage_effect_magnitude_stage_signal": lambda root: pf.stage_signal_rows(root),
    "stage_effect_geometry_stage_signal": lambda root: pf.stage_signal_rows(root),
    "core_geometry_anchor": lambda root: pf.atlas_readings(root),
    "stage_effect_geometry_anchor": lambda root: pf.atlas_readings(root),
    "nd_type_i_phase5": lambda root: pf.atlas_readings(root),
    "rank_floor_delta0.1": lambda root: pf.atlas_readings(root),
    "rank_floor_delta0.05": lambda root: pf.atlas_readings(root),
    "corrected_permutation": lambda root: pf.atlas_readings(root),
    "adni_design": lambda root: pf.atlas_readings(root),
    "design_guidance_power": lambda root: pf.power_rows(root),
    "nd_power_phase5": lambda root: pf.estimator_rows(root),
    "nd_power_phase6": lambda root: pf.estimator_rows(root),
    "diagnosis_delta0.1": lambda root: pf.diagnosis_readings(root),
    "diagnosis_delta0.05": lambda root: pf.diagnosis_readings(root),
}
#: (column, value) of rows every loader of the source needs, and the source's plotted-value column.
NEEDED: dict[str, tuple[str, Any, str]] = {
    "operating": ("trajectory_mode", "none", "rejection_rate"),
    "stage_signal": ("trajectory_mode", "none", "median_cv_auroc"),
    "anchor": ("role", "anchor", "rate"),
    "rank_floor": ("mode", "anchor", "rate_floor_2"),
    "corrected_permutation": ("mode", "anchor", "rate"),
    "adni_design": ("mode", "anchor", "rate"),
    "design_guidance_power": ("figure", "magnitude", "rejection_rate"),
    "nd_power": ("trajectory_mode", "none", "rejection_rate"),
    "diagnosis": ("mode", "anchor", "rate"),
}


def _needed(key: str) -> tuple[str, Any, str]:
    if key.startswith(("core_magnitude", "stage_effect_magnitude")):
        rate = "median_cv_auroc" if key.endswith("stage_signal") else "rejection_rate"
        return ("trajectory_mode", "magnitude", rate)
    for suffix in ("operating", "stage_signal", "anchor"):
        if key.endswith(suffix):
            return NEEDED[suffix]
    for prefix in ("rank_floor", "nd_power", "diagnosis"):
        if key.startswith(prefix):
            return NEEDED[prefix]
    if key == "nd_type_i_phase5":
        return NEEDED["anchor"]
    return NEEDED[key]


def _copy_sources(root: Path) -> None:
    for source in pf.SOURCES.values():
        if source.columns:
            target = root / source.path
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT / source.path, target)


def test_every_table_source_has_a_loader() -> None:
    assert set(LOADERS) == {key for key, source in pf.SOURCES.items() if source.columns}


def test_every_adapter_loads_its_committed_file() -> None:
    atlas = pf.atlas_readings(ROOT)
    assert len(atlas) == 18 * 2 * 3
    n1200, uses = pf.specificity_n1200(ROOT)
    n80, _ = pf.specificity_n80(ROOT)
    tables = {
        "atlas": atlas,
        "n1200": n1200,
        "n80": n80,
        "power": pf.power_rows(ROOT),
        "stage": pf.stage_signal_rows(ROOT),
        "estimator": pf.estimator_rows(ROOT),
        "diagnosis": pf.diagnosis_readings(ROOT),
    }
    sizes = {"n1200": 15, "n80": 15, "power": 54, "stage": 16, "estimator": 40, "diagnosis": 24}
    for name, size in sizes.items():
        assert len(tables[name]) == size, name
    for name, table in tables.items():
        assert np.isfinite(table["rate"].to_numpy(dtype=float)).all(), name
    assert {"rate", "se", "n_replicates", "source", "measurement"} <= set(atlas.columns)
    assert {"bound", "verdict"} <= set(atlas.columns)


@pytest.mark.parametrize("key", sorted(LOADERS))
@pytest.mark.parametrize("defect", ["missing_file", "missing_column", "missing_reading", "non_finite"])
def test_adapter_refuses_a_defective_source(tmp_path: Path, key: str, defect: str) -> None:
    _copy_sources(tmp_path)
    source = pf.SOURCES[key]
    path = tmp_path / source.path
    column, value, rate_column = _needed(key)
    frame = pd.read_csv(path)
    if defect == "missing_file":
        path.unlink()
        match = "source file is missing"
    elif defect == "missing_column":
        frame.drop(columns=[source.columns[-1]]).to_csv(path, index=False)
        match = f"missing column\\(s\\) \\['{source.columns[-1]}'\\]"
    elif defect == "missing_reading":
        frame[frame[column] != value].to_csv(path, index=False)
        match = "missing expected reading"
    else:
        frame.loc[frame[column] == value, rate_column] = np.nan
        frame.to_csv(path, index=False)
        match = "non-finite"
    with pytest.raises(pf.PresentationFigureError, match=match) as error:
        LOADERS[key](tmp_path)
    assert source.path in str(error.value)


def test_a_defective_source_fails_its_figure_without_output(tmp_path: Path) -> None:
    root = tmp_path / "root"
    _copy_sources(root)
    path = root / pf.SOURCES["adni_design"].path
    pd.read_csv(path).drop(columns=["bound"]).to_csv(path, index=False)
    out = tmp_path / "out"
    report = pf.render(out, ["size_atlas"], root=root, dpi=20)
    assert pf.SOURCES["adni_design"].path in report.failures["size_atlas"]
    assert "bound" in report.failures["size_atlas"]
    assert not list(out.glob("size_atlas*"))


# --------------------------------------------------------------------------- measurement tags (D5)


def _summary(root: Path, *, count: int, mismatches: int) -> None:
    path = root / pf.FLOOR_CHECK_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"count": count, "mismatches": mismatches, "upper_bound_95": 0.02}), encoding="utf-8")


def test_phase5_tag_follows_the_floor_check(tmp_path: Path) -> None:
    phase5, type_i, phase6 = (pf.SOURCES[k] for k in ("nd_power_phase5", "nd_type_i_phase5", "nd_power_phase6"))
    assert pf.resolve_measurement(phase5, tmp_path).floor == "pre-floor"
    assert "CV rank ≥ 1 (pre-floor)" in pf.resolve_measurement(phase5, tmp_path).tag()
    _summary(tmp_path, count=152, mismatches=0)
    for source in (phase5, type_i):
        measurement = pf.resolve_measurement(source, tmp_path)
        assert measurement.floor == "checked"
        assert "CV rank ≥ 2 by floor check (0/152 ranks change)" in measurement.tag()
    assert pf.resolve_measurement(phase6, tmp_path).floor == "pre-floor"
    _summary(tmp_path, count=152, mismatches=1)
    assert pf.resolve_measurement(phase5, tmp_path).floor == "pre-floor"
    assert pf.resolve_measurement(type_i, tmp_path).floor == "pre-floor"


def test_committed_tags_of_the_phase5_and_phase6_sources() -> None:
    check = pf.load_floor_check(ROOT)
    expected = "checked" if check is not None and check.agrees else "pre-floor"
    assert pf.resolve_measurement(pf.SOURCES["nd_power_phase5"], ROOT).floor == expected
    assert pf.resolve_measurement(pf.SOURCES["nd_type_i_phase5"], ROOT).floor == expected
    assert pf.resolve_measurement(pf.SOURCES["nd_power_phase6"], ROOT).floor == "pre-floor"
    for key, source in pf.SOURCES.items():
        if source.measurement is not None and not key.startswith(("nd_power", "nd_type_i")):
            assert pf.resolve_measurement(source, ROOT).floor == "current", key


# --------------------------------------------------------------------------- figures


def test_mode_geometry_colours_b_by_its_target_statistic(built: Any) -> None:
    (_, fig), = built["mode_geometry"].outputs
    lines = {line.get_gid(): line for ax in fig.axes for line in ax.get_lines()}
    for mode, color in (("magnitude", "#5a4fcf"), ("orientation", "#e0712c"), ("shape", "#1a9b7f"),
                        ("translation", pf.INK_MUTED)):
        assert to_hex(lines[f"path-B-{mode}"].get_color()) == color
        assert to_hex(lines[f"path-A-{mode}"].get_color()) == pf.INK_MUTED
    chips = {t.get_gid(): t.get_text() for t in fig.findobj(Text) if (t.get_gid() or "").startswith("chip-")}
    assert chips == {"chip-magnitude": "reads: size (δ)", "chip-orientation": "reads: orientation (θ)",
                     "chip-shape": "reads: shape", "chip-translation": "reads: none"}


def test_mode_simulated_is_deterministic_and_writes_six_figures(built: Any) -> None:
    first = pf.simulated_mode_data(seed=4, n_samples=SMALL_N)
    second = pf.simulated_mode_data(seed=4, n_samples=SMALL_N)
    assert list(first) == list(pf.SIMULATED_MODES)
    for mode in pf.SIMULATED_MODES:
        pd.testing.assert_frame_equal(first[mode]["means"], second[mode]["means"])
        pd.testing.assert_frame_equal(first[mode]["samples"], second[mode]["samples"])
        assert len(first[mode]["means"]) == 6
    stems = [stem for stem, _ in built["mode_simulated"].outputs]
    assert stems == ["mode_simulated", *(pf.SIMULATED_STEMS[m] for m in pf.SIMULATED_MODES)]
    note = pf.DISPLAY_NOTE
    for _, fig in built["mode_simulated"].outputs:
        assert any(note in t.get_text().replace("\n", " ") or "display projection" in t.get_text().lower()
                   for t in fig.findobj(Text))


def test_size_atlas_plots_d6_selection(built: Any) -> None:
    readings = built["size_atlas"].data
    sections = readings.drop_duplicates("row_index").groupby("section", sort=False).size()
    assert sections.tolist() == [6, 3, 1, 2, 5, 1]
    assert set(readings["control"]) == {"anchor", "translation"}
    assert (readings.groupby("row_index").size() == 6).all()
    adni = _csv(pf.SOURCES["adni_design"].path)
    rows = readings[readings["source"] == "adni_design"].sort_values(["control", "statistic"])
    expected = adni.assign(control=adni["mode"]).sort_values(["control", "statistic"])
    np.testing.assert_allclose(rows["rate"], expected["rate"])
    corrected = _csv(pf.SOURCES["corrected_permutation"].path)
    df = corrected[(corrected["tier"] == "reading") & (corrected["arm"] == "df")]
    plotted = readings[readings["source"] == "corrected_permutation"]
    assert sorted(plotted["rate"]) == sorted(df["rate"])
    assert set(plotted["bound"].round(6)) == {0.07924}
    anchors = readings[(readings["source"] == "core_geometry_anchor") & (readings["control"] == "anchor")]
    assert anchors["verdict_recorded"].all() and set(anchors["bound"].round(6)) == {0.069494}
    assert (readings["verdict"] == "holds").all()
    (_, fig), = built["size_atlas"].outputs
    for statistic in pf.STATISTICS:
        for control in ("anchor", "translation"):
            xs = [x for ax in fig.axes for c in ax.collections
                  if c.get_gid() == f"atlas-{statistic}-{control}-holds" for x in c.get_offsets()[:, 0]]
            want = readings[(readings["statistic"] == statistic) & (readings["control"] == control)]["rate"]
            assert sorted(xs) == pytest.approx(sorted(want))


def test_size_atlas_marks_an_exceedance(built: Any) -> None:
    readings = built["size_atlas"].data.copy()
    index = readings[(readings["statistic"] == "angle") & (readings["control"] == "anchor")].index[0]
    readings.loc[index, ["rate", "verdict"]] = [0.09, "exceeds"]
    with plt.rc_context(pf.style_rc()):
        fig = pf.draw_size_atlas(readings, pf.BuildContext(root=ROOT), built["size_atlas"].uses)
    hits = [c for ax in fig.axes for c in ax.collections if c.get_gid() == "atlas-angle-anchor-exceeds"]
    assert len(hits) == 1 and hits[0].get_offsets()[:, 0].tolist() == [0.09]
    labels = [t.get_text() for t in fig.findobj(Text) if t.get_gid() == "atlas-exceed-label-angle"]
    assert labels == ["exceeds: 0.090"]
    plt.close(fig)


def _specificity_image(result: Any) -> Any:
    (_, fig), = result.outputs
    (image,) = [im for ax in fig.axes for im in ax.get_images() if im.get_gid() == "specificity"]
    return image


def test_specificity_cells_equal_their_source_and_share_the_scale(built: Any) -> None:
    phase5 = _csv(pf.SOURCES["nd_power_phase5"].path)
    phase5["evaluation.orientation_estimator"] = phase5["evaluation.orientation_estimator"].fillna("pc1")
    estimator = phase5["evaluation.orientation_estimator"]
    nd = phase5[(estimator == "net_displacement") & (phase5["phase"] == "power_design")]
    image = _specificity_image(built["specificity_n1200"])
    for i, mode in enumerate(pf.SPECIFICITY_ROWS):
        effect = 0.0 if mode == "none" else 1.0
        for j, statistic in enumerate(pf.STATISTICS):
            cell = (nd["trajectory_mode"] == mode) & (nd["effect_size"] == effect)
            want = nd[cell & (nd["statistic"] == statistic)]
            assert image.get_array()[i, j] == pytest.approx(float(want["rejection_rate"].iloc[0]))
    small = _specificity_image(built["specificity_n80"])
    tables = {name: _csv(f"{pf.SLICE}/core_{name}/report/design_point_operating.csv")
              for name in ("magnitude", "geometry")}
    for i, mode in enumerate(pf.SPECIFICITY_ROWS):
        table = tables["magnitude" if mode == "magnitude" else "geometry"]
        effect = 0.0 if mode == "none" else pf.TOP_EFFECT[mode]
        for j, statistic in enumerate(pf.STATISTICS):
            want = table[(table["generator.n_samples"] == 80) & table["evaluation.integration_params.layers"].notna()
                         & (table["trajectory_mode"] == mode) & np.isclose(table["effect_size"], effect)
                         & (table["statistic"] == statistic)]
            assert small.get_array()[i, j] == pytest.approx(float(want["rejection_rate"].iloc[0]))
    assert (image.norm.vmin, image.norm.vmax) == (small.norm.vmin, small.norm.vmax) == (0.0, 1.0)
    assert image.get_cmap().name == small.get_cmap().name
    diagonals = {p.get_gid() for ax in image.figure.axes for p in ax.patches}
    assert {"diagonal-magnitude", "diagonal-orientation", "diagonal-shape"} <= diagonals


def test_power_vs_n_overview_equals_the_top_effect_diagonals(built: Any) -> None:
    result = built["power_vs_n"]
    assert [stem for stem, _ in result.outputs] == ["power_vs_n", "power_vs_n_delta", "power_vs_n_angle",
                                                    "power_vs_n_shape"]
    table = _csv(pf.SOURCES["design_guidance_power"].path)
    overview = result.outputs[0][1]
    lines = {line.get_gid(): line for ax in overview.axes for line in ax.get_lines() if line.get_gid()}
    for figure, mode, statistic, effects in pf.POWER_DIAGONALS:
        for blocks in (2, 3):
            want = table[(table["figure"] == figure) & (table["n_blocks"] == blocks) & (table["statistic"] == statistic)
                         & np.isclose(table["effect_size"], max(effects))].sort_values("generator.n_samples")
            line = lines[f"power-{statistic}-{blocks}"]
            assert list(line.get_xdata()) == list(want["generator.n_samples"])
            assert list(line.get_ydata()) == pytest.approx(list(want["rejection_rate"]))


def test_power_vs_stage_signal_has_three_own_axes_and_the_core_delta(built: Any) -> None:
    (_, fig), = built["power_vs_stage_signal"].outputs
    assert len(fig.axes) == 3
    boxes = [tuple(np.round(ax.get_position().bounds, 4)) for ax in fig.axes]
    assert len(set(boxes)) == 3
    rows = built["power_vs_stage_signal"].data
    core = {name: _csv(f"{pf.SLICE}/core_{name}/report/design_point_operating.csv") for name in ("magnitude",
                                                                                                    "geometry")}
    for mode in ("magnitude", "orientation", "shape"):
        table = core["magnitude" if mode == "magnitude" else "geometry"]
        want = table[(table["generator.n_samples"] == 200) & table["evaluation.integration_params.layers"].notna()
                     & (table["trajectory_mode"] == mode) & np.isclose(table["effect_size"], pf.TOP_EFFECT[mode])
                     & (table["statistic"] == pf.MODE_TARGET[mode])]
        got = rows[(rows["mode"] == mode) & (rows["delta"] == 2.0)]
        assert float(got["rate"].iloc[0]) == pytest.approx(float(want["rejection_rate"].iloc[0]))
        assert float(got["median_rank"].iloc[0]) == float(want["median_selected_lv"].iloc[0])
    assert sorted(rows["delta"].unique()) == [0.25, 0.5, 1.0, 2.0]


def test_estimator_comparison_top_effect_values(built: Any) -> None:
    rows = built["estimator_comparison"].data
    top = rows[(rows["mode"] == "orientation") & (rows["effect_size"] == 1.0)]
    values = {(d, e): float(r) for d, e, r in zip(top["design"], top["estimator"], top["rate"])}
    assert values == {("phase6", "pc1"): 0.300, ("phase6", "net_displacement"): 0.648,
                      ("phase5", "pc1"): 0.840, ("phase5", "net_displacement"): 0.924}
    (_, fig), = built["estimator_comparison"].outputs
    styles = {line.get_gid(): line for ax in fig.axes for line in ax.get_lines() if line.get_gid()}
    pc1, nd = styles["estimator-phase6-orientation-pc1"], styles["estimator-phase6-orientation-net_displacement"]
    assert to_hex(pc1.get_color()) == to_hex(nd.get_color()) == pf.STAT_COLORS["angle"]
    assert pc1.get_linestyle() != nd.get_linestyle()
    assert to_hex(pc1.get_markerfacecolor()) == pf.SURFACE and to_hex(nd.get_markerfacecolor()) != pf.SURFACE


def test_covariates_null_pooled_rates(built: Any) -> None:
    pooled = built["covariates_null"].data["pooled"].set_index(["covariate_columns", "correction"])
    assert pooled.loc[(1, "uncorrected"), "rate"] == pytest.approx(0.0545)
    assert pooled.loc[(6, "uncorrected"), "rate"] == pytest.approx(0.066)
    assert pooled.loc[(6, "df"), "rate"] == pytest.approx(0.0545)
    assert set(pooled["n_replicates"]) == {2000}
    assert pooled["bound"].iloc[0] == pytest.approx(0.05 + 2 * math.sqrt(0.05 * 0.95 / 2000))
    assert round(float(pooled["bound"].iloc[0]), 4) == 0.0597
