"""``scripts/design_guidance_figure.py`` on synthetic design-point operating tables.

The tables are built through ``resolve_operating_by_design_point`` from
hand-made records, so the block coordinate reaches the script in the CSV
spelling the engine actually writes (a tuple string, or empty for all three
blocks).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import pandas as pd
import pytest

matplotlib.use("Agg")

from motco.simulations.grid import SimulationReplicateResult
from motco.simulations.study import resolve_operating_by_design_point, write_report_csvs
from motco.simulations.study.enumerate import DESIGN_AXIS_MARKER, DESIGN_PHASE, DESIGN_POINT_KEY
from motco.simulations.study.report import ReportFrames

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import design_guidance_figure as dgf
finally:
    sys.path.pop(0)

N = "generator.n_samples"
STAGES = "generator.n_stages"
LAYERS = "evaluation.integration_params.layers"
TWO = ("methylation", "expression")
BASELINE = {N: 200, STAGES: 3, LAYERS: TWO}
R = 10


def _rejections(mode: str, effect: float, n: int, stages: int, blocks: int) -> int:
    """A deterministic rejection count out of ``R`` that differs between the block counts."""

    if mode == "none" or mode == "translation":
        return 1 if blocks == 2 else 2
    count = float(effect) * float(n) / 20.0 * (1.0 + 0.5 * (blocks - 2)) * (float(stages) / 3.0)
    return min(R, int(count + 1e-9))


def _rate(mode: str, effect: float, n: int, stages: int, blocks: int) -> float:
    return _rejections(mode, effect, n, stages, blocks) / R


def _record(point: dict, mode: str, effect: float, index: int, rejected: bool) -> SimulationReplicateResult:
    phase = "power_primary" if point == BASELINE else DESIGN_PHASE
    p = 0.01 if rejected else 0.5
    return SimulationReplicateResult(
        cell_id=f"{phase}-{mode}-{effect}-{point[N]}-{point[STAGES]}-{len(point[LAYERS] or ())}",
        phase=phase,
        replicate_index=index,
        replicate_seed=index,
        generator_seed=index,
        evaluation_seed=3,
        parameter_signature="sig",
        status="completed",
        p_values={"delta": p, "angle": p, "shape": p},
        cell_metadata={
            "trajectory_mode": mode,
            "effect_size": effect,
            DESIGN_POINT_KEY: dict(point),
            "varied_axis": DESIGN_AXIS_MARKER if phase == DESIGN_PHASE else None,
        },
        integration_metadata={"integration_method": "pls", "component_selection": "cv", "selected_lv": 2},
        null_summary={"angle": {"q95": 50.0}},
        config_spectrum={"version": 1, "pooled": {"relative_eigengap": 0.1}, "groups": {}},
    )


def _records(
    modes: dict[str, tuple[float, ...]], stages_present: tuple[int, ...] = (3, 4)
) -> list[SimulationReplicateResult]:
    records = []
    for n in (80, 200):
        for stages in stages_present:
            for layers in (TWO, None):
                point = {N: n, STAGES: stages, LAYERS: layers}
                blocks = 3 if layers is None else len(layers)
                cells = [("none", 0.0)] + [(mode, e) for mode, effects in modes.items() for e in effects]
                for mode, effect in cells:
                    n_reject = _rejections(mode, effect, n, stages, blocks)
                    records.extend(_record(point, mode, effect, i, i < n_reject) for i in range(R))
    return records


def _write_report(records: list[SimulationReplicateResult], out_dir: Path) -> Path:
    operating = resolve_operating_by_design_point(records)
    frames = ReportFrames(
        specificity_matrix=pd.DataFrame(), power_curves=pd.DataFrame(), type_i_table=pd.DataFrame(),
        design_point_operating=operating,
    )
    write_report_csvs(frames, out_dir)
    return out_dir


@pytest.fixture
def reports(tmp_path: Path) -> tuple[Path, Path]:
    magnitude = _write_report(_records({"magnitude": (0.05, 0.1)}), tmp_path / "magnitude" / "report")
    geometry = _write_report(
        _records({"orientation": (0.5, 1.0), "shape": (0.5, 1.0), "translation": (0.5, 1.0)}),
        tmp_path / "geometry" / "report",
    )
    return magnitude, geometry


def test_block_count_reads_the_csv_spelling() -> None:
    assert dgf.block_count(float("nan")) == 3
    assert dgf.block_count("") == 3
    assert dgf.block_count("('methylation', 'expression')") == 2
    assert dgf.block_count("('methylation',)") == 1
    with pytest.raises(dgf.DesignGuidanceError):
        dgf.block_count("not a list")


def test_figures_and_csv_are_written(reports: tuple[Path, Path], tmp_path: Path) -> None:
    magnitude, geometry = reports
    out = tmp_path / "figures"
    assert dgf.main(["--magnitude", str(magnitude), "--geometry", str(geometry), "--out-dir", str(out)]) == 0
    for name in ("magnitude", "orientation", "shape", "control"):
        path = out / f"design_guidance_{name}.png"
        assert path.exists() and path.stat().st_size > 0
    assert (out / "design_guidance_power.csv").exists()


def test_csv_rows_are_exactly_the_plotted_rows(reports: tuple[Path, Path], tmp_path: Path) -> None:
    magnitude, geometry = reports
    paths = dgf.run(magnitude, geometry, tmp_path / "figures")
    csv = pd.read_csv(paths["csv"])
    frames = dgf.build_frames(
        dgf.load_operating_table(magnitude, "magnitude"), dgf.load_operating_table(geometry, "geometry")
    )
    plotted = pd.concat(frames.values(), ignore_index=True)
    keys = ["figure", "profile", "trajectory_mode", "statistic", "effect_size", N, STAGES, "n_blocks", "rejection_rate"]
    expected = sorted(map(tuple, plotted[keys].to_records(index=False).tolist()))
    assert sorted(map(tuple, csv[keys].to_records(index=False).tolist())) == expected
    assert len(csv) == len(plotted)
    # Mode figures: every nonzero effect of the target statistic only; control: anchor + translation top.
    assert set(csv.loc[csv["figure"] == "magnitude", "statistic"]) == {"delta"}
    assert set(csv.loc[csv["figure"] == "orientation", "statistic"]) == {"angle"}
    assert set(csv.loc[csv["figure"] == "shape", "statistic"]) == {"shape"}
    assert (csv.loc[csv["figure"] != "control", "effect_size"] > 0).all()
    control = csv[csv["figure"] == "control"]
    assert set(control["statistic"]) == {"delta", "angle", "shape"}
    assert set(zip(control["trajectory_mode"], control["effect_size"])) == {("none", 0.0), ("translation", 1.0)}
    assert set(control["profile"]) == {"geometry"}


def test_two_and_three_block_columns_are_never_merged(reports: tuple[Path, Path], tmp_path: Path) -> None:
    magnitude, geometry = reports
    csv = pd.read_csv(dgf.run(magnitude, geometry, tmp_path / "figures")["csv"])
    assert set(csv["n_blocks"]) == {2, 3}
    grouped = csv.groupby(["figure", "statistic", "trajectory_mode", "effect_size", N, STAGES])
    distinct = 0
    for (figure, _, mode, effect, n, stages), group in grouped:
        # One row per block count at every (figure, statistic, mode, effect, n, stages) ...
        assert sorted(group["n_blocks"]) == [2, 3], (figure, mode, effect, n, stages)
        by_blocks = group.set_index("n_blocks")["rejection_rate"]
        # ... each carrying its own column's rate, never an average of the two.
        for blocks in (2, 3):
            assert by_blocks.loc[blocks] == pytest.approx(_rate(mode, effect, n, stages, blocks)), (figure, mode, n)
        distinct += int(by_blocks.loc[2] != by_blocks.loc[3])
    assert distinct > 0  # the construction makes the block columns differ (anchors: 0.1 vs 0.2)


def test_missing_axis_is_refused_without_writing(reports: tuple[Path, Path], tmp_path: Path, capsys) -> None:
    magnitude, geometry = reports
    table = pd.read_csv(geometry / "design_point_operating.csv").drop(columns=[STAGES])
    table.to_csv(geometry / "design_point_operating.csv", index=False)
    out = tmp_path / "figures"
    assert dgf.main(["--magnitude", str(magnitude), "--geometry", str(geometry), "--out-dir", str(out)]) == 1
    assert STAGES in capsys.readouterr().err
    assert not out.exists()

    # A missing report directory is refused the same way.
    assert dgf.main(["--magnitude", str(magnitude), "--geometry", str(tmp_path / "nope"), "--out-dir", str(out)]) == 1
    assert not out.exists()


# --- grade label and the single-stage-count layout ------------------------------

GEOMETRY_MODES = {"orientation": (0.5, 1.0), "shape": (0.5, 1.0), "translation": (0.5, 1.0)}
PRE_CHANGE_TITLES = {
    "magnitude": "magnitude / delta: power vs sample size (pilot; ±1 MC SE)",
    "orientation": "orientation / angle: power vs sample size (pilot; ±1 MC SE)",
    "shape": "shape / shape: power vs sample size (pilot; ±1 MC SE)",
    "control": "Size controls vs sample size (pilot; α = 0.05, ±1 MC SE)",
}


@pytest.fixture
def single_stage_reports(tmp_path: Path) -> tuple[Path, Path]:
    magnitude = _write_report(_records({"magnitude": (0.05, 0.1)}, (3,)), tmp_path / "magnitude" / "report")
    geometry = _write_report(_records(GEOMETRY_MODES, (3,)), tmp_path / "geometry" / "report")
    return magnitude, geometry


def _frames(magnitude: Path, geometry: Path) -> dict[str, pd.DataFrame]:
    return dgf.build_frames(
        dgf.load_operating_table(magnitude, "magnitude"), dgf.load_operating_table(geometry, "geometry")
    )


def test_single_stage_count_renders_one_panel_row(single_stage_reports: tuple[Path, Path], tmp_path: Path) -> None:
    magnitude, geometry = single_stage_reports
    frames = _frames(magnitude, geometry)
    for figure, mode, statistic, _ in dgf.MODE_FIGURES:
        fig = dgf.build_mode_figure(frames[figure], mode=mode, statistic=statistic)
        assert len(fig.axes) == 2  # one stage count x two block counts
        assert [ax.get_title() for ax in fig.axes] == ["3 stages, 2 blocks", "3 stages, 3 blocks"]
        matplotlib.pyplot.close(fig)
    control = dgf.build_control_figure(frames["control"], alpha=0.05, se_multiplier=2.0)
    assert len(control.axes) == 3 * 2  # three statistics x the same two panels
    matplotlib.pyplot.close(control)

    paths = dgf.run(magnitude, geometry, tmp_path / "figures")
    csv = pd.read_csv(paths["csv"])
    plotted = pd.concat(frames.values(), ignore_index=True)
    assert len(csv) == len(plotted)
    assert set(csv[STAGES]) == {3}
    keys = ["figure", "trajectory_mode", "statistic", "effect_size", N, STAGES, "n_blocks", "rejection_rate"]
    assert sorted(map(tuple, csv[keys].to_records(index=False).tolist())) == sorted(
        map(tuple, plotted[keys].to_records(index=False).tolist())
    )


def test_grade_label_is_in_every_title_and_defaults_to_the_pilot(reports: tuple[Path, Path], tmp_path: Path) -> None:
    magnitude, geometry = reports
    frames = _frames(magnitude, geometry)
    for figure, mode, statistic, _ in dgf.MODE_FIGURES:
        default = dgf.build_mode_figure(frames[figure], mode=mode, statistic=statistic)
        assert default.get_suptitle() == PRE_CHANGE_TITLES[figure]
        paper = dgf.build_mode_figure(frames[figure], mode=mode, statistic=statistic, grade="paper")
        assert paper.get_suptitle() == PRE_CHANGE_TITLES[figure].replace("(pilot;", "(paper;")
        matplotlib.pyplot.close(default)
        matplotlib.pyplot.close(paper)
    default = dgf.build_control_figure(frames["control"], alpha=0.05, se_multiplier=2.0)
    assert default.get_suptitle() == PRE_CHANGE_TITLES["control"]
    paper = dgf.build_control_figure(frames["control"], alpha=0.05, se_multiplier=2.0, grade="paper")
    assert paper.get_suptitle() == PRE_CHANGE_TITLES["control"].replace("(pilot;", "(paper;")
    matplotlib.pyplot.close(default)
    matplotlib.pyplot.close(paper)

    # The CLI accepts the label and still writes the four figures and the CSV.
    out = tmp_path / "figures"
    argv = ["--magnitude", str(magnitude), "--geometry", str(geometry), "--out-dir", str(out), "--grade", "paper"]
    assert dgf.main(argv) == 0
    assert {p.name for p in out.iterdir()} == {
        "design_guidance_magnitude.png",
        "design_guidance_orientation.png",
        "design_guidance_shape.png",
        "design_guidance_control.png",
        "design_guidance_power.csv",
    }
    assert dgf.build_parser().parse_args(["--magnitude", "m", "--geometry", "g", "--out-dir", "o"]).grade == "pilot"
