#!/usr/bin/env python3
"""Design-guidance figures from the two core design-guidance reports.

Reads ``design_point_operating.csv`` from the reports of the magnitude and the
geometry core profiles (``examples/trajectory_power_study/design_guidance_core_*_pilot.json``
or the paper-grade ``design_guidance_core_*_study.json``) and writes, under
``--out-dir``:

- ``design_guidance_magnitude.png`` (magnitude / ``delta``),
  ``design_guidance_orientation.png`` (orientation / ``angle``) and
  ``design_guidance_shape.png`` (shape / ``shape``): one panel per (stage
  count, block count), sample size on a log axis, one line per nonzero effect
  with ±1 Monte Carlo SE bars;
- ``design_guidance_control.png``: per statistic and panel, the zero-effect
  anchor's and translation's (at its largest effect) rejection rate against
  sample size, with the bound ``alpha + k·sqrt(alpha(1 − alpha)/R)`` drawn for
  the anchor's available replicate count ``R`` at every sample size;
- ``design_guidance_power.csv``: exactly the operating-table rows plotted,
  labelled with the figure, the source profile and the block count.

The engine's ``design_point_power.png`` plots two numeric axes and would mix
the two- and three-block columns on this grid. Here the block coordinate is read
in its CSV spelling — the measured-block list, or empty for all three blocks —
and the two block counts are never merged. The control figure reads the
geometry profile's anchor (the magnitude profile's anchor evaluates the same
datasets and is the report's reproduction check, not a second measurement).

``--grade`` (default ``pilot``) is the grade label written into every figure
title; the panels are derived from the stage counts and block counts present,
so a report with one stage count renders one panel row.

A report lacking any of the three design axes is refused before anything is
written: the script exits 1 and names the missing axis.

Example:

    python scripts/design_guidance_figure.py \\
        --magnitude results/design-guidance-slice-<date>/core_magnitude/report \\
        --geometry results/design-guidance-slice-<date>/core_geometry/report \\
        --out-dir results/design-guidance-slice-<date>/figures --grade paper
"""

from __future__ import annotations

import argparse
import ast
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import ScalarFormatter  # noqa: E402

N_AXIS = "generator.n_samples"
STAGES_AXIS = "generator.n_stages"
LAYERS_AXIS = "evaluation.integration_params.layers"
#: Design axes every core report must carry.
REQUIRED_AXES: tuple[str, ...] = (N_AXIS, STAGES_AXIS, LAYERS_AXIS)
#: Blocks the generator produces; an empty ``layers`` coordinate measures them all.
ALL_BLOCKS = 3
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
#: (figure, trajectory mode, target statistic, source report) of the mode figures.
MODE_FIGURES: tuple[tuple[str, str, str, str], ...] = (
    ("magnitude", "magnitude", "delta", "magnitude"),
    ("orientation", "orientation", "angle", "geometry"),
    ("shape", "shape", "shape", "geometry"),
)
ANCHOR_MODE = "none"
CONTROL_MODE = "translation"
OPERATING_TABLE = "design_point_operating.csv"
CSV_NAME = "design_guidance_power.csv"
#: Grade label the titles carry when none is given (the pilot's).
DEFAULT_GRADE = "pilot"


class DesignGuidanceError(ValueError):
    """Raised when a report cannot feed the design-guidance figure."""


def block_count(value: Any) -> int:
    """Number of measured blocks from the CSV spelling of the ``layers`` coordinate.

    Empty (``NaN``) means the key was absent, so every block was measured; a
    list spelling names the measured blocks.
    """

    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ALL_BLOCKS
    text = str(value).strip()
    if not text:
        return ALL_BLOCKS
    try:
        parsed = ast.literal_eval(text)
    except (ValueError, SyntaxError) as exc:
        raise DesignGuidanceError(f"cannot read the block coordinate {value!r} as a list of blocks.") from exc
    if isinstance(parsed, str):
        return 1
    return len(parsed)


def load_operating_table(report_dir: Path, profile: str) -> pd.DataFrame:
    """Read one report's operating table and label its block count; refuse a missing axis."""

    path = Path(report_dir) / OPERATING_TABLE
    if not path.exists():
        raise DesignGuidanceError(f"{profile} report: {path} does not exist.")
    frame = pd.read_csv(path)
    missing = [axis for axis in REQUIRED_AXES if axis not in frame.columns]
    if missing:
        raise DesignGuidanceError(
            f"{profile} report: {OPERATING_TABLE} lacks the design axis {missing}; "
            f"the design-guidance figure needs every one of {list(REQUIRED_AXES)}."
        )
    frame = frame.copy()
    frame[N_AXIS] = pd.to_numeric(frame[N_AXIS])
    frame[STAGES_AXIS] = pd.to_numeric(frame[STAGES_AXIS])
    frame["effect_size"] = pd.to_numeric(frame["effect_size"])
    frame["n_blocks"] = [block_count(value) for value in frame[LAYERS_AXIS]]
    frame["profile"] = profile
    return frame


def mode_rows(table: pd.DataFrame, mode: str, statistic: str) -> pd.DataFrame:
    """The rows a mode figure plots: the mode's nonzero effects at its target statistic."""

    sub = table[(table["trajectory_mode"] == mode) & (table["statistic"] == statistic) & (table["effect_size"] > 0.0)]
    return sub.sort_values([STAGES_AXIS, "n_blocks", "effect_size", N_AXIS]).reset_index(drop=True)


def control_rows(table: pd.DataFrame) -> pd.DataFrame:
    """The rows the control figure plots: the anchor and translation at its largest effect, every statistic."""

    anchor = table[(table["trajectory_mode"] == ANCHOR_MODE) & (table["effect_size"] == 0.0)]
    translation = table[table["trajectory_mode"] == CONTROL_MODE]
    if not translation.empty:
        translation = translation[translation["effect_size"] == translation["effect_size"].max()]
    rows = pd.concat([anchor, translation])
    rows = rows[rows["statistic"].isin(STATISTICS)]
    return rows.sort_values(["statistic", STAGES_AXIS, "n_blocks", "trajectory_mode", N_AXIS]).reset_index(drop=True)


def build_frames(magnitude: pd.DataFrame, geometry: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Every figure's rows, keyed by figure name, each labelled with its ``figure``."""

    tables = {"magnitude": magnitude, "geometry": geometry}
    frames = {
        figure: mode_rows(tables[profile], mode, statistic).assign(figure=figure)
        for figure, mode, statistic, profile in MODE_FIGURES
    }
    frames["control"] = control_rows(geometry).assign(figure="control")
    return frames


def _panels(rows: pd.DataFrame) -> tuple[list[int], list[int]]:
    return sorted(int(s) for s in rows[STAGES_AXIS].unique()), sorted(int(b) for b in rows["n_blocks"].unique())


def _log_n_axis(ax: Any, n_values: list[float]) -> None:
    ax.set_xscale("log")
    ax.set_xlim(min(n_values) / 1.4, max(n_values) * 1.4)
    ax.set_xticks(n_values)
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.minorticks_off()


def _empty_figure(out_path: Path, title: str) -> Path:
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.axis("off")
    ax.text(0.5, 0.5, title, ha="center", va="center", fontsize=9)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def render_mode_figure(
    rows: pd.DataFrame, out_path: Path, *, mode: str, statistic: str, grade: str = DEFAULT_GRADE
) -> Path:
    """Panels for stages × blocks; rate vs n on a log axis, one line per effect, ±1 SE."""

    out_path = Path(out_path)
    if rows.empty:
        return _empty_figure(out_path, f"{mode}/{statistic}: no rows")
    fig = build_mode_figure(rows, mode=mode, statistic=statistic, grade=grade)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def build_mode_figure(rows: pd.DataFrame, *, mode: str, statistic: str, grade: str = DEFAULT_GRADE) -> Any:
    """The mode figure (one panel per stage count × block count present), titled with the grade label."""

    stages, blocks = _panels(rows)
    n_values = sorted(float(n) for n in rows[N_AXIS].unique())
    fig, axes = plt.subplots(
        len(stages),
        len(blocks),
        figsize=(3.6 * len(blocks), 3.0 * len(stages)),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    for i, n_stages in enumerate(stages):
        for j, n_blocks in enumerate(blocks):
            ax = axes[i][j]
            panel = rows[(rows[STAGES_AXIS] == n_stages) & (rows["n_blocks"] == n_blocks)]
            for effect in sorted(panel["effect_size"].unique()):
                line = panel[panel["effect_size"] == effect].sort_values(N_AXIS)
                ax.errorbar(
                    line[N_AXIS].to_numpy(dtype=float),
                    line["rejection_rate"].to_numpy(dtype=float),
                    yerr=line["monte_carlo_se"].to_numpy(dtype=float),
                    marker="o",
                    capsize=2,
                    label=f"e = {effect:g}",
                )
            _log_n_axis(ax, n_values)
            ax.set_title(f"{n_stages} stages, {n_blocks} blocks", fontsize=9)
            ax.set_ylim(-0.02, 1.02)
            ax.grid(True, alpha=0.3)
            if i == len(stages) - 1:
                ax.set_xlabel("n_samples")
            if j == 0:
                ax.set_ylabel(f"{mode}/{statistic} rejection rate")
    axes[0][0].legend(fontsize=7, title="effect", title_fontsize=7)
    fig.suptitle(f"{mode} / {statistic}: power vs sample size ({grade}; ±1 MC SE)", fontsize=10)
    fig.tight_layout()
    return fig


def type_i_bound(alpha: float, se_multiplier: float, n_available: float) -> float:
    """``alpha + k·sqrt(alpha(1 − alpha)/R)``, ``NaN`` without replicates."""

    if not n_available or n_available <= 0:
        return float("nan")
    return alpha + se_multiplier * math.sqrt(alpha * (1.0 - alpha) / float(n_available))


def render_control_figure(
    rows: pd.DataFrame, out_path: Path, *, alpha: float, se_multiplier: float, grade: str = DEFAULT_GRADE
) -> Path:
    """Rows = statistics, columns = stages × blocks; anchor and translation vs n with the anchor bound."""

    out_path = Path(out_path)
    if rows.empty:
        return _empty_figure(out_path, "control: no rows")
    fig = build_control_figure(rows, alpha=alpha, se_multiplier=se_multiplier, grade=grade)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def build_control_figure(
    rows: pd.DataFrame, *, alpha: float, se_multiplier: float, grade: str = DEFAULT_GRADE
) -> Any:
    """The control figure (statistics × panels present), titled with the grade label."""

    stages, blocks = _panels(rows)
    panels = [(s, b) for s in stages for b in blocks]
    statistics = [statistic for statistic in STATISTICS if statistic in set(rows["statistic"])]
    n_values = sorted(float(n) for n in rows[N_AXIS].unique())
    fig, axes = plt.subplots(
        len(statistics), len(panels), figsize=(3.2 * len(panels), 2.6 * len(statistics)), squeeze=False, sharex=True
    )
    top = rows.loc[rows["trajectory_mode"] == CONTROL_MODE, "effect_size"]
    control_label = f"translation (e = {float(top.max()):g})" if not top.empty else "translation"
    for i, statistic in enumerate(statistics):
        for j, (n_stages, n_blocks) in enumerate(panels):
            ax = axes[i][j]
            panel = rows[
                (rows["statistic"] == statistic) & (rows[STAGES_AXIS] == n_stages) & (rows["n_blocks"] == n_blocks)
            ]
            anchor = panel[panel["trajectory_mode"] == ANCHOR_MODE].sort_values(N_AXIS)
            control = panel[panel["trajectory_mode"] == CONTROL_MODE].sort_values(N_AXIS)
            for block, label, marker in ((anchor, "anchor (e = 0)", "o"), (control, control_label, "s")):
                if block.empty:
                    continue
                ax.errorbar(
                    block[N_AXIS].to_numpy(dtype=float),
                    block["rejection_rate"].to_numpy(dtype=float),
                    yerr=block["monte_carlo_se"].to_numpy(dtype=float),
                    marker=marker,
                    capsize=2,
                    label=label,
                )
            if not anchor.empty:
                bound = [type_i_bound(alpha, se_multiplier, r) for r in anchor["n_available"].to_numpy(dtype=float)]
                ax.plot(
                    anchor[N_AXIS].to_numpy(dtype=float),
                    bound,
                    linestyle="--",
                    color="gray",
                    label=f"α + {se_multiplier:g}·SE bound",
                )
            ax.axhline(alpha, color="black", linewidth=0.6, alpha=0.5)
            _log_n_axis(ax, n_values)
            ax.set_ylim(-0.01, 0.35)
            ax.grid(True, alpha=0.3)
            if i == 0:
                ax.set_title(f"{n_stages} stages, {n_blocks} blocks", fontsize=9)
            if i == len(statistics) - 1:
                ax.set_xlabel("n_samples")
            if j == 0:
                ax.set_ylabel(f"{statistic} rejection rate")
    axes[0][0].legend(fontsize=6)
    fig.suptitle(f"Size controls vs sample size ({grade}; α = {alpha:g}, ±1 MC SE)", fontsize=10)
    fig.tight_layout()
    return fig


def plotted_table(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """The CSV: every plotted row, mode figures first, then the control rows."""

    order = [figure for figure, _, _, _ in MODE_FIGURES] + ["control"]
    table = pd.concat([frames[figure] for figure in order], ignore_index=True)
    leading = ["figure", "profile", "n_blocks"]
    return table[leading + [column for column in table.columns if column not in leading]]


def run(
    magnitude_dir: Path,
    geometry_dir: Path,
    out_dir: Path,
    *,
    alpha: float = 0.05,
    se_multiplier: float = 2.0,
    grade: str = DEFAULT_GRADE,
) -> dict[str, Path]:
    """Validate both reports, then write the four figures and the CSV. Nothing is written on a refusal."""

    magnitude = load_operating_table(Path(magnitude_dir), "magnitude")
    geometry = load_operating_table(Path(geometry_dir), "geometry")
    frames = build_frames(magnitude, geometry)
    table = plotted_table(frames)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for figure, mode, statistic, _ in MODE_FIGURES:
        paths[figure] = render_mode_figure(
            frames[figure], out_dir / f"design_guidance_{figure}.png", mode=mode, statistic=statistic, grade=grade
        )
    paths["control"] = render_control_figure(
        frames["control"],
        out_dir / "design_guidance_control.png",
        alpha=alpha,
        se_multiplier=se_multiplier,
        grade=grade,
    )
    paths["csv"] = out_dir / CSV_NAME
    table.to_csv(paths["csv"], index=False)
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="design_guidance_figure", description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--magnitude", type=Path, required=True, help="Report directory of the core magnitude profile.")
    parser.add_argument("--geometry", type=Path, required=True, help="Report directory of the core geometry profile.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Where to write the figures and the CSV.")
    parser.add_argument("--alpha", type=float, default=0.05, help="Nominal size of the anchor bound (default 0.05).")
    parser.add_argument(
        "--se-multiplier", type=float, default=2.0, help="k in alpha + k*sqrt(alpha(1-alpha)/R) (default 2)."
    )
    parser.add_argument(
        "--grade", default=DEFAULT_GRADE, help="Grade label in every figure title (default 'pilot'; e.g. 'paper')."
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        paths = run(
            args.magnitude,
            args.geometry,
            args.out_dir,
            alpha=args.alpha,
            se_multiplier=args.se_multiplier,
            grade=args.grade,
        )
    except DesignGuidanceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for key, path in paths.items():
        print(f"  {key}: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
