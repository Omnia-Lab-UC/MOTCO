#!/usr/bin/env python3
"""Standalone, provenance-tagged presentation figures of MOTCO's simulation evidence.

The set behind ``openspec/changes/simulation-presentation-figures/``: nine
figures for talks, each carrying one message, built from committed aggregates
under ``results/`` (never the gitignored replicate records; no simulation study,
no cluster, no R). The one exception is ``mode_simulated``, which generates one
seeded dataset per trajectory mode with the numpy generator.

==========================  ====================================================
``mode_geometry``           what magnitude, orientation, shape and translation do
                            to group B's trajectory, and which statistic reads each
``mode_simulated``          the same modes realized by the generator, in a display
                            projection (combined file + one file per mode)
``size_atlas``              every zero-effect and translation control rate,
                            n = 80 to n = 1,501, against its own α + k·SE bound
``specificity_n1200``       the mode × statistic response at n = 1,200
``specificity_n80``         the same matrix at n = 80, on the same colour scale
``power_vs_n``              target power against n, two blocks against three
                            (overview + one detail file per statistic)
``power_vs_stage_signal``   power, CV-selected rank and stage AUROC as the joint
                            stage effect δ weakens
``estimator_comparison``    PC1 against net displacement at both designs
``covariates_null``         covariate-adjusted ``delta`` size against covariate
                            columns, uncorrected against ``df``-corrected
==========================  ====================================================

Every source is a registry entry (path, run, measurement, required columns); its
adapter checks the columns and every expected reading and fails naming the file
and the missing item, so a figure never silently drops a reading. Each figure
carries its measurement tag and source runs in a footer (``--no-footer`` removes
it) and in ``manifest.json`` / ``captions.md`` beside the figures, with every
source's SHA-256, the git revision and whether the worktree was clean. Every
figure is written as PNG (300 dpi) and PDF on a 13.333 × 7.5 in canvas; writes
are atomic, so a failing figure leaves no file and does not touch another's.

Example:

    python scripts/presentation_figures.py --out-dir results/presentation-figures-<date>
    python scripts/presentation_figures.py --out-dir /tmp/figs --figure size_atlas --figure covariates_null
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import logging
import math
import os
import subprocess
import sys
import textwrap
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, MaxNLocator, NullFormatter, NullLocator  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import design_guidance_figure as dgf  # noqa: E402

sys.path.pop(0)

REPO_ROOT = Path(__file__).resolve().parents[1]

# --------------------------------------------------------------------------- style (design D3)

SURFACE = "#fcfcfb"
INK = "#1d1d1b"
INK_SECONDARY = "#57564f"
INK_MUTED = "#8a8980"
#: Hairline grid, one step off the surface.
GRID = "#e4e3de"
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
#: One validated colour per statistic (``validate_palette.js --mode light --pairs all`` against SURFACE).
STAT_COLORS: dict[str, str] = {"delta": "#5a4fcf", "angle": "#e0712c", "shape": "#1a9b7f"}
#: Controls (anchor, translation) and group A.
CONTROL_COLOR = INK_MUTED
PALETTE_VALIDATION: dict[str, Any] = {
    "date": "2026-11-02",
    "tool": "dataviz skill scripts/validate_palette.js",
    "arguments": "--mode light --pairs all",
    "surface": SURFACE,
    "colors": dict(STAT_COLORS),
    "result": "pass",
    "lightness_band": "pass (all inside L 0.43-0.77)",
    "chroma_floor": "pass (all >= 0.1)",
    "worst_cvd_delta_e": {"value": 10.1, "pair": ["shape", "angle"], "deficiency": "protan"},
    "worst_normal_delta_e": {"value": 24.7, "pair": ["shape", "angle"]},
    "contrast": "pass (all >= 3:1 against the surface)",
    "rejected_triplets": "three (lightness band; CVD dE 4.3; warn band 6-8)",
}
#: Neutral slate ramp of the heatmaps (light → dark); every cell also prints its value.
HEATMAP_RAMP: tuple[str, str] = ("#eef0f2", "#2b3440")
STAT_LABELS: dict[str, str] = {"delta": "size (δ)", "angle": "orientation (θ)", "shape": "shape"}
STAT_SYMBOLS: dict[str, str] = {"delta": "δ", "angle": "θ", "shape": "shape"}
MODE_TARGET: dict[str, str] = {"magnitude": "delta", "orientation": "angle", "shape": "shape"}
CANVAS: tuple[float, float] = (13.333, 7.5)
DPI = 300
TYPE_SIZES: dict[str, float] = {"tick": 15, "label": 18, "title": 18, "annotation": 14, "footer": 12}
FONT_FAMILY: tuple[str, ...] = ("Inter", "Source Sans 3", "DejaVu Sans")
ALPHA = 0.05
SE_MULTIPLIER = 2.0


def mode_color(mode: str) -> str:
    """A mode takes its target statistic's colour; the controls take neutral ink."""

    return STAT_COLORS[MODE_TARGET[mode]] if mode in MODE_TARGET else CONTROL_COLOR


def style_rc() -> dict[str, Any]:
    """The rc parameters every figure is built and saved under."""

    return {
        "font.family": "sans-serif",
        "font.sans-serif": list(FONT_FAMILY),
        "font.size": TYPE_SIZES["tick"],
        "axes.labelsize": TYPE_SIZES["label"],
        "axes.titlesize": TYPE_SIZES["title"],
        "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "xtick.labelsize": TYPE_SIZES["tick"],
        "ytick.labelsize": TYPE_SIZES["tick"],
        "legend.fontsize": TYPE_SIZES["annotation"],
        "legend.frameon": False,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": INK_MUTED,
        "axes.labelcolor": INK,
        "axes.titlecolor": INK,
        "text.color": INK,
        "xtick.color": INK_SECONDARY,
        "ytick.color": INK_SECONDARY,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": GRID,
        "grid.linewidth": 1.0,
        "grid.linestyle": "-",
        "lines.linewidth": 2.5,
        "lines.markersize": 9,
        "pdf.fonttype": 42,
        "figure.dpi": 100,
    }


logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# --------------------------------------------------------------------------- measurement and sources (D2, D5)


class PresentationFigureError(ValueError):
    """Raised when a figure cannot be built from its sources."""


@dataclass(frozen=True)
class Measurement:
    """What a source measured: rank rule, estimator, magnitude, permutation correction, R × permutations."""

    estimator: str
    correction: str
    replicates: int
    permutations: int
    #: ``current`` (floor 2), ``pre-floor`` (floor 1) or ``checked`` (floor 1, shown equivalent to floor 2).
    floor: str = "current"
    magnitude: str = "joint"
    integration: str = "PLS"
    covariates: str | None = None
    check: str | None = None

    def rank_rule(self) -> str:
        if self.floor == "current":
            return "CV rank ≥ 2"
        if self.floor == "checked":
            return f"CV rank ≥ 2 by floor check ({self.check})"
        return "CV rank ≥ 1 (pre-floor)"

    def fields(self) -> dict[str, str]:
        """The tag's parts by name, in tag order."""

        parts = {
            "rank": f"{self.integration}, {self.rank_rule()}",
            "estimator": self.estimator,
            "magnitude": f"{self.magnitude} magnitude",
        }
        if self.covariates:
            parts["covariates"] = self.covariates
        parts["correction"] = self.correction
        parts["size"] = f"{self.replicates} × {self.permutations}"
        return parts

    def tag(self) -> str:
        return " · ".join(self.fields().values())

    def short(self) -> str:
        """The in-figure label of a non-current measurement (empty for the current one)."""

        if self.floor == "current":
            return ""
        if self.floor == "checked":
            return f"pre-floor run; floor check: {self.check}"
        return "pre-floor (CV rank ≥ 1)"

    def as_dict(self) -> dict[str, Any]:
        return {
            "integration": self.integration,
            "rank_rule": self.rank_rule(),
            "floor": self.floor,
            "orientation_estimator": self.estimator,
            "magnitude": self.magnitude,
            "covariates": self.covariates,
            "permutation_correction": self.correction,
            "replicates": self.replicates,
            "permutations": self.permutations,
            "tag": self.tag(),
        }


@dataclass(frozen=True)
class Source:
    """One committed source table."""

    key: str
    path: str
    run: str
    measurement: Measurement | None
    columns: tuple[str, ...] = ()
    #: A pre-floor source whose tag follows the committed floor-equivalence check (D5).
    floor_check: bool = False
    #: Resolved against the repository rather than ``--root`` (package data).
    repo_file: bool = False

    @property
    def run_name(self) -> str:
        """The run directory under ``results/`` (the footer's short source name)."""

        parts = self.path.split("/")
        return parts[1] if parts[0] == "results" and len(parts) > 2 else self.run


@dataclass(frozen=True)
class FloorCheck:
    """The committed floor-equivalence check's summary."""

    path: str
    count: int
    mismatches: int
    upper_bound: float

    @property
    def agrees(self) -> bool:
        return self.count > 0 and self.mismatches == 0

    def citation(self) -> str:
        return f"{self.mismatches}/{self.count} ranks change"


N_AXIS = dgf.N_AXIS
LAYERS_AXIS = dgf.LAYERS_AXIS
DELTA_AXES: tuple[str, ...] = ("generator.delta_methyl", "generator.delta_expr", "generator.delta_protein")
ESTIMATOR_AXIS = "evaluation.orientation_estimator"
METHOD_AXIS = "evaluation.integration_method"
FORCED_AXIS = "evaluation.integration_params.forced_components"
N_COMPONENTS_AXIS = "evaluation.integration_params.n_components"
OPERATING_COLUMNS: tuple[str, ...] = (
    "phase",
    "trajectory_mode",
    "effect_size",
    "statistic",
    "n_replicates",
    "n_available",
    "rejection_rate",
    "monte_carlo_se",
    "median_selected_lv",
)
STAGE_SIGNAL_COLUMNS: tuple[str, ...] = ("trajectory_mode", "effect_size", "n_auroc_available", "median_cv_auroc")
ANCHOR_COLUMNS: tuple[str, ...] = (
    "role",
    "phase",
    "trajectory_mode",
    "effect_size",
    "statistic",
    "rate",
    "monte_carlo_se",
    "n_available",
    "bound",
    "classification",
)
RANK_FLOOR_COLUMNS: tuple[str, ...] = (
    "mode",
    "effect",
    "statistic",
    "replicates_floor_2",
    "rate_floor_2",
    "se_floor_2",
    "bound",
    "verdict",
)
SIZE_READING_COLUMNS: tuple[str, ...] = (
    "tier",
    "run",
    "mode",
    "effect",
    "statistic",
    "arm",
    "replicates",
    "rate",
    "bound",
    "verdict",
)
ADNI_COLUMNS: tuple[str, ...] = (
    "mode",
    "effect",
    "statistic",
    "replicates",
    "rate",
    "monte_carlo_se",
    "bound",
    "verdict",
)
DIAGNOSIS_COLUMNS: tuple[str, ...] = (
    "run",
    "mode",
    "effect",
    "arm",
    "statistic",
    "replicates",
    "rejections",
    "rate",
    "bound",
    "verdict",
)
POWER_COLUMNS: tuple[str, ...] = (
    "figure",
    "n_blocks",
    N_AXIS,
    "trajectory_mode",
    "effect_size",
    "statistic",
    "n_available",
    "rejection_rate",
    "monte_carlo_se",
)

CURRENT = Measurement(estimator="net displacement", correction="uncorrected RRPP", replicates=500, permutations=999)
BOTH_ESTIMATORS = "PC1 and net displacement on the same datasets"
PHASE5 = replace(CURRENT, floor="pre-floor", estimator=BOTH_ESTIMATORS)
PHASE6 = replace(CURRENT, floor="pre-floor", estimator=BOTH_ESTIMATORS)
CASE_STUDY = replace(CURRENT, covariates="6 covariate columns", correction="df-corrected RRPP")
ADNI = replace(CURRENT, covariates="6 covariate columns")
DIAGNOSIS = replace(CURRENT, covariates="1 / 3 / 6 covariate columns", correction="uncorrected vs df-corrected RRPP")

SLICE = "results/design-guidance-slice-2026-10-29"
SLICE_RUN = "design-guidance-slice-2026-10-29"
FLOOR_CHECK_PATH = "results/rank-floor-equivalence-2026-11-02/summary.json"
FLOOR_CHECK_RUN = "rank-floor-equivalence-2026-11-02"
REFERENCE_DATA = "src/motco/simulations/data/intersim_reference.npz"


def _slice_sources() -> list[Source]:
    sources = []
    for profile, axes in (
        ("core_magnitude", (N_AXIS, LAYERS_AXIS)),
        ("core_geometry", (N_AXIS, LAYERS_AXIS)),
        ("stage_effect_magnitude", DELTA_AXES),
        ("stage_effect_geometry", DELTA_AXES),
    ):
        run = f"{SLICE_RUN} · {profile}"
        report = f"{SLICE}/{profile}/report"
        sources.append(
            Source(
                f"{profile}_operating", f"{report}/design_point_operating.csv", run, CURRENT, OPERATING_COLUMNS + axes
            )
        )
        sources.append(
            Source(
                f"{profile}_stage_signal", f"{report}/design_point_stage_signal.csv", run, CURRENT,
                STAGE_SIGNAL_COLUMNS + axes,
            )
        )
        if profile.endswith("geometry"):
            sources.append(
                Source(f"{profile}_anchor", f"{report}/anchor_decision.csv", run, CURRENT, ANCHOR_COLUMNS + axes)
            )
    return sources


SOURCES: dict[str, Source] = {
    source.key: source
    for source in (
        *_slice_sources(),
        Source(
            "design_guidance_power", f"{SLICE}/figures/design_guidance_power.csv", f"{SLICE_RUN} · figures", CURRENT,
            POWER_COLUMNS,
        ),
        Source(
            "nd_power_phase5",
            "results/net-displacement-power-2026-10-14/phase5/report/design_point_operating.csv",
            "net-displacement-power-2026-10-14 · phase5",
            PHASE5,
            OPERATING_COLUMNS + (ESTIMATOR_AXIS,),
            floor_check=True,
        ),
        Source(
            "nd_power_phase6",
            "results/net-displacement-power-2026-10-14/phase6/report/design_point_operating.csv",
            "net-displacement-power-2026-10-14 · phase6",
            PHASE6,
            OPERATING_COLUMNS + (ESTIMATOR_AXIS,),
        ),
        Source(
            "nd_type_i_phase5",
            "results/net-displacement-type-i-2026-10-13/phase5/report/anchor_decision.csv",
            "net-displacement-type-i-2026-10-13 · phase5",
            replace(PHASE5, estimator="net displacement"),
            ANCHOR_COLUMNS + (METHOD_AXIS, FORCED_AXIS, N_COMPONENTS_AXIS, ESTIMATOR_AXIS),
            floor_check=True,
        ),
        *(
            Source(
                f"rank_floor_delta{delta}", f"results/rank-floor-2026-10-23/delta{delta}/rank_floor_rates.csv",
                f"rank-floor-2026-10-23 · delta{delta} (floor_2 arm)", CURRENT, RANK_FLOOR_COLUMNS,
            )
            for delta in ("0.1", "0.05")
        ),
        Source(
            "corrected_permutation", "results/corrected-permutation-2026-10-30/size_readings.csv",
            "corrected-permutation-2026-10-30 (df arm)", CASE_STUDY, SIZE_READING_COLUMNS,
        ),
        Source(
            "adni_design", "results/adni-design-2026-11-02/run/size_readings.csv", "adni-design-2026-11-02 · run",
            ADNI, ADNI_COLUMNS,
        ),
        *(
            Source(
                f"diagnosis_delta{delta}", f"results/case-study-diagnosis-2026-10-30/delta{delta}/diagnosis_rates.csv",
                f"case-study-diagnosis-2026-10-30 · delta{delta}", DIAGNOSIS, DIAGNOSIS_COLUMNS,
            )
            for delta in ("0.1", "0.05")
        ),
        Source("floor_check", FLOOR_CHECK_PATH, FLOOR_CHECK_RUN, None),
        Source("intersim_reference", REFERENCE_DATA, "numpy generator reference data", None, repo_file=True),
    )
}


def load_floor_check(root: Path) -> FloorCheck | None:
    """The committed floor-equivalence summary, or ``None`` when it is absent."""

    path = Path(root) / FLOOR_CHECK_PATH
    if not path.is_file():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    try:
        return FloorCheck(FLOOR_CHECK_PATH, int(data["count"]), int(data["mismatches"]), float(data["upper_bound_95"]))
    except KeyError as exc:
        raise PresentationFigureError(f"{FLOOR_CHECK_PATH}: missing key {exc.args[0]!r}.") from exc


def resolve_measurement(source: Source, root: Path) -> Measurement:
    """The source's measurement, with a pre-floor tag lifted only by an agreeing committed floor check."""

    if source.measurement is None:
        raise PresentationFigureError(f"{source.path}: carries no measurement.")
    if not source.floor_check:
        return source.measurement
    check = load_floor_check(root)
    if check is None or not check.agrees:
        return replace(source.measurement, floor="pre-floor", check=None)
    return replace(source.measurement, floor="checked", check=check.citation())


def read_source(source: Source, root: Path) -> pd.DataFrame:
    """Read a source table and check its columns; refuse a missing file or column, naming it."""

    path = Path(root) / source.path
    if not path.is_file():
        raise PresentationFigureError(f"{source.path}: source file is missing.")
    frame = pd.read_csv(path)
    missing = [column for column in source.columns if column not in frame.columns]
    if missing:
        raise PresentationFigureError(f"{source.path}: missing column(s) {missing}.")
    if LAYERS_AXIS in frame.columns:
        frame["n_blocks"] = [dgf.block_count(value) for value in frame[LAYERS_AXIS]]
    if ESTIMATOR_AXIS in frame.columns:
        frame[ESTIMATOR_AXIS] = frame[ESTIMATOR_AXIS].fillna("pc1")
    return frame


def _describe(criteria: dict[str, Any]) -> str:
    return ", ".join(f"{key}={value!r}" for key, value in criteria.items())


def select_one(frame: pd.DataFrame, criteria: dict[str, Any], source: Source) -> pd.Series:
    """The single row matching every criterion (``None`` matches empty); refuse none or several."""

    mask = np.ones(len(frame), dtype=bool)
    for column, value in criteria.items():
        series = frame[column]
        if value is None:
            hit = series.isna().to_numpy()
        elif isinstance(value, float):
            numeric = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
            hit = np.isclose(numeric, value, rtol=0.0, atol=1e-9)
        else:
            hit = (series == value).to_numpy()
        mask &= hit
    hits = frame[mask]
    if hits.empty:
        raise PresentationFigureError(f"{source.path}: missing expected reading ({_describe(criteria)}).")
    if len(hits) > 1:
        raise PresentationFigureError(f"{source.path}: {len(hits)} rows match ({_describe(criteria)}); expected one.")
    return hits.iloc[0]


def finite(value: Any, source: Source, what: str) -> float:
    """``value`` as a finite float; refuse a missing or non-finite plotted value, naming it."""

    number = float(value) if value is not None else float("nan")
    if not math.isfinite(number):
        raise PresentationFigureError(f"{source.path}: non-finite {what}.")
    return number


def binomial_se(rate: float, replicates: float) -> float:
    return math.sqrt(rate * (1.0 - rate) / replicates) if replicates > 0 else float("nan")


@dataclass(frozen=True)
class Use:
    """A source as one figure uses it: an optional element label and the measurement it plots."""

    source: Source
    measurement: Measurement | None
    label: str = ""


def _use(key: str, root: Path, label: str = "", **override: Any) -> Use:
    source = SOURCES[key]
    measurement = resolve_measurement(source, root)
    return Use(source, replace(measurement, **override) if override else measurement, label)


# --------------------------------------------------------------------------- adapters (tidy frames, D2)


def _tidy(
    source: Source, measurement: Measurement, *, rate: Any, se: Any, n: Any, bound: Any = None, verdict: Any = None,
    what: str, **keys: Any,
) -> dict[str, Any]:
    row = {
        **keys,
        "rate": finite(rate, source, f"rate of {what}"),
        "se": finite(se, source, f"SE of {what}"),
        "n_replicates": int(finite(n, source, f"replicates of {what}")),
        "bound": None if bound is None else finite(bound, source, f"bound of {what}"),
        "verdict": verdict,
        "source": source.key,
        "run": source.run,
        "measurement": measurement.tag(),
    }
    return row


def operating_rates(key: str, root: Path, cells: Iterable[dict[str, Any]], **override: Any) -> pd.DataFrame:
    """One tidy row per expected operating-table cell (rate, MC SE, replicates, median CV rank)."""

    source = SOURCES[key]
    frame = read_source(source, root)
    measurement = _use(key, root, **override).measurement
    assert measurement is not None
    rows = []
    for criteria in cells:
        hit = select_one(frame, criteria, source)
        what = _describe(criteria)
        row = _tidy(
            source, measurement, rate=hit["rejection_rate"], se=hit["monte_carlo_se"], n=hit["n_available"],
            what=what, **criteria,
        )
        row["median_rank"] = float(hit["median_selected_lv"])
        rows.append(row)
    return pd.DataFrame(rows)


def stage_signal(key: str, root: Path, cells: Iterable[dict[str, Any]]) -> pd.DataFrame:
    """One row per expected stage-signal cell: the median held-out CV stage AUROC."""

    source = SOURCES[key]
    frame = read_source(source, root)
    rows = []
    for criteria in cells:
        hit = select_one(frame, criteria, source)
        auroc = finite(hit["median_cv_auroc"], source, f"median CV AUROC of {_describe(criteria)}")
        rows.append({**criteria, "median_auroc": auroc, "source": source.key, "run": source.run})
    return pd.DataFrame(rows)


def anchor_table_readings(key: str, root: Path, design: dict[str, Any], **labels: Any) -> list[dict[str, Any]]:
    """Anchor and translation (e = 1.0) readings of an ``anchor_decision.csv`` design column.

    The anchor carries its recorded bound and classification. The table records
    none for translation, so its bound is the same ``alpha + 2·SE_alpha`` at its
    own replicate count and its verdict is derived against it.
    """

    source = SOURCES[key]
    frame = read_source(source, root)
    frame = frame[frame["phase"].fillna("") != "type_i_baseline"]
    measurement = _use(key, root).measurement
    assert measurement is not None
    rows = []
    for control, criteria in (
        ("anchor", {"role": "anchor", "trajectory_mode": "none", "effect_size": 0.0}),
        ("translation", {"trajectory_mode": "translation", "effect_size": 1.0}),
    ):
        for statistic in STATISTICS:
            full = {**design, **criteria, "statistic": statistic}
            hit = select_one(frame, full, source)
            what = _describe(full)
            if control == "anchor":
                bound, verdict, recorded = hit["bound"], str(hit["classification"]), True
            else:
                bound = dgf.type_i_bound(ALPHA, SE_MULTIPLIER, float(hit["n_available"]))
                verdict = "holds" if float(hit["rate"]) <= bound else "exceeds"
                recorded = False
            rows.append(
                {
                    **_tidy(
                        source, measurement, rate=hit["rate"], se=hit["monte_carlo_se"], n=hit["n_available"],
                        bound=bound, verdict=verdict, what=what, control=control, statistic=statistic,
                    ),
                    "verdict_recorded": recorded,
                    **labels,
                }
            )
    return rows


def rank_floor_readings(key: str, root: Path, **labels: Any) -> list[dict[str, Any]]:
    """The floor-2 arm's anchor and translation readings of a ``rank_floor_rates.csv``."""

    source = SOURCES[key]
    frame = read_source(source, root)
    measurement = _use(key, root).measurement
    assert measurement is not None
    rows = []
    for control, effect in (("anchor", 0.0), ("translation", 1.0)):
        for statistic in STATISTICS:
            criteria = {"mode": control, "effect": effect, "statistic": statistic}
            hit = select_one(frame, criteria, source)
            rows.append(
                {
                    **_tidy(
                        source, measurement, rate=hit["rate_floor_2"], se=hit["se_floor_2"],
                        n=hit["replicates_floor_2"], bound=hit["bound"], verdict=str(hit["verdict"]),
                        what=_describe(criteria), control=control, statistic=statistic,
                    ),
                    "verdict_recorded": True,
                    **labels,
                }
            )
    return rows


def size_reading_rows(key: str, root: Path, *, criteria: dict[str, Any], **labels: Any) -> list[dict[str, Any]]:
    """Anchor and translation readings of a ``size_readings.csv`` (corrected permutation or ADNI)."""

    source = SOURCES[key]
    frame = read_source(source, root)
    measurement = _use(key, root).measurement
    assert measurement is not None
    rows = []
    for control, effect in (("anchor", 0.0), ("translation", 1.0)):
        for statistic in STATISTICS:
            full = {**criteria, "mode": control, "effect": effect, "statistic": statistic}
            hit = select_one(frame, full, source)
            rate = finite(hit["rate"], source, f"rate of {_describe(full)}")
            replicates = finite(hit["replicates"], source, f"replicates of {_describe(full)}")
            se = hit["monte_carlo_se"] if "monte_carlo_se" in hit.index else binomial_se(rate, replicates)
            rows.append(
                {
                    **_tidy(
                        source, measurement, rate=rate, se=se, n=replicates, bound=hit["bound"],
                        verdict=str(hit["verdict"]), what=_describe(full), control=control, statistic=statistic,
                    ),
                    "verdict_recorded": True,
                    **labels,
                }
            )
    return rows


def diagnosis_readings(root: Path) -> pd.DataFrame:
    """Covariate-adjusted ``delta`` readings: 1/3/6 columns × uncorrected/df × 2 strengths × 2 controls."""

    rows = []
    for delta in ("0.1", "0.05"):
        key = f"diagnosis_delta{delta}"
        source = SOURCES[key]
        frame = read_source(source, root)
        measurement = _use(key, root).measurement
        assert measurement is not None
        for noise, columns in ((0, 1), (2, 3), (5, 6)):
            for correction, suffix in (("uncorrected", ""), ("df", "_corrected")):
                for control, effect in (("anchor", 0.0), ("translation", 1.0)):
                    criteria = {
                        "run": f"delta{delta}",
                        "mode": control,
                        "effect": effect,
                        "arm": f"z_noise{noise}{suffix}",
                        "statistic": "delta",
                    }
                    hit = select_one(frame, criteria, source)
                    what = _describe(criteria)
                    rate = finite(hit["rate"], source, f"rate of {what}")
                    rows.append(
                        {
                            **_tidy(
                                source, measurement, rate=rate, se=binomial_se(rate, float(hit["replicates"])),
                                n=hit["replicates"], bound=hit["bound"], verdict=str(hit["verdict"]), what=what,
                                control=control, statistic="delta",
                            ),
                            "rejections": int(finite(hit["rejections"], source, f"rejections of {what}")),
                            "stage_effect": float(delta),
                            "covariate_columns": columns,
                            "correction": correction,
                        }
                    )
    return pd.DataFrame(rows)


def pooled_diagnosis(readings: pd.DataFrame) -> pd.DataFrame:
    """Pooled rate and bound per (covariate columns, correction) over the four readings."""

    pooled = readings.groupby(["covariate_columns", "correction"], as_index=False)[["rejections", "n_replicates"]].sum()
    pooled["rate"] = pooled["rejections"] / pooled["n_replicates"]
    pooled["bound"] = [dgf.type_i_bound(ALPHA, SE_MULTIPLIER, r) for r in pooled["n_replicates"]]
    return pooled


def power_rows(root: Path) -> pd.DataFrame:
    """Diagonal power rows of ``design_guidance_power.csv``: every (figure, blocks, effect, n) expected."""

    source = SOURCES["design_guidance_power"]
    frame = read_source(source, root)
    measurement = _use("design_guidance_power", root).measurement
    assert measurement is not None
    rows = []
    for figure, mode, statistic, effects in POWER_DIAGONALS:
        for n_blocks in (2, 3):
            for effect in effects:
                for n in POWER_N:
                    criteria = {
                        "figure": figure,
                        "trajectory_mode": mode,
                        "statistic": statistic,
                        "n_blocks": n_blocks,
                        "effect_size": effect,
                        N_AXIS: n,
                    }
                    hit = select_one(frame, criteria, source)
                    rows.append(
                        _tidy(
                            source, measurement, rate=hit["rejection_rate"], se=hit["monte_carlo_se"],
                            n=hit["n_available"], what=_describe(criteria), mode=mode, statistic=statistic,
                            n_blocks=n_blocks, effect=effect, n_samples=n,
                        )
                    )
    return pd.DataFrame(rows)


#: (figure, mode, target statistic, effects) of the design-guidance diagonals.
POWER_DIAGONALS: tuple[tuple[str, str, str, tuple[float, ...]], ...] = (
    ("magnitude", "magnitude", "delta", (0.02, 0.05, 0.1)),
    ("orientation", "orientation", "angle", (0.25, 0.5, 1.0)),
    ("shape", "shape", "shape", (0.25, 0.5, 1.0)),
)
POWER_N: tuple[int, ...] = (80, 200, 500)
#: Top effect of each mode in the slice (magnitude on its own axis).
TOP_EFFECT: dict[str, float] = {"magnitude": 0.1, "orientation": 1.0, "shape": 1.0, "translation": 1.0}
STAGE_DELTAS: tuple[float, ...] = (0.25, 0.5, 1.0)
CORE_DELTA = 2.0


# --------------------------------------------------------------------------- figure plumbing


@dataclass
class BuildContext:
    root: Path = REPO_ROOT
    seed: int = 0
    footer: bool = True
    titles: bool = False
    #: Samples per simulated dataset (reduced in tests).
    simulated_n: int = 1000


@dataclass
class FigureResult:
    """A built figure: its output files (stem → figure), the sources it used and the plotted data."""

    outputs: list[tuple[str, Figure]]
    uses: list[Use]
    data: pd.DataFrame | dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class FigureSpec:
    id: str
    title: str
    message: str
    caveats: tuple[str, ...]
    build: Callable[[BuildContext], FigureResult]


def new_figure(**gridspec: Any) -> Figure:
    return plt.figure(figsize=CANVAS, layout="constrained", **gridspec)


FOOTER_WIDTH = 150


def footer_lines(uses: Sequence[Use], note: str | None = None) -> list[str]:
    """The measurement tag and source runs: the parts every source shares, then each distinct rest."""

    lines = textwrap.wrap(note, width=FOOTER_WIDTH) if note else []
    measured = [use for use in uses if use.measurement is not None]
    if not measured:
        return lines
    fields = [use.measurement.fields() for use in measured if use.measurement is not None]
    common = [key for key, value in fields[0].items() if all(other.get(key) == value for other in fields)]
    groups: dict[tuple[str, str], list[str]] = {}
    for use, parts in zip(measured, fields):
        rest = " · ".join(value for key, value in parts.items() if key not in common)
        runs = groups.setdefault((use.label, rest), [])
        if use.source.run_name not in runs:
            runs.append(use.source.run_name)
    shared = " · ".join(fields[0][key] for key in common)
    if len({rest for _, rest in groups}) == 1:
        rest = next(iter(groups))[1]
        runs = list(dict.fromkeys(run for group in groups.values() for run in group))
        tag = " · ".join(part for part in (shared, rest) if part)
        text = f"Measurement: {tag}  —  source: {'; '.join(runs)}"
        return lines + textwrap.wrap(text, width=FOOTER_WIDTH, subsequent_indent="    ")
    pieces = [f"{label + ': ' if label else ''}{rest or 'as above'} ({', '.join(runs)})" for (label, rest), runs in
              groups.items()]
    lines += textwrap.wrap(f"Measurement: {shared}" if shared else "Measurement:", width=FOOTER_WIDTH)
    return lines + textwrap.wrap("  |  ".join(pieces), width=FOOTER_WIDTH, initial_indent="    ",
                                 subsequent_indent="    ")


def finish(
    fig: Figure, ctx: BuildContext, spec: FigureSpec, uses: Sequence[Use], *, note: str | None = None
) -> Figure:
    """Add the optional title and the footer, reserving their space in the layout."""

    bottom = 0.0
    if ctx.footer:
        lines = footer_lines(uses, note)
        if lines:
            height_pt = len(lines) * TYPE_SIZES["footer"] * 1.4 + 12
            bottom = height_pt / (CANVAS[1] * 72)
            fig.text(
                0.012, 0.012, "\n".join(lines), ha="left", va="bottom", fontsize=TYPE_SIZES["footer"],
                color=INK_SECONDARY, linespacing=1.3, gid="footer",
            )
    if ctx.titles:
        fig.suptitle(spec.title, x=0.012, ha="left", fontsize=22, fontweight="semibold", color=INK)
    engine = fig.get_layout_engine()
    if engine is not None and bottom:
        engine.set(rect=(0.0, bottom, 1.0, 1.0 - bottom))
    return fig


def value_grid(ax: Any, axis: str = "y") -> None:
    ax.grid(True, axis=axis, color=GRID, linewidth=1.0)
    ax.set_axisbelow(True)


def log_axis(ax: Any, values: Sequence[float], *, labels: Sequence[str] | None = None, pad: float = 1.35) -> None:
    ax.set_xscale("log")
    ax.set_xlim(min(values) / pad, max(values) * pad)
    ax.xaxis.set_major_locator(FixedLocator(list(values)))
    ax.set_xticklabels(list(labels) if labels is not None else [f"{v:g}" for v in values])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())


def spread(values: Sequence[float], gap: float, top: float | None = None) -> list[float]:
    """Label positions at least ``gap`` apart, moved as little as possible and kept at or below ``top``."""

    out = [float(v) for v in values]
    order = sorted(range(len(out)), key=lambda i: out[i])
    last: float | None = None
    for index in order:
        if last is not None and out[index] - last < gap:
            out[index] = last + gap
        last = out[index]
    if top is not None and order and out[order[-1]] > top:
        ceiling = top
        for index in reversed(order):
            out[index] = min(out[index], ceiling)
            ceiling = out[index] - gap
    return out


def end_labels(
    ax: Any, points: Sequence[tuple[float, float, str]], *, gap: float, dx: float, log: bool = False,
    top: float | None = None,
) -> None:
    """Direct labels at line ends, spread vertically with a hairline leader when moved."""

    if not points:
        return
    ys = spread([p[1] for p in points], gap, top)
    for (x, y, text), y_label in zip(points, ys):
        x_label = x * dx if log else x + dx
        leader = {"arrowstyle": "-", "color": INK_MUTED, "lw": 0.8, "shrinkA": 0, "shrinkB": 4}
        ax.annotate(
            text, xy=(x, y), xytext=(x_label, y_label), textcoords="data", ha="left", va="center",
            fontsize=TYPE_SIZES["annotation"], color=INK, arrowprops=leader if abs(y_label - y) > 1e-9 else None,
        )


def reference_line(ax: Any, value: float, text: str, *, axis: str = "y", color: str = INK_MUTED) -> None:
    """A solid hairline reference (α, a bound) with an ink label."""

    if axis == "y":
        ax.axhline(value, color=color, linewidth=1.2, zorder=1)
        ax.annotate(
            text, xy=(1.0, value), xycoords=("axes fraction", "data"), xytext=(-4, 4), textcoords="offset points",
            ha="right", va="bottom", fontsize=TYPE_SIZES["annotation"], color=INK_SECONDARY,
        )
    else:
        ax.axvline(value, color=color, linewidth=1.2, zorder=1)


# --------------------------------------------------------------------------- figure 1: mode geometry

#: Group A's constructed three-stage trajectory.
SCHEMATIC_A = np.array([[0.0, 0.0], [1.0, 0.75], [2.2, 1.0]])


def schematic_paths() -> dict[str, np.ndarray]:
    """Group B's constructed path per mode (design D4)."""

    a = SCHEMATIC_A
    angle = math.radians(35.0)
    rotation = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
    bent = a.copy()
    bent[1] = a[1] + np.array([0.35, -0.95])
    return {
        "magnitude": a[0] + 1.6 * (a - a[0]),
        "orientation": a[0] + (a - a[0]) @ rotation.T,
        "shape": bent,
        "translation": a + np.array([0.45, -0.95]),
    }


def _path(ax: Any, xy: np.ndarray, *, color: str, group: str, gid: str, arrows: bool = True) -> None:
    filled = group == "B"
    ax.plot(
        xy[:, 0], xy[:, 1], color=color, linewidth=3.0 if filled else 2.0, marker="o",
        markersize=10 if filled else 9, markerfacecolor=color if filled else SURFACE, markeredgecolor=color,
        markeredgewidth=2.0, zorder=3 if filled else 2, gid=gid, solid_capstyle="round",
    )
    if arrows:
        for start, end in zip(xy[:-1], xy[1:]):
            middle = start + 0.55 * (end - start)
            ax.annotate(
                "", xy=middle, xytext=start + 0.45 * (end - start),
                arrowprops={"arrowstyle": "-|>,head_length=0.7,head_width=0.4", "color": color, "lw": 0,
                            "shrinkA": 0, "shrinkB": 0},
                zorder=4,
            )


def build_mode_geometry(ctx: BuildContext) -> FigureResult:
    spec = FIGURES["mode_geometry"]
    fig = new_figure()
    axes = fig.subplots(1, 4)
    paths = schematic_paths()
    reads = {"magnitude": "size (δ)", "orientation": "orientation (θ)", "shape": "shape", "translation": "none"}
    titles = {"magnitude": "Magnitude", "orientation": "Orientation", "shape": "Shape", "translation": "Translation"}
    every = np.vstack([SCHEMATIC_A, *paths.values()])
    lo, hi = every.min(axis=0) - 0.45, every.max(axis=0) + 0.45
    for ax, (mode, b) in zip(axes, paths.items()):
        color = mode_color(mode)
        _path(ax, SCHEMATIC_A, color=CONTROL_COLOR, group="A", gid=f"path-A-{mode}")
        _path(ax, b, color=color, group="B", gid=f"path-B-{mode}")
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_title(titles[mode], fontsize=TYPE_SIZES["title"], fontweight="semibold")
        ax.text(
            0.5, -0.04, f"reads: {reads[mode]}", transform=ax.transAxes, ha="center", va="top",
            fontsize=TYPE_SIZES["label"], color=INK, gid=f"chip-{mode}",
            bbox={"boxstyle": "round,pad=0.45,rounding_size=0.8", "facecolor": SURFACE, "edgecolor": color,
                  "linewidth": 2.5},
        )
        for index, (x, y) in enumerate(SCHEMATIC_A, start=1):
            ax.annotate(
                f"{index}", xy=(x, y), xytext=(4, -10), textcoords="offset points", ha="left", va="top",
                fontsize=TYPE_SIZES["annotation"], color=INK_SECONDARY,
            )
    handles = [
        Line2D([], [], color=CONTROL_COLOR, linewidth=2.0, marker="o", markerfacecolor=SURFACE, markersize=9,
               markeredgewidth=2.0, label="group A (reference trajectory)"),
        Line2D([], [], color=INK_SECONDARY, linewidth=3.0, marker="o", markersize=10,
               label="group B (coloured by the statistic that reads the mode)"),
    ]
    fig.legend(handles=handles, loc="outside upper center", ncol=2, fontsize=TYPE_SIZES["annotation"])
    finish(fig, ctx, spec, [], note="Schematic: stage means 1 → 2 → 3 of two groups on constructed coordinates "
                                    "(no data, no measurement).")
    return FigureResult([("mode_geometry", fig)], [], {"paths": paths})


# --------------------------------------------------------------------------- figure 2: simulated modes

SIMULATED_MODES: tuple[str, ...] = ("none", "magnitude", "orientation", "shape", "translation")
SIMULATED_STEMS: dict[str, str] = {
    "none": "mode_simulated_anchor",
    "magnitude": "mode_simulated_magnitude",
    "orientation": "mode_simulated_orientation",
    "shape": "mode_simulated_shape",
    "translation": "mode_simulated_translation",
}
SIMULATED_TITLES: dict[str, str] = {
    "none": "Zero effect (anchor)",
    "magnitude": "Magnitude",
    "orientation": "Orientation",
    "shape": "Shape",
    "translation": "Translation",
}
DISPLAY_NOTE = "Display projection (PCA of the six group × stage means per panel), not the measurement space."


def display_axes(means: np.ndarray, groups: Sequence[str], mode: str) -> np.ndarray:
    """Two orthonormal display axes: the first two principal axes of the cell means.

    Translation's offset is orthogonal to the stage geometry, so its panel keeps
    the first axis and replaces the second with the group offset (mean B − mean A
    over stages) orthogonalized against it.
    """

    centred = means - means.mean(axis=0)
    _, _, components = np.linalg.svd(centred, full_matrices=False)
    first, second = components[0], components[1]
    if mode == "translation":
        labels = np.asarray(groups)
        offset = means[labels == "B"].mean(axis=0) - means[labels == "A"].mean(axis=0)
        offset = offset - (offset @ first) * first
        second = offset / np.linalg.norm(offset)
    return np.vstack([first, second])


def simulated_mode_data(seed: int, n_samples: int) -> dict[str, dict[str, pd.DataFrame]]:
    """Per mode: sample scores and group × stage means in a 2-D display projection (deterministic)."""

    from motco.simulations.evaluation import SimulationEvaluationParams, integrate_semisynthetic_dataset
    from motco.simulations.showcase import generate_showcase_datasets

    datasets = generate_showcase_datasets(
        modes=SIMULATED_MODES, seed=seed, n_samples=n_samples, n_stages=3, effect_size=1.0
    )
    params = SimulationEvaluationParams(integration_method="concat", group_col="group", stage_col="stage")
    out: dict[str, dict[str, pd.DataFrame]] = {}
    for mode, dataset in datasets.items():
        matrix = integrate_semisynthetic_dataset(dataset, params).matrix.to_numpy(dtype=float)
        metadata = dataset.metadata.reset_index(drop=True)
        cells = metadata.groupby(["group", "stage"]).indices
        keys = sorted(cells)
        means = np.vstack([matrix[cells[key]].mean(axis=0) for key in keys])
        axes = display_axes(means, [group for group, _ in keys], mode)
        centre = means.mean(axis=0)
        scores = (matrix - centre) @ axes.T
        samples = pd.DataFrame({"x": scores[:, 0], "y": scores[:, 1], "group": metadata["group"],
                                "stage": metadata["stage"]})
        reference = samples[samples["group"] == "A"]
        stages = sorted(metadata["stage"].unique())
        span = reference[reference["stage"] == stages[-1]][["x", "y"]].mean() - reference[
            reference["stage"] == stages[0]][["x", "y"]].mean()
        if span["x"] < 0:
            samples["x"] = -samples["x"]
        flip_y = (samples[samples["group"] == "B"]["y"].mean() < reference["y"].mean()) if mode == "translation" \
            else span["y"] < 0
        if flip_y:
            samples["y"] = -samples["y"]
        cell_means = samples.groupby(["group", "stage"], as_index=False)[["x", "y"]].mean()
        out[mode] = {"samples": samples, "means": cell_means}
    return out


def _simulated_panel(ax: Any, mode: str, data: dict[str, pd.DataFrame], *, note: bool) -> None:
    color = mode_color(mode)
    samples, means = data["samples"], data["means"]
    for group in ("A", "B"):
        cloud = samples[samples["group"] == group]
        tint = CONTROL_COLOR if group == "A" else color
        ax.scatter(
            cloud["x"], cloud["y"], s=14, facecolors="none" if group == "A" else tint, edgecolors=tint,
            linewidths=0.6, alpha=0.16, zorder=1, rasterized=True,
        )
    for group in ("A", "B"):
        path = means[means["group"] == group].sort_values("stage")[["x", "y"]].to_numpy()
        _path(ax, path, color=CONTROL_COLOR if group == "A" else color, group=group, gid=f"means-{group}-{mode}")
    ax.set_title(SIMULATED_TITLES[mode], fontsize=TYPE_SIZES["title"], fontweight="semibold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("display axis 1", fontsize=TYPE_SIZES["annotation"], color=INK_SECONDARY)
    second = "group offset axis" if mode == "translation" else "display axis 2"
    ax.set_ylabel(second, fontsize=TYPE_SIZES["annotation"], color=INK_SECONDARY)
    if note:
        ax.text(
            0.99, 0.02, DISPLAY_NOTE, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=TYPE_SIZES["annotation"], color=INK_SECONDARY,
        )


def _simulated_handles(mode: str | None = None) -> list[Line2D]:
    b_color = INK_SECONDARY if mode is None else mode_color(mode)
    return [
        Line2D([], [], color=CONTROL_COLOR, linewidth=2.0, marker="o", markerfacecolor=SURFACE, markersize=9,
               markeredgewidth=2.0, label="group A stage means (open)"),
        Line2D([], [], color=b_color, linewidth=3.0, marker="o", markersize=10, label="group B stage means (filled)"),
    ]


def build_mode_simulated(ctx: BuildContext) -> FigureResult:
    spec = FIGURES["mode_simulated"]
    data = simulated_mode_data(ctx.seed, ctx.simulated_n)
    use = Use(SOURCES["intersim_reference"], None)
    outputs: list[tuple[str, Figure]] = []
    combined = new_figure()
    axes = combined.subplots(2, 3).ravel()
    for ax, mode in zip(axes, SIMULATED_MODES):
        _simulated_panel(ax, mode, data[mode], note=False)
    key = axes[-1]
    key.axis("off")
    key.legend(handles=_simulated_handles(), loc="upper left", fontsize=TYPE_SIZES["annotation"])
    key.text(
        0.0, 0.42,
        textwrap.fill(DISPLAY_NOTE, 34) + f"\n\nOne dataset per mode: seed {ctx.seed}, n = {ctx.simulated_n:,},"
        "\n3 stages, effect e = 1.0; faint dots are samples.",
        transform=key.transAxes, ha="left", va="top", fontsize=TYPE_SIZES["annotation"], color=INK_SECONDARY,
    )
    note = (
        f"Generated: generate_showcase_datasets(seed={ctx.seed}, n_samples={ctx.simulated_n}, n_stages=3, "
        "effect_size=1.0), numpy generator; display projection only, no measurement."
    )
    outputs.append(("mode_simulated", finish(combined, ctx, spec, [], note=note)))
    for mode in SIMULATED_MODES:
        fig = new_figure()
        ax = fig.subplots(1, 1)
        _simulated_panel(ax, mode, data[mode], note=True)
        ax.legend(handles=_simulated_handles(mode), loc="upper left", fontsize=TYPE_SIZES["annotation"])
        outputs.append((SIMULATED_STEMS[mode], finish(fig, ctx, spec, [], note=note)))
    return FigureResult(outputs, [use], data)


# --------------------------------------------------------------------------- figure 3: size atlas (D6)


def atlas_readings(root: Path) -> pd.DataFrame:
    """Every size-atlas reading in D6's order: section, row, control, statistic, rate, bound, verdict."""

    rows: list[dict[str, Any]] = []

    def add(section: str, label: str, readings: list[dict[str, Any]]) -> None:
        for reading in readings:
            rows.append({"section": section, "row": label, **reading})

    section = "Balanced · 3 stages · no covariates"
    for n in (80, 200, 500):
        for blocks in (2, 3):
            add(section, f"n = {n} · {blocks} blocks",
                anchor_table_readings("core_geometry_anchor", root, {N_AXIS: n, "n_blocks": blocks}))
    section = "Weaker stage signal · n = 200"
    for delta in STAGE_DELTAS:
        add(section, f"δ = {delta:g} · 2 blocks",
            anchor_table_readings("stage_effect_geometry_anchor", root, {axis: delta for axis in DELTA_AXES}))
    section = "n = 1,200 · 4 stages · 3 blocks"
    phase5 = resolve_measurement(SOURCES["nd_type_i_phase5"], root).short()
    add(section, phase5 or "n = 1,200", anchor_table_readings(
        "nd_type_i_phase5", root,
        {METHOD_AXIS: "pls", FORCED_AXIS: None, N_COMPONENTS_AXIS: None, ESTIMATOR_AXIS: "net_displacement"},
    ))
    section = "SEA-AD table (n = 80) · weak stage signal"
    for delta in ("0.1", "0.05"):
        add(section, f"δ = {delta} · 2 blocks", rank_floor_readings(f"rank_floor_delta{delta}", root))
    section = "SEA-AD configuration · 6 covariates · df"
    for label, run in (
        ("astrocyte (n = 80) · δ = 0.1", "type-i-astrocyte-delta0.1"),
        ("astrocyte (n = 80) · δ = 0.05", "type-i-astrocyte-delta0.05"),
        ("microglia (n = 76) · δ = 0.1", "type-i-microglia-delta0.1"),
        ("microglia (n = 76) · δ = 0.05", "type-i-microglia-delta0.05"),
        ("astrocyte · δ = 2 (power run)", "power"),
    ):
        add(section, label, size_reading_rows(
            "corrected_permutation", root, criteria={"tier": "reading", "arm": "df", "run": run}))
    section = "ADNI design · 6 covariates"
    add(section, "n = 1,501 · uncorrected", size_reading_rows("adni_design", root, criteria={}))
    frame = pd.DataFrame(rows)
    frame["section_index"] = pd.factorize(frame["section"])[0]
    frame["row_index"] = pd.factorize(frame["section"] + "|" + frame["row"])[0]
    return frame


ATLAS_KEYS: tuple[str, ...] = (
    "core_geometry_anchor",
    "stage_effect_geometry_anchor",
    "nd_type_i_phase5",
    "rank_floor_delta0.1",
    "rank_floor_delta0.05",
    "corrected_permutation",
    "adni_design",
)


def _atlas_layout(readings: pd.DataFrame) -> tuple[dict[str, float], dict[int, float], float]:
    """y of each section header and row (one slot each, top to bottom)."""

    headers: dict[str, float] = {}
    rows: dict[int, float] = {}
    y = 0.0
    for section in pd.unique(readings["section"]):
        headers[section] = y
        y += 1.0
        for row_index in pd.unique(readings.loc[readings["section"] == section, "row_index"]):
            rows[int(row_index)] = y
            y += 1.0
        y += 0.25
    return headers, rows, y


def draw_size_atlas(readings: pd.DataFrame, ctx: BuildContext, uses: Sequence[Use]) -> Figure:
    """The sectioned dot plot: a label column and one panel per statistic."""

    spec = FIGURES["size_atlas"]
    fig = new_figure()
    grid = fig.add_gridspec(1, 4, width_ratios=[1.7, 1, 1, 1], wspace=0.09)
    label_ax = fig.add_subplot(grid[0, 0])
    panels = [fig.add_subplot(grid[0, i + 1], sharey=label_ax) for i in range(3)]
    headers, row_y, height = _atlas_layout(readings)
    for ax in (label_ax, *panels):
        ax.set_ylim(height - 0.25, -0.7)
    label_ax.axis("off")
    for section, y in headers.items():
        label_ax.text(0.0, y, section, ha="left", va="center", fontsize=TYPE_SIZES["annotation"],
                      fontweight="semibold", color=INK)
    for row_index, y in row_y.items():
        label = readings.loc[readings["row_index"] == row_index, "row"].iloc[0]
        label_ax.text(0.04, y, label, ha="left", va="center", fontsize=TYPE_SIZES["annotation"], color=INK_SECONDARY)
    rate_max = max(0.1, float(readings["rate"].max()) + 0.01, float(readings["bound"].max()) + 0.01)
    offsets = {"anchor": -0.17, "translation": 0.17}
    for ax, statistic in zip(panels, STATISTICS):
        color = STAT_COLORS[statistic]
        sub = readings[readings["statistic"] == statistic]
        ax.axvline(ALPHA, color=INK_MUTED, linewidth=1.2, zorder=1)
        bounds = sub.groupby("row_index")["bound"].max()
        ax.scatter(
            bounds.to_numpy(dtype=float), [row_y[int(i)] for i in bounds.index], marker="|", s=260,
            color=INK_SECONDARY, linewidths=2.0, zorder=2, gid=f"atlas-bound-{statistic}",
        )
        for control in ("anchor", "translation"):
            for exceeds in (False, True):
                part = sub[(sub["control"] == control) & ((sub["verdict"] == "exceeds") == exceeds)]
                if part.empty:
                    continue
                ys = [row_y[int(i)] + offsets[control] for i in part["row_index"]]
                ax.scatter(
                    part["rate"].to_numpy(dtype=float), ys, marker="D" if exceeds else "o",
                    s=150 if exceeds else 80, facecolors=color if control == "anchor" else SURFACE,
                    edgecolors=INK if exceeds else color, linewidths=2.0, zorder=3,
                    gid=f"atlas-{statistic}-{control}-{'exceeds' if exceeds else 'holds'}",
                )
                if exceeds:
                    for rate, y in zip(part["rate"], ys):
                        ax.annotate(
                            f"exceeds: {rate:.3f}", xy=(rate, y), xytext=(8, 0), textcoords="offset points",
                            ha="left", va="center", fontsize=TYPE_SIZES["annotation"], color=INK,
                            gid=f"atlas-exceed-label-{statistic}",
                        )
        ax.set_xlim(0.0, rate_max)
        ticks = [t for t in (0.0, 0.05, 0.10, 0.15) if t <= rate_max + 1e-9]
        ax.set_xticks(ticks, ["0" if t == 0 else f"{t:.2f}" for t in ticks])
        ax.set_title(STAT_LABELS[statistic], fontsize=TYPE_SIZES["title"], fontweight="semibold")
        ax.set_xlabel("rejection rate", fontsize=TYPE_SIZES["annotation"])
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.spines["left"].set_visible(False)
        value_grid(ax, "x")
        for y in headers.values():
            ax.axhline(y - 0.45, color=GRID, linewidth=1.0, zorder=0)
    handles = [
        Line2D([], [], linestyle="none", marker="o", markersize=9, markerfacecolor=INK_SECONDARY,
               markeredgecolor=INK_SECONDARY, label="zero-effect anchor"),
        Line2D([], [], linestyle="none", marker="o", markersize=9, markerfacecolor=SURFACE,
               markeredgecolor=INK_SECONDARY, markeredgewidth=2.0, label="translation (e = 1)"),
        Line2D([], [], linestyle="none", marker="|", markersize=16, markeredgewidth=2.0, color=INK_SECONDARY,
               label="bound α + k·SE"),
        Line2D([], [], color=INK_MUTED, linewidth=1.2, label="α = 0.05"),
        Line2D([], [], linestyle="none", marker="D", markersize=10, markerfacecolor=SURFACE, markeredgecolor=INK,
               markeredgewidth=2.0, label="exceeds its bound"),
    ]
    fig.legend(handles=handles, loc="outside upper right", ncol=5, fontsize=TYPE_SIZES["annotation"])
    return finish(fig, ctx, spec, uses)


def build_size_atlas(ctx: BuildContext) -> FigureResult:
    readings = atlas_readings(ctx.root)
    uses = [_use(key, ctx.root) for key in ATLAS_KEYS]
    uses += _floor_check_uses(ctx.root, uses)
    return FigureResult([("size_atlas", draw_size_atlas(readings, ctx, uses))], uses, readings)


def _floor_check_uses(root: Path, uses: Sequence[Use]) -> list[Use]:
    """The floor-equivalence summary as a source whenever a figure plots a floor-check source."""

    if any(use.source.floor_check for use in uses) and load_floor_check(root) is not None:
        return [Use(SOURCES["floor_check"], None)]
    return []


# --------------------------------------------------------------------------- figures 4–5: specificity

SPECIFICITY_ROWS: tuple[str, ...] = ("none", "magnitude", "orientation", "shape", "translation")
HEATMAP_NORM: tuple[float, float] = (0.0, 1.0)


def heatmap_cmap() -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("motco_slate", list(HEATMAP_RAMP))


def specificity_n1200(root: Path) -> tuple[pd.DataFrame, list[Use]]:
    cells = []
    for mode in SPECIFICITY_ROWS:
        effect = 0.0 if mode == "none" else 1.0
        for statistic in STATISTICS:
            cells.append({ESTIMATOR_AXIS: "net_displacement", "phase": "power_design", "trajectory_mode": mode,
                          "effect_size": effect, "statistic": statistic})
    use = _use("nd_power_phase5", root, estimator="net displacement")
    frame = operating_rates("nd_power_phase5", root, cells, estimator="net displacement")
    return frame, [use, *_floor_check_uses(root, [use])]


def specificity_n80(root: Path) -> tuple[pd.DataFrame, list[Use]]:
    frames = []
    for mode in SPECIFICITY_ROWS:
        key = "core_magnitude_operating" if mode == "magnitude" else "core_geometry_operating"
        effect = 0.0 if mode == "none" else TOP_EFFECT[mode]
        cells = [{N_AXIS: 80, "n_blocks": 2, "trajectory_mode": mode, "effect_size": effect, "statistic": statistic}
                 for statistic in STATISTICS]
        frames.append(operating_rates(key, root, cells))
    uses = [_use("core_magnitude_operating", root), _use("core_geometry_operating", root)]
    return pd.concat(frames, ignore_index=True), uses


def specificity_matrix(frame: pd.DataFrame) -> np.ndarray:
    return np.array(
        [[float(frame[(frame["trajectory_mode"] == mode) & (frame["statistic"] == s)]["rate"].iloc[0])
          for s in STATISTICS] for mode in SPECIFICITY_ROWS]
    )


def draw_specificity(frame: pd.DataFrame, ctx: BuildContext, spec: FigureSpec, uses: Sequence[Use], *,
                     design: str) -> Figure:
    fig = new_figure()
    grid = fig.add_gridspec(1, 4, width_ratios=[0.55, 2.2, 0.06, 0.55])
    for i in (0, 3):
        fig.add_subplot(grid[0, i]).axis("off")
    ax = fig.add_subplot(grid[0, 1])
    cax = fig.add_subplot(grid[0, 2])
    matrix = specificity_matrix(frame)
    image = ax.imshow(matrix, cmap=heatmap_cmap(), norm=Normalize(*HEATMAP_NORM), aspect="auto", gid="specificity")
    for i, mode in enumerate(SPECIFICITY_ROWS):
        for j, statistic in enumerate(STATISTICS):
            value = matrix[i, j]
            red, green, blue = to_rgb(heatmap_cmap()(value))
            light = 0.2126 * red + 0.7152 * green + 0.0722 * blue > 0.5
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=TYPE_SIZES["label"],
                    color=INK if light else SURFACE)
            if MODE_TARGET.get(mode) == statistic:
                ax.add_patch(Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, fill=False, edgecolor=STAT_COLORS[statistic],
                                       linewidth=4.0, gid=f"diagonal-{mode}"))
    effects = frame.drop_duplicates("trajectory_mode").set_index("trajectory_mode")["effect_size"]
    labels = ["zero effect (anchor)" if m == "none" else f"{m} · e = {float(effects[m]):g}" for m in SPECIFICITY_ROWS]
    ax.set_yticks(range(len(SPECIFICITY_ROWS)), labels)
    ax.set_xticks(range(len(STATISTICS)), [STAT_LABELS[s] for s in STATISTICS])
    ax.xaxis.tick_top()
    ax.tick_params(length=0, labelsize=TYPE_SIZES["tick"] + 1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel(design, fontsize=TYPE_SIZES["annotation"], color=INK_SECONDARY)
    bar = fig.colorbar(image, cax=cax)
    bar.set_label("rejection rate", fontsize=TYPE_SIZES["annotation"])
    bar.outline.set_visible(False)
    return finish(fig, ctx, spec, uses)


def build_specificity_n1200(ctx: BuildContext) -> FigureResult:
    frame, uses = specificity_n1200(ctx.root)
    spec = FIGURES["specificity_n1200"]
    note = uses[0].measurement.short() if uses[0].measurement is not None else ""
    design = "n = 1,200 · 4 stages · 3 blocks · net displacement" + (f" · {note}" if note else "")
    return FigureResult([("specificity_n1200", draw_specificity(frame, ctx, spec, uses, design=design))], uses, frame)


def build_specificity_n80(ctx: BuildContext) -> FigureResult:
    frame, uses = specificity_n80(ctx.root)
    spec = FIGURES["specificity_n80"]
    design = "n = 80 · 3 stages · 2 blocks (methylation + expression) · net displacement"
    return FigureResult([("specificity_n80", draw_specificity(frame, ctx, spec, uses, design=design))], uses, frame)


# --------------------------------------------------------------------------- figure 6: power vs n

BLOCK_STYLES: dict[int, dict[str, Any]] = {
    2: {"linestyle": "-", "label": "2 blocks (methylation + expression)"},
    3: {"linestyle": (0, (5, 3)), "label": "3 blocks (+ protein)"},
}
EFFECT_MARKERS: tuple[str, ...] = ("o", "s", "^")
MODE_NAMES: dict[str, str] = {"magnitude": "Magnitude", "orientation": "Orientation", "shape": "Shape"}


def build_power_vs_n(ctx: BuildContext) -> FigureResult:
    spec = FIGURES["power_vs_n"]
    rows = power_rows(ctx.root)
    uses = [_use("design_guidance_power", ctx.root)]
    outputs = []
    fig = new_figure()
    axes = fig.subplots(1, 3, sharey=True)
    for ax, (_, mode, statistic, effects) in zip(axes, POWER_DIAGONALS):
        color = STAT_COLORS[statistic]
        top = max(effects)
        labels = []
        for n_blocks, style in BLOCK_STYLES.items():
            line = rows[(rows["mode"] == mode) & (rows["n_blocks"] == n_blocks) & (rows["effect"] == top)]
            line = line.sort_values("n_samples")
            ax.errorbar(
                line["n_samples"], line["rate"], yerr=line["se"], color=color, linestyle=style["linestyle"],
                marker="o", markersize=9, markerfacecolor=color if n_blocks == 2 else SURFACE, markeredgewidth=2.0,
                capsize=4, elinewidth=1.5, label=style["label"],
            )[0].set_gid(f"power-{statistic}-{n_blocks}")
            labels.append((float(line["n_samples"].iloc[-1]), float(line["rate"].iloc[-1]), f"{n_blocks} blocks"))
        log_axis(ax, POWER_N, pad=1.3)
        ax.set_xlim(POWER_N[0] / 1.3, POWER_N[-1] * 2.4)
        end_labels(ax, labels, gap=0.07, dx=1.12, log=True, top=1.0)
        ax.set_ylim(0.0, 1.04)
        ax.set_title(f"{MODE_NAMES[mode]} → {STAT_SYMBOLS[statistic]} (e = {top:g})")
        ax.set_xlabel("samples (n)")
        value_grid(ax)
    axes[0].set_ylabel("power (rejection rate)")
    handles = [Line2D([], [], color=INK_SECONDARY, linestyle=s["linestyle"], marker="o", markersize=9,
                      markerfacecolor=INK_SECONDARY if b == 2 else SURFACE, markeredgewidth=2.0, label=s["label"])
               for b, s in BLOCK_STYLES.items()]
    fig.legend(handles=handles, loc="outside upper right", ncol=2)
    outputs.append(("power_vs_n", finish(fig, ctx, spec, uses)))
    for _, mode, statistic, effects in POWER_DIAGONALS:
        fig = new_figure()
        axes = fig.subplots(1, 2, sharey=True)
        color = STAT_COLORS[statistic]
        for ax, n_blocks in zip(axes, (2, 3)):
            labels = []
            for effect, marker in zip(effects, EFFECT_MARKERS):
                line = rows[(rows["mode"] == mode) & (rows["n_blocks"] == n_blocks) & (rows["effect"] == effect)]
                line = line.sort_values("n_samples")
                ax.errorbar(
                    line["n_samples"], line["rate"], yerr=line["se"], color=color, linestyle="-", marker=marker,
                    markersize=10, capsize=4, elinewidth=1.5, label=f"e = {effect:g}",
                )[0].set_gid(f"detail-{statistic}-{n_blocks}-{effect:g}")
                labels.append((float(line["n_samples"].iloc[-1]), float(line["rate"].iloc[-1]), f"e = {effect:g}"))
            log_axis(ax, POWER_N, pad=1.3)
            ax.set_xlim(POWER_N[0] / 1.3, POWER_N[-1] * 2.4)
            end_labels(ax, labels, gap=0.07, dx=1.12, log=True, top=1.0)
            ax.set_ylim(0.0, 1.04)
            ax.set_title(BLOCK_STYLES[n_blocks]["label"])
            ax.set_xlabel("samples (n)")
            value_grid(ax)
        axes[0].set_ylabel(f"{MODE_NAMES[mode]} mode: power of {STAT_LABELS[statistic]}")
        axes[0].legend(loc="upper left", title="group effect", title_fontsize=TYPE_SIZES["annotation"])
        outputs.append((f"power_vs_n_{statistic}", finish(fig, ctx, spec, uses)))
    return FigureResult(outputs, uses, rows)


# --------------------------------------------------------------------------- figure 7: power vs stage signal

STAGE_MODES: tuple[str, ...] = ("magnitude", "orientation", "shape")


def stage_signal_rows(root: Path) -> pd.DataFrame:
    """Per (δ, mode) at n = 200, 2 blocks: target power, median CV rank, median stage AUROC (anchor: rank, AUROC)."""

    rows = []
    for delta in (*STAGE_DELTAS, CORE_DELTA):
        if delta == CORE_DELTA:
            design: dict[str, Any] = {N_AXIS: 200, "n_blocks": 2}
            profile = "core"
        else:
            design = {axis: delta for axis in DELTA_AXES}
            profile = "stage_effect"
        for mode in (*STAGE_MODES, "none"):
            family = "magnitude" if mode == "magnitude" else "geometry"
            effect = 0.0 if mode == "none" else TOP_EFFECT[mode]
            statistic = MODE_TARGET.get(mode, "delta")
            cell = {**design, "trajectory_mode": mode, "effect_size": effect}
            rate = operating_rates(f"{profile}_{family}_operating", root, [{**cell, "statistic": statistic}])
            auroc = stage_signal(f"{profile}_{family}_stage_signal", root, [cell])
            rows.append(
                {
                    "delta": delta,
                    "mode": mode,
                    "statistic": statistic,
                    "effect": effect,
                    "rate": float(rate["rate"].iloc[0]),
                    "se": float(rate["se"].iloc[0]),
                    "median_rank": float(rate["median_rank"].iloc[0]),
                    "median_auroc": float(auroc["median_auroc"].iloc[0]),
                    "source": rate["source"].iloc[0],
                }
            )
    return pd.DataFrame(rows)


def build_power_vs_stage_signal(ctx: BuildContext) -> FigureResult:
    spec = FIGURES["power_vs_stage_signal"]
    rows = stage_signal_rows(ctx.root)
    keys = [f"{p}_{f}_{t}" for p in ("stage_effect", "core") for f in ("magnitude", "geometry")
            for t in ("operating", "stage_signal")]
    uses = [_use(key, ctx.root) for key in keys]
    fig = new_figure()
    axes = fig.subplots(1, 3, sharex=True)
    deltas = (*STAGE_DELTAS, CORE_DELTA)
    for ax, (column, title) in zip(axes, (
        ("rate", "power at the top effect"),
        ("median_rank", "median CV-selected rank"),
        ("median_auroc", "median held-out stage AUROC"),
    )):
        labels = []
        for mode in (*STAGE_MODES, "none"):
            if column == "rate" and mode == "none":
                continue
            line = rows[rows["mode"] == mode].sort_values("delta")
            color = mode_color(mode)
            kwargs: dict[str, Any] = {"color": color, "marker": "o", "markersize": 9, "markeredgewidth": 2.0}
            if mode == "none":
                kwargs.update(markerfacecolor=SURFACE, linewidth=2.0)
            if column == "rate":
                data_line = ax.errorbar(line["delta"], line[column], yerr=line["se"], capsize=4, elinewidth=1.5,
                                        **kwargs)[0]
            else:
                (data_line,) = ax.plot(line["delta"], line[column], **kwargs)
            data_line.set_gid(f"stage-{column}-{mode}")
            name = "zero effect" if mode == "none" else mode
            labels.append((float(line["delta"].iloc[-1]), float(line[column].iloc[-1]), name))
        log_axis(ax, deltas, labels=["0.25", "0.5", "1", "2"], pad=1.3)
        if column == "rate":
            ax.set_xlim(deltas[0] / 1.3, deltas[-1] * 2.9)
            ax.set_ylim(0.0, 1.04)
            end_labels(ax, labels, gap=0.07, dx=1.12, log=True, top=1.0)
        elif column == "median_rank":
            ax.set_ylim(0.0, float(rows["median_rank"].max()) + 2.0)
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        else:
            ax.set_ylim(0.5, 1.02)
        ax.set_title(title)
        ax.set_xlabel("joint stage effect δ")
        value_grid(ax)
    handles = [Line2D([], [], color=mode_color(m), marker="o", markersize=9, label=f"{m} (e = {TOP_EFFECT[m]:g})")
               for m in STAGE_MODES]
    handles.append(Line2D([], [], color=CONTROL_COLOR, linewidth=2.0, marker="o", markersize=9,
                          markerfacecolor=SURFACE, markeredgewidth=2.0, label="zero effect"))
    fig.legend(handles=handles, loc="outside upper right", ncol=4)
    note = ("n = 200 · 3 stages · 2 blocks; δ is each omic's stage effect (2 is the design-guidance default); "
            "top effects: magnitude e = 0.1, orientation and shape e = 1.")
    return FigureResult([("power_vs_stage_signal", finish(fig, ctx, spec, uses, note=note))], uses, rows)


# --------------------------------------------------------------------------- figure 8: estimator comparison

ESTIMATOR_STYLES: dict[str, dict[str, Any]] = {
    "pc1": {"phase": "power_primary", "linestyle": (0, (5, 3)), "filled": False, "label": "PC1"},
    "net_displacement": {"phase": "power_design", "linestyle": "-", "filled": True, "label": "net displacement"},
}
COMPARISON_EFFECTS: tuple[float, ...] = (0.0, 0.25, 0.5, 0.75, 1.0)


def estimator_rows(root: Path) -> pd.DataFrame:
    frames = []
    for design, key in (("phase6", "nd_power_phase6"), ("phase5", "nd_power_phase5")):
        for estimator, style in ESTIMATOR_STYLES.items():
            label = "PC1" if estimator == "pc1" else "net displacement"
            for mode in ("orientation", "shape"):
                cells = [
                    {ESTIMATOR_AXIS: estimator, "phase": style["phase"],
                     "trajectory_mode": "none" if effect == 0.0 else mode, "effect_size": effect, "statistic": "angle"}
                    for effect in COMPARISON_EFFECTS
                ]
                frame = operating_rates(key, root, cells, estimator=label)
                frame["design"], frame["estimator"], frame["mode"] = design, estimator, mode
                frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def build_estimator_comparison(ctx: BuildContext) -> FigureResult:
    spec = FIGURES["estimator_comparison"]
    rows = estimator_rows(ctx.root)
    uses = [_use("nd_power_phase6", ctx.root, label="n = 80 panel"),
            _use("nd_power_phase5", ctx.root, label="n = 1,200 panel")]
    uses += _floor_check_uses(ctx.root, uses)
    fig = new_figure()
    axes = fig.subplots(1, 2, sharey=True)
    titles = {"phase6": "n = 80 · SEA-AD table · 3 stages · 2 blocks", "phase5": "n = 1,200 · 4 stages · 3 blocks"}
    tags = {
        "phase6": resolve_measurement(SOURCES["nd_power_phase6"], ctx.root).short(),
        "phase5": resolve_measurement(SOURCES["nd_power_phase5"], ctx.root).short(),
    }
    short = {"pc1": "PC1", "net_displacement": "ND"}
    for ax, design in zip(axes, ("phase6", "phase5")):
        labels = []
        for mode in ("orientation", "shape"):
            color = mode_color(mode)
            for estimator, style in ESTIMATOR_STYLES.items():
                line = rows[(rows["design"] == design) & (rows["mode"] == mode) & (rows["estimator"] == estimator)]
                line = line.sort_values("effect_size")
                ax.plot(
                    line["effect_size"], line["rate"], color=color, linestyle=style["linestyle"], marker="o",
                    markersize=9, markerfacecolor=color if style["filled"] else SURFACE, markeredgewidth=2.0,
                    gid=f"estimator-{design}-{mode}-{estimator}",
                )
                end = float(line["rate"].iloc[-1])
                labels.append((1.0, end, f"{short[estimator]} {end:.3f}"))
        end_labels(ax, labels, gap=0.065, dx=0.05)
        ax.axhline(ALPHA, color=INK_MUTED, linewidth=1.2, zorder=1)
        ax.set_xlim(-0.04, 1.42)
        ax.set_xticks(list(COMPARISON_EFFECTS), [f"{e:g}" for e in COMPARISON_EFFECTS])
        ax.set_ylim(0.0, 1.04)
        ax.set_xlabel("group effect e")
        ax.set_title(titles[design] + (f"\n{tags[design]}" if tags[design] else ""), linespacing=1.4)
        value_grid(ax)
    axes[0].set_ylabel("θ (orientation) rejection rate")
    handles = [
        Line2D([], [], color=mode_color("orientation"), linewidth=2.5, label="orientation mode: target"),
        Line2D([], [], color=mode_color("shape"), linewidth=2.5, label="shape mode: leak into θ"),
        Line2D([], [], color=INK_SECONDARY, linestyle=ESTIMATOR_STYLES["pc1"]["linestyle"], marker="o",
               markerfacecolor=SURFACE, markeredgewidth=2.0, label="PC1 (dashed, open)"),
        Line2D([], [], color=INK_SECONDARY, marker="o", label="net displacement (solid, filled)"),
        Line2D([], [], color=INK_MUTED, linewidth=1.2, label="α = 0.05"),
    ]
    fig.legend(handles=handles, loc="outside upper center", ncol=3)
    return FigureResult([("estimator_comparison", finish(fig, ctx, spec, uses))], uses, rows)


# --------------------------------------------------------------------------- figure 9: covariates and the null

CORRECTION_STYLES: dict[str, dict[str, Any]] = {
    "uncorrected": {"linestyle": (0, (5, 3)), "filled": False, "offset": -0.07, "label": "uncorrected"},
    "df": {"linestyle": "-", "filled": True, "offset": 0.07, "label": "df-corrected"},
}
COVARIATE_COLUMNS: tuple[int, ...] = (1, 3, 6)


def build_covariates_null(ctx: BuildContext) -> FigureResult:
    spec = FIGURES["covariates_null"]
    readings = diagnosis_readings(ctx.root)
    pooled = pooled_diagnosis(readings)
    uses = [_use("diagnosis_delta0.1", ctx.root), _use("diagnosis_delta0.05", ctx.root)]
    color = STAT_COLORS["delta"]
    fig = new_figure()
    grid = fig.add_gridspec(1, 3, width_ratios=[0.25, 1.6, 0.25])
    for i in (0, 2):
        fig.add_subplot(grid[0, i]).axis("off")
    ax = fig.add_subplot(grid[0, 1])
    positions = {columns: index for index, columns in enumerate(COVARIATE_COLUMNS)}
    jitter = {(0.1, "anchor"): -0.03, (0.1, "translation"): -0.01, (0.05, "anchor"): 0.01, (0.05, "translation"): 0.03}
    for correction, style in CORRECTION_STYLES.items():
        part = readings[readings["correction"] == correction]
        xs = [positions[c] + style["offset"] + jitter[(s, k)] for c, s, k in
              zip(part["covariate_columns"], part["stage_effect"], part["control"])]
        ax.scatter(xs, part["rate"], s=55, facecolors=color if style["filled"] else SURFACE, edgecolors=color,
                   linewidths=1.5, alpha=0.4, zorder=2, gid=f"readings-{correction}")
        line = pooled[pooled["correction"] == correction].sort_values("covariate_columns")
        xs = [positions[c] + style["offset"] for c in line["covariate_columns"]]
        ax.plot(xs, line["rate"], color=color, linestyle=style["linestyle"], linewidth=3.0, marker="o",
                markersize=12, markerfacecolor=color if style["filled"] else SURFACE, markeredgewidth=2.5, zorder=3,
                gid=f"pooled-{correction}", label=f"{style['label']} (pooled over 4 readings)")
        for x, rate in zip(xs, line["rate"]):
            ax.annotate(f"{rate:.4f}".rstrip("0"), xy=(x, rate), xytext=(-12 if correction == "uncorrected" else 12, 0),
                        textcoords="offset points", ha="right" if correction == "uncorrected" else "left",
                        va="center", fontsize=TYPE_SIZES["annotation"], color=INK)
    bound = float(pooled["bound"].iloc[0])
    reference_line(ax, bound, f"pooled bound α + 2·SE = {bound:.4f}", color=INK_SECONDARY)
    reference_line(ax, ALPHA, "α = 0.05")
    ax.set_xticks(range(len(COVARIATE_COLUMNS)), [str(c) for c in COVARIATE_COLUMNS])
    ax.set_xlim(-0.5, len(COVARIATE_COLUMNS) - 0.5)
    ax.set_xlabel("covariate columns in the model (z + noise columns)")
    ax.set_ylabel("δ (size) rejection rate under the null")
    low = min(0.03, float(readings["rate"].min()) - 0.005)
    ax.set_ylim(low, max(0.08, float(readings["rate"].max()) + 0.005))
    value_grid(ax)
    ax.legend(loc="upper left")
    ax.set_title("SEA-AD astrocyte table · δ = 0.1 and 0.05 · anchor and translation", fontsize=TYPE_SIZES["label"])
    data = {"readings": readings, "pooled": pooled}
    return FigureResult([("covariates_null", finish(fig, ctx, spec, uses))], uses, data)


# --------------------------------------------------------------------------- registry of figures

_SPECS: tuple[FigureSpec, ...] = (
    FigureSpec(
        "mode_geometry",
        "What each trajectory mode does",
        "Each simulated mode changes group B's trajectory in one geometric way, read by one statistic: magnitude "
        "scales it (size δ), orientation rotates it (θ), shape bends it (shape), and translation moves it without "
        "changing any of the three.",
        (
            "Constructed coordinates, not data or a fitted space.",
            "Magnitude scales by 1.6 and orientation rotates by 35° about stage 1; shape moves the interior stage "
            "only, so it needs at least three stages.",
        ),
        build_mode_geometry,
    ),
    FigureSpec(
        "mode_simulated",
        "The four modes as the generator realizes them",
        "The numpy generator realizes each mode on a shared baseline: group B's stage-mean path is scaled, "
        "rotated, bent or shifted against group A's, and coincides with it at zero effect.",
        (
            "Display projection: per panel, the first two principal axes of the six group × stage means of the "
            "concatenated standardized blocks; translation's offset is orthogonal to the stage geometry, so its "
            "panel's second axis is the group offset. It is not the measurement space (the double-CV PLS latent "
            "space).",
            "One dataset per mode from generate_showcase_datasets (seed given in the manifest, n = 1,000, 3 stages, "
            "e = 1.0, surgery_censoring='clamp' as in the showcase); the modes share group A's baseline draw but are "
            "not otherwise identical.",
            "Axes are flipped so group A runs from stage 1 toward positive values (translation: B above A); the "
            "sign carries nothing.",
        ),
        build_mode_simulated,
    ),
    FigureSpec(
        "size_atlas",
        "Size holds from n = 80 to n = 1,501",
        "Every zero-effect anchor and translation control at every committed design, from n = 80 to n = 1,501 and "
        "with up to six covariate columns, rejects at or below its own α + k·SE bound.",
        (
            "Bounds are those the sources record: α + 2·SE at R = 500 (0.0695) except the corrected-permutation "
            "readings, judged per reading at α + 3·SE (0.0792).",
            "The harness tables record no verdict for translation; its bound is the anchor's form at its own R and "
            "its verdict is derived against it.",
            "Production arm per source: floor 2 (rank-floor check), df (corrected permutation), adjusted (ADNI), the "
            "net-displacement CV column (n = 1,200).",
            "Excluded: case-study-type-i (superseded by the df configuration) and covariate-type-i (pre-floor), whose "
            "uncorrected excess is what covariates_null shows; the latent-adjustment arms (a sensitivity analysis); "
            "the PC1 columns; the Phase 6 Type I profile (pre-floor, no equivalence check; the same table at the "
            "current measurement appears through the rank-floor and corrected-permutation rows).",
        ),
        build_size_atlas,
    ),
    FigureSpec(
        "specificity_n1200",
        "At n = 1,200 every statistic detects its mode",
        "At n = 1,200 each mode is detected by its target statistic; magnitude moves size alone, while rotating "
        "or bending a trajectory also changes its size (and rotation its shape).",
        (
            "Phase 5 net-displacement power run (2026-10-14), net-displacement column; rows at e = 1.0, the anchor at "
            "e = 0.",
            "Magnitude e = 1.0 here is not comparable with e = 0.1 at n = 80; the message is the off-diagonal pattern.",
        ),
        build_specificity_n1200,
    ),
    FigureSpec(
        "specificity_n80",
        "At n = 80 the statistics mix",
        "At n = 80 (the SEA-AD scale) only magnitude stays specific: orientation moves all three statistics "
        "about equally and shape moves size more than shape, so one significant statistic does not name the mode.",
        (
            "Design-guidance slice core profiles at n = 80, 3 stages, methylation + expression; magnitude at "
            "e = 0.1, orientation, shape and translation at e = 1.0, the anchor at e = 0.",
            "Same colour scale (0 to 1) as specificity_n1200.",
        ),
        build_specificity_n80,
    ),
    FigureSpec(
        "power_vs_n",
        "Power against sample size",
        "Power rises with n for every statistic and with the third block: size (δ) is well powered by n = 200 at "
        "e = 0.1, shape needs n ≈ 500 with two blocks, and orientation stays below 0.75 even at n = 500.",
        (
            "Design-guidance slice (2026-10-29), 3 stages, 500 × 999 per cell; ±1 Monte Carlo SE.",
            "Orientation responds little to e at this rank (2–3); see the detail figure.",
        ),
        build_power_vs_n,
    ),
    FigureSpec(
        "power_vs_stage_signal",
        "A weak stage signal inflates the rank",
        "As the joint stage effect weakens from δ = 2 to 0.25, the CV-selected rank inflates (median 2 to 11) and "
        "size and shape lose their power while the stage stays almost perfectly classified; orientation peaks at "
        "intermediate δ.",
        (
            "n = 200, 3 stages, two blocks; δ = 0.25–1.0 from the slice's stage-effect arm, δ = 2 from the slice "
            "core at the same design.",
            "Stage AUROC is the median held-out CV AUROC of the rank-selection fits; it stays near 1 in this range "
            "(the pilot's δ = 0.05 and 0.1 reach 0.55–0.73, bracketing SEA-AD).",
        ),
        build_power_vs_stage_signal,
    ),
    FigureSpec(
        "estimator_comparison",
        "Net displacement measures orientation better",
        "On identical datasets, net displacement doubles orientation power at n = 80 (0.300 → 0.648 at e = 1) and "
        "lifts it at n = 1,200 (0.840 → 0.924), while removing the shape mode's leak into θ.",
        (
            "Both designs predate the rank floor; the n = 1,200 panel carries the floor check's tag, the n = 80 panel "
            "stays pre-floor (the comparison is between estimators, not a current-measurement claim).",
            "Points at e = 0 are the shared zero-effect anchor.",
        ),
        build_estimator_comparison,
    ),
    FigureSpec(
        "covariates_null",
        "Covariates and the permutation null",
        "Each covariate column shrinks the residual df and inflates the uncorrected δ null (pooled 0.0545 at one "
        "column, 0.066 at six, above the 0.0597 bound); the df-corrected permutation holds it (0.0545 at six).",
        (
            "Case-study diagnosis (2026-10-30), covariate-blind basis (z_noise arms), δ = 0.1 and 0.05 × anchor and "
            "translation; pooled over 2,000 replicates per point.",
            "Faint points are the four individual readings (R = 500 each, bound 0.0695).",
        ),
        build_covariates_null,
    ),
)
FIGURES: dict[str, FigureSpec] = {spec.id: spec for spec in _SPECS}
FIGURE_IDS: tuple[str, ...] = tuple(FIGURES)


# --------------------------------------------------------------------------- writer, manifest and captions


class UnknownFigureError(PresentationFigureError):
    """Raised for a figure name that is not in the set."""


def validate_selection(names: Sequence[str] | None) -> list[str]:
    """The selected figure ids in set order; refuse an unknown name before anything is written."""

    if not names:
        return list(FIGURE_IDS)
    unknown = [name for name in names if name not in FIGURES]
    if unknown:
        raise UnknownFigureError(f"unknown figure(s) {unknown}; valid names: {', '.join(FIGURE_IDS)}.")
    return [figure for figure in FIGURE_IDS if figure in set(names)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_state(root: Path = REPO_ROOT) -> dict[str, Any]:
    """HEAD and whether the worktree is clean (tracked changes or untracked files make it dirty)."""

    try:
        revision = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(root), "status", "--porcelain"], check=True, capture_output=True, text=True
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return {"git_revision": None, "worktree_clean": None}
    return {"git_revision": revision, "worktree_clean": status.strip() == ""}


def _save_atomic(outputs: Sequence[tuple[str, Figure]], out_dir: Path, dpi: int) -> list[str]:
    """Write every output as PNG and PDF to temporary files, then move them all into place."""

    pending: list[tuple[Path, Path]] = []
    try:
        for stem, fig in outputs:
            for extension in ("png", "pdf"):
                final = out_dir / f"{stem}.{extension}"
                temporary = out_dir / f".{stem}.{extension}.tmp"
                pending.append((temporary, final))
                metadata = {"CreationDate": None} if extension == "pdf" else None
                fig.savefig(temporary, format=extension, dpi=dpi, facecolor=SURFACE, metadata=metadata)
    except BaseException:
        for temporary, _ in pending:
            temporary.unlink(missing_ok=True)
        raise
    for temporary, final in pending:
        os.replace(temporary, final)
    return [final.name for _, final in pending]


def _write_text_atomic(path: Path, text: str) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def _source_entries(uses: Sequence[Use], root: Path) -> list[dict[str, Any]]:
    entries: dict[str, dict[str, Any]] = {}
    for use in uses:
        base = REPO_ROOT if use.source.repo_file else Path(root)
        entry = entries.setdefault(
            use.source.path,
            {"path": use.source.path, "sha256": sha256(base / use.source.path), "run": use.source.run,
             "measurement": []},
        )
        tag = use.measurement.tag() if use.measurement is not None else _unmeasured(use.source)
        label = f"{use.label}: {tag}" if use.label else tag
        if label not in entry["measurement"]:
            entry["measurement"].append(label)
    for entry in entries.values():
        if len(entry["measurement"]) == 1:
            entry["measurement"] = entry["measurement"][0]
    return list(entries.values())


def _unmeasured(source: Source) -> str:
    if source.key == "floor_check":
        return "latent-space fits at floor 1 and floor 2 on rebuilt datasets; no RRPP"
    return "generator input (no measurement)"


@dataclass
class RenderReport:
    out_dir: Path
    manifest: dict[str, Any]
    failures: dict[str, str]


def render(
    out_dir: str | Path,
    figures: Sequence[str] | None = None,
    *,
    root: str | Path = REPO_ROOT,
    seed: int = 0,
    footer: bool = True,
    titles: bool = False,
    dpi: int = DPI,
    simulated_n: int = 1000,
) -> RenderReport:
    """Render the selected figures (all by default), then write the manifest and the captions."""

    selected = validate_selection(figures)
    state = git_state()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ctx = BuildContext(root=Path(root), seed=seed, footer=footer, titles=titles, simulated_n=simulated_n)
    entries: list[dict[str, Any]] = []
    failures: dict[str, str] = {}
    for figure_id in selected:
        spec = FIGURES[figure_id]
        result: FigureResult | None = None
        try:
            with plt.rc_context(style_rc()):
                result = spec.build(ctx)
                files = _save_atomic(result.outputs, out_dir, dpi)
            caveats = list(spec.caveats)
            if figure_id == "mode_simulated":
                caveats.append(f"Rendered with seed {seed} and n = {simulated_n:,} per dataset.")
            entries.append(
                {
                    "id": figure_id,
                    "title": spec.title,
                    "message": spec.message,
                    "files": files,
                    "sources": _source_entries(result.uses, ctx.root),
                    "caveats": caveats,
                }
            )
        except Exception as exc:  # one figure's failure must not touch another's files
            failures[figure_id] = f"{type(exc).__name__}: {exc}"
        finally:
            if result is not None:
                for _, fig in result.outputs:
                    plt.close(fig)
            plt.close("all")
    import motco

    manifest = {
        "set": {
            **state,
            "motco_version": motco.__version__,
            "matplotlib_version": matplotlib.__version__,
            "rendered_at": _dt.datetime.now().astimezone().isoformat(timespec="seconds"),
            "seed": seed,
            "footer": footer,
            "titles": titles,
            "dpi": dpi,
            "canvas_in": list(CANVAS),
            "palette": {"surface": SURFACE, "ink": [INK, INK_SECONDARY, INK_MUTED], "statistics": dict(STAT_COLORS),
                        "heatmap_ramp": list(HEATMAP_RAMP), "validation": PALETTE_VALIDATION},
        },
        "figures": entries,
        "failed": [{"id": key, "error": value} for key, value in failures.items()],
    }
    _write_text_atomic(out_dir / "manifest.json", json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    _write_text_atomic(out_dir / "captions.md", captions(manifest))
    return RenderReport(out_dir, manifest, failures)


def captions(manifest: dict[str, Any]) -> str:
    """Speaker notes: per figure, its title, message, measurement, caveats, files and sources."""

    state = manifest["set"]
    clean = {True: "clean worktree", False: "DIRTY worktree", None: "worktree state unknown"}[state["worktree_clean"]]
    lines = [
        "# Presentation figures: speaker notes",
        "",
        f"Rendered {state['rendered_at']} at {state['git_revision'] or 'an unknown revision'} ({clean}), "
        f"MOTCO {state['motco_version']}, matplotlib {state['matplotlib_version']}, seed {state['seed']}. "
        "Regenerate with `python scripts/presentation_figures.py --out-dir <dir>`; `manifest.json` holds the same "
        "content with every source's SHA-256.",
        "",
    ]
    for number, entry in enumerate(manifest["figures"], start=1):
        lines += [f"## {number}. {entry['title']} (`{entry['id']}`)", "", f"**Message.** {entry['message']}", ""]
        if not entry["sources"]:
            lines += ["**Measurement.** None: a constructed schematic, no data and no source file.", ""]
        else:
            lines.append("**Measurement and sources.**")
            lines.append("")
            for source in entry["sources"]:
                tags = source["measurement"] if isinstance(source["measurement"], list) else [source["measurement"]]
                for tag in tags:
                    lines.append(f"- {tag} (`{source['path']}`, {source['run']})")
            lines.append("")
        lines.append("**Caveats.**")
        lines.append("")
        lines += [f"- {caveat}" for caveat in entry["caveats"]]
        lines += ["", "**Files.** " + ", ".join(f"`{name}`" for name in entry["files"]), ""]
    if manifest["failed"]:
        lines += ["## Failed figures", ""]
        lines += [f"- `{failure['id']}`: {failure['error']}" for failure in manifest["failed"]]
        lines.append("")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="presentation_figures", description=(__doc__ or "").split("\n\n")[0])
    parser.add_argument("--out-dir", type=Path, required=True, help="Where to write the figures and the manifest.")
    parser.add_argument(
        "--figure", action="append", metavar="NAME",
        help=f"Render only this figure (repeatable; default all): {', '.join(FIGURE_IDS)}.",
    )
    parser.add_argument("--no-footer", action="store_true", help="Leave the measurement footer out of the figures.")
    parser.add_argument("--titles", action="store_true", help="Add each figure's title inside the figure.")
    parser.add_argument("--seed", type=int, default=0, help="Seed of the simulated mode figure (default 0).")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        validate_selection(args.figure)
    except UnknownFigureError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    report = render(args.out_dir, args.figure, seed=args.seed, footer=not args.no_footer, titles=args.titles)
    for entry in report.manifest["figures"]:
        print(f"  {entry['id']}: {', '.join(entry['files'])}")
    for figure_id, error in report.failures.items():
        print(f"error: {figure_id}: {error}", file=sys.stderr)
    clean = report.manifest["set"]["worktree_clean"]
    if clean is False:
        print("note: the worktree is dirty; the manifest records it.", file=sys.stderr)
    return 1 if report.failures else 0


if __name__ == "__main__":
    sys.exit(main())
