#!/usr/bin/env python3
"""Stage-effect sample-size scaling: whether a weak stage effect costs the trajectory test more than samples buy back.

The driver behind ``openspec/changes/stage-effect-sample-size-scaling/`` (spec
``stage-effect-scaling-validation``, roadmap open item 14). The generator's
stage program is δ times a fixed pattern over a noise covariance that does not
depend on δ, so a design's information about the program, and about a plant
defined relative to it, grows as δ√n. This driver reads designs with half the
stage effect and four times the samples against committed designs with the same
δ√n, cell by cell, under a rule fixed before any unit.

**The unit** is the univariate baseline's (``univariate_baseline_check.evaluate_unit``,
imported unchanged): one dataset, one cross-validated production fit (stage
PLS on methylation, expression and proteomics, 658 features, rank floor 2, net
displacement), the trajectory statistics with their null distributions, then
the split test, the interaction scan and the false-discovery reading at the
unit's RRPP seed, and ``pca_matched`` as the one recorded space. The unit is
handed a parent-shaped ``stage_span`` settings block (``unit_settings``: the
parent's settings with only the designs and cells replaced), so its records
equal the parent driver's for the same generator parameters, cell and seeds.
Judged arms: ``production`` (min p ≤ 0.05 / 3), ``fdr`` and ``total``; every
other arm is recorded and decides nothing.

**Profiles** (base seed 2043, the parent's; three stages, equal groups;
controls at 500 replicates, other cells at 200, 999 permutations):

- ``saturated`` — the parent's thirteen stage-span cells on five designs:
  the rerun ``n200_delta0.5`` (the parent's own design name and seeds, the
  reproduction gate) and four new designs, each matched in δ√n to a committed
  one: ``n320_delta1`` and ``n1280_delta0.5`` against ``n80`` (δ√n = 17.9),
  ``n800_delta1`` against ``n200`` (28.3), ``n800_delta0.25`` against
  ``n200_delta0.5`` (7.1). 16,000 units.
- ``band`` — designs below the AUROC ceiling, built from ``calibrate``'s δ_c:
  ``band_n320`` (δ_c, n = 320), ``band_n1280_half`` (δ_c / 2, n = 1,280, the
  matched design) and ``band_n1280`` (δ_c, n = 1,280). Eleven cells under the
  planted scale: the anchor, translation e = 1.0, magnitude e ∈ {0.5, 1, 2},
  in-span rotation and bend e ∈ {0.25, 0.5, 1}. 8,400 units.

``calibrate`` fits production, with no permutation, on 40 anchor datasets per
δ ∈ {0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1} at n = 320 on three blocks, from
the latent-space comparison's calibration seed stream
(``latent_space_comparison.calibration_seed`` at base seed 2043), and adopts the
δ whose median held-out AUROC is nearest that comparison's target, 0.636 (a tie
to the smaller δ). The median must lie in its band, 0.60–0.67, or the ``band``
settings are refused, naming the ladder's readings. The comparison's own
``calibration_unit`` is not used: it fits two blocks, and the band measures
three.

**Parent.** The univariate baseline's committed ``stage_span`` aggregates
(``results/univariate-baseline-2026-11-23/stage_span/``: ``settings.json``,
``rate_table.csv``, ``cell_table.csv``) are bound by sha256 in both profiles'
settings and a changed file is refused by name in ``run``, ``merge`` and
``decide``. The ``band`` settings also bind ``calibration/calibration.json``.

Subcommands (outputs under ``--out-dir``):

- ``calibrate`` — ``calibration/calibration.json`` (settings, per-δ AUROC and
  rank quantiles, the adopted δ) and ``calibration/calibration_units.csv``.
- ``run --profile saturated|band`` — (a shard of) the units into
  ``<profile>/records/shard_<i>_of_<N>.jsonl``; ``<profile>/settings.json``
  holds the settings and their hash, a resume or merge under other settings
  exits naming the differing keys, and a resume skips complete units.
- ``merge --profile …`` — the parent's tables for the profile's designs
  (``rate_table.csv``, ``paired_table.csv``, ``cell_table.csv``), a
  ``rank_table.csv`` (per design, cell and selected rank: each judged arm's and
  statistic's counts), ``accounting.json`` and, for ``saturated``,
  ``reproduction.json`` (the rerun design against the committed rate table,
  exact on available and rejected counts for every cell and arm, every
  differing row with both counts).
- ``decide`` — ``stage_effect_scaling_decision.json`` (rule version 1): the
  gate ``reproduction_failed``; size of ``production``, ``fdr`` and ``total``
  at the anchors and translations of the new designs in two families
  (``saturated``, m = 8; ``band``, m = 6) under the stage-span paper grade's
  count-plus-ceiling rule (``judge_rates``, imported unchanged),
  ``production_size_excess`` when production fails either; the 53 pair-cells
  (the eleven signal cells of each saturated pair, the reference counts read
  from the committed rate table, and the nine of the band pair
  ``band_n1280_half`` against ``band_n320``) labelled ``weak_costs`` |
  ``weak_gains`` | ``unresolved`` by Fisher's exact test, one-sided in each
  direction at 0.025 / 53 with a margin of 0.10; the same labels on ``fdr``,
  the control: a ``production`` label counts unless ``fdr`` carries the same
  label there; the verdict ``weak_stage_costs`` | ``weak_stage_gains`` |
  ``mixed`` | ``scaling_unrefuted`` by which counted labels occur. Beside it,
  never changing it: the parent's characterization grade per design (the four
  committed designs graded from the committed table), the parent's labels of
  ``production`` against the admissible per-feature arms over the new designs'
  71 signal cells (``univariate_baseline_check.label_cells``), and descriptive
  tables (the index of every design by δ, n, δ√n, held-out AUROC, rank and
  capture; every arm's rate per cell; the stretch rates against e·δ·√n; the
  rank strata of ``weak_costs`` pair-cells; the calibration's ladder).

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

    python scripts/stage_effect_scaling_check.py calibrate --out-dir /tmp/ses --workers 7
    # smoke
    python scripts/stage_effect_scaling_check.py run --profile saturated --reps 2 --control-reps 2 \\
        --out-dir /tmp/ses --workers 7
    python scripts/stage_effect_scaling_check.py run --profile band --reps 2 --control-reps 2 \\
        --out-dir /tmp/ses --workers 7
    python scripts/stage_effect_scaling_check.py merge --profile saturated --out-dir /tmp/ses
    python scripts/stage_effect_scaling_check.py merge --profile band --out-dir /tmp/ses
    python scripts/stage_effect_scaling_check.py decide --out-dir /tmp/ses

    # one SLURM array task (scripts/stage_effect_scaling_array.sbatch), then merge and decide
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
import interaction_scan_check as isc  # noqa: E402  (the size observation of one rate row)
import latent_space_comparison as lsc  # noqa: E402  (seeds, settings binding, calibration target and band, CSV)
import span_split_check as ssc  # noqa: E402  (parent loading, records, paths)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)
import univariate_baseline_check as ubc  # noqa: E402  (the unit, its arms, its tables, the labels and the grade)

sys.path.pop(0)

DRIVER = "stage_effect_scaling_check"
DRIVER_VERSION = 1
RULE_VERSION = 1
PROFILES: tuple[str, ...] = ("saturated", "band")

BASE_SEED = ubc.STAGE_SPAN_BASE_SEED
DEFAULT_PERMUTATIONS = ubc.DEFAULT_PERMUTATIONS
DEFAULT_REPLICATES: dict[str, int] = {"control": 500, "other": 200}

REPO = Path(__file__).resolve().parents[1]
PARENT_DIR = REPO / "results" / "univariate-baseline-2026-11-23" / "stage_span"
PARENT_FILE_NAMES: tuple[str, ...] = ("settings.json", "rate_table.csv", "cell_table.csv")
PARENT_SOURCE = "univariate baseline 2026-11-23, stage_span profile"

#: The parent's design rerun as the reproduction gate, under its own name and seeds.
RERUN_DESIGN = "n200_delta0.5"
#: design -> (sample size, every per-omic stage effect δ); the rerun first, then the new designs.
SATURATED_DESIGNS: dict[str, tuple[int, float]] = {
    RERUN_DESIGN: ubc.STAGE_SPAN_DESIGNS[RERUN_DESIGN],
    "n320_delta1": (320, 1.0),
    "n1280_delta0.5": (1280, 0.5),
    "n800_delta1": (800, 1.0),
    "n800_delta0.25": (800, 0.25),
}
#: (weaker-staged design, its committed match with the same δ√n).
SATURATED_PAIRS: tuple[tuple[str, str], ...] = (
    ("n320_delta1", "n80"),
    ("n1280_delta0.5", "n80"),
    ("n800_delta1", "n200"),
    ("n800_delta0.25", "n200_delta0.5"),
)
#: band design -> (sample size, its δ as a multiple of the calibrated δ_c).
BAND_DESIGNS: dict[str, tuple[int, float]] = {
    "band_n320": (320, 1.0),
    "band_n1280_half": (1280, 0.5),
    "band_n1280": (1280, 1.0),
}
BAND_PAIRS: tuple[tuple[str, str], ...] = (("band_n1280_half", "band_n320"),)
PLANTED = ubc.PLANTED
#: (cell, trajectory mode, effect, φ or None, support or None, span scale); every cell under the planted scale.
BAND_CELLS: tuple[tuple[str, str, float, float | None, int | None, str | None], ...] = (
    ("anchor", "none", 0.0, None, None, PLANTED),
    ("translation", "translation", 1.0, None, None, PLANTED),
    *((f"magnitude_e{effect:g}", "magnitude", effect, None, None, PLANTED) for effect in (0.5, 1.0, 2.0)),
    *((f"orient_e{effect:g}_phi0", "orientation", effect, 0.0, None, PLANTED) for effect in (0.25, 0.5, 1.0)),
    *((f"shape_e{effect:g}", "shape", effect, None, None, PLANTED) for effect in (0.25, 0.5, 1.0)),
)

CALIBRATION_LADDER: tuple[float, ...] = (0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1)
CALIBRATION_N = 320
CALIBRATION_REPLICATES = lsc.CALIBRATION_REPLICATES
CALIBRATION_TARGET = lsc.CALIBRATION_TARGET
CALIBRATION_BAND: tuple[float, float] = lsc.CALIBRATION_BAND
CALIBRATION_FILE = "calibration/calibration.json"
CALIBRATION_UNIT_COLUMNS = lsc.CALIBRATION_UNIT_COLUMNS
#: Distances to the target equal to this many decimals are a tie (taken by the smaller δ).
_TIE_DECIMALS = 12

CONTROL_CELLS: tuple[str, ...] = ubc.CONTROL_CELLS
JUDGED_ARMS: tuple[str, ...] = ("production", "fdr", "total")
PER_FEATURE_ARMS: tuple[str, ...] = ("fdr", "total")
CONTROL_ARM = "fdr"
LABEL_ALPHA = 0.025
MARGIN = ubc.MARGIN
_MARGIN_TOLERANCE = 1e-9
PAIR_CELLS_EXPECTED = 53
MAP_CELLS_EXPECTED = 71
PAIR_LABELS: tuple[str, ...] = ("weak_costs", "weak_gains", "unresolved")
STRATA_ARMS: tuple[str, ...] = ("production", *ubc.STATISTICS, "fdr", "total")
STRETCH_ARMS: tuple[str, ...] = ("production", *ubc.STATISTICS, "fdr")

VERDICTS: tuple[str, ...] = (
    "reproduction_failed", "production_size_excess",
    "weak_stage_costs", "weak_stage_gains", "mixed", "scaling_unrefuted",
)
_LIMIT = (
    " The reading holds at 658 simulated features with Gaussian noise, equal groups and three stages; the band "
    "designs match the case studies' held-out AUROC, not their cell tables, covariates or rank, and no real-noise "
    "reading exists on methylation, expression and proteomics (roadmap open item 2). No verdict changes what "
    "`motco analyze` runs."
)
CONSEQUENCES: dict[str, str] = {
    "reproduction_failed": (
        "Nothing is read: the rerun of n200_delta0.5 does not reproduce the committed counts, so no pair-cell "
        "compares with the committed table. Find the drift before any rerun; a mismatch is diagnosed, never waived."
    ),
    "production_size_excess": (
        "Nothing is read from the pairs: the production omnibus does not hold size at the failing controls of the "
        "new designs, so a cost or gain could not be told from excess rejections." + _LIMIT
    ),
    "weak_stage_costs": (
        "The documentation states that a weak stage effect is not bought back by sample size at the matched rate, "
        "with the plant families and the size of the loss, and guidance is indexed by the stage effect and the "
        "sample size separately. If the costs sit in the strata where the rank inflates, the rank rule under a "
        "fixed measured dimension is reopened (roadmap, On hold)." + _LIMIT
    ),
    "weak_stage_gains": (
        "The documentation states where a weakly staged design with a matched sample reads better than its "
        "match." + _LIMIT
    ),
    "mixed": "The documentation states the map of costs and gains by plant family." + _LIMIT,
    "scaling_unrefuted": (
        "The documentation may state that no cost beyond δ√n was resolved at these designs, with the margin of "
        "0.10: a difference below it is not excluded, so the absence of a label is not agreement. A design may be "
        "placed by n·(δ / 2)² as an approximate effective sample size." + _LIMIT
    ),
}
DECISION_NAME = "stage_effect_scaling_decision.json"


# ── Shared helpers ───────────────────────────────────────────────────────────


def profile_dir(out_dir: Path, profile: str) -> Path:
    return Path(out_dir) / profile


def _ladder_readings(steps: Iterable[Mapping[str, Any]]) -> str:
    return "; ".join(
        f"δ = {float(step['delta']):g}: median AUROC "
        + ("unavailable" if step.get("auroc_median") is None else f"{float(step['auroc_median']):.3f}")
        for step in steps
    )


def design_delta(unit_settings: Mapping[str, Any], design: str) -> float:
    return float(unit_settings["units"]["designs"][design]["generator"]["delta_methyl"])


def design_size(unit_settings: Mapping[str, Any], design: str) -> int:
    return int(unit_settings["units"]["designs"][design]["generator"]["n_samples"])


def delta_sqrt_n(unit_settings: Mapping[str, Any], design: str) -> float:
    return design_delta(unit_settings, design) * math.sqrt(design_size(unit_settings, design))


def restricted(unit_settings: Mapping[str, Any], designs: Sequence[str]) -> dict[str, Any]:
    """The unit settings read over ``designs`` only (the parent's tables, labels and grade follow ``designs``)."""

    return {**unit_settings, "designs": list(designs)}


# ── Settings ─────────────────────────────────────────────────────────────────


def unit_settings(
    designs: Mapping[str, tuple[int, float]],
    cells: Sequence[Sequence[Any]],
    *,
    replicates: int,
    control_replicates: int,
    permutations: int,
) -> dict[str, Any]:
    """The parent's ``stage_span`` settings with only the designs and cells replaced.

    Every other key (rules, arms, CV knobs, rank floor, estimator, recorded space, base seed and seed functions) is
    the one ``univariate_baseline_check.stage_span_settings`` writes, so ``evaluate_unit`` treats the units as the
    parent's.
    """

    settings = ubc.stage_span_settings(
        designs=[RERUN_DESIGN], replicates=replicates, control_replicates=control_replicates,
        permutations=permutations, base_seed=BASE_SEED, recorded_spaces=ubc.DEFAULT_RECORDED_SPACES,
    )
    settings["units"]["designs"] = {
        name: {"generator": ubc.stage_span_generator(n_samples, delta), "layers": None}
        for name, (n_samples, delta) in designs.items()
    }
    settings["units"]["cells"] = [list(cell) for cell in cells]
    settings["designs"] = list(designs)
    return settings


def _parent_binding(parent_dir: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """The parent's bound files and settings hash, and its committed settings (hash-checked)."""

    parent_dir = Path(parent_dir)
    committed = ssc.load_parent(parent_dir)
    binding = {
        "parent_dir": ssc._display_path(parent_dir),
        "parent_files": ubc._bind(parent_dir / name for name in PARENT_FILE_NAMES),
        "parent_settings_hash": committed["hash"],
        "parent_source": PARENT_SOURCE,
    }
    return binding, committed["settings"]


def _common(profile: str, binding: Mapping[str, Any], units: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "driver": DRIVER,
        "driver_version": DRIVER_VERSION,
        "profile": profile,
        **binding,
        "judged_arms": list(JUDGED_ARMS),
        "control_arm": CONTROL_ARM,
        "unit": "univariate_baseline_check.evaluate_unit on unit_settings (the parent's stage_span settings with "
                "the designs and cells replaced), pca_matched the one recorded space",
        "unit_settings": dict(units),
    }


def _check_designs(designs: Sequence[str], known: Mapping[str, Any], profile: str) -> list[str]:
    unknown = sorted(set(designs) - set(known))
    if unknown:
        raise SystemExit(f"Unknown {profile} design(s) {unknown}; expected a subset of {list(known)}.")
    return [name for name in known if name in designs]


def saturated_settings(
    *,
    replicates: int = DEFAULT_REPLICATES["other"],
    control_replicates: int = DEFAULT_REPLICATES["control"],
    permutations: int = DEFAULT_PERMUTATIONS,
    designs: Sequence[str] = tuple(SATURATED_DESIGNS),
    parent_dir: Path = PARENT_DIR,
) -> dict[str, Any]:
    """The rerun design and the four matched designs on the parent's thirteen cells, the parent bound."""

    ordered = _check_designs(designs, SATURATED_DESIGNS, "saturated")
    binding, committed = _parent_binding(parent_dir)
    if [list(cell) for cell in committed["units"]["cells"]] != [list(cell) for cell in ubc.STAGE_SPAN_CELLS]:
        raise SystemExit(f"{parent_dir}/settings.json has other cells than the parent driver's; the pairs would "
                         "not compare the same plants.")
    for _, match in SATURATED_PAIRS:
        if match not in committed["units"]["designs"]:
            raise SystemExit(f"{parent_dir}/settings.json has no design {match!r} to match.")
    units = unit_settings({name: SATURATED_DESIGNS[name] for name in ordered}, ubc.STAGE_SPAN_CELLS,
                          replicates=replicates, control_replicates=control_replicates, permutations=permutations)
    if RERUN_DESIGN in ordered and (units["units"]["designs"][RERUN_DESIGN]
                                    != committed["units"]["designs"][RERUN_DESIGN]):
        raise SystemExit(f"The rerun design {RERUN_DESIGN!r} differs from the parent's; it would not reproduce.")
    return {
        **_common("saturated", binding, units),
        "rerun_design": RERUN_DESIGN,
        "pairs": [[weaker, match] for weaker, match in SATURATED_PAIRS if weaker in ordered],
        "pair_reference": "the match's counts are read from the committed rate table",
    }


def load_calibration(path: Path) -> dict[str, Any]:
    """A calibration written by ``calibrate``, refused when missing, corrupt or adopting no δ."""

    path = Path(path)
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `calibrate` first. The band profile is built from it.")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("hash") != lsc.settings_hash(payload.get("settings", {})):
        raise SystemExit(f"{path} is corrupt: its hash does not match its settings.")
    if "adopted" not in payload:
        raise SystemExit(f"{path} adopts no δ.")
    return payload


def band_settings(
    out_dir: Path,
    *,
    replicates: int = DEFAULT_REPLICATES["other"],
    control_replicates: int = DEFAULT_REPLICATES["control"],
    permutations: int = DEFAULT_PERMUTATIONS,
    designs: Sequence[str] = tuple(BAND_DESIGNS),
    parent_dir: Path = PARENT_DIR,
) -> dict[str, Any]:
    """The three band designs from ``<out_dir>/calibration/calibration.json``; refused outside the band."""

    ordered = _check_designs(designs, BAND_DESIGNS, "band")
    path = Path(out_dir) / CALIBRATION_FILE
    payload = load_calibration(path)
    adopted = payload["adopted"]
    if not adopted.get("in_band"):
        low, high = payload["settings"]["band"]
        raise SystemExit(
            f"{path}: the ladder value nearest the target {payload['settings']['target']} has median AUROC "
            f"{adopted.get('auroc_median')}, outside the band {low}–{high}; refusing to build the band settings. "
            f"Ladder: {_ladder_readings(payload['steps'])}. Amend the ladder in notes.md before any band unit."
        )
    delta_c = float(adopted["delta"])
    binding, _ = _parent_binding(parent_dir)
    units = unit_settings(
        {name: (BAND_DESIGNS[name][0], BAND_DESIGNS[name][1] * delta_c) for name in ordered}, BAND_CELLS,
        replicates=replicates, control_replicates=control_replicates, permutations=permutations,
    )
    return {
        **_common("band", binding, units),
        "calibration": {
            "file": CALIBRATION_FILE,
            "sha256": ssc.sha256_file(path),
            "settings_hash": payload["hash"],
            "delta": delta_c,
            "auroc_median": adopted["auroc_median"],
            "target": payload["settings"]["target"],
            "band": list(payload["settings"]["band"]),
        },
        "pairs": [[weaker, match] for weaker, match in BAND_PAIRS if weaker in ordered and match in ordered],
        "pair_reference": "both designs of the pair are this profile's",
    }


def check_bound_files(settings: Mapping[str, Any], out_dir: Path) -> None:
    """Refuse, by name, a changed parent file or, for ``band``, a missing or changed calibration."""

    ubc.check_parent_files(settings)
    calibration = settings.get("calibration")
    if calibration is not None:
        path = Path(out_dir) / calibration["file"]
        if not path.exists():
            raise SystemExit(f"{path} does not exist; the band profile is bound to its calibration.")
        if ssc.sha256_file(path) != calibration["sha256"]:
            raise SystemExit(f"{path} has changed since the band settings were written (sha256 differs); refusing.")


# ── Calibration ──────────────────────────────────────────────────────────────


def calibration_settings(
    *, replicates: int = CALIBRATION_REPLICATES, ladder: Sequence[float] = CALIBRATION_LADDER
) -> dict[str, Any]:
    return {
        "driver": DRIVER,
        "driver_version": DRIVER_VERSION,
        "n_samples": CALIBRATION_N,
        "generator": ubc.stage_span_generator(CALIBRATION_N, float(ladder[0])),
        "anchor": list(ubc.STAGE_SPAN_CELLS[0]),
        "ladder": [float(delta) for delta in ladder],
        "replicates": int(replicates),
        "target": CALIBRATION_TARGET,
        "band": list(CALIBRATION_BAND),
        "base_seed": BASE_SEED,
        "seed_function": "latent_space_comparison.calibration_seed",
        "pls_cv_params": dict(lsc.PLS_CV_PARAMS),
        "min_components": lsc.MIN_COMPONENTS,
        "orientation_estimator": lsc.ORIENTATION_ESTIMATOR,
        "fit": "univariate_baseline_check.stage_span_params at zero permutations (the units' production fit, "
               "three blocks); no permutation",
        "rule": "the ladder value whose median held-out AUROC is nearest the target, a tie to the smaller δ; "
                "the median must lie in the band",
    }


def calibration_unit(delta: float, replicate: int, settings: Mapping[str, Any]) -> dict[str, Any]:
    """One anchor dataset at ``delta``: the production fit only, no permutation."""

    from motco.simulations.evaluation import integrate_semisynthetic_dataset
    from motco.simulations.reference import load_reference
    from motco.simulations.semisynthetic import generate_semisynthetic_trajectory

    start = perf_counter()
    seed = lsc.calibration_seed(int(settings["base_seed"]), int(replicate))
    entry: dict[str, Any] = {"delta": float(delta), "replicate": int(replicate), "generator_seed": seed}
    try:
        _, mode, effect, phi, support, scale = settings["anchor"]
        params = lsc._generator_params(
            settings["generator"], seed=seed, trajectory_mode=mode, group_effect_size=float(effect),
            span_alignment=phi, span_support=support, span_scale=scale,
            **{key: float(delta) for key in ubc.STAGE_EFFECT_KEYS},
        )
        dataset = generate_semisynthetic_trajectory(params, reference=load_reference())
        fit = {"units": {key: settings[key] for key in ("pls_cv_params", "min_components", "orientation_estimator")},
               "permutations": 0}
        latent = integrate_semisynthetic_dataset(dataset, ubc.stage_span_params(fit, seed=0))
        entry.update({"status": "completed", "auroc": float(latent.metadata["cv_mean_auroc"]),
                      "rank": int(latent.metadata["selected_lv"])})
    except Exception as exc:  # a failed fit is recorded, never fatal
        entry.update({"status": "failed", "error_type": type(exc).__name__, "error_message": str(exc)})
    entry["seconds"] = perf_counter() - start
    return entry


def summarize_step(delta: float, units: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """The comparison's step summary with the rank's 10 % and 90 % quantiles added."""

    ranks = [int(unit["rank"]) for unit in units if unit.get("status") == "completed"]
    return {**lsc.summarize_step(float(delta), units),
            "rank_q10": ubc._quantile(ranks, 0.1), "rank_q90": ubc._quantile(ranks, 0.9)}


def adopt_calibration(
    steps: Sequence[Mapping[str, Any]],
    *,
    band: Sequence[float] = CALIBRATION_BAND,
    target: float = CALIBRATION_TARGET,
) -> dict[str, Any]:
    """The step whose median AUROC is nearest ``target`` (a tie to the smaller δ), flagged in or out of ``band``."""

    ordered = sorted((step for step in steps if step.get("auroc_median") is not None), key=lambda s: float(s["delta"]))
    if not ordered:
        raise SystemExit("No calibration step produced an AUROC; nothing to adopt.")
    nearest = min(ordered, key=lambda s: round(abs(float(s["auroc_median"]) - float(target)), _TIE_DECIMALS))
    auroc = float(nearest["auroc_median"])
    return {
        "delta": float(nearest["delta"]),
        "auroc_median": auroc,
        "distance": abs(auroc - float(target)),
        "in_band": bool(float(band[0]) <= auroc <= float(band[1])),
        "rule": "nearest_target",
    }


def calibrate(out_dir: Path, settings: Mapping[str, Any], *, workers: int = 1) -> dict[str, Any]:
    """Every ladder step's fits into ``calibration/calibration.json`` and ``calibration_units.csv``."""

    path = Path(out_dir) / CALIBRATION_FILE
    if path.exists():
        raise SystemExit(f"{path} exists; a band profile may be bound to it. Use a fresh --out-dir.")
    path.parent.mkdir(parents=True, exist_ok=True)
    jobs = [(float(delta), replicate) for delta in settings["ladder"] for replicate in range(settings["replicates"])]
    if workers <= 1:
        units = [calibration_unit(delta, replicate, settings) for delta, replicate in jobs]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            units = list(pool.map(calibration_unit, *zip(*jobs), [settings] * len(jobs)))
    steps = [summarize_step(float(delta), [u for u in units if u["delta"] == float(delta)])
             for delta in settings["ladder"]]
    payload = {
        "settings": dict(settings),
        "hash": lsc.settings_hash(settings),
        "steps": steps,
        "adopted": adopt_calibration(steps, band=settings["band"], target=float(settings["target"])),
        "failures": sum(1 for u in units if u.get("status") != "completed"),
        "core_hours": sum(float(u["seconds"]) for u in units) / 3600.0,
    }
    lsc._write_json(path, payload)
    lsc._write_csv(path.parent / "calibration_units.csv", units, CALIBRATION_UNIT_COLUMNS)
    return payload


# ── Units and the run ────────────────────────────────────────────────────────


def profile_units(settings: Mapping[str, Any]) -> list[lsc.Unit]:
    """Design-major, then cell, then replicate; the controls at the control count (the parent's order)."""

    return ubc.profile_units(settings["unit_settings"])


def evaluate_unit(unit: lsc.Unit, units: Mapping[str, Any]) -> dict[str, Any]:
    """The parent's unit, unchanged, on the parent-shaped ``unit_settings``."""

    return ubc.evaluate_unit(unit, units)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    index: int,
    count: int,
    *,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units into its file; returns how many were evaluated."""

    profile = settings["profile"]
    check_bound_files(settings, out_dir)
    directory = profile_dir(out_dir, profile)
    lsc.ensure_settings(directory, settings)
    units_settings = settings["unit_settings"]
    path = lsc.shard_path(directory, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = ssc.read_records([path], ubc.record_key) if path.exists() else {}
    # Rewritten with complete lines only, so an interrupted unit is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(done):
            handle.write(lsc._line(done[key]))
    units = profile_units(settings)
    pending = [unit for position, unit in enumerate(units)
               if position % count == index and ubc.unit_key(unit) not in done]
    if limit is not None:
        pending = pending[:limit]

    def append(position: int, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{profile} {index}/{count}] {position}/{len(pending)} {ubc.record_key(record)} {record['status']}",
              file=sys.stderr, flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            append(position, evaluate_unit(unit, units_settings))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(evaluate_unit, unit, units_settings) for unit in pending]
            for position, future in enumerate(as_completed(futures), start=1):
                append(position, future.result())
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


RANK_COLUMNS: tuple[str, ...] = ("design", "cell", "rank", "arm", "count", "rejections")


def rank_rows(records: Sequence[Mapping[str, Any]], units: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Per design, cell and selected rank: each strata arm's available and rejected counts."""

    rows = []
    for fields, group in ubc._cell_groups(records, units):
        completed = [r for r in group if r.get("status") == "completed"]
        by_rank = ssc._grouped(completed, lambda r: (int(r["production"]["rank"]),))
        for (rank,), members in sorted(by_rank.items()):
            for arm in STRATA_ARMS:
                rate = ssc._rate_row([ubc._arm_rejection(r, arm) for r in members])
                rows.append({"design": fields["design"], "cell": fields["cell"], "rank": rank, "arm": arm,
                             "count": rate["count"], "rejections": rate["rejections"]})
    return rows


def rerun_reproduction(
    rates: Sequence[Mapping[str, Any]], committed: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]
) -> dict[str, Any]:
    """The rerun design's available and rejected counts against the committed rate table, every cell and arm."""

    base = {"design": RERUN_DESIGN, "source": PARENT_SOURCE, "files": dict(settings["parent_files"])}
    if RERUN_DESIGN not in settings["unit_settings"]["designs"]:
        return {**base, "passed": False, "rows_compared": 0, "differing": [], "unresolved_rows": [],
                "reason": f"the rerun design {RERUN_DESIGN!r} is not in this run"}
    mine = {(str(r["cell"]), str(r["arm"])): r for r in rates if str(r["design"]) == RERUN_DESIGN}
    theirs = {(str(r["cell"]), str(r["arm"])): r for r in committed if str(r["design"]) == RERUN_DESIGN}
    differing, unresolved = [], []
    compared = 0
    for key in [*mine, *(key for key in theirs if key not in mine)]:
        if key not in mine or key not in theirs:
            unresolved.append({"cell": key[0], "arm": key[1],
                               "missing_from": "this run" if key not in mine else "the committed table"})
            continue
        compared += 1
        wanted, got = ubc._counts(theirs[key]), ubc._counts(mine[key])
        if got != wanted:
            assert wanted is not None
            differing.append({"cell": key[0], "arm": key[1],
                              "committed": {"count": wanted[0], "rejections": wanted[1]},
                              "reproduced": None if got is None else {"count": got[0], "rejections": got[1]}})
    return {**base, "passed": compared > 0 and not differing and not unresolved, "rows_compared": compared,
            "differing": differing, "unresolved_rows": unresolved}


def merge(out_dir: Path, profile: str) -> dict[str, Any]:
    """Merge one profile's shards and write its tables, accounting and (``saturated``) reproduction check."""

    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    check_bound_files(settings, out_dir)
    units = settings["unit_settings"]
    unit_hash = lsc.settings_hash(units)
    expected_keys = [ubc.unit_key(unit) for unit in profile_units(settings)]
    raw: list[dict[str, Any]] = []
    for path in sorted((directory / "records").glob("shard_*.jsonl")):
        raw.extend(ssc.read_records([path], ubc.record_key).values())
    stale = sorted({ubc.record_key(r) for r in raw if r.get("settings_hash") != unit_hash})
    if stale:
        raise SystemExit(f"{len(stale)} record(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    for record in raw:
        seen[ubc.record_key(record)] = seen.get(ubc.record_key(record), 0) + 1
    by_key = {ubc.record_key(r): r for r in raw}
    records = [by_key[key] for key in expected_keys if key in by_key]
    with (directory / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    summary = {
        **ubc.accounting(records, units, expected=len(expected_keys),
                         duplicates=sum(count - 1 for count in seen.values())),
        "profile": profile,
        "settings_hash": lsc.settings_hash(settings),
        "unit_settings_hash": unit_hash,
    }
    rates = ubc.rate_rows(records, units)
    lsc._write_csv(directory / "rate_table.csv", rates, ubc.RATE_COLUMNS)
    lsc._write_csv(directory / "paired_table.csv", ubc.paired_rows(records, units), ubc.PAIRED_COLUMNS)
    lsc._write_csv(directory / "cell_table.csv", ubc.cell_rows(records, units), ubc.CELL_COLUMNS)
    lsc._write_csv(directory / "rank_table.csv", rank_rows(records, units), RANK_COLUMNS)
    lsc._write_json(directory / "accounting.json", summary)
    result = None
    if profile == "saturated":
        committed = lsc._read_csv(ssc._resolve(settings["parent_dir"]) / "rate_table.csv")
        result = rerun_reproduction(rates, committed, settings)
        lsc._write_json(directory / "reproduction.json", result)
    return {
        "profile": profile,
        **{key: summary[key] for key in ("complete_units", "missing_units", "failed_units", "recorded_core_hours")},
        "reproduction_passed": None if result is None else result["passed"],
        "reproduction_differing": None if result is None else len(result["differing"]),
    }


# ── Decision (pure functions) ────────────────────────────────────────────────


def _units(block: Mapping[str, Any]) -> dict[str, Any]:
    """A block's parent-shaped settings: ``unit_settings`` of a profile, the committed settings themselves."""

    settings = block["settings"]
    return settings["unit_settings"] if "unit_settings" in settings else settings


def new_designs(saturated: Mapping[str, Any]) -> list[str]:
    return [design for design in _units(saturated)["designs"] if design != RERUN_DESIGN]


def _rates(block: Mapping[str, Any]) -> dict[tuple, Mapping[str, Any]]:
    return ubc._index(block["rates"], "design", "cell", "arm")


def family_of(mode: str, alignment: Any) -> str:
    """The plant family: the trajectory mode, or ``out_of_span`` for a rotation off the stage plane."""

    return "out_of_span" if alignment is not None and float(alignment) > 0 else str(mode)


def size_reading(saturated: Mapping[str, Any], band: Mapping[str, Any]) -> dict[str, Any]:
    """Each judged arm's anchors and translations of the new designs in two ``judge_rates`` families."""

    families = {"saturated": (saturated, new_designs(saturated)), "band": (band, list(_units(band)["designs"]))}
    out: dict[str, Any] = {}
    for arm in JUDGED_ARMS:
        judged: dict[str, Any] = {}
        for family, (block, designs) in families.items():
            rates = _rates(block)
            judged[family] = grade.judge_rates([
                isc._size_observation(rates.get((design, cell, arm)), family=family, design=design, cell=cell,
                                      arm=arm)
                for design in designs
                for cell in CONTROL_CELLS
            ])
        out[arm] = {"admissible": all(family["holds"] for family in judged.values()), "families": judged}
    return out


def fisher_p_values(weaker: tuple[int, int], match: tuple[int, int]) -> dict[str, float]:
    """One-sided exact p-values on the 2 × 2 table of (rejected, not rejected) in the weaker design and its match.

    ``p_costs`` is small when the weaker-staged design rejects less often, ``p_gains`` when it rejects more often.
    """

    import numpy as np
    from scipy.stats import fisher_exact

    (n_weak, r_weak), (n_match, r_match) = weaker, match
    table = np.array([[r_weak, n_weak - r_weak], [r_match, n_match - r_match]], dtype=np.int64)
    return {"p_costs": float(fisher_exact(table, alternative="less").pvalue),
            "p_gains": float(fisher_exact(table, alternative="greater").pvalue)}


def exact_label(weaker: Mapping[str, Any] | None, match: Mapping[str, Any] | None, *, level: float) -> dict[str, Any]:
    """``weak_costs`` | ``weak_gains`` | ``unresolved`` for one arm of one pair-cell (rows of a rate table)."""

    weak_counts, match_counts = ubc._counts(weaker), ubc._counts(match)
    if weak_counts is None or match_counts is None or weak_counts[0] <= 0 or match_counts[0] <= 0:
        return {"label": "unresolved", "unavailable": True,
                "weaker": None if weak_counts is None else {"count": weak_counts[0], "rejections": weak_counts[1]},
                "match": None if match_counts is None else {"count": match_counts[0], "rejections": match_counts[1]}}
    p = fisher_p_values(weak_counts, match_counts)
    rate_weak, rate_match = weak_counts[1] / weak_counts[0], match_counts[1] / match_counts[0]
    difference = rate_weak - rate_match
    costs = p["p_costs"] <= level and -difference >= MARGIN - _MARGIN_TOLERANCE
    gains = p["p_gains"] <= level and difference >= MARGIN - _MARGIN_TOLERANCE
    return {
        "label": "weak_costs" if costs else ("weak_gains" if gains else "unresolved"),
        "unavailable": False,
        "weaker": {"count": weak_counts[0], "rejections": weak_counts[1], "rate": rate_weak},
        "match": {"count": match_counts[0], "rejections": match_counts[1], "rate": rate_match},
        "difference": difference,
        **p,
    }


def _pairs(saturated: Mapping[str, Any], band: Mapping[str, Any], committed: Mapping[str, Any]) -> list[tuple]:
    """(profile, weaker, match, weaker's block, match's block) for every pair, in order."""

    return [
        *(("saturated", weaker, match, saturated, committed) for weaker, match in saturated["settings"]["pairs"]),
        *(("band", weaker, match, band, band) for weaker, match in band["settings"]["pairs"]),
    ]


def pair_cell_labels(
    saturated: Mapping[str, Any], band: Mapping[str, Any], committed: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], float | None]:
    """Every pair-cell's ``production`` and ``fdr`` labels, and whether ``production``'s counts toward the verdict."""

    plan = []
    for profile, weaker, match, weak_block, match_block in _pairs(saturated, band, committed):
        units = _units(weak_block)
        for cell in ubc._cell_names(units):
            if cell not in CONTROL_CELLS:
                plan.append((profile, weaker, match, weak_block, match_block, cell))
    m = len(plan)
    level = LABEL_ALPHA / m if m else None
    entries = []
    for profile, weaker, match, weak_block, match_block, cell in plan:
        meta = ubc._cell_meta({"settings": _units(weak_block)}, weaker, cell)
        weak_rates, match_rates = _rates(weak_block), _rates(match_block)
        assert level is not None
        arms = {arm: exact_label(weak_rates.get((weaker, cell, arm)), match_rates.get((match, cell, arm)),
                                 level=level)
                for arm in ("production", CONTROL_ARM)}
        label, control = arms["production"]["label"], arms[CONTROL_ARM]["label"]
        entries.append({
            "pair": f"{weaker}~{match}", "profile": profile, "weaker": weaker, "match": match, "cell": cell, **meta,
            "family": family_of(meta["trajectory_mode"], meta["span_alignment"]),
            "label": label, "control_label": control,
            "counted": label != "unresolved" and control != label,
            "production": arms["production"], CONTROL_ARM: arms[CONTROL_ARM],
        })
    return entries, level


def design_grades(
    saturated: Mapping[str, Any], band: Mapping[str, Any], committed: Mapping[str, Any]
) -> list[dict[str, Any]]:
    """The parent's characterization grade per design, one ``judge_rates`` family each; committed designs too."""

    plan = [
        *((committed, design, "committed") for design in _units(committed)["designs"]),
        *((saturated, design, "this run, saturated profile") for design in new_designs(saturated)),
        *((band, design, "this run, band profile") for design in _units(band)["designs"]),
    ]
    out = []
    for block, design, source in plan:
        units = _units(block)
        graded = ubc.characterization_grade({"settings": restricted(units, [design]), "rates": block["rates"]})
        out.append({"design": design, "source": source, "delta": design_delta(units, design),
                    "n_samples": design_size(units, design), "delta_sqrt_n": delta_sqrt_n(units, design), **graded})
    return out


def map_reading(saturated: Mapping[str, Any], band: Mapping[str, Any], admissible: Sequence[str]) -> dict[str, Any]:
    """The parent's labels of ``production`` against the admissible per-feature arms on the new designs' cells."""

    views = {
        "stage_span": {"settings": restricted(_units(saturated), new_designs(saturated)),
                       "paired": saturated["paired"]},
        "surgery": {"settings": _units(band), "paired": band["paired"]},
    }
    profile_of = {"stage_span": "saturated", "surgery": "band"}
    m = sum(len(ubc.signal_cells(view["settings"])) for view in views.values())
    # The parent's rule reads its two profiles under their names; the two views stand in for them.
    labels = ubc.label_cells(views, list(admissible), m=m)
    tally: dict[tuple[str, str, str], dict[str, int]] = {}
    for entry in labels:
        entry["profile"] = profile_of[entry["profile"]]
        entry["family"] = family_of(entry["trajectory_mode"], entry["span_alignment"])
        counts = tally.setdefault((entry["profile"], entry["family"], entry["design"]),
                                  {label: 0 for label in ubc.LABELS})
        counts[entry["label"]] += 1
    out: dict[str, Any] = {
        "m_cells": m,
        "m_cells_expected": MAP_CELLS_EXPECTED,
        "admissible_arms": list(admissible),
        "levels": {
            "production_leads": ubc.LABEL_ALPHA / m if m else None,
            "baseline_leads": ubc.LABEL_ALPHA / (m * len(admissible)) if m and admissible else None,
        },
        "table": [{"profile": profile, "family": family, "design": design, **counts}
                  for (profile, family, design), counts in tally.items()],
        "labels": labels,
    }
    if not admissible:
        out["reason"] = "no admissible per-feature arm: neither fdr nor total holds size in both families"
    return out


def _index_row(block: Mapping[str, Any], design: str, source: str) -> dict[str, Any]:
    units = _units(block)
    anchor = ubc._index(block["cells"], "design", "cell").get((design, "anchor")) or {}
    return {
        "design": design, "source": source, "delta": design_delta(units, design),
        "n_samples": design_size(units, design), "delta_sqrt_n": delta_sqrt_n(units, design),
        **{key: ubc._finite(anchor.get(key)) for key in (
            "auroc_median", "rank_q10", "rank_median", "rank_q90", "capture_pooled_stage_median",
            "capture_group_difference_median", "completed")},
    }


def descriptive_tables(
    saturated: Mapping[str, Any],
    band: Mapping[str, Any],
    committed: Mapping[str, Any],
    calibration: Mapping[str, Any] | None,
    pair_cells: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """The index of every design, every arm's rate per cell, the stretch rates, the rank strata, the calibration."""

    this_run = {"saturated": saturated, "band": band}
    index = [_index_row(committed, design, "committed") for design in _units(committed)["designs"]]
    for profile, block in this_run.items():
        index.extend(_index_row(block, design, f"this run, {profile} profile")
                     for design in _units(block)["designs"])
    cell_rates = []
    for profile, block in this_run.items():
        units, rates = _units(block), _rates(block)
        for design in units["designs"]:
            for cell in ubc._cell_names(units):
                cell_rates.append({
                    "profile": profile, "design": design, "cell": cell,
                    **ubc._cell_meta({"settings": units}, design, cell),
                    "rates": {arm: ubc._rate_of(rates.get((design, cell, arm))) for arm in ubc.ARMS},
                })
    stretch = []
    sources = [(committed, design, "committed") for design in _units(committed)["designs"]]
    sources += [(saturated, design, "this run, saturated profile") for design in new_designs(saturated)]
    sources += [(band, design, "this run, band profile") for design in _units(band)["designs"]]
    for block, design, source in sources:
        units, rates = _units(block), _rates(block)
        for cell in ubc._cell_names(units):
            meta = ubc._cell_meta({"settings": units}, design, cell)
            if meta["trajectory_mode"] != "magnitude":
                continue
            effect = float(meta["effect_size"])
            stretch.append({
                "design": design, "cell": cell, "source": source, "effect_size": effect,
                "delta": design_delta(units, design), "n_samples": design_size(units, design),
                "e_delta_sqrt_n": effect * delta_sqrt_n(units, design),
                "rates": {arm: ubc._rate_of(rates.get((design, cell, arm))) for arm in STRETCH_ARMS},
            })
    strata = []
    for entry in pair_cells:
        if entry["label"] != "weak_costs":
            continue
        block = this_run[entry["profile"]]
        rows = [row for row in block.get("ranks", [])
                if (str(row["design"]), str(row["cell"])) == (entry["weaker"], entry["cell"])]
        by_rank: dict[int, dict[str, Any]] = {}
        for row in rows:
            count, rejections = int(row["count"]), int(row["rejections"])
            by_rank.setdefault(int(row["rank"]), {})[str(row["arm"])] = {
                "count": count, "rejections": rejections, "rate": rejections / count if count else None}
        strata.append({"pair": entry["pair"], "cell": entry["cell"], "design": entry["weaker"],
                       "counted": entry["counted"],
                       "ranks": [{"rank": rank, "arms": arms} for rank, arms in sorted(by_rank.items())]})
    ladder = None
    if calibration is not None:
        ladder = {"steps": calibration.get("steps"), "adopted": calibration.get("adopted"),
                  "failures": calibration.get("failures")}
    return {"index": index, "cell_rates": cell_rates, "stretch": stretch, "rank_strata": strata,
            "calibration": ladder}


def decision(
    saturated: Mapping[str, Any],
    band: Mapping[str, Any],
    committed: Mapping[str, Any],
    calibration: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Apply rule version 1.

    ``saturated`` and ``band`` hold ``settings``, ``accounting``, ``rates``, ``paired``, ``cells`` and ``ranks``
    (rows), ``saturated`` also ``reproduction``; ``committed`` holds the parent's ``settings``, ``rates`` and
    ``cells``. Every criterion is recorded whatever the verdict.
    """

    size = size_reading(saturated, band)
    admissible = [arm for arm in PER_FEATURE_ARMS if size[arm]["admissible"]]
    pair_cells, level = pair_cell_labels(saturated, band, committed)
    reproduction = saturated.get("reproduction") or {}
    costs = [entry for entry in pair_cells if entry["counted"] and entry["label"] == "weak_costs"]
    gains = [entry for entry in pair_cells if entry["counted"] and entry["label"] == "weak_gains"]
    failing: dict[str, Any] = {}
    if not reproduction.get("passed"):
        verdict = "reproduction_failed"
        failing["reproduction"] = {key: reproduction.get(key) for key in ("differing", "unresolved_rows", "reason")
                                   if reproduction.get(key)}
    elif not size["production"]["admissible"]:
        verdict = "production_size_excess"
        failing["production_size"] = {family: ssc._failing_observations(judged)
                                      for family, judged in size["production"]["families"].items()
                                      if not judged["holds"]}
    elif costs and gains:
        verdict = "mixed"
    elif costs:
        verdict = "weak_stage_costs"
    elif gains:
        verdict = "weak_stage_gains"
    else:
        verdict = "scaling_unrefuted"

    def keys(entries: Iterable[Mapping[str, Any]]) -> list[list[str]]:
        return [[entry["pair"], entry["cell"]] for entry in entries]

    blocks = {"saturated": saturated, "band": band}
    return {
        "rule_version": RULE_VERSION,
        "alpha": ubc.ALPHA,
        "omnibus_level": ubc.ALPHA / ubc.OMNIBUS_DIVISOR,
        "margin": MARGIN,
        "label_level": level,
        "m_pair_cells": len(pair_cells),
        "m_pair_cells_expected": PAIR_CELLS_EXPECTED,
        "control_arm": CONTROL_ARM,
        "judged_arms": list(JUDGED_ARMS),
        "settings_hashes": {name: lsc.settings_hash(block["settings"]) for name, block in blocks.items()},
        "parent_settings_hash": saturated["settings"].get("parent_settings_hash"),
        "calibration": band["settings"].get("calibration"),
        "complete": all(block["accounting"]["missing_units"] == 0 for block in blocks.values()),
        "failed_units": {name: block["accounting"]["failed_units"] for name, block in blocks.items()},
        "reproduction": {key: reproduction.get(key) for key in ("passed", "rows_compared")},
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "consequences": dict(CONSEQUENCES),
        "failing": failing,
        "named": {"weak_costs": keys(costs), "weak_gains": keys(gains)},
        "inadmissible_per_feature_arms": {
            arm: {family: ssc._failing_observations(judged) for family, judged in size[arm]["families"].items()
                  if not judged["holds"]}
            for arm in PER_FEATURE_ARMS if not size[arm]["admissible"]
        },
        "criteria": {
            "size": size,
            "pair_cells": pair_cells,
            "label_cells": {
                "counted_weak_costs": keys(costs),
                "counted_weak_gains": keys(gains),
                "shared_with_control": keys(e for e in pair_cells if e["label"] != "unresolved" and not e["counted"]),
                "unresolved": keys(e for e in pair_cells if e["label"] == "unresolved"),
                "unavailable": keys(e for e in pair_cells if e["production"]["unavailable"]),
            },
        },
        "grades": design_grades(saturated, band, committed),
        "map": map_reading(saturated, band, admissible),
        "descriptive": descriptive_tables(saturated, band, committed, calibration, pair_cells),
    }


def _profile_tables(out_dir: Path, profile: str) -> dict[str, Any]:
    directory = profile_dir(out_dir, profile)
    settings = lsc.load_settings(directory)
    check_bound_files(settings, out_dir)
    path = directory / "accounting.json"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge --profile {profile}` first.")
    summary = json.loads(path.read_text(encoding="utf-8"))
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{path} was merged under other settings; run `merge` again.")
    block: dict[str, Any] = {
        "settings": settings,
        "accounting": summary,
        "rates": ubc._read_csv(directory / "rate_table.csv"),
        "paired": ubc._read_csv(directory / "paired_table.csv"),
        "cells": ubc._read_csv(directory / "cell_table.csv"),
        "ranks": ubc._read_csv(directory / "rank_table.csv"),
    }
    if profile == "saturated":
        block["reproduction"] = json.loads((directory / "reproduction.json").read_text(encoding="utf-8"))
    return block


def committed_block(settings: Mapping[str, Any]) -> dict[str, Any]:
    """The parent's committed settings, rate table and cell table from the directory the settings bind."""

    directory = ssc._resolve(settings["parent_dir"])
    stored = ssc.load_parent(directory)
    if stored["hash"] != settings["parent_settings_hash"]:
        raise SystemExit(f"{directory}/settings.json is not the parent the settings were written against.")
    return {"settings": stored["settings"], "rates": ubc._read_csv(directory / "rate_table.csv"),
            "cells": ubc._read_csv(directory / "cell_table.csv")}


def decide(out_dir: Path) -> dict[str, Any]:
    """Read both merged profiles, the committed parent and the calibration; write the decision."""

    saturated = _profile_tables(out_dir, "saturated")
    band = _profile_tables(out_dir, "band")
    if saturated["settings"]["parent_files"] != band["settings"]["parent_files"]:
        raise SystemExit("The two profiles were written against different parent files; refusing.")
    calibration = load_calibration(Path(out_dir) / CALIBRATION_FILE)
    result = decision(saturated, band, committed_block(saturated["settings"]), calibration)
    lsc._write_json(Path(out_dir) / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    calibrate_parser = sub.add_parser("calibrate", help="Calibrate the band's δ_c to the held-out AUROC target")
    calibrate_parser.add_argument("--out-dir", type=Path, required=True)
    calibrate_parser.add_argument("--reps", type=int, default=CALIBRATION_REPLICATES,
                                  help=f"Anchor fits per ladder value (default {CALIBRATION_REPLICATES})")
    calibrate_parser.add_argument("--workers", type=int, default=1)
    run = sub.add_parser("run", help="Evaluate (a shard of) one profile's units")
    run.add_argument("--profile", choices=PROFILES, required=True)
    run.add_argument("--out-dir", type=Path, required=True)
    run.add_argument("--shard", type=str, default="0/1", help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1)
    run.add_argument("--perms", type=int, default=DEFAULT_PERMUTATIONS)
    run.add_argument("--reps", type=int, default=DEFAULT_REPLICATES["other"],
                     help=f"Replicates of every non-control cell (default {DEFAULT_REPLICATES['other']})")
    run.add_argument("--control-reps", type=int, default=DEFAULT_REPLICATES["control"],
                     help=f"Replicates of the anchor and translation (default {DEFAULT_REPLICATES['control']})")
    run.add_argument("--designs", type=str, default=None, help="Comma-separated designs; default: all")
    merge_parser = sub.add_parser("merge", help="Merge one profile's shards")
    merge_parser.add_argument("--profile", choices=PROFILES, required=True)
    merge_parser.add_argument("--out-dir", type=Path, required=True)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME} from both profiles")
    decide_parser.add_argument("--out-dir", type=Path, required=True)
    return parser


def settings_from_args(args: argparse.Namespace) -> dict[str, Any]:
    chosen = None if args.designs is None else [name for name in args.designs.split(",") if name]
    common = {"replicates": args.reps, "control_replicates": args.control_reps, "permutations": args.perms}
    if args.profile == "saturated":
        return saturated_settings(designs=chosen or tuple(SATURATED_DESIGNS), **common)
    return band_settings(args.out_dir, designs=chosen or tuple(BAND_DESIGNS), **common)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "calibrate":
        payload = calibrate(args.out_dir, calibration_settings(replicates=args.reps), workers=args.workers)
        print(json.dumps({"adopted": payload["adopted"], "failures": payload["failures"],
                          "steps": payload["steps"]}, indent=1))
        return 0
    if args.command == "run":
        index, count = lsc.parse_shard(args.shard)
        evaluated = run_shard(args.out_dir, settings_from_args(args), index, count, limit=args.limit,
                              workers=args.workers)
        print(f"evaluated {evaluated} {args.profile} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.profile), indent=1))
        return 0
    result = decide(args.out_dir)
    print(json.dumps({
        "verdict": result["verdict"],
        "reproduction": result["reproduction"],
        "admissible_per_feature_arms": result["map"]["admissible_arms"],
        "counted": {key: len(value) for key, value in result["criteria"]["label_cells"].items()},
        "grades": {entry["design"]: entry["grade"] for entry in result["grades"]},
        "complete": result["complete"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
