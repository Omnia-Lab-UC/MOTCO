#!/usr/bin/env python3
"""Latent-rank scaling probe for the orientation→shape response (numpy generator, no R).

The readiness item-2 follow-up behind ``unify-shape-reflection-policy``. The
geometry audit ruled out reflection as the mechanism by which an
``orientation``-mode surgery — a pure rotation of the group-B trajectory, which
changes no shape in the population — produces a rejected ``shape`` test in the
PLS latent space. What remains is the **rank-limited-projection** account: the
CV-selected latent space keeps only ~3 components, and a rotation that leaves
the population configuration unchanged need not leave its low-rank projection
unchanged.

The probe holds the generated data fixed (matched seeds) and varies only the
retained latent rank, via the diagnostic ``forced_components`` override. For
each rank it records:

- the observed group-vs-group Procrustes ``shape`` distance in that latent space,
- whether it clears that replicate's *own* RRPP shape null (rejection rate),
- the pooled relative eigengap of the latent stage-mean configuration
  (``config_spectrum``, the audit's resolvability covariate).

The reference line is the **population value**: the joint standardized
population ``shape`` distance for the same replicates, which is the response an
unlimited-rank measurement would return. **If the latent response decays toward
that population value as rank grows, the response is a projection artifact** —
the predeclaration the Phase-5 report contract needs.

The ``none`` mode is carried alongside as a null control: its latent shape
response should not depend on rank in the same way.

Runs serial (``n_jobs=1``) for the same reason as
``scripts/geometry_specificity_probe.py``: pool spawn overhead dwarfs the work
at these permutation counts.

Example (smoke):

    python scripts/latent_rank_probe.py --reps 3 --perms 19 --n-samples 60 \
        --ranks 3 4 --out-dir /tmp/latent_rank

Example (pilot design point):

    python scripts/latent_rank_probe.py --out-dir results/latent-rank-probe
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from statistics import fmean

import numpy as np

from motco.simulations.evaluation import (
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
)
from motco.simulations.reference import load_reference
from motco.simulations.semisynthetic import (
    SemiSyntheticTrajectoryParams,
    generate_semisynthetic_trajectory,
)

#: Rank ladder: 3 is the CV-selected production point, the rest span toward the
#: rank of the joint differential program. Clamped to feasibility per dataset.
DEFAULT_RANKS: tuple[int, ...] = (3, 4, 6, 9, 12)

#: Modes carried: the target and a null control.
DEFAULT_MODES: tuple[str, ...] = ("orientation", "none")

_CSV_COLUMNS = (
    "mode",
    "forced_components",
    "n_replicates",
    "observed_shape_mean",
    "observed_shape_sd",
    "rejection_rate",
    "pooled_eigengap_mean",
    "population_shape_mean",
    "excess_over_population",
)


@dataclass(frozen=True)
class RankCell:
    """Per-(mode, rank) summary over replicates."""

    mode: str
    forced_components: int
    n_replicates: int
    observed_shape_mean: float
    observed_shape_sd: float
    rejection_rate: float
    pooled_eigengap_mean: float
    population_shape_mean: float

    @property
    def excess_over_population(self) -> float:
        """Latent response minus the population value it should decay toward."""

        return self.observed_shape_mean - self.population_shape_mean

    def as_row(self) -> dict[str, object]:
        return {
            "mode": self.mode,
            "forced_components": self.forced_components,
            "n_replicates": self.n_replicates,
            "observed_shape_mean": self.observed_shape_mean,
            "observed_shape_sd": self.observed_shape_sd,
            "rejection_rate": self.rejection_rate,
            "pooled_eigengap_mean": self.pooled_eigengap_mean,
            "population_shape_mean": self.population_shape_mean,
            "excess_over_population": self.excess_over_population,
        }


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Latent-rank scaling probe (no R).")
    p.add_argument("--reps", type=int, default=20, help="Replicates per (mode, rank) (default: 20)")
    p.add_argument("--perms", type=int, default=199, help="RRPP permutations (default: 199)")
    p.add_argument("--n-samples", type=int, default=300, help="Samples per dataset (default: 300)")
    p.add_argument("--n-stages", type=int, default=4, help="Stages per trajectory (default: 4)")
    p.add_argument("--effect-size", type=float, default=1.0, help="Group effect size (default: 1.0)")
    p.add_argument("--p-dmp", type=float, default=0.2, help="Per-stage methylation DMP prob (default: 0.2)")
    p.add_argument("--alpha", type=float, default=0.05, help="Rejection threshold (default: 0.05)")
    p.add_argument("--base-seed", type=int, default=400, help="Base seed (default: 400)")
    p.add_argument(
        "--surgery-censoring",
        choices=("clamp", "error"),
        default="clamp",
        help=(
            "Pool-limited surgery policy. Defaults to 'clamp', matching the pilot config this "
            "probe re-measures; 'error' refuses a partial surgery instead (default: clamp)"
        ),
    )
    p.add_argument(
        "--ranks",
        type=int,
        nargs="+",
        default=list(DEFAULT_RANKS),
        help=f"Forced latent ranks (default: {' '.join(str(r) for r in DEFAULT_RANKS)})",
    )
    p.add_argument(
        "--modes",
        type=str,
        nargs="+",
        default=list(DEFAULT_MODES),
        help=f"Trajectory modes (default: {' '.join(DEFAULT_MODES)})",
    )
    p.add_argument(
        "--out-dir",
        type=str,
        default="results/latent-rank-probe",
        help="Directory for the CSV + markdown summary (default: results/latent-rank-probe)",
    )
    return p


def _population_joint_shape(result) -> float | None:
    """Joint standardized population ``shape`` — the unlimited-rank reference."""

    checkpoints = (result.realized_geometry or {}).get("checkpoints", {})
    joint = (checkpoints.get("population_standardized") or {}).get("joint")
    if not joint:
        return None
    value = joint.get("shape")
    return float(value) if value is not None and np.isfinite(value) else None


def _pooled_eigengap(result) -> float | None:
    pooled = (result.config_spectrum or {}).get("pooled") or {}
    value = pooled.get("relative_eigengap")
    return float(value) if value is not None and np.isfinite(value) else None


def run_rank_cell(
    mode: str,
    forced_components: int,
    *,
    n_replicates: int,
    n_samples: int,
    n_stages: int,
    effect_size: float,
    p_dmp: float,
    permutations: int,
    alpha: float,
    base_seed: int,
    surgery_censoring: str,
    reference,
) -> RankCell | None:
    """Evaluate one (mode, rank) cell over matched-seed replicates.

    Seeds are ``base_seed + rep`` regardless of rank, so every rank in the
    ladder measures the *same* generated datasets: the only thing that varies
    across a row of the table is how much of the latent space is retained.
    Returns ``None`` when the rank was infeasible for every replicate.
    """

    observed: list[float] = []
    eigengaps: list[float] = []
    population: list[float] = []
    rejections = 0

    for rep in range(n_replicates):
        params = SemiSyntheticTrajectoryParams(
            seed=base_seed + rep,
            trajectory_mode=mode,  # type: ignore[arg-type]
            n_samples=n_samples,
            n_stages=n_stages,
            group_effect_size=effect_size,
            p_dmp=p_dmp,
            surgery_censoring=surgery_censoring,  # type: ignore[arg-type]
            # The probe predates the joint default; keep its magnitude construction.
            magnitude_kind="all",
        )
        dataset = generate_semisynthetic_trajectory(params, reference=reference)
        try:
            result = evaluate_semisynthetic_trajectory(
                dataset,
                SimulationEvaluationParams(
                    integration_method="pls",
                    integration_params={"forced_components": forced_components},
                    permutations=permutations,
                    n_jobs=1,
                    seed=base_seed + rep,
                    orientation_estimator="pc1",  # committed probe measured PC1 angle
                ),
            )
        except Exception as exc:  # infeasible rank for this sample, or a failed fit
            print(f"      rep {rep}: skipped ({exc})", file=sys.stderr)
            continue

        shape = result.pair_statistics.get("shape")
        if shape is None or not np.isfinite(shape):
            continue
        observed.append(float(shape))
        gap = _pooled_eigengap(result)
        if gap is not None:
            eigengaps.append(gap)
        pop = _population_joint_shape(result)
        if pop is not None:
            population.append(pop)
        p = result.p_values.get("shape")
        if p is not None and np.isfinite(p) and p < alpha:
            rejections += 1

    if not observed:
        return None
    return RankCell(
        mode=mode,
        forced_components=forced_components,
        n_replicates=len(observed),
        observed_shape_mean=fmean(observed),
        observed_shape_sd=float(np.std(observed, ddof=1)) if len(observed) > 1 else 0.0,
        rejection_rate=rejections / len(observed),
        pooled_eigengap_mean=fmean(eigengaps) if eigengaps else float("nan"),
        population_shape_mean=fmean(population) if population else float("nan"),
    )


def describe_decay(cells: list[RankCell], mode: str) -> str:
    """State the decay measurement and its direction for one mode.

    Reported on the *observed* response rather than on ``|excess|``: the latent
    and population values live in differently scaled spaces, so the response is
    not expected to land exactly on the population value, and a response that
    decays through it would read as "growing away" on an absolute-gap measure.
    What the projection-artifact account predicts is a monotone decay with rank,
    which is what this states.
    """

    row = [cell for cell in sorted(cells, key=lambda c: c.forced_components) if cell.mode == mode]
    if len(row) < 2:
        return f"`{mode}`: too few feasible ranks to measure a trend."

    first, last = row[0], row[-1]
    start, end = first.observed_shape_mean, last.observed_shape_mean
    fraction = (end - start) / start if abs(start) > 1e-12 else float("nan")
    if end < start:
        direction = "decays"
    elif end > start:
        direction = "grows"
    else:
        direction = "is flat"
    monotone = all(b.observed_shape_mean <= a.observed_shape_mean for a, b in zip(row, row[1:]))
    shape_of_trend = "monotonically" if monotone else "non-monotonically"

    crossing = ""
    for previous, current in zip(row, row[1:]):
        gap_before, gap_after = previous.excess_over_population, current.excess_over_population
        if np.isfinite(gap_before) and np.isfinite(gap_after) and gap_before > 0 >= gap_after:
            crossing = (
                f" It crosses the population value {first.population_shape_mean:.4f} between rank "
                f"{previous.forced_components} and rank {current.forced_components}."
            )
            break

    return (
        f"`{mode}`: observed latent shape {direction} {shape_of_trend} with rank — "
        f"{start:.4f} at rank {first.forced_components} → {end:.4f} at rank "
        f"{last.forced_components} ({fraction:+.1%}). Rejection rate "
        f"{first.rejection_rate:.2f} → {last.rejection_rate:.2f}.{crossing}"
    )


def render_markdown(cells: list[RankCell], modes: list[str], args: argparse.Namespace) -> str:
    lines = [
        "# Latent-rank scaling probe",
        "",
        f"Run: {args.reps} reps/cell, {args.perms} perms, `n_samples={args.n_samples}`, "
        f"`n_stages={args.n_stages}`, `effect_size={args.effect_size}`, `p_dmp={args.p_dmp}`, "
        f"`alpha={args.alpha}`, `base_seed={args.base_seed}`, "
        f"`surgery_censoring={args.surgery_censoring}`, serial (`n_jobs=1`).",
        "",
        "Matched seeds: every rank in a row measures the same generated datasets; only the "
        "retained latent rank varies (`integration_params[\"forced_components\"]`).",
        "",
        "## Per-rank latent shape response",
        "",
        "| mode | rank | reps | observed shape (sd) | reject | pooled eigengap | population shape | excess |",
        "|------|------|------|---------------------|--------|-----------------|------------------|--------|",
    ]
    for cell in sorted(cells, key=lambda c: (c.mode, c.forced_components)):
        lines.append(
            f"| `{cell.mode}` | {cell.forced_components} | {cell.n_replicates} | "
            f"{cell.observed_shape_mean:.4f} ({cell.observed_shape_sd:.4f}) | "
            f"{cell.rejection_rate:.2f} | {cell.pooled_eigengap_mean:.3f} | "
            f"{cell.population_shape_mean:.4f} | {cell.excess_over_population:+.4f} |"
        )
    lines += [
        "",
        "`excess` is the latent response minus the joint standardized population value — "
        "the quantity that must shrink with rank if the response is a projection artifact.",
        "",
        "## Decay measurement",
        "",
    ]
    lines += [f"- {describe_decay(cells, mode)}" for mode in modes]
    lines += [
        "",
        "Reading: a **decaying** `orientation` response — falling toward (and possibly through) the "
        "population value as rank grows, with the rejection rate collapsing alongside it — supports "
        "the rank-limited-projection account: the orientation→shape response is an artifact of "
        "measuring a rotated trajectory in a rank-limited latent space, not a shape difference in "
        "the population. A response flat in rank does not support it. The `none` control shows how "
        "much of any decay is generic to rank rather than specific to the orientation surgery.",
        "",
        "Caveat: the population column is measured in the standardized observed space, not in the "
        "latent space, so the two are not on a common scale — the crossing rank is indicative, and "
        "the load-bearing evidence is the decay and the rejection-rate collapse, not exact "
        "agreement with the population value.",
        "",
    ]
    return "\n".join(lines)


def write_outputs(cells: list[RankCell], markdown: str, out_dir: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "latent_rank_probe.csv"
    md_path = out_dir / "latent_rank_probe.md"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(_CSV_COLUMNS))
        writer.writeheader()
        for cell in sorted(cells, key=lambda c: (c.mode, c.forced_components)):
            writer.writerow(cell.as_row())
    md_path.write_text(markdown, encoding="utf-8")
    return {"csv": csv_path, "markdown": md_path}


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    reference = load_reference()
    modes = list(dict.fromkeys(str(mode) for mode in args.modes))

    # Feasibility clamp: PLS cannot retain more components than
    # ``min(n_features, n_samples)``, and the override rejects rather than
    # clamps, so ranks past the bound are dropped from the ladder up front.
    max_rank = min(reference.n_cpg + reference.n_gene + reference.n_protein, args.n_samples)
    requested = sorted({int(rank) for rank in args.ranks})
    ranks = [rank for rank in requested if 2 <= rank <= max_rank]
    dropped = [rank for rank in requested if rank not in ranks]
    if dropped:
        print(
            f"Dropping infeasible rank(s) {dropped}: outside [2, {max_rank}] at this design point.",
            file=sys.stderr,
        )
    if not ranks:
        print(f"No feasible rank in {requested} for max rank {max_rank}.", file=sys.stderr)
        return 1

    cells: list[RankCell] = []
    total = len(modes) * len(ranks)
    step = 0
    for mode in modes:
        for rank in ranks:
            step += 1
            print(f"[{step}/{total}] {mode} @ rank {rank}...", file=sys.stderr)
            cell = run_rank_cell(
                mode,
                rank,
                n_replicates=args.reps,
                n_samples=args.n_samples,
                n_stages=args.n_stages,
                effect_size=args.effect_size,
                p_dmp=args.p_dmp,
                permutations=args.perms,
                alpha=args.alpha,
                base_seed=args.base_seed,
                surgery_censoring=args.surgery_censoring,
                reference=reference,
            )
            if cell is None:
                # Infeasible for every replicate: clamped out of the ladder.
                print(f"      rank {rank} infeasible at this design point", file=sys.stderr)
                continue
            cells.append(cell)

    if not cells:
        print("No feasible (mode, rank) cell produced a shape statistic.", file=sys.stderr)
        return 1

    markdown = render_markdown(cells, modes, args)
    paths = write_outputs(cells, markdown, Path(args.out_dir))
    print(f"Wrote {paths['csv']} and {paths['markdown']}", file=sys.stderr)
    print(markdown)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
