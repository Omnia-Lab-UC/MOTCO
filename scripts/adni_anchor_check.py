#!/usr/bin/env python3
"""The ADNI real-noise anchor: size of the trajectory tests on ADNI's blocks under permuted labels (no R).

The driver behind ``openspec/changes/adni-real-noise-anchor/`` (spec
``cohort-real-noise-validation``, roadmap open items 1 and 2). It reads the
``adni_sexdiff`` export in place, subsamples participants within sex × stage
cells into design points, permutes the sex labels within stage so that the
null is true on real noise, and reads the rejection rates of the production
scores and of ``predictive_orthonormal`` under a decision fixed before any unit
runs. No analysis runs under the real labels, and no package code changes.

**Contract** (design D2). A :class:`CohortContract` names the identifier, group
and stage columns, the two groups, the stage order, the blocks and their files,
the covariates with their kinds and levels, and the rank and AUROC the cohort's
fit must reproduce. :func:`load_cohort` reads the export through it and raises
``plasmode_check``'s named contract errors: the manifest is the ``export``
step's, records ``motco_transform: identity`` and the sha256 of the four
participant-level files; each file hashes to it; the metadata has its three
columns, no repeated participant, the two groups, the declared stages, no empty
cell and the manifest's participant count and cell table; each block lists the
metadata's participants in order, has the manifest's feature count and is
finite; the covariates are complete and at their declared levels.
:data:`ADNI_CONTRACT` is the first contract; a later cohort's driver imports
these functions with its own.

**Design points** (design D3): ``n80``, ``n200`` and ``n500`` draw that many
participants without replacement within each sex × stage cell, at the cohort's
proportions rounded by largest remainder (ties in the table's order);
``cohort`` holds all 1,501. 500 replicates × 999 permutations per point. A
subsample missing a level of a categorical covariate is redrawn from the next
attempt's seed (counted). Each unit permutes the sex labels within stage
(``plasmode_check.permute_within_stage``); a draw equal to the real labels is
redrawn (counted).

**Seeds** (design D4): sha256 of ``[driver name, 2039, point, replicate,
purpose(, attempt)]``, truncated as ``plasmode_check._seed`` does; purposes
``subsample``, ``labels``, ``rrpp``, ``folds``. Every arm of a unit takes the
unit's ``rrpp`` seed.

**Arms** (designs D5, D6), each one ``analyze_trajectories`` call at the case
study's configuration (five covariates, rank floor 2, net displacement, the
harness CV defaults) with the ``df`` residual correction:

- operating: ``pls_cv`` (the cross-validated production fit; at ``cohort``
  forced at the pre-flight's rank), ``predictive_orthonormal`` on the same fit,
  and both under within-stage draws (``__within_stage``);
- fixed rank: the same four combinations forced at ranks 8 and 19
  (``pls_rank8``, ``predictive_orthonormal_rank19__within_stage``, …); at
  ``cohort`` the rank-8 arms are the operating arms and are not run twice;
- control: ``pls_cv__uncorrected``, ``pls_cv`` without the correction.

A subsampled unit holds thirteen arms and makes one cross-validation; a cohort
unit holds nine and makes none.

**Checks on every unit** (design D7): ``rank`` (each forced arm reports its
rank), ``same_fit`` (per fit rank, every arm's model equals the first one's to
1e-8), ``measured_dimension`` (each predictive arm measures min(rank, 2)
dimensions), ``same_statistics`` (each within-stage arm and the control arm
observe exactly their counterpart's statistics) and ``cell_table``. A failing
unit is ``inconsistent`` and leaves every rate.

**Fit readings** (design D12), never decisive: per fit rank the in-sample
stage R² and the out-of-fold explained share on four stage-stratified folds,
and per coordinate system and fit rank the residual dispersion by stage, each
through the imported helpers with stage positions as labels.

**Parents** (sha256 in the settings; every subcommand refuses a changed file by
name): the within-stage run's committed SEA-AD anchor tables
``results/within-stage-permutation-2026-11-16/plasmode/rate_table.csv`` and
``cell_table.csv`` (the cross-cohort table). The export's five files are bound
the same way, and ``merge`` and ``decide`` re-hash them.

Subcommands (``--export-dir`` required on each, no default):

- ``preflight`` — writes ``settings.json`` and ``preflight.json``: on cohort
  replicates 0 and 1 the cross-validated fit must select the contract's rank and
  AUROC (ADNI: 8 and 0.599514) and equal the forced call exactly.
- ``run`` — evaluates (a shard of) the units into
  ``records/shard_<i>_of_<N>.jsonl``; round-robin shards, resume of complete
  units, refusal of other settings naming the keys; waits on a passed
  pre-flight.
- ``merge`` — ``merged.jsonl`` (keep off the repository), ``rate_table.csv``,
  ``paired_table.csv``, ``point_table.csv`` and ``accounting.json``.
- ``decide`` — ``adni_anchor_decision.json`` (rule version 1): the verdict
  ``gate_failed`` | ``production_excess`` | ``excess_persists`` |
  ``excess_scoped`` | ``excess_absent`` on the operating arms' size scopes
  (design D8, D9), the fixed-rank grade ``reference_exceeds`` | ``absent`` |
  ``reproduced`` | ``reproduced_unrestored`` (design D10), the scoping rule
  S1–S5 (design D11) and the descriptive readings (design D12).

Records hold seeds and statistics only: no participant identifier, label
vector, subsample membership or fold assignment. The export is never copied;
there is no sbatch array. With ``--workers`` above 1, pin BLAS in the shell.

Examples (the smoke: four units per design point):

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../adni_sexdiff/results/motco_input
    OUT=$HOME/motco-runs/adni-anchor-smoke
    python scripts/adni_anchor_check.py preflight --export-dir $EXPORT --out-dir $OUT --reps 4
    python scripts/adni_anchor_check.py run --export-dir $EXPORT --out-dir $OUT --reps 4 --workers 7
    python scripts/adni_anchor_check.py merge --export-dir $EXPORT --out-dir $OUT
    python scripts/adni_anchor_check.py decide --export-dir $EXPORT --out-dir $OUT
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter
from typing import Any, Iterable, Mapping, NamedTuple, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import latent_space_comparison as lsc  # noqa: E402  (settings hashing and files, JSON and CSV helpers)
import plasmode_check as pc  # noqa: E402  (contract errors, file hashing, the within-stage label permutation)
import predictive_size_excess_check as pse  # noqa: E402  (folds, out-of-fold indicators, dispersion by stage)
import predictive_subspace_check as psc  # noqa: E402  (arm records, stage R², model deviation, checks)
import production_latent_space_check as plsc  # noqa: E402  (parents, the exact paired tests)
import stage_span_paper_grade_decision as grade  # noqa: E402  (the count-plus-ceiling rule)

sys.path.pop(0)

DRIVER_NAME = "adni-real-noise-anchor"
DRIVER_VERSION = 1
RULE_VERSION = 1
BASE_SEED = 2039
PURPOSES: tuple[str, ...] = ("subsample", "labels", "rrpp", "folds")
MAX_ATTEMPTS = 100

#: Subsampled design points (participants); the full cohort is the last point.
POINT_SIZES: tuple[int, ...] = (80, 200, 500)
COHORT_POINT = "cohort"
REPLICATES = 500
PERMUTATIONS = 999
FIXED_RANKS: tuple[int, ...] = (8, 19)
PREFLIGHT_REPLICATES: tuple[int, ...] = (0, 1)

#: The case study's configuration (``adni_sexdiff``, MOTCO 0.10.0) with the ``df`` correction.
CV: dict[str, Any] = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 5, "max_components": 20, "random_state": 1203}
MIN_COMPONENTS = 2
ORIENTATION_ESTIMATOR = "net_displacement"
RESIDUAL_CORRECTION = "df"
MISSING_COVARIATES = "error"
FOLDS = 4

STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
ALPHA = 0.05
REJECTION = "p <= alpha"
CONSISTENCY_TOLERANCE = psc.CONSISTENCY_TOLERANCE
CONSISTENCY_CHECKS: tuple[str, ...] = ("rank", "same_fit", "measured_dimension", "same_statistics", "cell_table")

SCORES = "scores"
PREDICTIVE = "predictive_orthonormal"
COORDINATE_SYSTEMS: tuple[str, ...] = (SCORES, PREDICTIVE)
WITHIN_STAGE = "__within_stage"
UNCORRECTED = "__uncorrected"
STRATA = "stage"
FIT_CV = "cv"
REFERENCE_ARM = "pls_cv"
CANDIDATE_ARM = PREDICTIVE
CONTROL_ARM = f"{REFERENCE_ARM}{UNCORRECTED}"
OPERATING_ARMS: tuple[str, ...] = (
    REFERENCE_ARM, PREDICTIVE, f"{REFERENCE_ARM}{WITHIN_STAGE}", f"{PREDICTIVE}{WITHIN_STAGE}",
)
WITHIN_STAGE_OPERATING_ARMS: tuple[str, ...] = OPERATING_ARMS[2:]

REPO = Path(__file__).resolve().parents[1]
#: The within-stage run's committed SEA-AD anchor tables (cells ``anchor`` and ``anchor_fresh``).
PARENT_DIR = REPO / "results" / "within-stage-permutation-2026-11-16" / "plasmode"
SEAAD_DESIGN = "plasmode"
SEAAD_CELLS: tuple[str, ...] = ("anchor", "anchor_fresh")
#: SEA-AD astrocytes: 10,000 joint features on 80 donors (the plasmode export, ``plasmode-2026-11-02``).
SEAAD_FEATURES_PER_DONOR = 10_000 / 80

PREFLIGHT_NAME = "preflight.json"
DECISION_NAME = "adni_anchor_decision.json"


# ── Contract (design D2) ─────────────────────────────────────────────────────


@dataclass(frozen=True)
class CohortContract:
    """A cohort export's declared layout and the fit its cohort analysis must reproduce."""

    name: str
    sample_col: str
    group_col: str
    stage_col: str
    groups: tuple[str, str]
    stage_order: tuple[str, ...]
    metadata_file: str
    covariates_file: str
    #: (block name, file), in analysis order.
    blocks: tuple[tuple[str, str], ...]
    #: (covariate, ``numeric`` | ``categorical``), in model order.
    covariates: tuple[tuple[str, str], ...]
    #: (categorical covariate, its declared levels as strings).
    covariate_levels: tuple[tuple[str, tuple[str, ...]], ...]
    #: The cohort analysis's CV-selected rank and mean held-out stage AUROC, to ``fidelity_decimals``.
    fidelity_rank: int
    fidelity_auroc: float
    fidelity_decimals: int = 6
    fidelity_source: str = ""
    manifest_file: str = "manifest.json"
    manifest_step: str = "export"
    transform: str = "identity"

    @property
    def block_names(self) -> list[str]:
        return [name for name, _ in self.blocks]

    @property
    def covariate_names(self) -> list[str]:
        return [name for name, _ in self.covariates]

    @property
    def categorical_covariates(self) -> list[str]:
        return [name for name, kind in self.covariates if kind == "categorical"]

    @property
    def levels(self) -> dict[str, tuple[str, ...]]:
        return {name: tuple(levels) for name, levels in self.covariate_levels}

    @property
    def participant_files(self) -> list[str]:
        return [self.metadata_file, self.covariates_file, *(file for _, file in self.blocks)]

    @property
    def files_read(self) -> list[str]:
        return [self.manifest_file, *self.participant_files]

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "sample_col": self.sample_col,
            "group_col": self.group_col,
            "stage_col": self.stage_col,
            "groups": list(self.groups),
            "stage_order": list(self.stage_order),
            "metadata_file": self.metadata_file,
            "covariates_file": self.covariates_file,
            "blocks": [[name, file] for name, file in self.blocks],
            "covariates": [[name, kind] for name, kind in self.covariates],
            "covariate_levels": {name: list(levels) for name, levels in self.covariate_levels},
            "fidelity": {"selected_lv": self.fidelity_rank, "cv_mean_auroc": self.fidelity_auroc,
                         "decimals": self.fidelity_decimals, "source": self.fidelity_source},
            "manifest_file": self.manifest_file,
            "manifest_step": self.manifest_step,
            "transform": self.transform,
        }


ADNI_CONTRACT = CohortContract(
    name="adni",
    sample_col="rid",
    group_col="sex",
    stage_col="stage",
    groups=("Female", "Male"),
    stage_order=("CN", "MCI", "AD"),
    metadata_file="metadata.csv",
    covariates_file="covariates.csv",
    blocks=(("p180", "p180.csv"), ("nmr", "nmr.csv")),
    covariates=(("age", "numeric"), ("education", "numeric"), ("apoe4", "categorical"), ("phase", "categorical"),
                ("statin", "categorical")),
    covariate_levels=(("apoe4", ("0", "1", "2")), ("phase", ("ADNI1", "ADNI2GO")), ("statin", ("0", "1"))),
    fidelity_rank=8,
    fidelity_auroc=0.599514,
    fidelity_source="adni_sexdiff results/motco_analysis/result.json (MOTCO 0.10.0): selected_lv 8, cv_mean_auroc "
                    "0.5995142920769246",
)
CONTRACTS: dict[str, CohortContract] = {"adni": ADNI_CONTRACT}


@dataclass
class Cohort:
    """The loaded export: blocks and metadata in the metadata file's row order, with the files' sha256."""

    directory: Path
    contract: CohortContract
    manifest: dict[str, Any]
    blocks: dict[str, pd.DataFrame]
    metadata: pd.DataFrame
    files: dict[str, str]
    _cache: dict[str, Any] = field(default_factory=dict, repr=False)

    @property
    def stages(self) -> np.ndarray:
        return self.metadata[self.contract.stage_col].astype(str).to_numpy()

    @property
    def real_labels(self) -> np.ndarray:
        return self.metadata[self.contract.group_col].astype(str).to_numpy()

    @property
    def n_features(self) -> int:
        return int(sum(frame.shape[1] for frame in self.blocks.values()))

    def table(self) -> list[list[int]]:
        """The real group × stage table (rows: groups, columns: stages, in declared order)."""

        return cell_table(self.real_labels, self.stages, self.contract)


def cell_table(labels: Any, stages: Any, contract: CohortContract) -> list[list[int]]:
    labels, stages = np.asarray(labels).astype(str), np.asarray(stages).astype(str)
    return [[int(((labels == group) & (stages == stage)).sum()) for stage in contract.stage_order]
            for group in contract.groups]


def _require(condition: bool, error: type[pc.ExportContractError], message: str) -> None:
    if not condition:
        raise error(message)


def _repeated(values: Sequence[str]) -> list[str]:
    counts = pd.Series(list(values), dtype=object).value_counts()
    return sorted(str(value) for value in counts[counts > 1].index)


def _manifest_counts(manifest: Mapping[str, Any], path: Path) -> dict[str, Any]:
    counts = manifest.get("counts")
    _require(isinstance(counts, Mapping), pc.ManifestError, f"{path} records no counts.")
    return dict(counts)  # type: ignore[arg-type]


def load_cohort(directory: Path, contract: CohortContract = ADNI_CONTRACT) -> Cohort:
    """Read the export through the contract; raise the named error of the first violation. Writes nothing."""

    directory = Path(directory)
    manifest_path = directory / contract.manifest_file
    _require(manifest_path.is_file(), pc.ManifestError, f"{manifest_path} does not exist.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    _require(manifest.get("step") == contract.manifest_step, pc.ManifestError,
             f"{manifest_path} is the {manifest.get('step')!r} step's manifest, not {contract.manifest_step!r}'s.")
    transform = (manifest.get("settings") or {}).get("motco_transform")
    _require(transform == contract.transform, pc.ManifestError,
             f"{manifest_path} records motco_transform {transform!r}, not {contract.transform!r}.")
    outputs = manifest.get("outputs") or {}
    for name in contract.participant_files:
        entry = outputs.get(name)
        _require(isinstance(entry, Mapping) and bool(entry.get("sha256")) and entry.get("participant_level") is True,
                 pc.ManifestError, f"{manifest_path} does not list {name} as a participant-level file with its sha256.")
    for name in contract.participant_files:
        _require((directory / name).is_file(), pc.ManifestError, f"{directory / name} does not exist.")
    files = {name: pc.sha256_file(directory / name) for name in contract.files_read}
    for name in contract.participant_files:
        recorded = outputs[name]["sha256"]
        _require(files[name] == recorded, pc.ManifestError,
                 f"{directory / name} hashes to {files[name]}, the manifest records {recorded}.")
    counts = _manifest_counts(manifest, manifest_path)

    sample, group, stage = contract.sample_col, contract.group_col, contract.stage_col
    meta_path = directory / contract.metadata_file
    metadata = pd.read_csv(meta_path, dtype={sample: str, group: str, stage: str})
    _require(list(metadata.columns) == [sample, group, stage], pc.ManifestError,
             f"{meta_path} has columns {list(metadata.columns)}, not {[sample, group, stage]}.")
    ids = metadata[sample].astype(str).tolist()
    repeated = _repeated(ids)
    _require(not repeated, pc.DonorMismatchError, f"{meta_path} repeats {len(repeated)} participant(s).")
    observed = sorted(set(metadata[group].astype(str)))
    _require(observed == sorted(contract.groups), pc.GroupCountError,
             f"{meta_path} holds groups {observed}; the contract needs exactly {list(contract.groups)}.")
    unknown = sorted(set(metadata[stage].astype(str)) - set(contract.stage_order))
    _require(not unknown, pc.UnknownStageError,
             f"{meta_path} holds stage label(s) {unknown} outside the declared order {list(contract.stage_order)}.")
    table = cell_table(metadata[group], metadata[stage], contract)
    empty = [f"{g}/{s}" for i, g in enumerate(contract.groups) for j, s in enumerate(contract.stage_order)
             if table[i][j] == 0]
    _require(not empty, pc.EmptyCellError, f"{meta_path} has no participant in cell(s) {empty}.")
    _require(int(counts.get("participants", -1)) == len(ids), pc.ManifestError,
             f"{meta_path} lists {len(ids)} participants; the manifest records {counts.get('participants')}.")
    recorded_table = counts.get("sex_stage") or {}
    expected = [[int((recorded_table.get(g) or {}).get(s, -1)) for s in contract.stage_order] for g in contract.groups]
    _require(expected == table, pc.ManifestError,
             f"{meta_path}'s cell table {table} differs from the manifest's {expected}.")

    blocks: dict[str, pd.DataFrame] = {}
    feature_counts = counts.get("features") or {}
    for name, file in contract.blocks:
        path = directory / file
        frame = pd.read_csv(path, index_col=0)
        frame.index = frame.index.astype(str)
        _require(frame.index.tolist() == ids, pc.DonorMismatchError,
                 f"{path} does not list {contract.metadata_file}'s participants in the same order "
                 f"({len(set(frame.index) ^ set(ids))} differ in membership).")
        kept = (feature_counts.get(name) or {}).get("kept")
        _require(kept is not None and int(kept) == frame.shape[1], pc.ManifestError,
                 f"{path} has {frame.shape[1]} features; the manifest records {kept}.")
        try:
            values = frame.to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise pc.BlockValueError(f"{path} holds non-numeric values.") from exc
        _require(bool(np.isfinite(values).all()), pc.BlockValueError, f"{path} holds missing or infinite values.")
        blocks[name] = frame

    cov_path = directory / contract.covariates_file
    covariates = pd.read_csv(cov_path, dtype={sample: str})
    _require(sample in covariates.columns, pc.ManifestError, f"{cov_path} lacks {sample!r}.")
    cov_ids = covariates[sample].astype(str).tolist()
    repeated = _repeated(cov_ids)
    _require(not repeated, pc.DonorMismatchError, f"{cov_path} repeats {len(repeated)} participant(s).")
    _require(set(cov_ids) == set(ids), pc.DonorMismatchError,
             f"{cov_path} and {contract.metadata_file} list different participants "
             f"({len(set(cov_ids) ^ set(ids))} differ).")
    absent = [name for name in contract.covariate_names if name not in covariates.columns]
    _require(not absent, pc.MissingCovariateError, f"{cov_path} lacks covariate column(s) {absent}.")
    levels = contract.levels
    for name, kind in contract.covariates:
        missing = int(covariates[name].isna().sum())
        _require(not missing, pc.MissingCovariateError, f"{cov_path}: covariate {name!r} is missing for "
                                                         f"{missing} participant(s).")
        if kind == "categorical":
            outside = sorted(set(covariates[name].astype(str)) - set(levels[name]))
            _require(not outside, pc.CovariateValueError,
                     f"{cov_path}: covariate {name!r} holds level(s) {outside} outside {list(levels[name])}.")
        else:
            numeric = pd.to_numeric(covariates[name], errors="coerce").to_numpy(dtype=float)
            _require(bool(np.isfinite(numeric).all()), pc.CovariateValueError,
                     f"{cov_path}: covariate {name!r} holds non-numeric or non-finite values.")
    merged = metadata.merge(covariates[[sample, *contract.covariate_names]], on=sample, how="inner",
                            validate="one_to_one")
    merged = merged[[sample, group, stage, *contract.covariate_names]].reset_index(drop=True)
    return Cohort(directory=directory, contract=contract, manifest=manifest, blocks=blocks, metadata=merged,
                  files=files)


def contract_exit(exc: pc.ExportContractError) -> SystemExit:
    return pc.contract_exit(exc)


def verify_export_files(directory: Path, settings: Mapping[str, Any]) -> None:
    """Re-hash the export's files named in the settings; refuse a difference, naming every file."""

    recorded = settings["export"]["files"]
    differing = []
    for name, digest in sorted(recorded.items()):
        path = Path(directory) / name
        if not path.is_file() or pc.sha256_file(path) != digest:
            differing.append(name)
    if differing:
        raise SystemExit(f"The export's file hashes differ from the settings' for: {', '.join(differing)}; refusing.")


def verify(settings: Mapping[str, Any], cohort: Cohort | None = None) -> None:
    """Refuse a changed parent file, or an export whose files differ from the settings', by name."""

    plsc.verify_parents(settings)
    if cohort is not None:
        plsc.verify_export(cohort, settings)


# ── Design points (design D3) ────────────────────────────────────────────────


def point_table(table: Sequence[Sequence[int]], n: int) -> list[list[int]]:
    """The table scaled to ``n`` and rounded by largest remainder, ties broken in the table's (row-major) order."""

    cells = [int(value) for row in table for value in row]
    total = sum(cells)
    if not 0 < int(n) <= total:
        raise ValueError(f"A design point of {n} participants needs 0 < n <= {total}.")
    floors = [value * int(n) // total for value in cells]
    remainders = [value * int(n) % total for value in cells]
    order = sorted(range(len(cells)), key=lambda index: (-remainders[index], index))
    for index in order[: int(n) - sum(floors)]:
        floors[index] += 1
    width = len(table[0])
    return [floors[row * width:(row + 1) * width] for row in range(len(table))]


def point_name(size: int) -> str:
    return f"n{int(size)}"


def design_points(cohort_table: Sequence[Sequence[int]], sizes: Sequence[int] = POINT_SIZES) -> list[dict[str, Any]]:
    total = int(sum(sum(row) for row in cohort_table))
    if list(sizes) != sorted(set(int(size) for size in sizes)) or any(int(size) >= total for size in sizes):
        raise ValueError(f"Point sizes must ascend strictly below the cohort's {total} participants: {list(sizes)}.")
    points = [
        {"name": point_name(size), "participants": int(size), "table": point_table(cohort_table, int(size)),
         "subsampled": True}
        for size in sizes
    ]
    points.append({"name": COHORT_POINT, "participants": total, "table": [list(map(int, row)) for row in cohort_table],
                   "subsampled": False})
    return points


def point_names(settings: Mapping[str, Any]) -> list[str]:
    return [str(point["name"]) for point in settings["points"]]


def point_spec(settings: Mapping[str, Any], name: str) -> dict[str, Any]:
    for point in settings["points"]:
        if point["name"] == name:
            return dict(point)
    raise KeyError(f"Unknown design point {name!r}.")


# ── Arms (designs D5, D6) ────────────────────────────────────────────────────


def fit_label(rank: int | None) -> str:
    return FIT_CV if rank is None else f"rank{int(rank)}"


def arm_name(
    coordinates: str, fit: str, draws: str = "unrestricted", correction: str | None = RESIDUAL_CORRECTION
) -> str:
    stem = (REFERENCE_ARM if coordinates == SCORES else PREDICTIVE) if fit == FIT_CV else (
        f"pls_{fit}" if coordinates == SCORES else f"{PREDICTIVE}_{fit}"
    )
    return stem + (WITHIN_STAGE if draws == "within_stage" else "") + ("" if correction else UNCORRECTED)


def arm_specs(fixed_ranks: Sequence[int] = FIXED_RANKS) -> dict[str, dict[str, Any]]:
    """Every arm's fit, coordinates, draws, correction and family, in evaluation order."""

    specs: dict[str, dict[str, Any]] = {}
    for fit, rank in ((FIT_CV, None), *((fit_label(rank), int(rank)) for rank in fixed_ranks)):
        family = "operating" if fit == FIT_CV else "fixed_rank"
        for coordinates in COORDINATE_SYSTEMS:
            for draws in ("unrestricted", "within_stage"):
                specs[arm_name(coordinates, fit, draws)] = {
                    "fit": fit, "fixed_rank": rank, "coordinates": coordinates, "draws": draws,
                    "residual_correction": RESIDUAL_CORRECTION, "family": family,
                }
    specs[CONTROL_ARM] = {"fit": FIT_CV, "fixed_rank": None, "coordinates": SCORES, "draws": "unrestricted",
                          "residual_correction": None, "family": "control"}
    return specs


def point_arms(settings: Mapping[str, Any], point: str) -> list[str]:
    """The arms a unit of the point holds: all at a subsampled point; at the cohort, none at the cohort's rank."""

    specs = settings["arm_specs"]
    if point_spec(settings, point)["subsampled"]:
        return list(specs)
    rank = int(settings["cohort_rank"])
    return [arm for arm, spec in specs.items() if spec["fixed_rank"] != rank]


def counterpart(arm: str) -> str | None:
    """The arm whose observed statistics this one must equal: the unrestricted, corrected arm of its fit."""

    if arm.endswith(WITHIN_STAGE):
        return arm.removesuffix(WITHIN_STAGE)
    if arm.endswith(UNCORRECTED):
        return arm.removesuffix(UNCORRECTED)
    return None


def scores_reference(arm: str, specs: Mapping[str, Mapping[str, Any]]) -> str | None:
    """A predictive arm's scores arm at the same fit and null."""

    spec = specs[arm]
    if spec["coordinates"] != PREDICTIVE:
        return None
    return arm_name(SCORES, spec["fit"], spec["draws"], spec["residual_correction"])


def paired_arms(settings: Mapping[str, Any]) -> list[tuple[str, str]]:
    specs = settings["arm_specs"]
    pairs = [(arm, counterpart(arm)) for arm in specs if counterpart(arm) is not None]
    pairs += [(arm, scores_reference(arm, specs)) for arm in specs if scores_reference(arm, specs) is not None]
    return [(arm, str(reference)) for arm, reference in pairs]


# ── Settings ─────────────────────────────────────────────────────────────────


def _seed(*parts: Any) -> int:
    encoded = json.dumps([DRIVER_NAME, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def build_settings(
    cohort: Cohort,
    *,
    replicates: int = REPLICATES,
    permutations: int = PERMUTATIONS,
    point_sizes: Sequence[int] | None = None,
    fixed_ranks: Sequence[int] | None = None,
    cv: Mapping[str, Any] | None = None,
    parent_dir: Path | None = None,
) -> dict[str, Any]:
    """The run's settings; unset arguments take the module's values at call time (ADNI's design)."""

    contract = cohort.contract
    point_sizes = POINT_SIZES if point_sizes is None else tuple(point_sizes)
    fixed_ranks = FIXED_RANKS if fixed_ranks is None else tuple(fixed_ranks)
    if int(replicates) < max(PREFLIGHT_REPLICATES) + 1:
        raise SystemExit(f"--reps must be at least {max(PREFLIGHT_REPLICATES) + 1} (the pre-flight's replicates).")
    if int(permutations) < 1:
        raise SystemExit("--perms must be at least 1.")
    specs = arm_specs(fixed_ranks)
    parent_dir = Path(PARENT_DIR if parent_dir is None else parent_dir)
    return {
        "driver": DRIVER_NAME,
        "driver_version": DRIVER_VERSION,
        "rule_version": RULE_VERSION,
        "contract": contract.to_dict(),
        "export": {"files": dict(cohort.files), "participants": int(len(cohort.metadata)),
                   "features": {name: int(frame.shape[1]) for name, frame in cohort.blocks.items()}},
        "parents": {
            "seaad_anchor_rates": plsc._parent(parent_dir / "rate_table.csv"),
            "seaad_anchor_cells": plsc._parent(parent_dir / "cell_table.csv"),
        },
        "base_seed": BASE_SEED,
        "seed_derivation": f"sha256 of [{DRIVER_NAME!r}, {BASE_SEED}, point, replicate, purpose(, attempt)], truncated "
                           "as plasmode_check._seed; purposes subsample, labels, rrpp, folds",
        "max_attempts": MAX_ATTEMPTS,
        "points": design_points(cohort.table(), point_sizes),
        "replicates": int(replicates),
        "permutations": int(permutations),
        "fixed_ranks": [int(rank) for rank in fixed_ranks],
        "cohort_rank": int(contract.fidelity_rank),
        "arms": list(specs),
        "arm_specs": specs,
        "analysis": {
            "cv": dict(CV if cv is None else cv),
            "min_components": MIN_COMPONENTS,
            "orientation_estimator": ORIENTATION_ESTIMATOR,
            "residual_correction": RESIDUAL_CORRECTION,
            "missing_covariates": MISSING_COVARIATES,
            "permutation_strata": STRATA,
            "latent_coordinates": PREDICTIVE,
            "n_jobs": 1,
        },
        "fit": "a subsampled unit's pls_cv cross-validates once and every other operating arm is forced at its rank; "
               "a cohort unit forces every arm at cohort_rank; fixed-rank arms are forced at their rank",
        "folds": FOLDS,
        "alpha": ALPHA,
        "rejection": REJECTION,
        "consistency_tolerance": CONSISTENCY_TOLERANCE,
        "preflight_replicates": list(PREFLIGHT_REPLICATES),
    }


class Unit(NamedTuple):
    point: str
    replicate: int


def unit_key(point: str, replicate: int) -> str:
    return f"{point}/{int(replicate)}"


def record_key(record: Mapping[str, Any]) -> str:
    return unit_key(record["point"], record["replicate"])


def enumerate_units(settings: Mapping[str, Any]) -> list[Unit]:
    """Replicate-major, so a partial run covers every point."""

    return [Unit(point, replicate) for replicate in range(int(settings["replicates"])) for point in
            point_names(settings)]


# ── Units: seeds, subsample and labels (designs D3, D4) ──────────────────────


def unit_seeds(settings: Mapping[str, Any], point: str, replicate: int) -> dict[str, int]:
    """The unit's first-attempt seeds (a redraw takes the next attempt's)."""

    base = int(settings["base_seed"])
    seeds = {}
    if point_spec(settings, point)["subsampled"]:
        seeds["subsample"] = _seed(base, point, int(replicate), "subsample", 0)
    seeds["labels"] = _seed(base, point, int(replicate), "labels", 0)
    seeds["rrpp"] = _seed(base, point, int(replicate), "rrpp")
    seeds["folds"] = _seed(base, point, int(replicate), "folds")
    return seeds


def _levels_present(frame: pd.DataFrame, contract: CohortContract) -> bool:
    levels = contract.levels
    return all(set(frame[name].astype(str)) == set(levels[name]) for name in contract.categorical_covariates)


def draw_subsample(
    cohort: Cohort, settings: Mapping[str, Any], point: str, replicate: int
) -> tuple[np.ndarray, int]:
    """Row indices into the cohort, in its row order, and how many draws were refused for an absent level."""

    contract = cohort.contract
    table = point_spec(settings, point)["table"]
    labels, stages = cohort.real_labels, cohort.stages
    base = int(settings["base_seed"])
    for attempt in range(int(settings["max_attempts"])):
        rng = np.random.default_rng(_seed(base, point, int(replicate), "subsample", attempt))
        chosen = []
        for g, group in enumerate(contract.groups):
            for s, stage in enumerate(contract.stage_order):
                rows = np.flatnonzero((labels == group) & (stages == stage))
                size = int(table[g][s])
                if size > rows.size:
                    raise ValueError(f"Cell {group}/{stage} holds {rows.size} participants; the point needs {size}.")
                chosen.append(rng.choice(rows, size=size, replace=False))
        rows = np.sort(np.concatenate(chosen))
        if _levels_present(cohort.metadata.iloc[rows], contract):
            return rows, attempt
    raise RuntimeError(f"{settings['max_attempts']} subsamples all missed a declared covariate level.")


def draw_labels(
    real: np.ndarray, stages: np.ndarray, settings: Mapping[str, Any], contract: CohortContract, point: str,
    replicate: int,
) -> tuple[np.ndarray, int]:
    """A within-stage permutation of ``real`` that differs from it, and how many draws were refused."""

    base = int(settings["base_seed"])
    for attempt in range(int(settings["max_attempts"])):
        rng = np.random.default_rng(_seed(base, point, int(replicate), "labels", attempt))
        drawn = pc.permute_within_stage(real, stages, contract.stage_order, rng)
        if not np.array_equal(drawn, np.asarray(real, dtype=object)):
            return drawn, attempt
    raise RuntimeError(f"{settings['max_attempts']} within-stage permutations all reproduced the real labels.")


def build_unit(cohort: Cohort, settings: Mapping[str, Any], point: str, replicate: int) -> dict[str, Any]:
    """The unit's sample (blocks and covariates on their rows) under its permuted labels; in memory only."""

    contract = cohort.contract
    subsampled = bool(point_spec(settings, point)["subsampled"])
    if subsampled:
        rows, subsample_redraws = draw_subsample(cohort, settings, point, replicate)
    else:
        rows, subsample_redraws = np.arange(len(cohort.metadata)), 0
    metadata = cohort.metadata.iloc[rows].reset_index(drop=True)
    real = metadata[contract.group_col].astype(str).to_numpy()
    stages = metadata[contract.stage_col].astype(str).to_numpy()
    labels, label_redraws = draw_labels(real, stages, settings, contract, point, replicate)
    metadata[contract.group_col] = np.asarray(labels).astype(str)
    blocks = {name: frame.iloc[rows] for name, frame in cohort.blocks.items()}
    return {
        "point": point,
        "replicate": int(replicate),
        "seeds": unit_seeds(settings, point, replicate),
        "subsample_redraws": int(subsample_redraws),
        "label_redraws": int(label_redraws),
        "rows": rows,
        "metadata": metadata,
        "blocks": blocks,
        "table": cell_table(metadata[contract.group_col], stages, contract),
    }


def make_dataset(blocks: Mapping[str, pd.DataFrame], metadata: pd.DataFrame, contract: CohortContract) -> Any:
    from motco.analysis import OmicsBlock, TrajectoryDataset

    return TrajectoryDataset(
        blocks=[OmicsBlock(name, blocks[name], contract.transform) for name in contract.block_names],
        metadata=metadata,
        sample_col=contract.sample_col,
        group_col=contract.group_col,
        stage_col=contract.stage_col,
        stage_order=list(contract.stage_order),
        group_order=list(contract.groups),
        covariate_cols=contract.covariate_names,
        categorical_covariates=contract.categorical_covariates,
        missing_covariates=MISSING_COVARIATES,
    )


def analysis_params(
    settings: Mapping[str, Any], seed: int, *, forced: int | None = None, coordinates: str = SCORES,
    draws: str = "unrestricted", correction: str | None = RESIDUAL_CORRECTION, permutations: int | None = None,
) -> Any:
    """The case-study parameters with the arm's settings; nothing else differs between arms."""

    from motco.analysis import TrajectoryAnalysisParams

    analysis = settings["analysis"]
    return TrajectoryAnalysisParams(
        permutations=int(settings["permutations"] if permutations is None else permutations),
        seed=int(seed),
        n_jobs=1,
        orientation_estimator=analysis["orientation_estimator"],
        min_components=int(analysis["min_components"]),
        residual_correction=correction,
        forced_components=None if forced is None else int(forced),
        latent_coordinates=None if coordinates == SCORES else coordinates,
        permutation_strata=STRATA if draws == "within_stage" else None,
        **analysis["cv"],
    )


def arm_params(arm: str, settings: Mapping[str, Any], seed: int, *, subsampled: bool, production_rank: int) -> Any:
    spec = settings["arm_specs"][arm]
    if spec["fit"] == FIT_CV:
        forced = None if (arm == REFERENCE_ARM and subsampled) else int(production_rank)
    else:
        forced = int(spec["fixed_rank"])
    return analysis_params(settings, seed, forced=forced, coordinates=spec["coordinates"], draws=spec["draws"],
                           correction=spec["residual_correction"])


def _analyse(dataset: Any, params: Any) -> Any:
    from motco.analysis import analyze_trajectories

    return analyze_trajectories(dataset, params)


def arm_record(arm: str, result: Any, seconds: float, spec: Mapping[str, Any], fit: str) -> dict[str, Any]:
    record = psc.analysis_arm_record(result, seconds)
    record.pop("span_capture", None)
    record.pop("population_latent", None)
    return {**record, **{key: spec[key] for key in ("coordinates", "draws", "residual_correction")}, "fit": fit,
            "fit_rank": spec["fit"]}


def failed_arm(arm: str, exc: Exception, seconds: float, spec: Mapping[str, Any]) -> dict[str, Any]:
    return {**psc._failed_arm(exc, seconds), **{key: spec[key] for key in ("coordinates", "draws",
                                                                            "residual_correction")},
            "fit_rank": spec["fit"]}


# ── Per-unit checks (design D7) ──────────────────────────────────────────────


def _statistics_deviation(left: Mapping[str, Any], right: Mapping[str, Any]) -> float:
    worst = 0.0
    for statistic in STATISTICS:
        a, b = left.get(statistic), right.get(statistic)
        if a is None and b is None:
            continue
        if a is None or b is None:
            return math.inf
        worst = max(worst, abs(float(a) - float(b)))
    return worst


def _named_check(deviations: Mapping[str, float], tolerance: float) -> dict[str, Any]:
    worst = max(deviations.values(), default=None)
    check = psc._check(worst, tolerance, applicable=bool(deviations))
    failing = {arm: value for arm, value in deviations.items()
               if not (math.isfinite(float(value)) and float(value) <= tolerance)}
    return {**check, "failing_arms": failing}


def check_rank(arms: Mapping[str, Mapping[str, Any]], forced: Mapping[str, int]) -> dict[str, Any]:
    """Each forced arm reports the rank it was forced at."""

    deviations = {arm: float(abs(int(arms[arm]["rank"]) - int(rank))) for arm, rank in forced.items()
                  if (arms.get(arm) or {}).get("status") == "completed"}
    return _named_check(deviations, 0.0)


def check_same_fit(models: Mapping[str, Any], order: Sequence[str]) -> dict[str, Any]:
    """Per fit rank, every arm's model against the first one fitted at that rank (rotations, Y loadings)."""

    first: dict[int, Any] = {}
    deviations: dict[str, float] = {}
    for arm in order:
        if arm not in models:
            continue
        model = models[arm]
        rank = int(np.asarray(model.x_rotations_).shape[1])
        if rank not in first:
            first[rank] = model
            continue
        deviations[arm] = psc._model_deviation(first[rank], model)
    return _named_check(deviations, CONSISTENCY_TOLERANCE)


def check_measured_dimension(arms: Mapping[str, Mapping[str, Any]], specs: Mapping[str, Mapping[str, Any]],
                             n_stages: int) -> dict[str, Any]:
    deviations = {
        arm: float(abs(int(entry["measured_dimension"]) - min(int(entry["rank"]), int(n_stages) - 1)))
        for arm, entry in arms.items()
        if entry.get("status") == "completed" and specs[arm]["coordinates"] == PREDICTIVE
    }
    return _named_check(deviations, 0.0)


def check_same_statistics(arms: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """Each within-stage arm and the control arm observe exactly their counterpart's statistics."""

    deviations = {}
    for arm, entry in arms.items():
        source = counterpart(arm)
        if source is None or entry.get("status") != "completed":
            continue
        if (arms.get(source) or {}).get("status") != "completed":
            continue
        deviations[arm] = _statistics_deviation(entry["pair_statistics"], arms[source]["pair_statistics"])
    return _named_check(deviations, 0.0)


def check_cell_table(expected: Sequence[Sequence[int]], realized: Mapping[str, Sequence[Sequence[int]]]) -> dict:
    """The realized group × stage tables (the analysed labels, each arm's record) against the point's."""

    target = np.asarray(expected, dtype=float)
    deviations = {
        source: float(np.max(np.abs(np.asarray(table, dtype=float) - target))) if np.shape(table) == target.shape
        else math.inf
        for source, table in realized.items()
    }
    return _named_check(deviations, 0.0)


def recorded_table(result: Any, contract: CohortContract) -> list[list[int]]:
    sizes = result.record["design"]["cell_sizes"]
    return [[int(sizes[group][stage]) for stage in contract.stage_order] for group in contract.groups]


# ── Fit readings (design D12) ────────────────────────────────────────────────


def fit_readings(
    *, inputs: Any, results: Mapping[str, Any], settings: Mapping[str, Any], fits: Sequence[str], folds_seed: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Per fit rank the in-sample stage R² and out-of-fold explained share; per system and fit the dispersion.

    Every helper receives each row's stage position in the declared order
    (``inputs.working``'s stage codes), so labels that sort out of order are
    never handed to a helper that sorts them.
    """

    from motco import analysis

    order = [str(index) for index in range(len(inputs.stage_levels))]
    strata = inputs.working[analysis._STAGE].astype(str).to_numpy()
    assignment = pse.stratified_folds(strata, order, int(folds_seed), folds=int(settings["folds"]))
    readings: dict[str, Any] = {}
    dispersion: dict[str, Any] = {coordinates: {} for coordinates in COORDINATE_SYSTEMS}
    for fit in fits:
        scores_arm = arm_name(SCORES, fit)
        entry: dict[str, Any] = {}
        result = results.get(scores_arm)
        if result is None:
            readings[fit] = {"error": f"the arm {scores_arm} did not complete"}
        else:
            model = result.latent_artifacts.model
            features = result.latent_artifacts.features.to_numpy(dtype=float)
            rank = int(np.asarray(model.x_rotations_).shape[1])
            entry = {"rank": rank, "in_sample_r2": psc.stage_r2(model, features, strata.tolist())}
            try:
                built = pse.out_of_fold_fitted(features, strata, assignment, rank=rank,
                                               scale=bool(getattr(model, "scale", True)))
                entry["out_of_fold_share"] = pse.explained_share(strata, built["predictions"])
            except Exception as exc:  # a reading, never decisive
                entry["out_of_fold_share"] = None
                entry["out_of_fold_error"] = f"{type(exc).__name__}: {exc}"
            readings[fit] = entry
        for coordinates in COORDINATE_SYSTEMS:
            source = results.get(arm_name(coordinates, fit))
            if source is not None:
                dispersion[coordinates][fit] = pse.stage_dispersion(source.latent_scores, inputs.model_reduced, strata,
                                                                    order)
    return readings, dispersion


# ── One unit ─────────────────────────────────────────────────────────────────


def _failed_unit(base: Mapping[str, Any], exc: Exception, start: float, timing: dict[str, float]) -> dict[str, Any]:
    return lsc._jsonable(psc._failed_unit(base, exc, start, timing))


def evaluate_unit(unit: Unit, settings: Mapping[str, Any], cohort: Cohort) -> dict[str, Any]:
    """Build the unit, fit production once, evaluate every arm; one record of seeds and statistics."""

    from motco import analysis

    start = perf_counter()
    contract = cohort.contract
    point = point_spec(settings, unit.point)
    subsampled = bool(point["subsampled"])
    specs = settings["arm_specs"]
    seeds = unit_seeds(settings, unit.point, unit.replicate)
    seed = int(seeds["rrpp"])
    n_stages = len(contract.stage_order)
    base: dict[str, Any] = {
        "point": unit.point,
        "replicate": int(unit.replicate),
        "participants": int(point["participants"]),
        "seeds": seeds,
        "rrpp_seed": seed,
        "settings_hash": lsc.settings_hash(settings),
    }
    timing: dict[str, float] = {}
    arms: dict[str, Any] = {}
    results: dict[str, Any] = {}
    try:
        mark = perf_counter()
        built = build_unit(cohort, settings, unit.point, unit.replicate)
        base.update({"subsample_redraws": built["subsample_redraws"], "label_redraws": built["label_redraws"],
                     "cell_table": built["table"]})
        dataset = make_dataset(built["blocks"], built["metadata"], contract)
        timing["build"] = perf_counter() - mark
        production_params = arm_params(REFERENCE_ARM, settings, seed, subsampled=subsampled,
                                       production_rank=int(settings["cohort_rank"]))
        mark = perf_counter()
        production = _analyse(dataset, production_params)
        timing["production"] = perf_counter() - mark
        latent = production.record["latent"]
        rank = int(latent["selected_lv"])
    except Exception as exc:  # recorded, never aborts the shard
        return _failed_unit(base, exc, start, timing)

    production_fit = "cross_validated" if subsampled else "forced_at_cohort_rank"
    base["production"] = {
        "rank": rank,
        "component_selection": latent.get("component_selection"),
        "cv_mean_auroc": latent.get("cv_mean_auroc"),
        "fit": production_fit,
    }
    arms[REFERENCE_ARM] = arm_record(REFERENCE_ARM, production, timing["production"], specs[REFERENCE_ARM],
                                     production_fit)
    results[REFERENCE_ARM] = production
    forced: dict[str, int] = {} if subsampled else {REFERENCE_ARM: rank}
    order = point_arms(settings, unit.point)
    for arm in order:
        if arm == REFERENCE_ARM:
            continue
        spec = specs[arm]
        mark = perf_counter()
        target = rank if spec["fit"] == FIT_CV else int(spec["fixed_rank"])
        forced[arm] = target
        fit = ("forced_at_cv_rank" if subsampled else "forced_at_cohort_rank") if spec["fit"] == FIT_CV \
            else "forced_at_fixed_rank"
        try:
            result = _analyse(dataset, arm_params(arm, settings, seed, subsampled=subsampled, production_rank=rank))
            arms[arm] = arm_record(arm, result, perf_counter() - mark, spec, fit)
            results[arm] = result
        except Exception as exc:  # recorded on the arm, the unit goes on
            arms[arm] = failed_arm(arm, exc, perf_counter() - mark, spec)

    realized = {"labels": built["table"]}
    realized.update({arm: recorded_table(result, contract) for arm, result in results.items()})
    checks = {
        "rank": check_rank(arms, forced),
        "same_fit": check_same_fit({arm: result.latent_artifacts.model for arm, result in results.items()}, order),
        "measured_dimension": check_measured_dimension(arms, specs, n_stages),
        "same_statistics": check_same_statistics(arms),
        "cell_table": check_cell_table(point["table"], realized),
    }

    mark = perf_counter()
    fits = [FIT_CV, *(fit_label(r) for r in settings["fixed_ranks"] if subsampled or int(r) != rank)]
    try:
        inputs = analysis._trajectory_inputs(dataset, production_params)
        readings, dispersion = fit_readings(inputs=inputs, results=results, settings=settings, fits=fits,
                                            folds_seed=int(seeds["folds"]))
    except Exception as exc:  # readings, never decisive
        readings, dispersion = {"error": f"{type(exc).__name__}: {exc}"}, {}
    timing["fit_readings"] = perf_counter() - mark
    failing = [name for name in CONSISTENCY_CHECKS if not checks[name]["passed"]]
    timing["unit"] = perf_counter() - start
    return lsc._jsonable({
        **base,
        "status": "inconsistent" if failing else "completed",
        "inconsistent_checks": failing,
        "consistency": checks,
        "fit_readings": readings,
        "dispersion": dispersion,
        "folds": {"seed": int(seeds["folds"]), "count": int(settings["folds"])},
        "arms": {arm: arms[arm] for arm in order if arm in arms},
        "runtime_seconds": timing,
    })


# ── Pre-flight (design D5) ───────────────────────────────────────────────────


def preflight_unit(replicate: int, settings: Mapping[str, Any], cohort: Cohort) -> dict[str, Any]:
    """A cohort unit's cross-validated ``pls_cv`` against the contract's fidelity and the forced call."""

    contract = cohort.contract
    built = build_unit(cohort, settings, COHORT_POINT, replicate)
    dataset = make_dataset(built["blocks"], built["metadata"], contract)
    seed = int(built["seeds"]["rrpp"])
    mark = perf_counter()
    cross_validated = _analyse(dataset, analysis_params(settings, seed))
    cv_seconds = perf_counter() - mark
    rank = int(settings["cohort_rank"])
    mark = perf_counter()
    forced = _analyse(dataset, analysis_params(settings, seed, forced=rank))
    forced_seconds = perf_counter() - mark
    latent = cross_validated.record["latent"]
    fidelity = settings["contract"]["fidelity"]
    auroc = latent.get("cv_mean_auroc")
    rounded = None if auroc is None else round(float(auroc), int(fidelity["decimals"]))
    fidelity_passed = int(latent["selected_lv"]) == int(fidelity["selected_lv"]) and rounded is not None and \
        math.isclose(rounded, float(fidelity["cv_mean_auroc"]), rel_tol=0.0, abs_tol=1e-12)
    differing = {
        key: sorted(statistic for statistic in set(cross_validated.record[key]) | set(forced.record[key])
                    if cross_validated.record[key].get(statistic) != forced.record[key].get(statistic))
        for key in ("statistics", "p_values")
    }
    forced_rank = int(forced.record["latent"]["selected_lv"])
    passed = bool(fidelity_passed) and forced_rank == rank and not any(differing.values())
    return lsc._jsonable({
        "unit": unit_key(COHORT_POINT, replicate),
        "rrpp_seed": seed,
        "label_redraws": built["label_redraws"],
        "fidelity": {"rank": int(latent["selected_lv"]), "cv_mean_auroc_rounded": rounded,
                     "expected": {"selected_lv": fidelity["selected_lv"], "cv_mean_auroc": fidelity["cv_mean_auroc"]},
                     "passed": bool(fidelity_passed)},
        "forced_rank": forced_rank,
        "statistics_differing": differing["statistics"],
        "p_values_differing": differing["p_values"],
        "cross_validated": {"statistics": cross_validated.record["statistics"],
                            "p_values": cross_validated.record["p_values"]},
        "forced": {"statistics": forced.record["statistics"], "p_values": forced.record["p_values"]},
        "seconds": {"cross_validated": cv_seconds, "forced": forced_seconds},
        "passed": passed,
    })


def preflight(out_dir: Path, settings: Mapping[str, Any], cohort: Cohort) -> dict[str, Any]:
    """Write the settings and ``preflight.json``; ``run`` waits on its ``passed``."""

    out_dir = Path(out_dir)
    if (out_dir / "settings.json").exists():
        verify(lsc.load_settings(out_dir))
    verify(settings, cohort)
    lsc.ensure_settings(out_dir, settings)
    units = []
    for replicate in settings["preflight_replicates"]:
        try:
            units.append(preflight_unit(int(replicate), settings, cohort))
        except Exception as exc:  # a pre-flight that raises fails
            units.append({"unit": unit_key(COHORT_POINT, int(replicate)), "passed": False,
                          "error": f"{type(exc).__name__}: {exc}"})
    result = {
        "settings_hash": lsc.settings_hash(settings),
        "rule": "the cross-validated pls_cv on the cohort selects the contract's rank and AUROC (rounded), and its "
                "statistics and p-values equal the call forced at that rank exactly",
        "expected": dict(settings["contract"]["fidelity"]),
        "units": units,
        "passed": bool(units) and all(unit["passed"] for unit in units),
    }
    lsc._write_json(out_dir / PREFLIGHT_NAME, result)
    return result


def load_preflight(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any] | None:
    path = Path(out_dir) / PREFLIGHT_NAME
    if not path.exists():
        return None
    record = json.loads(path.read_text(encoding="utf-8"))
    return record if record.get("settings_hash") == lsc.settings_hash(settings) else None


# ── Shards, resume and the run ───────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> list[dict[str, Any]]:
    """Every parsable unit line; a line torn by an interrupted write is dropped, so its unit is re-run."""

    return psc.read_records(paths)


_WORKER_COHORT: Cohort | None = None


def _init_worker(directory: str, contract: CohortContract, files: Mapping[str, str]) -> None:
    global _WORKER_COHORT
    cohort = load_cohort(Path(directory), contract)
    if cohort.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_COHORT = cohort


def _evaluate_pending(unit: Unit, settings: Mapping[str, Any]) -> dict[str, Any]:
    assert _WORKER_COHORT is not None
    return evaluate_unit(unit, settings, _WORKER_COHORT)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    cohort: Cohort,
    index: int,
    count: int,
    *,
    limit: int | None = None,
    workers: int = 1,
) -> int:
    """Evaluate this shard's missing units into its file; returns how many were evaluated."""

    out_dir = Path(out_dir)
    if (out_dir / "settings.json").exists():
        verify(lsc.load_settings(out_dir))
    verify(settings, cohort)
    lsc.ensure_settings(out_dir, settings)
    checked = load_preflight(out_dir, settings)
    if checked is None:
        raise SystemExit(f"{out_dir / PREFLIGHT_NAME} is missing or was written under other settings; run "
                         "`preflight` first.")
    if not checked.get("passed"):
        raise SystemExit(f"The pre-flight failed ({out_dir / PREFLIGHT_NAME}); the run stops.")
    path = lsc.shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    done = {record_key(record): record for record in (read_records([path]) if path.exists() else [])}
    # Rewritten with complete units only, so a torn line is re-run, not duplicated.
    with path.open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(done[key]) for key in sorted(done)))
    pending = [
        unit for position, unit in enumerate(enumerate_units(settings))
        if position % count == index and unit_key(*unit) not in done
    ]
    if limit is not None:
        pending = pending[:limit]

    def write(position: int, unit: Unit, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(lsc._line(record))
            handle.flush()
        print(f"[{index}/{count}] {position}/{len(pending)} {unit_key(*unit)} {record['status']}", file=sys.stderr,
              flush=True)

    if workers <= 1:
        for position, unit in enumerate(pending, start=1):
            write(position, unit, evaluate_unit(unit, settings, cohort))
        return len(pending)
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker,
                             initargs=(str(cohort.directory), cohort.contract, dict(cohort.files))) as pool:
        futures = {pool.submit(_evaluate_pending, unit, settings): unit for unit in pending}
        for position, future in enumerate(as_completed(futures), start=1):
            write(position, futures[future], future.result())
    return len(pending)


# ── Merge: rate, paired and point tables, accounting ─────────────────────────


RATE_COLUMNS: tuple[str, ...] = (
    "point", "participants", "arm", "family", "fit", "coordinates", "draws", "residual_correction", "statistic",
    "rate", "se", "count", "rejections", "failures", "inconsistent",
)
PAIRED_COLUMNS: tuple[str, ...] = (
    "point", "statistic", "arm", "reference_arm", "comparison", "n_paired", "rate_reference", "rate_arm",
    "difference", "b", "c", "p_loss", "p_gain", "p_two_sided",
)
POINT_COLUMNS: tuple[str, ...] = ("point", "quantity", "arm", "fit", "coordinates", "statistic", "key", "value")


def _rejects(p_value: Any, settings: Mapping[str, Any]) -> bool | None:
    return plsc._rejects(p_value, float(settings["alpha"]), str(settings["rejection"]))


def _arm(record: Mapping[str, Any], arm: str) -> Mapping[str, Any] | None:
    """The arm's completed entry on a completed (consistent) unit, else ``None``."""

    if record.get("status") != "completed":
        return None
    entry = (record.get("arms") or {}).get(arm)
    return entry if entry is not None and entry.get("status") == "completed" else None


def _p(record: Mapping[str, Any], arm: str, statistic: str) -> Any:
    entry = _arm(record, arm)
    return None if entry is None else (entry.get("p_values") or {}).get(statistic)


def _by_point(records: Iterable[Mapping[str, Any]]) -> dict[str, list[Mapping[str, Any]]]:
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for record in records:
        grouped.setdefault(str(record["point"]), []).append(record)
    return grouped


def rate_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Rejections (p ≤ α) per point × arm × statistic; failed and inconsistent units excluded."""

    grouped = _by_point(records)
    specs = settings["arm_specs"]
    rows = []
    for point in point_names(settings):
        group = grouped.get(point, [])
        inconsistent = sum(1 for record in group if record.get("status") == "inconsistent")
        participants = point_spec(settings, point)["participants"]
        for arm in point_arms(settings, point):
            failures = sum(1 for record in group if record.get("status") == "failed"
                           or (record.get("status") == "completed" and _arm(record, arm) is None))
            spec = specs[arm]
            for statistic in STATISTICS:
                decisions = [_rejects(_p(record, arm, statistic), settings) for record in group]
                available = [bool(value) for value in decisions if value is not None]
                count, hits = len(available), sum(available)
                rate = hits / count if count else None
                rows.append({
                    "point": point, "participants": participants, "arm": arm, "family": spec["family"],
                    "fit": spec["fit"], "coordinates": spec["coordinates"], "draws": spec["draws"],
                    "residual_correction": spec["residual_correction"] or "none", "statistic": statistic,
                    "rate": rate, "se": None if rate is None else math.sqrt(rate * (1.0 - rate) / count),
                    "count": count, "rejections": hits, "failures": failures, "inconsistent": inconsistent,
                })
    return rows


def paired_table(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Each within-stage and the control arm against its counterpart, each predictive arm against the scores."""

    grouped = _by_point(records)
    rows = []
    for point in point_names(settings):
        present = set(point_arms(settings, point))
        group = sorted(grouped.get(point, []), key=lambda record: int(record["replicate"]))
        for statistic in STATISTICS:
            for arm, reference in paired_arms(settings):
                if arm not in present or reference not in present:
                    continue
                pairs = [(_rejects(_p(record, reference, statistic), settings),
                          _rejects(_p(record, arm, statistic), settings)) for record in group]
                both = [(ref, cand) for ref, cand in pairs if ref is not None and cand is not None]
                n = len(both)
                b = sum(1 for ref, cand in both if ref and not cand)
                c = sum(1 for ref, cand in both if cand and not ref)
                rows.append({
                    "point": point, "statistic": statistic, "arm": arm, "reference_arm": reference,
                    "comparison": "same_measurement" if counterpart(arm) == reference else "coordinates",
                    "n_paired": n,
                    "rate_reference": sum(ref for ref, _ in both) / n if n else None,
                    "rate_arm": sum(cand for _, cand in both) / n if n else None,
                    "difference": (c - b) / n if n else None,
                    "b": b, "c": c, **plsc.paired_p_values(b, c),
                })
    return rows


def _median(values: Iterable[Any]) -> float | None:
    return plsc._quantile(values, 0.5)


def point_rows(records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], n_stages: int,
               n_features: int) -> list[dict[str, Any]]:
    """The point table in long form: one quantity per row; blank keys are not applicable."""

    grouped = _by_point(records)
    rows: list[dict[str, Any]] = []

    def add(point: str, quantity: str, value: Any, **keys: Any) -> None:
        rows.append({"point": point, "quantity": quantity, **keys, "value": value})

    floor_rank = int(n_stages) - 1
    for point in point_names(settings):
        spec = point_spec(settings, point)
        units = [record for record in grouped.get(point, []) if record.get("status") == "completed"]
        everything = grouped.get(point, [])
        add(point, "participants", spec["participants"])
        add(point, "features_per_participant", n_features / int(spec["participants"]))
        add(point, "units_completed", len(units))
        ranks = [int(record["production"]["rank"]) for record in units]
        for rank in sorted(set(ranks)):
            add(point, "cv_rank_units", ranks.count(rank), key=str(rank))
        add(point, "cv_rank_share_k_minus_1", ranks.count(floor_rank) / len(ranks) if ranks else None)
        add(point, "cv_rank_median", _median(ranks))
        add(point, "cv_auroc_median", _median(record["production"].get("cv_mean_auroc") for record in units))
        fits = [FIT_CV, *(fit_label(rank) for rank in settings["fixed_ranks"])]
        for fit in fits:
            readings = [(record.get("fit_readings") or {}).get(fit) or {} for record in units]
            if not any(readings):
                continue
            add(point, "in_sample_r2_median", _median(entry.get("in_sample_r2") for entry in readings), fit=fit)
            add(point, "out_of_fold_share_median", _median(entry.get("out_of_fold_share") for entry in readings),
                fit=fit)
            for coordinates in COORDINATE_SYSTEMS:
                entries = [((record.get("dispersion") or {}).get(coordinates) or {}).get(fit) for record in units]
                entries = [entry for entry in entries if entry]
                if not entries:
                    continue
                for key in (*(f"stage_{index}" for index in range(int(n_stages))), "max_min_ratio"):
                    add(point, "dispersion_median", _median(entry.get(key) for entry in entries), fit=fit,
                        coordinates=coordinates, key=key)
        for arm in point_arms(settings, point):
            add(point, "angle_null_q95_median", _median((_arm(record, arm) or {}).get("angle_null_q95")
                                                        for record in units), arm=arm)
        classes = (("k_minus_1", [r for r in units if int(r["production"]["rank"]) == floor_rank]),
                   ("above_k_minus_1", [r for r in units if int(r["production"]["rank"]) > floor_rank]))
        for arm in OPERATING_ARMS:
            for label, members in classes:
                for statistic in STATISTICS:
                    decisions = [_rejects(_p(record, arm, statistic), settings) for record in members]
                    available = [bool(value) for value in decisions if value is not None]
                    add(point, "rank_class_count", len(available), arm=arm, statistic=statistic, key=label)
                    add(point, "rank_class_rejections", sum(available), arm=arm, statistic=statistic, key=label)
        add(point, "subsample_redraws", sum(int(record.get("subsample_redraws") or 0) for record in everything))
        add(point, "label_redraws", sum(int(record.get("label_redraws") or 0) for record in everything))
    return rows


def accounting(
    records: Sequence[Mapping[str, Any]], settings: Mapping[str, Any], *, expected: Mapping[str, int], duplicates: int
) -> dict[str, Any]:
    """Units, failures and inconsistent units by point, arm failures, redraws, core-hours and per-arm timing."""

    by_point: dict[str, dict[str, int]] = {}
    inconsistent_units: list[dict[str, Any]] = []
    failed_units: list[dict[str, Any]] = []
    arm_failures: dict[str, dict[str, int]] = {}
    arm_seconds: dict[str, list[float]] = {}
    unit_seconds: dict[str, list[float]] = {}
    for record in records:
        point = str(record["point"])
        counts = by_point.setdefault(point, {"complete": 0, "failed": 0, "inconsistent": 0})
        counts["complete"] += 1
        timing = record.get("runtime_seconds") or {}
        if timing.get("unit") is not None:
            unit_seconds.setdefault(point, []).append(float(timing["unit"]))
        status = record.get("status")
        if status == "failed":
            counts["failed"] += 1
            failed_units.append({"unit": record_key(record), "error_type": record.get("error_type"),
                                 "error_message": record.get("error_message")})
            continue
        if status == "inconsistent":
            counts["inconsistent"] += 1
            inconsistent_units.append({
                "unit": record_key(record),
                "checks": {name: {"deviation": record["consistency"][name].get("deviation"),
                                  "failing_arms": record["consistency"][name].get("failing_arms")}
                           for name in record["inconsistent_checks"]},
            })
        for arm, entry in (record.get("arms") or {}).items():
            if entry.get("status") != "completed":
                failures = arm_failures.setdefault(str(arm), {})
                failures[point] = failures.get(point, 0) + 1
            elif entry.get("seconds") is not None:
                arm_seconds.setdefault(str(arm), []).append(float(entry["seconds"]))
    points = point_names(settings)
    total_expected = sum(int(expected.get(point, 0)) for point in points)
    complete = sum(counts["complete"] for counts in by_point.values())
    all_seconds = [value for values in unit_seconds.values() for value in values]
    return {
        "settings_hash": lsc.settings_hash(settings),
        "expected_units": total_expected,
        "complete_units": complete,
        "missing_units": total_expected - complete,
        "failed_units": sum(counts["failed"] for counts in by_point.values()),
        "inconsistent_units": sum(counts["inconsistent"] for counts in by_point.values()),
        "by_point": {
            point: {
                "expected": int(expected.get(point, 0)),
                **by_point.get(point, {"complete": 0, "failed": 0, "inconsistent": 0}),
                "usable": by_point.get(point, {}).get("complete", 0) - by_point.get(point, {}).get("failed", 0)
                - by_point.get(point, {}).get("inconsistent", 0),
                "unit_seconds_median": _median(unit_seconds.get(point, [])),
            }
            for point in points
        },
        "failed_unit_errors": failed_units,
        "inconsistent_unit_checks": inconsistent_units,
        "arm_failures_by_point": arm_failures,
        "duplicate_records": int(duplicates),
        "recorded_core_hours": sum(all_seconds) / 3600.0,
        "arm_seconds_median": {arm: _median(arm_seconds.get(arm, [])) for arm in settings["arms"]},
        "arm_seconds_total_hours": {arm: sum(arm_seconds.get(arm, [])) / 3600.0 for arm in settings["arms"]},
    }


def merge(out_dir: Path, export_dir: Path) -> dict[str, Any]:
    """Merge every shard and write the rate, paired and point tables and the accounting."""

    out_dir = Path(out_dir)
    settings = lsc.load_settings(out_dir)
    verify(settings)
    verify_export_files(export_dir, settings)
    raw = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    stale = sorted({record_key(r) for r in raw if r.get("settings_hash") != lsc.settings_hash(settings)})
    if stale:
        raise SystemExit(f"{len(stale)} unit(s) were written under other settings, e.g. {stale[:3]}; refusing.")
    seen: dict[str, int] = {}
    by_key: dict[str, dict[str, Any]] = {}
    for record in raw:
        key = record_key(record)
        seen[key] = seen.get(key, 0) + 1
        by_key[key] = record
    units = enumerate_units(settings)
    records = [by_key[unit_key(*unit)] for unit in units if unit_key(*unit) in by_key]
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        handle.write("".join(lsc._line(record) for record in records))
    contract = settings["contract"]
    n_features = sum(int(value) for value in settings["export"]["features"].values())
    lsc._write_csv(out_dir / "rate_table.csv", rate_table(records, settings), RATE_COLUMNS)
    lsc._write_csv(out_dir / "paired_table.csv", paired_table(records, settings), PAIRED_COLUMNS)
    lsc._write_csv(out_dir / "point_table.csv",
                   point_rows(records, settings, len(contract["stage_order"]), n_features), POINT_COLUMNS)
    expected: dict[str, int] = {}
    for unit in units:
        expected[unit.point] = expected.get(unit.point, 0) + 1
    summary = accounting(records, settings, expected=expected, duplicates=sum(n - 1 for n in seen.values()))
    lsc._write_json(out_dir / "accounting.json", summary)
    return {key: summary[key] for key in ("expected_units", "complete_units", "missing_units", "failed_units",
                                          "inconsistent_units", "duplicate_records", "recorded_core_hours")}


# ── Decision (pure functions; designs D8–D12) ────────────────────────────────


VERDICTS: tuple[str, ...] = ("gate_failed", "production_excess", "excess_persists", "excess_scoped", "excess_absent")
CONSEQUENCES: dict[str, str] = {
    "gate_failed": (
        "No reading is interpreted. The cause is fixed and the affected units are rerun under the same rule."
    ),
    "production_excess": (
        "Production's unrestricted null with df exceeds size on ADNI's noise at the named readings. No reading of "
        "predictive_orthonormal is interpreted at a point where production fails. A diagnosis of those points comes "
        "before any other reading on ADNI. If the cohort is among them, the roadmap states that a rejection at ADNI's "
        "design lacks a size-controlled null on its own noise; the case study's primary non-rejections stand, since "
        "an excess cannot produce them."
    ),
    "excess_persists": (
        "The excess appears at ADNI's own design, at 0.2 features per participant, so item 19's prediction for ADNI "
        "is refuted. The veto rests on two cohorts and the second reopening condition of open item 1 is closed for "
        "ADNI. predictive_orthonormal stays an opt-in diagnostic."
    ),
    "excess_scoped": (
        "predictive_orthonormal holds size on ADNI's blocks from the named point up to the cohort and exceeds below "
        "it. The size part of open item 1's second condition is met for ADNI inside that scope only. It stays an "
        "opt-in diagnostic; the gain and specificity readings at ADNI's design are the next change and are confined "
        "to the scope."
    ),
    "excess_absent": (
        "predictive_orthonormal holds size on ADNI's blocks at every design point under RRPP's unrestricted draws: "
        "SEA-AD's excess does not transfer to ADNI. The size part of open item 1's second condition is met for ADNI. "
        "It stays an opt-in diagnostic until the gain and specificity readings at ADNI's design, the next change. "
        "The fixed-rank grade says whether the hold is owed to ADNI's noise or to the ranks its CV selects."
    ),
}
STANDING = (
    "In every verdict the SEA-AD anchor's reading stands for SEA-AD's design, and keep_production stands."
)
GRADES: tuple[str, ...] = ("reference_exceeds", "absent", "reproduced", "reproduced_unrestored")
GRADE_STATEMENTS: dict[str, str] = {
    "reference_exceeds": (
        "An over-fitted rank mis-sizes the scores themselves on ADNI's noise. The predictive arms are reported "
        "without a grade, and the roadmap records that production's size depends on a rank its CV does not select "
        "here."
    ),
    "absent": (
        "No excess at either fixed rank. An over-fitted stage fit is not sufficient on ADNI's noise: SEA-AD's excess "
        "needs a property of its blocks this ladder lacks. The lead is a feature-subsampling ladder on SEA-AD's own "
        "export, then ROSMAP."
    ),
    "reproduced": (
        "ADNI's noise gives the predictive subspace the excess once the fit is over-fitted, and the within-stage "
        "draws remove it. Item 19's mechanism is not specific to SEA-AD's blocks. The report tabulates the excess "
        "against the in-sample stage fit."
    ),
    "reproduced_unrestored": (
        "The excess appears at an over-fitted rank and the within-stage draws do not remove it at the named "
        "readings. Unequal stage dispersion is not the whole mechanism on ADNI's noise; the dispersion table and the "
        "out-of-fold shares are the leads."
    ),
}
SCOPING_RULE: dict[str, str] = {
    "S1": "A size reading speaks for its export and its design points. The points are declared before the run. "
          "Nothing is interpolated between them or extended beyond them.",
    "S2": "An option's scope on a cohort is its size scope (design D8). An option with no scope there is not "
          "admissible on that cohort.",
    "S3": "Cohorts do not vouch for one another. A hold on one export does not lift a failure on another and does "
          "not extend to a cohort without its own reading. The SEA-AD anchor's reading stands for SEA-AD's design.",
    "S4": "Size is one of three readings. Reopening an option at a design also needs a gain over production at that "
          "design in simulation and specificity there, each under a decision declared before its reading. Reopening "
          "means that a separate change may propose its paper-grade run; only that run can change a default.",
    "S5": "Covariates describe; they do not extend. Features per participant, the selected ranks and the in-sample "
          "stage fit are reported beside a scope. They license no reading at a design that was not read.",
}
SIZE_SCOPE_NOTE = (
    "A scope is the size part of open item 1's reopening condition only; the gain over production in simulation at "
    "ADNI's design and specificity there are still missing, each under its own predeclared decision."
)
DESCRIPTIVE_NOTE = (
    "Features per participant, the selected ranks and the in-sample stage fit describe the design points; they "
    "license no reading at another design (S5). Every table under 'descriptive' is a reading beside the verdict and "
    "the grade and never changes them."
)

#: A reading: (point, arm, statistic).
Reading = tuple[str, str, str]


def index_rates(rows: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str, str], Mapping[str, Any]]:
    return {(str(row["point"]), str(row["arm"]), str(row["statistic"])): row for row in rows}


def index_paired(rows: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str, str, str], Mapping[str, Any]]:
    return {(str(row["point"]), str(row["statistic"]), str(row["arm"]), str(row["reference_arm"])): row
            for row in rows}


def observation(rates: Mapping[tuple[str, str, str], Mapping[str, Any]], reading: Reading) -> dict[str, Any]:
    point, arm, statistic = reading
    entry = {"point": point, "arm": arm, "statistic": statistic}
    row = rates.get(reading)
    if row is None or row.get("count") in (None, "") or int(float(row["count"])) <= 0:
        return {**entry, "status": "unavailable"}
    return {**entry, "n_available": int(float(row["count"])), "n_rejected": int(float(row["rejections"]))}


def size_sets(points: Sequence[str]) -> list[list[str]]:
    """The upper sets of the ladder, each reaching the cohort: all points, then from the second on, …"""

    return [list(points[start:]) for start in range(len(points))]


def judge_family(rates: Mapping[tuple[str, str, str], Mapping[str, Any]], readings: Sequence[Reading]) -> dict:
    judged = grade.judge_rates([observation(rates, reading) for reading in readings])
    judged["failing_readings"] = [row for row in judged["observations"] if row.get("status") != "holds"]
    return judged


def size_scope(rates: Mapping[tuple[str, str, str], Mapping[str, Any]], arm: str, points: Sequence[str]) -> dict:
    """The arm's families on the ladder's upper sets and its scope: the first set that holds, or none."""

    families = []
    scope: list[str] | None = None
    for members in size_sets(points):
        judged = judge_family(rates, [(point, arm, statistic) for point in members for statistic in STATISTICS])
        families.append({"points": members, **judged})
        if scope is None and judged["holds"]:
            scope = members
    return {
        "arm": arm,
        "scope_start": None if scope is None else scope[0],
        "scope": scope,
        "holds_everywhere": scope is not None and scope[0] == points[0],
        "families": families,
    }


def gate(preflight_record: Mapping[str, Any] | None, summary: Mapping[str, Any]) -> dict[str, Any]:
    failures = []
    if preflight_record is None:
        failures.append("the pre-flight record is missing or was written under other settings")
    elif not preflight_record.get("passed"):
        failures.append("the pre-flight did not pass")
    inconsistent = int(summary.get("inconsistent_units", 0) or 0)
    if inconsistent:
        failures.append(f"{inconsistent} inconsistent unit(s)")
    return {"passed": not failures, "failures": failures}


def verdict_of(gated: Mapping[str, Any], production: Mapping[str, Any], candidate: Mapping[str, Any],
               points: Sequence[str]) -> str:
    if not gated["passed"]:
        return "gate_failed"
    if production["scope_start"] != points[0]:
        return "production_excess"
    if candidate["scope_start"] is None:
        return "excess_persists"
    if candidate["scope_start"] != points[0]:
        return "excess_scoped"
    return "excess_absent"


def fixed_rank_readings(settings: Mapping[str, Any], coordinates: str, draws: str) -> list[Reading]:
    """The family of one fixed-rank combination: every fixed rank at every subsampled point, and at the cohort the
    ranks other than the cohort's (its readings are the operating arms')."""

    readings = []
    for rank in settings["fixed_ranks"]:
        arm = arm_name(coordinates, fit_label(int(rank)), draws)
        for point in settings["points"]:
            if not point["subsampled"] and int(rank) == int(settings["cohort_rank"]):
                continue
            readings.extend((str(point["name"]), arm, statistic) for statistic in STATISTICS)
    return readings


def _r2_lookup(point_table_rows: Iterable[Mapping[str, Any]]) -> dict[tuple[str, str], float | None]:
    return {(str(row["point"]), str(row["fit"])): plsc._number(row["value"]) for row in point_table_rows
            if str(row["quantity"]) == "in_sample_r2_median"}


def fixed_rank_grade(rates: Mapping[tuple[str, str, str], Mapping[str, Any]], settings: Mapping[str, Any],
                     point_table_rows: Sequence[Mapping[str, Any]] = ()) -> dict[str, Any]:
    """Design D10: the four combinations' 21-reading families and the first grade that applies."""

    specs = settings["arm_specs"]
    r2 = _r2_lookup(point_table_rows)
    families: dict[str, Any] = {}
    for coordinates in COORDINATE_SYSTEMS:
        for draws in ("unrestricted", "within_stage"):
            judged = judge_family(rates, fixed_rank_readings(settings, coordinates, draws))
            for row in judged["observations"]:
                fit = specs[row["arm"]]["fit"]
                row.update({"fit": fit, "in_sample_r2_median": r2.get((row["point"], fit))})
            for row in judged["failing_readings"]:
                fit = specs[row["arm"]]["fit"]
                row.update({"fit": fit, "in_sample_r2_median": r2.get((row["point"], fit))})
            families[f"{coordinates}__{draws}"] = judged
    scores, candidate = families[f"{SCORES}__unrestricted"], families[f"{PREDICTIVE}__unrestricted"]
    restricted = families[f"{PREDICTIVE}__within_stage"]
    if not scores["holds"]:
        name = "reference_exceeds"
    elif candidate["holds"]:
        name = "absent"
    elif restricted["holds"]:
        name = "reproduced"
    else:
        name = "reproduced_unrestored"
    named = (candidate["failing_readings"] if name in ("reproduced", "reproduced_unrestored")
             else scores["failing_readings"] if name == "reference_exceeds" else [])
    return {
        "grade": name,
        "statement": GRADE_STATEMENTS[name],
        "statements": dict(GRADE_STATEMENTS),
        "named_readings": named,
        "graded_arms": [] if name == "reference_exceeds" else [f"{PREDICTIVE}__unrestricted",
                                                               f"{PREDICTIVE}__within_stage"],
        "families": families,
        "rule": "one family per combination of coordinates and null: three statistics at every fixed rank on every "
                "subsampled point and at the ranks other than the cohort's on the cohort, judged by the "
                "count-plus-ceiling rule; it never changes the verdict",
    }


def indicative(settings: Mapping[str, Any], summary: Mapping[str, Any]) -> dict[str, Any]:
    reasons = []
    if int(settings["replicates"]) < REPLICATES:
        reasons.append(f"{settings['replicates']} replicates per point (paper grade: {REPLICATES})")
    if int(settings["permutations"]) < PERMUTATIONS:
        reasons.append(f"{settings['permutations']} permutations (paper grade: {PERMUTATIONS})")
    for point, counts in (summary.get("by_point") or {}).items():
        if int(counts.get("usable", 0)) < REPLICATES:
            reasons.append(f"{counts.get('usable', 0)} usable units at {point}")
    return {"indicative": bool(reasons), "reasons": reasons}


def decision(
    settings: Mapping[str, Any],
    rates_rows: Sequence[Mapping[str, Any]],
    paired_rows: Sequence[Mapping[str, Any]],
    point_table_rows: Sequence[Mapping[str, Any]],
    summary: Mapping[str, Any],
    preflight_record: Mapping[str, Any] | None,
    *,
    seaad_rates: Sequence[Mapping[str, Any]] = (),
    seaad_cells: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """The predeclared verdict (rule version 1), the fixed-rank grade and the descriptive readings beside them."""

    rates = index_rates(rates_rows)
    points = point_names(settings)
    gated = gate(preflight_record, summary)
    production = size_scope(rates, REFERENCE_ARM, points)
    candidate = size_scope(rates, CANDIDATE_ARM, points)
    within = {arm: size_scope(rates, arm, points) for arm in WITHIN_STAGE_OPERATING_ARMS}
    verdict = verdict_of(gated, production, candidate, points)
    named: list[dict[str, Any]] = []
    if verdict == "production_excess":
        named = production["families"][0]["failing_readings"]
    elif verdict in ("excess_persists", "excess_scoped"):
        named = candidate["families"][0]["failing_readings"]
    graded = fixed_rank_grade(rates, settings, point_table_rows)
    return {
        "rule_version": RULE_VERSION,
        "driver_version": DRIVER_VERSION,
        "verdict": verdict,
        "consequence": f"{CONSEQUENCES[verdict]} {STANDING}",
        "standing": STANDING,
        "named_readings": named,
        "scope_start": candidate["scope_start"] if verdict == "excess_scoped" else None,
        "consequences": dict(CONSEQUENCES),
        "gate": gated,
        "indicative": indicative(settings, summary),
        "operating": {
            "reference": production,
            "candidate": candidate,
            "candidate_interpreted_at": [] if verdict == "gate_failed" else (
                points if verdict != "production_excess" else
                [point for point in points if point in (production["scope"] or [])]),
        },
        "within_stage_operating": {
            arm: {**judged, "named_readings": [] if judged["holds_everywhere"] else judged["families"][0][
                "failing_readings"]}
            for arm, judged in within.items()
        },
        "fixed_rank": graded,
        "scoping_rule": dict(SCOPING_RULE),
        "size_scope_note": SIZE_SCOPE_NOTE,
        "descriptive_note": DESCRIPTIVE_NOTE,
        "alpha": float(settings["alpha"]),
        "rejection": str(settings["rejection"]),
        "points": [{key: point[key] for key in ("name", "participants", "table", "subsampled")}
                   for point in settings["points"]],
        "descriptive": descriptive_tables(settings, rates_rows, paired_rows, point_table_rows,
                                          seaad_rates=seaad_rates, seaad_cells=seaad_cells),
    }


def _point_value(rows: Iterable[Mapping[str, Any]], point: str, quantity: str, **keys: str) -> float | None:
    for row in rows:
        if str(row["point"]) == point and str(row["quantity"]) == quantity and all(
                str(row.get(key) or "") == value for key, value in keys.items()):
            return plsc._number(row["value"])
    return None


def descriptive_tables(
    settings: Mapping[str, Any],
    rates_rows: Sequence[Mapping[str, Any]],
    paired_rows: Sequence[Mapping[str, Any]],
    point_table_rows: Sequence[Mapping[str, Any]],
    *,
    seaad_rates: Sequence[Mapping[str, Any]] = (),
    seaad_cells: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Design D12: rates with their status, rank classes, the control arm, paired tests, points, cross-cohort."""

    rates = index_rates(rates_rows)
    points = point_names(settings)
    every_arm = {
        arm: judge_family(rates, [(point, arm, statistic) for point in points if arm in point_arms(settings, point)
                                  for statistic in STATISTICS])["observations"]
        for arm in settings["arms"]
    }
    rank_classes = []
    for point in points:
        for arm in OPERATING_ARMS:
            for label in ("k_minus_1", "above_k_minus_1"):
                for statistic in STATISTICS:
                    count = _point_value(point_table_rows, point, "rank_class_count", arm=arm, statistic=statistic,
                                         key=label)
                    hits = _point_value(point_table_rows, point, "rank_class_rejections", arm=arm,
                                        statistic=statistic, key=label)
                    rank_classes.append({"point": point, "arm": arm, "rank_class": label, "statistic": statistic,
                                         "count": count, "rejections": hits,
                                         "rate": hits / count if count and hits is not None else None})
    control = [
        {**{key: row.get(key) for key in PAIRED_COLUMNS}}
        for row in paired_rows if str(row["arm"]) == CONTROL_ARM
    ]
    cross = cross_cohort_table(settings, rates, point_table_rows, seaad_rates, seaad_cells)
    return {
        "rates_with_status": every_arm,
        "rank_classes": rank_classes,
        "control_arm": control,
        "paired": [{key: row.get(key) for key in PAIRED_COLUMNS} for row in paired_rows],
        "points": [{key: row.get(key) for key in POINT_COLUMNS} for row in point_table_rows],
        "angle_null_q95_median": [
            {"point": str(row["point"]), "arm": str(row["arm"]), "value": plsc._number(row["value"])}
            for row in point_table_rows if str(row["quantity"]) == "angle_null_q95_median"
        ],
        "cross_cohort": cross,
    }


def cross_cohort_table(
    settings: Mapping[str, Any],
    rates: Mapping[tuple[str, str, str], Mapping[str, Any]],
    point_table_rows: Sequence[Mapping[str, Any]],
    seaad_rates: Sequence[Mapping[str, Any]],
    seaad_cells: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """SEA-AD's committed anchor counts beside each ADNI point's, for the operating arms."""

    rows = []
    reference_cell = {
        str(row["cell"]): row for row in seaad_cells
        if str(row.get("design")) == SEAAD_DESIGN and str(row.get("arm")) == REFERENCE_ARM
    }
    for arm in OPERATING_ARMS:
        for cell in SEAAD_CELLS:
            counts = {
                str(row["statistic"]): {"count": int(float(row["count"])), "rejections": int(float(row["rejections"]))}
                for row in seaad_rates
                if str(row.get("design")) == SEAAD_DESIGN and str(row.get("cell")) == cell and str(row["arm"]) == arm
            }
            source = reference_cell.get(cell) or {}
            rows.append({
                "cohort": "seaad", "point": cell, "arm": arm, "features_per_participant": SEAAD_FEATURES_PER_DONOR,
                "fit_rank": plsc._number(source.get("rank_median")),
                "in_sample_r2": plsc._number(source.get("stage_r2_median")),
                **{statistic: counts.get(statistic) for statistic in STATISTICS},
            })
        for point in point_names(settings):
            observed = {statistic: observation(rates, (point, arm, statistic)) for statistic in STATISTICS}
            rows.append({
                "cohort": settings["contract"]["name"], "point": point, "arm": arm,
                "features_per_participant": _point_value(point_table_rows, point, "features_per_participant"),
                "fit_rank": _point_value(point_table_rows, point, "cv_rank_median"),
                "in_sample_r2": _point_value(point_table_rows, point, "in_sample_r2_median", fit=FIT_CV),
                **{statistic: None if entry.get("status") == "unavailable" else
                   {"count": entry["n_available"], "rejections": entry["n_rejected"]}
                   for statistic, entry in observed.items()},
            })
    return rows


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run `merge` first.")
    return json.loads(path.read_text(encoding="utf-8"))


_TABLE_FILES: dict[str, str] = {"rates": "rate_table.csv", "paired": "paired_table.csv", "points": "point_table.csv"}


def decide(out_dir: Path, export_dir: Path) -> dict[str, Any]:
    """Read the merged tables, the accounting and the pre-flight; write ``adni_anchor_decision.json``."""

    out_dir = Path(out_dir)
    settings = lsc.load_settings(out_dir)
    verify(settings)
    verify_export_files(export_dir, settings)
    summary = _read_json(out_dir / "accounting.json")
    if summary.get("settings_hash") != lsc.settings_hash(settings):
        raise SystemExit(f"{out_dir}/accounting.json was merged under other settings; run `merge` again.")
    tables = {}
    for name, filename in _TABLE_FILES.items():
        path = out_dir / filename
        if not path.exists():
            raise SystemExit(f"{path} does not exist; run `merge` first.")
        tables[name] = lsc._read_csv(path)
    checked = load_preflight(out_dir, settings)
    seaad_rates = lsc._read_csv(plsc._resolve(settings["parents"]["seaad_anchor_rates"]["path"]))
    seaad_cells = lsc._read_csv(plsc._resolve(settings["parents"]["seaad_anchor_cells"]["path"]))
    result = decision(settings, tables["rates"], tables["paired"], tables["points"], summary, checked,
                      seaad_rates=seaad_rates, seaad_cells=seaad_cells)
    result["complete"] = {key: summary[key] for key in ("expected_units", "complete_units", "missing_units",
                                                        "failed_units", "inconsistent_units")}
    result["preflight"] = None if checked is None else {
        "passed": checked.get("passed"),
        "units": [{key: unit.get(key) for key in ("unit", "passed", "fidelity", "forced_rank", "statistics_differing",
                                                  "p_values_differing", "error")} for unit in checked.get("units", [])],
    }
    inputs = {filename: plsc.sha256_file(out_dir / filename)
              for filename in ("settings.json", *_TABLE_FILES.values(), "accounting.json")}
    if (out_dir / PREFLIGHT_NAME).exists():
        inputs[PREFLIGHT_NAME] = plsc.sha256_file(out_dir / PREFLIGHT_NAME)
    result["inputs"] = inputs
    lsc._write_json(out_dir / DECISION_NAME, result)
    return result


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    def common(command: argparse.ArgumentParser) -> None:
        command.add_argument("--export-dir", type=Path, required=True,
                             help="The cohort export, read in place (ADNI: adni_sexdiff/results/motco_input)")
        command.add_argument("--out-dir", type=Path, required=True)

    def units(command: argparse.ArgumentParser) -> None:
        command.add_argument("--reps", type=int, default=REPLICATES,
                             help=f"Replicates per point (default {REPLICATES})")
        command.add_argument("--perms", type=int, default=PERMUTATIONS,
                             help=f"Permutations (default {PERMUTATIONS})")

    pre = sub.add_parser("preflight", help="Write the settings and check the cohort's fit against the contract")
    common(pre)
    units(pre)
    run = sub.add_parser("run", help="Evaluate (a shard of) the units")
    common(run)
    units(run)
    run.add_argument("--shard", type=lsc.parse_shard, default=(0, 1), help="i/N, round-robin over the units")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    merge_parser = sub.add_parser("merge", help="Merge the shards and write the tables")
    common(merge_parser)
    decide_parser = sub.add_parser("decide", help=f"Write {DECISION_NAME}")
    common(decide_parser)
    return parser


def _load(directory: Path) -> Cohort:
    try:
        return load_cohort(directory, CONTRACTS["adni"])
    except pc.ExportContractError as exc:
        raise contract_exit(exc) from exc


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command in ("preflight", "run"):
        cohort = _load(args.export_dir)
        settings = build_settings(cohort, replicates=args.reps, permutations=args.perms)
        if args.command == "preflight":
            result = preflight(args.out_dir, settings, cohort)
            print(json.dumps({"passed": result["passed"], "units": [
                {key: unit.get(key) for key in ("unit", "passed", "fidelity", "forced_rank", "statistics_differing",
                                                "p_values_differing", "error")}
                for unit in result["units"]
            ]}, indent=1))
            return 0 if result["passed"] else 1
        index, count = args.shard
        evaluated = run_shard(args.out_dir, settings, cohort, index, count, limit=args.limit, workers=args.workers)
        print(f"evaluated {evaluated} unit(s)", file=sys.stderr)
        return 0
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, args.export_dir), indent=1))
        return 0
    result = decide(args.out_dir, args.export_dir)
    print(json.dumps({
        "verdict": result["verdict"],
        "scope_start": result["scope_start"],
        "reference_scope": result["operating"]["reference"]["scope_start"],
        "candidate_scope": result["operating"]["candidate"]["scope_start"],
        "within_stage_scopes": {arm: entry["scope_start"] for arm, entry in result["within_stage_operating"].items()},
        "fixed_rank_grade": result["fixed_rank"]["grade"],
        "gate": result["gate"],
        "indicative": result["indicative"],
        "complete": result["complete"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
