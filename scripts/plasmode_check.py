#!/usr/bin/env python3
"""Plasmode attribution demonstration on the SEA-AD astrocyte export (no R).

The check behind ``openspec/changes/plasmode-attribution-demonstration/`` (spec
``plasmode-attribution-demonstration``, roadmap item 12): does ``motco
analyze``, at the case study's configuration, detect a group × stage
orientation difference planted in the real exported blocks, and does its
covariate-adjusted attribution trace the detection back to the planted
features? A predeclared rule either lifts Phase 6's withheld attribution leg
(``attribution_leg_lifted``) or reports it failed (``attribution_leg_failed``),
naming every failing criterion and cell.

**Input contract** (``--export-dir``, read in place, never copied into MOTCO):

- ``manifest.json`` declares exactly the blocks ``rna`` and ``atac`` (each a
  donors × features CSV with ``donor_id`` as its first column, ``identity``
  transform and the recorded ``n_features``), the metadata file with its
  ``sample_col``/``group_col``/``stage_col``, the groups ``Female``/``Male``
  and the stage order ``0``/``1``/``2``, and the covariates file;
- ``metadata.csv`` lists every donor once with a group among the two declared
  and a stage in the declared order, every group × stage cell observed;
- ``covariates.csv`` lists the same donors with ``age_at_death``, ``pmi``,
  ``rin``, ``years_education``, ``n_rna`` and ``n_atac``, none missing. The
  model's covariates are the case-study primary's six, in its order:
  ``age_at_death``, ``pmi``, ``rin``, ``years_education``, ``log10_n_rna`` and
  ``log10_n_atac``, the last two the base-10 logarithms of ``n_rna`` and
  ``n_atac`` (``sead_sexdiff`` ``ANALYSIS_LOG10_COLUMNS`` at ``2d427b0``);
- every block lists the metadata's donors; rows are taken in ``donor_id``
  order, as ``sead_sexdiff`` analyses them.

A violation raises a named :class:`ExportContractError` before any analysis.
The sha256 of the five files read is part of the settings, so a resume on
another export is refused.

**Fidelity gate** (``fidelity``): the latent space is stage-supervised and
covariate-blind, so it does not depend on the group labels. One analysis on a
within-stage permutation of the labels (stream ``fidelity``), with no
permutations and no attribution, must select rank 5 with mean held-out stage
AUROC 0.635556 to six decimals (the primary's, SEA-AD report §3, §10). Every
anchor unit is held to the same check.

**Null-making step.** Each replicate permutes the group labels within each
stage (cell table kept, each donor's blocks and covariates on its own row).
A draw that equals the real labels is refused and redrawn from the next seed
of the replicate's stream; the redraw is counted. No analysis ever runs under
the real labels, and no record holds a donor identifier or a label vector.

**Plant.** In the donors of the second declared group (``Male``) under the
permuted labels, each feature ``f`` of a planted set ``K`` (10 per block, drawn
per replicate, uniform without replacement, with uniform signs ``s_f``) is
shifted by ``e · c(stage) · σ_f · s_f``, with ``σ_f`` the export's pooled
within-stage SD (ddof = number of stages) and ``c`` = 0, 0.5, 1 at stages 0,
1, 2. At ``e = 0`` the blocks are the export's frames, untouched.

**Seeds** derive from (base seed 2034, stream ``paper`` | ``calibration`` |
``fidelity``, replicate, purpose ``labels`` | ``plant`` | ``rrpp`` |
``bootstrap``) and never from the strength, so the cells of a replicate share
the labels, ``K``, the signs, the RRPP seed and the bootstrap seed.

**Calibration** (``calibrate``): the ladder e ∈ {0.5, 1.0, 1.5, 2.0, 3.0},
40 units per step at 99 permutations with attribution off (stream
``calibration``). ``e*`` is the first step whose ``angle`` rejection rate is at
least 0.90; otherwise the cap (3.0, ``cap_reached``). The descriptive cell is
the step below ``e*`` (``e*``/2 at the lowest step). Every iteration is appended
to ``calibration.json``; ``run --calibration`` names the adopted one.

**Units and arms.** The anchor (e = 0), ``e*`` and the descriptive step, 500
replicates each, 999 permutations, B = 200, top-k 20. Each unit is analysed
once by :func:`motco.analysis.analyze_trajectories` with the six covariates,
covariate-adjusted attribution, rank floor 2, net displacement,
``residual_correction="df"`` and the harness CV defaults (``random_state``
1203): the ``adjusted`` (primary) arm. The ``unadjusted`` arm attributes on the
analysis's ``latent_artifacts`` without covariates, with the same bootstrap
seed (descriptive). Each unit first checks that the exposed estimator
reproduces the latent scores and that the adjusted means project to the stage
means (1e-8 relative); a failed check raises and the unit is not recorded.

**Metrics** per arm, block (transitions ``s->s+1`` and ``principal``) and
component (observed, PLS-captured, residual): the planted-set precision
|top-k ∩ K| / min(k, |K|), the mean bootstrap sign stability and top-k
frequency over all of K, and the share of recovered planted features whose
observed sign is the planted sign.

**Decision** (pure, :func:`plasmode_decision`, read from the committed
summary):

- (i) size: at the anchor, the rate (p ≤ 0.05) of each of ``delta``,
  ``angle`` and ``shape`` ≤ α + 2·√(α(1 − α)/R);
- (ii) detection: at ``e*``, the ``angle`` rate ≥ 0.80;
- (iii) recovery and stability: on the ``e*`` units whose ``angle`` rejects,
  primary arm, observed principal contrast, mean − 2·SE ≥ 0.50 (planted-set
  precision), 0.90 (sign stability over K) and 0.50 (top-k frequency over K).

The unadjusted arm, the captured and residual components, the transition
blocks, the unconditional readings, ``delta``/``shape`` at ``e*`` and the
descriptive cell are reported beside the verdict and never change it. A
verdict below 500 replicates per cell, 999 permutations or B = 200, or at a
strength not adopted from calibration, is indicative.

Outputs (``--out-dir``): ``fidelity.json``; ``calibration.json``;
``settings.json``; ``records/shard_<i>_of_<N>.jsonl`` (one line per unit;
resume skips recorded units); after ``merge``: ``merged.jsonl``,
``plasmode_units.csv`` and ``plasmode_metrics.csv`` (per unit, kept on the run
machine), ``plasmode_summary.csv`` and ``merge.json`` (aggregates); after
``decide``: ``plasmode_decision.json``.

With ``--workers`` above 1, pin BLAS to one thread per process
(``OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1``).

Examples:

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
    EXPORT=../sead_sexdiff/results/motco_input/cohort_30

    python scripts/plasmode_check.py fidelity --export-dir $EXPORT --out-dir results/pm
    python scripts/plasmode_check.py calibrate --export-dir $EXPORT --fidelity results/pm/fidelity.json \\
        --out-dir results/pm --workers 8

    # smoke at an explicit strength (indicative)
    python scripts/plasmode_check.py run --export-dir $EXPORT --fidelity results/pm/fidelity.json \\
        --e-star 1.0 --reps 20 --perms 99 --bootstrap 50 --out-dir results/pm-smoke --workers 8

    # paper grade, sharded, then merge and decide
    python scripts/plasmode_check.py run --export-dir $EXPORT --fidelity results/pm/fidelity.json \\
        --calibration results/pm/calibration.json --index 3 --count 30 --out-dir results/pm/run
    python scripts/plasmode_check.py merge --out-dir results/pm/run
    python scripts/plasmode_check.py decide --out-dir results/pm/run
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import adjusted_attribution_check as aa  # noqa: E402  (top-k rule and consistency error)
import covariate_type_i_check as d10  # noqa: E402  (settings files, JSON and CSV helpers)

sys.path.pop(0)

SCRIPT_VERSION = 1
CHECK = "plasmode"
BASE_SEED = 2034
STREAMS: tuple[str, ...] = ("paper", "calibration", "fidelity")
PURPOSES: tuple[str, ...] = ("labels", "plant", "rrpp", "bootstrap")
#: The calibration ladder (design D5), in within-stage SD units at the last stage.
LADDER: tuple[float, ...] = (0.5, 1.0, 1.5, 2.0, 3.0)
ADOPTION_RATE = 0.90
CALIBRATION_REPLICATES = 40
CALIBRATION_PERMUTATIONS = 99
DEFAULT_REPLICATES = 500
DEFAULT_PERMUTATIONS = 999
DEFAULT_BOOTSTRAP = 200
DEFAULT_TOP_K = 20
PLANTED_PER_BLOCK = 10
ZERO_TOLERANCE = 1e-12
PROJECTION_RTOL = 1e-8
MAX_LABEL_ATTEMPTS = 100

#: The export's declared layout (design D2).
BLOCKS: tuple[str, ...] = ("rna", "atac")
GROUPS: tuple[str, str] = ("Female", "Male")
STAGE_ORDER: tuple[str, ...] = ("0", "1", "2")
#: The case-study primary's covariates, in model order (``sead_sexdiff`` ``ANALYSIS_COVARIATES``).
COVARIATES: tuple[str, ...] = ("age_at_death", "pmi", "rin", "years_education", "log10_n_rna", "log10_n_atac")
#: Derived column → exported count it is the base-10 logarithm of (``ANALYSIS_LOG10_COLUMNS``).
LOG10_COLUMNS: dict[str, str] = {"log10_n_rna": "n_rna", "log10_n_atac": "n_atac"}
DERIVATION_SOURCE = "sead_sexdiff src/sead_sexdiff/data.py ANALYSIS_LOG10_COLUMNS at 2d427b0"
MISSING_COVARIATES = "error"
#: The harness CV defaults the primary used (``random_state`` 1203).
CV: dict[str, Any] = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 5, "max_components": 20, "random_state": 1203}
MIN_COMPONENTS = 2
ORIENTATION_ESTIMATOR = "net_displacement"
RESIDUAL_CORRECTION = "df"
#: The primary's label-free latent space (design D6): CV-selected rank and mean held-out stage AUROC.
FIDELITY_GATE: dict[str, Any] = {"selected_lv": 5, "cv_mean_auroc": 0.635556, "decimals": 6}

ANCHOR = 0.0
ARMS: tuple[str, ...] = ("adjusted", "unadjusted")
PRIMARY_ARM = "adjusted"
COMPONENTS: tuple[str, ...] = ("observed", "pls_captured", "residual")
PRINCIPAL = "principal"
STATISTICS: tuple[str, ...] = ("delta", "angle", "shape")
ALPHA = 0.05
SE_MULTIPLIER = 2.0
#: The roadmap's power floor (``adni_design_check.POWER_FLOOR``).
DETECTION_FLOOR = 0.80
#: Criterion (iii)'s thresholds on mean − 2·SE (design D9); part of the settings hash.
THRESHOLDS: dict[str, float] = {"planted_precision": 0.50, "sign_stability": 0.90, "top_k_frequency": 0.50}
RULE_VERSION = 1
METRICS: tuple[str, ...] = (
    "planted_precision",
    "sign_stability",
    "top_k_frequency",
    "sign_agreement",
    "valid_replicates",
)
SUBSETS: tuple[str, ...] = ("all", "angle_rejects")
VERDICTS: tuple[str, ...] = ("attribution_leg_lifted", "attribution_leg_failed")
CONSEQUENCES: dict[str, str] = {
    "attribution_leg_lifted": (
        "the Phase 6 G1 attribution leg is recorded as demonstrated on planted data (SEA-AD astrocyte cohort_30); "
        "the paper cites the report; no case-study reading changes"
    ),
    "attribution_leg_failed": (
        "the Phase 6 G1 attribution leg is reported failed, naming the failing criteria; Phase 6 is not reopened "
        "and no case-study reading changes"
    ),
}


# ── Input contract ───────────────────────────────────────────────────────────


class ExportContractError(ValueError):
    """The export violates the input contract; nothing is analysed."""

    violation = "contract"


class ManifestError(ExportContractError):
    violation = "manifest"


class DonorMismatchError(ExportContractError):
    violation = "donor_mismatch"


class UnknownStageError(ExportContractError):
    violation = "unknown_stage"


class GroupCountError(ExportContractError):
    violation = "group_count"


class EmptyCellError(ExportContractError):
    violation = "empty_cell"


class MissingCovariateError(ExportContractError):
    violation = "missing_covariate"


class CovariateValueError(ExportContractError):
    violation = "covariate_value"


class BlockValueError(ExportContractError):
    violation = "block_value"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass
class Export:
    """The loaded export: blocks and metadata in ``donor_id`` order, with the files' sha256."""

    directory: Path
    manifest: dict[str, Any]
    blocks: dict[str, pd.DataFrame]
    metadata: pd.DataFrame
    files: dict[str, str]
    sample_col: str
    group_col: str
    stage_col: str
    groups: tuple[str, str] = GROUPS
    stage_order: tuple[str, ...] = STAGE_ORDER
    _cache: dict[str, Any] = field(default_factory=dict, repr=False)

    @property
    def stages(self) -> np.ndarray:
        return self.metadata[self.stage_col].astype(str).to_numpy()

    @property
    def real_labels(self) -> np.ndarray:
        return self.metadata[self.group_col].astype(str).to_numpy()

    def group_stage_sizes(self, labels: np.ndarray | None = None) -> list[list[int]]:
        labels = self.real_labels if labels is None else labels
        stages = self.stages
        return [
            [int(((labels == group) & (stages == stage)).sum()) for stage in self.stage_order] for group in self.groups
        ]

    def sigma(self) -> dict[str, np.ndarray]:
        """Each block's pooled within-stage SD (ddof = number of stages), label-free."""

        if "sigma" not in self._cache:
            stages = self.stages
            self._cache["sigma"] = {
                name: pooled_within_stage_sd(frame, stages, self.stage_order) for name, frame in self.blocks.items()
            }
        return self._cache["sigma"]

    def geometry(self) -> dict[str, np.ndarray]:
        """The analysed scale (pooled SD, ddof 0) and pooled displacement in joint feature order."""

        if "geometry" not in self._cache:
            values = np.concatenate([self.blocks[name].to_numpy(dtype=float) for name in BLOCKS], axis=1)
            scale = values.std(axis=0)
            scale[scale < 1e-12] = 1.0
            standardized = (values - values.mean(axis=0)) / scale
            stages = self.stages
            displacement = standardized[stages == self.stage_order[-1]].mean(axis=0) - standardized[
                stages == self.stage_order[0]
            ].mean(axis=0)
            self._cache["geometry"] = {"scale": scale, "displacement": displacement}
        return self._cache["geometry"]


def _require(condition: bool, error: type[ExportContractError], message: str) -> None:
    if not condition:
        raise error(message)


def _repeated(values: Sequence[str]) -> list[str]:
    counts = pd.Series(list(values), dtype=object).value_counts()
    return sorted(str(value) for value in counts[counts > 1].index)


def load_export(directory: Path) -> Export:
    """Read the export through the input contract; raise the named error of the first violation."""

    directory = Path(directory)
    manifest_path = directory / "manifest.json"
    _require(manifest_path.is_file(), ManifestError, f"{manifest_path} does not exist.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    blocks_entry = manifest.get("blocks", {})
    _require(
        sorted(blocks_entry) == sorted(BLOCKS),
        ManifestError,
        f"{manifest_path} declares blocks {sorted(blocks_entry)}; the contract needs exactly {list(BLOCKS)}.",
    )
    for name in BLOCKS:
        entry = blocks_entry[name]
        _require(
            entry.get("motco_transform") == "identity",
            ManifestError,
            f"{manifest_path}: block {name!r} declares transform {entry.get('motco_transform')!r}, not 'identity'.",
        )
        _require(bool(entry.get("file")), ManifestError, f"{manifest_path}: block {name!r} names no file.")
    meta_entry = manifest.get("metadata", {})
    for key in ("file", "sample_col", "group_col", "stage_col"):
        _require(bool(meta_entry.get(key)), ManifestError, f"{manifest_path}: metadata entry lacks {key!r}.")
    _require(
        list(meta_entry.get("group_labels", [])) == list(GROUPS),
        ManifestError,
        f"{manifest_path} declares groups {meta_entry.get('group_labels')}; the contract needs {list(GROUPS)}.",
    )
    _require(
        [str(stage) for stage in meta_entry.get("stage_order", [])] == list(STAGE_ORDER),
        ManifestError,
        f"{manifest_path} declares stage order {meta_entry.get('stage_order')}; the contract needs "
        f"{list(STAGE_ORDER)}.",
    )
    covariates_file = manifest.get("covariates", {}).get("file")
    _require(bool(covariates_file), ManifestError, f"{manifest_path}: covariates entry names no file.")
    sample_col, group_col, stage_col = (str(meta_entry[key]) for key in ("sample_col", "group_col", "stage_col"))

    names = {
        "manifest": "manifest.json",
        "metadata": str(meta_entry["file"]),
        "covariates": str(covariates_file),
        **{name: str(blocks_entry[name]["file"]) for name in BLOCKS},
    }
    for role, name in names.items():
        _require((directory / name).is_file(), ManifestError, f"{directory / name} ({role}) does not exist.")
    files = {name: sha256_file(directory / name) for name in sorted(set(names.values()))}

    metadata = pd.read_csv(directory / names["metadata"], dtype={sample_col: str, group_col: str, stage_col: str})
    absent = [column for column in (sample_col, group_col, stage_col) if column not in metadata.columns]
    _require(not absent, ManifestError, f"{directory / names['metadata']} lacks column(s) {absent}.")
    donors = metadata[sample_col].astype(str).tolist()
    repeated = _repeated(donors)
    _require(not repeated, DonorMismatchError, f"{directory / names['metadata']} repeats donor(s) {repeated}.")
    observed_groups = sorted(set(metadata[group_col].astype(str)))
    _require(
        observed_groups == sorted(GROUPS),
        GroupCountError,
        f"{directory / names['metadata']} holds groups {observed_groups}; the contract needs exactly the two "
        f"declared groups {list(GROUPS)}.",
    )
    unknown = sorted(set(metadata[stage_col].astype(str)) - set(STAGE_ORDER))
    _require(
        not unknown,
        UnknownStageError,
        f"{directory / names['metadata']} holds stage label(s) {unknown} outside the declared order "
        f"{list(STAGE_ORDER)}.",
    )
    empty = [
        f"{group}/{stage}"
        for group in GROUPS
        for stage in STAGE_ORDER
        if not ((metadata[group_col] == group) & (metadata[stage_col].astype(str) == stage)).any()
    ]
    _require(not empty, EmptyCellError, f"{directory / names['metadata']} has no donor in cell(s) {empty}.")

    blocks: dict[str, pd.DataFrame] = {}
    for name in BLOCKS:
        path = directory / names[name]
        frame = pd.read_csv(path, index_col=0)
        frame.index = frame.index.astype(str)
        repeated = _repeated(frame.index.tolist())
        _require(not repeated, DonorMismatchError, f"{path} repeats donor(s) {repeated}.")
        _require(
            set(frame.index) == set(donors),
            DonorMismatchError,
            f"{path} and {names['metadata']} list different donors (only in the block: "
            f"{sorted(set(frame.index) - set(donors))}; only in the metadata: "
            f"{sorted(set(donors) - set(frame.index))}).",
        )
        expected = blocks_entry[name].get("n_features")
        _require(
            expected is None or int(expected) == frame.shape[1],
            ManifestError,
            f"{manifest_path} records n_features = {expected} for {name} but {path} has {frame.shape[1]}.",
        )
        try:
            values = frame.to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise BlockValueError(f"{path} holds non-numeric values.") from exc
        _require(bool(np.isfinite(values).all()), BlockValueError, f"{path} holds missing or infinite values.")
        blocks[name] = frame

    covariates = pd.read_csv(directory / names["covariates"], dtype={sample_col: str})
    _require(sample_col in covariates.columns, ManifestError, f"{names['covariates']} lacks {sample_col!r}.")
    cov_donors = covariates[sample_col].astype(str).tolist()
    repeated = _repeated(cov_donors)
    _require(not repeated, DonorMismatchError, f"{names['covariates']} repeats donor(s) {repeated}.")
    _require(
        set(cov_donors) == set(donors),
        DonorMismatchError,
        f"{names['covariates']} and {names['metadata']} list different donors (only in the covariates: "
        f"{sorted(set(cov_donors) - set(donors))}; only in the metadata: {sorted(set(donors) - set(cov_donors))}).",
    )
    needed = [name for name in COVARIATES if name not in LOG10_COLUMNS] + list(LOG10_COLUMNS.values())
    absent = [name for name in needed if name not in covariates.columns]
    _require(not absent, MissingCovariateError, f"{names['covariates']} lacks covariate column(s) {absent}.")
    for name in needed:
        missing = covariates.loc[covariates[name].isna(), sample_col].astype(str).tolist()
        _require(
            not missing,
            MissingCovariateError,
            f"{names['covariates']}: covariate {name!r} is missing for donor(s) {sorted(missing)}.",
        )
    merged = metadata[[sample_col, group_col, stage_col]].merge(
        covariates[[sample_col, *needed]], on=sample_col, how="inner", validate="one_to_one"
    )
    for derived, source in LOG10_COLUMNS.items():
        counts = merged[source].astype(float)
        bad = merged.loc[counts <= 0, sample_col].astype(str).tolist()
        _require(not bad, CovariateValueError, f"{names['covariates']}: {source} must be positive; donor(s) {bad}.")
        merged[derived] = np.log10(counts)
    merged = merged[[sample_col, group_col, stage_col, *COVARIATES]].sort_values(sample_col).reset_index(drop=True)
    merged[stage_col] = merged[stage_col].astype(str)
    order = merged[sample_col].tolist()
    blocks = {name: frame if frame.index.tolist() == order else frame.loc[order] for name, frame in blocks.items()}
    return Export(
        directory=directory,
        manifest=manifest,
        blocks=blocks,
        metadata=merged,
        files=files,
        sample_col=sample_col,
        group_col=group_col,
        stage_col=stage_col,
    )


def contract_exit(exc: ExportContractError) -> SystemExit:
    return SystemExit(f"export contract violation ({exc.violation}): {exc}")


# ── Seeds and the null-making step ───────────────────────────────────────────


def _seed(*parts: Any) -> int:
    encoded = json.dumps([CHECK, *parts], separators=(",", ":")).encode("utf-8")
    return int(hashlib.sha256(encoded).hexdigest()[:8], 16) & 0x7FFFFFFF


def labels_seed(base_seed: int, stream: str, replicate: int, attempt: int) -> int:
    return _seed(int(base_seed), stream, int(replicate), "labels", int(attempt))


def unit_seeds(base_seed: int, stream: str, replicate: int) -> dict[str, int]:
    """The replicate's seeds; ``labels`` is the first draw's (a redraw takes the next attempt's)."""

    if stream not in STREAMS:
        raise ValueError(f"Unknown stream {stream!r}; allowed: {list(STREAMS)}.")
    seeds = {"labels": labels_seed(base_seed, stream, replicate, 0)}
    seeds.update({purpose: _seed(int(base_seed), stream, int(replicate), purpose) for purpose in PURPOSES[1:]})
    return seeds


def permute_within_stage(
    labels: np.ndarray, stages: np.ndarray, stage_order: Sequence[str], rng: np.random.Generator
) -> np.ndarray:
    """Permute the labels among the donors of each stage, stage by stage in declared order."""

    permuted = np.array(labels, dtype=object)
    for stage in stage_order:
        rows = np.flatnonzero(stages == stage)
        permuted[rows] = np.asarray(labels, dtype=object)[rows][rng.permutation(rows.size)]
    return permuted


def draw_labels(
    real: np.ndarray, stages: np.ndarray, stage_order: Sequence[str], base_seed: int, stream: str, replicate: int
) -> tuple[np.ndarray, int]:
    """A within-stage permutation that differs from the real labels, and how many draws were refused."""

    for attempt in range(MAX_LABEL_ATTEMPTS):
        rng = np.random.default_rng(labels_seed(base_seed, stream, replicate, attempt))
        drawn = permute_within_stage(real, stages, stage_order, rng)
        if not np.array_equal(drawn, np.asarray(real, dtype=object)):
            return drawn, attempt
    raise RuntimeError(f"{MAX_LABEL_ATTEMPTS} within-stage permutations all reproduced the real labels.")


# ── Plant ────────────────────────────────────────────────────────────────────


def pooled_within_stage_sd(frame: pd.DataFrame, stages: np.ndarray, stage_order: Sequence[str]) -> np.ndarray:
    """Per feature, √(Σ_s Σ_{i ∈ s} (x_i − x̄_s)² / (n − S)), S the number of stages."""

    values = frame.to_numpy(dtype=float)
    residuals = values.copy()
    for stage in stage_order:
        rows = stages == stage
        residuals[rows] -= values[rows].mean(axis=0)
    return np.sqrt((residuals**2).sum(axis=0) / (values.shape[0] - len(stage_order)))


def stage_weights(stage_order: Sequence[str]) -> dict[str, float]:
    """``c(stage)``: the stage's position, 0 at the first stage to 1 at the last."""

    last = len(stage_order) - 1
    return {str(stage): index / last for index, stage in enumerate(stage_order)}


def draw_plant(feature_names: Mapping[str, Sequence[str]], seed: int, per_block: int) -> dict[str, dict[str, Any]]:
    """Per block, ``per_block`` features (uniform without replacement, sorted) and their uniform signs."""

    rng = np.random.default_rng(seed)
    plant: dict[str, dict[str, Any]] = {}
    for name in BLOCKS:
        names = list(feature_names[name])
        indices = np.sort(rng.choice(len(names), size=per_block, replace=False))
        signs = rng.choice(np.array([-1, 1]), size=per_block)
        plant[name] = {
            "indices": [int(index) for index in indices],
            "features": [str(names[index]) for index in indices],
            "signs": [int(sign) for sign in signs],
        }
    return plant


def plant_blocks(
    export: Export, labels: np.ndarray, plant: Mapping[str, Mapping[str, Any]], effect: float
) -> dict[str, pd.DataFrame]:
    """The export's blocks with ``e · c(stage) · σ_f · s_f`` added to the planted group's donors on K.

    At ``effect == 0`` the export's frames are returned untouched; otherwise
    only the planted group's donors at stages with ``c > 0`` change, on K only.
    """

    if float(effect) == 0.0:
        return dict(export.blocks)
    weights = stage_weights(export.stage_order)
    c = np.array([weights[stage] for stage in export.stages], dtype=float)
    rows = np.flatnonzero((np.asarray(labels) == export.groups[1]) & (c > 0))
    sigma = export.sigma()
    planted: dict[str, pd.DataFrame] = {}
    for name, frame in export.blocks.items():
        indices = np.asarray(plant[name]["indices"], dtype=int)
        signs = np.asarray(plant[name]["signs"], dtype=float)
        values = frame.to_numpy(dtype=float, copy=True)
        shift = float(effect) * c[rows][:, None] * sigma[name][indices][None, :] * signs[None, :]
        values[np.ix_(rows, indices)] = values[np.ix_(rows, indices)] + shift
        planted[name] = pd.DataFrame(values, index=frame.index, columns=frame.columns)
    return planted


def planted_vectors(
    feature_names: Sequence[str], plant: Mapping[str, Mapping[str, Any]]
) -> tuple[np.ndarray, np.ndarray]:
    """The planted mask and the planted sign (0 off K) over joint feature names ``block__feature``."""

    signs = {
        f"{name}__{feature}": int(sign)
        for name in BLOCKS
        for feature, sign in zip(plant[name]["features"], plant[name]["signs"])
    }
    sign_vector = np.array([signs.get(str(name), 0) for name in feature_names], dtype=int)
    if int((sign_vector != 0).sum()) != len(signs):
        raise ValueError("The planted features are not all among the analysed features.")
    return sign_vector != 0, sign_vector


def plant_cosine(export: Export, plant: Mapping[str, Mapping[str, Any]]) -> float | None:
    """Cosine of the plant direction with the export's pooled displacement, in analysed coordinates."""

    geometry = export.geometry()
    sigma = export.sigma()
    direction = np.zeros_like(geometry["displacement"])
    offset = 0
    for name in BLOCKS:
        indices = np.asarray(plant[name]["indices"], dtype=int)
        direction[offset + indices] = sigma[name][indices] * np.asarray(plant[name]["signs"], dtype=float)
        offset += export.blocks[name].shape[1]
    direction = direction / geometry["scale"]
    norms = float(np.linalg.norm(direction) * np.linalg.norm(geometry["displacement"]))
    return float(direction @ geometry["displacement"] / norms) if norms > 0 else None


def build_unit(effect: float, replicate: int, settings: Mapping[str, Any], export: Export) -> dict[str, Any]:
    """One unit's labels, plant, blocks and metadata; everything drawn depends on the replicate only."""

    base_seed, stream = int(settings["base_seed"]), str(settings["stream"])
    seeds = unit_seeds(base_seed, stream, replicate)
    labels, redraws = draw_labels(export.real_labels, export.stages, export.stage_order, base_seed, stream, replicate)
    plant = draw_plant(
        {name: list(frame.columns.astype(str)) for name, frame in export.blocks.items()},
        seeds["plant"],
        int(settings["plant"]["per_block"]),
    )
    metadata = export.metadata.copy()
    metadata[export.group_col] = labels.astype(str)
    return {
        "effect": float(effect),
        "replicate": int(replicate),
        "seeds": seeds,
        "label_redraws": redraws,
        "plant": plant,
        "blocks": plant_blocks(export, labels, plant, effect),
        "metadata": metadata,
    }


# ── Settings ─────────────────────────────────────────────────────────────────


def export_settings(export: Export) -> dict[str, Any]:
    return {
        "cohort": export.manifest.get("cohort"),
        "files": dict(export.files),
        "n_donors": int(export.metadata.shape[0]),
        "blocks": {name: int(export.blocks[name].shape[1]) for name in BLOCKS},
        "sample_col": export.sample_col,
        "group_col": export.group_col,
        "stage_col": export.stage_col,
        "groups": list(export.groups),
        "stage_order": list(export.stage_order),
        "group_stage_sizes": export.group_stage_sizes(),
    }


def base_settings(export: Export, *, base_seed: int = BASE_SEED) -> dict[str, Any]:
    """Every constant and the export's hashes; the run and calibration settings add their units."""

    return {
        "script_version": SCRIPT_VERSION,
        "check": CHECK,
        "base_seed": int(base_seed),
        "export": export_settings(export),
        "covariates": {
            "declared": list(COVARIATES),
            "derived": dict(LOG10_COLUMNS),
            "derivation": "base-10 logarithm of the exported count",
            "derivation_source": DERIVATION_SOURCE,
            "missing_policy": MISSING_COVARIATES,
        },
        "cv": dict(CV),
        "min_components": MIN_COMPONENTS,
        "orientation_estimator": ORIENTATION_ESTIMATOR,
        "residual_correction": RESIDUAL_CORRECTION,
        "fidelity_gate": dict(FIDELITY_GATE),
        "null": {
            "kind": "group labels permuted within each stage; blocks and covariates stay on their donors' rows",
            "guard": "a draw equal to the real labels is refused and redrawn from the next seed of the stream",
        },
        "plant": {
            "per_block": PLANTED_PER_BLOCK,
            "group": GROUPS[1],
            "stage_weights": stage_weights(STAGE_ORDER),
            "shift": "e * c(stage) * sigma_f * s_f",
            "sigma": "pooled within-stage SD of the export (ddof = number of stages)",
            "signs": "uniform on {-1, +1}",
            "selection": "uniform without replacement per block, per replicate",
        },
        "seeds": {"streams": list(STREAMS), "purposes": list(PURPOSES), "derivation": "sha256 of [check, base seed, "
                  "stream, replicate, purpose]; never the strength"},
        "ladder": list(LADDER),
        "calibration_rule": {
            "replicates": CALIBRATION_REPLICATES,
            "permutations": CALIBRATION_PERMUTATIONS,
            "attribution": False,
            "adoption_rate": ADOPTION_RATE,
            "cap": max(LADDER),
            "descriptive": "the ladder step below e*, or e*/2 at the lowest step",
        },
        "rule": {
            "version": RULE_VERSION,
            "alpha": ALPHA,
            "se_multiplier": SE_MULTIPLIER,
            "size": "anchor rate of each of delta, angle and shape <= alpha + k*sqrt(alpha(1 - alpha)/R)",
            "detection_floor": DETECTION_FLOOR,
            "detection": "angle rate at e* >= detection_floor",
            "thresholds": dict(THRESHOLDS),
            "recovery": (
                "on the e* units whose angle rejects, adjusted arm, observed principal contrast: mean - k*SE >= "
                "threshold for planted_precision, sign_stability and top_k_frequency"
            ),
            "verdicts": list(VERDICTS),
            "indicative": (
                f"fewer than {DEFAULT_REPLICATES} replicates in any cell, {DEFAULT_PERMUTATIONS} permutations or "
                f"{DEFAULT_BOOTSTRAP} bootstrap replicates, or e* not adopted from calibration"
            ),
            "descriptive": (
                "the unadjusted arm, the PLS-captured and residual components, the transition blocks, the readings "
                "not conditioned on rejection, delta and shape at e*, and the descriptive cell; none changes the "
                "verdict"
            ),
        },
    }


def check_fidelity(path: Path, export: Export) -> dict[str, Any]:
    """The gate's result, refused unless it passed on this export, the gate and the configuration."""

    if not Path(path).exists():
        raise SystemExit(f"{path} does not exist; run fidelity first.")
    stored = json.loads(Path(path).read_text(encoding="utf-8"))
    if not stored.get("passed"):
        raise SystemExit(f"{path} did not pass the fidelity gate.")
    if stored.get("export_files") != export.files:
        raise SystemExit(f"{path} was written on another export (file hashes differ).")
    if stored.get("gate") != FIDELITY_GATE or stored.get("cv") != CV or stored.get("min_components") != MIN_COMPONENTS:
        raise SystemExit(f"{path} was written under another gate or configuration.")
    return {
        "selected_lv": int(stored["selected_lv"]),
        "cv_mean_auroc": round(float(stored["cv_mean_auroc"]), int(FIDELITY_GATE["decimals"])),
    }


def build_settings(
    export: Export,
    fidelity: Mapping[str, Any],
    *,
    replicates: int = DEFAULT_REPLICATES,
    permutations: int = DEFAULT_PERMUTATIONS,
    bootstrap: int = DEFAULT_BOOTSTRAP,
    top_k: int = DEFAULT_TOP_K,
    base_seed: int = BASE_SEED,
    calibration: Path | None = None,
    e_star: float | None = None,
    descriptive: float | None = None,
) -> dict[str, Any]:
    """The paper-grade (or smoke) run: the anchor, ``e*`` and the descriptive step at the given scale.

    ``e*`` and the descriptive step come from ``calibration`` (the adopted
    iteration is named), or are given explicitly (smoke and timing; the
    decision is then indicative).
    """

    if (calibration is None) == (e_star is None):
        raise SystemExit("Give exactly one of --calibration and --e-star.")
    settings = base_settings(export, base_seed=base_seed)
    if calibration is not None:
        adopted = adopted_calibration(Path(calibration), settings)
        strength = {
            "source": "calibration",
            "e_star": adopted["e_star"],
            "descriptive": adopted["descriptive"],
            "outcome": adopted["outcome"],
        }
        settings["adopted_calibration"] = adopted
    else:
        assert e_star is not None
        if float(e_star) <= 0:
            raise SystemExit("--e-star must be positive.")
        strength = {
            "source": "explicit",
            "e_star": float(e_star),
            "descriptive": None if descriptive is None else float(descriptive),
            "outcome": None,
        }
        settings["adopted_calibration"] = None
    settings.update(
        {
            "stream": "paper",
            "replicates": int(replicates),
            "permutations": int(permutations),
            "attribution": {
                "bootstrap_replicates": int(bootstrap),
                "top_k": int(top_k),
                "zero_tolerance": ZERO_TOLERANCE,
                "bootstrap_seed": "derived per replicate (shared by both arms and every cell)",
            },
            "strength": strength,
            "cells": cells_from_strength(strength),
            "fidelity": dict(fidelity),
        }
    )
    return settings


settings_hash = d10.settings_hash
ensure_settings = d10.ensure_settings
load_settings = d10.load_settings
shard_path = d10.shard_path


def cells_from_strength(strength: Mapping[str, Any]) -> list[float]:
    cells = [ANCHOR, float(strength["e_star"])]
    if strength.get("descriptive") is not None:
        cells.append(float(strength["descriptive"]))
    return cells


def cells(settings: Mapping[str, Any]) -> list[float]:
    return [float(effect) for effect in settings["cells"]]


def cell_label(effect: float) -> str:
    return "anchor" if float(effect) == ANCHOR else f"e={float(effect):g}"


def unit_key(effect: float, replicate: int) -> str:
    return json.dumps([float(effect), int(replicate)])


def units(settings: Mapping[str, Any]) -> list[tuple[float, int]]:
    """Every (effect, replicate), replicate-major, so a shard holds the matched cells of its replicates."""

    return [(effect, replicate) for replicate in range(int(settings["replicates"])) for effect in cells(settings)]


# ── Unit evaluation ──────────────────────────────────────────────────────────


ConsistencyError = aa.ConsistencyError
top_k_indices = aa.top_k_indices


def model_transform(model: Any, frame: pd.DataFrame) -> np.ndarray:
    return np.asarray(model.transform(frame.to_numpy(dtype=float)), dtype=float)


def check_consistency(result: Any) -> None:
    """The exposed estimator reproduces the latent scores; with attribution, the adjusted means project
    to the analysis's stage means (relative tolerance 1e-8)."""

    artifacts = result.latent_artifacts
    scores = result.latent_scores.to_numpy(dtype=float)
    transformed = model_transform(artifacts.model, artifacts.features)
    if not np.allclose(transformed, scores, rtol=0.0, atol=1e-10 * max(1.0, float(np.abs(scores).max()))):
        raise ConsistencyError("latent_scores: the exposed estimator does not reproduce the latent scores.")
    if result.attribution is None:
        return
    means = result.attribution.means
    projected = model_transform(artifacts.model, means[list(artifacts.features.columns)])
    stage_means = result.stage_means.drop(columns=["group", "stage"]).to_numpy(dtype=float)
    labels_match = list(means["group"]) == list(result.stage_means["group"]) and list(means["stage"]) == list(
        result.stage_means["stage"]
    )
    scale = max(1.0, float(np.abs(stage_means).max()))
    if not labels_match or not np.allclose(projected, stage_means, rtol=PROJECTION_RTOL, atol=PROJECTION_RTOL * scale):
        raise ConsistencyError("projection: the adjusted means do not project to the analysis's stage means.")


def check_anchor_latent(result: Any, fidelity: Mapping[str, Any]) -> None:
    """An anchor unit's latent space is the fidelity gate's (rank, and AUROC to the gate's decimals)."""

    latent = result.record["latent"]
    rank = int(latent["selected_lv"])
    auroc = latent.get("cv_mean_auroc")
    rounded = None if auroc is None else round(float(auroc), int(FIDELITY_GATE["decimals"]))
    if rank != int(fidelity["selected_lv"]) or rounded is None or not math.isclose(
        rounded, float(fidelity["cv_mean_auroc"]), rel_tol=0.0, abs_tol=1e-12
    ):
        raise ConsistencyError(
            f"anchor_latent: rank {rank} and AUROC {rounded} differ from the fidelity gate's "
            f"{fidelity['selected_lv']} and {fidelity['cv_mean_auroc']}."
        )


def analysis_params(
    settings: Mapping[str, Any], seeds: Mapping[str, int] | None, *, permutations: int, n_jobs: int = 1
) -> Any:
    from motco.analysis import TrajectoryAnalysisParams
    from motco.simulations.evaluation import AttributionDiagnosticSettings

    attribution = settings.get("attribution")
    return TrajectoryAnalysisParams(
        permutations=int(permutations),
        seed=None if seeds is None else seeds["rrpp"],
        n_jobs=n_jobs,
        orientation_estimator=settings["orientation_estimator"],
        min_components=int(settings["min_components"]),
        residual_correction=settings["residual_correction"],
        attribution=(
            AttributionDiagnosticSettings()
            if attribution is None or seeds is None
            else AttributionDiagnosticSettings(
                enabled=True,
                bootstrap_replicates=int(attribution["bootstrap_replicates"]),
                bootstrap_seed=seeds["bootstrap"],
                top_k=int(attribution["top_k"]),
                zero_tolerance=float(attribution["zero_tolerance"]),
            )
        ),
        **settings["cv"],
    )


def analyse(
    blocks: Mapping[str, pd.DataFrame], metadata: pd.DataFrame, export: Export, params: Any
) -> tuple[Any, Any]:
    """The case-study analysis of one assignment: the dataset and :func:`analyze_trajectories`' result."""

    from motco.analysis import OmicsBlock, TrajectoryDataset, analyze_trajectories

    dataset = TrajectoryDataset(
        blocks=[OmicsBlock(name, blocks[name], "identity") for name in BLOCKS],
        metadata=metadata,
        sample_col=export.sample_col,
        group_col=export.group_col,
        stage_col=export.stage_col,
        stage_order=list(export.stage_order),
        covariate_cols=list(COVARIATES),
        missing_covariates=MISSING_COVARIATES,
    )
    return dataset, analyze_trajectories(dataset, params)


def contrast_metrics(
    contrast: np.ndarray | None,
    *,
    planted: np.ndarray,
    signs: np.ndarray,
    bootstrap: Any,
    feature_names: Sequence[str],
    top_k: int,
) -> dict[str, Any]:
    """Recovery and stability of the planted set on one contrast.

    The stability readings come from the bootstrap and stay defined when the
    observed contrast is unavailable or misses K entirely.
    """

    record: dict[str, Any] = {name: None for name in METRICS}
    record["available"] = contrast is not None
    record["planted_size"] = int(planted.sum())
    record["valid_replicates"] = int(bootstrap.valid_replicates)
    sign = np.asarray(bootstrap.sign_stability, dtype=float)[planted]
    frequency = np.asarray(bootstrap.top_k_selection_frequency, dtype=float)[planted]
    record["sign_stability"] = float(np.nanmean(sign)) if np.isfinite(sign).any() else None
    record["top_k_frequency"] = float(np.nanmean(frequency)) if np.isfinite(frequency).any() else None
    if contrast is None or not planted.any():
        return record
    selected = top_k_indices(contrast, feature_names, top_k)
    hits = selected[planted[selected]]
    record["planted_precision"] = float(len(hits) / min(top_k, int(planted.sum())))
    if len(hits):
        record["sign_agreement"] = float(np.mean(np.sign(contrast[hits]) == signs[hits]))
    return record


def _blocks(result: Any) -> dict[str, Any]:
    blocks: dict[str, Any] = {record.transition_id: record for record in result.transitions}
    blocks[result.principal_orientation.component_id] = result.principal_orientation
    return blocks


def arm_metrics(result: Any, *, planted: np.ndarray, signs: np.ndarray, top_k: int) -> dict[str, Any]:
    """Per block and component, the planted-set metrics of one arm's attribution."""

    names = list(result.feature_names)
    bootstrap = {(record.transition_id, record.component): record for record in result.bootstrap_records}
    return {
        block_id: {
            component: contrast_metrics(
                getattr(block, component).directional_contrast,
                planted=planted,
                signs=signs,
                bootstrap=bootstrap[(block_id, component)],
                feature_names=names,
                top_k=top_k,
            )
            for component in COMPONENTS
        }
        for block_id, block in _blocks(result).items()
    }


def analysis_summary(result: Any) -> dict[str, Any]:
    record = result.record
    summary: dict[str, Any] = {
        "statistics": record["statistics"],
        "p_values": record["p_values"],
        "selected_lv": record["latent"]["selected_lv"],
        "cv_mean_auroc": record["latent"].get("cv_mean_auroc"),
    }
    if result.attribution is not None:
        principal = result.attribution.principal_orientation
        summary.update(
            {
                "degenerate": bool(principal.degenerate),
                "degeneracy_sources": list(principal.degeneracy_sources),
                "observed_eigengaps": list(principal.observed_eigengaps),
                "rank_deficient_replicates": result.attribution.config.rank_deficient_replicates,
            }
        )
    return summary


def evaluate_built(
    built: Mapping[str, Any], settings: Mapping[str, Any], export: Export, *, n_jobs: int = 1
) -> dict[str, Any]:
    """One analysis of a built unit, its consistency checks, and both arms' metrics when attribution is on."""

    from motco.stats.attribution import analyze_orientation_attribution

    params = analysis_params(settings, built["seeds"], permutations=int(settings["permutations"]), n_jobs=n_jobs)
    dataset, result = analyse(built["blocks"], built["metadata"], export, params)
    check_consistency(result)
    if float(built["effect"]) == ANCHOR and settings.get("fidelity") is not None:
        check_anchor_latent(result, settings["fidelity"])
    evaluated: dict[str, Any] = {"analysis": analysis_summary(result)}
    if result.attribution is None:
        return evaluated
    attribution = settings["attribution"]
    artifacts = result.latent_artifacts
    labels = pd.DataFrame(
        {
            "group": dataset.metadata[export.group_col].astype(str).to_numpy(),
            "stage": dataset.metadata[export.stage_col].astype(str).to_numpy(),
        },
        index=artifacts.features.index,
    )
    unadjusted = analyze_orientation_attribution(
        artifacts.features,
        labels,
        artifacts.model,
        groups=list(dataset.groups),
        stages=list(dataset.stages),
        bootstrap_replicates=int(attribution["bootstrap_replicates"]),
        bootstrap_seed=built["seeds"]["bootstrap"],
        top_k=int(attribution["top_k"]),
        zero_tolerance=float(attribution["zero_tolerance"]),
    )
    planted, signs = planted_vectors(list(artifacts.features.columns), built["plant"])
    top_k = int(attribution["top_k"])
    evaluated["arms"] = {
        arm: {
            "blocks": arm_metrics(arm_result, planted=planted, signs=signs, top_k=top_k),
            "mean_source": arm_result.config.mean_source,
        }
        for arm, arm_result in (("adjusted", result.attribution), ("unadjusted", unadjusted))
    }
    return evaluated


def evaluate_unit(
    effect: float, replicate: int, settings: Mapping[str, Any], export: Export, *, n_jobs: int = 1
) -> dict[str, Any]:
    """One record line: the unit, its seeds and plant, the analysis and the metrics. No donor-level value."""

    start = perf_counter()
    built = build_unit(effect, replicate, settings, export)
    evaluated = evaluate_built(built, settings, export, n_jobs=n_jobs)
    return d10._jsonable(
        {
            "key": unit_key(effect, replicate),
            "effect": float(effect),
            "replicate": int(replicate),
            "stream": settings["stream"],
            "settings_hash": settings_hash(settings),
            "seeds": built["seeds"],
            "label_redraws": built["label_redraws"],
            "planted": {
                name: {"features": built["plant"][name]["features"], "signs": built["plant"][name]["signs"]}
                for name in BLOCKS
            },
            "plant_cosine": plant_cosine(export, built["plant"]),
            **evaluated,
            "runtime_seconds": perf_counter() - start,
        }
    )


# ── Fidelity ─────────────────────────────────────────────────────────────────


class FidelityError(RuntimeError):
    """The gate's analysis does not reproduce the primary's latent space."""


def fidelity_analysis(export: Export, *, n_jobs: int = 1) -> Any:
    """The gate's analysis: a within-stage permutation (stream ``fidelity``), no permutations, no attribution."""

    labels, _ = draw_labels(export.real_labels, export.stages, export.stage_order, BASE_SEED, "fidelity", 0)
    metadata = export.metadata.copy()
    metadata[export.group_col] = labels.astype(str)
    settings = {**base_settings(export), "attribution": None}
    _, result = analyse(export.blocks, metadata, export, analysis_params(settings, None, permutations=0, n_jobs=n_jobs))
    return result


def run_fidelity(export: Export, out_dir: Path, *, n_jobs: int = 1, analysis: Callable[..., Any] | None = None) -> dict:
    """Run the gate; write ``fidelity.json`` on a match, raise :class:`FidelityError` naming both values otherwise."""

    result = (analysis or fidelity_analysis)(export, n_jobs=n_jobs)
    latent = result.record["latent"]
    rank = int(latent["selected_lv"])
    auroc = latent.get("cv_mean_auroc")
    decimals = int(FIDELITY_GATE["decimals"])
    rounded = None if auroc is None else round(float(auroc), decimals)
    rank_ok = rank == int(FIDELITY_GATE["selected_lv"])
    auroc_ok = rounded is not None and math.isclose(
        rounded, float(FIDELITY_GATE["cv_mean_auroc"]), rel_tol=0.0, abs_tol=1e-12
    )
    if auroc is None or not (rank_ok and auroc_ok):
        raise FidelityError(
            f"fidelity gate failed: rank {rank} (gate {FIDELITY_GATE['selected_lv']}), AUROC "
            f"{'n/a' if rounded is None else format(rounded, f'.{decimals}f')} (gate "
            f"{FIDELITY_GATE['cv_mean_auroc']:.{decimals}f}); no calibration or unit may run."
        )
    payload = {
        "check": CHECK,
        "script_version": SCRIPT_VERSION,
        "passed": True,
        "gate": dict(FIDELITY_GATE),
        "selected_lv": rank,
        "cv_mean_auroc": float(auroc),
        "cv_mean_auroc_rounded": rounded,
        "component_selection": latent.get("component_selection"),
        "cv_repeats_completed": latent.get("cv_repeats_completed"),
        "cv": dict(CV),
        "min_components": MIN_COMPONENTS,
        "export_files": dict(export.files),
        "export_cohort": export.manifest.get("cohort"),
        "assignment": {"stream": "fidelity", "replicate": 0, "labels_seed": labels_seed(BASE_SEED, "fidelity", 0, 0)},
        "label_free": "the latent space is stage-supervised and covariate-blind; no statistic is read",
        "software": result.record.get("software"),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "fidelity.json").write_text(
        json.dumps(d10._jsonable(payload), indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    return payload


# ── Shards and resume ────────────────────────────────────────────────────────


def read_records(paths: Iterable[Path]) -> dict[str, dict[str, Any]]:
    """Records keyed by unit; torn lines (an interrupted write) and partial records are re-run."""

    records: dict[str, dict[str, Any]] = {}
    for path in paths:
        with path.open(encoding="utf-8") as handle:
            for line in handle:
                try:
                    record = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if set(record.get("arms", {})) == set(ARMS):
                    records[str(record["key"])] = record
    return records


def pending_units(out_dir: Path, settings: Mapping[str, Any], index: int, count: int) -> list[tuple[float, int]]:
    """This shard's units not yet recorded in any shard file; rewrites this shard's file with its complete lines."""

    path = shard_path(out_dir, index, count)
    path.parent.mkdir(parents=True, exist_ok=True)
    own = read_records([path]) if path.exists() else {}
    with path.open("w", encoding="utf-8") as handle:
        for key in sorted(own):
            handle.write(json.dumps(own[key], sort_keys=True, allow_nan=False) + "\n")
    others = [other for other in sorted(path.parent.glob("shard_*.jsonl")) if other != path]
    done = {**read_records(others), **own}
    return [
        unit
        for position, unit in enumerate(units(settings))
        if position % count == index and unit_key(*unit) not in done
    ]


_WORKER_EXPORT: Export | None = None


def _init_worker(directory: str, files: Mapping[str, str]) -> None:
    global _WORKER_EXPORT
    export = load_export(Path(directory))
    if export.files != dict(files):
        raise SystemExit(f"{directory} changed after the run started (file hashes differ).")
    _WORKER_EXPORT = export


def _evaluate_pending(unit: tuple[float, int], settings: Mapping[str, Any]) -> dict[str, Any]:
    assert _WORKER_EXPORT is not None
    return evaluate_unit(*unit, settings, _WORKER_EXPORT)


def run_shard(
    out_dir: Path,
    settings: Mapping[str, Any],
    export: Export,
    index: int,
    count: int,
    *,
    n_jobs: int = 1,
    workers: int = 1,
    limit: int | None = None,
) -> int:
    """Evaluate this shard's missing units, appending one line each; returns how many were evaluated."""

    if export.files != settings["export"]["files"]:
        raise SystemExit("The export's file hashes differ from the settings'; refusing to append.")
    pending = pending_units(out_dir, settings, index, count)
    if limit is not None:
        pending = pending[:limit]
    path = shard_path(out_dir, index, count)

    def append(position: int, record: Mapping[str, Any]) -> None:
        with path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
            handle.flush()
        print(
            f"[shard {index}/{count}] {position}/{len(pending)} {record['key']} {record['runtime_seconds']:.1f}s",
            file=sys.stderr,
            flush=True,
        )

    if workers > 1 and pending:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker, initargs=(str(export.directory), dict(export.files))
        ) as pool:
            for position, record in enumerate(pool.map(_evaluate_pending, pending, [settings] * len(pending)), 1):
                append(position, record)
    else:
        for position, unit in enumerate(pending, start=1):
            append(position, evaluate_unit(*unit, settings, export, n_jobs=n_jobs))
    return len(pending)


# ── Merge ────────────────────────────────────────────────────────────────────


UNIT_COLUMNS: tuple[str, ...] = (
    "effect",
    "replicate",
    *(f"{statistic}_statistic" for statistic in STATISTICS),
    *(f"{statistic}_p_value" for statistic in STATISTICS),
    "selected_lv",
    "cv_mean_auroc",
    "degenerate",
    "plant_cosine",
    "label_redraws",
    "rank_deficient_replicates",
)
METRIC_COLUMNS: tuple[str, ...] = (
    "effect",
    "replicate",
    "arm",
    "block",
    "component",
    "available",
    "planted_size",
    *METRICS,
    "angle_p_value",
)
SUMMARY_COLUMNS: tuple[str, ...] = (
    "effect",
    "cell",
    "arm",
    "block",
    "component",
    "subset",
    "quantity",
    "n",
    "mean",
    "se",
    "role",
)
#: Per-unit analysis readings summarized under the ``analysis`` pseudo-arm.
UNIT_QUANTITIES: tuple[str, ...] = (
    *(f"reject_{statistic}" for statistic in STATISTICS),
    "selected_lv",
    "cv_mean_auroc",
    "degenerate",
    "plant_cosine",
    "label_redraws",
)


def unit_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(records):
        record = records[key]
        analysis = record["analysis"]
        rows.append(
            {
                "effect": float(record["effect"]),
                "replicate": int(record["replicate"]),
                **{f"{s}_statistic": analysis["statistics"].get(s) for s in STATISTICS},
                **{f"{s}_p_value": analysis["p_values"].get(s) for s in STATISTICS},
                "selected_lv": analysis["selected_lv"],
                "cv_mean_auroc": analysis.get("cv_mean_auroc"),
                "degenerate": analysis.get("degenerate"),
                "plant_cosine": record.get("plant_cosine"),
                "label_redraws": record.get("label_redraws"),
                "rank_deficient_replicates": analysis.get("rank_deficient_replicates"),
            }
        )
    return rows


def metric_rows(records: Mapping[str, Mapping[str, Any]]) -> list[dict[str, Any]]:
    """One row per unit × arm × block × component, sorted by unit key."""

    rows: list[dict[str, Any]] = []
    for key in sorted(records):
        record = records[key]
        for arm in ARMS:
            for block, components in record["arms"][arm]["blocks"].items():
                for component in COMPONENTS:
                    values = components[component]
                    rows.append(
                        {
                            "effect": float(record["effect"]),
                            "replicate": int(record["replicate"]),
                            "arm": arm,
                            "block": block,
                            "component": component,
                            "available": bool(values["available"]),
                            "planted_size": values["planted_size"],
                            **{name: values.get(name) for name in METRICS},
                            "angle_p_value": record["analysis"]["p_values"].get("angle"),
                        }
                    )
    return rows


def _mean_se(values: Iterable[Any]) -> tuple[int, float | None, float | None]:
    finite = [float(value) for value in values if value is not None and value != "" and math.isfinite(float(value))]
    if not finite:
        return 0, None, None
    mean = float(np.mean(finite))
    se = float(np.std(finite, ddof=1) / math.sqrt(len(finite))) if len(finite) > 1 else None
    return len(finite), mean, se


def _rejects(p_value: Any, alpha: float) -> bool:
    return p_value is not None and p_value != "" and float(p_value) <= alpha


def _role(
    effect: float, arm: str, block: str, component: str, subset: str, quantity: str, settings: Mapping[str, Any]
) -> str:
    """What the decision reads a summary row for (empty for the descriptive rows)."""

    principal_observed = arm == PRIMARY_ARM and block == PRINCIPAL and component == "observed"
    if principal_observed and effect == ANCHOR and subset == "all" and quantity == "planted_precision":
        return "false_recovery_reference"
    e_star = float(settings["strength"]["e_star"])
    if principal_observed and effect == e_star and subset == "angle_rejects" and quantity in THRESHOLDS:
        return "criterion_iii"
    if arm == "analysis" and effect == ANCHOR and quantity in {f"reject_{s}" for s in STATISTICS}:
        return "criterion_i"
    if arm == "analysis" and effect == e_star and quantity == "reject_angle":
        return "criterion_ii"
    return ""


def summary_rows(
    units_table: Sequence[Mapping[str, Any]],
    metrics_table: Sequence[Mapping[str, Any]],
    settings: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Per cell: the analysis readings, then per arm × block × component × subset every metric's mean and SE.

    The subsets are every unit (``all``) and the units whose ``angle`` p-value
    is at most α (``angle_rejects``).
    """

    alpha = float(settings["rule"]["alpha"])
    out: list[dict[str, Any]] = []
    for effect in cells(settings):
        cell_units = [row for row in units_table if float(row["effect"]) == effect]
        base = {"effect": effect, "cell": cell_label(effect), "block": "", "component": "", "subset": "all"}
        for quantity in UNIT_QUANTITIES:
            if quantity.startswith("reject_"):
                statistic = quantity.removeprefix("reject_")
                values = [
                    float(_rejects(row[f"{statistic}_p_value"], alpha))
                    for row in cell_units
                    if row.get(f"{statistic}_p_value") not in (None, "")
                ]
            elif quantity == "degenerate":
                values = [
                    float(bool(row["degenerate"])) for row in cell_units if row.get("degenerate") not in (None, "")
                ]
            else:
                values = [row.get(quantity) for row in cell_units]
            n, mean, se = _mean_se(values)
            out.append(
                {**base, "arm": "analysis", "quantity": quantity, "n": n, "mean": mean, "se": se,
                 "role": _role(effect, "analysis", "", "", "all", quantity, settings)}
            )
        cell_metrics = [row for row in metrics_table if float(row["effect"]) == effect]
        keys = sorted({(str(row["arm"]), str(row["block"]), str(row["component"])) for row in cell_metrics})
        keys.sort(key=lambda key: (ARMS.index(key[0]), key[1] != PRINCIPAL, key[1], COMPONENTS.index(key[2])))
        for arm, block, component in keys:
            rows = [
                row for row in cell_metrics if (row["arm"], row["block"], row["component"]) == (arm, block, component)
            ]
            for subset in SUBSETS:
                chosen = rows if subset == "all" else [row for row in rows if _rejects(row.get("angle_p_value"), alpha)]
                for quantity in METRICS:
                    n, mean, se = _mean_se(row.get(quantity) for row in chosen)
                    out.append(
                        {
                            "effect": effect,
                            "cell": cell_label(effect),
                            "arm": arm,
                            "block": block,
                            "component": component,
                            "subset": subset,
                            "quantity": quantity,
                            "n": n,
                            "mean": mean,
                            "se": se,
                            "role": _role(effect, arm, block, component, subset, quantity, settings),
                        }
                    )
    return out


def merge(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """Merge every shard; write the per-unit tables (local), the per-cell summary and ``merge.json``."""

    records = read_records(sorted((out_dir / "records").glob("shard_*.jsonl")))
    stale = sorted(key for key, record in records.items() if record.get("settings_hash") != settings_hash(settings))
    if stale:
        raise SystemExit(f"{len(stale)} record(s) were written under other settings, e.g. {stale[:3]}.")
    expected = [unit_key(*unit) for unit in units(settings)]
    unexpected = sorted(set(records) - set(expected))
    if unexpected:
        raise SystemExit(f"{len(unexpected)} record(s) are not units of these settings, e.g. {unexpected[:3]}.")
    with (out_dir / "merged.jsonl").open("w", encoding="utf-8") as handle:
        for key in sorted(records):
            handle.write(json.dumps(records[key], sort_keys=True, allow_nan=False) + "\n")
    units_table = unit_rows(records)
    metrics_table = metric_rows(records)
    d10._write_csv(out_dir / "plasmode_units.csv", units_table, UNIT_COLUMNS)
    d10._write_csv(out_dir / "plasmode_metrics.csv", metrics_table, METRIC_COLUMNS)
    summary = summary_rows(units_table, metrics_table, settings)
    d10._write_csv(out_dir / "plasmode_summary.csv", summary, SUMMARY_COLUMNS)
    reference = next((row for row in summary if row["role"] == "false_recovery_reference"), None)
    payload = {
        "check": CHECK,
        "settings_hash": settings_hash(settings),
        "units": len(records),
        "expected": len(expected),
        "missing": [key for key in expected if key not in records][:20],
        "units_by_cell": {
            cell_label(effect): sum(1 for record in records.values() if float(record["effect"]) == effect)
            for effect in cells(settings)
        },
        "label_redraws": sum(int(record.get("label_redraws") or 0) for record in records.values()),
        "false_recovery_reference": None
        if reference is None
        else {key: reference[key] for key in ("cell", "arm", "block", "component", "subset", "n", "mean", "se")},
        "local_files": {
            name: sha256_file(out_dir / name)
            for name in ("merged.jsonl", "plasmode_units.csv", "plasmode_metrics.csv")
        },
    }
    (out_dir / "merge.json").write_text(
        json.dumps(d10._jsonable(payload), indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8"
    )
    return payload


# ── Decision ─────────────────────────────────────────────────────────────────


def _value(raw: Any) -> float | None:
    return None if raw is None or raw == "" or (isinstance(raw, float) and math.isnan(raw)) else float(raw)


def _lookup(summary: Iterable[Mapping[str, Any]]) -> dict[tuple, Mapping[str, Any]]:
    return {
        (float(row["effect"]), str(row["arm"]), str(row["block"] or ""), str(row["component"] or ""),
         str(row["subset"]), str(row["quantity"])): row
        for row in summary
    }


def size_bound(alpha: float, replicates: int, k: float) -> float:
    return alpha + k * math.sqrt(alpha * (1 - alpha) / replicates)


def _reading(
    criterion: str, effect: float, quantity: str, row: Mapping[str, Any] | None, *, direction: str, bound: float | None,
    k: float = 0.0,
) -> dict[str, Any]:
    """One bound reading: ``upper`` holds when mean + k·SE ≤ bound, ``lower`` when mean − k·SE ≥ bound."""

    n = 0 if row is None else int(row["n"])
    mean = None if row is None else _value(row["mean"])
    se = None if row is None else _value(row["se"])
    entry = {
        "criterion": criterion,
        "cell": cell_label(effect),
        "quantity": quantity,
        "n": n,
        "mean": mean,
        "se": se,
        "k": k,
        "bound": bound,
        "direction": direction,
    }
    if mean is None or bound is None or (k > 0 and se is None):
        return {**entry, "value": None, "verdict": "unavailable"}
    value = mean + k * (se or 0.0) if direction == "upper" else mean - k * (se or 0.0)
    holds = value <= bound if direction == "upper" else value >= bound
    return {**entry, "value": value, "verdict": "holds" if holds else "fails"}


def recovery_readings(
    lookup: Mapping[tuple, Mapping[str, Any]],
    settings: Mapping[str, Any],
    effect: float,
    *,
    arm: str = PRIMARY_ARM,
    block: str = PRINCIPAL,
    component: str = "observed",
    subset: str = "angle_rejects",
) -> list[dict[str, Any]]:
    """Criterion (iii)'s three bound readings on one arm, block, component and subset."""

    k = float(settings["rule"]["se_multiplier"])
    return [
        _reading(
            "recovery",
            effect,
            quantity,
            lookup.get((effect, arm, block, component, subset, quantity)),
            direction="lower",
            bound=float(threshold),
            k=k,
        )
        for quantity, threshold in settings["rule"]["thresholds"].items()
    ]


def rate_readings(lookup: Mapping[tuple, Mapping[str, Any]], effect: float) -> dict[str, Any]:
    """Each statistic's rejection rate at one cell: units with a p-value, rate and its SE."""

    readings: dict[str, Any] = {}
    for statistic in STATISTICS:
        row = lookup.get((effect, "analysis", "", "", "all", f"reject_{statistic}"))
        readings[statistic] = None if row is None else {
            "n": int(row["n"]), "rate": _value(row["mean"]), "se": _value(row["se"])
        }
    return readings


def verdict_from_summary(summary: Iterable[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """Criteria (i), (ii) and (iii) and the verdict; an unavailable reading counts as failing."""

    lookup = _lookup(summary)
    rule = settings["rule"]
    alpha, k = float(rule["alpha"]), float(rule["se_multiplier"])
    e_star = float(settings["strength"]["e_star"])
    size = []
    for statistic in STATISTICS:
        row = lookup.get((ANCHOR, "analysis", "", "", "all", f"reject_{statistic}"))
        replicates = 0 if row is None else int(row["n"])
        size.append(
            _reading(
                "size",
                ANCHOR,
                statistic,
                row,
                direction="upper",
                bound=size_bound(alpha, replicates, k) if replicates > 0 else None,
            )
        )
    detection = [
        _reading(
            "detection",
            e_star,
            "angle",
            lookup.get((e_star, "analysis", "", "", "all", "reject_angle")),
            direction="lower",
            bound=float(rule["detection_floor"]),
        )
    ]
    recovery = recovery_readings(lookup, settings, e_star)
    criteria = [*size, *detection, *recovery]
    failing = [
        {key: entry[key] for key in ("criterion", "cell", "quantity", "value", "bound", "verdict")}
        for entry in criteria
        if entry["verdict"] != "holds"
    ]
    verdict = "attribution_leg_failed" if failing else "attribution_leg_lifted"
    return {
        "verdict": verdict,
        "consequence": CONSEQUENCES[verdict],
        "failing": failing,
        "criteria": {"size": size, "detection": detection, "recovery": recovery},
    }


def plasmode_decision(summary: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    """The verdict (:func:`verdict_from_summary`) with every descriptive reading beside it."""

    verdict = verdict_from_summary(summary, settings)
    lookup = _lookup(summary)
    e_star = float(settings["strength"]["e_star"])
    descriptive_effect = settings["strength"].get("descriptive")
    blocks = sorted({str(row["block"]) for row in summary if str(row["block"] or "") not in ("", PRINCIPAL)})
    replicates = {
        cell_label(effect): int((lookup.get((effect, "analysis", "", "", "all", "selected_lv")) or {"n": 0})["n"])
        for effect in cells(settings)
    }
    reasons = []
    short = {label: count for label, count in replicates.items() if count < DEFAULT_REPLICATES}
    if short:
        reasons.append(f"fewer than {DEFAULT_REPLICATES} replicates: {short}")
    if int(settings["permutations"]) < DEFAULT_PERMUTATIONS:
        reasons.append(f"{settings['permutations']} permutations < {DEFAULT_PERMUTATIONS}")
    if int(settings["attribution"]["bootstrap_replicates"]) < DEFAULT_BOOTSTRAP:
        reasons.append(f"B = {settings['attribution']['bootstrap_replicates']} < {DEFAULT_BOOTSTRAP}")
    if settings["strength"]["source"] != "calibration":
        reasons.append("e* was given explicitly, not adopted from calibration")

    def cell_readings(effect: float) -> dict[str, Any]:
        return {
            "rates": rate_readings(lookup, effect),
            "recovery_angle_rejects": recovery_readings(lookup, settings, effect),
            "recovery_all": recovery_readings(lookup, settings, effect, subset="all"),
        }

    def mean_of(effect: float, quantity: str, **where: str) -> dict[str, Any] | None:
        key = (effect, where.get("arm", "analysis"), where.get("block", ""), where.get("component", ""),
               where.get("subset", "all"), quantity)
        row = lookup.get(key)
        return None if row is None else {"n": int(row["n"]), "mean": _value(row["mean"]), "se": _value(row["se"])}

    principal_observed = {"arm": PRIMARY_ARM, "block": PRINCIPAL, "component": "observed"}
    descriptive = {
        "unadjusted_arm": recovery_readings(lookup, settings, e_star, arm="unadjusted"),
        "components": {
            component: recovery_readings(lookup, settings, e_star, component=component)
            for component in ("pls_captured", "residual")
        },
        "transitions": {block: recovery_readings(lookup, settings, e_star, block=block) for block in blocks},
        "unconditional": recovery_readings(lookup, settings, e_star, subset="all"),
        "e_star_rates": rate_readings(lookup, e_star),
        "descriptive_cell": None
        if descriptive_effect is None
        else {"cell": cell_label(float(descriptive_effect)), **cell_readings(float(descriptive_effect))},
        "false_recovery_reference": mean_of(ANCHOR, "planted_precision", **principal_observed),
        "sign_agreement": {
            subset: mean_of(e_star, "sign_agreement", **principal_observed, subset=subset) for subset in SUBSETS
        },
        "per_cell": {
            cell_label(effect): {
                quantity: mean_of(effect, quantity)
                for quantity in ("selected_lv", "cv_mean_auroc", "degenerate", "plant_cosine", "label_redraws")
            }
            for effect in cells(settings)
        },
    }
    return {
        "check": CHECK,
        "settings_hash": settings_hash(settings),
        "adopted_calibration": settings.get("adopted_calibration"),
        "strength": dict(settings["strength"]),
        "verdict": verdict["verdict"],
        "consequence": verdict["consequence"],
        "failing": verdict["failing"],
        "indicative": bool(reasons),
        "indicative_reasons": reasons,
        "replicates": replicates,
        "permutations": int(settings["permutations"]),
        "bootstrap_replicates": int(settings["attribution"]["bootstrap_replicates"]),
        "top_k": int(settings["attribution"]["top_k"]),
        "alpha": float(settings["rule"]["alpha"]),
        "se_multiplier": float(settings["rule"]["se_multiplier"]),
        "detection_floor": float(settings["rule"]["detection_floor"]),
        "thresholds": dict(settings["rule"]["thresholds"]),
        "criteria": verdict["criteria"],
        "descriptive": descriptive,
    }


def read_summary(path: Path) -> list[dict[str, Any]]:
    frame = pd.read_csv(path, dtype={"block": str, "component": str, "role": str}, keep_default_na=False)
    return [{str(key): value for key, value in row.items()} for row in frame.to_dict(orient="records")]


def decide(out_dir: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    path = out_dir / "plasmode_summary.csv"
    if not path.exists():
        raise SystemExit(f"{path} does not exist; run merge first.")
    decision = plasmode_decision(read_summary(path), settings)
    target = out_dir / "plasmode_decision.json"
    target.write_text(json.dumps(d10._jsonable(decision), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return decision


# ── Calibration ──────────────────────────────────────────────────────────────


def calibration_settings(export: Export, fidelity: Mapping[str, Any], *, base_seed: int = BASE_SEED) -> dict:
    """The calibration's settings: the ``calibration`` stream, 40 units per step, 99 permutations, no attribution."""

    settings = base_settings(export, base_seed=base_seed)
    settings.update(
        {
            "stream": "calibration",
            "replicates": CALIBRATION_REPLICATES,
            "permutations": CALIBRATION_PERMUTATIONS,
            "attribution": None,
            "fidelity": dict(fidelity),
            "cells": None,
        }
    )
    return settings


def step_settings(settings: Mapping[str, Any], effect: float) -> dict[str, Any]:
    return {**settings, "cells": [float(effect)]}


def _calibration_unit(unit: tuple[float, int], settings: Mapping[str, Any]) -> dict[str, Any]:
    assert _WORKER_EXPORT is not None
    return calibration_unit(*unit, settings, _WORKER_EXPORT)


def calibration_unit(effect: float, replicate: int, settings: Mapping[str, Any], export: Export) -> dict[str, Any]:
    built = build_unit(effect, replicate, settings, export)
    params = analysis_params(settings, built["seeds"], permutations=int(settings["permutations"]))
    _, result = analyse(built["blocks"], built["metadata"], export, params)
    check_consistency(result)
    return {"effect": float(effect), "replicate": int(replicate), "analysis": analysis_summary(result)}


def calibration_entry(evaluated: Sequence[Mapping[str, Any]], settings: Mapping[str, Any]) -> dict[str, Any]:
    alpha = float(settings["rule"]["alpha"])
    rates = {}
    for statistic in STATISTICS:
        values = [
            _rejects(record["analysis"]["p_values"].get(statistic), alpha)
            for record in evaluated
            if record["analysis"]["p_values"].get(statistic) is not None
        ]
        rates[statistic] = {"n": len(values), "rejections": int(sum(values)),
                            "rate": float(np.mean(values)) if values else None}
    ranks: dict[str, int] = {}
    for record in evaluated:
        rank = str(record["analysis"]["selected_lv"])
        ranks[rank] = ranks.get(rank, 0) + 1
    aurocs = [record["analysis"].get("cv_mean_auroc") for record in evaluated]
    n, mean, se = _mean_se(aurocs)
    return {
        "effect": float(settings["cells"][0]),
        "rate": rates["angle"]["rate"],
        "rates": rates,
        "units": len(evaluated),
        "selected_rank": dict(sorted(ranks.items(), key=lambda item: int(item[0]))),
        "cv_mean_auroc": {"n": n, "mean": mean, "se": se},
    }


def calibrate_once(settings: Mapping[str, Any], export: Export, *, workers: int = 1) -> dict[str, Any]:
    pending = [(float(settings["cells"][0]), replicate) for replicate in range(int(settings["replicates"]))]
    if workers > 1:
        with ProcessPoolExecutor(
            max_workers=workers, initializer=_init_worker, initargs=(str(export.directory), dict(export.files))
        ) as pool:
            evaluated = list(pool.map(_calibration_unit, pending, [settings] * len(pending)))
    else:
        evaluated = [calibration_unit(*unit, settings, export) for unit in pending]
    return calibration_entry(evaluated, settings)


def calibration_outcome(effect: float, rate: float | None) -> tuple[str, float | None, float | None]:
    """``adopt`` / ``cap_reached`` with (e*, descriptive), or ``climb``."""

    position = LADDER.index(float(effect))
    below = LADDER[position - 1] if position > 0 else float(effect) / 2
    if rate is not None and rate >= ADOPTION_RATE:
        return "adopt", float(effect), below
    if position == len(LADDER) - 1:
        return "cap_reached", float(effect), below
    return "climb", None, None


def calibrate(
    out_dir: Path,
    settings: Mapping[str, Any],
    export: Export,
    *,
    workers: int = 1,
    run_once: Callable[..., Mapping[str, Any]] = calibrate_once,
) -> dict[str, Any]:
    """Climb the ladder, appending every iteration to ``calibration.json``; resume after the last iteration."""

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "calibration.json"
    history: list[dict[str, Any]] = json.loads(path.read_text(encoding="utf-8")) if path.exists() else []
    digest = settings_hash(settings)
    if any(entry.get("calibration_hash") != digest for entry in history):
        raise SystemExit(f"{path} was written under other calibration settings or another export.")
    if history and history[-1]["outcome"] != "climb":
        return history[-1]
    start = LADDER.index(float(history[-1]["effect"])) + 1 if history else 0
    for effect in LADDER[start:]:
        entry = dict(run_once(step_settings(settings, effect), export, workers=workers))
        outcome, e_star, descriptive = calibration_outcome(effect, entry.get("rate"))
        entry.update(
            {
                "effect": float(effect),
                "iteration": len(history),
                "calibration_hash": digest,
                "base_seed": int(settings["base_seed"]),
                "stream": settings["stream"],
                "replicates": int(settings["replicates"]),
                "permutations": int(settings["permutations"]),
                "adoption_rate": ADOPTION_RATE,
                "export_files": dict(settings["export"]["files"]),
                "outcome": outcome,
                "e_star": e_star,
                "descriptive": descriptive,
            }
        )
        history.append(d10._jsonable(entry))
        path.write_text(json.dumps(history, indent=1, allow_nan=False) + "\n", encoding="utf-8")
        print(json.dumps({key: entry[key] for key in ("effect", "rate", "outcome", "e_star", "descriptive")}),
              file=sys.stderr, flush=True)
        if outcome != "climb":
            return entry
    raise AssertionError("the ladder ended without an outcome")


def adopted_calibration(path: Path, settings: Mapping[str, Any]) -> dict[str, Any]:
    """The terminal calibration iteration, refused if it was run on another export or did not finish."""

    if not path.exists():
        raise SystemExit(f"{path} does not exist; run calibrate first or give --e-star.")
    history = json.loads(path.read_text(encoding="utf-8"))
    if not history or history[-1]["outcome"] == "climb":
        raise SystemExit(f"{path} holds no finished calibration (the last iteration climbs).")
    entry = history[-1]
    if entry.get("export_files") != settings["export"]["files"]:
        raise SystemExit(f"{path} was calibrated on another export (file hashes differ).")
    if int(entry["base_seed"]) != int(settings["base_seed"]):
        raise SystemExit(f"{path} was calibrated under base seed {entry['base_seed']}.")
    digest = hashlib.sha256(json.dumps(entry, sort_keys=True).encode("utf-8")).hexdigest()
    return {
        "iteration": int(entry["iteration"]),
        "outcome": entry["outcome"],
        "e_star": float(entry["e_star"]),
        "descriptive": float(entry["descriptive"]),
        "rate": entry["rate"],
        "entry_sha256": digest,
    }


# ── CLI ──────────────────────────────────────────────────────────────────────


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plasmode attribution demonstration on the SEA-AD astrocyte export.")
    sub = parser.add_subparsers(dest="command", required=True)

    def export_arguments(command: argparse.ArgumentParser) -> None:
        command.add_argument("--export-dir", type=Path, required=True, help="sead_sexdiff export directory (cohort_30)")
        command.add_argument("--out-dir", type=Path, required=True, help="Output directory")

    fidelity = sub.add_parser("fidelity", help="Check the latent space reproduces the primary's (rank 5, AUROC)")
    export_arguments(fidelity)
    fidelity.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV")
    calibrate_command = sub.add_parser("calibrate", help="Climb the strength ladder on angle detection")
    export_arguments(calibrate_command)
    calibrate_command.add_argument("--fidelity", type=Path, required=True, help="fidelity.json of a passed gate")
    calibrate_command.add_argument("--base-seed", type=int, default=BASE_SEED, help=f"Base seed (default {BASE_SEED})")
    calibrate_command.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel")
    run = sub.add_parser("run", help="Evaluate (a shard of) the units")
    export_arguments(run)
    run.add_argument("--fidelity", type=Path, required=True, help="fidelity.json of a passed gate")
    strength = run.add_mutually_exclusive_group(required=True)
    strength.add_argument("--calibration", type=Path, help="calibration.json whose terminal iteration is adopted")
    strength.add_argument("--e-star", type=float, help="Explicit e* (smoke and timing; indicative)")
    run.add_argument("--descriptive", type=float, default=None, help="Explicit descriptive step (with --e-star)")
    run.add_argument("--reps", type=int, default=DEFAULT_REPLICATES, help=f"Replicates (default {DEFAULT_REPLICATES})")
    run.add_argument("--perms", type=int, default=DEFAULT_PERMUTATIONS, help=f"Permutations ({DEFAULT_PERMUTATIONS})")
    run.add_argument("--bootstrap", type=int, default=DEFAULT_BOOTSTRAP, help=f"Bootstrap B ({DEFAULT_BOOTSTRAP})")
    run.add_argument("--top-k", type=int, default=DEFAULT_TOP_K, help=f"Top-k (default {DEFAULT_TOP_K})")
    run.add_argument("--base-seed", type=int, default=BASE_SEED, help=f"Base seed (default {BASE_SEED})")
    run.add_argument("--index", type=int, default=0, help="Shard index in [0, count)")
    run.add_argument("--count", type=int, default=1, help="Number of shards")
    run.add_argument("--n-jobs", type=int, default=1, help="Jobs for CV and RRPP (not part of the settings)")
    run.add_argument("--workers", type=int, default=1, help="Units evaluated in parallel (not part of the settings)")
    run.add_argument("--limit", type=int, default=None, help="Evaluate at most this many pending units")
    for name, text in (("merge", "Merge every shard into the summary"), ("decide", "Write the decision")):
        command = sub.add_parser(name, help=f"{text} (reads settings.json)")
        command.add_argument("--out-dir", type=Path, required=True, help="Output directory")
    return parser


def _load(directory: Path) -> Export:
    try:
        return load_export(directory)
    except ExportContractError as exc:
        raise contract_exit(exc) from exc


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "fidelity":
        try:
            payload = run_fidelity(_load(args.export_dir), args.out_dir, n_jobs=args.n_jobs)
        except FidelityError as exc:
            raise SystemExit(str(exc)) from exc
        print(json.dumps({key: payload[key] for key in ("passed", "selected_lv", "cv_mean_auroc")}, indent=1))
        return 0
    if args.command == "calibrate":
        export = _load(args.export_dir)
        settings = calibration_settings(export, check_fidelity(args.fidelity, export), base_seed=args.base_seed)
        entry = calibrate(args.out_dir, settings, export, workers=args.workers)
        print(json.dumps({key: entry[key] for key in ("effect", "rate", "outcome", "e_star", "descriptive")}, indent=1))
        return 0
    if args.command == "run":
        if args.count < 1 or not 0 <= args.index < args.count:
            raise SystemExit(f"--index must be in [0, --count), got {args.index}/{args.count}.")
        if args.descriptive is not None and args.e_star is None:
            raise SystemExit("--descriptive goes with --e-star; a calibration names its own.")
        export = _load(args.export_dir)
        settings = build_settings(
            export,
            check_fidelity(args.fidelity, export),
            replicates=args.reps,
            permutations=args.perms,
            bootstrap=args.bootstrap,
            top_k=args.top_k,
            base_seed=args.base_seed,
            calibration=args.calibration,
            e_star=args.e_star,
            descriptive=args.descriptive,
        )
        ensure_settings(args.out_dir, settings)
        evaluated = run_shard(
            args.out_dir, settings, export, args.index, args.count,
            n_jobs=args.n_jobs, workers=args.workers, limit=args.limit,
        )
        print(f"evaluated {evaluated} unit(s) into {shard_path(args.out_dir, args.index, args.count)}", file=sys.stderr)
        return 0
    settings = load_settings(args.out_dir)
    if args.command == "merge":
        print(json.dumps(merge(args.out_dir, settings), indent=1))
        return 0
    decision = decide(args.out_dir, settings)
    print(json.dumps({key: decision[key] for key in ("verdict", "failing", "indicative", "replicates")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
