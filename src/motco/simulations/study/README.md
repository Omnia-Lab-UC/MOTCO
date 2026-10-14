# MOTCO trajectory power study

This package (`motco.simulations.study`) is the engine behind the
reproducible Type I error / power study for the MOTCO trajectory test
(`delta`, `angle`, `shape`). It enumerates a grid of semi-synthetic
datasets from the numpy generator (cached InterSIM reference data, no R
at runtime), runs the full MOTCO pipeline on each replicate, and
produces summary tables, figures, and an acceptance-target report.

This README is the operational handbook: enough for anyone
to run and replicate the study end-to-end, locally or on
a SLURM cluster.

---

## 1. What the study computes

For each cell in the grid × replicate, the runner:

1. Generates aligned methylation, expression, and proteomics matrices
   via the numpy generator and cached InterSIM reference data (no R).
2. Injects a group trajectory difference as feature-set surgery on the
   per-stage differential indicators, in one of the modes below.
3. Runs the MOTCO trajectory pipeline (integration → design →
   `estimate_difference` → RRPP).
4. Records per-statistic p-values (`delta`, `angle`, `shape`) plus a
   pre-registered combined rule.

**Trajectory modes**

| Mode          | Role                                                                  |
|---------------|-----------------------------------------------------------------------|
| `none`        | Type I baseline (identical groups). Always added by the enumerator.   |
| `translation` | Negative control — constant location offset (no geometry change).     |
| `magnitude`   | Power probe for `delta` (uniformly scales every step).                |
| `orientation` | Power probe for `angle` (global feature permutation = rotation).      |
| `shape`       | Power probe for `shape` (perturbs interior-stage overlaps).           |

Reports cover three views:

- **Type I table** — rejection rates on `none` cells (target ≈ α).
- **Specificity matrix** — off-diagonal rejection rates: each mode ×
  each statistic. Diagonal entries are power; off-diagonals should
  stay near α.
- **Power curves** — per-statistic rejection rate vs `effect_size` per
  mode.

Optional **acceptance targets** (Type I control, power monotonicity,
specificity) are evaluated against Monte Carlo uncertainty and saved
to a non-gating report.

---

## 2. Module layout

```
src/motco/simulations/study/
├── __init__.py        Public API re-exports
├── config.py          StudyConfig dataclasses + load/dump
├── enumerate.py       StudyConfig → SimulationGrid (cells)
├── sharding.py        Partition grid into shards, run_shard()
├── merge.py           Combine shard_*.jsonl → merged.jsonl
├── summary.py         Per-statistic + combined-rule summaries
├── report.py          Specificity matrix, power curves, Type I table, plots
└── targets.py         Acceptance-target evaluation
```

Driver scripts live at the repo root under `scripts/`:

- `scripts/run_study_shard.py` — one shard per invocation (one cluster array task).
- `scripts/motco_study.py` — `merge` + `report` subcommands for post-processing.
- `scripts/motco_study_array.sbatch` — SLURM array template.

Example configs live under `examples/trajectory_power_study/`:

- `smoke.json` — tiny end-to-end smoke (~minutes locally, not paper-grade).

---

## 3. Pipeline

```
config (YAML or JSON)
    │
    ▼  scripts/run_study_shard.py        ← one cluster array task per shard
    │     writes shard_<i>.jsonl (signature-guarded, resumable)
    ▼  scripts/motco_study.py merge      ← shard_*.jsonl → merged.jsonl
    │
    ▼  scripts/motco_study.py report     ← CSVs + PNGs + acceptance report
```

Two orthogonal parallelism layers:

1. **Across shards** — one SLURM array task per shard
   (`$SLURM_ARRAY_TASK_ID`).
2. **Within a shard** — `--n-jobs $SLURM_CPUS_PER_TASK` forwards to
   RRPP's permutation loop via the evaluation params.

---

## 4. Configuration

Configs are YAML or JSON; both go through `load_study_config()`. The
schema lives in `config.py` (`StudyConfig`). Required top-level keys:

| Field              | Purpose                                                                |
|--------------------|------------------------------------------------------------------------|
| `generator`        | Baseline numpy-generator params (sizing, `n_stages`, `p_dmp`, per-omic `delta_*`, perturbation). `seed` is required. |
| `evaluation`       | Integration method, RRPP permutations, `n_jobs`, eval seed.            |
| `trajectory_modes` | Modes enumerated in the power grid. `none` is always added.            |
| `effect_sizes`     | Non-negative effect-size sweep (per mode).                             |
| `axes`             | Optional OFAT axes. Keys must be namespaced `generator.*` or `evaluation.*`, or the one nested form `evaluation.integration_params.<key>` (see below). |
| `design_grid`      | Optional design points, each with its own anchored power grid. **Crossed** form `{"axes": {...}}`: every combination of the values is a point, and each axis must list its baseline value. **Explicit** form `{"points": [{axis: value, ...}, ...]}`: the points are listed directly (for non-rectangular columns such as a forced-rank PLS, a PCA and a concat column in one study); every point names the same axes and the baseline point must be listed. |
| `n_replicates`     | Replicates per cell.                                                   |
| `base_seed`        | Deterministic seed root for replicates.                                |
| `alpha`            | Significance level for rejection rates.                                |
| `acceptance`       | Pre-specified targets: `type_i`, `power`, `specificity` (each `type_i`/`specificity` entry takes an optional `bound`: `two_sided` \| `upper`, see *Acceptance target bound form* below), optional `gate`, optional `design_point` (advisory design-point rule over a `design_grid`), optional `rank_decision` (advisory retained-rank rule over the rank axis), optional `anchor_decision` (advisory Type I rule over every design point's zero-effect anchor). |
| `report_contract`  | Optional declared reporting/execution rules: `driver_component` (`observed` \| `pls_captured` \| `residual`; which attribution component `driver_report.csv` and the attribution figure present), `cross_replicate_driver_agreement` (only `descriptive` is legal — cross-replicate Jaccard/sign agreement is never a stability claim), `n_jobs_override` (`forbid` \| `warn`, default `warn`; `forbid` makes the shard runner refuse a differing `--n-jobs`). When declared, `report` also writes `report_contract.json` and `driver_report.csv`; when absent, every output is byte-identical to a pre-contract run. The Phase 5 profile declares `observed` / `descriptive` / `forbid`. |
| `metadata`         | Free-form provenance (name, intent, notes).                            |

Validation enforces:

- `trajectory_modes` ⊆ `{none, translation, magnitude, orientation, shape}`.
- `effect_sizes` are non-negative.
- `axes` keys use a known namespace prefix (`generator.*` or `evaluation.*`).
- `design_grid.axes` keys use the same namespaces, are not also OFAT axes, and
  each lists the baseline value; `acceptance.design_point.prefer` names only
  design-grid axes.
- `design_grid.points` (the explicit form) is not declared together with
  `axes`; every point names exactly the first point's axes; no two points are
  equal after value normalization (`"pc1"` and `null` are the same
  `evaluation.orientation_estimator`); one point equals the baseline. Errors
  name the offending point (`design_grid.points[i]`) or axis. `axes` is derived
  from the points, so every axis-based rule above applies unchanged.
- `acceptance.anchor_decision` has `0 < alpha < 1`, a positive
  `se_multiplier`, and a non-empty, duplicate-free `statistics` ⊆
  {`delta`, `angle`, `shape`}.
- `acceptance.rank_decision.axis` is a declared design-grid axis whose values
  include `null`; its target and protected modes are in `trajectory_modes`; both
  SE multipliers are non-negative.
- `0 < alpha < 1`.
- **The root mapping contains only the keys in the table above.** An unknown
  top-level key (a misspelled `report_contract`, say) fails to load with an
  error naming it instead of being silently ignored; every sub-block already
  rejected unknown fields.
- `report_contract`, when present, names a known `driver_component`, uses
  `cross_replicate_driver_agreement: descriptive`, and a known `n_jobs_override`.
- `acceptance.type_i[i].bound` and `acceptance.specificity[i].bound`, when
  present, are `two_sided` or `upper`; any other value fails to load with an
  error naming the entry and the value.

### Acceptance target bound form

Each Type I control (`acceptance.type_i`) and specificity
(`acceptance.specificity`) target judges its rejection rate in one of two forms.
`k` is the target's `se_tolerance` and `R` is the cell's available replicates.

| `bound` | Met when | SE |
|---|---|---|
| `two_sided` (default) | `\|rate − α\| ≤ k · SE_emp` | empirical, `√(rate(1 − rate)/R)` |
| `upper` | `rate ≤ α + k · √(α(1 − α)/R)` | α-based |

`upper` is the Phase 4 gate's control and Type I inflation bound
(`phase4.type_i_inflation_bound`, called directly), so an acceptance report
that uses it agrees with the gate. Under `two_sided` a *conservative* control is
flagged as a failure. A rate of 0 has an empirical SE of 0, so its tolerance
collapses to 0 and anything other than exactly α fails. Under `upper`, a
control that rejects less often than α is met.

Only an `upper` target adds observations: `bound_kind: "upper"`, `alpha_se`
(the α-based SE), `bound` (the upper limit), and a **signed** `deviation`
(`rate − α`, negative is conservative). Its rationale reads
`rate = 0.0240 <= alpha + 2.00·SE_alpha (0.0695)`. A `two_sided` target writes
exactly the observations and rationale it wrote before the key existed, with
`deviation = |rate − α|`. A report whose targets are all `two_sided` therefore
gains no `obs_bound_kind` / `obs_alpha_se` column and regenerates
byte-identically. In a mixed report those two columns are blank on the
two-sided rows, and `obs_bound_kind` marks the rows where `obs_deviation` is
signed. `dump_study_config` writes the resolved `bound` on every target.

No committed config sets `upper`. Flipping the default would change the
regenerated `acceptance_report.*` of every committed `results/` tree, so it would
require re-issuing every committed report.

### Nested evaluation axes

Both `axes` and `design_grid.axes` accept, besides top-level fields, exactly one
nested form: `evaluation.integration_params.<key>`. The value is set inside a
copy of the evaluation integration-parameter mapping (every other key
unchanged), so it enters the cell's evaluation parameters and therefore the
parameter signature. JSON `null` is a legal value meaning **key absent**: applying
it removes the key rather than storing `None`, so the baseline column's cells are
byte-identical (ids and signatures) to a config without the axis. The baseline
value of such an axis is the mapping's current value, or `null` when the key is
not set. Any other nesting (`generator.a.b`, `evaluation.attribution.x`) is
rejected by name.

The one use today is the retained PLS rank,
`evaluation.integration_params.forced_components`: `null` is the cross-validated
column, every other value fits the pooled PLS model at that fixed rank
(`component_selection = "forced"`, `selected_lv = <rank>` in the record's
integration metadata). Because the axis lives in the evaluation namespace,
design columns that differ only in it share the primary matched-seed family
**and** identical generator parameters — at every replicate index they evaluate
the same generated dataset under different measurement settings. The
duplicate-dataset guard accepts them on purpose (it keys on evaluation identity
too), and the report states the pairing wherever such columns are compared.

The second use is block selection, `evaluation.integration_params.layers`: a
list of omic blocks to measure (`null` = all three). Columns that differ only
in it evaluate the same generated dataset on different blocks — the Phase 6
three-block column is paired with the two-block baseline this way.

### Explicit group × stage sizes

`generator.group_stage_sizes` takes the generator's exact size table as nested
JSON lists (`[[11, 10, 28], [9, 10, 12]]`, rows = groups, columns = stages). The
loader normalizes it, in the baseline block and as an `axes` or
`design_grid.axes` value, to a tuple of int tuples, so equal tables give equal
cell ids, signatures and matched seeds whatever their JSON spelling. A config
that sets it must leave `n_samples`, `stage_sample_prop` and `group_ratio` out
(the generator rejects the combination). Hash payloads, and dumped configs, omit
the field while it is unset, so every config committed before it existed keeps
its cell ids, signatures and matched seeds
(`tests/test_committed_config_identity.py`).

See `examples/trajectory_power_study/smoke.json` for a complete,
minimal example.

### Scaling smoke → paper-grade

The smoke config is meant to finish in minutes; for a real study you
typically want:

| Field                       | Smoke | Paper-grade (typical) |
|-----------------------------|-------|-----------------------|
| `generator.n_samples`       | 60    | 200–600               |
| `evaluation.permutations`   | 49    | 999 or 4999           |
| `n_replicates`              | 8     | 500–1000              |
| `effect_sizes`              | 4 pts | 5–8 pts incl. `0.0`   |

Keep `0.0` as the first effect size — it anchors the within-mode null
and is what the specificity/Type I checks read.

**Phase 5 cost.** The committed paper-grade profile
`examples/trajectory_power_study/phase5_power_study.json` (n = 1200, four
stages, `p_dmp = 0.1`, pooled PLS with double-CV rank, 999 permutations,
attribution on the orientation cells) is 19 cells × 500 replicates = 9,500
units. A one-replicate rehearsal on 2026-09-09 (19 units as 19 concurrent
single-threaded processes on a 24-core workstation, BLAS pinned) measured a
median **43 s per unit** (min 30, max 51; CV fit + 999 RRPP permutations +
attribution) — the attribution units cost the same as the others (median 42.9 s
vs 43.0 s), so the bootstrap is immaterial at n = 1200. Budget roughly 120
workstation core-hours, or about 150 EPYC 7662 core-hours if the ladder's
per-core ratio (57 s per CV unit at 199 permutations) carries over; 100
single-CPU shards finish in about 1.5–2 h wall.

---

## 5. Local smoke run

From the repo root, with the project virtualenv active:

```bash
uv venv && source .venv/bin/activate
uv sync --extra test
```

Run all shards locally:

```bash
for i in 0 1 2 3; do
  python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/smoke.json \
    --out-dir /tmp/motco-smoke \
    --shard-index $i --n-shards 4 \
    --error-policy record
done

python scripts/motco_study.py merge  --out-dir /tmp/motco-smoke
python scripts/motco_study.py report \
    --config examples/trajectory_power_study/smoke.json \
    --out-dir /tmp/motco-smoke
```

Outputs (under `/tmp/motco-smoke/report/`):

- `specificity_matrix.csv` / `.png`
- `power_curves.csv` / `.png`
- `type_i_table.csv` / `type_i.png`
- `acceptance_report.csv` / `.json`

Run the smoke first on any new environment — it doubles as a sanity
check that the package (and its cached reference data) is wired up
correctly. No R is required at runtime.

---

## 6. SLURM cluster run

### 6.1 Prerequisites on the cluster

- Python (matching `pyproject.toml`) and `uv` (or `pip`) available on
  the compute nodes.
- **No R needed.** Generation runs on the numpy generator and the
  cached reference data (`src/motco/simulations/data/intersim_reference.npz`),
  which ships in the repo. R is only ever needed to *regenerate* that
  cache (see `export_reference.R`), not to run the study.
- The repo cloned and the virtualenv built once on a login node:
  `uv venv && uv sync --extra test`. The sbatch template activates
  `.venv/bin/activate` from the project root.

If your cluster uses `module load`/conda, add the appropriate `module
load python` (and any `conda activate`) lines before the `source
.venv/bin/activate` step in the sbatch script. No R module is required.

### 6.2 Choose `N_SHARDS`

Total replicates = `n_replicates × n_cells`, where `n_cells` is the
sum of (modes × effect_sizes) for the power grid plus one cell per
OFAT axis value. Aim for **~2–6 h per shard**:

```
N_SHARDS ≈ ceil(total_replicate_seconds / target_shard_seconds)
```

A short probe submission (`--array=0-0 --time=2:00:00`) gives a
realistic per-shard wallclock and memory profile. Use that to pick
both `N_SHARDS` and `--cpus-per-task`.

### 6.3 Edit the sbatch template

Open `scripts/motco_study_array.sbatch` and adjust:

- `#SBATCH --time`, `--cpus-per-task`, `--mem` — set from your probe.
- Partition / QOS / account flags (`-p`, `--qos`, `-A`) — cluster-specific.
- The `#SBATCH --array=0-N_SHARDS_MINUS_ONE` line is a placeholder;
  override with `--array=` on the command line at submit time.
- If your cluster needs `module load` / conda, add those calls before
  the `.venv/bin/activate` line.

The runner forwards `--n-jobs $SLURM_CPUS_PER_TASK` so RRPP saturates
the CPUs allocated to each task.

### 6.4 Submit

```bash
mkdir -p logs results

sbatch \
  --array=0-63 \
  --export=ALL,STUDY_CONFIG=$(pwd)/examples/trajectory_power_study/phase5_power_study.json,STUDY_OUT=$(pwd)/results,N_SHARDS=64 \
  scripts/motco_study_array.sbatch
```

Leave `STUDY_N_JOBS` unset. A config whose `report_contract` sets
`n_jobs_override: forbid` (the Phase 5 profile does) makes every array task
exit 2 before enumerating anything if `STUDY_N_JOBS` differs from its
`evaluation.n_jobs`; other configs warn and run at the overridden signature.

Required environment variables (via `--export`):

| Variable       | Meaning                                                         |
|----------------|-----------------------------------------------------------------|
| `STUDY_CONFIG` | Absolute path to the study config (YAML or JSON).               |
| `STUDY_OUT`    | Output directory for `shard_<i>.jsonl` files.                   |
| `N_SHARDS`     | Total shards. Must match the `--array` width.                   |

### 6.5 Merge + report after completion

```bash
python scripts/motco_study.py merge  --out-dir results
python scripts/motco_study.py report \
    --config examples/trajectory_power_study/phase5_power_study.json \
    --out-dir results
```

---

## 7. Failure handling and resumption

- **`--error-policy record`** (default in the sbatch template) lets a
  shard continue past a failed replicate; the failure is captured as
  a row with `status="failed"` in the shard JSONL. Use `raise` only
  when you want a single failure to abort the shard.
- **Signature-guarded resumability** — each `shard_<i>.jsonl` records
  the parameter signature it was produced with. Re-running the same
  shard skips replicates already completed for that signature and
  fills only what's missing or failed.
- **Resubmit only failed array tasks**:

  ```bash
  sbatch --array=7,12,40 \
         --export=ALL,STUDY_CONFIG=...,STUDY_OUT=...,N_SHARDS=64 \
         scripts/motco_study_array.sbatch
  ```

- **Forced overwrite** — pass `--overwrite` to `run_study_shard.py` to
  discard an existing shard JSONL before running. Use with care; this
  loses any completed replicates in that shard.
- **Worker-count lock** — `--n-jobs` enters the cell parameter signature
  and changes the realized permutation draws. When the config declares
  `report_contract.n_jobs_override: forbid`, `run_study_shard.py` refuses a
  `--n-jobs` that differs from `evaluation.n_jobs`: it prints both values
  and the signature consequence to stderr and exits **2** before enumerating
  or writing anything (an equal value is accepted). Under `warn`, or with no
  contract, it warns and proceeds at the overridden signature, as before.

---

## 8. Outputs

After `merge` + `report`, the output directory looks like:

```
results/
├── shard_0.jsonl            (per-shard raw records, signature-guarded)
├── shard_1.jsonl
├── ...
├── merged.jsonl             (deduplicated by (cell, replicate))
└── report/
    ├── specificity_matrix.csv     mode × statistic rejection rates
    ├── specificity_matrix.png
    ├── power_curves.csv           rejection rate vs effect_size, per mode × statistic
    ├── power_curves.png
    ├── type_i_table.csv           per-statistic + combined-rule on null cells
    ├── type_i.png
    ├── config_spectrum.csv        recorded latent eigengaps per cell × configuration
    ├── eigengap_stratified_power.csv  orientation power by eigengap tercile
    ├── continuity_resolved_orientation.csv  only when the grid sweeps baseline continuity
    ├── design_point_operating.csv only when the config declares a design_grid
    ├── design_point_power.png     (same condition)
    ├── design_point_decision.json only when acceptance.design_point is declared
    ├── design_point_decision.csv  (same condition)
    ├── rank_ladder.png            only when the design grid declares the retained-rank axis
    ├── rank_decision.json         only when acceptance.rank_decision is declared
    ├── rank_decision.csv          (same condition)
    ├── anchor_decision.json       only when acceptance.anchor_decision is declared
    ├── anchor_decision.csv        (same condition)
    ├── driver_report.csv          only when the config declares a report_contract: the declared
    │                              attribution component per (mode, effect, transition) — precision/recall
    │                              vs truth, selected count, within-replicate bootstrap stability,
    │                              replicate accounting; no cross-replicate agreement columns
    ├── report_contract.json       (same condition) resolved contract: driver component, the descriptive
    │                              status of cross-replicate agreement, the uniform n_jobs the records
    │                              carry, and the shared zero-effect anchor (cell_id, resolves_modes,
    │                              counted_as: 1) with one plain-language statement per item
    ├── acceptance_report.csv      acceptance target evaluation
    └── acceptance_report.json
```

With `acceptance.gate.enabled` the `phase4_*` tables and figures are written
as well (see "Phase 4 outputs" in the examples README); under a report
contract the attribution figure plots within-replicate bootstrap stability
only, for the declared component, while `phase4_attribution.csv` keeps every
component and the cross-replicate columns.

### Per-replicate record fields

Every JSONL record is one `SimulationReplicateResult`. Beyond the cell identity,
seeds, `p_values`, and `pair_statistics`, two fields carry the per-replicate
diagnostics:

- `realized_geometry` / `integration_metadata` / `attribution_diagnostics` — the
  Phase 4 diagnostics.
- `null_summary` — a compact per-statistic description of *that replicate's own*
  RRPP permutation null: `count` (retained draws, non-finite excluded), `mean`,
  `sd`, and `q50`/`q90`/`q95`/`q99`. It is written on every run with
  `permutations > 0`, independently of `include_null_distributions` (which
  controls only whether the full draw vectors are *returned*, never persisted).

  `null_summary` makes the pair *(observed statistic, its own critical value)*
  computable from a single record, which is what any pivotality or calibration
  question needs. It is **not** part of `parameter_signature` — it summarizes
  draws that already happened and changes no generation, integration, or
  permutation behavior — so records written before it existed still load (with an
  empty summary) and remain resumable. A record from a `permutations = 0` run is
  told apart from a pre-field record by `runtime_metadata["permutations"]`.

- `config_spectrum` — the eigenspectrum of the centered **stage-mean
  configuration in the evaluated latent space**, written on every run. It carries
  `pooled` (stages averaged across groups) and `groups` (one entry per group
  level), each with `n_points`, `n_dimensions`, `total_variance`, the normalized
  `spectrum`, and the **relative eigengap** `(l1 − l2) / Σl`. With
  `permutations > 0` it also carries `permutation_pooled_eigengap` — `count`,
  `mean`, `sd`, `q05`/`q50`/`q95` of the pooled eigengap over the permutation
  draws — so a replicate can be located against its own permutation
  distribution. Full per-permutation spectra are never retained.

  The eigengap is the covariate the [geometry
  audit](../../../../docs/reports/geometry-audit-2026-09-01.md) found predicts how
  wide a replicate's own `angle` permutation null will be, and therefore whether
  its orientation is resolvable at all. It is **recorded and reported only** — it
  enters no statistic, no p-value, and no decision rule. A configuration with
  zero total variance records `relative_eigengap: null` and an empty `spectrum`
  rather than a non-finite float; a two-stage configuration's eigengap is
  identically `1.0` and uninformative by construction.

  Unlike `null_summary`, the spectrum **is** signature-bearing:
  `parameter_signature` includes `config_spectrum_version`. A shard written
  before this field existed therefore **refuses resume** with the usual
  signature-mismatch error, rather than producing a merged set in which only some
  records carry the covariate — a missing covariate and a degenerate one are
  indistinguishable at analysis time. Already-merged pre-change record sets still
  *load* (the block is empty) and their committed reports are unaffected; only
  resuming *into* an old shard is refused.

### Eigengap reporting

`motco_study.py report` writes up to three spectrum tables beside the existing
ones, from recorded values only — no dataset is regenerated and no spectrum
recomputed:

- `config_spectrum.csv` — one row per (cell, configuration ∈ {`pooled`, group
  levels}) with `n_replicates` / `n_recorded` / `n_available` and the eigengap
  mean, sd, quartiles, min and max. `n_recorded` counts records carrying the
  block at all and `n_available` those whose eigengap is defined, so pre-field
  records and degenerate configurations stay separately visible.
- `eigengap_stratified_power.csv` — for orientation-mode power cells, rejection
  rates within **within-cell terciles** of the recorded pooled eigengap, with
  per-stratum counts and Monte Carlo SEs. A cell whose records carry no spectrum
  is emitted with `status = unavailable` rather than dropped.
- `continuity_resolved_orientation.csv` — **only when the merged set spans more
  than one `generator.baseline_continuity` value.** One row per (continuity,
  mode, effect, statistic) with the rejection rate, the recorded pooled-eigengap
  distribution (mean and terciles), and the dispersion of
  `null_summary["angle"]["q95"]`. The eigengap is the point of the table: power
  that rises along ρ should be traceable to a configuration that acquired a
  dominant axis, which is the observable that carries a continuity-conditioned
  orientation claim to real data — not the knob itself. A study that holds the
  axis fixed writes no such file, and records predating the axis are excluded
  rather than folded into the ρ = 0 bin. When the grid is a `design_grid`, rows
  are additionally keyed on every other design coordinate (emitted as columns),
  so a ρ × `n_samples` grid never pools three operating points into one row.

### Design-point reporting

A config with a `design_grid` enumerates one anchored power grid per design
point (phase `power_design`, metadata `design_point = {axis: value}`,
`varied_axis = "design_grid"`); the baseline point is the primary grid, stamped
with its coordinates. Design cells join the primary matched-seed family and are
invisible to every baseline reader (power curves, specificity matrix, Type I
table, Phase 4 gate, acceptance targets). The report then writes:

- `design_point_operating.csv` — one row per (design coordinates, mode, effect,
  statistic), anchors included as `none` at `0.0`: rejection rate ± MC SE, the
  recorded pooled-eigengap distribution, the `angle` null-width (`q95`)
  dispersion, and the selected-dimensionality distribution.
- `design_point_power.png` — the target statistic's power at the top effect per
  point against the first design axis, one line per value of the second,
  annotated with the median eigengap.
- `design_point_decision.json` / `.csv` — when `acceptance.design_point`
  declares the rule (`trajectory_mode`, `statistic`, `min_power_at_top`,
  `confirmation_se_threshold`, `prefer`): per column `meets` (`rate − k·SE ≥
  floor`), `marginal` (point estimate only), `fails`, or `unavailable`, with
  the anchor's Type I rates beside it, and a verdict — the first `meets` column
  in preference order, or `revise_claim`. Advisory only.

The committed Phase 5 design-point pilot,
`examples/trajectory_power_study/phase5_design_point_pilot.json`, crosses
ρ ∈ {0.0, 0.5, 0.8} with `n_samples` ∈ {300, 600, 1200} at `p_dmp = 0.1`.

Every row of `design_point_operating.csv` also carries `component_selection`
(`cv` / `forced`, or `mixed` — which cannot arise from a well-formed study and
flags a configuration error), read from the records' integration metadata.

### Retained-rank ladder reporting

When the design grid declares `evaluation.integration_params.forced_components`
(see *Nested evaluation axes*), the report additionally writes:

- `rank_ladder.png` — per trajectory mode (the zero-effect anchor first), each
  statistic's rejection rate at the mode's top effect against the retained rank
  with MC error bars; the cross-validated column is drawn at its recorded median
  selected rank with a star marker and a `CV` label so it is never read as a
  forced rung. Built from `design_point_operating.csv` alone; not written for
  grids without the rank axis.
- `rank_decision.json` / `.csv` — when `acceptance.rank_decision` declares the
  rule (`axis`, `target` pair, `protected` pairs, `type_i_bound`,
  `gain_se_multiplier`, `loss_se_multiplier`). With the `null` (CV) column as
  reference, a forced rank *qualifies* when (a) its target power at the top
  effect exceeds the reference by more than the gain multiplier × pooled MC SE,
  (b) its zero-effect anchor is within `alpha + se_tolerance·sqrt(alpha(1−alpha)/n)`
  on every statistic, and (c) no protected power at the top effect falls below
  the reference by more than the loss multiplier × pooled SE. Verdict `keep_cv`
  when no rank qualifies, else `adopt_fixed_rank` naming the **smallest**
  qualifying rank (parsimony predeclared); per column, every criterion's status
  with the rates, SEs, and thresholds it used, and the failing pairs/statistics
  by name. Advisory only: it never feeds the Phase 4 gate or the acceptance
  targets. The JSON states that columns are paired on identical datasets.

The production guard that refuses forced-rank records
(`assert_production_component_selection`) admits a forced record only when its
design point declares the rank axis with the same rank the integration recorded;
any other forced record is still rejected. Forced columns stay invisible to the
power curves, specificity matrix, Type I table, gate, and acceptance targets.

The committed Phase 5 latent-rank ladder,
`examples/trajectory_power_study/phase5_latent_rank_ladder.json`, holds the
chosen design point fixed (ρ = 0, n = 1200, `p_dmp = 0.1`) and varies only the
rank over {`null`, 3, 4, 6, 9, 12} for all four modes.

### Measurement columns and the anchor Type I rule

`evaluation.orientation_estimator` and `evaluation.integration_method` are
ordinary evaluation axes, OFAT or design-grid. Columns that differ only in
evaluation axes share generator parameters, the matched-seed family and the RRPP
seed (`evaluation.seed`), so they measure the same datasets with the same
permutation draws in different spaces; a `null` estimator column writes records
identical to a config without the axis. `concat` and `pca` columns record no
PLS component selection: `design_point_operating.csv` reports their rates,
eigengap and `angle` null width and leaves `component_selection` and the
selected-rank fields empty, and the rank ladder ignores them. Mixed columns that
are not a rectangle use the explicit `design_grid.points` form.

`acceptance.anchor_decision` (`alpha`, `se_multiplier`, `statistics`, `name`)
writes `anchor_decision.json` / `.csv`. Per design point — the baseline column
included — and per named statistic, the zero-effect anchor's rate over its `R`
available replicates is compared with `alpha + se_multiplier·sqrt(alpha(1 − alpha)/R)`
(0.0695 at R = 500, α = 0.05, k = 2): `holds` at or below, `exceeds` above,
`unavailable` without replicates. A point is `size_holds` when every statistic
holds, `excess` when any exceeds, `unavailable` otherwise. Every other cell at
the point (and the `type_i_baseline` cells, under the baseline point) is listed
under `other_controls` with its rates and no classification. Advisory: it never
feeds the Phase 4 gate or the acceptance targets, and a config without it
writes nothing new.

The committed **net-displacement Type I profiles**,
`examples/trajectory_power_study/net_displacement_type_i_phase6.json` and
`net_displacement_type_i_phase5.json`, are null-only (the shared anchor and the
translation control at e = 1.0) across five explicit columns — `pls_cv × pc1`
(baseline), `pls_cv × net_displacement` (the lead arm), PLS forced to rank
2(K − 1), PCA at rank 2(K − 1) and `concat`, each × net displacement — at
500 × 999, 12 cells / 6,000 units per design. Each sits in its parent study's
seed family with the same evaluation parameters, so its baseline column
reproduces the parent's records exactly. After both are reported,
`scripts/net_displacement_follow_up.py --phase6 … --phase5 … --out follow_up_decision.json`
reads the two `anchor_decision.json` files and selects the next change:
`paper_grade_power_comparison` when the lead arm is `size_holds` at both
designs, `diagnose_anchor_excess` otherwise, with a per-design localization
(`every_column` / `some_columns` / `none` of the net-displacement columns).

### Angle-pivotality diagnostic

`scripts/angle_null_pivotality.py` reads a merged record set and reports whether
each replicate's null moved with its own observed statistic:

```bash
# 5 cells x 100 replicates = 500 work units; ~26 minutes on 16 cores.
for i in $(seq 0 15); do
  python scripts/run_study_shard.py \
      --config examples/trajectory_power_study/angle_pivotality_diagnostic.json \
      --out-dir results/angle-pivotality-2026-09-01 \
      --shard-index "$i" --n-shards 16 --error-policy record &
done
wait
python scripts/motco_study.py merge --out-dir results/angle-pivotality-2026-09-01
python scripts/angle_null_pivotality.py \
    --merged  results/angle-pivotality-2026-09-01/merged.jsonl \
    --out-dir results/angle-pivotality-2026-09-01/report
```

Do **not** pass `--n-jobs`: it is part of the cell parameter signature, so
overriding it changes the permutation draws and breaks resume.

It writes `pivotality_association.csv` (correlation and slope of the null's
mean/sd/q95 against the observed statistic, with a Fisher-z interval),
`pivotality_rejection_split.csv` (observed statistic and critical value split by
rejection outcome), `pivotality_standardized.csv` (as-specified vs
cross-replicate standardized rejection rate, controls included), and
`pivotality_spectrum.csv` (per cell and statistic, the Spearman and log–log
Pearson association between the *recorded* pooled eigengap and the width of that
replicate's own null, with Fisher-z intervals; `status = unavailable` for record
sets predating the spectrum field).

The standardized counterfactual is a **diagnostic, not a deployable test** — it
borrows a reference `z` distribution from the null-control cells. Within-replicate
studentization is a no-op: it rescales both sides of the comparison and leaves
the p-value unchanged.

Note that `enumerate_study` adds two cells the profile does not request — a
`none` Type I baseline and a `translation` negative control. They are kept on
purpose: they are the null-control reference the standardized counterfactual
calibrates against.

The 2026-09-01 run of this profile is written up in
`docs/reports/angle-null-pivotality-2026-09-01.md`.

Interpretation notes:

- **Diagonal of the specificity matrix** is power at the largest
  `effect_size` for the matching statistic; off-diagonals are
  specificity (should sit near α).
- **Type I table** reads the `none` cells. Both per-statistic and the
  combined rule should respect α within Monte Carlo uncertainty.
- **Acceptance report** is non-gating — it's a structured record of
  whether pre-specified targets were met given the SE of the
  rejection-rate estimates.

---

## 9. Programmatic API

The CLI scripts are thin wrappers around the public API exposed in
`motco.simulations.study.__init__`. You can drive everything from
Python if you prefer:

```python
from pathlib import Path
from motco.simulations.study import (
    load_study_config, enumerate_study, run_shard,
    summarize_study, summarize_combined_rule,
    build_specificity_matrix, build_power_curves, build_type_i_table,
    evaluate_targets, write_report_csvs, write_target_report,
)

config = load_study_config("examples/trajectory_power_study/smoke.json")
grid   = enumerate_study(config)

# Shard 0 of 4, in-process:
records = run_shard(grid, shard_index=0, n_shards=4,
                    out_dir=Path("/tmp/motco-smoke"),
                    error_policy="record")
```

---

## 10. Troubleshooting

- **Reference cache missing** (`ReferenceCacheMissingError`) — the
  committed `data/intersim_reference.npz` is absent from the install.
  Reinstall the package, or regenerate it once in an R environment with
  InterSIM: `Rscript src/motco/simulations/export_reference.R
  --output-dir <dir>` then `build_cache_from_export(<dir>)`. No R is
  needed for normal runs.
- **Shard wallclock spills past `#SBATCH --time`** — the shard is
  resumable: just resubmit the failed array task ids. Then either
  lower the per-shard load (raise `N_SHARDS`) or raise `--time`.
- **Memory blowups under high `n_jobs`** — RRPP's parallel workers
  each hold a copy of the permuted residuals. Lower
  `--cpus-per-task` and re-probe.
- **Reports disagree with what you expect on a null cell** — check
  that `effect_sizes` starts at `0.0` and that `none` survived
  enumeration; both are required for the Type I view.
- **Re-running an old shard skips everything / produces no new
  records** — the parameter signature includes a seed-derivation
  version. After changes to derivation logic or the generator surface
  (e.g. the move to the numpy generator), shard files produced before
  the change have a stale signature and are re-executed automatically
  on the next `run_shard`. No manual deletion needed.
