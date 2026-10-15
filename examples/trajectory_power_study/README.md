# Trajectory power study

A reproducible, sharded study that characterizes the Type I error and power of the
MOTCO trajectory test (delta, angle, shape) under semi-synthetic InterSIM datasets.

## Workflow

```
study config (YAML/JSON)
    │
    ▼  enumerate_study → SimulationGrid (Type I + power cells)
    │
    ▼  shard runner: scripts/run_study_shard.py
    │     (one cluster array task per shard, writes shard_<i>.jsonl)
    │
    ▼  merge: python scripts/motco_study.py merge --out-dir <dir>
    │     (combines shards into merged.jsonl, dedup by (cell, replicate))
    │
    ▼  report: python scripts/motco_study.py report --config <cfg> --out-dir <dir>
          (per-statistic + combined-rule summaries → specificity matrix,
           power curves, Type I table; CSV + PNG; acceptance-target report)
```

## Local smoke run

```bash
# Generate a few shards locally (no cluster):
python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/smoke.json \
    --out-dir /tmp/motco-smoke \
    --shard-index 0 --n-shards 4 --error-policy record
python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/smoke.json \
    --out-dir /tmp/motco-smoke \
    --shard-index 1 --n-shards 4 --error-policy record
# (repeat for shards 2 and 3, or run with --n-shards 1 to do them all locally)

python scripts/motco_study.py merge --out-dir /tmp/motco-smoke
python scripts/motco_study.py report \
    --config examples/trajectory_power_study/smoke.json \
    --out-dir /tmp/motco-smoke
```

Outputs land under `/tmp/motco-smoke/report/`:

- `specificity_matrix.csv` / `.png` — mode × statistic rejection rates
- `power_curves.csv` / `.png` — per-statistic rejection rate vs effect size
- `type_i_table.csv` / `type_i.png` — null-cell per-statistic + combined-rule rates
- `acceptance_report.csv` / `.json` — pre-specified targets evaluated against
  observed Monte Carlo uncertainty (non-gating)

## SLURM cluster run

```bash
# Submit an array of size N_SHARDS:
sbatch \
    --array=0-63 \
    --export=ALL,STUDY_CONFIG=$(pwd)/examples/trajectory_power_study/smoke.json,STUDY_OUT=$(pwd)/results,N_SHARDS=64 \
    scripts/motco_study_array.sbatch

# After completion:
python scripts/motco_study.py merge  --out-dir results
python scripts/motco_study.py report --config examples/trajectory_power_study/smoke.json --out-dir results
```

Failed array tasks can be resubmitted with `--array=7,12,40` — the shard-resume
guard (parameter signature) skips already-completed replicates.

## Config quick reference

| Field             | Purpose                                                    |
|-------------------|------------------------------------------------------------|
| `intersim`        | Baseline InterSIM params (R generator)                     |
| `generator`       | Baseline semi-synthetic perturbation params                |
| `evaluation`      | Integration method, RRPP permutations, n_jobs              |
| `trajectory_modes`| Power-grid modes (e.g. `magnitude`, `orientation`, …)      |
| `effect_sizes`    | Power-grid effect sizes                                    |
| `axes`            | OFAT axes, namespaced `generator.` / `evaluation.`, or the nested `evaluation.integration_params.<key>` (`null` = key absent; see the ladder section) |
| `design_grid`     | Crossed design points: `{"axes": {...}}`, every axis listing its baseline value (see below) |
| `n_replicates`    | Replicates per cell                                        |
| `base_seed`       | Deterministic seed root                                    |
| `alpha`           | Significance level for rejection rates                     |
| `acceptance`      | Pre-specified Type I, power, and specificity targets       |
| `acceptance.gate` | Phase 4 gate parameters (see below); omit for pre-Phase-4 configs |
| `acceptance.design_point` | Advisory design-point rule over a `design_grid` (see the design-point pilot) |
| `acceptance.rank_decision` | Advisory retained-rank rule over the rank axis (see the latent-rank ladder) |
| `acceptance.anchor_decision` | Advisory Type I rule over every design point's zero-effect anchor (see the net-displacement studies) |
| `acceptance.column_decision` | Advisory candidate-vs-baseline column rule: target gain/no-loss, anchor, mandatory controls, protected pairs (see the net-displacement studies) |
| `attribution`     | Which cells get orientation-attribution diagnostics        |
| `matched_seeds`   | Opt-in matched generator seeds across primary cells        |
| `generator.surgery_censoring` | Pool-limited-surgery policy; leave at the `"error"` default (see below) |
| `generator.group_stage_sizes` | Exact group × stage sample sizes, `[[A per stage], [B per stage]]`; replaces `n_samples`/`stage_sample_prop`/`group_ratio` (which must then be omitted). Absent = proportional sizing, byte-identical to before (see the Phase 6 small-n study) |
| `evaluation.integration_params.layers` | Omic blocks to measure, e.g. `["methylation", "expression"]`; absent = all three, byte-identical to before. Usable as a design-grid axis (`null` = all three) whose columns share data |
| `generator.magnitude_kind` | Magnitude construction (`joint` default, `all`, `extremes`); every committed config names it, and a new config must choose it deliberately (see below) |
| `evaluation.orientation_estimator` | Orientation estimator for `angle` (`net_displacement` default, `pc1`); every committed config pins `"pc1"`, and a new config must choose it deliberately (see below). In a design grid `null` names the default and `"pc1"` names PC1 |
| `report_contract` | Declared reporting/execution rules (`driver_component`, `cross_replicate_driver_agreement: descriptive`, `n_jobs_override: forbid`/`warn`); adds `report_contract.json` + `driver_report.csv` and, under `forbid`, makes the runner refuse `--n-jobs` (see the Phase 5 paper-grade study) |

Unknown top-level keys are rejected by name, so a misspelled block cannot be
silently ignored.

**Which config is the paper-grade one?** `phase5_power_study.json`. The older
`study.json` is **superseded**: it predates every decision since Phase 4
(n = 300, `p_dmp = 0.2`, `n_jobs = -1`, `"surgery_censoring": "clamp"`, no
matched seeds, no gate, and a specificity target on orientation → `shape` that
the predeclared cross-talk fails by construction). It stays as the historical
record and must not be run as Phase 5.

`none` is always present as the Type I baseline (enforced by enumeration);
`translation` is added explicitly as a second negative control.

### `generator.surgery_censoring` — why the pre-Phase-5 configs set `"clamp"`

`orientation`, `translation`, and `shape` with `shape_kind="relocate"` draw
their surgery from a finite pool of CpGs, so a large `group_effect_size` can
request more sites than the pool holds. The generator's default policy,
`"error"`, refuses to realize a partial surgery, and enumeration rejects any
cell whose requested effect exceeds the expected pool headroom before compute
is spent.

**Every config in this directory older than the Phase 5 pilots predates that
policy** and was run under the old silent clamping, so each one carries an explicit
`"surgery_censoring": "clamp"` to stay loadable and enumerable as the record of
what was actually run. At `p_dmp = 0.2` with four stages the axis saturates at
roughly `e ≈ 0.56` for `orientation` and `e ≈ 0.29` for `translation` — above
those, distinct requested effects produced near-identical realized datasets
(see [the geometry audit](../../docs/reports/geometry-audit-2026-09-01.md),
finding F2). Their top cells are therefore **not** independent power
measurements, and `realized_surgery.csv` in the study report flags the affected
pairs.

**A new config must not copy the flag.** Leave `surgery_censoring` at its
default and choose an effect axis that respects the headroom — lower the axis
top, lower `p_dmp`, or use fewer stages. Enumeration reports the saturating
effect for the offending cell, which is the number to design against.

### `generator.magnitude_kind` — why the historical configs pin `"all"`

`magnitude_kind` defaults to `"joint"`, which scales every omic's δ by the same
`1 + e`. It became the default after the
[Phase 5 magnitude re-measurement](../../docs/reports/phase5-magnitude-remeasurement-2026-10-08.md)
showed that it meets the magnitude controls and that the former default, `"all"`
(methylation δ alone), does not. Ten configs in this directory predate the flip and
named no kind, so they ran under `"all"`: `angle_pivotality_diagnostic.json`,
`orientation_signfix_rerun_100x199.json`, `phase4_pilot_100x199.json`,
`phase4_smoke.json`, `phase5_design_point_pilot.json`,
`phase5_latent_rank_ladder.json`, `phase5_power_study.json`, `pilot_50x199.json`,
`smoke.json` and `study.json`. Each now carries an explicit
`"magnitude_kind": "all"` so it still describes what was run. The two magnitude
re-measurement profiles name `"joint"`.

**The pin is behavior-neutral.** `parameter_signature` hashes the resolved
generator parameters, so an explicit `"all"` serializes exactly like the old
implicit default. Every cell id, parameter signature, matched seed and shard
assignment is unchanged, and committed results stay resumable.
`tests/test_magnitude_default_pins.py` checks this against
`tests/data/historical_config_signatures.json`, a fixture captured on the
pre-flip code with the configs unedited.

**Committed provenance hashes.** Pinning changed these files' bytes, so their
sha256 no longer matches the `config_sha256` recorded in committed
`PROVENANCE.txt` files. Those hashes refer to the config as it stood when each
run executed, at the code revision its `PROVENANCE.txt` names; the provenance
files are frozen outputs and were not rewritten. The pin is the only edit.

**A new config must not copy the pin.** Choose `magnitude_kind` deliberately:
`"joint"` is the construction the paper describes as a size change; `"all"`
reproduces only the historical methylation-only construction. The same test
requires every committed config here to name the key, so a config that omits it
fails.

### `evaluation.orientation_estimator` — why every committed config pins `"pc1"`

`orientation_estimator` defaults to `"net_displacement"`: each group's orientation
is the unit vector from its first to its last stage mean. It became the harness
default after the
[net-displacement power comparison](../../docs/reports/net-displacement-power-2026-10-14.md)
selected it (orientation/`angle` power 0.300 → 0.648 at Phase 6 and 0.840 → 0.924
at Phase 5, anchor size held, shape → `angle` leak removed). Every config in this
directory was run under PC1, the former default, so each now carries an explicit
`"orientation_estimator": "pc1"` in its `evaluation` block, and the four
net-displacement profiles name their baseline column `"pc1"` rather than `null`
(a `null` there would now name the new default and silently swap which column is
the primary one).

**The pin is behavior-neutral.** `"pc1"` normalizes to `None`, which is left out
of every hash payload, so every cell id, parameter signature, matched seed and
shard assignment is unchanged, records are byte-identical, and committed results
stay resumable. `tests/test_orientation_default_pins.py` checks all 19 configs
against `tests/data/orientation_default_identity.json`, a fixture captured on the
pre-flip code with the configs unedited, and asserts that each baseline
evaluation resolves to PC1.

**Committed provenance hashes.** As with the `magnitude_kind` pin, these files'
bytes changed, so their sha256 no longer matches the `config_sha256` in committed
`PROVENANCE.txt` files. Those hashes refer to each config as it stood at the run
revision its `PROVENANCE.txt` names; the provenance files are frozen and were not
rewritten.

**A new config must not copy the pin.** Choose `orientation_estimator`
deliberately: `"net_displacement"` (or omitting the key) is the estimator the
paper describes; `"pc1"` reproduces only the historical measurement. The same
test requires every committed config here to name the key and forbids a `null`
estimator in a committed design grid.

## Phase 4 pilot

`phase4_pilot_100x199.json` is the committed Phase 4 medium pilot: pooled PLS on
M-value methylation, `n_samples=300`, four stages, 100 replicates per cell, 199
RRPP permutations, the four established modes, and effects `0.00`–`1.00` by
`0.25`. It supersedes `pilot_50x199.json`, which predates the corrected shape
estimator and realized-geometry diagnostics and is retained as historical
evidence only. `phase4_smoke.json` exercises the identical code paths at
development scale; its numbers carry no scientific meaning.

```bash
# Development-scale smoke over every Phase 4 path (PLS, matched seeds,
# attribution, geometry, gate report):
python scripts/run_study_shard.py \
    --config examples/trajectory_power_study/phase4_smoke.json \
    --out-dir /tmp/motco-phase4-smoke \
    --shard-index 0 --n-shards 1 --error-policy record
python scripts/motco_study.py merge  --out-dir /tmp/motco-phase4-smoke
python scripts/motco_study.py report \
    --config examples/trajectory_power_study/phase4_smoke.json \
    --out-dir /tmp/motco-phase4-smoke

# The medium pilot (1900 work units) on SLURM:
sbatch --array=0-63 \
    --export=ALL,STUDY_CONFIG=$(pwd)/examples/trajectory_power_study/phase4_pilot_100x199.json,STUDY_OUT=$(pwd)/results/phase4,N_SHARDS=64 \
    scripts/motco_study_array.sbatch
python scripts/motco_study.py merge  --out-dir results/phase4
python scripts/motco_study.py report \
    --config examples/trajectory_power_study/phase4_pilot_100x199.json \
    --out-dir results/phase4
```

Do **not** pass `--n-jobs` for a Phase 4 run. RRPP seeds one RNG stream per
worker, so the worker count changes the realized permutation draws, and `n_jobs`
is part of `evaluation_params` and therefore of each cell's parameter signature.
Overriding it makes completed replicates unresumable and the run irreproducible
from the committed config alone. Parallelize across shards instead; the sbatch
script only forwards `--n-jobs` when `STUDY_N_JOBS` is explicitly set.

Write the Phase 4 run to a **new** output directory. Its parameter signatures
include the bumped seed-derivation version and the diagnostic schema versions,
so July shards can neither be resumed into nor overwritten by it.

## Phase 5 design-point pilot

`phase5_design_point_pilot.json` is the committed Phase 5 design-point pilot
(readiness item 4). It keeps the Phase 4 integration and evaluation settings but
changes the baseline to `p_dmp = 0.1` — so the `0.25`–`1.00` effect axis is
realizable without censoring for every pool-limited mode at every design point,
including the ρ = 0 stress-test endpoint — and declares a **crossed design
grid**: baseline continuity ρ ∈ {0.0, 0.5, 0.8} × `n_samples` ∈ {300, 600,
1200}. Every design point gets its own zero-effect anchor plus `orientation`
(the estimand in question) and `translation` (the negative control that binds
the surgery headroom) at effects `0.25`, `0.50`, `1.00`; 100 replicates and 199
permutations per cell; 6,500 work units, with the `n = 1200` columns
dominating compute. `magnitude` and `shape` are left to the Phase 5 study at
the chosen design point — both reached power 1.00 at n = 300 in Phase 4.

**Run 2026-09-08** — see the
[findings report](../../docs/reports/phase5-design-point-pilot-2026-09-08.md) and
the committed outputs under `results/phase5-design-point-2026-09-08/`
(`report/` and `PROVENANCE.txt`). All 6,500 units completed with zero failures
and no censored surgery on a SLURM array of 100 single-CPU shards (AMD EPYC
7662): 70 core-hours of recorded unit runtime, 70 minutes wall. Verdict:
**chosen**, ρ = 0 at n = 1200 (orientation `angle` power 0.88, 1·SE lower
bound 0.85).

Measured per-unit cost with BLAS pinned to one thread (PLS double CV + 199
RRPP permutations): 12.5 s / 18.5 s / 27 s per unit at n = 300 / 600 / 1200 on
a desktop core (Ryzen 9 7900X), 28 s / — / 56 s on an EPYC 7662 core — about
35 desktop-core-hours or 70 EPYC-core-hours for the whole pilot. Running many
shards on one workstation costs more per unit than the isolated figures
(memory bandwidth: 25 s / 69 s / 124 s with 22 concurrent shards on 12 cores),
so a 22-shard local run takes roughly 6 hours. **Always pin BLAS**
(`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`): with the
default thread pool each shard spawns ~32 threads, 12 parallel shards drove
the load average past 100, and the 2026-09-04 rehearsal measured under that
oversubscription over-estimated the cost by about 7×.

```bash
# As executed on 2026-09-08 (partition/resource flags are cluster-specific).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
sbatch -p 512x1024 --cpus-per-task=1 --mem=2G --time=6:00:00 --array=0-99 \
    --export=ALL,STUDY_CONFIG=$(pwd)/examples/trajectory_power_study/phase5_design_point_pilot.json,STUDY_OUT=$(pwd)/results/phase5-design-point-2026-09-08,N_SHARDS=100 \
    scripts/motco_study_array.sbatch
python scripts/motco_study.py merge  --out-dir results/phase5-design-point-2026-09-08
python scripts/motco_study.py report \
    --config examples/trajectory_power_study/phase5_design_point_pilot.json \
    --out-dir results/phase5-design-point-2026-09-08
```

The report adds, beside the usual outputs:

- `design_point_operating.csv` — one row per (design point, mode, effect,
  statistic): rejection rate ± MC SE, the recorded pooled-eigengap distribution,
  the `angle` null-width (`q95`) dispersion, and the selected-dimensionality
  distribution. The baseline column and each point's anchor (`none` at `0.0`)
  are included.
- `design_point_power.png` — orientation `angle` power at the top effect vs
  `n_samples`, one line per ρ, annotated with the median eigengap.
- `continuity_resolved_orientation.csv` — resolved on ρ **and** every other
  design coordinate, so rows never pool across `n_samples`.
- `design_point_decision.json` / `.csv` — the predeclared rule
  (`acceptance.design_point`) evaluated per column: `meets` when
  `rate − k·SE ≥ floor`, `marginal` when only the point estimate clears it,
  `fails` otherwise. The verdict names the first `meets` column in the declared
  preference order (`n_samples` ascending, then ρ ascending) or is
  `revise_claim`. It is advisory: it never feeds the Phase 4 gate or the
  acceptance targets.

## Phase 5 latent-rank ladder

`phase5_latent_rank_ladder.json` is the committed Phase 5 latent-rank ladder
(readiness item 3). It derives from the design-point pilot and holds the chosen
design point fixed — ρ = 0, `n_samples = 1200`, four stages, `p_dmp = 0.1`,
default fail-loud censoring — and varies **only the retained PLS rank** through
a one-axis design grid over the nested evaluation parameter
`evaluation.integration_params.forced_components` ∈ {`null`, 3, 4, 6, 9, 12}.
`null` means "key absent": that column runs the production stage-supervised
double CV and is the reference; every other column fits the pooled PLS model at
that fixed rank (`component_selection = "forced"`). All four modes
(`magnitude`, `orientation`, `shape`, `translation`) at effects
`0.25`/`0.50`/`1.00` plus one zero-effect anchor per column; 100 replicates and
199 permutations; 6 × 13 = 78 power cells plus the two baseline Type I
controls, **8,000 work units**.

**Same data, different measurement.** The rank axis is an evaluation-namespace
axis, so every column shares the primary matched-seed family and identical
generator parameters with the baseline: at every replicate index the six
columns evaluate the *same* generated dataset at different ranks. Differences
between columns are measurement differences, not sampling differences (the
duplicate-dataset guard keys on evaluation identity and accepts this by
design). Forced 3 is included so CV-selection variability (recorded range 2–4,
median 3 = `n_stages − 1`) is separated from the rank itself, and because
"fixed `n_stages − 1`" is itself a candidate group-blind rule. 12 is the top
rung because the 2026-09-03 latent-rank probe saw no change from 9 to 12.

**Predeclared decision** (`acceptance.rank_decision`, advisory): against the
`null` column, a forced rank qualifies when it raises orientation `angle`
power at the top effect by more than 2 pooled MC SEs, keeps every statistic's
anchor within `alpha + 2·SE`, and lowers neither magnitude `delta` nor shape
`shape` power at the top effect by more than 2 pooled SEs. Verdict `keep_cv`,
or `adopt_fixed_rank` with the **smallest** qualifying rank. Only group-blind
rules are candidates — the geometry audit's caution stands: group-aware sizing
would void the fixed-latent-space RRPP conditioning.

**Run 2026-09-08** — see the
[findings report](../../docs/reports/latent-rank-ladder-2026-09-08.md) and the
committed outputs under `results/phase5-latent-rank-2026-09-08/` (`report/` and
`PROVENANCE.txt`). All 8,000 units completed with zero failures and no censored
surgery on a SLURM array of 100 single-CPU shards (20 min wall; CV units median
57 s, forced units 0.2 s; 24.5 core-hours). Verdict: **`keep_cv`** — no fixed
rank raised orientation `angle` power (0.79–0.86 vs 0.85 at CV) and every rank
above 3 lost shape `shape` power (1.00 → 0.57). Phase 5 keeps stage-supervised
double CV, which selects rank 3 = `n_stages − 1` there.

Cost: the CV column measured 55.8 s per unit on an EPYC 7662 core in the
design-point pilot (n = 1200, 199 permutations); forced columns skip the double
CV and are cheaper. Budget ~120 EPYC core-hours; 100 single-CPU shards finish
in roughly 1.5 h wall. **Pin BLAS** and **do not pass `--n-jobs`** (it enters the
parameter signature and changes the permutation draws).

```bash
# SLURM (partition/resource flags are cluster-specific).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RUN=results/phase5-latent-rank-$(date -u +%F)
sbatch -p 512x1024 --cpus-per-task=1 --mem=2G --time=6:00:00 --array=0-99 \
    --export=ALL,STUDY_CONFIG=$(pwd)/examples/trajectory_power_study/phase5_latent_rank_ladder.json,STUDY_OUT=$(pwd)/$RUN,N_SHARDS=100 \
    scripts/motco_study_array.sbatch
python scripts/motco_study.py merge  --out-dir $RUN
python scripts/motco_study.py report \
    --config examples/trajectory_power_study/phase5_latent_rank_ladder.json \
    --out-dir $RUN

# Locally, in K resumable shards (each shard is one process; pin BLAS first):
for i in $(seq 0 $((K-1))); do
  python scripts/run_study_shard.py \
      --config examples/trajectory_power_study/phase5_latent_rank_ladder.json \
      --out-dir $RUN --shard-index $i --n-shards $K --error-policy record &
done; wait
```

The report adds, beside the design-point outputs above:

- `design_point_operating.csv` gains `component_selection` (`cv`/`forced`) per
  row; `median_selected_lv` equals the forced rank on forced rows.
- `rank_ladder.png` — per mode (anchor first), each statistic's rate at the top
  effect vs rank with MC error bars; the CV column sits at its median selected
  rank with a star marker and a `CV` label.
- `rank_decision.json` / `.csv` — the verdict, the chosen rank, and every
  column's per-criterion status with the rates, SEs, and thresholds used.

### `design_grid` — crossed design points

`axes` varies one factor at a time off the baseline; `design_grid.axes` is
**crossed**. Every combination of the declared values is a *design point*, and
each non-baseline point enumerates the full power grid (one zero-effect anchor
when `0.0` is among the effect sizes, plus every mode × nonzero effect) with
the point's coordinates applied to the baseline generator/evaluation
parameters. Rules:

- every axis MUST list the baseline value — the baseline point is served by the
  primary cells (which are stamped with their coordinates) and never
  re-emitted;
- an axis is either crossed or OFAT, never both;
- design cells carry phase `power_design`, `design_point = {axis: value}`, and
  `varied_axis = "design_grid"`, so the primary power curves, specificity
  matrix, Type I table, gate, and acceptance targets read the baseline column
  only, exactly as they would without the grid;
- design cells join the primary matched-seed family, so columns are paired at
  the same replicate index; the duplicate-dataset guard and the surgery-headroom
  check apply to every design cell using its own parameters.

Configs without `design_grid` enumerate and report exactly as before.

### Matched seeds

`matched_seeds.enabled` makes every **primary** power cell draw its generator
seed from one shared `(seed family, replicate index)` key, so at a given
replicate index every mode and effect starts from the same generated reference
and requested-effect comparisons are paired. Negative-control, Type I, and OFAT
cells keep their own seed families and stay independent draws; design-grid
cells join the primary family (see above). Persistence keys remain
`(cell_id, replicate_index)`.

At `group_effect_size = 0` the generator returns group B's baseline unchanged
and consumes no extra randomness, so per-mode zero-effect cells inside one
family would be byte-identical datasets. With
`matched_seeds.shared_zero_effect_anchor`, enumeration therefore emits **one**
mode-agnostic zero-effect primary cell and every mode's power curve resolves its
`0.00` point from that shared anchor; the report flags those rows with
`from_shared_anchor` so the modes' nulls are not read as independent evidence.
Enumeration also asserts that no two primary cells would generate identical data
at the same replicate index.

Configs without `matched_seeds` keep the pre-Phase-4 per-cell seed derivation
unchanged.

### Attribution diagnostics

`attribution.enabled` turns on bounded orientation-attribution diagnostics for
the cells the selector names — in the Phase 4 pilot, every nonzero primary
`orientation` cell, with 100 frozen-model bootstrap replicates and `top_k=20`.
Eligibility is resolved **during enumeration**, so it enters the cell's
parameter signature and never depends on an observed p-value; conditioning on
significance would bias the reported stability.

Attribution requires `evaluation.integration_method: "pls"` — it conditions on
the fitted PLS estimator, which `concat` and `snf` do not produce. The exact
estimator and standardized joint matrix that produced the trajectory scores are
reused; no second fit or component selection happens.

Each eligible replicate persists a compact, versioned record: effective
settings, ordered transitions with observed / PLS-captured / residual path
lengths and retention, the top-k signed features per transition and component
(standardized **and** original units, with the unit basis labeled so M-value
methylation is not read as beta values), bootstrap sign and selection stability,
and precision/recall against generator truth. Fitted estimators, full
standardized matrices, bootstrap matrices, and unrestricted feature tables are
never persisted. A cell that was never selected records
`attribution_status: "not_requested"`; an eligible replicate whose attribution
fails records `"failed"` with a reason and still contributes its trajectory
measurement.

Generator truth is defined as the features whose group-stage differential *mean
change* differs between groups, including CpG→gene→protein propagated effects —
so a real downstream driver is not scored as a false positive.

### Phase 4 gate

`acceptance.gate` carries every threshold the gate applies; the study code
hard-codes none of its own, so a gate is re-specified by editing the config.

| Gate field | Meaning |
|---|---|
| `alpha` | Significance level for control checks |
| `control_se_tolerance` | `k` in the one-sided bound `alpha + k·sqrt(alpha(1-alpha)/n)` |
| `monotonicity_se_tolerance` | A downward power step is tolerated when it is at most this many combined MC SEs |
| `min_power_at_top` | Default power floor at the largest effect |
| `confirmation_se_threshold` | An exceedance smaller than this many MC SEs counts as *marginal* |
| `max_marginal_exceedances` | How many marginal control exceedances may be tolerated (default 1) |
| `control_modes` | Modes whose cells are Type I controls at every effect level |
| `rules` | `mandatory_power` / `mandatory_control` / `descriptive` mode-statistic pairs |
| `require_complete_records` | Whether incomplete records/diagnostics block `proceed` |

The mandatory rules are:

- **Type I inflation** — the `none` baseline and *every* `translation` effect
  level (translation is a location-only offset at any effect) must keep each
  available statistic at or below the one-sided bound.
- **Power** — each `mandatory_power` pair must reach its floor at the top effect
  and be non-decreasing within `monotonicity_se_tolerance` combined MC SEs.
- **Control** — each `mandatory_control` off-diagonal pair (magnitude's `angle`
  and `shape`) is checked against the same inflation bound at its largest effect.
- **Completeness** — every expected work unit resolved, every completed PLS
  record carrying selected-component and realized-geometry metadata, and every
  eligible orientation record carrying valid attribution diagnostics or a
  recorded failure.

`descriptive` pairs — orientation's and shape's off-diagonals — are reported
against realized geometry but never gate: Phase 2 established that both
constructions are genuinely mixed after biological propagation and joint
preprocessing, so demanding purity there would fail a correct estimator.

Gate multiplicity is predeclared. Each control cell contributes one one-sided
test per statistic, so at a true rate of exactly alpha each exceeds its bound
with probability around `0.023`. One marginal exceedance therefore must not
decide the phase: when exactly one control statistic exceeds its bound by less
than `confirmation_se_threshold` MC SEs and nothing else mandatory fails, the
report emits `indeterminate` and names the confirmation re-run. Two or more
exceedances, or any exceedance of at least one SE, is `hold`. Study execution
never aborts because a scientific gate fails.

### Phase 4 outputs

Alongside the existing specificity, Type I, and power outputs, `report` writes:

| Output | Contents |
|---|---|
| `phase4_operating.csv` | Rejection rate and MC SE by mode, effect, and statistic, with `from_shared_anchor` |
| `phase4_geometry.csv` / `.png` | Realized geometry by mode, effect, checkpoint, scope, statistic, and path length, with a `measurement_space` label and unavailable counts |
| `phase4_pls_selection.csv` / `phase4_selected_components.png` | Selected component counts, effective CV settings, AUROC, and missing-diagnostic counts |
| `phase4_attribution.csv` / `phase4_attribution_stability.png` | Availability, observed-versus-captured retention, cross-replicate top-k Jaccard and sign agreement, bootstrap stability, and truth recovery |
| `phase4_localization.csv` | First checkpoint at which each off-diagonal response becomes material |
| `phase4_gate.csv`, `phase4_gate_decision.json` | Every gate observation, and the `proceed` / `hold` / `indeterminate` decision |

Localization compares each checkpoint only against its **own** zero-effect null,
on a scale-free quantity (`delta` divided by that checkpoint's path length;
`angle` and `shape` already dimensionless). Raw distances are never compared
across the standardized feature space and the PLS latent space. Its labels —
construction-present, sampling/preprocessing-associated, projection-associated —
describe *where* a response first appears, not what caused it, and gate nothing.

## Phase 5 paper-grade study

`phase5_power_study.json` is the committed paper-grade Phase 5 configuration —
the study the paper's operating-characteristic claims are read from. **It ran
2026-09-10** (job 880599 on `ing`, 100 shards, 3 h 04 m wall, 182.7 core-hours,
9,500/9,500 units, 0 failures, 0 censored surgeries); the Phase 4 gate returned
**HOLD** on the two magnitude mandatory controls, with all three mandatory power
diagonals and all 21 Type I checks met. Results are in
`results/phase5-2026-09-10/` and the findings report is
[`docs/reports/phase5-paper-grade-2026-09-10.md`](../../docs/reports/phase5-paper-grade-2026-09-10.md).
It
supersedes `study.json` (see the config quick reference) and derives from the
latent-rank ladder's cross-validated column (`metadata.derives_from`): identical
`generator` and `evaluation.integration_params`, no `design_grid` (the ladder
returned `keep_cv`, so the retained rank is the stage-supervised double-CV
choice), and paper-grade Monte Carlo precision.

| | |
|---|---|
| design point | ρ = 0, n = 1200 (300 per group-stage cell), four stages, `p_dmp = 0.1`, default fail-loud `surgery_censoring` |
| measurement | pooled PLS on M-value methylation, double-CV rank (`cv1_splits` 3, `cv2_splits` 4, 5 repeats, ≤ 20 components) |
| grid | modes magnitude / orientation / shape / translation × effects 0 / 0.25 / 0.50 / 0.75 / 1.00; every point under headroom (orientation saturates ≈ 1.69, translation ≈ 2.00) |
| Monte Carlo | 500 replicates × 999 permutations; `n_jobs` 1 |
| cells / units | 2 Type I controls + 1 shared zero-effect anchor + 16 power cells = 19 cells, 9,500 units |
| seeds | matched seeds, family `phase5-primary` (independent of both pilots), one shared zero-effect anchor |
| attribution | nonzero orientation primary cells; 100 bootstraps, top-20, seed 0 |
| gate | Phase 4 roles — mandatory power magnitude/`delta`, orientation/`angle`, shape/`shape` (floor 0.80); mandatory control magnitude/`angle`, magnitude/`shape`; descriptive orientation/`delta`, orientation/`shape`, shape/`delta`, shape/`angle`; controls `none`, translation |
| acceptance | `type_i` 0.05 ± 2 SE; `power` 0.80 on the three diagonals; `specificity` = exactly the gate's mandatory controls (translation × 3, magnitude/`angle`, magnitude/`shape`) — orientation → `shape` is predeclared cross-talk, not a target |
| report contract | `driver_component: observed`, `cross_replicate_driver_agreement: descriptive`, `n_jobs_override: forbid` |

### Report contract

The `report_contract` block turns the readiness-item-5 rules into enforced or
echoed outputs (`docs/phase5-readiness.md` §5):

- `report/driver_report.csv` — the declared component only, one row per
  (mode, effect, transition): precision/recall vs generator truth, mean
  selected count, **within-replicate** bootstrap stability, replicate
  accounting. No `top_k_jaccard` / `sign_agreement` columns; those stay in
  `phase4_attribution.csv` (which keeps all three components) as descriptive
  quantities, never a stability claim.
- `phase4_attribution_stability.png` — retitled "Within-replicate bootstrap
  stability (observed component)" and drawn from the bootstrap series only.
- `report/report_contract.json` — the resolved contract plus the uniform
  `n_jobs` the records carry (a mixed set is refused) and the shared
  zero-effect anchor's `cell_id`, `resolves_modes`, and `counted_as: 1`.
- `scripts/run_study_shard.py` exits 2 before enumeration when `--n-jobs`
  (`STUDY_N_JOBS`) differs from `evaluation.n_jobs`.

The dated findings report must follow
[`phase5_report_template.md`](phase5_report_template.md), which carries a
section for every contract item, the `PROVENANCE.txt` field list, the two
design-point hand-offs (eigengap and `angle` null width beside every
orientation result; the orientation surgery's non-ρ-invariant realized
contrast), and reproduction commands without `--n-jobs`. The gate's outputs
keep their `phase4_*` names — they name the mechanism, not the phase.

### Cost

A one-replicate rehearsal (2026-09-09, 19 units as 19 concurrent
single-threaded processes on a 24-core workstation, BLAS pinned) measured a
median **43 s per unit** (min 30, max 51) for the CV fit + 999 permutations +
attribution at n = 1200; attribution units cost the same as the others
(median 42.9 s vs 43.0 s). Budget ≈ 120 workstation core-hours, ≈ 150 EPYC
7662 core-hours by the ladder's per-core ratio; 100 single-CPU shards finish
in roughly 1.5–2 h wall.

### Running it

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RUN=results/phase5-$(date -u +%F)

# SLURM (partition/resource flags are cluster-specific). Do NOT set STUDY_N_JOBS:
# the contract forbids it and every array task would exit 2.
sbatch -p <partition> --cpus-per-task=1 --mem=2G --time=6:00:00 --array=0-99 \
    --export=ALL,STUDY_CONFIG=$(pwd)/examples/trajectory_power_study/phase5_power_study.json,STUDY_OUT=$(pwd)/$RUN,N_SHARDS=100 \
    scripts/motco_study_array.sbatch

# Locally, in K resumable shards (one process each):
for i in $(seq 0 $((K-1))); do
  uv run python scripts/run_study_shard.py \
      --config examples/trajectory_power_study/phase5_power_study.json \
      --out-dir $RUN --shard-index "$i" --n-shards K --error-policy record &
done; wait

uv run python scripts/motco_study.py merge  --out-dir $RUN
uv run python scripts/motco_study.py report \
    --config examples/trajectory_power_study/phase5_power_study.json --out-dir $RUN
```

Commit `report/` and a hand-written `PROVENANCE.txt` (fields listed in the
template); the shard and merged JSONL stay gitignored.

## Phase 5 magnitude re-measurement (`joint` construction)

The Phase 5 paper-grade run held the gate on its two magnitude mandatory
controls (`magnitude`→`angle` 0.276, `magnitude`→`shape` 0.742 at e = 1.00
against a bound near 0.07). The
[magnitude-construction diagnostic](../../docs/reports/magnitude-construction-diagnostic-2026-09-10.md)
traced that to the construction, not the statistics: `magnitude_kind='all'`
scales methylation's δ alone, which is size-pure *within* every omic block and
impure only as the geometry of concatenating one block that grew with two that
did not. `magnitude_kind='joint'` scales all three per-omic δ by the same
`1 + e`, so group B's native-space trajectory is exactly `1 + e` times group A's
in every block — size-pure in the joint standardized space as well (population
`angle` ~1e-6, `shape` ~1e-15, the anchor's floor). The re-measurement runs the
magnitude mode under `joint` on a three-rung ladder; `'all'` remains the default
and no historical config changes.

**Rung 1 — analytic bracket** (`results/magnitude-axis-bracket-2026-10-08/`).
`scripts/magnitude_construction_diagnostic.py --bracket` evaluates the
population-standardized joint geometry of both constructions over e = 0…1
(steps of 0.0025) at the Phase 5 design point (minutes on a workstation, no
sampling, no cluster, no R) and writes `effect_axis_bracket.csv`. It relates the
two effect axes — joint's realized size is 2.47–2.58 × the `all` construction's
at the same nominal e — and gave the pilot its first grid: the `all`
construction's `delta` power is already 1.000 at its first Phase 5 point
(e = 0.25, realized size 5.95), and joint realizes 3.23 at e = 0.05 and 6.34 at
e = 0.10, so the bracket placed the rise between those. The bracket is
population geometry and cannot see the RRPP null width; the pilot can.

**Rung 2 — pilot** `phase5_magnitude_pilot.json`: 50 × 199, gate disabled.
Submission 1 ran 0 / 0.02 / 0.05 / 0.10 / 0.25 / 1.00 and measured `delta`
power 1.000 from e = 0.02 (realized PLS-latent joint `delta` ≈ 68 × e against a
`delta` null 95th percentile near 0.40, so the rise sits near e ≈ 0.006);
submission 2 extended the grid with 0.0025 / 0.005 / 0.01, resuming the
completed units in the same run directory. 11 cells (2 Type I baselines + 1
shared zero-effect anchor + 8 magnitude cells), 550 units, ≈ 12 core-hours. Its
note (`results/phase5-magnitude-pilot-2026-10-08/NOTES.md`) records the grid
chosen for paper grade and why.

**Rung 3 — paper grade** `phase5_magnitude_remeasurement.json`: 500 × 999,
`n_jobs = 1` with overrides forbidden, the Phase 5 `report_contract`, and the
Phase 4 gate reduced to exactly the three magnitude rules (`magnitude`/`delta`
mandatory power; `magnitude`/`angle` and `magnitude`/`shape` mandatory control;
`control_modes = ["none"]`). Effect grid **0 / 0.0025 / 0.005 / 0.01 / 0.02 /
0.25 / 1.00** on joint's native axis (e is the per-omic size ratio minus one,
exactly; it is *not* rescaled to the `all` curve — the realized joint `delta`
recorded per replicate relates the two axes): the three smallest points bracket
the measured rise, 0.02 is its saturation, 0.25 is the `all` construction's
first Phase 5 point, and 1.00 is the control stress point (2.47 × the production
construction's realized size). 9 cells, 4,500 units; the 2026-10-08 run took
106.7 core-hours and 2 h 06 m wall on 100 single-CPU shards (node contention;
per-unit median 88.7 s against Phase 5's 61.5 s).

**Run 2026-10-08** — see the
[addendum](../../docs/reports/phase5-magnitude-remeasurement-2026-10-08.md) and
the committed outputs under `results/phase5-magnitude-2026-10-08/` (`report/`
and `PROVENANCE.txt`). 4,500 / 4,500 units, 0 failures. Reduced gate
**PROCEED**: magnitude/`delta` 1.000 at the top with the rise resolved at
0.130 / 0.384 / 0.836 / 0.996 over e = 0.0025–0.02; magnitude/`angle` 0.024 and
magnitude/`shape` 0.000 at e = 1.00 against the 0.0695 bound (Phase 5's `all`
construction: 0.276 / 0.742); the anchor reproduces the Phase 5 anchor at all
500 replicate indices including every p-value (`report/anchor_reproduction.csv`,
`scripts/anchor_reproduction.py`). The withheld magnitude-specificity claim is
lifted for `joint`; `all` remains the default until a separate change flips it.

Both profiles copy `generator` (except `magnitude_kind`),
`evaluation.integration_params`, `base_seed` 600, and the matched-seed family
`phase5-primary` from `phase5_power_study.json`, so the shared zero-effect
anchor and every group A baseline are **byte-identical to the Phase 5 records
at the same replicate index** (the anchor's reproduction of Phase 5's anchor is
stated in the addendum as a reproducibility check). Orientation, shape, and
translation cells are not re-run: matched seeds make them identical to the
committed Phase 5 records and their findings stand.

```bash
# Bracket (workstation).
python scripts/magnitude_construction_diagnostic.py --bracket \
    --out-dir results/magnitude-axis-bracket-$(date -u +%F)

# Pilot, then paper grade, on SLURM (partition/resource flags are cluster-specific).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RUN=results/phase5-magnitude-pilot-$(date -u +%F)    # or results/phase5-magnitude-$(date -u +%F)
CFG=examples/trajectory_power_study/phase5_magnitude_pilot.json   # or phase5_magnitude_remeasurement.json
sbatch -p 512x1024 --cpus-per-task=1 --mem=2G --time=6:00:00 --array=0-99 \
    --export=ALL,STUDY_CONFIG=$(pwd)/$CFG,STUDY_OUT=$(pwd)/$RUN,N_SHARDS=100 \
    scripts/motco_study_array.sbatch
python scripts/motco_study.py merge  --out-dir $RUN
python scripts/motco_study.py report --config $CFG --out-dir $RUN
```

The findings go in a dated addendum,
`docs/reports/phase5-magnitude-remeasurement-<date>.md`, following
`phase5_report_template.md` section for section with the non-magnitude sections
marked not applicable; it lifts the withheld magnitude-specificity claim or
reports it as failed on the construction's own terms.

## Phase 6 small-n operating study (SEA-AD design)

Phase 6 applies MOTCO to SEA-AD MTG astrocytes: donor pseudobulks of RNA and
ATAC, Sex as the group, three merged ADNC stages (Not AD + Low / Intermediate /
High), 80 donors at ≥30 nuclei (F 11/10/28, M 9/10/12). Every operating
characteristic measured so far used three blocks, four balanced stages and
n = 300–1200, so this study measures `delta`/`angle`/`shape` Type I and power
at the cohort's own design before any SEA-AD result is interpreted. The
simulation is a structural analogue, not an RNA/ATAC emulator.

Two keys make the design expressible; both are byte-identical when absent:

- `generator.group_stage_sizes` — the exact donor table (group A = F, group
  B = M, so the surgery lands on the smaller group: the conservative choice).
- `evaluation.integration_params.layers = ["methylation", "expression"]` — all
  three InterSIM blocks are generated and two are measured. Methylation stays
  because every surgery acts on methylation indicators; it is the regulatory
  analogue of ATAC.

**Pilot** `phase6_small_n_pilot.json`: baseline column only, 100 × 199, gate
disabled, shared log-spaced effect axis 0 / 0.005 / 0.01 / 0.02 / 0.05 / 0.10 /
0.25 / 0.50 / 1.00 for all four modes, `magnitude_kind: joint`, new matched-seed
family `phase6-small-n` (base seed 800). 35 cells (2 Type I baselines + 1 shared
zero-effect anchor + 4 × 8), 3,500 units. A local smoke at 9 permutations ran
at about 15 s per unit on a workstation; the pilot records the cluster median.
It brackets where each mode's target statistic rises at n = 80.

**Pilot run 2026-10-09** (`results/phase6-small-n-pilot-2026-10-09/`, `NOTES.md`
and `PROVENANCE.txt`): 3,500 / 3,500 units, 0 failures, median 8.0 s per unit
(8.6 core-hours). Type I is controlled. Magnitude/`delta` rises over e =
0.02–0.25 (0.09 → 0.93) and saturates at 0.50. Orientation/`angle` and
shape/`shape` reach only 0.24 and 0.18 at e = 1.00, the relocation
construction's ceiling (the fraction is clamped at 1). At three stages the
stage-mean eigengap is small (anchor median 0.11), so the `angle` null is very
wide, and the stage-supervised PLS space attenuates the rotation.

**Paper grade**, split per the pilot (spec split rule), both at 500 × 999 with
the Phase 5 report contract and the Phase 4 gate rules as advisory targets:

- `phase6_small_n_magnitude.json` — magnitude, effects 0 / 0.02 / 0.05 / 0.10 /
  0.25 / 0.50 / 1.00; reduced magnitude gate. 30 cells, 15,000 units.
- `phase6_small_n_study.json` — orientation, shape, translation, effects 0 /
  0.25 / 0.50 / 0.75 / 1.00; orientation attribution. 54 cells, 27,000 units.

Both copy the pilot's generator, evaluation (permutations raised to 999), base
seed and matched-seed family, so all three profiles share their zero-effect
anchor datasets. Each crosses `generator.group_stage_sizes` (≥30 nuclei, ≥50
nuclei `[[10, 9, 25], [9, 9, 12]]`, n = 74) with
`evaluation.integration_params.layers` (two blocks, `null` = three) into four
columns. The design grid always crosses its axes, so the ≥50-nuclei
three-block column runs too, although neither contrast needs it. The findings
report reads cohort size as baseline vs ≥50 nuclei at two blocks and block count
as baseline vs three blocks at n = 80, and ends with a per-statistic
interpretability statement for the case study. Planning cost: about 10 s per
unit, about 117 core-hours, 140 with a 20% margin.

**Run 2026-10-09** — see the
[findings report](../../docs/reports/phase6-small-n-2026-10-09.md) and
`results/phase6-small-n-magnitude-2026-10-09/` and
`results/phase6-small-n-2026-10-09/` (`report/`, `PROVENANCE.txt`).
42,000 / 42,000 units, 0 failures, 150 recorded core-hours. Node n4 ran out of
memory mid-run; its 12 stalled shards were cancelled and resumed elsewhere
without loss. Type I is controlled (0.010–0.040).

- **Magnitude** (gate PROCEED): `delta` 0.206 / 0.542 / 0.956 at e = 0.05 /
  0.10 / 0.25, with `angle` and `shape` at the floor.
- **Orientation/shape** (gate HOLD): `angle` 0.300 and `shape` 0.228 at the
  construction maximum. Both constructions are detected more often by other
  statistics.
- **Third block** (same datasets): raises power for every mode.
- **≥50-nuclei cohort**: leaves magnitude unchanged. It lowers orientation at
  e = 1.00 through CV rank selection.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RUN=results/phase6-small-n-pilot-$(date -u +%F)
CFG=examples/trajectory_power_study/phase6_small_n_pilot.json
sbatch -p 512x1024 --cpus-per-task=1 --mem=2G --time=6:00:00 --array=0-99 \
    --export=ALL,STUDY_CONFIG=$(pwd)/$CFG,STUDY_OUT=$(pwd)/$RUN,N_SHARDS=100 \
    scripts/motco_study_array.sbatch
python scripts/motco_study.py merge  --out-dir $RUN
python scripts/motco_study.py report --config $CFG --out-dir $RUN
```

## Net-displacement estimator studies

Two pairs of committed profiles compare the `net_displacement` orientation
estimator (opt-in when they ran, the harness default since 2026-10-15) with PC1 at the Phase 6 and Phase 5 designs. Each profile
is an explicit `design_grid.points` grid over `evaluation.integration_method` ×
`evaluation.orientation_estimator` in its parent study's seed family, so the
baseline column re-computes the parent's records exactly and every column
measures the same datasets with the same permutation draws.

- **Type I** (`net_displacement_type_i_phase6.json`,
  `net_displacement_type_i_phase5.json`): null-only, anchor + translation at
  e = 1.0, five columns, 500 × 999, judged by `acceptance.anchor_decision`;
  `scripts/net_displacement_follow_up.py` picks the next change. Run 2026-10-13:
  size holds in every column at both designs; verdict
  `paper_grade_power_comparison`
  ([report](../../docs/reports/net-displacement-type-i-2026-10-13.md)).
- **Power** (`net_displacement_power_phase6.json`,
  `net_displacement_power_phase5.json`): `pls_cv × pc1` against
  `pls_cv × net_displacement`, four modes at e ∈ {0.25, 0.5, 0.75, 1.0},
  500 × 999, 18,000 units per design, judged by `acceptance.column_decision`
  (`gain` at Phase 6, `no_loss` at Phase 5);
  `scripts/net_displacement_power_follow_up.py` combines the two verdicts and
  `scripts/net_displacement_power_checks.py` verifies the shared statistics and
  the parent reproduction. Run 2026-10-14: the candidate qualifies at both
  designs (orientation/`angle` at e = 1.0: 0.300 → 0.648 at Phase 6,
  0.840 → 0.924 at Phase 5; shape → `angle` leak 0.138 → 0.054 and
  0.796 → 0.022); verdict `adopt_net_displacement`
  ([report](../../docs/reports/net-displacement-power-2026-10-14.md)).

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
RUN=results/net-displacement-power-$(date -u +%F)
for d in phase6 phase5; do
  CFG=examples/trajectory_power_study/net_displacement_power_$d.json
  sbatch -p 512x1024 --cpus-per-task=1 --mem=4G --time=12:00:00 --array=0-99 \
      --export=ALL,STUDY_CONFIG=$(pwd)/$CFG,STUDY_OUT=$(pwd)/$RUN/$d,N_SHARDS=100 \
      scripts/motco_study_array.sbatch
done
# After completion:
for d in phase6 phase5; do
  python scripts/motco_study.py merge  --out-dir $RUN/$d
  python scripts/motco_study.py report \
      --config examples/trajectory_power_study/net_displacement_power_$d.json --out-dir $RUN/$d
done
# Then the follow-up and checks; see the study README, "The column decision rule".
```

Do not pass `--n-jobs`; parallelize across shards. Planning cost on `ing`:
about 13 s per unit at Phase 6 (64 core-hours) and 62 s at Phase 5
(319 core-hours).
