# Simulations

InterSIM bridge utilities and semi-synthetic trajectory generators for multi-omics simulation studies.

## R dependency

The bridge is optional and requires `Rscript` plus the R `InterSIM` package:

```r
install.packages(
  "InterSIM",
  repos = c("https://cran.r-universe.dev", "https://cloud.r-project.org")
)
```

Check availability before running a simulation:

```python
from motco.simulations import check_intersim_available

availability = check_intersim_available()
if not availability.available:
    print(availability.message)
```

## Example

```python
from motco.simulations import InterSIMParams, run_intersim

result = run_intersim(
    InterSIMParams(
        seed=1203,
        n_sample=100,
        cluster_sample_prop=(0.3, 0.3, 0.4),
        delta_methyl=1.0,
        delta_expr=1.0,
        delta_protein=1.0,
        p_dmp=0.1,
    )
)

methylation = result.methylation
expression = result.expression
proteomics = result.proteomics
clusters = result.clusters
```

## Semi-synthetic trajectory generation

The semi-synthetic trajectory generator samples directly from the numpy-native
InterSIM reimplementation (no R at runtime) and returns a MOTCO-ready dataset with:

- aligned methylation, gene expression, and proteomics matrices
- sample metadata containing `sample_id`, `group`, and `stage`
- truth metadata recording trajectory mode, per-stage/group differential
  indicators, per-omic δ, and the generator seed

Group A is a random baseline trajectory; group B is a deterministic transform of
group A's per-stage differential indicators selected by `trajectory_mode`.

```python
from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    generate_semisynthetic_trajectory,
)

dataset = generate_semisynthetic_trajectory(
    SemiSyntheticTrajectoryParams(
        seed=99,
        trajectory_mode="magnitude",
        n_samples=120,
        n_stages=3,
        group_effect_size=0.2,
        group_ratio=0.5,
    ),
)

sample_metadata = dataset.metadata
truth = dataset.truth
```

Supported trajectory modes (all governed by the unified `group_effect_size` knob,
where `0` is the null for every mode):

| Mode | Injected group-specific pattern |
|------|---------------------------------|
| `none` | Identical to baseline; useful for Type I error scenarios |
| `translation` | Constant observed-space location offset |
| `magnitude` | Scales δ (size); see the `magnitude_kind` table below |
| `orientation` | One global per-omic feature permutation (rotation) |
| `shape` | Permutes interior stages only (bend); requires at least three stages |

`magnitude_kind` selects how `magnitude` scales group B (the `motco simulate --magnitude-kind`
flag has the same default):

| `magnitude_kind` | Construction |
|------------------|--------------|
| `'joint'` (default) | Scales every omic's δ (methylation, expression, proteomics) by the same `1 + e`. Size-pure in the joint standardized space, not only within each block. This is the construction the Phase 5 magnitude re-measurement validated. |
| `'all'` | Scales the methylation δ alone at every stage. Size-pure within the methylation block but rotates the concatenated trajectory. |
| `'extremes'` | Scales only the endpoint stages' methylation indicators. |

`'extremes'` and the shape mode's `shape_kind='magnitude'` remain methylation-only scalings; the
joint default does not extend to them.

**Explicit group × stage sizes.** `group_stage_sizes` replaces proportional sizing
(`n_samples` × `stage_sample_prop` × `group_ratio`) with an exact table: one tuple of per-stage
counts per group, rows ordered as `group_labels` and columns as stages. Every group-stage cell gets
exactly the stated number of samples, and truth metadata records the realized table under
`truth["group_stage_sizes"]`. Use it when the imbalance varies by stage, as in a real cohort:

```python
SemiSyntheticTrajectoryParams(
    seed=1,
    n_stages=3,
    group_stage_sizes=((11, 10, 28), (9, 10, 12)),  # 80 samples, the SEA-AD MTG design
)
```

The table must have two rows of `n_stages` integers, each at least 1. Setting it together with a
non-default `n_samples`, `stage_sample_prop` or `group_ratio` is an error that names the
conflicting settings. When the table is absent (the default `None`), generation, truth metadata,
cell ids and parameter signatures are byte-identical to the generator before the field existed:
hash payloads omit the field while it is `None`.

**Migration.** Before `default-joint-magnitude-construction` the default was `'all'`. Code that
relied on the default and wants the old construction must pass `magnitude_kind="all"` (or
`--magnitude-kind all`). The committed historical study configs under
`examples/trajectory_power_study/` already pin `"all"`, so their cells, signatures and datasets are
unchanged.

### Stage-span construction

The modes above are **feature-set surgeries** (the default, `construction=None`; `"surgery"`
normalizes to it). `construction="stage_span"` is an opt-in, measurement-defined alternative: it
edits group B's exact null population stage means inside the span of B's own stage configuration,
so the planted difference lies where the stage-supervised PLS measurement looks. Use it to ask how
sensitive the test is to a difference the measurement can see; the surgeries remain the
biologically grounded constructions.

```python
SemiSyntheticTrajectoryParams(
    seed=1,
    trajectory_mode="orientation",
    group_effect_size=0.5,          # a 45° rotation
    construction="stage_span",
    span_layers=("methylation", "expression"),  # must equal the measured blocks
)
```

The edit is made over `span_layers` (default: all three blocks) in analysed units (methylation in
M-values), each feature divided by its null population pooled SD: the within-cell variance (the
reference covariance diagonal) plus the between-cell variance of the null means at the realized
cell proportions. With `z_s` B's standardized stage means, centroid `c`, path length `L`, mean step
`L̄ = L / (S − 1)` and `u` the unit net displacement:

| Mode | Stage-span edit | Intended population geometry |
|------|-----------------|------------------------------|
| `magnitude` | `c + (1 + e)(z_s − c)` | path length ×(1 + e); `angle` = `shape` = 0 |
| `orientation` | rigid rotation about `c` by θ = 90°·e in the plane of `u` and the next principal direction `v` (e ≤ 2) | net-displacement `angle` = θ; `delta` = `shape` = 0 |
| `shape` | the stage the `shape` surgery perturbs moves e·L̄ off the chord between its neighbours, then the path is rescaled to `L` and recentred on `c` | `shape` > 0; `delta` = `angle` = 0 |
| `translation` | every stage moves e·L̄ along `u` | all three statistics 0 (the in-span negative control) |

The plant is mapped back to analysed units and added to group B's cell means before sampling
(`generate_omics(..., cell_shift_*=)`); blocks outside `span_layers` are untouched. The construction
draws no random numbers, so a stage-span dataset equals the `none` dataset at the same seed plus the
plant: group A and B's sampling noise are shared with the surgery datasets of the same seed family.
`orientation` and `shape` need a second in-span direction, so fewer than three stages, or a
configuration whose second singular value is at most `1e-8` times its first, is an error;
`span_layers` without `stage_span`, and unknown, duplicated or empty layers, are rejected.

Truth metadata records `construction`, `span_layers`, the effect parameters under `transform`
(`scale`, `rotation_degrees`, `shape_stage`/`displacement`/`rescale`, or `offset_length`) and a
`stage_span` block with the configuration's `path_length` and `mean_step`, the
`intended_geometry` (net-displacement `delta`, `angle`, `shape` and path lengths in the
null-standardized coordinates) and the `plant_norm` per stage. The population trajectories include
the plant. Stage-span cells draw from no feature pool, so they have no surgery headroom and no
realized-surgery summary. While both settings are unset every dataset, truth record, cell id and
signature is byte-identical to the generator before they existed.

## Evaluation harness

The evaluation harness runs one `SemiSyntheticTrajectoryDataset` through MOTCO integration and trajectory testing. It is the per-replicate layer used before larger Type I error or power grids.

```python
from motco.simulations import (
    SimulationEvaluationParams,
    evaluate_semisynthetic_trajectory,
)

evaluation = evaluate_semisynthetic_trajectory(
    dataset,
    SimulationEvaluationParams(
        integration_method="concat",
        permutations=0,
    ),
)

observed_delta = evaluation.pair_statistics["delta"]
truth = evaluation.truth_metadata
```

Supported integration methods:

| Method | Behavior |
|--------|----------|
| `concat` | Column-binds methylation, expression, and proteomics matrices after deterministic per-feature standardization by default |
| `snf` | Builds per-omic affinity matrices, fuses them with SNF, and uses spectral embedding as the trajectory outcome matrix |

**Block selection.** `integration_params["layers"]` names the omic blocks the evaluation measures,
for every integration method: e.g. `{"layers": ["methylation", "expression"]}`. All three blocks are
still generated; only the selected ones enter pooled preprocessing, integration, every
realized-geometry checkpoint (per-block scopes for the selected blocks only, joint scopes over them)
and orientation attribution. The selection is normalized to canonical order (methylation,
expression, proteomics) and must be non-empty, duplicate-free and known; integration metadata
records it as `integration_params.layers`. When the key is absent, all three blocks are measured and
the evaluation result is byte-identical to the harness before the key existed (no `layers` metadata
entry is added).

**Orientation estimator.** `orientation_estimator` selects how each group's orientation is measured
for the observed `angle`, every RRPP draw and every realized-geometry checkpoint `angle`. The default
is `"net_displacement"` (the unit vector from the first to the last stage mean), adopted after the
[net-displacement power comparison](../reports/net-displacement-power-2026-10-14.md); it is part of the
cell id and parameter signature and is recorded in `runtime_metadata["orientation_estimator"]`. Pass
`orientation_estimator="pc1"` to opt out to the principal-axis estimator: `"pc1"` normalizes to
`None`, which is left out of every hash payload and record, so a PC1 evaluation is byte-identical to
the harness before the selection existed. Every committed study config and historical diagnostic pins
`"pc1"`. Only `angle` depends on the selection. The stats primitives (`estimate_difference`, `RRPP`)
keep `"pc1"` as their default.

**Rank floor.** `min_components` is the smallest rank the stage-supervised double CV may select for
`pls` integration: the candidates are `min_components … max_components − 1`, so every outer-fold
selection, the modal rank and the pooled fit are at or above it. The default is `2`, so `angle` and
`shape` are measured in a space of at least two dimensions (at rank 1 two orientations are parallel or
antiparallel by construction and `angle`'s null is a two-point distribution). The default is hashed
into the cell id and parameter signature and recorded in `runtime_metadata["min_components"]` and in
the latent metadata's `integration_params.min_components`. Pass `min_components=1` for the unfloored
rule: `1` normalizes to `None`, which is left out of every hash payload and record, so a floor-1
evaluation is byte-identical to the harness before the floor existed. Every committed study config
pins `"min_components": 1` beside `"pc1"`, as do the historical diagnostics. `max_components` is
clamped to `[min_components + 1, min(n_features, n_samples)]`; a sample that cannot satisfy the floor
raises an error naming the floor and the counts. The diagnostic `forced_components` override keeps its
own range (`[2, max]`) and ignores the floor; `pca` and `pls_within_group` have fixed ranks by
construction. In a study config `evaluation.min_components` is a field and a design axis: `null`
names the evaluator's default and `1` the unfloored rule, as distinct values. The floor guarantees two
dimensions, not a second axis with signal; it was validated by `scripts/rank_floor_check.py` (the 2026-10-23 run
returned `adopt_floor` at both calibrated strengths, [report](../reports/rank-floor-2026-10-23.md)).

Set `permutations=0` for observed statistics only. When `permutations > 0`, the harness runs RRPP and computes upper-tail empirical p-values with plus-one correction:

```text
p = (1 + count(null >= observed)) / (1 + n_permutations)
```

The result includes observed `delta`, `angle`, and `shape` matrices, scalar two-group pair statistics, optional p-values, latent matrix metadata, generator truth metadata, runtime metadata, group/stage levels, and the trajectory contrast. Shape pair statistics and p-values are reported as unavailable for datasets with fewer than three stages.

### Realized-geometry decomposition

Generated semi-synthetic datasets also expose exact group-stage population means in integration units: methylation M-values and native expression/proteomics units. Evaluation uses one pooled fitted preprocessing artifact for both diagnostics and PLS, then reports:

| Checkpoint | Scopes | Interpretation |
|---|---|---|
| `population_native` | each omic | Exact construction and biological-cascade geometry without sampling noise |
| `population_standardized` | each omic and joint | Exact means under the replicate's pooled fitted feature scaling |
| `observed_standardized` | each omic and joint | Geometry actually presented to integration, including finite sampling |
| `pls_latent` | joint | Geometry retained in the fitted production PLS representation |
| `population_latent` | joint | Exact population means seen through the replicate's fitted PLS basis (only with `span_capture`) |

Every applicable scope records both group path lengths, pairwise `delta`, `angle`, and `shape`, and explicit availability flags. Joint native-space geometry is intentionally omitted because the unstandardized omic blocks use incompatible units. Raw distances across feature and latent checkpoints are not scale-equivalent; use the decomposition to locate changes in behavior rather than subtracting distances across spaces.

**Span capture.** `SimulationEvaluationParams(span_capture=True)` (opt-in; `False` normalizes to the
hash-omitted `None`) asks a `pls` evaluation how much of the population trajectories its fitted basis
sees. It reuses the evaluation's own fitted model and preprocessor — nothing is refitted, so the latent
matrix, rank, statistics and p-values are identical with and without it — and adds the
`population_latent` checkpoint above plus a `span_capture` block. The block is computed in the model's
input coordinates with an orthonormal (QR) basis of its rotations:

| Field | Meaning |
|---|---|
| `group_difference` | Share of Σ_s‖d_s‖² inside the span, d_s the population B − A difference at stage s minus its mean over stages; `null` when the difference is the same at every stage (the anchor, a translation), which none of the three statistics can see |
| `pooled_stage` | Share of the centred average stage means inside the span |
| `rank`, `n_features`, `random_share` | The fitted rank, the measured feature count and `rank / n_features`, the expected share of a random direction |

The shares are energy (squared-norm) shares. For a single direction the plasmode check's
`latent_capture.py` reported the norm share ‖Qᵀd‖/‖d‖, the square root of this one. Other integration
methods, and datasets without population trajectories, record neither the checkpoint nor the block.

**Stage-span guard.** A dataset generated with `construction="stage_span"` may only be measured on its
span layers: the harness raises, before any preprocessing or fit, when `integration_params["layers"]`
names other blocks, and the error names both block sets.

See the [Phase 2 findings](../reports/realized-generator-geometry-2026-08-05.md) for the revised orientation and shape construction contracts and the matched-seed characterization command.

## Grid orchestration

The grid orchestration layer enumerates parameter cells, runs local replicates through the evaluation harness, persists one JSONL row per replicate, resumes completed work, and summarizes rejection rates for Type I error or power studies.

```python
from pathlib import Path

from motco.simulations import (
    SemiSyntheticTrajectoryParams,
    SimulationEvaluationParams,
    SimulationRunConfig,
    enumerate_type_i_grid,
    run_simulation_grid,
    summarize_rejection_rates,
)

grid = enumerate_type_i_grid(
    baseline_generator_params=SemiSyntheticTrajectoryParams(seed=2, n_samples=60),
    evaluation_params=SimulationEvaluationParams(integration_method="concat", permutations=99),
    axes={
        "generator.n_samples": [60, 120],
        "generator.group_ratio": [0.5, 0.7],
    },
    n_replicates=3,
    base_seed=2026,
)

records = run_simulation_grid(
    grid,
    config=SimulationRunConfig(output_path=Path("simulation-results.jsonl")),
)
summaries = summarize_rejection_rates(records, alpha=0.05)
```

Each `SimulationCell` stores a stable `cell_id`, phase, `SemiSyntheticTrajectoryParams`, `SimulationEvaluationParams`, replicate count, base seed, and metadata such as the varied axis. Axis names use explicit namespaces: `generator.<field>` or `evaluation.<field>`.

Initial persistence is JSON Lines. Each row records cell and replicate IDs, deterministic seeds, a parameter signature, status, p-values, pair statistics, realized-geometry diagnostics, truth metadata, runtime metadata, cell metadata, and optional error details. The diagnostic schema version participates in the parameter signature. With `resume=True`, completed rows with matching parameter signatures are skipped. A matching cell/replicate with a different parameter signature raises unless `overwrite=True`; legacy rows remain readable but have no realized-geometry payload.

`summarize_rejection_rates` groups completed replicate rows by cell and statistic, then reports available replicate count, rejection count, rejection rate, Monte Carlo standard error, and unavailable replicate count. Missing statistics, such as shape p-values for two-stage datasets, remain unavailable rather than being counted as non-significant.

## API

::: motco.simulations.InterSIMParams

::: motco.simulations.InterSIMResult

::: motco.simulations.InterSIMAvailability

::: motco.simulations.SemiSyntheticTrajectoryParams

::: motco.simulations.SemiSyntheticTrajectoryDataset

::: motco.simulations.SimulationEvaluationParams

::: motco.simulations.SimulationEvaluationResult

::: motco.simulations.SimulationTrajectoryDesign

::: motco.simulations.LatentIntegrationResult

::: motco.simulations.SimulationCell

::: motco.simulations.SimulationGrid

::: motco.simulations.SimulationReplicateResult

::: motco.simulations.SimulationRunConfig

::: motco.simulations.SimulationSummaryResult

::: motco.simulations.check_intersim_available

::: motco.simulations.run_intersim

::: motco.simulations.generate_semisynthetic_trajectory

::: motco.simulations.integrate_semisynthetic_dataset

::: motco.simulations.build_simulation_trajectory_design

::: motco.simulations.evaluate_semisynthetic_trajectory

::: motco.simulations.enumerate_type_i_grid

::: motco.simulations.enumerate_power_grid

::: motco.simulations.run_simulation_replicate

::: motco.simulations.run_simulation_grid

::: motco.simulations.read_replicate_results

::: motco.simulations.append_replicate_results

::: motco.simulations.summarize_rejection_rates

## Phase 4 diagnostics and study gate

Phase 4 adds two optional, additive layers on top of the existing harness: a
bounded orientation-attribution diagnostic computed from the *same* fitted PLS
estimator that produced the trajectory scores, and a predeclared gate that turns
the study's operating characteristics into an explicit `proceed` / `hold` /
`indeterminate` decision for the paper-grade study.

Both are off by default. `SimulationEvaluationParams.attribution` is a disabled
`AttributionDiagnosticSettings`, and a study configuration without an
`attribution` or `matched_seeds` block behaves exactly as before. Enabling
attribution requires `integration_method="pls"`; `concat` and `snf` produce no
fitted estimator to condition on and are rejected with a descriptive error.

Every gate threshold lives in the study configuration's `acceptance.gate` block,
never in the summary or report code. See
`examples/trajectory_power_study/README.md` for the config fields, the matched-
seed and shared-anchor semantics, the gate rules, and the full output inventory.

::: motco.simulations.AttributionDiagnosticSettings

::: motco.simulations.attribution_diagnostics.compute_attribution_diagnostics

::: motco.simulations.attribution_diagnostics.derive_truth_driver_features

::: motco.simulations.study.AttributionSelector

::: motco.simulations.study.MatchedSeedPolicy

::: motco.simulations.study.Phase4GateConfig

::: motco.simulations.study.GateRule

::: motco.simulations.study.enumerate_study

::: motco.simulations.study.summarize_realized_geometry

## Angle-pivotality diagnostic

::: motco.simulations.pivotality.association_table

::: motco.simulations.pivotality.rejection_split_table

::: motco.simulations.pivotality.standardized_counterfactual_table

::: motco.simulations.pivotality.replicate_z

::: motco.simulations.pivotality.write_pivotality_tables

::: motco.simulations.study.summarize_pls_selection

::: motco.simulations.study.summarize_attribution

::: motco.simulations.study.localize_off_diagonal

::: motco.simulations.study.type_i_inflation_bound

::: motco.simulations.study.evaluate_phase4_gate

::: motco.simulations.study.Phase4GateDecision

::: motco.simulations.study.build_phase4_frames

::: motco.simulations.study.write_phase4_report

::: motco.simulations.study.render_phase4_figures
