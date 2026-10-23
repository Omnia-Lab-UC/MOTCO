# trajectory-power-study Specification

## Purpose
Define a declaratively-configured, cluster-executable study that characterizes the Type I error and power of the MOTCO trajectory test via per-statistic operating characteristics, with negative controls, pre-specified acceptance targets, and paper-ready reporting (specificity matrix, Type I tables, power curves).

## Requirements

### Requirement: Study is defined by a declarative configuration

MOTCO SHALL provide a declarative, file-based study configuration that fully determines an enumerated grid of simulation cells, so that a study is reproducible from its config alone.

The configuration MUST capture baseline InterSIM/generator/evaluation parameters, the set of trajectory modes, the set of effect sizes, the one-factor-at-a-time axes with their values, an optional design grid, the per-cell replicate count, the base seed, and the pre-specified acceptance targets. A design grid is declared in exactly one of two forms: a **crossed** form naming axes with their values, where every combination of values is a design point; or an **explicit** form listing the design points directly, where every listed point names the same set of axes and assigns each a value. In both forms the baseline design point MUST be present (in the crossed form every axis lists its baseline value; in the explicit form the baseline point is one of the listed points), every axis MUST use the same `generator.` / `evaluation.` namespace rules as the one-factor-at-a-time axes, and no axis may also be a one-factor-at-a-time axis. Declaring both forms at once, listing a point that omits or adds an axis, or listing the same point twice MUST be rejected at load time with an error naming the offending point or axis. Both kinds of axis MAY name a nested evaluation integration parameter (`evaluation.integration_params.<key>`); for such an axis the baseline value is the mapping's current value or `null` when the key is absent, `null` is a legal axis value meaning "key absent", and the baseline-value rule applies to that resolved value.

#### Scenario: Config enumerates a deterministic grid

- **WHEN** a study configuration is loaded and enumerated
- **THEN** it produces a `SimulationGrid` of Type I and power cells with stable, deterministic cell identifiers
- **AND** enumerating the same configuration twice yields identical cell identifiers and parameter signatures

#### Scenario: Config records negative-control modes

- **WHEN** a study configuration is enumerated
- **THEN** the resulting Type I cells include both the `none` (no group effect) and `translation` (location-only group effect) trajectory modes as negative controls

#### Scenario: Invalid configuration is rejected

- **WHEN** a study configuration omits a required field or specifies an unknown trajectory mode, a negative replicate count, or an unknown axis namespace
- **THEN** loading the configuration raises a clear validation error identifying the offending field

#### Scenario: Design grid without its baseline column is rejected

- **WHEN** a study configuration declares a design grid whose values for some axis omit that axis's baseline value, or names an axis with an unknown namespace, or names an axis also declared as a one-factor-at-a-time axis
- **THEN** loading the configuration raises a clear validation error identifying the offending axis

#### Scenario: Configuration without a design grid is unchanged

- **WHEN** a study configuration declares no design grid
- **THEN** enumeration emits exactly the cells it emitted before the design grid existed, with identical identifiers and parameter signatures

#### Scenario: Nested evaluation axis resolves its baseline from the integration parameters

- **WHEN** a study configuration declares the design-grid axis `evaluation.integration_params.forced_components` with values `[null, 4, 6]` and its baseline evaluation integration parameters do not set that key
- **THEN** the configuration loads, the axis's baseline value is `null`, and the `null` value survives a configuration dump and reload unchanged

#### Scenario: Nested evaluation axis without its baseline value is rejected

- **WHEN** a study configuration declares `evaluation.integration_params.forced_components` with values `[4, 6]` while the baseline does not set the key
- **THEN** loading the configuration raises a clear validation error naming the axis and the missing baseline value `null`

#### Scenario: Explicit design points are accepted

- **WHEN** a study configuration declares a design grid as an explicit list of points over the axes `evaluation.integration_method`, `evaluation.integration_params.forced_components`, `evaluation.integration_params.n_components` and `evaluation.orientation_estimator`, one of which is the baseline point
- **THEN** the configuration loads, its design points are exactly the listed points in listed order, and the list survives a configuration dump and reload unchanged

#### Scenario: Explicit design points without the baseline point are rejected

- **WHEN** a study configuration declares an explicit list of design points none of which equals the baseline point
- **THEN** loading the configuration raises a clear validation error naming the missing baseline point

#### Scenario: Malformed explicit design points are rejected

- **WHEN** a study configuration declares both the crossed and the explicit form, or lists a point that omits an axis another point names, or lists two points with identical coordinates
- **THEN** loading the configuration raises a clear validation error naming the offending point or axis

### Requirement: Study executes as resumable per-shard work units

MOTCO SHALL execute an enumerated study as `(cell, replicate)` work units partitioned across a fixed number of shards, so that the study can run as parallel cluster tasks without coordination.

Each shard MUST persist its own JSONL output file and MUST be independently resumable using the existing parameter-signature guard, so re-running a shard skips already-completed replicates and never appends duplicates.

#### Scenario: Work is partitioned deterministically across shards

- **WHEN** a study is run with `n_shards = N` and a given shard index `i`
- **THEN** the shard executes exactly the `(cell, replicate)` units assigned to index `i` by a deterministic partition of all units
- **AND** the union of units across all `N` shards equals the full set of units with no overlaps

#### Scenario: Shard resumes without duplicating completed work

- **WHEN** a shard is re-run and its JSONL output already contains completed records with matching parameter signatures
- **THEN** the shard skips those replicates and appends only missing or failed ones

#### Scenario: Shard records failures without aborting the study

- **WHEN** a replicate within a shard fails and the configured error policy is to record
- **THEN** the shard writes a failed replicate record with error details and continues with remaining units

### Requirement: Shards merge into a single deduplicated result set

MOTCO SHALL merge per-shard JSONL outputs into a single result set, deduplicating by `(cell_id, replicate_index)` and validating parameter-signature consistency across shards.

#### Scenario: Merge combines all shards

- **WHEN** per-shard JSONL files are merged
- **THEN** the merged result contains exactly one record per `(cell_id, replicate_index)` across all shards

#### Scenario: Merge detects inconsistent shards

- **WHEN** two shards contain the same `(cell_id, replicate_index)` with different parameter signatures
- **THEN** the merge raises a clear validation error rather than silently choosing one record

### Requirement: Study characterizes each statistic independently

MOTCO SHALL characterize the trajectory test using per-statistic operating characteristics, reporting the rejection rate of each of the `delta`, `angle`, and `shape` statistics as its own marginal quantity, without multiplicity correction across statistics.

#### Scenario: Per-statistic rejection rates are reported per cell

- **WHEN** merged results are summarized for a cell
- **THEN** the summary reports, for each of `delta`, `angle`, and `shape`, the completed replicate count, rejection count, rejection rate, and Monte Carlo standard error at the configured alpha

#### Scenario: Unavailable statistic is not counted as a rejection

- **WHEN** a statistic is unavailable for a cell (for example `shape` with fewer than three stages)
- **THEN** the summary reports it as unavailable rather than treating it as non-significant

### Requirement: Study reports a combined-rule Type I result

MOTCO SHALL report, as a secondary result, the Type I error rate of the combined decision rule that rejects when any of the three statistics is significant at the configured alpha, computed only over null cells.

#### Scenario: Combined-rule false-positive rate is computed on null cells

- **WHEN** a null cell (negative control) is summarized under the combined rule
- **THEN** a replicate counts as a rejection if any available statistic's p-value is below alpha
- **AND** the reported rate is the fraction of such replicates with its Monte Carlo standard error

### Requirement: Study produces paper-ready reports

MOTCO SHALL produce, from the merged and summarized results, a mode × statistic specificity matrix, Type I tables, and power-curve data, written as CSV and as figures.

#### Scenario: Specificity matrix is produced

- **WHEN** reporting runs on summarized results
- **THEN** it produces a matrix indexed by trajectory mode and statistic whose entries are rejection rates with Monte Carlo standard errors

#### Scenario: Power curves are produced

- **WHEN** reporting runs on power cells
- **THEN** it produces, for each trajectory mode and statistic, rejection rate as a function of effect size, suitable for plotting as a curve with error bars

#### Scenario: Type I table is produced

- **WHEN** reporting runs on null cells across the configured axes
- **THEN** it produces a table of per-statistic and combined-rule rejection rates with Monte Carlo standard errors

### Requirement: Study evaluates results against pre-specified targets

MOTCO SHALL evaluate summarized results against the acceptance targets declared in the configuration and report, per target, whether it is met given Monte Carlo uncertainty.

A Type I control target or a specificity target MUST declare its bound form as one of two values, defaulting to the two-sided form when the key is absent. Under the **two-sided** form the target is met when the empirical rejection rate lies within `k` empirical Monte Carlo standard errors of alpha on either side. Under the **upper** form the target is met when the empirical rejection rate is at or below `alpha + k · sqrt(alpha · (1 − alpha) / R)`, with `R` the number of available replicates — the same bound the Phase 4 gate applies to its control and Type I inflation rules — so that a control rejecting less often than alpha is met. `k` is the target's standard-error tolerance in both forms. The report MUST record the bound form, the alpha-based standard error and the signed deviation from alpha only for a target declared with the upper form; a target in the default form MUST produce exactly the observations and rationale it produced before the bound form existed, so that every committed acceptance report regenerates byte-identically. Power monotonicity targets have no bound form.

#### Scenario: Type I control target is evaluated

- **WHEN** a null cell is evaluated against a Type I control target at alpha
- **THEN** the report indicates whether the empirical rejection rate is within the target's tolerance (expressed in Monte Carlo standard errors) of alpha

#### Scenario: Power monotonicity target is evaluated

- **WHEN** a power mode's diagonal statistic is evaluated against a monotonicity target
- **THEN** the report indicates whether the rejection rate is non-decreasing in effect size and reaches the target's minimum power at the largest effect size

#### Scenario: Specificity target is evaluated

- **WHEN** an off-diagonal mode × statistic combination is evaluated against a specificity target
- **THEN** the report indicates whether its rejection rate stays within tolerance of alpha

#### Scenario: Conservative control is met under the upper bound

- **WHEN** a specificity target declared with the upper bound at alpha `0.05` and tolerance `2` evaluates a control cell whose rejection rate is `0.000` over `500` available replicates
- **THEN** the target is met, the recorded bound is `0.0695`, the recorded deviation is `−0.050`, and the recorded bound form is upper

#### Scenario: Zero-rejection control is not met under the two-sided bound

- **WHEN** the same control cell is evaluated against a specificity target in the default two-sided form
- **THEN** the target is not met, because the empirical standard error of a zero rate is zero and the tolerance collapses to zero

#### Scenario: Inflated control is not met under the upper bound

- **WHEN** a Type I control target declared with the upper bound at alpha `0.05` and tolerance `2` evaluates a null cell whose rejection rate is `0.080` over `500` available replicates
- **THEN** the target is not met and the rationale states the rate exceeds the alpha-based bound `0.0695`

#### Scenario: Default-form targets regenerate byte-identically

- **WHEN** a study whose acceptance targets declare no bound form is reported before and after the bound form exists
- **THEN** its acceptance report CSV and JSON are byte-identical, with no added column or key

#### Scenario: Unknown bound form is rejected

- **WHEN** a study configuration declares a Type I control or specificity target with a bound form other than the two-sided or upper value
- **THEN** loading the configuration raises a clear validation error naming the target entry and the offending value

### Requirement: Study runs on the numpy generator without an R runtime dependency
The trajectory power study SHALL generate every replicate through the numpy generator and cached reference data, so that study execution (including cluster shards) requires no `Rscript` or R `InterSIM` package.

#### Scenario: Shards run without R
- **WHEN** a study shard executes its `(cell, replicate)` units
- **THEN** each replicate is generated from the numpy generator and cached reference data, with no R invocation

#### Scenario: Negative-control modes are retained under the new semantics
- **WHEN** a study configuration is enumerated
- **THEN** the Type I cells still include `none` (no group effect) and `translation` (location-only group effect) as negative controls, now defined by the feature-surgery generator

### Requirement: Acceptance targets are re-specified for the new mode semantics
Because the trajectory modes are redefined, the study's pre-specified acceptance targets SHALL be reset to reflect the operating characteristics of the feature-surgery modes, and prior results SHALL be treated as superseded.

#### Scenario: Acceptance targets reference the new modes
- **WHEN** a study configuration's acceptance targets are evaluated
- **THEN** the per-statistic power and specificity targets correspond to the feature-surgery `magnitude`/`orientation`/`shape` modes

#### Scenario: Specificity demonstration is supported by indicator truth
- **WHEN** a replicate is summarized
- **THEN** the per-stage/group differential indicators emitted by the generator are available to confirm that the injected mode predominantly moves its matching statistic

### Requirement: Study provides a fixed Phase 4 medium-pilot profile

The study SHALL provide a version-controlled Phase 4 configuration using the numpy generator, pooled PLS integration with M-value methylation, 300 samples, four stages, 100 replicates per cell, 199 RRPP permutations, trajectory modes `magnitude`, `orientation`, `shape`, and `translation`, and effect sizes `0.00`, `0.25`, `0.50`, `0.75`, and `1.00`. The configuration MUST explicitly record PLS cross-validation parameters, seed policy, diagnostic settings, study intent, and an acceptance block holding every gate parameter: alpha, tolerance multipliers, minimum power at the top effect, confirmation-rule thresholds, and which mode/statistic pairs are mandatory versus descriptive.

#### Scenario: Phase 4 config is loaded

- **WHEN** the committed Phase 4 configuration is loaded
- **THEN** it deterministically enumerates the expected Type I controls and primary power cells with the declared sample size, stage count, replicate count, permutations, modes, and effects
- **AND** it requires no R runtime dependency

#### Scenario: Primary cells use matched seeds

- **WHEN** primary cells are enumerated for different modes or nonzero effects at the same replicate index
- **THEN** their generated datasets use the same matched replicate seed under a versioned seed-pairing policy
- **AND** the policy is included in parameter signatures so incompatible legacy shards cannot resume into the Phase 4 run

#### Scenario: Zero-effect cells are not duplicated across modes

- **WHEN** the Phase 4 configuration is enumerated
- **THEN** exactly one mode-agnostic zero-effect primary cell is emitted, inside the matched-seed family, and every trajectory mode's power curve resolves its `0.00` point from that shared anchor
- **AND** no two enumerated primary cells generate identical datasets at the same replicate index

#### Scenario: Gate parameters come from the configuration

- **WHEN** Phase 4 gate evaluation runs
- **THEN** every threshold it applies is read from the configuration's acceptance block
- **AND** no mandatory gate threshold is hard-coded in summary or report code

#### Scenario: Orientation diagnostics are selected without significance conditioning

- **WHEN** a primary orientation cell has a nonzero requested effect
- **THEN** its evaluation enables 100 frozen-model attribution bootstrap replicates with `top_k=20`
- **AND** eligibility does not depend on the observed global orientation p-value

### Requirement: Replicate persistence includes Phase 4 diagnostics

The study SHALL persist JSON-safe PLS integration metadata and optional attribution diagnostics beside existing p-values, observed statistics, realized geometry, truth metadata, and runtime metadata. Readers MUST remain able to load legacy records that lack the additive fields, while resume validation MUST reject records created with a different diagnostic configuration or schema version.

#### Scenario: Phase 4 replicate round-trips

- **WHEN** a completed Phase 4 replicate is serialized and loaded
- **THEN** its selected PLS component count, effective cross-validation metadata, realized-geometry diagnostics, and applicable attribution diagnostics are preserved

#### Scenario: Ineligible cell is persisted

- **WHEN** a cell is not selected for attribution diagnostics
- **THEN** its record explicitly identifies attribution as not requested rather than failed or unavailable

#### Scenario: Legacy record is loaded

- **WHEN** a reader loads a record written before Phase 4 diagnostic fields existed
- **THEN** the record loads with empty diagnostic fields
- **AND** it cannot satisfy resume for a Phase 4 cell with a different parameter signature

### Requirement: Phase 4 reporting relates operating characteristics to realized geometry

The study SHALL report per-statistic rejection rates and Monte Carlo uncertainty together with summaries of
every applicable realized-geometry checkpoint. It MUST interpret off-diagonal rejection in light of the
corresponding pre-integration geometry and MUST NOT label a response as estimator cross-talk solely because
the requested mode name differs from the responding statistic.

Materiality SHALL be determined in units of each checkpoint's own zero-effect null dispersion for the
statistic being judged, so that `delta`, `angle`, and `shape` are commensurable despite their
incommensurable native scales. A response counts as material at a checkpoint when its excess over the
zero-effect null at that same checkpoint is at least the configured number of null-dispersion units. A
single absolute threshold applied to all three statistics MUST NOT be used, because it renders a statistic
whose entire response range is smaller than that threshold structurally unclassifiable — reporting
`not_material` for a response the corresponding test rejects at high rate.

Where the zero-effect null has no usable dispersion — an analytically null construction leaves a population
checkpoint's dispersion at zero or at floating-point dust — materiality SHALL fall back to whether the
excess exceeds that dust tolerance, which is the limit of the dispersion rule as the null variance goes to
zero: a null with no variance is exceeded by any real difference. The system MUST NOT divide by a dust-valued
dispersion, and MUST NOT report such a checkpoint as unclassifiable merely because its null is degenerate.
Each emitted row MUST record whether the dispersion path or the fallback determined its verdict.

Changing the materiality rule MUST NOT silently change how an already-reported run's localization renders:
a run reported under an earlier rule MUST either continue to reproduce its recorded classification, or be
re-issued explicitly with the rule it was judged under recorded alongside it.

#### Scenario: Geometry-aware operating report is produced

- **WHEN** merged Phase 4 records are reported
- **THEN** structured tables summarize path lengths, `delta`, `angle`, and `shape` by mode, effect,
  checkpoint, and scope beside the rejection-rate tables
- **AND** the report identifies the first checkpoint where each material off-diagonal response appears

#### Scenario: Materiality is comparable across statistics

- **WHEN** an off-diagonal response is judged at a checkpoint whose zero-effect null has usable dispersion
- **THEN** its excess over the zero-effect null is expressed in that checkpoint's null-dispersion units for
  that statistic, and the reported row records both the excess and the dispersion it was scaled by

#### Scenario: A degenerate null does not silence a checkpoint

- **WHEN** the zero-effect null's dispersion at a checkpoint is zero or floating-point dust, as it is for an
  analytically null construction at a population checkpoint
- **THEN** the response is judged material when its excess exceeds the dust tolerance, no division by the
  degenerate dispersion occurs, and the row records that the fallback path determined the verdict

#### Scenario: A small-scale statistic remains classifiable

- **WHEN** a statistic's observed response and its null are both far smaller in native units than the
  responses of the other statistics, and the corresponding test rejects well above its nominal level
- **THEN** the localization output classifies the response at a checkpoint rather than reporting
  `not_material`

#### Scenario: A previously reported run does not change silently

- **WHEN** a run reported under an earlier materiality rule is regenerated from its merged records
- **THEN** its localization output either reproduces the recorded classification, or the regeneration is an
  explicit re-issue that records which rule each classification was produced under

#### Scenario: Measurement spaces are compared

- **WHEN** feature-space and PLS-latent checkpoints appear in one report
- **THEN** the report labels their measurement spaces and compares trends or retention without treating raw
  distances as scale-equivalent

#### Scenario: PLS and attribution stability are reported

- **WHEN** eligible PLS orientation replicates are available
- **THEN** the report summarizes selected component counts, attribution availability, observed-versus-captured
  retention, cross-replicate top-k selection and sign agreement, bootstrap stability, and generator-truth
  recovery by effect and transition

### Requirement: Phase 4 emits an explicit gate for the paper-grade study

The Phase 4 report SHALL evaluate predeclared mandatory gates and emit `proceed`, `hold`, or `indeterminate` for Phase 5. It MUST NOT emit `proceed` unless every mandatory gate is met with complete eligible diagnostics.

#### Scenario: Type I inflation gate is evaluated

- **WHEN** the `none` baseline and every translation cell at each enumerated effect level are complete
- **THEN** each available statistic is checked against the predeclared one-sided Monte Carlo inflation tolerance at alpha `0.05`

#### Scenario: A single control statistic is marginally above its bound

- **WHEN** exactly one control statistic exceeds its inflation bound by less than one Monte Carlo standard error and no other mandatory gate fails
- **THEN** the report emits `indeterminate`, names the cell and statistic, and states the confirmation re-run required before a Phase 5 decision
- **AND WHEN** two or more control statistics exceed their bounds, or any exceedance is at least one Monte Carlo standard error
- **THEN** the report emits `hold`

#### Scenario: Magnitude gate is evaluated

- **WHEN** magnitude cells are complete
- **THEN** `delta` power is checked for an uncertainty-tolerant non-decreasing curve and power of at least `0.80` at effect `1.00`
- **AND** magnitude `angle` and `shape` rejection rates are checked against the Type I inflation tolerance

#### Scenario: Orientation and shape gates are evaluated

- **WHEN** orientation and shape cells are complete
- **THEN** orientation `angle` and shape `shape` power are checked for uncertainty-tolerant non-decreasing curves and power of at least `0.80` at effect `1.00`
- **AND** off-diagonal responses are reported against realized geometry without imposing a false purity requirement on the mixed constructions
- **AND** the shape gate constrains only the diagonal `shape` response, reported against the validated invariance contract, so it tests detectability of the constructed bend rather than shape-specific purity

#### Scenario: Diagnostic completeness gate is evaluated

- **WHEN** all expected work units have been merged
- **THEN** every completed PLS replicate is required to contain selected-component and realized-geometry metadata
- **AND** every eligible orientation replicate is required to contain valid attribution diagnostics or a recorded failure reason

#### Scenario: A mandatory gate fails or lacks evidence

- **WHEN** any mandatory gate fails
- **THEN** the report emits `hold` with the failed criteria and supporting observations
- **AND WHEN** expected records or required diagnostics are incomplete
- **THEN** the report emits `indeterminate` and does not recommend the paper-grade study

### Requirement: Phase 4 findings are versioned and reproducible

The study SHALL produce a dated findings report tied to the exact configuration, code revision, parameter signatures, software versions, record counts, failure counts, and reproduction commands. The prior July pilot and its outputs MUST remain unchanged and be identified as superseded for Phase 4 gate decisions.

#### Scenario: Findings report is committed

- **WHEN** the Phase 4 run and reporting complete
- **THEN** a versioned report records the gate decision, scientific interpretation, limitations, and exact shard, merge, and report commands
- **AND** all claims in the report can be traced to structured CSV or JSON outputs

### Requirement: Reporting stratifies orientation power by configuration eigengap

The study report SHALL summarize the recorded pooled and per-group eigengaps per cell beside the existing geometry summaries, and for orientation-mode power cells SHALL report rejection rates stratified by the recorded pooled eigengap, with per-stratum replicate counts and Monte Carlo uncertainty. Stratification MUST read recorded values only — the report MUST NOT regenerate datasets or recompute spectra from raw data.

#### Scenario: Eigengap summaries join the cell tables

- **WHEN** merged records carrying the spectrum block are reported
- **THEN** per-cell tables include summaries of the recorded pooled and per-group eigengaps

#### Scenario: Orientation power is stratified

- **WHEN** an orientation-mode power cell with recorded spectra is reported
- **THEN** the report includes rejection rates by eigengap stratum with per-stratum counts and uncertainty
- **AND** the strata are defined from the recorded values, not from regenerated data

#### Scenario: Legacy records degrade gracefully

- **WHEN** a report runs over merged records lacking the spectrum block
- **THEN** the stratified table is reported as unavailable for those cells and every pre-existing table is produced unchanged

### Requirement: Enumeration rejects cells that exceed surgery headroom
Study enumeration SHALL validate, at configuration time and before any compute is spent, that every
cell's requested effect can be realized without censoring in expectation: for each pool-limited
surgery mode, the expected requested surgery size (derived from the cell's generator parameters,
e.g. `p_dmp`, `n_stages`, and baseline continuity) must not exceed the expected available pool,
with the expectation computed under the cell's own baseline continuity value (the
continuity-adjusted stage-active union, reducing to the independence union at zero continuity).
Enumeration of a config containing a cell that would saturate SHALL fail with an error naming the
cell, its requested effect, and the effect at which the surgery saturates. Configs that explicitly
opt in to the clamping policy are exempt from the error but MUST still be enumerable.

#### Scenario: Censored cell fails enumeration
- **WHEN** a study config requests an orientation or translation effect beyond the expected
  headroom of its generator parameters under the default (fail-loud) censoring policy
- **THEN** enumeration raises a configuration error naming the offending cell, the requested
  effect, and the saturating effect, before any replicate is run

#### Scenario: Phase-4 pilot grid is caught
- **WHEN** the Phase-4 pilot configuration (orientation and translation at effects 0.75 and 1.00
  with `p_dmp = 0.2`, four stages) is enumerated under the default policy
- **THEN** enumeration fails on the censored cells rather than producing a grid whose top cells
  regenerate near-identical datasets

#### Scenario: Headroom-respecting grids enumerate unchanged
- **WHEN** every cell's requested effect is within the expected surgery headroom
- **THEN** enumeration succeeds and produces the same cells, seeds, and signatures it would have
  apart from the generator's new policy parameter

#### Scenario: Headroom is evaluated per continuity value
- **WHEN** a study config sweeps baseline continuity as a generator axis alongside a pool-limited
  surgery mode
- **THEN** each enumerated cell's headroom check uses that cell's own continuity value, so an
  effect rejected at zero continuity may enumerate at a higher continuity where the expected pool
  is larger

### Requirement: Reports annotate realized effects and duplicate constructions
Study reports SHALL report the realized surgery per power cell (nominal size, mean realized size,
censored fraction) and SHALL flag any pair of power cells whose realized constructions are
identical in distribution — in particular matched-seed cells whose replicates generated identical
datasets — so no two such cells are presented as independent power measurements.

#### Scenario: Realized surgery appears in the study report
- **WHEN** a study report is built from merged records that carry surgery truth metadata
- **THEN** the report includes, per power cell, the nominal surgery size, mean realized surgery
  size, and censored-replicate fraction

#### Scenario: Duplicate realized constructions are flagged
- **WHEN** two power cells in the same seed family share identical realized datasets on more than
  5% of replicate pairs
- **THEN** the report flags the pair as a duplicated construction rather than reporting the two
  cells as independent measurements

### Requirement: Reporting resolves orientation operating characteristics by baseline continuity
When a study's enumerated grid varies baseline continuity, study reporting SHALL present the
orientation-relevant operating characteristics as functions of the continuity axis: per-cell
rejection rates for each statistic, the distribution of the recorded pooled configuration eigengap,
and the dispersion of the per-replicate `angle` null width, each resolved by continuity value. The
report MUST make the linking observable explicit — the eigengap is the quantity expected to carry a
continuity-conditioned orientation claim to real data — so a reader can trace power differences
along the axis to the recorded geometry rather than to the knob itself. When the grid also varies
other design coordinates, the continuity-resolved view MUST be resolved on those coordinates as well
and MUST NOT pool records that differ in any design coordinate other than continuity.

#### Scenario: Continuity-resolved orientation table is produced
- **WHEN** a merged result set contains cells that differ only in baseline continuity
- **THEN** the report includes a per-continuity-value summary of orientation rejection rates,
  eigengap distribution summaries, and `angle` null-width dispersion, computed from the persisted
  records alone

#### Scenario: Continuity-resolved table does not pool across other design coordinates
- **WHEN** a merged result set contains design-grid cells that vary baseline continuity and at least one other design coordinate
- **THEN** every row of the continuity-resolved view is identified by continuity together with the values of every other design coordinate
- **AND** no row aggregates records from cells that differ in a design coordinate other than continuity

#### Scenario: Studies without a continuity axis are unaffected
- **WHEN** a study's grid holds baseline continuity fixed
- **THEN** reporting produces its existing outputs unchanged, without a continuity-resolved view

### Requirement: Crossed design grid enumerates one matched, anchored power grid per design point

When a study configuration declares a design grid, enumeration SHALL cross the declared axes into design points and, for every design point other than the baseline, emit the full power grid — one mode-agnostic zero-effect anchor cell plus one cell per (trajectory mode × nonzero effect size) — with the design point's coordinates applied to the baseline generator and evaluation parameters. The baseline design point MUST be served by the existing primary power cells and shared anchor rather than re-emitted. Every design-point cell MUST carry its design coordinates and its distinct phase in cell metadata so that readers of the primary power curves, the Phase 4 gate, and the acceptance targets do not see design-point cells unless they opt in. Design-point cells MUST share the primary matched-seed family so that cells across design points are paired at the same replicate index, and the existing duplicate-dataset guard and surgery-headroom check MUST apply to every design-point cell using that cell's own parameters.

#### Scenario: Design points cross the declared axes

- **WHEN** a configuration declares a design grid with axes `A` (values a₁…aₖ, including the baseline) and `B` (values b₁…bₘ, including the baseline)
- **THEN** enumeration emits design-point power grids for every (aᵢ, bⱼ) pair except the baseline pair
- **AND** each such grid contains exactly one zero-effect anchor cell and one cell per (mode, nonzero effect)
- **AND** each cell's metadata records the values of both `A` and `B`

#### Scenario: Design-point cells are invisible to baseline readers

- **WHEN** a merged result set from a design-grid study is reported
- **THEN** the primary power curves, specificity matrix, Type I table, Phase 4 gate, and acceptance-target evaluation are computed from the baseline column only and equal what a study without the design grid would report

#### Scenario: Design-point cells share the matched-seed family

- **WHEN** matched seeds are enabled and a design grid is declared
- **THEN** every design-point power cell and anchor draws the same generator seed as the primary cells at the same replicate index
- **AND** enumeration still rejects any two cells in that family whose generator parameters would produce identical datasets

#### Scenario: Headroom is enforced at every design point

- **WHEN** a design point's coordinates (for example a lower baseline continuity) leave a pool-limited surgery's requested effect above the expected headroom at that point
- **THEN** enumeration fails before any compute is spent, naming the offending cell and its saturating effect
- **AND** a design point whose coordinates provide sufficient headroom for the same effect enumerates normally

### Requirement: Reporting resolves operating characteristics by design point

When a merged result set contains design-point cells, study reporting SHALL produce a design-point operating table with one row per (design point, trajectory mode, effect size, statistic), including the baseline column and each design point's zero-effect anchor, carrying the rejection rate with its Monte Carlo standard error, the distribution of the recorded pooled configuration eigengap, the dispersion of the per-replicate `angle` null width, the distribution of the selected latent dimensionality, and the component-selection mode (cross-validated or forced) of the records in that row — all computed from persisted records alone. Studies without design-point cells MUST produce no such table.

#### Scenario: Design-point operating table is produced

- **WHEN** a merged result set from a design-grid study is reported
- **THEN** the report includes a design-point operating table whose rows are identified by every design coordinate, the trajectory mode, the effect size, and the statistic
- **AND** each row carries rejection rate, Monte Carlo standard error, eigengap distribution summaries, `angle` null-width dispersion, selected-dimensionality summaries, and the component-selection mode

#### Scenario: Design-point nulls are reported

- **WHEN** the design-point operating table is built
- **THEN** every design point contributes rows for the `none` mode at effect `0.0` from its zero-effect anchor, for each statistic

#### Scenario: Forced-rank rows are marked as such

- **WHEN** a design column's records were evaluated with a forced retained rank
- **THEN** every row of that column reports the component-selection mode as forced and its selected-dimensionality summary equals the forced rank
- **AND** rows from the cross-validated baseline column report the mode as cross-validated

#### Scenario: Studies without a design grid produce no design-point table

- **WHEN** a merged result set contains no design-point cells
- **THEN** reporting produces its existing outputs unchanged and writes no design-point operating table

### Requirement: Design-point decision is predeclared in the configuration and evaluated in the report

A study configuration MAY declare a design-point decision rule naming a target trajectory mode and statistic, a minimum power at the top enumerated effect, a confirmation standard-error threshold, and a preference order over the design-grid axes. When declared, the report SHALL evaluate the rule per design point using the target statistic's rejection rate at the largest effect enumerated in that column, classify each design point as meeting the floor with Monte Carlo confirmation, meeting it marginally, or failing it, and SHALL record a decision: the first confirmed design point in the declared preference order, or a `revise_claim` verdict when none is confirmed. The decision output MUST list every design point's classification and the zero-effect anchor's rejection rate for every statistic at that point, and MUST NOT feed the Phase 4 gate or the acceptance-target report.

#### Scenario: A confirmed design point is chosen in preference order

- **WHEN** at least one design point's target rejection rate minus the confirmation threshold times its Monte Carlo standard error is at or above the floor
- **THEN** the decision names the first such design point in the declared preference order
- **AND** lists every design point's classification and anchor rejection rates

#### Scenario: No design point is confirmed

- **WHEN** no design point meets the floor with confirmation
- **THEN** the decision records the `revise_claim` verdict together with every design point's classification, so the readiness worklist can revise the claim rather than the Monte Carlo size

#### Scenario: Decision thresholds come from the configuration

- **WHEN** the design-point decision is evaluated
- **THEN** the floor, confirmation threshold, target pair, and preference order are read from the configuration's acceptance block
- **AND** no threshold is hard-coded in report code

### Requirement: Study provides a fixed Phase 5 design-point pilot profile

The study SHALL provide a version-controlled Phase 5 design-point pilot configuration using the numpy generator, pooled PLS integration with M-value methylation, four stages, differential-site density `p_dmp = 0.1`, the default fail-loud surgery-censoring policy, a design grid crossing baseline continuity `0.0`, `0.5`, `0.8` with sample sizes `300`, `600`, `1200`, trajectory modes `orientation` and `translation`, effect sizes `0.00`, `0.25`, `0.50`, and `1.00`, 100 replicates per cell, 199 RRPP permutations, matched seeds with a shared zero-effect anchor per design point, attribution disabled, the Phase 4 gate disabled, and a design-point decision rule targeting orientation's `angle` statistic at a 0.80 floor with preference for the smallest sample size and then the lowest continuity. The configuration MUST enumerate without any censored surgery at every design point.

#### Scenario: Phase 5 pilot config is loaded

- **WHEN** the committed Phase 5 design-point pilot configuration is loaded
- **THEN** it deterministically enumerates the baseline Type I controls, the baseline primary power grid, and eight further design-point power grids with the declared modes, effects, replicate count, and permutations
- **AND** every enumerated pool-limited cell fits its expected surgery headroom
- **AND** it requires no R runtime dependency

#### Scenario: Pilot does not copy the historical clamp flag

- **WHEN** the committed Phase 5 design-point pilot configuration is inspected
- **THEN** it does not set `surgery_censoring` to `clamp`

### Requirement: Phase 5 design-point findings are versioned and reproducible

The study SHALL produce a dated design-point findings report tied to the exact configuration, code revision, parameter signatures, software versions, record counts, failure counts, and reproduction commands, recording the design-point decision, the continuity-conditional interpretation with the recorded eigengap as the linking observable, the per-design-point Type I behaviour, and the hand-off of the retained-rank question to the latent-dimensionality readiness item. The Phase 5 readiness worklist MUST record the resulting design point or the revised claim.

#### Scenario: Findings report is committed

- **WHEN** the Phase 5 design-point pilot run and reporting complete
- **THEN** a versioned report records the decision verdict, scientific interpretation, limitations, and exact shard, merge, and report commands
- **AND** all claims in the report can be traced to structured CSV or JSON outputs
- **AND** the readiness worklist's design-point item states the chosen design point or the revised claim

### Requirement: Evaluation-only design columns share datasets by design

When two design points differ only in evaluation-namespace coordinates, enumeration SHALL place them in the same matched-seed family with identical generator parameters, so that at every replicate index they evaluate the same generated dataset under different measurement settings, and the duplicate-dataset guard MUST accept them. The report MUST state, wherever such columns are compared, that the comparison is paired on identical datasets and differs only in measurement.

#### Scenario: Rank columns evaluate identical datasets

- **WHEN** a design grid declares only the axis `evaluation.integration_params.forced_components`
- **THEN** every design column's power cells and anchor share the generator parameters and replicate seeds of the baseline column's corresponding cells
- **AND** enumeration succeeds without a duplicate-dataset error

#### Scenario: Generator-identical columns that also share evaluation settings are still rejected

- **WHEN** two cells in one matched-seed family have identical generator parameters and identical evaluation parameters
- **THEN** enumeration fails with the existing duplicate-dataset error

### Requirement: Report renders the retained-rank ladder

When a study's design grid declares the retained-rank axis, the report SHALL render a rank-ladder figure showing, per trajectory mode at the largest enumerated effect and for the zero-effect anchor, each statistic's rejection rate against the retained rank with Monte Carlo error bars, with the cross-validated baseline column placed at its recorded median selected rank and marked as cross-validated. The figure MUST be built from the design-point operating table alone, and studies without the rank axis MUST produce no such figure.

#### Scenario: Rank-ladder figure is produced

- **WHEN** a merged result set from a study declaring `evaluation.integration_params.forced_components` as a design axis is reported
- **THEN** a rank-ladder figure is written beside the design-point operating table
- **AND** the cross-validated column is visibly distinguished from the forced-rank columns

#### Scenario: Studies without the rank axis produce no ladder figure

- **WHEN** a merged result set from a design-grid study that does not declare the rank axis is reported
- **THEN** no rank-ladder figure is written and every existing output is unchanged

### Requirement: Retained-rank decision is predeclared in the configuration and evaluated against the baseline column

A study configuration MAY declare a retained-rank decision rule naming the rank axis, a target (trajectory mode, statistic) pair, one or more protected (trajectory mode, statistic) pairs, a Type I bound for the zero-effect anchor, and separate standard-error multipliers for the target gain and the protected loss. When declared, the report SHALL evaluate the rule with the cross-validated column (rank value `null`) as the reference: a forced-rank column *qualifies* when (a) its target rejection rate at the largest enumerated effect exceeds the reference rate by more than the gain multiplier times the pooled Monte Carlo standard error of the difference, (b) its zero-effect anchor's rejection rate is within the Type I bound for every statistic, and (c) none of its protected rejection rates at the largest enumerated effect falls below the reference rate by more than the loss multiplier times the pooled standard error of the difference. The decision MUST be `keep_cv` when no column qualifies and `adopt_fixed_rank` naming the smallest qualifying rank otherwise. The output MUST list every column's per-criterion status and the quantities each criterion was computed from, MUST read every threshold from the configuration, and MUST NOT feed the Phase 4 gate or the acceptance-target report.

#### Scenario: No rank qualifies

- **WHEN** every forced-rank column fails at least one of the three criteria
- **THEN** the decision records `keep_cv` together with every column's per-criterion status

#### Scenario: Smallest qualifying rank is adopted

- **WHEN** ranks 6 and 9 both satisfy all three criteria and rank 4 does not
- **THEN** the decision records `adopt_fixed_rank` with rank 6

#### Scenario: A protected statistic vetoes an otherwise better rank

- **WHEN** a forced-rank column clears the target gain and the Type I bound but lowers a protected rejection rate below the reference by more than the loss multiplier times the standard error
- **THEN** that column does not qualify and the protected pair is named as the failing criterion

#### Scenario: Anchor inflation vetoes a rank

- **WHEN** a forced-rank column clears the target gain but its zero-effect anchor exceeds the Type I bound on any statistic
- **THEN** that column does not qualify and the inflated statistic is named

#### Scenario: Decision thresholds come from the configuration

- **WHEN** the retained-rank decision is evaluated
- **THEN** the axis, target pair, protected pairs, Type I bound, and both multipliers are read from the configuration's acceptance block
- **AND** declaring the rule without the named rank axis in the design grid is rejected at configuration load

### Requirement: Study provides a fixed Phase 5 latent-rank ladder profile

The study SHALL provide a version-controlled Phase 5 latent-rank ladder configuration using the numpy generator, pooled PLS integration with M-value methylation, four stages, `n_samples = 1200`, baseline continuity `0.0`, differential-site density `p_dmp = 0.1`, the default fail-loud surgery-censoring policy, trajectory modes `magnitude`, `orientation`, `shape`, and `translation`, effect sizes `0.00`, `0.25`, `0.50`, and `1.00`, a design grid over `evaluation.integration_params.forced_components` with values `null`, `3`, `4`, `6`, `9`, and `12`, 100 replicates per cell, 199 RRPP permutations, matched seeds with a shared zero-effect anchor per design point, attribution disabled, the Phase 4 gate disabled, and a retained-rank decision rule targeting orientation's `angle` statistic while protecting magnitude's `delta` and shape's `shape`. The configuration MUST enumerate without any censored surgery.

#### Scenario: Ladder config is loaded

- **WHEN** the committed Phase 5 latent-rank ladder configuration is loaded
- **THEN** it deterministically enumerates the baseline Type I controls, the baseline primary power grid, and five further design-point power grids with the declared modes, effects, replicate count, and permutations
- **AND** every enumerated pool-limited cell fits its expected surgery headroom
- **AND** every design column shares generator parameters and seeds with the baseline column
- **AND** it requires no R runtime dependency

#### Scenario: Ladder does not copy the historical clamp flag

- **WHEN** the committed Phase 5 latent-rank ladder configuration is inspected
- **THEN** it does not set `surgery_censoring` to `clamp`

### Requirement: Phase 5 latent-rank findings are versioned and the readiness worklist records the committed rank rule

The study SHALL produce a dated latent-rank findings report tied to the exact configuration, code revision, parameter signatures, software versions, record counts, failure counts, and reproduction commands, recording the retained-rank decision, each statistic's response to rank per mode beside the recorded eigengap and `angle` null-width dispersion, whether the orientation→shape response decays with rank at the chosen design point, the re-measured magnitude and shape operating characteristics at the chosen design point under the cross-validated rank, and the group-blind statement of the committed rank rule. The Phase 5 readiness worklist MUST record the committed rank rule as the integration configuration Phase 5 runs with, and the roadmap's Phase 5 planned baseline MUST name it.

#### Scenario: Findings report is committed

- **WHEN** the Phase 5 latent-rank ladder run and reporting complete
- **THEN** a versioned report records the decision verdict, scientific interpretation, limitations, and exact shard, merge, and report commands
- **AND** all claims in the report can be traced to structured CSV or JSON outputs
- **AND** the readiness worklist's latent-dimensionality item states the committed rank rule and the roadmap's Phase 5 planned baseline names it

### Requirement: Study configuration declares a report contract and rejects unknown top-level keys

The declarative configuration SHALL accept an optional `report_contract` block declaring which attribution component driver reports present (`observed`, `pls_captured`, or `residual`), that cross-replicate driver agreement is descriptive only, and whether a worker-count override at run time is forbidden or merely warned about. Loading a configuration MUST fail with a message naming the offending key when the root mapping contains a key the schema does not define or when the contract contains an unknown field or value. The contract MUST survive a dump-and-reload round trip. Configurations without the block MUST load and report exactly as before.

#### Scenario: Contract is loaded

- **WHEN** a configuration declares `report_contract` with `driver_component` `observed`, `cross_replicate_driver_agreement` `descriptive`, and `n_jobs_override` `forbid`
- **THEN** the loaded configuration exposes those values and writing it back and reloading it yields the same contract

#### Scenario: Unknown key is rejected

- **WHEN** a configuration contains an unknown root key, an unknown `report_contract` field, or an `n_jobs_override`/`driver_component` value outside the declared set
- **THEN** loading fails and the error names the key or value

#### Scenario: Absent contract is inert

- **WHEN** every committed configuration under the study examples directory is loaded and enumerated
- **THEN** each loads without error and its cell identities and parameter signatures are unchanged from the recorded snapshot

### Requirement: Report echoes the resolved contract and identifies the shared zero-effect anchor as one measurement

When a configuration declares a report contract, reporting SHALL write a machine-readable contract echo containing the declared driver component, the descriptive-only status of cross-replicate driver agreement, the worker count the merged records were produced with (which MUST be uniform across records), and the shared zero-effect anchor's cell identifier together with the trajectory modes whose zero-effect power point it resolves and the statement that it is counted as one measurement. Operating-characteristic rows derived from the anchor MUST remain flagged as anchor-derived.

#### Scenario: Contract echo is written

- **WHEN** reporting runs on merged records for a configuration with a report contract and a shared zero-effect anchor
- **THEN** the report directory contains a contract echo naming the driver component, the anchor cell, the modes it resolves, and the uniform worker count

#### Scenario: Non-uniform worker count is refused

- **WHEN** merged records carry more than one worker-count value
- **THEN** reporting fails naming the values found rather than echoing one of them

### Requirement: Driver reporting presents the declared component and makes no cross-replicate stability claim

When a configuration declares a report contract, reporting SHALL write a driver table restricted to the declared attribution component, per trajectory mode, effect size, and transition, carrying truth precision and recall, mean selected count, within-replicate bootstrap stability, and replicate accounting, and containing no cross-replicate agreement column. The attribution figure SHALL plot within-replicate bootstrap stability only and SHALL be titled as such. The complete per-component attribution table, including cross-replicate agreement, MUST continue to be written unchanged, and configurations without a contract MUST produce exactly the prior outputs.

#### Scenario: Driver table is restricted to the declared component

- **WHEN** reporting runs with a contract declaring the `observed` component on records that carry attribution diagnostics for all three components
- **THEN** the driver table contains rows for the `observed` component only and has no cross-replicate agreement column
- **AND** the full attribution table still contains every component and the cross-replicate columns

#### Scenario: Legacy outputs are byte-identical without a contract

- **WHEN** reporting re-runs on the committed Phase 4, design-point, and latent-rank fixtures, whose configurations declare no contract
- **THEN** every existing report file is byte-identical to the committed output and no driver table or contract echo is written

### Requirement: Shard runner enforces the configured worker count when the contract forbids overrides

When the configuration's report contract sets the worker-count override policy to forbid, the shard runner SHALL refuse a command-line worker count that differs from the configured value, exiting non-zero with a message naming both values and stating that the override changes the permutation draws and the cell parameter signature. A command-line value equal to the configured value SHALL be accepted. Under the warn policy, or with no contract, the runner SHALL keep its existing warn-and-proceed behavior. The SLURM array template MUST document that a forbidding configuration refuses the optional worker-count environment variable.

#### Scenario: Differing override is refused

- **WHEN** the shard runner is invoked with a worker count differing from a forbidding configuration
- **THEN** it exits non-zero before enumerating or running any unit and names both values

#### Scenario: Equal or absent override proceeds

- **WHEN** the shard runner is invoked without a worker count, or with one equal to the configured value, under a forbidding configuration
- **THEN** it runs the shard normally

### Requirement: Study provides a fixed Phase 5 paper-grade profile

The study SHALL provide a version-controlled Phase 5 paper-grade configuration using the numpy generator, pooled PLS integration with M-value methylation and cross-validated component selection, four stages, `n_samples = 1200`, baseline continuity `0.0`, `p_dmp = 0.1`, the default fail-loud surgery-censoring policy, trajectory modes `magnitude`, `orientation`, `shape`, and `translation`, effect sizes `0.00`, `0.25`, `0.50`, `0.75`, and `1.00`, 500 replicates per cell, 999 RRPP permutations, a worker count of one with overrides forbidden, matched seeds with one shared zero-effect anchor in a family distinct from both Phase 5 pilots, attribution enabled on nonzero-effect orientation primary cells, the Phase 4 gate enabled with mandatory power rules on the three diagonal statistics, mandatory control rules on magnitude's off-diagonals, and descriptive rules on the remaining off-diagonals, acceptance targets whose specificity entries are exactly the gate's mandatory controls, and a report contract declaring the `observed` driver component. The configuration MUST enumerate without any censored surgery and MUST derive its generator and integration parameters from the committed latent-rank ladder profile.

#### Scenario: Paper-grade config is loaded

- **WHEN** the committed Phase 5 paper-grade configuration is loaded
- **THEN** it deterministically enumerates the Type I controls, one shared zero-effect anchor, and sixteen nonzero power cells at 500 replicates and 999 permutations, every pool-limited cell fits its expected surgery headroom, and it requires no R runtime dependency

#### Scenario: Paper-grade config does not copy the historical clamp flag or the failing target

- **WHEN** the committed Phase 5 paper-grade configuration is inspected
- **THEN** it does not set `surgery_censoring` to `clamp` and its specificity targets contain no entry for orientation's `shape` statistic

### Requirement: Phase 5 findings follow a committed report template

The study SHALL commit a findings-report template beside the Phase 5 profile that the dated Phase 5 report MUST follow. The template SHALL require: configuration and provenance including the provenance file's field list; unit and failure accounting; the gate decision; a Type I section stating which cell each null claim reads and that the shared anchor is one measurement; per-mode power with the recorded eigengap distribution and `angle` null-width dispersion beside every orientation result; the predeclared orientation-to-`shape` cross-talk statement; a driver section limited to the declared component with within-replicate bootstrap stability and an explicit statement that no cross-replicate driver-stability claim is made and why; construction limitations including the n-conditional orientation claim and the non-ρ-invariant realized orientation contrast; and reproduction commands that pass no worker-count override. The readiness worklist MUST record item 5 as resolved by the committed contract, profile, and template.

#### Scenario: Template names every contract item

- **WHEN** the committed template is inspected
- **THEN** it contains a section for each required element above, and the contract items it names agree with the committed profile's report contract

#### Scenario: Readiness worklist records item 5

- **WHEN** the readiness worklist is inspected after this change
- **THEN** item 5 is marked resolved with references to the profile, the template, and the contract fields that enforce each of its four points

### Requirement: Study provides fixed Phase 5 magnitude re-measurement profiles

The study SHALL provide two version-controlled Phase 5 magnitude re-measurement configurations — a pilot and a paper-grade profile — that re-measure the magnitude mode under the `joint` magnitude kind at the Phase 5 paper-grade design point. Both MUST derive their generator parameters (other than `magnitude_kind`), integration parameters, base seed, and matched-seed family from the committed Phase 5 paper-grade profile, so that the shared zero-effect anchor and every group A baseline are byte-identical to the Phase 5 records at the same replicate index. Both MUST set `trajectory_modes` to `magnitude` only, `generator.magnitude_kind` to `joint`, attribution disabled, the default fail-loud surgery-censoring policy, and effect sizes consisting of `0.00`, `1.00`, and at least two further nonzero values strictly below `0.25` chosen from the committed analytic effect-axis bracket so the `delta` power curve is measured below saturation. The pilot MUST use 50 replicates per cell and 199 RRPP permutations with the Phase 4 gate disabled. The paper-grade profile MUST use 500 replicates per cell, 999 RRPP permutations, a worker count of one with overrides forbidden, the same report contract as the Phase 5 paper-grade profile, and the Phase 4 gate enabled with exactly three rules — mandatory power on magnitude's `delta`, mandatory control on magnitude's `angle` and `shape` — with `none` as the only control mode, and acceptance targets whose specificity entries are exactly the gate's mandatory controls.

#### Scenario: Pilot config is loaded

- **WHEN** the committed Phase 5 magnitude pilot configuration is loaded
- **THEN** it deterministically enumerates one shared zero-effect anchor and one nonzero power cell per nonzero effect, all at 50 replicates and 199 permutations under the joint magnitude kind, and requires no R runtime dependency

#### Scenario: Paper-grade re-measurement config is loaded

- **WHEN** the committed Phase 5 magnitude paper-grade configuration is loaded
- **THEN** it deterministically enumerates one shared zero-effect anchor and one nonzero power cell per nonzero effect, all at 500 replicates and 999 permutations under the joint magnitude kind, with the gate's rules naming only the magnitude mode and the control modes naming only `none`

#### Scenario: Anchor reproduces the Phase 5 anchor

- **WHEN** the re-measurement's zero-effect anchor cell is generated at any replicate index
- **THEN** its dataset is identical to the Phase 5 paper-grade anchor's dataset at the same replicate index

#### Scenario: Profiles do not touch the Phase 5 artefacts

- **WHEN** the two re-measurement profiles are added
- **THEN** the committed Phase 5 paper-grade configuration, its results directory, its report, the gate role definitions, the acceptance-target semantics, the report contract, and the report template are unchanged

### Requirement: Phase 5 magnitude re-measurement findings are versioned as an addendum

The study SHALL produce a dated magnitude re-measurement findings report, tied to the exact configuration, code revision, parameter signatures, software versions, record counts, failure counts, and reproduction commands, that follows the committed Phase 5 report template section for section with non-magnitude sections marked not applicable. It MUST record the reduced gate's decision, the `delta` power curve on the joint construction's native effect axis beside the realized joint `delta` per cell, the two control rates at every effect against the α + 2·SE bound, whether the anchor records reproduced the Phase 5 anchor, and the construction's biological framing (every tier of the methylation→expression→protein cascade scales together). It MUST either lift the Phase 5 exit review's withheld magnitude-specificity claim or report the claim as failed on the construction's own terms, and MUST explain any deviation from a predeclared target as a method or claim revision rather than a Monte Carlo sample-size question. The roadmap MUST record the outcome.

#### Scenario: Addendum report is committed

- **WHEN** the paper-grade magnitude re-measurement run and reporting complete
- **THEN** a versioned addendum records the gate verdict, the power curve with realized sizes, the control rates, the anchor reproduction check, limitations, and exact shard, merge, and report commands, and every claim traces to a committed CSV or JSON output

#### Scenario: Specificity claim is resolved, not deferred

- **WHEN** the addendum is written
- **THEN** it states explicitly whether the withheld magnitude-specificity claim is lifted, and if not, which control failed and what revision that implies

### Requirement: Study configurations accept explicit group × stage sizes and block selection
The declarative study configuration SHALL accept the generator's explicit group × stage size table and the evaluation's layer selection. Each SHALL be accepted both in the baseline blocks and as an axis of a design grid. A size table given as nested JSON lists SHALL be normalized so that equal tables produce equal parameter signatures, cell identifiers, and matched seeds whatever their JSON spelling. Design-grid columns that differ only in the layer selection SHALL share generator parameters and matched seeds, as for every evaluation-namespace axis. A configuration that does not use either key SHALL enumerate exactly as before.

#### Scenario: Size table round-trips through a configuration
- **WHEN** a configuration sets `generator.group_stage_sizes` to `[[11, 10, 28], [9, 10, 12]]`
- **THEN** the configuration loads, enumerates, and generates datasets with those exact cell sizes, and reloading the same file yields identical cell identifiers and parameter signatures

#### Scenario: Size table as a design-grid axis
- **WHEN** a design grid declares `generator.group_stage_sizes` with two tables
- **THEN** enumeration emits one design point per table, each with its own anchored power grid

#### Scenario: Layer selection columns share data
- **WHEN** two design-grid columns differ only in `evaluation.integration_params.layers`
- **THEN** at every replicate index they evaluate the same generated dataset

#### Scenario: Committed configurations are unchanged
- **WHEN** any configuration committed before this change is enumerated
- **THEN** its cell identifiers, parameter signatures, and matched seeds are identical to those enumerated before this change

### Requirement: Phase 6 small-n study profiles
The study SHALL provide version-controlled Phase 6 small-n configurations, a pilot and a paper-grade profile, that reproduce the SEA-AD MTG astrocyte case-study design. Both SHALL use the numpy generator. Both SHALL use three stages and the explicit size table `((11, 10, 28), (9, 10, 12))` from the ≥30-nuclei cohort. Both SHALL use the layer selection methylation + expression as the two-block analogue of ATAC + RNA. Both SHALL use pooled PLS integration with M-value methylation and cross-validated component selection, the `joint` magnitude construction, the fail-loud surgery-censoring policy, matched seeds with a shared zero-effect anchor, and, across the profiles of the paper-grade study, the trajectory modes `magnitude`, `orientation`, `shape`, and `translation`. The paper-grade profile SHALL declare a design grid that crosses the ≥30-nuclei and ≥50-nuclei tables (`((10, 9, 25), (9, 9, 12))`) with the two-block and three-block layer selections. Its columns therefore include the baseline, the ≥50-nuclei table on two blocks, and the baseline table on all three blocks. The pilot SHALL measure the baseline column only, and SHALL exist to bracket the effect axes, which every column shares. The paper-grade profile SHALL fix its effect axis from recorded pilot evidence and SHALL record that evidence in its metadata. If the pilot shows that no single effect axis resolves the rise of every mode's target statistic, the paper-grade study SHALL be split into one profile per axis, as the Phase 5 magnitude re-measurement was. The split profiles SHALL share base seed, matched-seed family, generator, and evaluation parameters, so that their shared zero-effect anchors are byte-identical. The paper-grade profile SHALL run at least 500 replicates per cell and 999 permutations, and SHALL declare the Phase 5 report contract and the Phase 4 gate rules as advisory acceptance targets.

#### Scenario: Profiles load and enumerate without censoring
- **WHEN** the pilot or paper-grade profile is loaded and enumerated
- **THEN** it validates, every cell passes the configuration-time surgery-headroom check, and every cell's generator produces 80 samples (or 74 in the ≥50-nuclei column) in the declared group × stage sizes

#### Scenario: Effect axes are evidence-based
- **WHEN** the paper-grade profile is read
- **THEN** its metadata names the pilot result directory and states, per trajectory mode, why each enumerated effect was chosen

### Requirement: Phase 6 small-n findings report
The study SHALL produce a dated findings report for the paper-grade Phase 6 small-n run. The report SHALL follow the committed Phase 5 report template section for section, and SHALL be tied to the exact configuration, code revision, parameter signatures, software versions, record counts, failure counts, and reproduction commands. The report SHALL state, for each of `delta`, `angle`, and `shape`, the Type I rate and the power at each enumerated effect at the SEA-AD design, with Monte Carlo uncertainty. It SHALL state the recorded pooled-eigengap distribution of the zero-effect anchor and the orientation power stratified by eigengap tercile. It SHALL state the change in each operating characteristic from the baseline to the ≥50-nuclei column and to the three-block column. It SHALL conclude with an explicit per-statistic interpretability statement for the case study: whether a rejection, and whether a non-rejection, of that statistic may be reported as evidence at n = 80.

#### Scenario: Report states case-study interpretability per statistic
- **WHEN** the findings report is complete
- **THEN** it contains, for each of `delta`, `angle`, and `shape`, a statement of whether that statistic's result on the SEA-AD cohort is interpretable, grounded in the measured Type I and power at the SEA-AD design

#### Scenario: Report separates sample size from block count
- **WHEN** the findings report compares design columns
- **THEN** it attributes operating-characteristic differences to cohort size (baseline vs ≥50-nuclei column) and to block count (baseline vs three-block column) separately

### Requirement: Explicit design points enumerate exactly like crossed design points

When a design grid is declared in the explicit form, enumeration SHALL treat every listed non-baseline point exactly as it treats a non-baseline point of a crossed grid: it emits that point's anchored power grid under the `power_design` phase, stamps the point's coordinates in cell metadata, places the cells in the primary matched-seed family, and applies the duplicate-dataset guard and the surgery-headroom check with the point's own parameters. The baseline point MUST be served by the primary cells and never re-emitted. A configuration whose explicit point list equals the full cross of some axes MUST enumerate the same cell identifiers and parameter signatures as the crossed declaration of those axes.

#### Scenario: Explicit points enumerate their own anchored grids

- **WHEN** a configuration lists the baseline point and four further design points
- **THEN** enumeration emits four `power_design` grids, each with one zero-effect anchor and one cell per (mode, nonzero effect), each carrying its point's coordinates

#### Scenario: Explicit and crossed forms agree on a rectangle

- **WHEN** one configuration declares axes `A` and `B` in the crossed form and another lists every (aᵢ, bⱼ) pair explicitly, with everything else equal
- **THEN** both enumerate identical cell identifiers, parameter signatures and replicate seeds

### Requirement: Orientation estimator and integration method are design axes

The evaluation-namespace axes `evaluation.orientation_estimator` and `evaluation.integration_method` SHALL be accepted as design-grid axes (in either form) and as one-factor-at-a-time axes. The baseline value of `evaluation.orientation_estimator` is the baseline evaluation's estimator: the evaluator's default (net displacement) unless the baseline evaluation sets it. In a configuration, `null` SHALL resolve to the evaluator's default, as it does for the nested `evaluation.integration_params.<key>` form, and `"pc1"` SHALL name the PC1 estimator; `null` and `"pc1"` are therefore distinct axis values. Cells at a design point SHALL evaluate with the point's estimator and integration method, every RRPP draw and realized-geometry checkpoint included, and cells at the baseline column SHALL write records byte-identical to a configuration without these axes. The design-point operating table SHALL carry both axes as columns and SHALL accept columns whose integration method records no PLS component selection or selected rank, reporting those fields as unavailable rather than failing. The production forced-rank guard SHALL admit a forced-rank column declared through the explicit form on the same terms as one declared through the crossed form: its records pass when the recorded rank equals the rank at their design point.

#### Scenario: Estimator axis is applied per design point

- **WHEN** a configuration whose baseline evaluation names `"pc1"` declares the design axis `evaluation.orientation_estimator` with values `["pc1", "net_displacement"]`
- **THEN** the net-displacement column's records carry `orientation_estimator: "net_displacement"` in their runtime metadata and their `angle` p-values come from that estimator
- **AND** the baseline column's records carry no estimator entry and equal the records of the same configuration without the axis

#### Scenario: Null names the evaluator's default

- **WHEN** a configuration whose baseline evaluation omits the estimator declares the design axis `evaluation.orientation_estimator` with values `[null, "pc1"]`
- **THEN** the `null` point is the baseline column and evaluates with the net-displacement estimator, and the `"pc1"` column's records carry no estimator entry

#### Scenario: Non-PLS design columns are reported

- **WHEN** a merged result set contains design-point records whose integration method is `concat` or `pca`
- **THEN** the design-point operating table lists their rows with the rejection rate, Monte Carlo standard error, eigengap and `angle` null-width covariates, and leaves the component-selection and selected-rank fields unavailable
- **AND** reporting completes without error

#### Scenario: Forced-rank column declared as an explicit point passes the guard

- **WHEN** an explicit design point sets `evaluation.integration_params.forced_components` to 4 and its records record a forced selection at rank 4
- **THEN** the production component-selection guard admits those records

### Requirement: Anchor Type I decision is predeclared in the configuration and evaluated in the report

A study configuration MAY declare an advisory anchor Type I decision rule naming a significance level, a standard-error multiplier, and the statistics it applies to. When declared, the report SHALL evaluate every design point's zero-effect anchor, the baseline column included, against the bound α + k·√(α(1 − α)/R), where R is the number of available replicates of that anchor and statistic, and SHALL classify each (design point, statistic) as `holds` when the anchor rejection rate is at or below the bound and `exceeds` otherwise, with `unavailable` when no replicate is available. The decision output MUST list, per design point, every statistic's rate, Monte Carlo standard error, replicate count, bound and classification, plus the rejection rates of every other negative-control cell at that point as descriptive context, and MUST record a per-point verdict: `size_holds` when every named statistic holds, `excess` when any named statistic exceeds, and `unavailable` when none exceeds but some named statistic has no available replicate. The rule MUST NOT feed the Phase 4 gate or the acceptance-target report, and a configuration without the rule MUST produce outputs byte-identical to before the rule existed.

#### Scenario: Anchor rate at the bound holds

- **WHEN** an anchor's `angle` rejection rate is 0.068 over 500 available replicates at α = 0.05 and k = 2
- **THEN** the bound is 0.0695, the classification is `holds`, and the point's verdict is `size_holds` when its other named statistics also hold

#### Scenario: Anchor rate above the bound is an excess

- **WHEN** an anchor's `angle` rejection rate is 0.080 over 500 available replicates at α = 0.05 and k = 2
- **THEN** the classification is `exceeds` and the point's verdict is `excess`, whatever the other statistics show

#### Scenario: Rule absent leaves every output unchanged

- **WHEN** a configuration declares no anchor decision rule
- **THEN** no anchor decision output is written and every other report output is byte-identical to before the rule existed

### Requirement: Study provides fixed net-displacement Type I profiles

The study SHALL provide two version-controlled net-displacement Type I profiles that re-measure the zero-effect anchor and the `translation` negative control at paper-grade size across five measurement columns on matched datasets. Both SHALL use the numpy generator, `p_dmp = 0.1`, baseline continuity `0.0`, the joint magnitude construction, the default fail-loud surgery-censoring policy, the trajectory mode `translation` only with effect sizes `0.0` and `1.0`, 500 replicates per cell, 999 RRPP permutations, matched seeds with a shared zero-effect anchor, attribution disabled, the Phase 4 gate disabled, Type I and translation specificity targets at α = 0.05 with a two-SE tolerance, and an anchor decision rule at α = 0.05 with k = 2 over `delta`, `angle` and `shape`. The design grid SHALL be declared in the explicit form over the axes `evaluation.integration_method`, `evaluation.integration_params.forced_components`, `evaluation.integration_params.n_components` and `evaluation.orientation_estimator`, with the points: the baseline (cross-validated stage PLS, PC1); cross-validated stage PLS with net displacement; stage PLS forced to rank 2(K − 1) with net displacement; PCA at rank 2(K − 1) with net displacement; and standardized concatenation with net displacement. The **Phase 6** profile SHALL use three stages, `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, the methylation and expression blocks, and the Phase 6 small-n study's base seed and matched-seed family, so that its baseline anchor evaluates the same datasets as the committed Phase 6 small-n run. The **Phase 5** profile SHALL use `n_samples = 1200`, four stages, `group_ratio = 0.5`, all three blocks, and the Phase 5 primary base seed and matched-seed family. Neither profile may set `surgery_censoring` to `clamp`.

#### Scenario: Phase 6 Type I profile is loaded

- **WHEN** the committed Phase 6 net-displacement Type I profile is loaded and enumerated
- **THEN** it deterministically enumerates two Type I baseline cells, the primary anchor and translation cells, and four further design columns of one anchor and one translation cell each, twelve cells in all, at 500 replicates and 999 permutations
- **AND** its baseline anchor cell draws the same generator seeds as the Phase 6 small-n study's shared anchor at every replicate index
- **AND** every pool-limited cell fits its expected surgery headroom

#### Scenario: Phase 5 Type I profile is loaded

- **WHEN** the committed Phase 5 net-displacement Type I profile is loaded and enumerated
- **THEN** it deterministically enumerates the same twelve-cell structure at `n_samples = 1200` and four stages
- **AND** its baseline anchor cell draws the same generator seeds as the Phase 5 primary family's shared anchor at every replicate index, and the first replicates' datasets equal those of the Phase 5 paper-grade anchor

#### Scenario: Committed profiles are identity-pinned

- **WHEN** either profile is enumerated
- **THEN** its cell identifiers, parameter signatures and matched seeds equal the committed identity fixture, and every previously committed configuration's fixture entry is unchanged

### Requirement: Follow-up selection is computed from both profiles' anchor decisions

The study SHALL provide a follow-up selection that reads the anchor decision outputs of the Phase 6 and Phase 5 net-displacement Type I profiles and writes a machine-readable decision. The lead arm is cross-validated stage PLS with net displacement. The verdict SHALL be `paper_grade_power_comparison` when the lead arm's verdict is `size_holds` at both designs and `diagnose_anchor_excess` otherwise. The output MUST also carry a descriptive localization: `every_column` when every net-displacement column exceeds at a design, `some_columns` when only some do, and `none` when none does, per design. The selection is advisory and MUST state the rule version and the two input files it read.

#### Scenario: Lead arm holds size at both designs

- **WHEN** the lead arm's anchor verdict is `size_holds` in both inputs
- **THEN** the verdict is `paper_grade_power_comparison`

#### Scenario: Lead arm exceeds at one design

- **WHEN** the lead arm's anchor verdict is `excess` in either input
- **THEN** the verdict is `diagnose_anchor_excess` and the localization names, per design, whether every net-displacement column, some, or none exceeded

### Requirement: Study produces a dated net-displacement Type I findings report

The study SHALL produce a dated findings report for the net-displacement Type I run, tied to the exact configurations, code revision, parameter signatures, software versions, record counts, failure counts and reproduction commands. The report SHALL present, per design and per column, the anchor rejection rate of every statistic with its Monte Carlo standard error against the bound, the translation control's rates, the recorded eigengap and `angle` null-width covariates, the anchor decision and the follow-up verdict, and SHALL state which follow-up change the verdict selects. It SHALL state that the baseline column of each profile re-evaluates the datasets of the corresponding committed study with the same evaluation parameters and RRPP seed, and report whether its baseline records reproduce that study's recorded p-values exactly. It MUST NOT present any power figure or any orientation power claim.

#### Scenario: Report is tied to its run

- **WHEN** the findings report is written
- **THEN** it names both configurations, the code revision, the software versions, the record and failure counts, the committed results directories, and the commands that reproduce the merge, the report and the follow-up selection

#### Scenario: Report names the follow-up

- **WHEN** the follow-up selection has been computed
- **THEN** the report quotes its verdict and names the next change it selects, and the roadmap entry for item 3 is updated to the branch taken

### Requirement: Column decision is predeclared in the configuration and evaluated against the baseline column

A study configuration MAY declare an advisory column decision rule that judges one **candidate** design column against the baseline column on the matched power grid. The rule SHALL name the candidate by design-axis coordinates that resolve, after the same value normalization the design grid applies, to exactly one non-baseline design point; a target mode × statistic pair with a criterion kind of `gain` or `no_loss`; a gain standard-error multiplier and a loss standard-error multiplier; a Type I bound (significance level and standard-error tolerance); a non-empty list of **mandatory control** mode × statistic pairs; and a list of **protected** mode × statistic pairs. Every named mode MUST be enumerated by the study, the target pair MUST NOT also be a control or protected pair, and the loader MUST reject a candidate that matches no design point, more than one, or the baseline.

When declared, the report SHALL measure, in the baseline column and the candidate column, each named mode's rejection rate at the largest enumerated effect of that mode and each column's zero-effect anchor rate, with Monte Carlo standard errors and available replicate counts, and SHALL evaluate four criteria for the candidate:

- **target** — under `gain`, the candidate's target rate exceeds the baseline's by more than the gain multiplier times the pooled standard error; under `no_loss`, the baseline's target rate exceeds the candidate's by at most the loss multiplier times the pooled standard error;
- **anchor** — the candidate's zero-effect rate is at or below α + k·√(α(1 − α)/R) on `delta`, `angle` and `shape`, with α and k the rule's Type I bound and R the available anchor replicates;
- **controls** — every mandatory control pair's rate in the candidate column, at that mode's largest effect, is at or below the same one-sided bound computed with that pair's available replicates;
- **protected** — no protected pair's rate in the candidate column falls below the baseline's by more than the loss multiplier times the pooled standard error.

The verdict SHALL be `candidate_qualifies` when every criterion passes, `candidate_fails` when any criterion fails, naming the failing criteria and pairs, `unavailable` when no criterion fails but one cannot be evaluated, `no_candidate` when no record carries the candidate's coordinates, and `no_design_grid` when no design-point record is present. The output MUST list every rate, standard error, replicate count, threshold and bound each criterion used, MUST report the same measures for every other design column as descriptive context without a verdict, and MUST be written as a JSON document and a CSV table. The rule MUST NOT feed the Phase 4 gate or the acceptance-target report, and a configuration without the rule MUST produce outputs byte-identical to before the rule existed.

#### Scenario: Candidate qualifies on a gain criterion

- **WHEN** the rule declares `gain` on orientation/`angle` with multiplier 2, the candidate's rate at the top orientation effect is 0.670 and the baseline's is 0.300 over 500 replicates each, the candidate's anchor is at or below 0.0695 on every statistic, every mandatory control is at or below its bound, and no protected pair loses more than 2 pooled SE
- **THEN** the target criterion passes (the gain 0.370 exceeds 2 × 0.0294) and the verdict is `candidate_qualifies`

#### Scenario: Mandatory control above the bound fails the candidate

- **WHEN** the rule names translation/`angle` as a mandatory control and the candidate's translation `angle` rate at the top translation effect is 0.080 over 500 available replicates at α = 0.05 and k = 2
- **THEN** the controls criterion fails naming translation/`angle` against the bound 0.0695, and the verdict is `candidate_fails` whatever the other criteria show

#### Scenario: Candidate qualifies on a no-loss criterion

- **WHEN** the rule declares `no_loss` on orientation/`angle` with loss multiplier 2, the baseline's rate at the top orientation effect is 0.860 and the candidate's is 0.840 over 500 replicates each, and the other three criteria pass
- **THEN** the target criterion passes (the loss 0.020 is within 2 pooled SE) and the verdict is `candidate_qualifies`

#### Scenario: Candidate coordinates must resolve to one non-baseline point

- **WHEN** a configuration declares a column decision whose candidate coordinates match no design point, or match the baseline point, or match more than one point
- **THEN** loading the configuration raises a clear validation error naming the candidate coordinates and the problem

#### Scenario: Rule absent leaves every output unchanged

- **WHEN** a configuration declares no column decision rule
- **THEN** no column decision output is written and every other report output is byte-identical to before the rule existed

### Requirement: Study provides fixed net-displacement power comparison profiles

The study SHALL provide two version-controlled net-displacement power comparison profiles that measure the four trajectory modes at paper-grade size across two measurement columns on matched datasets. Both SHALL use the numpy generator, `p_dmp = 0.1`, baseline continuity `0.0`, the joint magnitude construction, the default fail-loud surgery-censoring policy, the trajectory modes `magnitude`, `orientation`, `shape` and `translation` with effect sizes `0.0`, `0.25`, `0.5`, `0.75` and `1.0`, 500 replicates per cell, 999 RRPP permutations, matched seeds with a shared zero-effect anchor, attribution disabled, the Phase 4 gate disabled, no power target, and Type I and specificity targets (translation × `delta`/`angle`/`shape`, magnitude × `angle`/`shape`) at α = 0.05 with a two-SE tolerance in the **upper** bound form. The design grid SHALL be declared in the explicit form over the axes `evaluation.integration_method` and `evaluation.orientation_estimator`, with the points: the baseline (cross-validated stage PLS, PC1) and cross-validated stage PLS with net displacement (the candidate). Both SHALL declare the column decision rule with the candidate at the net-displacement stage-PLS point, target orientation/`angle`, mandatory controls translation/`angle`, magnitude/`angle` and shape/`angle`, protected magnitude/`delta` and shape/`shape`, Type I bound α = 0.05 with k = 2, and multipliers 2. The **Phase 6** profile SHALL use three stages, `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, the methylation and expression blocks, the Phase 6 small-n study's base seed and matched-seed family, and the `gain` criterion. The **Phase 5** profile SHALL use `n_samples = 1200`, four stages, `group_ratio = 0.5`, all three blocks, the Phase 5 primary base seed and matched-seed family, and the `no_loss` criterion. Neither profile may set `surgery_censoring` to `clamp`.

#### Scenario: Phase 6 comparison profile is loaded

- **WHEN** the committed Phase 6 net-displacement power comparison profile is loaded and enumerated
- **THEN** it deterministically enumerates two Type I baseline cells, the primary anchor and sixteen primary power cells, and one further design column of one anchor and sixteen power cells, thirty-six cells in all, at 500 replicates and 999 permutations
- **AND** its baseline cells draw the same generator seeds as the Phase 6 small-n study's cells of the same mode and effect at every replicate index
- **AND** every pool-limited cell fits its expected surgery headroom

#### Scenario: Phase 5 comparison profile is loaded

- **WHEN** the committed Phase 5 net-displacement power comparison profile is loaded and enumerated
- **THEN** it deterministically enumerates the same thirty-six-cell structure at `n_samples = 1200` and four stages
- **AND** its baseline cells draw the same generator seeds as the Phase 5 primary family's cells of the same mode and effect at every replicate index, and the first replicates' orientation datasets equal those generated from the Phase 5 paper-grade profile at the same effect

#### Scenario: Committed profiles are identity-pinned

- **WHEN** either profile is enumerated
- **THEN** its cell identifiers, parameter signatures and matched seeds equal the committed identity fixture, and every previously committed configuration's fixture entry is unchanged

### Requirement: Comparison selection is computed from both profiles' column decisions

The study SHALL provide a comparison selection that reads the column decision outputs of the Phase 6 and Phase 5 net-displacement power comparison profiles and writes a machine-readable decision. The verdict SHALL be `adopt_net_displacement` when the candidate's verdict is `candidate_qualifies` at both designs and `keep_pc1` otherwise. The output MUST name, per design, the candidate's verdict and the failing criteria and pairs when it does not qualify, and MUST state that the selection is advisory, the rule version, and the two input files it read with their content hashes.

#### Scenario: Candidate qualifies at both designs

- **WHEN** the candidate's verdict is `candidate_qualifies` in both inputs
- **THEN** the verdict is `adopt_net_displacement`

#### Scenario: Candidate fails at one design

- **WHEN** the candidate's verdict is `candidate_fails` or `unavailable` in either input
- **THEN** the verdict is `keep_pc1` and the output names, per design, the failing criteria and pairs

### Requirement: Study produces a dated net-displacement power comparison findings report

The study SHALL produce a dated findings report for the net-displacement power comparison run, tied to the exact configurations, code revision, parameter signatures, software versions, record counts, failure counts and reproduction commands. The report SHALL present, per design: each mode's `angle` rejection rate at every effect under both estimators with Monte Carlo standard errors, with the probe's concatenation ceiling cited as descriptive context; the anchor and every mandatory control against the bound; the recorded eigengap and `angle` null-width covariates per column; the column decision and the comparison verdict. It SHALL state that `delta` and `shape` cannot differ between the two stage-PLS columns, because the estimator enters only the orientation statistic, and report whether their p-values are identical record for record. It SHALL state that the baseline column re-evaluates the datasets of the corresponding committed studies with the same evaluation parameters and RRPP seed, and report whether its records reproduce those studies' recorded p-values exactly wherever a parent record of the same mode and effect exists. When the verdict is `adopt_net_displacement`, the report SHALL state the revised case-study interpretability statement for `angle` that the measured operating characteristics support and SHALL name the production-default change it selects; otherwise it SHALL state that the Phase 6 interpretability statement and the production default are unchanged.

#### Scenario: Report is tied to its run

- **WHEN** the findings report is written
- **THEN** it names both configurations, the code revision, the software versions, the record and failure counts, the committed results directories, and the commands that reproduce the merge, the report and the comparison selection

#### Scenario: Report names the branch taken

- **WHEN** the comparison selection has been computed
- **THEN** the report quotes its verdict, states the interpretability consequence for the case study, and the roadmap entry for item 2 is updated to the branch taken

### Requirement: Committed configurations pin the PC1 orientation estimator

Every committed study configuration whose recorded results were produced under the PC1 estimator SHALL name `"pc1"` as `evaluation.orientation_estimator` in its baseline evaluation block, and every committed configuration that declares `evaluation.orientation_estimator` as a design axis SHALL name its PC1 column `"pc1"` rather than `null`, so that changing the evaluator's default alters none of its enumerated cells' parameter signatures, cell identifiers, matched seeds, or persisted records. Every committed diagnostic that constructs evaluation parameters directly, and whose committed results report `angle` under PC1, SHALL likewise select `"pc1"` explicitly. The pins record historical behavior and MUST be documented as such, so that new configurations choose an estimator deliberately rather than copy them.

#### Scenario: Committed config identity is unchanged by the default flip

- **WHEN** each committed study configuration is loaded and enumerated after the default changes
- **THEN** its set of cell identifiers, every cell's parameter signature, and every cell's matched seeds equal those recorded for that configuration before the change

#### Scenario: Committed design-grid baseline column stays PC1

- **WHEN** a committed net-displacement profile is loaded after the default changes
- **THEN** its baseline column is the `"pc1"` point, its phases and cell identifiers are unchanged, and its net-displacement columns are unchanged

#### Scenario: A new config that omits the estimator gets the default

- **WHEN** a configuration that names no `evaluation.orientation_estimator` is loaded
- **THEN** every cell evaluates with the net-displacement estimator and its records state that estimator

### Requirement: Committed configurations pin the unfloored rank rule

Every committed study configuration whose recorded results predate the rank floor SHALL name `1` as `evaluation.min_components` in its baseline evaluation block, and every committed diagnostic that constructs evaluation or analysis parameters directly and whose results predate the floor SHALL pin the floor of 1, so the default alters none of their cell identifiers, signatures, matched seeds or records. The pins record historical behavior and MUST be documented as such.

#### Scenario: Committed config identity is unchanged by the default floor

- **WHEN** each committed study configuration is loaded and enumerated after the default floor exists
- **THEN** its set of cell identifiers, every cell's parameter signature and every cell's matched seeds equal those recorded for that configuration before the change

#### Scenario: A committed config names the floor deliberately

- **WHEN** a committed configuration under `examples/trajectory_power_study/` is read
- **THEN** its evaluation block names `min_components`, and the pin is documented as historical

#### Scenario: A new config that omits the floor gets the default

- **WHEN** a configuration that names no `evaluation.min_components` is loaded
- **THEN** every PLS cell selects its rank from candidates at or above 2 and its records state the floor

### Requirement: Rank floor is a design axis

`evaluation.min_components` SHALL be accepted as a design-grid axis (in either form) and as a one-factor-at-a-time axis. In a configuration, `null` SHALL name the evaluator's default and `1` the unfloored rule, as distinct values. Design columns that differ only in the floor SHALL share generator parameters and matched seeds. The floor-1 column's records SHALL carry no floor entry and SHALL equal the records of the same configuration without the axis.

#### Scenario: Floor axis is applied per design point

- **WHEN** a configuration whose baseline evaluation names a floor of 1 declares the design axis `evaluation.min_components` with values `[1, 2]`
- **THEN** the floor-2 column's records carry `min_components: 2` in their runtime metadata and select ranks of at least 2 on the same generated datasets as the baseline column
- **AND** the baseline column's records carry no floor entry and equal the records of the same configuration without the axis

#### Scenario: Null names the evaluator's default

- **WHEN** a configuration whose baseline evaluation omits the floor declares the design axis `evaluation.min_components` with values `[null, 1]`
- **THEN** the `null` point is the baseline column and evaluates at a floor of 2, and the `1` column's records carry no floor entry
