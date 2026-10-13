# Spec Delta

## MODIFIED Requirements

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

## ADDED Requirements

### Requirement: Explicit design points enumerate exactly like crossed design points

When a design grid is declared in the explicit form, enumeration SHALL treat every listed non-baseline point exactly as it treats a non-baseline point of a crossed grid: it emits that point's anchored power grid under the `power_design` phase, stamps the point's coordinates in cell metadata, places the cells in the primary matched-seed family, and applies the duplicate-dataset guard and the surgery-headroom check with the point's own parameters. The baseline point MUST be served by the primary cells and never re-emitted. A configuration whose explicit point list equals the full cross of some axes MUST enumerate the same cell identifiers and parameter signatures as the crossed declaration of those axes.

#### Scenario: Explicit points enumerate their own anchored grids

- **WHEN** a configuration lists the baseline point and four further design points
- **THEN** enumeration emits four `power_design` grids, each with one zero-effect anchor and one cell per (mode, nonzero effect), each carrying its point's coordinates

#### Scenario: Explicit and crossed forms agree on a rectangle

- **WHEN** one configuration declares axes `A` and `B` in the crossed form and another lists every (aᵢ, bⱼ) pair explicitly, with everything else equal
- **THEN** both enumerate identical cell identifiers, parameter signatures and replicate seeds

### Requirement: Orientation estimator and integration method are design axes

The evaluation-namespace axes `evaluation.orientation_estimator` and `evaluation.integration_method` SHALL be accepted as design-grid axes (in either form) and as one-factor-at-a-time axes. The baseline value of `evaluation.orientation_estimator` is `null` (the PC1 estimator) unless the baseline evaluation sets it, and `"pc1"` is the same value as `null`. Cells at a design point SHALL evaluate with the point's estimator and integration method, every RRPP draw and realized-geometry checkpoint included, and cells at the baseline column SHALL write records byte-identical to a configuration without these axes. The design-point operating table SHALL carry both axes as columns and SHALL accept columns whose integration method records no PLS component selection or selected rank, reporting those fields as unavailable rather than failing. The production forced-rank guard SHALL admit a forced-rank column declared through the explicit form on the same terms as one declared through the crossed form: its records pass when the recorded rank equals the rank at their design point.

#### Scenario: Estimator axis is applied per design point

- **WHEN** a configuration declares the design axis `evaluation.orientation_estimator` with values `[null, "net_displacement"]`
- **THEN** the net-displacement column's records carry `orientation_estimator: "net_displacement"` in their runtime metadata and their `angle` p-values come from that estimator
- **AND** the baseline column's records carry no estimator entry and equal the records of the same configuration without the axis

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
