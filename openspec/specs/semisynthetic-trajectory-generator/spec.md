# semisynthetic-trajectory-generator Specification

## Purpose

Generate MOTCO-ready semi-synthetic trajectory datasets with the numpy generator, sampling from the cached InterSIM reference (no R at runtime), using the clusters-as-stages assumption, reproducible group assignment, trajectory effect injection, and explicit truth metadata.

## Requirements

### Requirement: Groups are assigned reproducibly within stages
The generator SHALL assign comparison group labels reproducibly within each stage according to configured group balance. Group balance SHALL be configurable in exactly one of two ways. The proportional way uses a total sample count, optional per-stage proportions, and one group ratio applied in every stage. The explicit way uses a group × stage sample-size table, one sequence of per-stage counts per group, ordered as the group labels and the stages. When the explicit table is given, every group-stage cell SHALL contain exactly the stated number of samples. The total sample count and stage counts SHALL be derived from the table. The generator SHALL reject a configuration that sets the explicit table together with any proportional setting that differs from its default. When the explicit table is absent, generation SHALL be identical to the generator's behavior before the table existed, at every seed.

#### Scenario: Two groups are assigned within every stage
- **WHEN** each stage has enough samples for two groups and no explicit size table is given
- **THEN** the generator assigns group labels within each stage according to the configured group ratio

#### Scenario: Same seed gives same group labels
- **WHEN** the same generator parameters and seed are used twice
- **THEN** the generated group labels are identical

#### Scenario: Insufficient stage size is rejected
- **WHEN** any stage has too few samples to assign both comparison groups
- **THEN** the generator raises a clear validation error

#### Scenario: Explicit group × stage sizes are realized exactly
- **WHEN** the caller supplies the size table `((11, 10, 28), (9, 10, 12))` with three stages
- **THEN** the generated metadata contains exactly 11, 10, and 28 samples of the first group and 9, 10, and 12 samples of the second group at stages 0, 1, and 2, and 80 samples in total
- **AND** truth metadata records the realized group × stage sizes

#### Scenario: Explicit table conflicts with proportional settings
- **WHEN** the caller supplies the size table together with a non-default `n_samples`, `stage_sample_prop`, or `group_ratio`
- **THEN** the generator raises a clear validation error naming the conflicting settings

#### Scenario: Malformed explicit table is rejected
- **WHEN** the size table has a number of rows other than two, a row length other than the number of stages, or any cell below one sample
- **THEN** the generator raises a clear validation error

#### Scenario: Absent table is byte-identical
- **WHEN** generator parameters omit the size table
- **THEN** the generated dataset, truth metadata, and parameter signature equal those produced before the table existed, at the same seed

### Requirement: Generator creates MOTCO-ready trajectory datasets from the numpy generator
MOTCO SHALL provide a semi-synthetic trajectory generator that builds datasets on top of the numpy omics generator and returns aligned omics matrices, sample metadata, and truth metadata, without an `InterSIMResult` or an R subprocess at runtime.

#### Scenario: Successful generation from the numpy generator
- **WHEN** a caller provides valid generator parameters
- **THEN** the generator returns methylation, gene expression, and proteomics matrices aligned to sample metadata, produced without invoking R

#### Scenario: Sample metadata contains trajectory design columns
- **WHEN** generation succeeds
- **THEN** sample metadata contains `sample_id`, `group`, `stage`, and `cluster` columns

#### Scenario: Truth metadata records generator parameters
- **WHEN** generation succeeds
- **THEN** truth metadata records trajectory mode, group effect size, per-stage/group differential indicators, per-omic effect sizes, and the group-assignment seed

### Requirement: Trajectory modes are feature-set surgery on methylation differential indicators
Under the default `surgery` construction, the generator SHALL define `none`, `translation`, `magnitude`, `orientation`, and `shape` as operations on group B's **per-stage methylation** differential-feature indicators only. Group A inherits a baseline set of per-stage methylation indicators (which need not form a continuous/straight trajectory). For both groups, the gene-expression and proteomics differential indicators SHALL be **derived from the (group-specific) methylation indicators** through the cached CpG→gene→protein incidence maps — the surgery never touches expression, proteomics, or the latent space directly. This keeps the simulated differences biologically grounded (methylation drives expression drives protein) and keeps the data realistic rather than tailored to MOTCO. The one surgery that acts on effect sizes rather than indicators is `magnitude`, whose per-omic δ scaling is selected by `magnitude_kind`. The opt-in `stage_span` construction defines the same four effect modes differently and is specified separately.

#### Scenario: Group B expression and protein indicators are derived from its methylation
- **WHEN** any non-null mode transforms group B's methylation indicators
- **THEN** group B's expression and proteomics differential indicators are re-derived from group B's methylation indicators via the incidence maps (not manipulated independently)

#### Scenario: None mode gives identical group trajectories
- **WHEN** `trajectory_mode` is `none` (or `group_effect_size` is 0)
- **THEN** group B uses the same methylation indicators and effects as group A, so the groups share an identical trajectory

#### Scenario: Translation mode adds an extra constant differential set
- **WHEN** `trajectory_mode` is `translation`
- **THEN** group B keeps group A's stage-changing methylation sites unchanged and additionally marks an extra set `U` of methylation sites — whose mapped genes are absent from the stage program — as differential at every group-B stage (and at none of group A's), producing a constant group offset that leaves the size, orientation, and shape statistics unchanged

#### Scenario: Magnitude mode scales the methylation effect
- **WHEN** `trajectory_mode` is `magnitude` and `magnitude_kind` is `all` or `extremes`
- **THEN** group B uses the same per-stage methylation indicators as group A but with a scaled methylation effect size `δ_methyl_B = (1 + e)·δ_methyl` (at every stage for `all`, through the endpoint indicators for `extremes`), leaving `δ_expr` and `δ_protein` at group A's values — a size/`delta` change that is exactly size-pure within each omic block but not in the joint space

#### Scenario: Joint magnitude mode scales every omic's effect together
- **WHEN** `trajectory_mode` is `magnitude` and `magnitude_kind` is `joint`
- **THEN** group B uses the same per-stage methylation indicators as group A with every per-omic effect size scaled by the same factor, `δ_methyl_B = (1 + e)·δ_methyl`, `δ_expr_B = (1 + e)·δ_expr`, `δ_protein_B = (1 + e)·δ_protein`, so group B's native-space stage-mean trajectory is exactly `(1 + e)` times group A's in every omic block — a size/`delta` change that is size-pure in the joint space as well

#### Scenario: Orientation mode relocates stage-changing sites consistently across stages
- **WHEN** `trajectory_mode` is `orientation`
- **THEN** a fraction `e` of group A's stage-changing methylation sites are relocated to different CpGs using a single relocation applied identically to every stage, so the per-stage on/off pattern is preserved on different feature axes (a rotation: orientation changes, with size and shape preserved in the linear limit)

#### Scenario: Shape mode perturbs a single interior stage
- **WHEN** `trajectory_mode` is `shape` and at least three stages are available
- **THEN** group B perturbs a single interior stage relative to group A — either by relocating a fraction `e` of that stage's methylation sites (`relocate`) or by scaling that stage's methylation effect (`magnitude`) — bending one interior vertex of the trajectory (a shape change, which may co-move size)

#### Scenario: Shape mode rejects fewer than three stages
- **WHEN** `trajectory_mode` is `shape` and fewer than three stages are available
- **THEN** the generator raises a clear validation error

### Requirement: Baseline indicators and the group transform are deterministic
The generator SHALL construct group A's baseline per-stage methylation indicators and group B's transform deterministically from the seed and parameters, recording both in truth metadata.

#### Scenario: Same seed gives the same indicators
- **WHEN** the same parameters and seed are used twice
- **THEN** the per-stage/group differential indicators are identical across runs

#### Scenario: Indicators are recorded as truth
- **WHEN** generation succeeds
- **THEN** truth metadata records the per-stage methylation (and derived expression/protein) indicators for both groups, the per-omic effect sizes, and any extra translation set

### Requirement: Generator provides a single-call convenience path
MOTCO SHALL provide a convenience API that generates a semi-synthetic trajectory dataset end to end from parameters, using the numpy generator and cached reference data.

#### Scenario: End-to-end generation succeeds without R
- **WHEN** a caller provides valid parameters
- **THEN** the convenience API returns a semi-synthetic trajectory dataset using only cached reference data, with no R invocation

### Requirement: Generator emits per-stage/group differential-indicator truth for characterization
The generator SHALL expose the differential-feature indicators for each stage and group so a consumer can characterize how MOTCO responds to an injected mode (which statistics move, and how much cross-talk there is).

#### Scenario: Indicator truth is consumable downstream
- **WHEN** a dataset is generated
- **THEN** the per-stage/group indicators are available in the dataset's truth structure for the showcase and study to characterize mode-to-statistic response (a descriptive specificity matrix, not a pass/fail gate)

### Requirement: Magnitude kind selects the construction, joint by default
The generator SHALL provide a `magnitude_kind` option with three selectable values: `joint` (the default) scales **every omic's** effect size (methylation, expression and proteomics) by the same factor `1 + e` at all stages; `all` scales group B's methylation effect alone at **all** stages; `extremes` scales the methylation effect only at the **extreme** stages (the first and last). `all` and `extremes` SHALL reproduce their existing methylation-only behavior byte for byte at every seed when selected explicitly. An unknown value SHALL be rejected with a clear validation error naming the allowed values. Every value SHALL be selectable through the generator parameters, the study configuration loader, and the `motco simulate` command line, and the command line's default SHALL be the generator's default. `extremes` and the shape mode's `magnitude` kind remain methylation-only scalings; the joint default does not extend to them.

#### Scenario: Joint magnitude is the default
- **WHEN** `trajectory_mode` is `magnitude` and `magnitude_kind` is unset
- **THEN** group B's methylation, expression, and proteomics effect sizes are each `(1 + e)` times group A's, and truth metadata records `magnitude_kind` as `joint`

#### Scenario: All-stages methylation-only magnitude is selectable
- **WHEN** `trajectory_mode` is `magnitude` and `magnitude_kind` is `all`
- **THEN** group B's methylation effect is scaled at every stage with expression and proteomics at group A's values, and the dataset and truth are identical to those the generator produced for an unset `magnitude_kind` before the default changed, at the same seed

#### Scenario: Extreme-stage magnitude scales only the endpoints
- **WHEN** `trajectory_mode` is `magnitude` and `magnitude_kind` is `extremes`
- **THEN** group B's methylation effect is scaled only at the first and last stages, leaving interior stages at the baseline effect

#### Scenario: Joint magnitude scales all three per-omic effects
- **WHEN** `trajectory_mode` is `magnitude` and `magnitude_kind` is `joint` with effect size `e`
- **THEN** group B's methylation, expression, and proteomics effect sizes are each `(1 + e)` times group A's, and group B's methylation indicators are identical to group A's

#### Scenario: Joint magnitude is exactly null at zero effect
- **WHEN** `trajectory_mode` is `magnitude`, `magnitude_kind` is `joint`, and `group_effect_size` is 0
- **THEN** the generated dataset and truth are identical to the `none` mode's at the same seed

#### Scenario: Joint magnitude is selectable from a study configuration
- **WHEN** a study configuration sets `generator.magnitude_kind` to `joint`
- **THEN** the configuration loads, enumerates, and generates with the joint construction, and a configuration naming any value outside `all`, `extremes`, `joint` is still refused at load time

#### Scenario: Command-line default follows the generator
- **WHEN** `motco simulate` is run in `magnitude` mode without `--magnitude-kind`
- **THEN** the written truth records `magnitude_kind` as `joint`, and passing `--magnitude-kind all` reproduces the methylation-only construction

#### Scenario: Magnitude variant is recorded as truth
- **WHEN** a `magnitude` dataset is generated
- **THEN** truth metadata records which `magnitude_kind` was used and, for `joint`, the common scale factor applied to each of the three per-omic effect sizes

### Requirement: Historical configurations pin the pre-default magnitude construction
Every committed study configuration that predates the `joint` default and does not already name a magnitude kind SHALL name `magnitude_kind` `all` explicitly in its generator block, so that changing the default alters none of its enumerated cells' parameter signatures, cell identifiers, matched seeds, or generated datasets. Every committed diagnostic that constructs generator parameters directly, and whose committed results were produced under the former default, SHALL likewise select `all` explicitly. The pins record historical behavior and MUST be documented as such, so that new configurations do not copy them.

#### Scenario: Historical config signatures are unchanged by the default flip
- **WHEN** each historical study configuration is loaded and enumerated after the default changes
- **THEN** every enumerated cell's parameter signature equals the signature recorded for that cell before the change

#### Scenario: Historical diagnostics regenerate their committed data
- **WHEN** a pinned historical diagnostic is re-run with its committed command after the default changes
- **THEN** it generates the same datasets it generated before the change

#### Scenario: New configuration without a kind gets the joint construction
- **WHEN** a study configuration that does not set `generator.magnitude_kind` requests the magnitude mode
- **THEN** it resolves to the `joint` construction

### Requirement: Pool-limited surgeries apply an explicit censoring policy
The generator SHALL apply an explicit, configurable censoring policy to every surgery whose size is
limited by a destination or candidate pool (`orientation` relocation, `translation` extra set,
`shape` with `shape_kind='relocate'`). The policy SHALL default to failing loudly: when the
requested surgery size exceeds the available pool, generation raises a validation error naming the
mode, the requested size, and the available pool size. An explicit opt-in policy value SHALL
preserve the previous clamping behavior (realize as much of the surgery as the pool allows). The
generator MUST NOT silently clamp under the default policy.

#### Scenario: Default policy fails loudly on a censored surgery
- **WHEN** a pool-limited surgery's requested size exceeds its available pool and no censoring
  policy is explicitly configured
- **THEN** generation raises a clear validation error identifying the trajectory mode, the
  requested surgery size, and the available pool size, and no dataset is returned

#### Scenario: Opt-in clamping preserves the previous behavior
- **WHEN** the censoring policy is explicitly set to clamp and a pool-limited surgery's requested
  size exceeds its available pool
- **THEN** the generator realizes the largest surgery the pool allows and generation succeeds,
  matching the pre-policy behavior replicate for replicate at the same seed

#### Scenario: Uncensored surgeries are unaffected by the policy
- **WHEN** a pool-limited surgery's requested size does not exceed its available pool
- **THEN** the realized surgery equals the requested size under either policy value, and the
  sampled dataset at a given seed is identical to the pre-policy generator's output

### Requirement: Truth metadata makes requested-vs-realized surgery explicit
For every pool-limited surgery, truth metadata SHALL record the nominal (requested) surgery size
beside the realized surgery size, and a flag stating whether the surgery was censored (realized
size smaller than nominal). This makes the requested-effect → realized-surgery relationship
auditable per record, without regenerating data.

#### Scenario: Nominal and realized sizes are recorded together
- **WHEN** a dataset with a pool-limited surgery is generated successfully
- **THEN** truth metadata contains the nominal surgery size, the realized surgery size, and a
  censored flag that is true exactly when realized < nominal

#### Scenario: Clamped generation is identifiable downstream
- **WHEN** a dataset is generated under the clamping policy and the clamp binds
- **THEN** the truth metadata's censored flag is true, so downstream summaries can identify the
  record without access to the generator

### Requirement: Baseline stage-program continuity is a generator axis
The generator SHALL expose a baseline continuity parameter ρ ∈ [0, 1) controlling how strongly
group A's per-stage methylation differential indicators persist across adjacent stages. Each CpG's
indicator sequence across stages SHALL follow a stationary first-order Markov chain whose per-stage
marginal is Bernoulli(`p_dmp`) at every ρ (first stage ~ Bernoulli(`p_dmp`); stay-differential
probability `p_dmp + ρ·(1 − p_dmp)`; become-differential probability `p_dmp·(1 − ρ)`), so that the
cross-stage indicator correlation is ρ^|t−s| while per-stage indicator counts, the per-omic effect
sizes, the cross-omic coupling derivation, and the meaning of `p_dmp` are unchanged along the axis.
At ρ = 0 the generator SHALL reproduce the independent per-stage draws of the pre-change generator
byte-identically at the same seed. Values outside [0, 1) SHALL be rejected with a clear validation
error, and truth metadata SHALL record the continuity value.

#### Scenario: Marginals are preserved along the axis
- **WHEN** datasets are generated at any ρ ∈ [0, 1)
- **THEN** each stage's methylation indicators are marginally Bernoulli(`p_dmp`), and the expected
  per-stage differential count does not depend on ρ

#### Scenario: Zero continuity reproduces the independent baseline exactly
- **WHEN** a dataset is generated with ρ = 0 and any seed
- **THEN** the baseline indicators, all derived indicators, and the sampled dataset are identical
  to the pre-change generator's output at that seed

#### Scenario: Continuity induces trending stage configurations
- **WHEN** ρ > 0
- **THEN** adjacent stages share differential programs with correlation ρ, so expected squared
  stage-mean distances grow with stage separation (proportional to 1 − ρ^|t−s|) instead of being
  equal for all pairs

#### Scenario: Continuity is recorded as truth
- **WHEN** generation succeeds
- **THEN** truth metadata records the baseline continuity value beside the other generator
  parameters

#### Scenario: Invalid continuity is rejected
- **WHEN** a caller passes a continuity value below 0 or at/above 1
- **THEN** generation raises a validation error naming the parameter and the allowed range

### Requirement: Expected surgery headroom accounts for baseline continuity
The generator's analytic expected-headroom computation for pool-limited surgeries SHALL use the
continuity-adjusted probability that a CpG is differential in at least one stage,
`1 − (1 − p_dmp)·(1 − p_dmp·(1 − ρ))^(n_stages − 1)`, which reduces to the independence union
probability at ρ = 0. Expected destination pools and saturating effects SHALL therefore reflect
that higher continuity shrinks the expected stage-active union (stage programs overlap more).

#### Scenario: Headroom at zero continuity is unchanged
- **WHEN** expected headroom is computed at ρ = 0
- **THEN** the expected pools and saturating effects equal the pre-change independence-based values

#### Scenario: Headroom tracks the continuity-adjusted union
- **WHEN** expected headroom is computed at ρ > 0 for an orientation or shape-relocate surgery
- **THEN** the expected stage-active fraction uses the continuity-adjusted union probability, and
  the reported saturating effect is correspondingly no smaller than at ρ = 0 for the same
  `p_dmp` and `n_stages`

### Requirement: Construction is an opt-in generator setting

The generator SHALL accept a `construction` setting, `surgery` (the default) or `stage_span`, and a `span_layers` setting naming the blocks the stage span is computed over. Both SHALL be absent by default, and `surgery` SHALL normalize to absent. While both are absent, every dataset, truth record, cell identifier, signature and matched seed SHALL equal those produced before the settings existed.

#### Scenario: Default construction is byte-identical
- **WHEN** a dataset is generated without `construction` and `span_layers`, or with `construction` set to `surgery`
- **THEN** it is byte-identical to the dataset generated before this change at the same parameters and seed

#### Scenario: Span layers are validated
- **WHEN** `construction` is `stage_span` and `span_layers` names an unknown, duplicated or empty set of blocks
- **THEN** the generator raises a validation error naming the offending value
- **AND** when `span_layers` is absent under `stage_span`, the span is computed over all three blocks

#### Scenario: Span layers without the stage-span construction are rejected
- **WHEN** `span_layers` is set while `construction` is absent or `surgery`
- **THEN** the generator raises a validation error

### Requirement: Stage-span edits are made in pooled-standardized joint coordinates

Under `stage_span`, the plant SHALL be defined on group B's exact null population stage means over `span_layers`, in analysed units (methylation M-values), each feature divided by its null population pooled standard deviation. That deviation includes the within-cell variance and the between-cell variance of the null means at the realized cell proportions. The plant SHALL be mapped back to analysed units and added to B's cell means before sampling. Blocks outside `span_layers` SHALL be unchanged.

#### Scenario: Plant is applied in analysed units before sampling
- **WHEN** a stage-span dataset is generated
- **THEN** group B's methylation cell means are shifted in M-value space before the inverse logit, and its expression and proteomics cell means in their native units

#### Scenario: Blocks outside the span are untouched
- **WHEN** `span_layers` names methylation and expression only
- **THEN** group B's proteomics samples equal those of the null dataset at the same seed

### Requirement: Stage-span modes edit the configuration inside its own span

With z_s B's standardized null stage means, centroid c, path length L, mean step L̄ = L/(S − 1) and u the unit net displacement z_{S−1} − z_0, the four modes SHALL be the edits below. v is the unit leading principal direction of the configuration after its u-components are removed. Mode `none`, or an effect of 0, SHALL yield the null dataset.

#### Scenario: Stage-span magnitude scales the configuration
- **WHEN** `trajectory_mode` is `magnitude` under `stage_span` with effect e
- **THEN** B's stage means become c + (1 + e)(z_s − c), so in the null-standardized coordinates B's path length is (1 + e) times A's, and `angle` and `shape` against A are 0

#### Scenario: Stage-span orientation rotates the configuration in its own plane
- **WHEN** `trajectory_mode` is `orientation` under `stage_span` with effect e, 0 ≤ e ≤ 2
- **THEN** B's configuration is rotated rigidly about c by θ = 90°·e in the plane spanned by u and v, so `angle` against A is θ, and `delta` and `shape` are 0

#### Scenario: Stage-span shape bends one interior stage at fixed size and direction
- **WHEN** `trajectory_mode` is `shape` under `stage_span` with effect e
- **THEN** the middle interior stage (the stage the `shape` surgery perturbs) moves by e·L̄ along the unit direction in the configuration's span that is perpendicular to the chord between its neighbours and points away from it
- **AND** the configuration is then rescaled to path length L and translated back to centroid c, so `delta` and the net-displacement `angle` against A are 0 and `shape` is positive for e > 0

#### Scenario: Stage-span translation offsets every stage along the progression
- **WHEN** `trajectory_mode` is `translation` under `stage_span` with effect e
- **THEN** every B stage mean moves by e·L̄·u, and `delta`, `angle` and `shape` against A are 0

#### Scenario: Degenerate configurations are rejected
- **WHEN** `orientation` or `shape` is requested under `stage_span` with fewer than three stages, or with a configuration whose second singular value is at most 1e-8 times its first
- **THEN** the generator raises a validation error stating that the stage span has no second direction
- **AND** an orientation effect above 2 is rejected

### Requirement: Stage-span construction consumes no randomness

The stage-span construction SHALL draw no random numbers, so at a given seed the group A samples, the baseline indicators and group B's sampling noise equal those of the `none` dataset.

#### Scenario: Stage-span dataset is the null dataset plus the plant
- **WHEN** a stage-span dataset and a `none` dataset are generated at the same seed and design
- **THEN** group A's samples are identical
- **AND** group B's samples in analysed units differ from the null dataset's by exactly the plant of their stage, to floating-point precision

### Requirement: Stage-span truth records the intended geometry

Truth metadata for a stage-span dataset SHALL record the construction, the span layers, the effect parameters (scale, rotation angle, bent stage and displacement, or offset length) and the intended population `delta`, `angle` and `shape` in the null-standardized coordinates. It SHALL also record the plant's norm per stage. Population trajectories SHALL include the plant.

#### Scenario: Truth and population trajectories carry the plant
- **WHEN** a stage-span dataset is generated
- **THEN** its truth names `construction: "stage_span"`, the span layers, the effect parameters and the intended geometry
- **AND** its population trajectories for group B equal the null means plus the plant in every span block
