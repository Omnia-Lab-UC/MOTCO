# Spec Delta

## MODIFIED Requirements

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

## ADDED Requirements

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
