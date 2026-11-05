# Spec Delta

## ADDED Requirements

### Requirement: Stage-span orientation accepts an alignment with the stage span

The generator SHALL accept `span_alignment`, an angle φ in degrees in [0, 90], for the stage-span `orientation`. The rotation SHALL be made in the plane of u and w_φ = cos φ·v + sin φ·o, where u and v are the existing in-span directions and o is a unit direction on features whose null stage means do not vary. φ = 0 SHALL normalize to absent.

#### Scenario: Intended geometry does not depend on alignment
- **WHEN** a stage-span orientation dataset is generated at effect e and any φ in [0, 90]
- **THEN** its truth records an intended `angle` of 90°·e and an intended `delta` and `shape` of 0, in the null-standardized coordinates

#### Scenario: Zero alignment is the in-span rotation
- **WHEN** a dataset is generated with `span_alignment` 0
- **THEN** its samples, truth, cell identifier and signature are identical to those generated without the setting at the same parameters and seed

#### Scenario: Full alignment leaves the pooled stage span
- **WHEN** a stage-span orientation is generated at φ = 90 with balanced cells
- **THEN** in the null-standardized coordinates the population B − A difference at every stage, minus its mean over stages, is orthogonal to the span of the pooled population stage means, to numerical tolerance

#### Scenario: Alignment is rejected outside its scope
- **WHEN** `span_alignment` is set outside [0, 90], or with a construction other than `stage_span`, or with a trajectory mode other than `orientation` or `none`
- **THEN** the generator raises a validation error naming the setting

### Requirement: The out-of-span direction has a declared support

The generator SHALL accept `span_support`, absent by default, for a stage-span orientation with φ > 0. Absent, o SHALL be a random unit direction over every span-layer feature whose null stage means do not vary. An integer k SHALL restrict o to k such features drawn at random. o SHALL have no component on any feature whose null stage means vary.

#### Scenario: Sparse support has exactly k features
- **WHEN** `span_support` is k
- **THEN** o has exactly k nonzero entries, all on features whose null stage means are constant across stages

#### Scenario: Too few unused features
- **WHEN** k exceeds the number of span-layer features whose null stage means do not vary
- **THEN** the generator raises a validation error naming k and the number available

#### Scenario: Support without alignment is rejected
- **WHEN** `span_support` is set while `span_alignment` is absent or 0, or is not a positive integer
- **THEN** the generator raises a validation error

### Requirement: The out-of-span direction leaves the null dataset unchanged

The out-of-span direction SHALL be drawn from a random stream separate from the generator's sampling stream and derived from the dataset seed. At a given seed the samples SHALL equal the null dataset's plus the plant. o SHALL be the same for every effect and every φ > 0 at that seed and support.

#### Scenario: Dataset is still the null dataset plus the plant
- **WHEN** a stage-span orientation dataset with φ > 0 and a `none` dataset are generated at the same seed and design
- **THEN** group A's samples are identical, and group B's differ by exactly the plant of their stage, to floating-point precision

#### Scenario: The direction is shared across cells
- **WHEN** two datasets are generated at the same seed and support, with different effects or different nonzero φ
- **THEN** their truth records the same out-of-span direction identifier and the same support

### Requirement: Alignment settings are recorded and hash-omitted

Truth SHALL record `span_alignment`, `span_support` and the support's size when the settings are present. While both are absent, every dataset, truth record, cell identifier, signature and matched seed SHALL equal those produced before the settings existed.

#### Scenario: Absent settings are byte-identical
- **WHEN** a dataset is generated without `span_alignment` and `span_support`
- **THEN** its samples, truth, cell identifier and signature are byte-identical to those produced before this change

#### Scenario: Truth names the alignment
- **WHEN** a dataset is generated with `span_alignment` 60 and `span_support` 20
- **THEN** its truth records φ = 60, a support of 20 features and the span layers they were drawn from
