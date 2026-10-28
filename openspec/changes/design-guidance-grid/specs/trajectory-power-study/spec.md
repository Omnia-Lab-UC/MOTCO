# Spec Delta

## ADDED Requirements

### Requirement: Design-guidance core profiles

The study SHALL provide two version-controlled design-guidance core pilot configurations, a magnitude profile and a geometry profile. Each crosses `n_samples` ∈ {80, 200, 500, 1000}, `n_stages` ∈ {3, 4} and the measured blocks (methylation + expression, or all three) in one design grid around the baseline of 200 samples, three stages and two blocks. Each runs at 100 replicates and 199 RRPP permutations.

#### Scenario: Core profiles enumerate the crossed grid

- **WHEN** each core profile is loaded and enumerated
- **THEN** it yields sixteen design columns, the baseline column included, each with one zero-effect anchor and one cell per (trajectory mode, nonzero effect)
- **AND** the magnitude profile's only mode is `magnitude` at effects 0.01, 0.02, 0.05, 0.1 and 0.25
- **AND** the geometry profile's modes are `orientation`, `shape` and `translation` at effects 0.25, 0.5 and 1.0

#### Scenario: Block columns share data

- **WHEN** two core columns differ only in the measured blocks
- **THEN** at every replicate index they evaluate the same generated dataset

#### Scenario: Every core cell fits the surgery headroom

- **WHEN** a core profile is enumerated
- **THEN** enumeration completes without a headroom failure at every design column

### Requirement: Design-guidance arm profiles

The study SHALL provide two version-controlled design-guidance arm pilot configurations, a magnitude profile and a geometry profile, using the core profiles' modes, effects, replicates and permutations. Each declares explicit design points, and each point moves exactly one generator setting away from the core baseline.

#### Scenario: Arms move one setting at a time

- **WHEN** each arm profile is loaded and enumerated
- **THEN** besides the baseline column it yields one column for each of: group ratio 0.65 and 0.8; stage proportions 0.5/0.3/0.2; `p_dmp` 0.05 and 0.15; stage effect 0.25, 0.5 and 1.0; baseline continuity 0.5 and 0.8
- **AND** every such column differs from the baseline in that one setting only

#### Scenario: The stage effect moves jointly

- **WHEN** an arm column sets the stage effect
- **THEN** its methylation, expression and protein stage effects all equal that value

#### Scenario: Every arm cell fits the surgery headroom

- **WHEN** an arm profile is enumerated
- **THEN** enumeration completes without a headroom failure at every design column

### Requirement: Design-guidance profiles name the production measurement

Every design-guidance profile SHALL name in its baseline blocks:
- the integration method `pls`;
- the orientation estimator `"net_displacement"`;
- the rank floor `2`;
- the magnitude kind `"joint"`.

It SHALL disable attribution and SHALL declare an advisory anchor decision. A committed identity fixture SHALL pin its cell identifiers, parameter signatures and matched seeds.

#### Scenario: Settings are named, not defaulted

- **WHEN** a design-guidance profile is read
- **THEN** its evaluation block names `orientation_estimator: "net_displacement"` and `min_components: 2`, and its generator block names `magnitude_kind: "joint"`

#### Scenario: Profile identity is pinned

- **WHEN** each design-guidance profile is loaded and enumerated
- **THEN** its cell identifiers, every cell's parameter signature and every cell's matched seeds equal those in the committed identity fixture

#### Scenario: Records state the production measurement

- **WHEN** a design-guidance cell is evaluated
- **THEN** its record states the net-displacement estimator in its runtime metadata and the rank floor 2 in its integration parameters

### Requirement: Reporting summarizes stage separability by design point

When a merged result set contains design-point cells, study reporting SHALL write a stage-signal table beside the design-point operating table, with one row per (design point, trajectory mode, effect size). Each row SHALL carry:
- the number of records with a recorded cross-validated stage AUROC;
- that AUROC's mean, 33rd percentile, median and 67th percentile.

The design-point operating table SHALL be unchanged.

#### Scenario: Stage-signal table is produced

- **WHEN** a merged result set from a design-grid study is reported
- **THEN** a stage-signal table is written whose rows are identified by every design coordinate, the trajectory mode and the effect size
- **AND** the baseline column and every design point's zero-effect anchor have rows

#### Scenario: Columns without a cross-validated rank are reported

- **WHEN** a design column's records carry no cross-validated stage AUROC, as for a forced rank, `concat` or `pca`
- **THEN** that column's rows report zero available records and empty summaries, and reporting completes without error

#### Scenario: Existing outputs are unchanged

- **WHEN** the same merged records of a design-grid study are reported with and without this change
- **THEN** every output written without it is byte-identical to the corresponding output written with it, and the stage-signal table is the only addition

#### Scenario: Studies without a design grid write no stage-signal table

- **WHEN** a merged result set contains no design-point cells
- **THEN** no stage-signal table is written

### Requirement: Design-guidance figure

The study SHALL provide a script that reads the design-point operating tables of the two core profiles and writes:
- for magnitude (`delta`), orientation (`angle`) and shape (`shape`), a figure with one panel per (stage count, block set) that plots the target statistic's rejection rate against sample size, one line per effect size, with Monte Carlo error bars;
- a control figure;
- the plotted values as CSV.

#### Scenario: Mode figures are written

- **WHEN** the script is run on the reports of both core profiles
- **THEN** it writes one figure each for magnitude/`delta`, orientation/`angle` and shape/`shape`, each with four panels (three or four stages × two or three blocks), sample size on the horizontal axis and one line per enumerated nonzero effect

#### Scenario: Control figure shows size

- **WHEN** the script is run on the reports of both core profiles
- **THEN** it writes a control figure plotting each statistic's rejection rate against sample size at the zero-effect anchor and at translation's largest effect, with the bound α + 2·√(α(1 − α)/R) drawn for the anchor's replicate count R

#### Scenario: Plotted values are written as data

- **WHEN** the script writes its figures
- **THEN** it also writes a CSV whose rows are exactly the operating-table rows it plotted

#### Scenario: A report without the core axes is refused

- **WHEN** the script is given a report whose design-point operating table lacks the sample-size, stage-count or block axis
- **THEN** it exits non-zero with a message naming the missing axis and writes no figure

### Requirement: Dated design-guidance pilot report

The study SHALL produce a dated findings report for the design-guidance pilot. It SHALL be tied to the exact configurations, code revision, parameter signatures, software versions, record counts, failure counts and reproduction commands. It SHALL report:
- per-column operating characteristics with Monte Carlo uncertainty;
- each column's stage signal and selected rank;
- for each factor of the originally planned simulation grid, what was run instead, deferred or dropped, and why.

#### Scenario: Report maps the planned grid

- **WHEN** the report is read
- **THEN** it lists each factor of the originally planned grid with its disposition (run, deferred to a named follow-up, or dropped) and the reason. The factors are:
  - sample size, group ratio, number of stages and stage proportions;
  - differential-feature proportion, shared stage shift and noise level;
  - group-specific mean shift, covariance scaling and correlation scaling;
  - feature filtering and integration method.

#### Scenario: Report reads the pilot as a pilot

- **WHEN** the report states a rejection rate
- **THEN** the rate carries its Monte Carlo standard error at the pilot's replicate count
- **AND** the report proposes the slice to rerun at paper grade without making paper-grade claims

#### Scenario: Report states every column's size

- **WHEN** the report is read
- **THEN** it gives every design column's zero-effect anchor rate per statistic and its advisory anchor-decision status

## MODIFIED Requirements

### Requirement: Committed configurations pin the PC1 orientation estimator

Every committed study configuration whose recorded results were produced under the PC1 estimator SHALL name `"pc1"` as `evaluation.orientation_estimator` in its baseline evaluation block. Every committed configuration that declares `evaluation.orientation_estimator` as a design axis SHALL name its PC1 column `"pc1"` rather than `null`. Together these ensure that changing the evaluator's default alters none of its enumerated cells' parameter signatures, cell identifiers, matched seeds, or persisted records. Every committed diagnostic that constructs evaluation parameters directly, and whose committed results report `angle` under PC1, SHALL likewise select `"pc1"` explicitly. The pins record historical behavior and MUST be documented as such, so that new configurations choose an estimator deliberately rather than copy them. A configuration committed after the default flip SHALL name its estimator explicitly. It is not held to the pre-flip identity fixture, and its identity SHALL be pinned by a fixture of its own.

#### Scenario: Committed config identity is unchanged by the default flip

- **WHEN** each study configuration committed before the default changed is loaded and enumerated after the default changes
- **THEN** its set of cell identifiers, every cell's parameter signature, and every cell's matched seeds equal those recorded for that configuration before the change

#### Scenario: Committed design-grid baseline column stays PC1

- **WHEN** a committed net-displacement profile is loaded after the default changes
- **THEN** its baseline column is the `"pc1"` point, its phases and cell identifiers are unchanged, and its net-displacement columns are unchanged

#### Scenario: A new config that omits the estimator gets the default

- **WHEN** a configuration that names no `evaluation.orientation_estimator` is loaded
- **THEN** every cell evaluates with the net-displacement estimator and its records state that estimator

#### Scenario: Every committed configuration is accounted for

- **WHEN** the committed study configurations are listed
- **THEN** each one is in the pre-flip identity fixture or in an identity fixture for configurations committed since, and none is in neither

### Requirement: Committed configurations pin the unfloored rank rule

Every committed study configuration whose recorded results predate the rank floor SHALL name `1` as `evaluation.min_components` in its baseline evaluation block. Every committed diagnostic that constructs evaluation or analysis parameters directly, and whose results predate the floor, SHALL pin the floor of 1. Together these ensure the default alters none of their cell identifiers, signatures, matched seeds or records. The pins record historical behavior and MUST be documented as such. A configuration committed after the floor SHALL name `evaluation.min_components` explicitly, with the floor it measures under, and is not held to the historical pin.

#### Scenario: Committed config identity is unchanged by the default floor

- **WHEN** each study configuration committed before the floor existed is loaded and enumerated after the default floor exists
- **THEN** its set of cell identifiers, every cell's parameter signature and every cell's matched seeds equal those recorded for that configuration before the change

#### Scenario: A committed config names the floor deliberately

- **WHEN** a committed configuration under `examples/trajectory_power_study/` is read
- **THEN** its evaluation block names `min_components`
- **AND** if the configuration predates the floor, the value is `1` and the pin is documented as historical

#### Scenario: A post-floor config names its floor

- **WHEN** a design-guidance profile is read
- **THEN** its evaluation block names `min_components: 2`, and the historical floor-1 check does not apply to it

#### Scenario: A new config that omits the floor gets the default

- **WHEN** a configuration that names no `evaluation.min_components` is loaded
- **THEN** every PLS cell selects its rank from candidates at or above 2 and its records state the floor
