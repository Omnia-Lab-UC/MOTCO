# candidate-simulation-validation Specification

## Purpose
Reads the measurement candidates (`predictive_orthonormal`, block weighting and sparse stage-PLS) against the production measurement on simulated data only, across sample size, stage effect, stage count, balance and stage dispersion. A predeclared decision, scoped by design family, states where each candidate is better than production and where it is worse.

## Requirements

### Requirement: Every arm reads one dataset at one seed

Each unit SHALL generate one dataset and evaluate every arm on it at the unit's RRPP seed with serial draws. `pls_cv` SHALL be the production measurement: stage-PLS with net displacement and rank floor 2. `predictive_orthonormal` SHALL be measured on the `pls_cv` fit without a second fit. `weighted` and `sparse_block` SHALL each be their own cross-validated fit. A failure SHALL be recorded for the unit without aborting the shard.

#### Scenario: The predictive arm shares production's fit
- **WHEN** a unit runs
- **THEN** the `predictive_orthonormal` arm reports the selected rank and held-out AUROC of the unit's `pls_cv` arm, and no cross-validation is run for it

#### Scenario: Arms share draws
- **WHEN** two arms of one unit are evaluated
- **THEN** both use the unit's RRPP seed, so their permutation indices are equal

#### Scenario: Failed unit
- **WHEN** generation, fitting or an arm fails for a unit
- **THEN** the failure is recorded with its error for that unit and arm, and the shard continues

### Requirement: The sparse arm runs on a named subset

`sparse_block` SHALL run on the designs `n200`, `n200_delta0.5` and `band_n320` only. On every other design its readings SHALL be recorded as not run, and no rate SHALL be reported for it there.

#### Scenario: A design outside the subset
- **WHEN** a unit of a design outside the subset runs
- **THEN** no sparse fit is made, and the merged tables mark `sparse_block` as not run for that design

### Requirement: Committed designs are rerun and must reproduce

The `committed`, `surgery` and `adni_design` profiles SHALL rerun committed designs with their committed seeds, cells and replicate counts. For every cell and statistic the `pls_cv` available and rejected counts SHALL equal the committed rate table exactly, and on `adni_design` the `predictive_orthonormal` counts SHALL too. The committed settings and tables SHALL be bound by sha256 and refused by name if changed.

#### Scenario: Reproduction holds
- **WHEN** every compared row of a profile matches its committed counts
- **THEN** the merge records the reproduction as passed with the number of rows compared

#### Scenario: Reproduction fails
- **WHEN** any compared row differs
- **THEN** the merge records the reproduction as failed and lists every differing row with both counts

#### Scenario: A changed parent is refused
- **WHEN** a bound parent file's sha256 differs from the one in the stored settings
- **THEN** the driver exits naming the file and runs no unit

### Requirement: New conditions use fresh datasets

The `conditions` profile SHALL hold five three-block designs with the univariate baseline's thirteen stage-span cells: four stages at n = 200 with δ = 2.0 and with δ = 0.5; the SEA-AD cell table with δ = 2.0 and with δ = 0.5 at three stages; and a first-group share of 0.8 at n = 200, δ = 2.0 and three stages. Controls SHALL run at 500 replicates and other cells at 200, at 999 permutations.

#### Scenario: Unit count
- **WHEN** the profile enumerates its units
- **THEN** it has 16,000 units

#### Scenario: New datasets
- **WHEN** a design's dataset seeds are derived
- **THEN** none equals a dataset seed of a committed run, and cells of one design and replicate share their null dataset

### Requirement: The dispersion profile is null only

The `dispersion` profile SHALL hold ten three-block designs that differ from an equal-dispersion design only in `stage_dispersion_scale`: the patterns `last_narrow` and `last_wide` and their stronger forms on the SEA-AD cell table at δ = 2.0 and at δ = 0.5, and the two stronger forms at balanced n = 200, δ = 2.0. Each design SHALL have the anchor and the translation at 500 replicates and 999 permutations.

#### Scenario: Unit count
- **WHEN** the profile enumerates its units
- **THEN** it has 10,000 units, none of them a plant of a trajectory difference

#### Scenario: Only the noise differs
- **WHEN** a dispersion design's generator parameters are compared with its equal-dispersion design's
- **THEN** they differ only in `stage_dispersion_scale`

#### Scenario: Within-stage draws are recorded
- **WHEN** a dispersion unit runs
- **THEN** `pls_cv` and `predictive_orthonormal` are also read under within-stage draws on the same fit, statistics and seed, and those readings do not enter the decision

### Requirement: The decision is scoped by design family

The decision SHALL group the designs into the families `core`, `weak_stage`, `four_stages`, `unbalanced` and `adni_design`, and SHALL label each candidate arm in each family as `size_excess`, `power_loss`, `specificity_excess`, `no_gain` or `qualifies`, the first that applies, or `not_run`. Every criterion SHALL be recorded for every arm and family whatever the label.

#### Scenario: Size
- **WHEN** an arm's controls in a family are judged
- **THEN** its `delta`, `angle` and `shape` counts at p < 0.05 on the anchors and translations are judged under the count-plus-ceiling rule, and a failure labels the family `size_excess`

#### Scenario: Protected power
- **WHEN** an arm's in-span target cells in a family are compared with `pls_cv` on shared units
- **THEN** a one-sided exact paired loss test is applied per cell under the count-plus-ceiling rule, and a failure labels the family `power_loss`

#### Scenario: Specificity at the core stage effect
- **WHEN** a family's design has δ = 2.0
- **THEN** the arm's off-target rates on its in-span target cells are judged under the count-plus-ceiling rule

#### Scenario: Specificity below the core stage effect
- **WHEN** a family's design has δ below 2.0
- **THEN** the arm's off-target rejections are judged by a one-sided exact paired excess test against `pls_cv` under the count-plus-ceiling rule, and production's absolute off-target rates are reported beside it

#### Scenario: Gain
- **WHEN** an arm holds size, protected power and specificity in a family
- **THEN** it is labelled `qualifies` if at least one target cell has a one-sided exact paired gain p-value at or below 0.05 divided by the family's target cells, and `no_gain` otherwise

#### Scenario: An arm absent from a family
- **WHEN** an arm has no reading in a family
- **THEN** the family is labelled `not_run` for that arm

### Requirement: Size under stage dispersion is judged against production's own reading

The decision SHALL judge each arm's controls on the dispersion designs under the count-plus-ceiling rule. If `pls_cv` itself fails that rule, the criterion for every candidate SHALL instead be a one-sided exact paired excess test against `pls_cv` under the count-plus-ceiling rule, and production's failure SHALL be reported as a finding of its own.

#### Scenario: Production holds
- **WHEN** `pls_cv` holds size on the dispersion designs
- **THEN** a candidate whose counts fail the rule is recorded as failing the dispersion criterion

#### Scenario: Production fails
- **WHEN** `pls_cv` fails size on the dispersion designs
- **THEN** the decision records `production_dispersion_excess`, and a candidate fails the criterion only if it rejects in significant excess of `pls_cv`

### Requirement: Incidental routes are labelled and do not decide

On the feature-set surgery cells and the out-of-span rotation cells, the decision SHALL label each candidate against `pls_cv` per target cell as `candidate_loses`, `candidate_gains` or `unresolved`, by one-sided exact paired tests with the familywise error held at 0.05 and a margin of 0.10. These labels SHALL be stated in the arm's consequence and SHALL NOT change a family label or a verdict.

#### Scenario: A loss on a surgery cell
- **WHEN** a candidate rejects a surgery's target statistic significantly less often than `pls_cv`, by more than the margin
- **THEN** the cell is labelled `candidate_loses` and named in the arm's consequence

#### Scenario: Labels are outside the verdict
- **WHEN** the decision is computed with the surgery and out-of-span rates changed
- **THEN** every family label and verdict is unchanged

### Requirement: Each arm and the run receive a verdict

Each candidate arm SHALL receive `qualifies_everywhere` when every family is run and qualifies and the dispersion criterion holds, `scoped_candidate` when at least one family qualifies, and `keep_production` otherwise. The run SHALL receive `propose_default_change` naming one arm when an arm qualifies everywhere, `scoped_candidates` when an arm is a scoped candidate, and `keep_production` otherwise. A failed reproduction SHALL return `reproduction_failed` and no other verdict.

#### Scenario: An arm with an unread family
- **WHEN** an arm qualifies in every family it was run in and one family is `not_run`
- **THEN** its verdict is `scoped_candidate`

#### Scenario: Two arms qualify everywhere
- **WHEN** more than one arm qualifies everywhere
- **THEN** the run names the arm with the fewest off-target readings above the bound, ties taken in the order `predictive_orthonormal`, `weighted`, `sparse_block`

#### Scenario: Reproduction gates the verdict
- **WHEN** any reproduction failed
- **THEN** the verdict is `reproduction_failed`, and the criteria are still recorded

#### Scenario: The rule is fixed before the run
- **WHEN** the decision is written
- **THEN** it carries its rule version, and the rule's text was committed before any array unit ran

### Requirement: The report states the map and its limits

A dated report SHALL state the verdict, each reproduction, and each arm's label and criteria by family, with the dispersion criterion and the incidental-route labels. It SHALL state that no real-noise reading enters the decision and cite the SEA-AD and ADNI readings as recorded limits. Every number it quotes SHALL be readable from the committed aggregates.

#### Scenario: Descriptive readings
- **WHEN** the report is written
- **THEN** it gives the recorded within-stage readings, and the selected rank, held-out AUROC and span capture by design

#### Scenario: Reader asks about many features per sample
- **WHEN** the reader asks whether a candidate holds size at many features per sample
- **THEN** the report states that no design has more than 8 features per sample and that the question is unread

#### Scenario: Reader asks what the real data said
- **WHEN** the reader looks for the real-noise readings of an arm
- **THEN** the report names each one with its verdict and states that it did not enter this decision

#### Scenario: Reader checks a number
- **WHEN** the reader takes a rate quoted in the report
- **THEN** it can be read from a committed table under the run's results directory
