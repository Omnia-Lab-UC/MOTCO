# Spec Delta

## ADDED Requirements

### Requirement: Design-guidance paper-grade slice profiles

The study SHALL provide two version-controlled paper-grade design-guidance profiles, a magnitude profile and a geometry profile. Each crosses `n_samples` ∈ {80, 200, 500}, three stages and the measured blocks (methylation + expression, or all three) in one design grid around the baseline of 200 samples, three stages and two blocks, runs 500 replicates at 999 RRPP permutations in the pilot's seed family, and satisfies "Design-guidance profiles name the production measurement".

#### Scenario: Slice profiles enumerate six columns

- **WHEN** each slice profile is loaded and enumerated
- **THEN** it yields six design columns, the baseline column included, each with one zero-effect anchor and one cell per (trajectory mode, nonzero effect)
- **AND** the magnitude profile's only mode is `magnitude` at effects 0.02, 0.05 and 0.1
- **AND** the geometry profile's modes are `orientation`, `shape` and `translation` at effects 0.25, 0.5 and 1.0
- **AND** every cell runs 500 replicates at 999 permutations

#### Scenario: Slice columns evaluate the pilot's datasets first

- **WHEN** a slice column cell in the primary or design-point phase is compared with the pilot core cell at the same sample size, three stages, the same blocks, mode and effect
- **THEN** their matched seeds at replicate indices 0 through 99 are identical
- **AND** the slice's Type I baseline cells are not held to this, because their seeds are keyed on their own identity

#### Scenario: Block columns share data

- **WHEN** two slice columns differ only in the measured blocks
- **THEN** at every replicate index they evaluate the same generated dataset

#### Scenario: Slice profiles predeclare the size rules only

- **WHEN** a slice profile is read
- **THEN** its acceptance block declares the advisory anchor decision at α 0.05 and multiplier 2 over `delta`, `angle` and `shape`, a Type I entry without a mode, and a column-scoped Type I entry naming `translation` and `angle` under the upper bound at α 0.05 and tolerance 2
- **AND** it declares no power target, no specificity target, no gate and no report contract

#### Scenario: Every slice cell fits the surgery headroom

- **WHEN** a slice profile is enumerated
- **THEN** enumeration completes without a headroom failure at every design column

### Requirement: Design-guidance stage-effect profiles

The study SHALL provide two version-controlled stage-effect profiles, a magnitude profile and a geometry profile, with the slice profiles' modes, effects, replicates, permutations, seed family and measurement. Each declares explicit design points that set the joint stage effect (methylation, expression and protein equal) to 0.25, 0.5 and 1.0 at 200 samples, three stages and two measured blocks, with 0.25 as the configuration baseline so that no column repeats a slice column.

#### Scenario: Stage-effect profiles enumerate three columns

- **WHEN** each stage-effect profile is loaded and enumerated
- **THEN** it yields three design columns at joint stage effects 0.25 (the baseline column), 0.5 and 1.0, each with one zero-effect anchor and one cell per (trajectory mode, nonzero effect)
- **AND** every other generator setting of every column equals the slice baseline's
- **AND** at every column the methylation, expression and protein stage effects all equal that column's value

#### Scenario: Stage-effect profiles are pinned and fit the headroom

- **WHEN** each stage-effect profile is loaded and enumerated
- **THEN** its cell identifiers, parameter signatures and matched seeds equal those in the committed post-floor identity fixture
- **AND** enumeration completes without a headroom failure

### Requirement: Anchor reproduction compares anchors by design column

The anchor reproduction tool SHALL offer a by-column mode that pairs a candidate run's zero-effect anchors with a reference run's on (design column, replicate index) and reports, per column, the number of pairs, the seed and statistic matches at the relative tolerance, whether p-values are comparable (equal permutation counts) and, when comparable, whether they match. Without the mode its output SHALL be unchanged.

#### Scenario: Sibling profiles reproduce each other per column

- **WHEN** the magnitude and geometry slice runs are compared by column
- **THEN** every column reports every replicate paired, with matching generator seeds, observed statistics and p-values

#### Scenario: Pilot datasets are reproduced at the first replicates

- **WHEN** a slice run at 999 permutations is compared by column against the pilot run at 199 permutations
- **THEN** every slice column pairs replicate indices 0 through 99 with matching seeds and observed statistics
- **AND** p-values are reported as not comparable and excluded from the match verdict

#### Scenario: A column without a counterpart is listed

- **WHEN** one run has an anchor at a design column the other run lacks
- **THEN** the summary lists that column with zero pairs, and the match verdict covers the paired columns only

#### Scenario: Default mode is unchanged

- **WHEN** the tool is run without the by-column mode on single-anchor runs
- **THEN** its CSV and summary equal those produced before the mode existed

### Requirement: Dated design-guidance paper-grade report

The study SHALL produce a dated findings report for the paper-grade slice, tied to the exact configurations, code revision, parameter signatures, software versions, record counts, failure counts and reproduction commands. It SHALL state every column's rejection rates with Monte Carlo uncertainty at 500 replicates, each column's stage signal and selected rank, the predeclared verdicts, the reproduction checks and the design-guidance statement the paper carries.

#### Scenario: Report states every column's size verdicts

- **WHEN** the report is read
- **THEN** for every design column it gives the zero-effect anchor rate per statistic with its anchor-decision classification against the bound at R = 500
- **AND** translation's `angle` rate at its largest effect with the column-scoped Type I verdict

#### Scenario: Report states the reproduction checks

- **WHEN** the report is read
- **THEN** it reports that the sibling profiles' anchors reproduce each other at every column with identical p-values
- **AND** that the slice reproduces the pilot's datasets at replicate indices 0 through 99 at every column
- **AND** it never pools pilot and slice replicates

#### Scenario: Report carries the guidance statement

- **WHEN** the report is read
- **THEN** for each trajectory component it states what a cohort of each sample size and block count can detect, with every `angle` reading accompanied by the column's selected rank
- **AND** it states the limitations: translation measured at three effects, block count confounded with feature count, the structural-analogue blocks, the saturated stage signal at the core stage effect, and the rank confound of orientation power

#### Scenario: Report covers the stage-effect arm only when it ran

- **WHEN** the stage-effect profiles were run
- **THEN** the report gives the same per-column content for every stage-effect column
- **AND** otherwise it states that the arm was not run and why

## MODIFIED Requirements

### Requirement: Study evaluates results against pre-specified targets

MOTCO SHALL evaluate summarized results against the acceptance targets declared in the configuration and report, per target, whether it is met given Monte Carlo uncertainty.

A Type I control target or a specificity target MUST declare its bound form as one of two values, defaulting to the two-sided form when the key is absent. Under the **two-sided** form the target is met when the empirical rejection rate lies within `k` empirical Monte Carlo standard errors of alpha on either side. Under the **upper** form the target is met when the empirical rejection rate is at or below `alpha + k · sqrt(alpha · (1 − alpha) / R)`, with `R` the number of available replicates — the same bound the Phase 4 gate applies to its control and Type I inflation rules — so that a control rejecting less often than alpha is met. `k` is the target's standard-error tolerance in both forms. The report MUST record the bound form, the alpha-based standard error and the signed deviation from alpha only for a target declared with the upper form; a target in the default form MUST produce exactly the observations and rationale it produced before the bound form existed, so that every committed acceptance report regenerates byte-identically. Power monotonicity targets have no bound form.

A Type I control target MAY additionally name a negative-control `trajectory_mode` (`none` or `translation`) and, only beside it, a `statistic`. Such a **column-scoped** entry is evaluated at every design column — the baseline column and every design point — on that control's cell at the column (the zero-effect anchor for `none`, the largest enumerated effect for `translation`), once per statistic (the named one, or every statistic), under the entry's bound form, and its observations MUST carry the column's design coordinates, whether it is the baseline, the cell, phase, mode and effect. It MUST NOT read the Type I baseline phase. An entry without a mode is evaluated on the Type I phases exactly as before, and a configuration whose entries name no mode MUST produce acceptance reports and configuration echoes byte-identical to before the fields existed. A mode outside the negative controls, a `statistic` without a mode, or an unknown statistic MUST be rejected at load with an error naming the entry.

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

#### Scenario: Column-scoped entry is evaluated at every design column

- **WHEN** a Type I entry naming `translation` and `angle` under the upper bound at alpha `0.05` and tolerance `2` evaluates a design-grid study with six columns
- **THEN** the acceptance report holds six evaluations for it, one per column, the baseline column included, each on that column's translation cell at its largest effect, each judged against `alpha + 2 · sqrt(alpha (1 − alpha) / R)` with `R` that cell's available replicates, and each carrying the column's design coordinates and whether it is the baseline

#### Scenario: Column-scoped entry without a statistic covers every statistic

- **WHEN** a Type I entry names `none` and no statistic
- **THEN** every design column's zero-effect anchor is evaluated once per statistic, and the Type I baseline phase's anchor is not evaluated by it

#### Scenario: Column without the control is reported unavailable

- **WHEN** a column has no completed replicate of the named control
- **THEN** the entry's evaluation for that column is neither met nor unmet and its rationale states that no replicate is available

#### Scenario: Entries without a mode are unchanged

- **WHEN** a study whose Type I entries name no mode is loaded, dumped and reported before and after the column-scoped fields exist
- **THEN** its acceptance report CSV and JSON and its configuration echo are byte-identical

#### Scenario: Invalid column-scoped entry is rejected

- **WHEN** a Type I entry names a mode that is not a negative control, names a statistic without a mode, or names an unknown statistic
- **THEN** loading the configuration raises a validation error naming the entry and the offending value

### Requirement: Design-guidance figure

The study SHALL provide a script that reads the design-point operating tables of the two core profiles and writes:
- for magnitude (`delta`), orientation (`angle`) and shape (`shape`), a figure with one panel per (stage count, block set) that plots the target statistic's rejection rate against sample size, one line per effect size, with Monte Carlo error bars;
- a control figure;
- the plotted values as CSV.

The script SHALL take a grade label, defaulting to the pilot's, that appears in the figure titles, and SHALL lay out one panel row per stage count present in the tables.

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

#### Scenario: A single stage count renders one panel row

- **WHEN** the script is run on reports whose operating tables carry one stage count
- **THEN** each mode figure has two panels (that stage count × two or three blocks) and the control figure lays out the same way, and the CSV rows are exactly the rows plotted

#### Scenario: Grade label is in the titles and defaults to the pilot's

- **WHEN** the script is run with a grade label
- **THEN** every figure title carries that label
- **AND** when it is run without one, the titles equal those written before the option existed
