# Spec Delta

## Purpose

Present MOTCO's simulation evidence (operating characteristics, specificity, design guidance and the measurement decisions behind them) as a curated set of standalone, provenance-tagged figures for talks, built from committed aggregates.

## ADDED Requirements

### Requirement: One command renders the figure set
MOTCO SHALL provide `scripts/presentation_figures.py`. It renders the nine named figures into an output directory given on the command line, either all of them or a named subset. Figures SHALL be independent, so a failure in one does not alter another's files.

#### Scenario: Render every figure
- **WHEN** the script runs with an output directory and no figure selection
- **THEN** it writes every file of the nine figures (`mode_geometry`, `mode_simulated`, `size_atlas`, `specificity_n1200`, `specificity_n80`, `power_vs_n`, `power_vs_stage_signal`, `estimator_comparison`, `covariates_null`), a manifest and a caption file, and exits 0

#### Scenario: Render a subset
- **WHEN** the script runs with a selection naming one figure
- **THEN** it writes only that figure's files and the manifest entries for it

#### Scenario: Unknown figure name
- **WHEN** the selection names a figure that is not in the set
- **THEN** the script exits non-zero, names the unknown figure and lists the valid names, before writing anything

### Requirement: Figures read committed aggregates only
Every figure except `mode_simulated` SHALL be built only from files tracked in git under `results/`. No figure SHALL read gitignored replicate records or run a simulation study. `mode_simulated` SHALL generate one seeded dataset per mode with the numpy generator and no R runtime.

#### Scenario: Clean clone
- **WHEN** the script runs in a fresh clone of the repository, with no gitignored files present
- **THEN** every figure renders

#### Scenario: Deterministic simulated figure
- **WHEN** `mode_simulated` is rendered twice with the same seed
- **THEN** the plotted coordinates are identical

### Requirement: Source tables are validated before plotting
Each source file SHALL be checked for its required columns and for the rows its figure plots. A missing file, a missing column, a missing expected row or a non-finite plotted rate SHALL fail that figure with an error naming the file and the missing item. A figure SHALL never silently drop a reading.

#### Scenario: Missing source file
- **WHEN** a source file the figure needs is absent
- **THEN** the figure fails with an error naming that path, and no partial figure file is left in the output directory

#### Scenario: Missing expected reading
- **WHEN** a source table lacks a reading the figure plots (for example, one design column of the size atlas)
- **THEN** the figure fails with an error naming the table and the missing reading

### Requirement: Every figure states its measurement
Each figure SHALL carry, in a footer and in the manifest, the measurement of every plotted source: the integration and rank rule, the orientation estimator, the magnitude construction, the permutation correction, and replicates × permutations. A figure that plots sources under different measurements SHALL label each element with its own measurement.

#### Scenario: Footer present by default
- **WHEN** a figure is rendered without options
- **THEN** its footer shows the measurement tag and the source run name

#### Scenario: Footer suppressed for slides
- **WHEN** a figure is rendered with the footer turned off
- **THEN** the figure has no footer, and the manifest still records the measurement tag and sources

#### Scenario: Mixed measurements
- **WHEN** a figure plots one source run under PC1 and another under net displacement
- **THEN** each plotted element is labelled with its estimator

### Requirement: Pre-floor sources are labelled unless the floor is shown not to bind
A source produced before the rank floor (floor 1) SHALL be tagged as pre-floor. It MAY be tagged as the current measurement only for a design whose floor-1 and floor-2 selected ranks have been checked to agree, and the tag SHALL then cite that check.

#### Scenario: Pre-floor source without a check
- **WHEN** a figure plots a floor-1 source for a design with no recorded floor-equivalence check
- **THEN** its measurement tag reads as pre-floor

#### Scenario: Pre-floor source with a check
- **WHEN** a figure plots a floor-1 source for a design whose recorded check found every sampled replicate's rank unchanged under floor 2
- **THEN** its measurement tag reads as the current measurement and cites the check, with its sample size

### Requirement: Statistics have fixed, validated colours
Each statistic SHALL have one colour in every figure: `delta`, `angle` and `shape`. The three colours SHALL pass a categorical palette validation against the figure surface for all pairs: lightness band, chroma floor, colour-vision-deficiency separation, normal-vision separation, and contrast of at least 3:1. Text SHALL use ink colours, never a statistic's colour.

#### Scenario: Same statistic, same colour
- **WHEN** `angle` appears in any two figures
- **THEN** it is drawn in the same colour in both

#### Scenario: Palette validation
- **WHEN** the palette is changed
- **THEN** it is re-validated for all pairs against the figure surface before any figure is committed, and the validation result is recorded in the design

### Requirement: Modes take their target statistic's colour
A trajectory mode with a target statistic SHALL be drawn in that statistic's colour: magnitude in `delta`'s, orientation in `angle`'s, shape in `shape`'s. The negative controls (the zero-effect anchor and translation) SHALL be drawn in neutral ink. Every other category (block count, estimator, permutation correction, effect size) SHALL be encoded by line style, marker or a direct label, never by a new hue.

#### Scenario: Orientation row
- **WHEN** a figure shows the orientation mode
- **THEN** it is drawn in `angle`'s colour

#### Scenario: Secondary category
- **WHEN** a figure compares two estimators for one statistic
- **THEN** the estimators differ in line style or marker, both drawn in that statistic's colour

### Requirement: Figures are slide-ready
Every figure SHALL be drawn on a light surface on a 16:9 canvas, with no text smaller than 12 pt at the rendered size and every axis labelled in plain language rather than internal field names. Each figure SHALL be written as PNG at 300 dpi and as PDF.

#### Scenario: Output formats
- **WHEN** a figure is rendered
- **THEN** a `.png` at 300 dpi and a `.pdf` of the same figure are written

#### Scenario: Plain-language axes
- **WHEN** a figure plots sample size
- **THEN** the axis reads as a sample size (for example, "samples (n)"), not `n_samples` or `generator.n_samples`

### Requirement: Manifest and captions record provenance
The script SHALL write `manifest.json` and `captions.md`. For each figure they SHALL record: its message (one sentence), its files, every source path with its SHA-256, its measurement tags and its caveats. They SHALL also record the git revision, whether the worktree was clean, the MOTCO version and the render date.

#### Scenario: Source hash
- **WHEN** a figure is rendered from a source file
- **THEN** the manifest records that file's path and SHA-256

#### Scenario: Dirty worktree
- **WHEN** the script runs with uncommitted changes
- **THEN** the manifest records the worktree as dirty

### Requirement: Mode geometry schematic
`mode_geometry` SHALL show, on constructed two-dimensional trajectories rather than data, how magnitude, orientation, shape and translation change group B's trajectory relative to group A's. It SHALL name the statistic that reads each mode, and show translation as moving none of the three.

#### Scenario: Four modes shown
- **WHEN** `mode_geometry` is rendered
- **THEN** it shows the magnitude, orientation, shape and translation constructions, each labelled with the statistic it targets, or with "none" for translation

### Requirement: Simulated mode figure
`mode_simulated` SHALL show the zero-effect anchor and the four modes as generated by the numpy generator from a fixed seed. Each mode SHALL be shown as group stage-mean trajectories in a display projection, labelled as display-only and not the measurement. It SHALL write one combined file and one file per mode.

#### Scenario: Combined and per-mode files
- **WHEN** `mode_simulated` is rendered
- **THEN** it writes one combined figure and one figure for each of the anchor and the four modes

### Requirement: Size atlas
`size_atlas` SHALL plot every zero-effect anchor and translation-control rejection rate, for each statistic, from each committed source design: n = 80 to n = 1,501, including the case-study and ADNI configurations. Each reading SHALL show its own α + 2·SE bound and the verdict its source recorded.

#### Scenario: Exceeding reading is visible
- **WHEN** a plotted reading's source verdict is an exceedance
- **THEN** it is drawn with a distinct marker and labelled, not hidden or dropped

#### Scenario: Production arm only
- **WHEN** a source records several analysis arms
- **THEN** the atlas plots the arm that is the production analysis for that configuration, and names it in the measurement tag

### Requirement: Specificity at n = 1,200
`specificity_n1200` SHALL show the mode × statistic rejection-rate matrix at n = 1,200, four stages and three blocks, under net displacement and joint magnitude. Rows are the anchor and the four modes at their top effect, and every cell shows its value.

#### Scenario: Source values
- **WHEN** `specificity_n1200` is rendered
- **THEN** each cell equals the net-displacement column's rejection rate in the source table for that mode, effect and statistic

### Requirement: Specificity at n = 80
`specificity_n80` SHALL show the same matrix layout at n = 80, three stages and two blocks, under the current measurement, at each mode's top effect. Its colour scale SHALL equal `specificity_n1200`'s, so the two can be shown side by side.

#### Scenario: Shared scale
- **WHEN** both specificity figures are rendered
- **THEN** they use the same colour scale, from 0 to 1

### Requirement: Power against sample size
`power_vs_n` SHALL plot each mode's target-statistic power against sample size (n = 80, 200 and 500) for two and three blocks, with ±1 Monte Carlo SE. It SHALL also write one detail figure per statistic that shows each effect size separately.

#### Scenario: Overview and details
- **WHEN** `power_vs_n` is rendered
- **THEN** it writes one overview figure and three detail figures, for `delta`, `angle` and `shape`

### Requirement: Power against stage signal
`power_vs_stage_signal` SHALL plot target-statistic power at n = 200 against the joint stage effect δ ∈ {0.25, 0.5, 1.0, 2.0}. It SHALL show beside it the median CV-selected rank and the median held-out stage AUROC at each δ, as separate panels with their own axes, never as a second y-axis.

#### Scenario: Rank shown separately
- **WHEN** `power_vs_stage_signal` is rendered
- **THEN** power, selected rank and stage AUROC each have their own panel and axis

### Requirement: Estimator comparison
`estimator_comparison` SHALL plot, at the Phase 6 (n = 80) and Phase 5 (n = 1,200) designs, orientation/`angle` power and the shape-mode `angle` response against effect size, for PC1 and net displacement on identical datasets.

#### Scenario: Source values at the top effect
- **WHEN** `estimator_comparison` is rendered
- **THEN** orientation/`angle` at e = 1.0 reads 0.300 (PC1) and 0.648 (net displacement) at Phase 6, and 0.840 and 0.924 at Phase 5, as in the source tables

### Requirement: Covariates and the permutation null
`covariates_null` SHALL plot the covariate-adjusted `delta` null rejection rate against the number of covariate columns (1, 3 and 6), for the uncorrected and the `df`-corrected permutation. Each control and stage strength SHALL be shown with its pooled rate and the bound.

#### Scenario: Pooled rates
- **WHEN** `covariates_null` is rendered
- **THEN** the pooled uncorrected `delta` rate reads 0.0545 at one column and 0.066 at six, and the `df`-corrected rate at six reads 0.0545, as computed from the source tables

### Requirement: Rendered set is committed with provenance
The rendered figure set used for a presentation SHALL be committed under `results/presentation-figures-<render date>/`, together with its manifest and captions. It SHALL be rendered from a clean worktree.

#### Scenario: Committed set
- **WHEN** a rendered set is committed
- **THEN** its manifest records a clean worktree and a git revision that contains the script that rendered it
