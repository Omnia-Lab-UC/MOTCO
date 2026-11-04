# Spec Delta

## MODIFIED Requirements

### Requirement: Visibility of group differences in the stage-supervised measurement is documented

The documentation SHALL state that the production stage-supervised PLS measurement sees a group difference only through its component inside the fitted basis's span. It SHALL state which differences that span contains, why a rotation of the stage program onto other features can cancel, and what a non-rejection therefore means. Each claim SHALL link to the study that measured it.

#### Scenario: Reader asks which differences the measurement can see
- **WHEN** the reader consults the measurement-visibility documentation
- **THEN** it explains that differences re-weighting the stage program both groups share lie inside the span and are visible
- **AND** it explains that at rank K − 1 with balanced cells the span is that of the pooled stage means, so an equal-length rotation of the program onto other features is perpendicular to it (the rhombus argument) and cancels

#### Scenario: Reader asks when the rhombus holds beyond two stages
- **WHEN** the reader consults the rhombus argument for three or more stages
- **THEN** the documentation states that the cancellation is exact when B's steps keep the lengths and mutual angles of A's and are orthogonal to every one of A's steps (B's program moved onto features A's does not use)
- **AND** it states that a rotation within A's own features is partly or wholly visible, the stage-span rotation being the wholly visible case

#### Scenario: Reader asks how out-of-span differences can still register
- **WHEN** the reader consults the same documentation
- **THEN** it lists the incidental routes (CV rank above K − 1, unequal cell sizes or feature scales, correlated noise, nonlinear units) and states that through them a rotation registers as a mix of `delta`, `angle` and `shape` rather than as `angle` alone
- **AND** it cites the orientation probe's measurement that the label-free `pca` and the two-sided `pls_within_group` spaces keep the full rotation, with the costs that probe recorded

#### Scenario: Reader compares a latent space with concat
- **WHEN** the reader consults the statement that every linear latent space is a projection of `concat`
- **THEN** the documentation states that a projection cannot see a difference `concat` does not contain, and that it can still have more power than `concat` by discarding noise directions, citing the probe's Phase 6 `shape` power (0.19 for `pls_cv` against 0.00 for `concat`)

#### Scenario: Reader interprets a non-rejection
- **WHEN** the reader interprets a non-rejection of `angle` or `shape`
- **THEN** the documentation states that it means no detectable difference inside the stage-supervised span, not the same direction or shape in feature space
- **AND** it cites the span-capture shares measured for every simulator construction

## ADDED Requirements

### Requirement: Each candidate space's visibility profile is documented

The measurement-visibility documentation SHALL report, for each space in the latent-space comparison, the span capture and the target power of the stage-span orientation as its alignment with the stage span moves from 0° to 90°. It SHALL report each space's size at the controls and the cost on protected statistics, at every design, each linked to the dated report.

#### Scenario: Reader chooses a measurement space for a question
- **WHEN** the reader consults the documentation to decide which space matches a question
- **THEN** it shows, for each space, how capture and orientation power change with alignment, its anchor and translation size, and its in-span `angle`, `shape` and `delta` power against the production space

#### Scenario: Reader asks whether blindness depends on the features touched
- **WHEN** the reader consults the sparse-versus-dense contrast at full alignment
- **THEN** the documentation states whether the comparison found capture and power to depend on the out-of-span direction's support at a fixed norm, with the measured difference and its standard error
