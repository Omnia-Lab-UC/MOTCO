# Spec Delta

## ADDED Requirements

### Requirement: The score-coordinate metric of the PLS measurement is documented

The documentation SHALL state that `delta`, `angle` and `shape` are measured in PLS score coordinates, `X·W·(PᵀW)⁻¹`, a sheared image of the fitted span rather than an orthonormal projection. It SHALL state what this does to size and to specificity, and cite the measured residuals.

#### Scenario: Reader asks what coordinates the statistics use
- **WHEN** the reader consults the measurement-visibility documentation
- **THEN** it states that the production statistics are computed in the PLS score coordinates, and that these coordinates do not preserve the standardized feature space's lengths and angles within the span

#### Scenario: Reader asks whether size is affected
- **WHEN** the reader asks whether the score coordinates affect Type I error
- **THEN** the documentation states that under the null both groups share one population trajectory, which any linear map keeps equal, so size is unaffected

#### Scenario: Reader asks whether specificity is affected
- **WHEN** the reader asks whether a shape-only or orientation-only difference can register in another statistic
- **THEN** the documentation states that a difference pure in the standardized feature space can move `delta`, `angle` or `shape` in score coordinates
- **AND** it cites the leak decomposition's score-coordinate residuals per column (planted-scale bend, e = 0.5), including that they are largest on the SEA-AD table and at four stages

#### Scenario: Reader asks what would remove the residual
- **WHEN** the reader asks how the score-coordinate route could be avoided
- **THEN** the documentation names measuring in orthonormal span coordinates as the candidate, states that it changes the production measurement and is not validated, and links the roadmap item that holds it
