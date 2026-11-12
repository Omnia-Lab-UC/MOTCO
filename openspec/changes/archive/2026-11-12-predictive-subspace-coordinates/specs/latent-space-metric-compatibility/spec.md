# Spec Delta

## ADDED Requirements

### Requirement: The measured dimension of the PLS measurement is documented

The documentation SHALL state that the production measurement uses the selected rank of the stage-PLS fit both to regularize the fit and as the number of dimensions the statistics are measured in. It SHALL state what measuring in the predictive subspace keeps and drops, report the pilot of the predictive-subspace coordinates against production, and state the production consequence.

#### Scenario: Reader asks how many dimensions the statistics use
- **WHEN** the reader consults the measurement-visibility documentation
- **THEN** it states that the production statistics are computed in as many dimensions as the cross-validation selects components, and cites the measured cost of dimensions beyond the stage program

#### Scenario: Reader asks what the predictive subspace is
- **WHEN** the reader asks how fewer dimensions can be measured on the same fit
- **THEN** the documentation states that the fitted stage labels move along at most the number of stages minus one directions of the score space, that measuring there keeps the fit and its rank, and that the harness and the real-data analysis offer three coordinate systems on that subspace as an opt-in setting

#### Scenario: Reader asks why a lower rank is not the same thing
- **WHEN** the reader asks whether forcing the fit to fewer components is equivalent
- **THEN** the documentation states that the later components aim the predictive directions, and cites the control arm's rates and capture beside the predictive arms'

#### Scenario: Reader asks what is given up
- **WHEN** the reader asks what the predictive subspace cannot see or report
- **THEN** the documentation states that a difference seen only through the dropped components is no longer read as it was, citing the surgery cell
- **AND** it states how the observed statistics compare with the population's through the same map, per stage effect

#### Scenario: Reader asks which measurement production uses
- **WHEN** the reader asks whether production and `motco analyze` measure in the predictive subspace
- **THEN** the documentation states the pilot's verdict and its size, power, specificity and gain readings, and that production measures in the fitted span's scores until a separate change adopts another measurement, which only a paper grade following `propose_paper_grade` permits
