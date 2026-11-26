# Spec Delta

## ADDED Requirements

### Requirement: The simulation-only reading of the measurement candidates is documented

The measurement-visibility documentation SHALL state, for `predictive_orthonormal`, block weighting and sparse stage-PLS, the candidate simulation validation's verdict and each arm's label by design family, with the measured rates. It SHALL state that the decision used simulated data only, and SHALL cite beside it each arm's real-noise readings and the feature count of the simulated designs as limits.

#### Scenario: Reader asks where a candidate beats production
- **WHEN** the reader asks at which designs a candidate is better than the production measurement
- **THEN** the documentation gives the arm's label by design family with the cells and rates behind each label

#### Scenario: Reader asks what a candidate gives up
- **WHEN** the reader asks what a candidate loses
- **THEN** the documentation names the families where it loses size, protected power or specificity, and the surgery and out-of-span cells where it rejects less often than production

#### Scenario: Reader asks about unequal noise by stage
- **WHEN** the reader asks whether the tests hold size when stages differ in noise
- **THEN** the documentation states production's and each candidate's reading on the dispersion designs, and what the within-stage draws read on the same units

#### Scenario: Reader asks how far the result carries
- **WHEN** the reader asks whether the map holds for a real cohort
- **THEN** the documentation states that it was read on 658 Gaussian features at 1 to 8 features per sample, names the SEA-AD and ADNI readings of each arm with their verdicts, and states that neither entered the decision
