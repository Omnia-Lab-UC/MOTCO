# Spec Delta

## ADDED Requirements

### Requirement: motco de exposes the orientation estimator selection

`motco de` SHALL accept `--orientation-estimator` with the values `pc1` and `net_displacement`. Its default SHALL be `pc1`, so every existing invocation produces output identical to the output produced before the flag existed. When `net_displacement` is selected, the estimate and every RRPP permutation SHALL use the net-displacement estimator, and the written JSON SHALL state the estimator used. The CLI reference SHALL state that `net_displacement` is the estimator the simulation harness and the published operating characteristics describe.

#### Scenario: Default de output is unchanged

- **WHEN** `motco de` runs without `--orientation-estimator`
- **THEN** its JSON output equals the output produced before the flag existed for the same inputs and seed

#### Scenario: Net displacement is selectable on the command line

- **WHEN** `motco de` runs with `--orientation-estimator net_displacement` and RRPP permutations
- **THEN** the observed and permuted `angles` are the net-displacement angles, and the JSON states `net_displacement` as the estimator

#### Scenario: Unknown estimator is rejected at the command line

- **WHEN** `motco de` runs with an estimator name other than `pc1` or `net_displacement`
- **THEN** the command exits non-zero listing the allowed values
