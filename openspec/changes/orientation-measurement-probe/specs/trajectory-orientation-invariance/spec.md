## ADDED Requirements

### Requirement: Net-displacement orientation estimator is an opt-in alternative

MOTCO SHALL offer a second orientation estimator, selected by an explicit option on trajectory estimation and on RRPP. The estimator returns the unit vector of each group's net displacement, from its first stage's LS mean to its last stage's LS mean, in stage order. The pairwise `angle` SHALL be the angle between the two groups' net-displacement unit vectors, in degrees in [0, 180]. The PC1 estimator SHALL remain the default. When the option is absent or names the PC1 estimator, every output SHALL be identical to the output produced before the option existed. Only `angle` depends on the option: `delta` and `shape` SHALL be identical under both estimators. When either group's net displacement has zero length within numerical tolerance, `angle` SHALL be reported as unavailable, as it is for the PC1 estimator, and never as a non-finite value.

#### Scenario: Default estimation is unchanged

- **WHEN** trajectory differences or RRPP are computed without selecting an estimator
- **THEN** the returned `deltas`, `angles`, `shapes`, and permutation distributions equal those produced before the option existed, for the same inputs and seed

#### Scenario: Net displacement is the first-to-last direction

- **WHEN** the net-displacement estimator is selected for two groups with stage means `a_1..a_K` and `b_1..b_K`
- **THEN** the pairwise `angle` equals the angle between `a_K − a_1` and `b_K − b_1`

#### Scenario: Two-stage estimators agree

- **WHEN** a trajectory has exactly two stages
- **THEN** the net-displacement and PC1 estimators return the same pairwise `angle` within numerical tolerance

#### Scenario: Net displacement is stable on an isotropic configuration

- **WHEN** both groups share an equilateral three-stage configuration with the same stage order
- **THEN** the net-displacement `angle` is zero within numerical tolerance

#### Scenario: Reversed progression is 180 degrees

- **WHEN** one group's stage order is the reverse of the other's over the same stage means
- **THEN** the net-displacement `angle` is 180 degrees within numerical tolerance

#### Scenario: Only angle depends on the estimator

- **WHEN** the same inputs are estimated under both estimators
- **THEN** `deltas` and `shapes` are identical, and only `angles` may differ

#### Scenario: RRPP uses the selected estimator throughout

- **WHEN** RRPP runs with the net-displacement estimator selected
- **THEN** the observed `angle` and every permuted `angle` are computed with the net-displacement estimator
