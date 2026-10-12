# trajectory-orientation-invariance Specification

## Purpose

Defines the production contract for MOTCO trajectory `angle` as principal-axis divergence — the angle between the groups' leading principal axes, each signed by net displacement, which equals a difference in direction of progression only for straight or two-stage trajectories: the estimator that summarizes a multi-stage trajectory's principal axis, the transformations that axis must ignore, and the null configurations that must never be reported as large orientation differences.

## Requirements

### Requirement: Orientation is the principal axis of the stage configuration, signed by progression

MOTCO SHALL estimate a trajectory's orientation as the leading principal axis of its centered stage configuration, with the sign resolved so the returned vector points along the trajectory's direction of progression — from its first stage toward its last. The estimator MUST return a direction, not an axis: reversing the stage order MUST reverse the returned vector.

#### Scenario: Orientation follows the direction of progression

- **WHEN** a trajectory's orientation is estimated
- **THEN** the returned vector has a non-negative inner product with the trajectory's net displacement from its first stage to its last

#### Scenario: Reversing stage order reverses orientation

- **WHEN** the same stage configuration is presented in reversed stage order
- **THEN** the returned orientation vector is negated

#### Scenario: Two-stage orientation is the transition direction

- **WHEN** a trajectory has exactly two stages
- **THEN** its orientation is the unit vector of the transition from the first stage to the second

### Requirement: Orientation is invariant to translation and uniform scale

MOTCO SHALL report a trajectory `angle` of zero, within documented numerical tolerance, between any two trajectories whose stage configurations differ only by translation, by uniform scale, or by both. The angle MUST NOT depend on where the trajectories sit relative to the coordinate origin.

#### Scenario: Translated trajectory has zero angle

- **WHEN** two trajectories have the same stage configuration and one is translated by a constant vector
- **THEN** their pairwise `angle` is zero within numerical tolerance

#### Scenario: Trajectories on opposite sides of the origin have zero angle

- **WHEN** two trajectories with the same stage configuration are translated so they lie on opposite sides of the coordinate origin
- **THEN** their pairwise `angle` is zero within numerical tolerance

#### Scenario: Uniformly scaled trajectory has zero angle

- **WHEN** two trajectories have the same stage configuration and one is uniformly scaled by a positive factor
- **THEN** their pairwise `angle` is zero within numerical tolerance

### Requirement: Bent trajectories do not produce spurious reversals

The orientation sign MUST remain stable for trajectories whose stage configuration bends — including configurations that depart from and return toward their starting region, where the first stage lies near the configuration centroid along the principal axis. Two trajectories that carry no orientation difference MUST NOT be reported as near-antiparallel.

#### Scenario: Bent identical trajectories report zero angle

- **WHEN** two copies of a bent multi-stage trajectory whose first stage projects near zero onto its own centered principal axis are compared
- **THEN** their pairwise `angle` is zero within numerical tolerance

#### Scenario: Sign is stable under small perturbation

- **WHEN** a bent trajectory's stages are perturbed by noise small relative to its extent, repeatedly
- **THEN** the returned orientation does not change sign across those perturbations

#### Scenario: Null configurations never approach antiparallel

- **WHEN** trajectories differing only by translation, uniform scale, or small noise are compared across repeated draws
- **THEN** no pairwise `angle` approaches 180 degrees

### Requirement: Genuine orientation differences are preserved

The estimator SHALL continue to report the true angle between trajectories that genuinely differ in principal-axis orientation, including differences exceeding 90 degrees, and SHALL reproduce the published reference outputs it was ported from. Documentation of the `angle` statistic SHALL state its estimand as principal-axis divergence: the angle between the groups' leading principal axes, each signed by net displacement — which equals the angle between directions of progression only for straight or two-stage trajectories.

#### Scenario: Known angle is recovered

- **WHEN** two trajectories are constructed with a known angle between their signed principal axes
- **THEN** the reported `angle` equals that angle within numerical tolerance

#### Scenario: Reference outputs are reproduced

- **WHEN** the committed reference datasets are analyzed
- **THEN** the reported `angle` values match the committed reference results, including values greater than 90 degrees

#### Scenario: Estimand is documented as principal-axis divergence

- **WHEN** a reader consults the orientation documentation for a multi-stage trajectory
- **THEN** it states that `angle` compares signed principal axes, not per-step transition directions, and that the two coincide only for straight or two-stage trajectories

### Requirement: The sign convention is documented as a deviation from the reference

The implementation SHALL record that the sign convention deviates deliberately from the reference supplement, which anchors on the raw first-stage row, and SHALL state why that anchor is unsuitable: it makes the orientation sign depend on the trajectory's position relative to the coordinate origin, so a pure translation can reverse it.

#### Scenario: Deviation is discoverable at the implementation

- **WHEN** a contributor reads the orientation estimator
- **THEN** it names the reference line it departs from, the reason, and the invariance the chosen anchor provides

### Requirement: Two-stage pairwise angle equals the transition-vector angle

For trajectories with exactly two stages, the pairwise `angle` statistic SHALL equal the angle between the two groups' unit transition vectors (first stage toward second stage). This is an identity, not an approximation: with two stages the centered stage configuration is rank one, its leading principal axis is the transition direction, and the net-displacement sign anchor orients it along the progression. Agreement MUST hold at machine precision, assessed on the cosine of the angle (the arccos map is ill-conditioned near 0 and 180 degrees, so the tolerance contract is stated on the cosine). Trajectories whose two stages coincide are the documented zero-displacement degeneracy and are outside this requirement.

#### Scenario: Synthetic two-stage angle is the transition-vector angle

- **WHEN** two two-stage trajectories with distinct stages are compared via the pairwise `angle` statistic
- **THEN** the cosine of the reported angle equals the inner product of the two unit transition vectors within machine precision

#### Scenario: Example1 fixture angles equal the direct-vector angles

- **WHEN** the committed two-stage reference dataset (`evo_649_sm_example1.csv`) is analyzed
- **THEN** every pairwise `angle` equals the angle computed directly from that pair's unit transition vectors, derived from the same fitted stage means, within machine precision on the cosine

#### Scenario: Committed R angles are the sign-anchor artifact of the direct-vector angles

- **WHEN** the reported two-stage angles are compared against the committed R reference outputs for the two-stage dataset
- **THEN** each committed angle equals either the direct-vector angle or its supplement (180 degrees minus it), and the supplement cases are attributable to the reference's raw first-stage sign anchor

### Requirement: Two-stage regression expectations are pinned to the progression convention

Regression tests that compare two-stage `angle` outputs against committed reference values MUST NOT accept both an angle and its supplement. The expected value SHALL be the direct transition-vector angle under the progression sign convention, so that a sign-flip regression in the orientation estimator fails the comparison rather than passing through a supplementary-angle acceptance.

#### Scenario: A sign flip is not absorbed by supplementary acceptance

- **WHEN** a two-stage regression comparison is evaluated against its expected angle
- **THEN** an output equal to the supplement of the expected angle fails the comparison

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
