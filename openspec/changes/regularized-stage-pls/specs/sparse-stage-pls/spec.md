# Spec Delta

## Purpose

Defines the sparse stage-PLS-DA estimator MOTCO may fit its production latent space with: a PLS-DA whose components load on a budget of features, chosen in one of two scopes, that equals ordinary PLS-DA at the full budget and exposes the same fitted surface.

## ADDED Requirements

### Requirement: Each component loads on a budget of features

The estimator SHALL fit PLS-DA components in which the X weight vector of every component is soft-thresholded to a budget of features: the features with the largest absolute weights are kept, every other weight is zero, and the kept weights are shrunk toward zero by the largest discarded absolute weight. The budget SHALL be read in one of two scopes:
- `joint`: one count over all features;
- `block`: one count per declared block, applied within each block.

In either scope every budget SHALL be at least one feature, and the count SHALL be derived from a keep fraction of the scope's feature count, rounded to the nearest integer and floored at one.

#### Scenario: Joint budget
- **WHEN** an estimator with scope `joint` and keep fraction 0.1 is fitted on 500 features
- **THEN** every component's weight vector has exactly 50 non-zero entries

#### Scenario: Block budget
- **WHEN** an estimator with scope `block` and keep fraction 0.1 is fitted on blocks of 367 and 131 features
- **THEN** every component's weight vector has 37 non-zero entries in the first block and 13 in the second

#### Scenario: Smallest block keeps one feature
- **WHEN** a block budget rounds to zero for some block
- **THEN** that block keeps one feature per component

### Requirement: The full budget is ordinary PLS-DA

With every keep fraction equal to 1, the estimator's scores, rotations, loadings and predictions SHALL equal those of the ordinary PLS-DA fit on the same rows within numerical tolerance (1e-6 relative), up to the signs of the components, under per-feature rescaling and without it.

#### Scenario: Full budget equals PLS
- **WHEN** an estimator with keep fraction 1 and an ordinary PLS-DA at the same rank are fitted on the same rows
- **THEN** their score matrices agree to 1e-6 up to column signs and their predictions agree to 1e-6

### Requirement: The estimator exposes the fitted PLS surface

The fitted estimator SHALL expose the X weights, X loadings, Y loadings, X rotations, X scores, input centre and input scale, the transformation of new rows to latent coordinates, the inverse transformation of latent coordinates to input rows and the prediction of the outcome, with the shapes and meanings of the ordinary PLS-DA estimator, so that every consumer of a fitted PLS-DA model in MOTCO accepts it unchanged. It SHALL also expose the resolved budget per component and scope and the support (the kept features) per component.

#### Scenario: Consumers accept the sparse model
- **WHEN** a fitted sparse estimator is handed to the latent map, the VIP computation, attribution and the span-capture measurement
- **THEN** each runs without change, and the latent map's basis spans the sparse rotations

### Requirement: The estimator is deterministic and agrees with the reference implementation

The fit SHALL be deterministic (no random state) and SHALL agree with the committed reference fixture produced once by the reference implementation of sparse PLS-DA (regression mode, no per-feature rescaling): the support of every component's weight vector equal, and the weights and scores equal to 1e-6 up to column signs. No external software SHALL be needed at test time.

#### Scenario: Reference fixture
- **WHEN** the estimator is fitted on the fixture's rows at its rank and budgets
- **THEN** every component's support equals the fixture's and the weights and scores agree to 1e-6 up to column signs

### Requirement: Malformed budgets are rejected before fitting

An unknown scope, a keep fraction outside (0, 1], a block scope without block membership, or a membership that does not cover every feature exactly once SHALL be rejected before any component is fitted, naming the offending value.

#### Scenario: Fraction out of range
- **WHEN** a keep fraction of 0 or 1.5 is given
- **THEN** the estimator is rejected naming the fraction, and nothing is fitted
