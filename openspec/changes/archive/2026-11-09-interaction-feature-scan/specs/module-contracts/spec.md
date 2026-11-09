# Spec Delta

## MODIFIED Requirements

### Requirement: interaction module exposes public symbols
`motco.stats.interaction` SHALL export `interaction_split_test`, `interaction_feature_scan` and their result types, and the flat `motco.stats` namespace SHALL return the same functions.

#### Scenario: Direct import from interaction submodule
- **WHEN** a caller does `from motco.stats.interaction import interaction_split_test`
- **THEN** the import succeeds and the function is callable

#### Scenario: Flat namespace returns the same object
- **WHEN** a caller does `from motco.stats import interaction_split_test`
- **THEN** the import succeeds and returns the same object as the direct submodule import

#### Scenario: The feature scan is exported from both paths
- **WHEN** a caller does `from motco.stats.interaction import interaction_feature_scan` and `from motco.stats import interaction_feature_scan`
- **THEN** both imports succeed and return the same callable object
