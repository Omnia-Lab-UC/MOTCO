# Spec Delta

## ADDED Requirements

### Requirement: interaction module exposes public symbols
`motco.stats.interaction` SHALL export `interaction_split_test` and its result type, and the flat `motco.stats` namespace SHALL return the same function.

#### Scenario: Direct import from interaction submodule
- **WHEN** a caller does `from motco.stats.interaction import interaction_split_test`
- **THEN** the import succeeds and the function is callable

#### Scenario: Flat namespace returns the same object
- **WHEN** a caller does `from motco.stats import interaction_split_test`
- **THEN** the import succeeds and returns the same object as the direct submodule import
