# Spec Delta

## ADDED Requirements

### Requirement: sparse_pls module exposes public symbols
`motco.stats.sparse_pls` SHALL export `SparsePLSDA`, `SparsityDesign`, `SPARSITY_SCOPES` and `DEFAULT_KEEP_FRACTIONS`, and the flat `motco.stats` namespace SHALL return the same objects.

#### Scenario: Direct import from the sparse_pls submodule
- **WHEN** a caller does `from motco.stats.sparse_pls import SparsePLSDA`
- **THEN** the import succeeds and the class is instantiable

#### Scenario: Flat namespace returns the same object
- **WHEN** a caller does `from motco.stats import SparsePLSDA`
- **THEN** the import succeeds and returns the same object as the direct submodule import
