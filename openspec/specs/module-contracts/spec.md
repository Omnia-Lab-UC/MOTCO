# module-contracts Specification

## Purpose
Defines the public import surface of the `motco.stats` subpackage: which symbols each of the `design`, `trajectory`, and `permutation` modules exports, that the flat `motco.stats` namespace keeps returning the same objects, and that the retired `sd` module stays unimportable.

## Requirements

### Requirement: design module exposes public symbols
`motco.stats.design` SHALL export `center_matrix`, `get_model_matrix`, and `build_ls_means`.

#### Scenario: Direct import from design submodule
- **WHEN** a caller does `from motco.stats.design import get_model_matrix`
- **THEN** the import succeeds and the function is callable

#### Scenario: Flat namespace unchanged
- **WHEN** a caller does `from motco.stats import get_model_matrix`
- **THEN** the import succeeds and returns the same object as the direct submodule import

### Requirement: trajectory module exposes public symbols
`motco.stats.trajectory` SHALL export `pair_difference`, `estimate_difference`, `estimate_betas`, and `get_observed_vectors`.

#### Scenario: Direct import from trajectory submodule
- **WHEN** a caller does `from motco.stats.trajectory import estimate_difference`
- **THEN** the import succeeds and the function is callable

### Requirement: permutation module exposes public symbols
`motco.stats.permutation` SHALL export `RRPP` and `residual_correction_scale`.

#### Scenario: Direct import from permutation submodule
- **WHEN** a caller does `from motco.stats.permutation import RRPP`
- **THEN** the import succeeds and the function is callable

#### Scenario: Residual-correction helper is importable
- **WHEN** a caller does `from motco.stats.permutation import residual_correction_scale`
- **THEN** the import succeeds and the function returns one multiplier per sample for a reduced model and a correction

### Requirement: sd module no longer exists
`motco.stats.sd` SHALL NOT be importable after this change.

#### Scenario: Direct import from sd fails
- **WHEN** a caller does `import motco.stats.sd`
- **THEN** an `ImportError` or `ModuleNotFoundError` is raised
