# Spec Delta

## MODIFIED Requirements

### Requirement: permutation module exposes public symbols
`motco.stats.permutation` SHALL export `RRPP` and `residual_correction_scale`.

#### Scenario: Direct import from permutation submodule
- **WHEN** a caller does `from motco.stats.permutation import RRPP`
- **THEN** the import succeeds and the function is callable

#### Scenario: Residual-correction helper is importable
- **WHEN** a caller does `from motco.stats.permutation import residual_correction_scale`
- **THEN** the import succeeds and the function returns one multiplier per sample for a reduced model and a correction
