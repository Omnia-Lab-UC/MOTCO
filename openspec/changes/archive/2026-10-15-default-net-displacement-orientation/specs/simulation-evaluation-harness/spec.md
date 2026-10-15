# Spec Delta

## REMOVED Requirements

### Requirement: Harness accepts an orientation estimator selection
**Reason**: The requirement made PC1 the default and promised that the default selection is invisible in records, signatures and cell ids. Both are superseded: net displacement is the default after the paper-grade comparison selected it, and it is the PC1 selection, not the default, that stays invisible.
**Migration**: Replaced by "Harness measures orientation with net displacement by default" below. Callers wanting PC1 pass `orientation_estimator="pc1"`; committed configurations are pinned per the trajectory-power-study requirement "Committed configurations pin the PC1 orientation estimator".

## ADDED Requirements

### Requirement: Harness measures orientation with net displacement by default

The evaluation parameters SHALL accept an orientation estimator selection, either the net-displacement estimator (the default) or the PC1 estimator. The selection SHALL be used for the observed `angle`, for every RRPP permutation, and for every realized-geometry checkpoint `angle`. When the PC1 estimator is selected, the selection SHALL be omitted from the persisted record, the parameter signature and the cell identity, so that PC1 evaluation results, signatures and cell ids are identical to those produced before the selection existed. When the selection holds its default, or names the net-displacement estimator explicitly, the persisted record SHALL state the estimator and the parameter signature and cell identity SHALL include it. An unknown estimator SHALL be rejected with a clear validation error.

#### Scenario: Net displacement is the default

- **WHEN** the evaluation parameters omit the estimator
- **THEN** the observed `angle`, the RRPP null and the realized-geometry `angle` values use the net-displacement estimator, and the record states `net_displacement` as the selected estimator

#### Scenario: PC1 estimator is invisible

- **WHEN** the evaluation parameters name the PC1 estimator
- **THEN** the evaluation result, its persisted record, its parameter signature and its cell id equal those produced before the selection existed

#### Scenario: Net-displacement estimator is recorded

- **WHEN** the evaluation parameters select the net-displacement estimator, by default or by name
- **THEN** the observed `angle`, the RRPP null and the realized-geometry `angle` values use it, and the record states the selected estimator

#### Scenario: Unknown estimator is rejected

- **WHEN** the evaluation parameters name an unknown estimator
- **THEN** the harness raises a clear validation error listing the allowed estimators
