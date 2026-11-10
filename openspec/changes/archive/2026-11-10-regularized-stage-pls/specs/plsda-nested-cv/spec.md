# Spec Delta

## ADDED Requirements

### Requirement: Joint selection of rank and keep fraction

When a sparsity design is given, the inner loop SHALL score every candidate pair of rank (the floor to the maximum less one) and keep fraction (the design's ladder) by the mean held-out AUROC across the inner folds, and select the pair with the highest mean; on an exact tie the smaller rank and, within a rank, the fuller budget SHALL win. The outer fold SHALL evaluate the selected pair, each repeat SHALL aggregate the rank by mode with the smaller-value tie-break and the fraction by mode with the fuller-budget tie-break, and the per-repeat refits SHALL be sparse models at the repeat's pair. The table SHALL gain a `keep_fraction` column. Without a sparsity design the function SHALL be byte-identical to its previous output, with no new column.

#### Scenario: Pair selection
- **WHEN** the inner folds give their highest mean AUROC to rank 3 at fraction 0.2
- **THEN** the outer fold is evaluated at rank 3 and fraction 0.2, and the repeat's row records both

#### Scenario: Tie to the fuller budget
- **WHEN** rank 2 at fractions 1.0 and 0.5 tie for the highest mean inner AUROC
- **THEN** fraction 1.0 is selected

#### Scenario: A ladder of one full budget reproduces production
- **WHEN** the sparsity design's ladder is the single fraction 1.0
- **THEN** the table's rank and AUROC columns and the refit models' scores equal those of the run without a sparsity design (within 1e-6 up to component signs)

#### Scenario: Absent design is byte-identical
- **WHEN** the function runs without a sparsity design
- **THEN** its table and models equal the committed identity fixture byte for byte
