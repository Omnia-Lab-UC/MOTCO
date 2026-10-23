# Spec Delta

## ADDED Requirements

### Requirement: Candidate floor option for motco plsr
`motco plsr` SHALL accept `--min-components N` (default `1`) and pass it to the cross-validation as the smallest number of latent variables it considers, so the candidates are `N … max_components − 1`. A value below 1, or not below `--max-components`, SHALL exit non-zero with a message naming both values. Without the flag the command's table, VIP and score outputs SHALL be unchanged.

#### Scenario: Floor restricts the CV table and the score model
- **WHEN** `motco plsr --min-components 2 --max-components 6 ...` runs with `--out-scores` and without `--n-components`
- **THEN** every `LV` in the table is at least 2 and the score matrix has at least 2 columns

#### Scenario: Default keeps the previous output
- **WHEN** `motco plsr` runs without `--min-components`
- **THEN** its table, VIP and score outputs equal those produced before the flag existed

#### Scenario: Infeasible floor is rejected
- **WHEN** `--min-components 6 --max-components 6` is given
- **THEN** the command exits non-zero and the message names both values
