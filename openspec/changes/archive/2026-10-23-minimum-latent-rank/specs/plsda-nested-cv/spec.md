# Spec Delta

## MODIFIED Requirements

### Requirement: Inner-loop selection by mean AUROC across folds
`plsda_doubleCV` SHALL select the number of latent variables for each outer fold by computing the mean inner-CV AUROC across the V inner folds for each candidate `n_LV ∈ {min_components, …, max_components - 1}`, where `min_components` is a keyword-only argument with default `1`, and then choosing the `n_LV` that maximizes this mean. Ties SHALL be broken by selecting the smaller `n_LV`. `min_components` MUST satisfy `1 ≤ min_components < max_components`; otherwise the function SHALL raise an error naming both values before any fold is fitted.

#### Scenario: Per-fold scores are averaged before argmax
- **WHEN** the inner CV produces AUROCs `{n_LV=1: [0.7, 0.6, 0.8], n_LV=2: [0.85, 0.85, 0.85], n_LV=3: [0.9, 0.5, 0.6]}` for a given outer fold
- **THEN** the selected `n_LV` for that outer fold is `2` (mean 0.85, the highest mean) — not `3` (single-fold max 0.9), and not `1` (lowest mean)

#### Scenario: Tie at the maximum mean prefers the smaller n_LV
- **WHEN** two or more candidate `n_LV` values share the maximum mean inner AUROC
- **THEN** the smaller `n_LV` is chosen

#### Scenario: Candidates below the floor are never evaluated
- **WHEN** `plsda_doubleCV` is called with `min_components=2` and `max_components=5`
- **THEN** only `n_LV ∈ {2, 3, 4}` are fitted in the inner CV, every outer fold's selection is at least 2, every `LV` in the table is at least 2, and every stored model has at least 2 components

#### Scenario: Floor at or above the largest candidate is rejected
- **WHEN** `plsda_doubleCV` is called with `min_components=5` and `max_components=5`, or with `min_components=0`
- **THEN** the function raises an error naming the floor and the maximum before any fold is fitted

#### Scenario: Floor of 1 is the previous rule
- **WHEN** `plsda_doubleCV` is called without `min_components`, or with `min_components=1`, on a fixed input and `random_state`
- **THEN** its table and models equal those produced before the argument existed

### Requirement: Public signature and table row count are preserved
The public signature of `plsda_doubleCV` SHALL remain unchanged in its existing positional and keyword arguments and in its return-type structure (`{"table": pd.DataFrame, "models": list[PLSRegression]}`). It MAY gain keyword-only arguments whose defaults reproduce the previous behaviour exactly: `scale` (default `True`), `nuisance` (default `None`), `strata` (default `None`) and `min_components` (default `1`). The number of rows in `table` SHALL equal `n_repeats`. The number of entries in `models["models"]` SHALL equal `n_repeats`.

#### Scenario: Row count matches n_repeats
- **WHEN** `plsda_doubleCV(..., n_repeats=R)` returns
- **THEN** `result["table"].shape[0] == R` and `len(result["models"]) == R`

#### Scenario: Column 1 is LV; column 2 is AUROC (positional contract)
- **WHEN** existing callers access `result["table"].iloc[:, 1]` and `result["table"].iloc[:, 2]`
- **THEN** they continue to receive the LV (integer) and AUROC (float) values respectively, with no change in semantics for positional readers

#### Scenario: Defaults reproduce the previous output
- **WHEN** `plsda_doubleCV` is called without `scale`, `nuisance`, `strata` and `min_components` on a fixed input and `random_state`
- **THEN** its table and models equal those produced before the arguments existed
