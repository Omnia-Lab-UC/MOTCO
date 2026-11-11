# plsda-nested-cv Specification

## Purpose
Define the canonical double-nested cross-validation contract for `plsda_doubleCV`: how the inner CV selects the number of latent variables (mean inner-fold AUROC, parsimony tie-break), how the outer CV aggregates per-repeat (mean and sample std across K outer folds), what is stored in `models["models"]` (full-data refits with the per-repeat modal n_LV), and the resulting `table` schema (`rep`, `LV`, `AUROC`, `AUROC_std`).
## Requirements
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

### Requirement: Per-repeat AUROC is the mean across outer folds
`plsda_doubleCV` SHALL report, for each repeat, the arithmetic mean of the K outer-fold test AUROCs as the `AUROC` value in the returned `table`. The function SHALL NOT report the maximum, minimum, or any single fold's AUROC as the canonical per-repeat value.

#### Scenario: Mean is reported, not max
- **WHEN** the K outer folds for one repeat produce test AUROCs `[0.80, 0.82, 0.78, 0.85, 0.81]`
- **THEN** the corresponding row in `table` reports `AUROC = 0.812` (the mean), not `0.85`

#### Scenario: Saturated folds do not pin the report
- **WHEN** one outer fold lands at AUROC = 1.0 and the others at 0.85, 0.86, 0.83, 0.84
- **THEN** the reported `AUROC` is the mean of all five values (~ 0.876), not 1.0

### Requirement: AUROC_std column reports outer-fold dispersion
The `table` returned by `plsda_doubleCV` SHALL include a column named `AUROC_std` at position 3 (zero-indexed), containing the sample standard deviation (`numpy.std(..., ddof=1)`) of the K outer-fold test AUROCs for each repeat. When K < 2, the value SHALL be `NaN`.

#### Scenario: Column shape and dtype
- **WHEN** `plsda_doubleCV(..., n_repeats=R)` returns successfully
- **THEN** `result["table"].shape == (R, 4)` and `result["table"].columns.tolist() == ["rep", "LV", "AUROC", "AUROC_std"]`

#### Scenario: AUROC_std is non-negative
- **WHEN** any repeat completes with K >= 2 outer folds
- **THEN** the corresponding `AUROC_std` value is finite and `>= 0`

### Requirement: LV column is the modal outer-fold selection per repeat
The `LV` column of the returned `table` SHALL contain, for each repeat, the integer mode of the K per-outer-fold selected `n_LV*` values. Ties SHALL be broken by selecting the smaller value.

#### Scenario: Modal LV with a clear majority
- **WHEN** the K outer folds select `n_LV*` values `[2, 2, 3, 2, 2]`
- **THEN** the reported `LV` is `2`

#### Scenario: Tied modes prefer the smaller LV
- **WHEN** the K outer folds select `n_LV*` values `[2, 2, 3, 3, 4]`
- **THEN** the reported `LV` is `2`

### Requirement: Per-repeat model is refit on full data with the selected LV
`models["models"][i]` SHALL be a `PLSRegression(n_components=LV_i, scale=<scale>, max_iter=1000)` fit on the full input arrays `(X, y_one_hot)`, where `LV_i` is the modal `n_LV*` recorded for repeat `i` in the table and `<scale>` is the function's `scale` argument (default `True`). When a nuisance design is given, the full input `X` SHALL first be residualized with coefficients estimated on all samples. The function SHALL NOT store a model that was fit on a single CV fold.

#### Scenario: Stored model is a full-data fit
- **WHEN** `plsda_doubleCV` returns and `model = result["models"][i]` for any `i`
- **THEN** `model.x_scores_.shape[0] == X.shape[0]` (the model was fit on every input sample)

#### Scenario: Stored model uses the per-repeat selected LV
- **WHEN** `result["table"].iloc[i, 1] == 3`
- **THEN** `result["models"][i].n_components == 3`

#### Scenario: Default scale is unchanged
- **WHEN** `plsda_doubleCV` is called without `scale`
- **THEN** every stored model was fit with `scale=True`, as before

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

### Requirement: Per-fold nuisance residualization in the nested CV
When a nuisance design (conditioning and covariate columns, row-aligned to `X`) is given, every outer rest set and inner training set SHALL estimate the covariate coefficients of `X` on an intercept, the conditioning and the covariate columns from its own rows, residualize its own rows, and residualize the held-out rows with the same coefficients and the held-out covariate columns only. Each design MUST have full column rank, or the function SHALL raise an error naming the repeat and fold.

#### Scenario: Held-out rows use training coefficients

- **WHEN** a nuisance design is given and one outer fold's rest set yields covariate coefficients `γ`
- **THEN** that fold's test rows are residualized as `X_test − C_test γ` and the inner training sets within the rest set estimate their own coefficients

#### Scenario: Rank-deficient fold is an error

- **WHEN** a training set's nuisance design has fewer independent columns than columns
- **THEN** the function raises an error naming the repeat and fold

### Requirement: Optional stratification labels for the CV folds
When `strata` (row-aligned labels) is given, the outer and inner folds SHALL be stratified on `strata` instead of `y`; the AUROC SHALL still score `y`. Without `strata`, the folds SHALL be stratified on `y` as before.

#### Scenario: Same partition, same folds

- **WHEN** `strata` partitions the rows exactly as `y` does, under different label names
- **THEN** the table equals the one produced without `strata`

#### Scenario: Conditioning cells keep every cell in training

- **WHEN** a group × stage nuisance design is given on rows sorted by group within stage, with `strata` set to its conditioning cells
- **THEN** no training design loses a conditioning column and the function completes, where stratifying on `y` alone raises naming the fold

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
