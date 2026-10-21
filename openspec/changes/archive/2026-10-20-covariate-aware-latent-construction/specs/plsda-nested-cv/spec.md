# Spec Delta

## MODIFIED Requirements

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
The public signature of `plsda_doubleCV` SHALL remain unchanged in its existing positional and keyword arguments and in its return-type structure (`{"table": pd.DataFrame, "models": list[PLSRegression]}`). It MAY gain keyword-only arguments whose defaults reproduce the previous behaviour exactly: `scale` (default `True`), `nuisance` (default `None`) and `strata` (default `None`). The number of rows in `table` SHALL equal `n_repeats`. The number of entries in `models["models"]` SHALL equal `n_repeats`.

#### Scenario: Row count matches n_repeats
- **WHEN** `plsda_doubleCV(..., n_repeats=R)` returns
- **THEN** `result["table"].shape[0] == R` and `len(result["models"]) == R`

#### Scenario: Column 1 is LV; column 2 is AUROC (positional contract)
- **WHEN** existing callers access `result["table"].iloc[:, 1]` and `result["table"].iloc[:, 2]`
- **THEN** they continue to receive the LV (integer) and AUROC (float) values respectively, with no change in semantics for positional readers

#### Scenario: Defaults reproduce the previous output
- **WHEN** `plsda_doubleCV` is called without `scale` and `nuisance` on a fixed input and `random_state`
- **THEN** its table and models equal those produced before the arguments existed

## ADDED Requirements

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

