# plsda-nested-cv Specification

## Purpose
Define the canonical double-nested cross-validation contract for `plsda_doubleCV`: how the inner CV selects the number of latent variables (mean inner-fold AUROC, parsimony tie-break), how the outer CV aggregates per-repeat (mean and sample std across K outer folds), what is stored in `models["models"]` (full-data refits with the per-repeat modal n_LV), and the resulting `table` schema (`rep`, `LV`, `AUROC`, `AUROC_std`).
## Requirements
### Requirement: Inner-loop selection by mean AUROC across folds
`plsda_doubleCV` SHALL select the number of latent variables for each outer fold by computing the mean inner-CV AUROC across the V inner folds for each candidate `n_LV ∈ {1, …, max_components - 1}`, and then choosing the `n_LV` that maximizes this mean. Ties SHALL be broken by selecting the smaller `n_LV`.

#### Scenario: Per-fold scores are averaged before argmax
- **WHEN** the inner CV produces AUROCs `{n_LV=1: [0.7, 0.6, 0.8], n_LV=2: [0.85, 0.85, 0.85], n_LV=3: [0.9, 0.5, 0.6]}` for a given outer fold
- **THEN** the selected `n_LV` for that outer fold is `2` (mean 0.85, the highest mean) — not `3` (single-fold max 0.9), and not `1` (lowest mean)

#### Scenario: Tie at the maximum mean prefers the smaller n_LV
- **WHEN** two or more candidate `n_LV` values share the maximum mean inner AUROC
- **THEN** the smaller `n_LV` is chosen

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
