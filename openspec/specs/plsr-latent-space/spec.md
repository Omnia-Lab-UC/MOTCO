# plsr-latent-space Specification

## Purpose
Defines the PLS-DA latent-space capability: the public fitting and transform functions in `stats/pls.py`, and the `motco plsr` multi-omics input, supervision-label, and score-output contracts that expose the fitted latent space to downstream trajectory analysis.

## Requirements

### Requirement: fit_plsda_transform function
`stats/pls.py` SHALL expose a public function `fit_plsda_transform(X, y, n_components) -> np.ndarray` that fits a `PLSRegression(n_components=n_components, scale=True)` on the full arrays `X` and `y`, and returns the resulting `x_scores_` matrix of shape `(n_samples, n_components)`.

#### Scenario: Returns correct shape
- **WHEN** called with X of shape (n, p) and n_components=k
- **THEN** the returned array has shape (n, k)

#### Scenario: Consistent with sklearn transform
- **WHEN** called with valid X, y, n_components
- **THEN** the returned scores match `PLSRegression(n_components).fit(X, y).x_scores_`

### Requirement: Multi-omics input mode for motco plsr
`motco plsr` SHALL accept `--input` (repeatable flag, one per omics layer CSV file) as an alternative to `--data` / `--x` / `--y`. When `--input` is provided, the command SHALL standardize each layer independently (z-score per feature; columns with std=0 are left unchanged) and concatenate them horizontally to form X.

#### Scenario: Two input files are concatenated
- **WHEN** `--input methyl.csv --input expr.csv` is provided and both files have n rows
- **THEN** X used for PLS has n rows and (methyl_cols + expr_cols) columns

#### Scenario: Each layer is standardized before concatenation
- **WHEN** input files with different feature scales are provided
- **THEN** each layer's columns have mean ≈ 0 and std ≈ 1 before concatenation

#### Scenario: --input and --data are mutually exclusive
- **WHEN** both `--input` and `--data` are provided
- **THEN** the command exits with a non-zero status and a descriptive error message

### Requirement: Metadata label column for multi-omics mode
When `--input` is used, `motco plsr` SHALL require `--metadata <csv>` and `--label-col <column>` to specify the supervision signal y. The label column is read from the metadata CSV and used as-is (string or numeric) as the outcome for PLS-DA.

#### Scenario: Missing --metadata when using --input
- **WHEN** `--input` is provided without `--metadata`
- **THEN** the command exits with a non-zero status and a descriptive error message

#### Scenario: Missing --label-col when using --input
- **WHEN** `--input` and `--metadata` are provided but `--label-col` is omitted
- **THEN** the command exits with a non-zero status and a descriptive error message

#### Scenario: Label column not found in metadata
- **WHEN** `--label-col foo` is specified but `foo` is not a column in the metadata CSV
- **THEN** the command exits with a non-zero status naming the missing column

### Requirement: Latent space score output
When `--out-scores <path>` is provided, `motco plsr` SHALL:
1. Determine the final number of components: use `--n-components` if specified, otherwise the modal value from the CV table's LV column.
2. Call `fit_plsda_transform(X, y, n_components)` on the full dataset.
3. Save the resulting score matrix as a CSV (no index, column names `lv_1, lv_2, …`).

#### Scenario: Scores file is created with correct shape
- **WHEN** `--out-scores scores.csv` is provided after a successful double-CV run
- **THEN** `scores.csv` exists and has n_samples rows and n_components columns

#### Scenario: --n-components overrides modal LV
- **WHEN** `--n-components 3` is specified
- **THEN** the score matrix has exactly 3 columns regardless of the CV table result

#### Scenario: --out-scores works with single-matrix mode
- **WHEN** `--data` or `--x`/`--y` is used together with `--out-scores`
- **THEN** scores are produced using the same modal-LV logic as multi-omics mode

### Requirement: fit_plsda_model function
`stats/pls.py` SHALL expose a public function `fit_plsda_model(X, y, n_components) -> PLSRegression` that fits a `PLSRegression(n_components=n_components, scale=True)` on the full arrays `X` and one-hot-encoded `y` (matching `fit_plsda_transform`/`plsda_doubleCV` encoding) and returns the fitted estimator. No cross-validation is performed.

#### Scenario: Returns a fitted projector
- **WHEN** called with valid `X`, `y`, and `n_components=k`
- **THEN** the return value is a fitted `PLSRegression` whose `.transform()` maps an `(m, p)` feature matrix to `(m, k)` scores

#### Scenario: fit_plsda_transform is consistent with fit_plsda_model
- **WHEN** `fit_plsda_transform(X, y, k)` and `fit_plsda_model(X, y, k).transform(X)` are computed
- **THEN** the two score matrices are equal

#### Scenario: Exported from the stats package
- **WHEN** a caller writes `from motco.stats import fit_plsda_model`
- **THEN** the import succeeds

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
