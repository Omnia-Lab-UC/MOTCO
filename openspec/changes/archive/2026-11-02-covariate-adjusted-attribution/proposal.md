# covariate-adjusted-attribution

## Why

Phase 6 closed on 2026-11-02 with its attribution leg withheld ([exit review](../../../docs/reports/phase6-exit-review-2026-11-02.md) §4 D6, §6). Both case-study primaries declare six covariate columns, and `motco analyze` cannot attribute a covariate-adjusted contrast. Two guards block it:
- `analyze_trajectories` refuses attribution whenever covariates are declared (`src/motco/analysis.py:519–523`).
- `analyze_orientation_attribution` accepts adjusted means only as a precomputed `mean_table`, and then requires `bootstrap_replicates = 0` (`src/motco/stats/attribution.py:368–372`).

So a claim from either primary could not be traced to *stable* drivers. Unadjusted attribution is no substitute: under a confounded covariate it estimates a different contrast. The covariate Type I check's unadjusted arm rejected `shape` at 0.364 on null data ([covariate Type I report](../../../docs/reports/covariate-type-i-2026-10-16.md) §2).

This change is roadmap item 11. Item 12, the plasmode demonstration on SEA-AD's real blocks, either lifts the withheld leg or reports it failed. It must run at the configuration in which a real claim would be read, so it cannot start until this change lands.

## What Changes

- **Covariate-adjusted cell means in the attribution primitive.**
  - `analyze_orientation_attribution` gains the keyword-only `covariates` (a row-aligned numeric design-column table) and `covariate_reference` (the covariate row at which the means are evaluated, by default the column means).
  - When `covariates` is given, every group × stage feature mean is the least-squares mean of the cell-means model with the covariates added. That is the model the analysis fits in the latent space, here fitted in feature space and evaluated at the reference row.
  - The observed transitions and the principal-orientation contrast do not depend on the reference row: it shifts every cell by the same vector. The means tables do depend on it.
  - The reconstructed (PLS-captured) and residual components follow through the same frozen model.
  - The configuration records `mean_source = "covariate_adjusted"`, the covariate columns and the reference row.
  - If the covariates are collinear within cells, or `covariates` is combined with a precomputed `mean_table`, the primitive raises an error before computing any driver values.
- **A bootstrap that refits the adjustment.** Each replicate draws the same stratified within-cell resample as the unadjusted bootstrap for the same seed. It then refits the covariate coefficients on that resample, evaluates the means at the fixed reference row, and contrasts them through the frozen preprocessing and PLS model. A resample whose covariates are collinear within cells is an invalid replicate. The number of such replicates is recorded.
- **Attribution with covariates in `motco analyze`.**
  - The covariate guard is removed. With covariates declared, attribution uses the adjusted means at the record's `covariates.reference_row`, together with the adjusted bootstrap.
  - Projecting the attribution's observed means through the fitted model reproduces the analysis's latent `stage_means`.
  - The compact attribution record gains an `adjustment` block only when covariates are declared. `motco analyze --attribution` prints the mean source and the rank-deficient replicate count.
  - Attribution combined with **latent adjustment** is still rejected before fitting, with its own message (design D6).
  - `TrajectoryAnalysisResult` exposes the fitted estimator and the standardized joint feature matrix **in memory only**, and never writes or pickles them. The validation needs them, and so will item 12.
- **Byte identity without covariates.** Without covariates the primitive, `run_attribution_diagnostics` and `analyze_trajectories` write every table and record exactly as in 0.13.0. A fixture captured on the unmodified tree pins this. No committed test covers attribution outputs byte for byte today, so the fixture is new.
- **Validation: `scripts/adjusted_attribution_check.py`.**
  - **Construction.** It runs at the case study's configuration, using the construction of `case_study_type_i_check`: the astrocyte `cohort_30` table (n = 80), methylation + expression blocks, a covariate `z` confounded with group × stage and injected into a fixed feature subset, five noise covariates (six design columns), rank floor 2, net displacement and `df`. The stage effect is the Phase 6 one (δ = 2.0), so the planted orientation has a trajectory to rotate. Seeds come from a fresh family (base seed 2033).
  - **Cells.** The anchor, and orientation at e = 0.5 and 1.0, each at 500 replicates × 999 permutations with 200 bootstrap replicates.
  - **Arms.** Each unit is analysed once by `analyze_trajectories`, with covariates and attribution. That is the `adjusted` arm. Its frozen model and features then feed two more attributions:
    - `oracle`: the injection removed exactly in the analysed coordinates;
    - `unadjusted`: the existing path, descriptive.
  - **Predeclared decision.** `adjusted_attribution_decision.json` returns `adjusted_attribution_validated`, `diagnose_adjusted_attribution` or `construction_uninformative`. It reads four criteria on the observed principal-orientation contrast:
    - contamination of the top-k drivers by injected features, against the oracle;
    - recovery of the generator's truth drivers, non-inferior to the oracle;
    - residual confounding along the injection direction, against the unadjusted arm;
    - bootstrap sign stability and top-k frequency, non-inferior to the oracle.
  - **Calibration.** A calibration step on the injection strength comes first, with every iteration recorded.
- **Run, report, release.** The sequence is:
  - a local smoke run and timing, then an approval gate;
  - the paper-grade run, of about 1,500 units at roughly 15 s each, so about 6–7 core-hours. The smoke timing decides between local workers and `ing` (design D11).
  - a dated report and `results/adjusted-attribution-<date>/`;
  - updates to `docs/roadmap.md` (item 11), `docs/api/analysis.md`, `AGENTS.md` and the version (0.14.0).
- **Out of scope.**
  - The plasmode demonstration on real SEA-AD blocks, which is item 12 and has its own proposal.
  - Attribution on a latent-adjusted basis.
  - Covariate interactions with group or stage.
  - Changes to the simulation harness's attribution, which has no covariates.
  - Pathway or module enrichment, which stays in the case-study repositories.
  - Any change in `sead_sexdiff` or `adni_sexdiff`.

## Capabilities

### New Capabilities

- `adjusted-attribution-validation`: the semi-synthetic check that covariate-adjusted attribution removes a confounded covariate's contribution, and recovers and stabilizes planted orientation drivers as well as an oracle that knows the injection. It covers:
  - the construction, its calibration, and the three arms with their consistency checks;
  - the metrics, records and merge;
  - the predeclared decision and its verdicts.

### Modified Capabilities

- `orientation-driver-attribution`: attribution accepts declared covariates. It then computes the cell means as covariate-adjusted least-squares means at a reference row and runs a bootstrap that refits the adjustment within each resample, and it records both. Without covariates its outputs are unchanged.
- `real-data-trajectory-analysis`: attribution is allowed with declared covariates and still rejected with latent adjustment. Its means project to the analysis's stage means. The result exposes the fitted estimator and the features in memory. `motco analyze` prints the attribution's mean source.

## Impact

- **Code:**
  - `src/motco/stats/attribution.py`: the `covariates` and `covariate_reference` keywords, adjusted means, the adjusted bootstrap, and the configuration and interpretation rows.
  - `src/motco/simulations/attribution_diagnostics.py`: pass-through and the `adjustment` block.
  - `src/motco/analysis.py`: the guards, the attribution call and `TrajectoryAnalysisResult.latent_artifacts`.
  - `src/motco/cli.py`: the summary line.
  - `src/motco/__init__.py`: version 0.14.0.
  - New `scripts/adjusted_attribution_check.py`, which imports `case_study_type_i_check` and `covariate_type_i_check`. `covariate_type_i_check.build_replicate` additionally returns the generated dataset for the truth drivers, and its records are unchanged.
  - A new `scripts/adjusted_attribution_array.sbatch` only if the run goes to `ing`.
- **Tests:**
  - new `tests/test_adjusted_attribution.py` (primitive and analysis);
  - new `tests/test_adjusted_attribution_check.py`;
  - CLI cases in `tests/test_analysis_cli.py`;
  - a new fixture `tests/data/attribution_no_covariate_identity.json` with its capture script;
  - the existing `tests/test_attribution.py`, `tests/test_attribution_diagnostics.py`, `tests/test_evaluation_attribution.py`, `tests/test_analysis_identity.py` and `tests/test_analysis.py` must pass unchanged, except the one test that pins the rejection.
- **API:** additive. A request for attribution with covariates used to be an error, and now returns a result. The rejection test (`tests/test_analysis.py::test_attribution_with_covariates_is_rejected_before_fitting`) is replaced, not deleted.
- **Results and docs:**
  - `results/adjusted-attribution-<date>/`: settings, calibration, PROVENANCE, metric tables and decision;
  - `docs/reports/adjusted-attribution-<date>.md`;
  - `docs/roadmap.md`, `docs/api/analysis.md` and `AGENTS.md`.
- **Downstream:**
  - Item 12's proposal starts after the release and pins 0.14.0.
  - Under `diagnose_adjusted_attribution` the withheld leg stays withheld, and a diagnosis is proposed before item 12.
- **Compute:** about 6–7 core-hours, re-timed at smoke.
