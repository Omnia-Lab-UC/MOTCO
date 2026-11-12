# Proposal

## Why

Production uses the CV-selected rank A of the stage-PLS fit for two jobs: it regularizes the fit, and it is the number of dimensions the trajectory statistics are measured in. Three results say the second job costs power on what the span sees: ranks above K − 1 lose `shape` power at the Phase 5 design (rank ladder, 2026-09-08); two near-empty added components halve magnitude `delta` (latent-space comparison, 0.74 → 0.37 at n = 200); and at δ = 2.0 every candidate of the production latent-space selection gained where it chose rank 2 and production's CV chose 3 (`angle` 0.46 → 0.90–1.00 on the SEA-AD table's 22.5° rotation).

Roadmap open item 1 named two ways to measure fewer dimensions: a stricter rank rule, and the predictive subspace of the production fit. An exploratory probe on 2026-11-12 (`notes.md`; not a result) separated them:
- **Lowering the fit rank is the wrong lever.** The later components aim the K − 1 stage directions. A forced rank of K − 1 and a tolerance rule gained `angle` only at δ = 2.0, and at δ = 0.25 a forced rank of 2 kept 0.31 of the pooled stage program against production's 0.67–0.77. A one-standard-error rule equals production wherever the held-out AUROC is at its ceiling, because the best candidate's fold SE is exactly zero.
- **The predictive subspace keeps the fit and drops the dimensions.** The stage labels are one-hot, so the Y loadings have rank K − 1 whatever A is, and the fitted stage labels move only along K − 1 directions of the score space. Measuring there gained on all three target statistics at δ = 2.0 (`delta` 0.23 → 0.53, `angle` 0.45 → 1.00, `shape` 0.68 → 0.88 on the SEA-AD table) and more at δ = 0.25, with no paired loss.
- **It has costs the probe could not size.** Off-target readings depend on the metric chosen on the subspace, one fit component too many leaks a pure rotation into `shape`, the weak design's `shape` control read 0.077 over 300 units, and nothing gained at the weak design.

The predictive subspace changes what production measures, so it must be read against production on matched data under a predeclared decision before any default moves. The probe's datasets are spent, and so are the production selection's δ = 2.0 cells, where the gains were first read.

## What Changes

- **Opt-in predictive-subspace coordinates in the harness.** `SimulationEvaluationParams.latent_coordinates` gains three values. Each measures the final model's scores in the row space of its Y loadings (dimension K − 1, or the fit rank when that is smaller):
  - `"predictive"`: the scores projected on an orthonormal basis of that row space. It equals the score measurement up to an orthogonal map whenever the fit rank is at most K − 1.
  - `"predictive_orthonormal"`: the model's input coordinates projected on an orthonormal basis of the same subspace in feature space.
  - `"predictive_fitted"`: the fitted stage indicators in the model's response units.

  Nothing else moves: the fitted model, the cross-validation, the selected or forced rank, the rank floor and the permutation draws are those of the score measurement. The three share one feature-space subspace, so their span-capture shares are equal. The values are rejected before any fit with a method other than `pls`, with attribution, and with fewer than three stages. While the setting is absent every committed identifier, signature, record and report is unchanged.
- **A prebuilt-fit path.** A predictive-subspace integration is built from an existing score-coordinate `pls` integration without refitting, equal to fitting it from scratch, so a validation unit fits one model and measures it several times.
- **The same opt-in setting in `motco.analysis`.** `TrajectoryAnalysisParams.latent_coordinates` takes the three predictive values, validated before fitting, omitted from `to_dict()` while unset and recorded when set. It is rejected with attribution and with the split interaction test. The plasmode profile needs it; it gets no CLI flag in this change.
- **A paired pilot.** `scripts/predictive_subspace_check.py` (`run --profile simulation|weak_size|plasmode` | `merge --profile …` | `decide`; sbatch array `scripts/predictive_subspace_array.sbatch`) fits production once per unit and reads six arms on the same dataset and RRPP seed:
  - `pls_cv`: production;
  - `predictive`, `predictive_orthonormal`, `predictive_fitted`: the candidates, built from the production fit;
  - `forced_k_minus_1`: the fit forced to K − 1 components, the control that lowers the fit rank;
  - `predictive_rank_plus_one`: the `predictive` coordinates of the fit forced to A + 1 components, a descriptive reading of how near the CV rank sits to over-aiming (simulation only).
- **Three profiles.**
  - **`simulation`**: a fresh seed family (base seed 2036) on the latent-space comparison's design generators and seed derivations, under the planted scale, at 100 × 199. Seven designs: n = 200 and the SEA-AD table at δ = 2.0; n = 200 at δ ∈ {1.0, 0.5, 0.25} and the SEA-AD table at δ = 0.5 (the stage-effect band, where the CV rank inflates and the probe's candidates separate); and the SEA-AD table at the committed weak δ = 0.11. Eleven cells each: anchor, translation, magnitude at e ∈ {0.05, 0.1, 0.25}, orientation at e ∈ {0.1, 0.25, 0.5} (φ = 0) and shape at e ∈ {0.1, 0.25, 0.5}. One descriptive cell adds the Phase 6 feature-set orientation surgery at e = 1.0 on the SEA-AD table, where extra components are the only route by which the plant is seen. 7,800 units.
  - **`weak_size`**: the comparison's `seaad_weak` anchor and translation at replicates 0–499 and 999 permutations (base seed 2035). The `pls_cv` arm must reproduce the span-split check's committed `production_table.csv` counts exactly. 1,000 units.
  - **`plasmode`**: the SEA-AD astrocyte export through the plasmode's contract, at the case-study configuration and 999 permutations: the anchor at the `paper` stream's replicates 0–499, and the production selection's real-noise span plant at e ∈ {0.25, 0.5, 1.0} at replicates 0–99. The `pls_cv` arm must reproduce the plasmode's committed anchor counts over 500, the production selection's committed plant rows, and the label-free fidelity values (rank 5, AUROC 0.635556). 800 units.
- **A predeclared decision.** `predictive_subspace_decision.json` (rule version 1) judges each candidate arm against `pls_cv`:
  - **S1 size:** the controls' `delta` / `angle` / `shape` in three families under the stage-span paper grade's count-plus-ceiling rule (`judge_rates`, imported unchanged): the simulation's anchors and translations (m = 42), `weak_size` (m = 6) and the plasmode anchor (m = 3). All three must hold.
  - **S2 protected power:** the 54 target cells of the six designs at δ ≥ 0.25 and the plasmode's three plant `delta` cells (m = 57), by a one-sided exact paired loss test under the count-plus-ceiling rule.
  - **S3 specificity:** the off-target rates of the target cells. At the δ = 2.0 designs (m = 36), where production's specificity is established, the bound rule. In the band (m = 72), where production's own off-target rates have never been measured, a one-sided exact paired excess test against `pls_cv` under the count-plus-ceiling rule.
  - **G gain:** among the 18 target cells of the δ = 2.0 designs, at least one one-sided exact paired gain p-value at or below 0.05/18. The band's and the plasmode's gains are recorded with the same test and do not decide.

  Per arm: `size_excess` | `power_loss` | `specificity_excess` | `no_gain` | `qualifies`, in that order. Overall: `propose_paper_grade` naming one qualifying arm (the one with the fewest off-target readings above the bound; ties in the order `predictive`, `predictive_orthonormal`, `predictive_fitted`) | `keep_production`; `reproduction_failed` gates everything. Descriptive beside it: the selected rank, held-out AUROC and capture per cell; each arm's observed statistic beside the population's through the same map; the control and rank-plus-one arms' rates; production's own off-target rates in the band; and the surgery cell.
- **Documentation and records.** A dated report `docs/reports/predictive-subspace-<date>.md`; aggregates under `results/predictive-subspace-<date>/`; a section of `docs/measurement-visibility.md` on the fit rank against the measured dimension; `docs/api/`; roadmap open item 1; `AGENTS.md`; `CHANGELOG.md`; MOTCO 0.24.0.

**Non-goals.**
- **The paper grade and any default flip.** On `propose_paper_grade` a separate change runs the named arm at 500 × 999, and only its verdict may flip the harness and `motco analyze` defaults, with pins for every committed configuration and the SEA-AD re-issue.
- **A rank rule.** The candidates keep production's CV. The stricter rules the probe read are not arms; a rule for the fit rank under a fixed measured dimension is a question only if an arm qualifies.
- **`motco analyze --latent-coordinates`.** The CLI flag follows a verdict; until then the setting is a typed field the pilot uses.
- **Attribution and the split interaction test under predictive coordinates.** Both are defined on the fitted span, and their counterpart on a K − 1-dimensional measurement needs its own definition.
- **Other regularizers of the same subspace** (shrinkage LDA, roadmap open item 4) and the oblique projection whose dropped scores are uncorrelated with the kept ones, which had the largest off-target readings in the probe.
- **Out-of-span questions**, which the split interaction test and the scan own.

## Capabilities

### New Capabilities
- `predictive-subspace-validation`: the paired pilot of the predictive-subspace coordinates against the production measurement. It covers the six arms on one fit, the three profiles and their reproduction of the committed parents, the unit record and its consistency checks, the predeclared decision with its descriptive tables, and the dated report.

### Modified Capabilities
- `simulation-evaluation-harness`: `latent_coordinates` gains the three predictive values with their definition, their rejection rules and the prebuilt-fit path; the single-pass requirement is restated to cover them.
- `real-data-trajectory-analysis`: the analysis parameters gain the opt-in `latent_coordinates` setting with its validation, rejections, record and invisibility while absent.
- `latent-space-metric-compatibility`: the documentation states that the production measurement uses the fit rank as its dimension, what the predictive subspace keeps and drops, and the pilot's verdict.

## Impact

- **Code:**
  - `src/motco/simulations/evaluation.py`: the three values, their validation, the coordinate step after the final fit, and the prebuilt-fit helper;
  - `src/motco/simulations/diagnostics.py`: a `LinearLatentMap` constructor for the predictive subspace;
  - `src/motco/analysis.py`: the typed field, its validation and record;
  - `src/motco/__init__.py` and `CHANGELOG.md`: 0.24.0.

  `motco.stats`, `motco.cli` and `grid._SIGNATURE_OMIT_WHEN_NONE` (which already omits `latent_coordinates`) are not touched.
- **Scripts:** `scripts/predictive_subspace_check.py` and `scripts/predictive_subspace_array.sbatch`. They import unchanged `latent_space_comparison` (designs, seeds, `unit_dataset`), `span_split_check` (the weak units), `production_latent_space_check` and `plasmode_check` (the export contract, label permutation, span plant, fidelity values, paired tests) and `stage_span_paper_grade_decision` (`judge_rates`).
- **Tests:** the predictive matrices span one subspace of the fitted span and have the rank of the Y loadings; `predictive` equals the scores up to an orthogonal map at a fit rank of K − 1; the fitted labels are preserved; statistics are invariant to the orthonormal basis; rejections; byte identity while absent in both entry points; prebuilt-fit equality; harness–analysis equivalence under the setting; driver units, seeds, parents, resume, merge, every reproduction outcome and every decision path.
- **Compute:** about 9,600 units with one production CV each. The estimate before the smoke is 100–140 core-hours on `ing` (simulation about 65, `weak_size` about 11, plasmode about 33). A local smoke replaces the estimate, and launch waits on the user's approval with a cap.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged. The setting is absent by default in both entry points, and the parents' committed files are read, hashed and refused by name if changed.
