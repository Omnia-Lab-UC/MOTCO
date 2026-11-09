# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **The models and draws already exist.**
  - `analyze_trajectories` builds the full (group × stage) and reduced (group + stage) models with covariates appended to both (`_covariate_adjusted_design`).
  - RRPP permutes the reduced residuals by Freedman–Lane, drawing `rng.permutation(n)` per draw from `np.random.default_rng(seed)` on its serial path.
  - `residual_correction_scale` gives the `df` and `leverage` multipliers.
  - The split test (`motco.stats.interaction`) reuses all of this and holds the column-space helpers (`_column_basis`, `_interaction_basis`) and the input validation a per-feature scan needs.
- **The features.** `analyze_trajectories` builds `X`, the per-block standardized features concatenated in declared order with columns `block__feature`. A per-feature F is invariant to each feature's location and scale, so the scan's results do not depend on the standardization. They do depend on the declared unit transform (for example B to M values), which is nonlinear.
- **The parent run.** The split-test validation (`results/span-split-2026-11-05/`, 0 failures) reran two older runs' units and reproduced both exactly.
  - Its simulation profile: 6,600 units at a median of 11.1 s, 25.6 core-hours, 498 features (methylation + expression).
  - Its plasmode profile: 1,500 units at a median of 80 s, 36.9 core-hours, 10,000 features, n = 80, six covariates, `df`.
  - Both profiles' settings, `split_table.csv` and trajectory or production tables are committed.
- **The plasmode's planted set is known per unit.** `plasmode_check.build_unit` returns the 20 planted feature names. The labels are permuted within stage, so every feature outside the planted set is null by randomization.
- **The latest release is 0.19.0.** Its `latent_coordinates` setting is opt-in and byte-identical when absent, so the parent's units still rebuild under it.

## Goals / Non-Goals

**Goals:**
- A per-feature interaction test whose draws, models and corrections are the analysis's, with a global p-value and familywise-controlled per-feature p-values.
- A cost of seconds per analysis at 10,000 features and 999 draws.
- A validation whose only new compute is the scan itself plus a rerun of committed units, paired unit by unit with the split test and the trajectory statistics.

**Non-Goals:**
- Choosing between the scan, the split test and the trajectory statistics for a case study. Their joint multiplicity policy belongs to the case study.
- Effect-size estimation beyond the per-stage group differences in the table.
- Any change to the split test's results or messages.

## Decisions

**D1. The per-feature statistic is the interaction F of one column.**
- **Definition.** For column y: SS_int = ‖(H_F − H_R)·y‖², RSS = ‖(I − H_F)·y‖², F = (SS_int/df_int)/(RSS/df_res). This is the split test's pseudo-F restricted to one coordinate, and the per-feature sums add up to the split test's `total` part.
- **Why a ratio.** It is the reason D1 of the split test gave: the Freedman–Lane residuals' scale shrinks by about n/(n − r_R), and a ratio cancels that. Per feature this also makes every column's null approximately F(df_int, df_res), whatever its scale, which is what makes a maximum over features meaningful.
- **Alternatives.**
  - **The per-feature SS.** It is not comparable across features of different variance and carries the residual shrinkage.
  - **A moderated F** (limma-style variance shrinkage). It is likely more powerful at n = 80, but its hyperparameters would have to be re-estimated inside every draw, and its permutation validity would need its own check. It is a follow-up if the measured power is short.
  - **A per-stage t contrast.** It needs a predeclared stage; the interaction F tests every stage pattern at once.

**D2. Features are combined by the maximum, with Westfall–Young step-down max-T.**
- **Why the maximum.** A sum of per-feature energy has noise growing like √p, so k planted features must carry energy of order √p to stand out. A maximum's threshold grows only like log p: the F bar at df (2, 68) rises from 3.1 for one feature to 14.7 at a Bonferroni level over 10,000. The back-of-envelope calculation in the proposal (planted-feature power 0.06 / 0.51 at e = 2.0 / 3.0) uses that bar.
- **Why step-down.** It gives per-feature adjusted p-values with strong familywise control under subset pivotality. It is uniformly at least as powerful as single-step max-T, and its top feature's adjusted p-value is the single-step global p-value. Permutation keeps the correlation between features, so the threshold is less conservative than Bonferroni.
- **The rule.** Order the non-degenerate features by observed F (descending, ties by column index). For draw b and rank j, u_bj = max over ranks ≥ j of F*_b. Then p̃_j = (1 + #{b: u_bj ≥ F_(j)})/(1 + B), and the adjusted p-value is the running maximum down the order.
- **Alternatives rejected on resolution.**
  - **minP** (combining per-feature permutation p-values) is degenerate at B = 999 and p = 10,000. Per-feature permutation p-values bottom out at 0.001, and the minimum over 10,000 features sits at that floor in nearly every draw. It needs B ≫ p/α = 200,000.
  - **Benjamini–Hochberg on permutation p-values** fails for the same reason: its first threshold is 0.05/10,000 = 5·10⁻⁶, below the 1/1000 resolution.
  - **BH on parametric F p-values** relies on F tails on real counts that nothing here validates.
- **Other alternatives.** Higher criticism or a top-k sum targets intermediate sparsity with a tuning choice. It is out of scope.
- **The hit level is fixed at 0.05** (`INTERACTION_SCAN_ALPHA`), the level the validation measures. The per-feature table carries every adjusted p-value, so a reader can apply a stricter level.

**D3. The permutation costs O(n · rank(F) · p) per draw and never forms a permuted matrix.**
- **The reduction.** Let W = [U_int, U_R] be an orthonormal basis of col(F), interaction directions first, as in the split test. For permuted residuals E[π]:
  - SS*_f = ‖U_intᵀ·E[π]_f‖²;
  - RSS*_f = ‖E_f‖² − ‖Wᵀ·E[π]_f‖², since a column's norm is invariant to row permutation.
- **The computation.** With W' = W[argsort(π)], Wᵀ·E[π] = W'ᵀ·E. Each draw is one (rank(F) × n)·(n × p) product. At n = 80, rank 12 and p = 10,000 that is about 10⁷ multiply-adds, so 999 draws take seconds.
- **Step-down counting.** The features are ordered once. Per draw, the suffix maximum of F* in that order is compared with the ordered observed F, and the counts accumulate. That is O(p) per draw.
- **The observed F** goes through the same code with the identity permutation.
- **Corrections** scale E's rows by s before any draw, as in the split test. Under `df`, s is constant, so every F* is unchanged.
- **Returned nulls.** The null maxima and their columns are always kept (2·B numbers). The B × p null F is kept only on request: it is 80 MB at 10,000 features, and tests need it, not the analysis.
- **Degenerate features.** A feature is degenerate when its full-model RSS is at most 10⁻¹² of its centred sum of squares, or when it is constant. Such features are excluded before ordering. RSS* is clipped at 0, as in the split test.

**D4. Where it lives.**
- **`motco.stats.interaction.interaction_feature_scan`**, with an `InteractionScanResult` type. It shares the module's validation and column-space helpers; the shared validation is factored out so that the split test's messages and results stay byte-identical (its tests are the check). It imports nothing from `motco.simulations`.
- **`analyze_trajectories`** calls it on `X` with the covariate-adjusted models.
- **The validation driver** calls it on the split test's own inputs.
- **The simulation harness and study configs** are left out.

**D5. The analysis record and table.**
- **Inputs.** `X`, the covariate-adjusted full and reduced models, the seed and permutations. The residual correction applies only when permutations run, as for RRPP and the split test.
- **The `interaction_scan` block** comes after `measurement` and any `span_split` block. It holds:
  - `n_features`, `n_degenerate`, `df_interaction`, `df_residual`;
  - `max_f` with its `feature` and `block`;
  - `p_value_global`;
  - `hit_alpha` (0.05) and `hits`: block, feature, F and adjusted p, ordered by F;
  - per block: features, degenerate features, hits, and the share of null maxima;
  - `permutations`, `seed`, `residual_correction` and `draws_match_rrpp` (true only with one job).
- **The per-feature table** (`TrajectoryAnalysisResult.interaction_scan`, written as `interaction_scan.csv`) has the columns `block`, `feature`, `ss_interaction`, `rss`, `f`, `p_value_adjusted`, `degenerate`, and `difference[<stage>]` for each declared stage.
  - The differences are (LS row of the second group at stage s − LS row of the first) · β, with β = `estimate_betas(model_full, X)`. They use the same LS-mean rows (covariate reference row included) that give `stage_means.csv` in latent coordinates, here in standardized feature units.
  - The feature names are split from `block__feature` as `pls_loadings.csv` does.
- **The resolved parameters** carry `interaction_scan` only when it is true. Without it, the record and every file are byte-identical to 0.19.0's.
- **Allowed with latent adjustment and attribution.** The scan uses neither the latent basis nor the residualized features: `X` and the RRPP models are the same under latent adjustment, and a test pins that.
- **The CLI** gains `--interaction-scan` and one printed line: the maximum F and its feature, the global p-value, and the hits at 0.05.

**D6. The validation reruns the split-test validation's units.**
- **Why that parent.** Its units already pair the production statistics and all three split parts on the same datasets and draws, and it reproduced both of its own parents exactly. Reproducing it ties the scan to the latent-space comparison and the plasmode in one step.
- **Simulation units.**
  - The settings are rebuilt through `span_split_check.simulation_settings` and must equal the committed `results/span-split-2026-11-05/simulation/settings.json`, whose sha256 is recorded.
  - Each unit composes the parent's own steps (`lsc.unit_dataset`, `production_params`, integrate, evaluate, `production_record`, `split_inputs`, `split_record`), then runs the scan on `split_inputs`' features and models at the unit's RRPP seed. The features are the fit's input coordinates, a per-feature affine map of `X`, which the scan ignores.
- **Plasmode units.**
  - The settings come through `span_split_check.plasmode_settings`, which must equal the committed `plasmode/settings.json`. The export hashes must match the plasmode's settings, and the fidelity gate must pass.
  - Each unit is `plasmode_check.build_unit`, then one `pc.analyse` with the parent's parameters plus `interaction_scan=True`, then the consistency and anchor latent checks.
  - The planted counts are computed from `build_unit`'s planted names against the scan table and recorded as counts only. The parent's uncorrected split rerun (the SS diagnostic) is not repeated.
- **Reproduction.** Per cell (and design), the split test's pseudo-F rows for the three parts and the trajectory statistics' rows must equal the parent's committed `split_table.csv` and `production_table.csv` / `trajectory_table.csv`, in count and rejections, under the parent's rule (p < α in the simulation, p ≤ α in the plasmode). A failure gates the decision (`reproduction_failed`).
- **The scan's own rejection rule** is adjusted p ≤ 0.05 in both profiles. With B = 999 that is exact: P(p ≤ 0.05) = 50/1000 under exchangeability.

**D7. The criteria.** They use the stage-span paper grade's `judge_rates` by import.
- **S1, simulation size (m = 6).** The global rate at the anchor and translation of `n200`, `seaad` and `seaad_weak` at R = 500. Translation is a constant offset, so it has no interaction and is a null for the scan.
- **S2, real-noise size (m = 1).** The global rate at the plasmode anchor.
- **S3, strong familywise control (m = 2).**
  - **What is measured.** At e = 2.0 and 3.0, the rate of units with any hit outside the planted set.
  - **Why it is needed.** Strong control is what licenses naming features, and step-down needs subset pivotality, which Freedman–Lane with covariates satisfies only approximately. The plasmode is the one place with real noise where the null set is known.
- **G, sparse gain.**
  - **The comparison.** At e = 3.0 the scan's global rate must exceed the largest of `delta`, `angle`, `shape` and the three split parts' rates on the same units by more than 2 pooled SE.
  - **Why e = 3.0.** It is the plasmode's predeclared cap, the cell where the gap was documented. The comparator is the best existing test, chosen after the fact among six, which is the conservative choice. e = 2.0 is descriptive.
- **Precedence:** `reproduction_failed`, `size_excess`, `fwer_excess`, `no_sparse_gain`, `scan_validated`.
- **The consequence states that no verdict reopens Phase 6.** The plasmode report already said that a new test needs its own predeclared rule, not a rescue of the old one.

**D8. Descriptive tables.**
- **Plasmode recovery.** Recall |hits ∩ K|/20, hit precision |hits ∩ K|/|hits| over units with hits, and top-20 precision |top-20 by F ∩ K|/20. They sit beside the plasmode's committed attribution planted-set precision (`plasmode_summary.csv`: adjusted arm 0.344, residual component 0.432 at e = 3.0), the same top-20 metric.
- **Block balance.** Under per-feature pivotality each non-degenerate feature is equally likely to give a draw's maximum. A block's share of null maxima should therefore equal its share of features (0.5 / 0.5 on the export). Correlated features lower a block's share and heavy tails raise it.
- **Simulation localization.** The scan beside the split parts and the production statistics, per design and cell, on the same units. The prediction: the scan trails the complement on the dense φ = 90 rotation, which spreads over hundreds of features.
- **Simulation false hits.**
  - **The null set.** A feature is population-null when its population B − A difference, from `dataset.population_trajectories` in integration units, is constant across stages within 10⁻⁹ of its largest absolute population mean. That test is invariant to the per-feature affine standardization.
  - **What is reported.** The rate of units with a hit on a population-null feature. It is reported where the cell has population-null features, and it is the generator-noise analogue of S3.

**D9. Documentation.**
- **The measurement-visibility page** gains a section on the scan: what a hit and a non-rejection mean, why the maximum, how it differs from attribution and the split test, and what it does not provide. It is filled from the report once the decision exists.
- **API pages:** `docs/api/sd.md` documents the primitive and `docs/api/analysis.md` the option, the block and the CSV.
- **Roadmap and repository guide:** roadmap item 16's thread 3 records the verdict, and `AGENTS.md` gains the function, the option and the driver.

## Risks / Trade-offs

- **[Subset pivotality is approximate under Freedman–Lane with covariates]** → S3 measures strong control on real noise. A failure is the verdict `fwer_excess`, and only the global test could then be reported.
- **[Heavy-tailed or discrete features could dominate the maximum]** → The block-balance table shows it. minP, the usual remedy, is infeasible at B = 999 (D2). Per-block maxima are a follow-up if the balance is materially off.
- **[Correlated features give clustered hits]** → The documentation states that a hit names a feature whose interaction is detected. Its correlated neighbours may be hits too, and a hit is not a causal driver.
- **[The plant is stylized: 20 features, 3 SDs, one group]** → G is a statement about the documented gap, not realistic power. e = 2.0 is reported, and a non-rejection is read as uninformative without power.
- **[The plain F is noisy at n = 80]** → If the measured power at e = 2.0 is low, a moderated statistic is the named follow-up.
- **[The simulation has 498 features, not 10,000]** → S2 and S3 on the real export cover the high-dimensional case, and S1 covers generator noise and the SEA-AD-like cell table.
- **[Size on the export is conditional on one data matrix]** → As in the plasmode; S1 gives unconditional size on generated data.
- **[Parent drift across 0.18.0 → 0.20.0]** → Reproduction gates the decision. The settings and parent tables are bound by sha256.
- **[More than one job changes RRPP's draws but not the scan's]** → `draws_match_rrpp` records it, and the validation runs one job per unit.
- **[Compute estimate]** → The local smoke times both profiles, and the projection in `notes.md` replaces the estimate before launch.

## Migration Plan

Opt-in throughout, so there is nothing to migrate.
- **Version:** MOTCO 0.20.0.
- **Rollback:** without the analysis setting every record and file equals 0.19.0's; the stats function is new.

## Open Questions

- **Multiplicity in a case study.** How a case study reads the scan beside the split test and the trajectory statistics (a fixed primary, Holm across the families, or the omnibus first). It belongs to `sead_sexdiff`'s change, before any real-label analysis, and changes nothing here.
- **A moderated statistic or per-block maxima.** Whether either is worth its own change depends on this run's power at e = 2.0 and its block balance.
