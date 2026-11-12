# Design

## Context

See proposal.md (Why) for the motivation, and `notes.md` for the probe's design and tables. The facts below shape the approach.

- **The algebra, checked numerically** (2026-11-12, sklearn 1.8.0; n = 80, p = 300, three stages, ranks 2 and 5):
  - `PLSRegression` fits `Y ≈ T·Qᵀ` with `T = x_scores_ = X_in·R`, `R = x_rotations_` and `Q = y_loadings_` (K × A). `predict` equals `T·Qᵀ` rescaled by the model's response scale and mean (1e-15).
  - The stage labels are one-hot, so the centred response has rank K − 1 and so has `Q` when A ≥ K − 1: its third singular value is exactly 0 at rank 5. `Q` had rank 2 on every one of the probe's 1,080 production fits.
  - Let `V` (A × rank(Q)) be an orthonormal basis of the row space of `Q`. Then `T·V·Vᵀ·Qᵀ = T·Qᵀ` (1e-15): the part of a score row orthogonal to that row space does not move the fitted labels.
  - The columns of `R·Qᵀ` lie in the span of `R·V` (1e-16), so the projected scores, the fitted labels and the input coordinates on `qr(R·V)` are three coordinate systems on one feature-space subspace.
  - At a fit rank of at most K − 1, `V` is square and orthogonal.
- **What exists.** `fit_pls_latent_space` (`evaluation.py`) is the core shared by the harness and `motco.analysis`. It takes `coordinates` (`"scores"` | `"orthonormal"`), replaces the matrix after the final fit and returns a `LinearLatentMap` whose `transform` measures in the same coordinates. `SimulationEvaluationParams.latent_coordinates` is hash-omitted while `None`, rejected with attribution and with non-`pls` methods, and `pls_orthonormal_from_pls` builds the alternative integration from a prebuilt fit. `motco.analysis` passes no coordinates.
- **What the earlier runs measured about coordinates.** Orthonormal coordinates of the full span lost power in 19 of 45 target cells against the scores (`power_loss`, 2026-11-09). In the probe the three coordinate systems on the predictive subspace had equal target power at every cell and differed only off target.
- **The parents.**
  - The latent-space comparison's `design_generators`, `generator_seed`, `rrpp_seed`, `unit_dataset` and `cell_spec` rebuild any design and cell from a settings dictionary and a base seed.
  - The span-split check committed `pls_cv` counts at 999 permutations for the comparison's anchor and translation at replicates 0–499 (`results/span-split-2026-11-05/simulation/production_table.csv`).
  - The plasmode committed its anchor rates over 500 replicates (`results/plasmode-2026-11-02/plasmode_summary.csv`: `delta` 0.042, `angle` 0.048, `shape` 0.026, rank 5 on every unit). The production selection committed the `pls_cv` rows of its real-noise span plant at replicates 0–99 (`results/production-latent-space-2026-11-10/plasmode/rate_table.csv`) and showed that its anchor p-values equal the plasmode's on the same replicates.
- **Costs.** The probe ran a production CV unit with eleven RRPP passes at 199 permutations in 19–27 s on the SEA-AD table and about 35 s at n = 200 and δ = 0.25. The production selection's plasmode `pls_cv` arm took 134 s per unit at 999 permutations. A predictive arm adds one RRPP pass and no fit.

## Goals / Non-Goals

**Goals:**
- Measure the trajectory statistics in the K − 1 stage-predictive directions of the production fit, as opt-in values of the existing coordinate setting that change nothing while absent.
- Read the three coordinate systems against production on matched data, fit and draws, on fresh datasets, at the stage effects where the probe saw the candidates separate and at the stage signal where real data sit.
- Decide, under a rule fixed before any unit, whether a paper grade is worth its cost and for which coordinate system.

**Non-Goals:**
- Changing the fit: the CV, the rank rule, the rank floor and the scaling stay as they are.
- Resolving rates near the bound at R = 100 in the simulation profile. The controls the probe flagged are read at R = 500.
- Four-stage designs. At the Phase 5 design the CV selects K − 1 (range 2–3), where `predictive` is the score measurement; the paper grade carries those columns.

## Decisions

**D1. The predictive subspace is the row space of the Y loadings, in score coordinates.**
- **Definition.** With the final model's `Q` and `V` as in Context, the subspace is `span(V)` in score space and `span(R·V)` in the model's input space. Its dimension is the numerical rank of `Q` (singular values above 1e-10 of the largest): K − 1 when the fit rank is at least K − 1, the fit rank otherwise.
- **Why this subspace.** It is the largest part of the fitted span the stage fit uses and the smallest that keeps the fitted labels exactly. It needs no new fit, no new tuning parameter and no group label.
- **Why not a lower fit rank.** A forced rank of K − 1 truncates the fit and loses the components that aim the directions (proposal, Why). It stays in the pilot as the control arm.
- **A degenerate fit.** If `Q` has rank below min(A, K − 1), the step raises naming the rank; the probe met none.

**D2. Three coordinate systems on that subspace, each a value of `latent_coordinates`.**
- **`predictive`: `T·V`.** The score metric restricted to the subspace. At a fit rank of at most K − 1 it is the score matrix times an orthogonal matrix, so every statistic and p-value equals production's. It is the candidate nearest to production.
- **`predictive_orthonormal`: `X_in·qr(R·V)`.** The standardized feature space's own lengths and angles within the subspace, the geometry the stage-span constructions define. Every orthonormal basis gives the same statistics (the orthonormal change's D2).
- **`predictive_fitted`: `T·Qᵀ`.** The fitted stage indicators in the model's response units: K columns of rank K − 1. Its metric depends on the response scaling, so on the stage prevalences.
- **Why all three.** The probe could not rank them: target power was equal at every cell, and the off-target rates under a 45° rotation at δ ≤ 1.0 read `shape` 0.50–0.80 for `predictive`, 0.00–0.25 for `predictive_orthonormal` and 0.00 for `predictive_fitted` on 20 units each. An arm is one RRPP pass, so the pilot measures all three and the decision names one.
- **Alternative: a separate `latent_subspace` parameter crossed with `latent_coordinates`.** It separates the subspace from its metric, but the fitted-label metric exists only on the predictive subspace, and a second hash-omitted parameter would need its own pins at a default flip. One setting names what the statistics are measured in.
- **Alternative: the oblique projection whose dropped scores are uncorrelated with the kept ones** (the orthogonal-PLS reading). It had the largest off-target rates in the probe (`shape` 0.80–1.00 in the same cells) and is not an arm.

**D3. Validity carries over by linearity; size is still read.**
- **The argument.** The fit is group-blind and fixed. Each predictive matrix is the score matrix times a fixed matrix, the LS-mean fit is linear, and RRPP's permutation indices do not depend on the outcome's columns. A predictive arm's draws are therefore production's draws times that matrix, and the exchangeability argument is production's.
- **Why read size anyway.** The statistics are not linear, and the subspace is where the fit's in-sample optimism concentrates. The probe's weak-design `shape` control read 0.077 over 300 units against production's 0.050.

**D4. The setting, its rejections and the analysis field.**
- **Harness.** The three values join the allowed values of `SimulationEvaluationParams.latent_coordinates`; `fit_pls_latent_space(coordinates=)` takes them. The metadata records the value, the latent matrix has rank(Q) columns (K for `predictive_fitted`), and the returned `LinearLatentMap` has basis `qr(R·V)` and a `transform` in the measured coordinates, so `population_latent` and span capture follow.
- **Composition.** The step is a column map of the final model's scores, so it follows a forced rank, block weighting, sparsity and the latent adjustment unchanged.
- **Rejections before any fit.** A method other than `pls`; attribution (its reconstruction is defined on the fitted span); fewer than three stages (the subspace has one dimension, where `angle` is 0° or 180° by construction, the case the rank floor exists to prevent).
- **Analysis.** `TrajectoryAnalysisParams.latent_coordinates` accepts the three predictive values, is validated in `__post_init__`, omitted from `to_dict()` while unset, passed to the shared core and recorded under `measurement`. It is rejected with attribution and with `span_split_test`: the split's span part is the fitted span, and whether it should follow a K − 1-dimensional measurement is undecided. `orthonormal` is not offered there; it lost its verdict and was never an analysis option.
- **No CLI flag.** The plasmode profile needs the field; a case study does not need the flag before a verdict. `AGENTS.md` records the exception.
- **Alternative: keep `motco.analysis` untouched and let the driver measure from `latent_artifacts`.** The driver would have to rebuild the covariate-adjusted models and the residual correction from private helpers. A typed field on the shared core keeps one measurement path and its equivalence test.

**D5. One production fit, six arms per unit.**
- **The unit.** It generates once, fits `pls` once in scores, builds the three predictive integrations from that fit without refitting, and fits the two forced-rank arms without a CV (`forced_k_minus_1` at K − 1; `predictive_rank_plus_one` at A + 1, then the `predictive` step). It evaluates every arm with `latent=` at the unit's RRPP seed, serially, so the draws match.
- **Why the production arm reproduces exactly.** It is the parents' evaluation, unchanged.
- **Consistency checks**, which reproduction cannot provide because it reads only `pls_cv`: rank(Q) equals min(A, K − 1); `T·V·Vᵀ·Qᵀ = T·Qᵀ` to 1e-8; the three predictive arms report equal span-capture shares; and at A ≤ K − 1 the `predictive` arm's p-values equal `pls_cv`'s.
- **On the plasmode** the unit calls `analyze_trajectories` once per arm with the same parameters and seeds. `pls_cv` runs the CV. Each predictive arm runs at `forced_components` = the unit's `pls_cv` rank, because its final fit is the same `fit_plsda_model` call as the CV's, so the CV is not repeated. Its record carries `pls_cv`'s AUROC, and `same_fit` checks that its fitted model is `pls_cv`'s. `forced_k_minus_1` runs at `forced_components = 2`. The fidelity values confirm the label-free fit. `predictive_rank_plus_one` is not run there.
    - *Amended 2026-11-12, before any array unit* (`notes.md`, Smoke). As first written, every predictive arm repeated the label-free CV, which took about 100 s per arm and three quarters of the plasmode's cost. On the smoke's units the forced call reproduced each arm's statistics and p-values exactly. No rule, unit, seed or reproduction changed.

**D6. The profiles.**
- **`simulation` is a fresh family.** The gains were read on the production selection's datasets (base seed 2035) and on the probe's (4242), so neither can confirm them. Base seed 2036 on the comparison's generators and seed derivations gives new datasets with nothing else changed.
- **The planted scale.** Under `span_scale="planted"` a plant is exact where the sample scaler looks, so an off-target change is attributable to the measurement and not to the construction. The probe ran the comparison's null scale, which is one reason its off-target readings in the band cannot be attributed.
- **The band.** n = 200 at δ ∈ {1.0, 0.5, 0.25} are the design-guidance slice's stage-effect columns, where the held-out AUROC stays at 0.997–1.000 while the median CV rank rises to 3, 4 and 11. The SEA-AD table at δ = 0.5 adds the case-study cell sizes. There a forced rank of K − 1 and the predictive subspace are different spaces (largest principal angle near 88° at δ = 0.25 in the probe).
- **The effect ladder** is one ladder for every design: the smallest effects keep the candidates below saturation at δ = 2.0 (the probe's 22.5° rotation was already at 1.00), and the largest keep production above the null in the band.
- **The surgery cell.** The Phase 6 orientation surgery at e = 1.0 is largely outside the stage span; production reads it through components beyond the third (`angle` 0.996 above three components, 0.280 otherwise). It measures what the lower dimension gives up and is descriptive.
- **`weak_size` reuses the orthonormal check's weak units.** They are controls, so no gain was ever read on them, and the span-split check's committed counts make the production arm's reproduction exact at R = 500 and 999 permutations, the grade the probe's `shape` reading needs.
- **`plasmode`.** The anchor runs at replicates 0–499 because the production selection's 0–99 hold 11 of the plasmode's 21 `delta` rejections, which put production itself above that family's ceiling. The plant cells reuse replicates 0–99: no predictive arm was read there, and production's committed rows make the reproduction exact.

**D7. The decision.**
- **S1 and the δ = 2.0 half of S3 are bounds; S2, G and the band half of S3 are paired.** Size and specificity are properties of an arm, so a candidate is held to its own rate wherever production's is established (the production selection's D4). S2 and G compare arms on shared data and draws, where the information is in the discordant units.
- **Why S3 is paired in the band.** Production's off-target rates under stage-span plants below δ = 2.0 have never been measured, and the probe read production's own `shape` at 0.50–0.60 under a pure rotation there (null scale, 20 units). An absolute bound could fail every arm whatever it changes, the situation the production selection amended for. The candidate must not exceed production: a one-sided exact test on the units only the candidate rejects against those only production rejects, under the count-plus-ceiling rule over m = 72. Production's absolute rates are reported beside it.
- **Why G reads δ = 2.0 only.** The roadmap's question is whether the gains read there are recovered on fresh data. One Bonferroni-protected gain among 18 cells answers it; the probe's paired counts (0/10 to 0/33) would pass it at R = 60. The band is where the payoff is larger, but no committed run has a gain there to confirm, so its gains are recorded with the same test and do not decide.
- **Why the weak design is in S1 only.** No arm had power there in any run. Its target cells are reported.
- **Naming one arm.** Among qualifying arms the pilot names the one with the fewest off-target readings above the bound over both S3 families, because specificity is the only axis on which the probe saw the three differ. Ties go to the arm nearest production: `predictive`, then `predictive_orthonormal`, then `predictive_fitted`.
- **Order.** Size before power before specificity before gain. Every criterion is recorded for every arm whatever the verdict.
- **Alternative: an absolute S3 bound everywhere.** Rejected for the reason above. If the planted scale makes production's band rates hold the bound, the paired test is the weaker rule there, and the descriptive table shows it; the paper grade can then bound it.

**D8. What is descriptive.**
- **`forced_k_minus_1`** shows what lowering the fit rank does on the same draws: its rates, its paired differences against `pls_cv` and `predictive`, and its `pooled_stage` capture.
- **`predictive_rank_plus_one`** shows how near the CV rank sits to over-aiming. In the probe's forced-rank ladder on the SEA-AD table at δ = 2.0, rank 3 read the rotation's off-target `shape` at 0.03 and rank 4 at 0.85, and power fell to zero once the fit interpolated the stage labels.
- **Estimation.** Per cell and arm, the median observed statistic beside the population's through the same map. In the probe the `predictive` arm's observed angle read 16° against 41° at δ = 0.25 and agreed within 3° at δ ≥ 1.0. A qualifying arm is a detector first.
- **Fit and capture.** The selected rank (median, q10, q90), the held-out AUROC, the `pooled_stage` and `group_difference` capture, and the in-sample share of the stage labels the fit explains.

**D9. Execution follows the recent validation drivers.**
- One script with `run --profile` / `merge --profile` / `decide`; `settings.json` bound by hash with the parents' sha256; round-robin shards; resume by complete unit key; JSONL records holding seeds and statistics only; an sbatch array with a profile variable and a dry-run variable, BLAS pinned to one thread, run on `ing` from a clean worktree at a pushed revision.
- Committed: aggregates, the decision, provenance and the report under `results/predictive-subspace-<date>/`. Records stay on the cluster. No donor id or label vector is written.

## Risks / Trade-offs

- **[The planted fixed point fails in the band, where it has never run]** → Failed units are recorded and counted per cell; the smoke runs every design; a design whose failures exceed a few units is amended before any array unit, with the amendment in `notes.md`.
- **[The weak `shape` control is a real excess]** → S1's `weak_size` family reads it at R = 500. `size_excess` is a valid outcome, and it would hold for the real-noise anchor too.
- **[All three arms hold S3 only because production leaks as much in the band]** → Production's absolute off-target rates are in the report, the named arm is the one with the fewest readings above the bound, and the paper grade bounds it.
- **[The CV rank sits one component from the leak]** → `predictive_rank_plus_one` measures it. A knife edge is a reason to hold a qualifying arm until a rank rule for a fixed measured dimension is studied.
- **[A predictive arm's reported statistics are attenuated]** → The descriptive table measures it. Measurement of the difference stays with roadmap open item 3.
- **[The gain does not reach the case studies]** → Expected: the probe gained nothing at δ = 0.11, where the span holds 0.08–0.10 of the stage program. The plasmode's plant cells are the only real-noise gain reading, and they are recorded, not required.
- **[Reproduction fails from environment drift]** → The gate stops the reading. A mismatch is diagnosed, never waived.
- **[Compute overrun]** → A local smoke measures per-unit time for every profile before launch. The levers are the band's SEA-AD design and the plasmode anchor's replicate count. The projection goes to the user with a cap.

## Migration Plan

Nothing migrates. The setting is opt-in and absent in both entry points.
- **On `propose_paper_grade`:** a separate change runs the named arm at 500 × 999 on the stage-span paper grade's columns and the band, and on its verdict another flips the defaults with pins for every committed configuration, adds the CLI flag, defines attribution and the split test on the measured subspace, and re-issues SEA-AD.
- **On `keep_production`:** the three values stay opt-in diagnostics, and the measurement-visibility page records what the measured dimension costs and why it was kept.
- **Rollback:** remove the values and the analysis field. Nothing reads them while absent.

## Open Questions

- The core-hour cap for the arrays is set from the smoke timing, at the approval gate.
