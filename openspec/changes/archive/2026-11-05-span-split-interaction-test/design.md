# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **The trajectory models.** `build_simulation_trajectory_design` gives the full model (group × stage) and the reduced model (group + stage). `analyze_trajectories` appends covariates to both. RRPP permutes the reduced model's residuals (Freedman–Lane), so the tested effect is the group × stage interaction, in every path.
- **RRPP's draws.** The serial path draws `idx = rng.permutation(n)` per permutation from `np.random.default_rng(seed)`. The parallel path spawns per-worker seeds, so its draws differ. `residual_correction_scale` gives the `df` and `leverage` multipliers.
- **The fitted span is already exposed.** `LinearLatentMap.from_pls(model)` gives the PLS input transform (`_x_mean`, `_x_std`) and `qr(x_rotations_)`. `analyze_trajectories` exposes the fitted model and its standardized features in `latent_artifacts`. The comparison driver holds the same through its `pls_cv` integration.
- **Two committed parent runs.**
  - The latent-space comparison (`results/latent-space-comparison-2026-11-04/`: settings, calibration, `rate_table.csv`, `capture_table.csv`) ran 100 × 199 at base seed 2035.
  - The plasmode (`results/plasmode-2026-11-02/`: settings with the export's sha256, `plasmode_summary.csv`) ran 500 × 999 at base seed 2034 on the SEA-AD export, which is available locally and on `ing`.
- **Costs.** A comparison unit spent about 20 s on data generation and the PLS CV fit, and 0.43 s on `pls_cv`'s RRPP at 199 permutations. `concat`'s RRPP took 16 s, because every draw re-estimates the trajectories in 498 dimensions.

## Goals / Non-Goals

**Goals:**
- A test of the group × stage interaction, with an exact decomposition into the fitted span and its complement, valid without refitting any basis.
- A permutation cost independent of the number of features, so 10,000 features cost what 498 do.
- Validation on units that reproduce committed runs, so its results pair with the five-arm comparison and the plasmode.

**Non-Goals:**
- Geometry of an out-of-span difference, and a localization finer than span against complement.
- Changing the production statistics, their coordinates or their defaults.
- A multiplicity policy across the three parts for real-data reading. That belongs to the case study's reading rules (Open Questions).

## Decisions

**D1. The statistic is the trace pseudo-F of the interaction.**
- **Definition.** Per part: SS_int = ‖(H_F − H_R)·Y_part‖², RSS = ‖(I − H_F)·Y_part‖² and F = (SS_int/df_int) / (RSS/df_res). This is the classical multivariate interaction sum of squares; with two groups df_int = K − 1.
- **Why a ratio.** The Freedman–Lane residuals have mean square RSS_R/n, not RSS_R/(n − r_R), so the permuted null's interaction energy is too small by about n/(n − r_R). Item 8 measured that shrinkage at the case-study configuration. A sum over p coordinates concentrates: its relative null SD is about √(2/(df_int·p_eff)). So a few-percent scale bias puts the observed energy several null SDs out, and the null becomes liberal in proportion to dimension. The pseudo-F computes both sums from the same permuted residuals, so the scale cancels, and it is the standard pivotal choice for Freedman–Lane.
- **Relation to span capture.** With balanced cells, SS_int = (n_cell/2)·Σ_s‖d_s − d̄‖². That is the same energy that span capture's `group_difference` measures on population means, up to a constant.
- **Alternatives.**
  - The energy alone with the `df` correction restores the null's mean, but only under homogeneous leverage. It is kept as a diagnostic (D10).
  - The unweighted LS-mean energy Σ_s‖d_s − d̄‖² matches span capture's definition. But it is not a sum of squares of the fitted model and needs the same correction.
  - Hotelling- or Wilks-type statistics need p < n. They are undefined at 10,000 features.
  - A per-feature max-T scan targets sparse differences. It is thread 3.

**D2. The split is an orthogonal projection in the PLS input coordinates.**
- **The basis** is `qr(x_rotations_)` in the model's input coordinates `(x − _x_mean)/_x_std`, the coordinates in which it is orthonormal and in which span capture is measured. For pooled-standardized features `_x_std` is a uniform factor, so `total` is the standardized feature space up to a scale, which the pseudo-F ignores.
- **Exact decomposition.** The parts use Y·Q·Qᵀ and Y·(I − Q·Qᵀ). Orthogonal projectors split both sums of squares exactly, and the observed in-span share is SS_span/SS_total.
- **Why not the score coordinates.** The scores `X·W·(PᵀW)⁻¹` are a sheared image of the span (item 15). Splitting in them would not decompose the energy.
- **Why the fitted span.** Real data has no population span. The reading is therefore defined against the fitted one, and the documentation says so (spec `latent-space-metric-compatibility`).

**D3. Validity rests on a fixed basis built without the group labels.**
- **The rule.** The basis comes from the features and the stage labels, never the groups, and stays fixed across permutations, exactly like the production statistics' basis.
- **Why it can hold exactly.** The cell-size-weighted pooled stage mean is uncorrelated with the within-stage group difference for any cell sizes: cov(w_A·m_A + w_B·m_B, m_B − m_A) = σ²/(n_A + n_B) − σ²/(n_A + n_B) = 0. A rank-(K − 1) span built from the pooled means is therefore independent of the interaction's noise under Gaussian noise.
- **Where it is empirical.** Components beyond K − 1 follow the total covariance, which contains the within-stage differences. So validity there is an empirical claim, checked by C1 at every design and by C2 on real noise.
- **Alternative: refitting the basis per permutation.** That is the price of a group-aware basis (thread 2), and it is unnecessary here.

**D4. Permutations run on Gram matrices.**
- **The reduction.** H_int·Ŷ_R = 0 and (I − H_F)·Ŷ_R = 0. So a permuted outcome's interaction and residual sums depend only on the permuted reduced residuals E_part[π].
- **The computation.** Each part's G = E_part·E_partᵀ (n × n) is computed once. Per draw, SS* = ⟨H_int, G[π, π]⟩ and RSS* = ⟨I − H_F, G[π, π]⟩. The cost is O(n²) per draw and part, whatever p is.
- **Numerics.** E_span = E·Q and E_comp = E − (E·Q)·Qᵀ are formed directly, which costs n·p·r once. Forming G_comp as G_total − G_span would cancel badly when the span holds most of the energy.
- **Corrections** scale G by s·sᵀ, with s from `residual_correction_scale`. Under `df`, s is constant, so F* is unchanged and SS* scales by κ².
- **Draws** use `np.random.default_rng(seed)` and one `rng.permutation(n)` per draw, RRPP's serial sequence. The test can therefore be paired with a trajectory analysis at the same seed.
- **Alternative: materializing Ŷ_R + E[π] and refitting per draw.** That is what `concat`'s RRPP does, and it costs 16 s per unit at 199 draws in 498 dimensions.

**D5. Where it lives.**
- **`motco.stats.interaction`** holds the test and imports nothing from `motco.simulations`. Callers map their features to input coordinates and pass the basis.
- **`analyze_trajectories`** calls it after the fit, with `latent_artifacts`' features and model and the covariate-adjusted models.
- **The validation driver** calls it with its `pls_cv` integration's linear map.
- **The simulation harness and study configs** are left out (non-goal). The driver and a cross-path test cover the semi-synthetic side.

**D6. The analysis record.**
- **`span_split` follows `measurement`.** It holds the parts (dimension, SS_int, RSS, df_int, df_res, F and the p-values), `span_share`, `residual_correction`, `permutations`, `seed` and `draws_match_rrpp` (true only with one job).
- **The resolved parameters** carry `span_split_test` only when it is true. That keeps absent-setting records byte-identical, as for `residual_correction`.
- **Latent adjustment is rejected before fitting.** Its basis is fitted on residualized features while the scores project the raw ones, and a split in that pair of coordinates has not been defined or validated.

**D7. The simulation profile: the comparison's units at mixed replicate counts.**
- **Units.** The same three designs, fourteen cells, weak δ (read from the committed calibration, sha256 checked) and seed derivations. Replicates 0–99 are therefore the comparison's datasets and RRPP seeds.
- **Replicate counts.**
  - The anchor and translation run 500 replicates. Size at every design is then measured at paper grade (bound 0.0695), including the weak design where `concat`'s `angle` was liberal.
  - The other twelve cells keep the comparison's 100, so their rates pair with its five arms.
- **Cost.** About 6,600 units at about 23 s each, so about 42 core-hours. Running all cells at 500 replicates would cost about 134 core-hours. The extra precision would serve only C3, whose expected margin is large.
- **Permutations.** 999 everywhere. RRPP's serial draws are nested: the first 199 of 999 are the 199-draw sequence. The production arm's p-values from the first 199 draws therefore reproduce the comparison exactly.

**D8. The plasmode profile: real noise at the case-study configuration.**
- **Units.** The plasmode's units: the export contract, within-stage label permutation, the 20-feature plant and its seeds. The cells are the anchor, e = 2.0 and e = 3.0 at 500 × 999. Each unit is one `motco analyze` with the six covariates, floor 2, net displacement, `df`, the split test on and attribution off.
- **Why.** It is the repository's only real-data null: real noise, covariates and cell table, at 10,000 features. It is where a high-dimensional energy test is most likely to misbehave. Its plant is exactly the difference the production measurement missed.
- **Reproduction.** Attribution draws from its own bootstrap seed, so turning it off leaves the RRPP draws unchanged, and the trajectory rejection rates must equal the committed summary.
- **Cost.** About 1,500 units without attribution, estimated at no more than 25 core-hours. The smoke replaces the estimate.

**D9. Reproduction gates the decision.**
- **Why exact reproduction.** The descriptive comparisons with the five arms and the plasmode are only meaningful on the same datasets, fits and draws.
- **What is checked.** The committed aggregates suffice: exact counts and rejections per design, cell and statistic for the comparison, and exact rates per cell and statistic for the plasmode.
- **On failure** the verdict is `reproduction_failed`, as in the case-study diagnosis, and nothing else is read.

**D10. The criteria.**
- **C1 and C2: size.** They reuse the stage-span paper grade's count-plus-ceiling implementation (`scripts/stage_span_paper_grade_decision.py`) by import.
  - m = 18 at R = 500 for C1. That allows one or two chance exceedances of the per-cell bound and none above the Bonferroni ceiling.
  - m = 3 for C2. That allows none.
- **C3: gain.** It compares the complement with `pls_cv`'s `angle`, the production target statistic for the rotation, on the same units, at both δ = 2.0 designs, by more than 2 pooled SE.
  - It does not compare with `concat`, which is a reference using a different statistic and is liberal at the weak design.
  - The weak design is not used, because stage-span plants there are as weak as the program they reshape.
- **Localization is descriptive, not a criterion.** A complement rejection of an in-span plant still means energy outside the fitted span, from the 3–7 % capture shortfall or estimation noise, which is the documented reading.
- **Precedence:** `reproduction_failed`, then `size_excess`, then `no_out_of_span_gain`, then `split_validated`.
- **The SS diagnostic.** The interaction-SS p-values, uncorrected and under `df`, are recorded by part and design. The prediction is that uncorrected energy is liberal and more so in the complement (high dimension), while `df` and F hold. That tests D1's rationale. It may bear on `concat`'s liberal `angle` at the weak design but does not test it, and the report must not claim it does.

**D11. Documentation.**
- **The measurement-visibility page** gains a section on testing outside the span: the parts, the fitted-span caveat, the reading of rejections and non-rejections, and the ratio. It is filled from the report once the decision exists.
- **API pages:** `docs/api/sd.md` documents the primitive and `docs/api/analysis.md` the option and the record block.
- **Roadmap and repository guide:** roadmap item 16 records the verdict, and `AGENTS.md` gains the module, the option and the driver.

## Risks / Trade-offs

- **[An energy test in 10,000 dimensions dilutes a 20-feature plant]** → The plasmode profile measures that power, and a weak result is reported as such. Thread 3 (sparse scan) is the targeted follow-up.
- **[Fitted components beyond K − 1 depend on the within-stage differences]** → C1 at every design and C2 on real noise check size. `seaad_weak` (unequal cells, a span fitted mostly to noise) and the real export (rank 5) are the hard cases.
- **[R = 100 in the power cells]** → They are descriptive except C3, whose expected margin (`concat` 1.00 / 0.90 against `pls_cv` 0.08 / 0.11) far exceeds 2 pooled SE (about 0.09).
- **[Three p-values per analysis invite selective reading]** → The documentation states which question each part answers. The case study's reading rules fix the multiplicity policy before any real-label analysis.
- **[Parent inputs could drift]** → The calibration file's and the export's sha256 are checked against the committed settings before any unit.
- **[More than one job changes RRPP's draws but not the split test's]** → `draws_match_rrpp` records it. The validation runs with one job per unit.
- **[The weak design shows no power anywhere]** → Its plants are relative to a weak program, as in the comparison. The fitted-span caveat is shown there through capture and the observed share, not through power.
- **[Compute estimate]** → The local smoke times both profiles, and the projection in `notes.md` replaces the estimate before launch.

## Migration Plan

Opt-in throughout, so there is nothing to migrate.
- **Version:** MOTCO 0.18.0.
- **Rollback:** without the analysis setting every record and table equals 0.17.0's, and the stats module is new.

## Open Questions

- How the case study adjusts for reading three parts: a fixed primary part, Holm, or no adjustment with the omnibus read first. It belongs to `sead_sexdiff`'s change, before any real-label analysis, and changes nothing here.
