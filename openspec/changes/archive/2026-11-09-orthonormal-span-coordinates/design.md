# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **The algebra, checked numerically.** sklearn's `PLSRegression` has the following properties, each confirmed to machine precision (2026-11-09, n = 80, p = 300, rank 4):
  - The weights `W = x_weights_` are orthonormal.
  - `PᵀW` is unit upper triangular, so `R = x_rotations_ = W·U`, with `U = (PᵀW)⁻¹` also unit upper triangular.
  - `qr(R)` therefore returns `W` up to column signs (1e-16).
  - `x_scores_ = transform(X) = X_in·W·U` (4e-15), where `X_in = (X − _x_mean)/_x_std`.
- **One production fit.** `fit_pls_latent_space` (`src/motco/simulations/evaluation.py:919`) is the core shared by the harness and `motco.analysis`. After the CV it fits one final model and returns `x_scores_` as the latent matrix, or `model.transform(X)` under latent adjustment. It hands the same model to attribution and to `LinearLatentMap.from_pls`, whose basis is already `qr(x_rotations_)`.
- **Who already uses the orthonormal basis.** `LinearLatentMap.basis` is used by:
  - the span-capture shares;
  - the leak check's `pop_span` checkpoint;
  - the split interaction test (`analysis._span_split_inputs`).

  `population_latent` and the leak check's `pop_lat` use `LinearLatentMap.transform`, which is `model.transform`, the scores.
- **The other linear spaces are already orthonormal.** `pca` measures `X·V`, `pls_within_group` uses a QR basis, and `concat` uses the identity. Only `pls`, and the first block of `pls_residual_pca`, are sheared.
- **Attribution's projection is oblique.** `stats.attribution._reconstruct_means` uses `inverse_transform(transform(m))`, which is `X_in·R·Pᵀ` rescaled. That is a projection onto the span of the loadings `P` along the null space of `Rᵀ`, not the orthogonal projection onto `span(W)`.
- **The parents.**
  - The stage-span paper grade (`results/stage-span-paper-grade-2026-11-04/`) committed per-cell `n_available` / `n_rejected` (p < α, as `spectrum.py` counts) for every column, cell and statistic, plus the leak check's per-unit path lengths.
  - The span-split check (`results/span-split-2026-11-05/simulation/production_table.csv`) committed `pls_cv` counts at 999 permutations for the comparison's anchor and translation at replicates 0–499.
- **Costs.** The paper grade ran 48,000 units in 229 core-hours, about 17 s per unit, dominated by data generation and the PLS CV fit. A `pls_cv` RRPP pass costs about 0.43 s at 199 permutations (latent-space comparison), so about 2 s at 999.

## Goals / Non-Goals

**Goals:**
- Measure the trajectory statistics on the fitted span with its own Euclidean metric, as one opt-in setting that changes nothing while absent.
- Validate it against the production scores on units where the score shear's leak was measured, paired on data, fit and permutation draws, under a predeclared decision.

**Non-Goals:**
- Changing the production default, `motco analyze`, attribution or study configurations (proposal, Non-goals).
- Changing the fitted span itself: the CV, rank rule, scaling and rank floor stay as they are.
- Re-reading committed results in the new coordinates beyond the two parents. That re-issue belongs to the default flip.

## Decisions

**D1. The coordinates are `X_in·Q` with `Q = qr(x_rotations_)`.**
- **One definition of the span's basis.** It is the basis `LinearLatentMap.from_pls` already exposes. Span capture, the leak check's `pop_span` and the split test therefore all measure the coordinates the statistics will be measured in.
- **The same span.** Since `R = W·U` with `U` unit upper triangular, `Q = W·diag(±1)`. The orthonormal coordinates are `X_in·W` up to column signs, and every statistic is invariant to the signs (D2).
- **The input coordinates are the model's.** `X_in` uses the model's own `_x_mean` and `_x_std`. The features arrive pooled-standardized, so `_x_std` differs from 1 only by the uniform factor `√(n/(n − 1))`. That factor scales `delta` and leaves `angle`, `shape` and every p-value unchanged. Under latent adjustment (`scale=False`), `X_in = X − _x_mean`, and the step is the same.
- **Alternatives.**
  - **`x_weights_` directly.** It is equal up to signs, but it relies on an sklearn implementation property that `qr(R)` does not need. It would also give the span two definitions in MOTCO.
  - **Whitened or decorrelated scores** (`T·(TᵀT)^(−1/2)`). This replaces one data-dependent metric with another and does not preserve the feature space's lengths.
  - **QR of the score matrix `T`.** This orthonormalizes sample columns, not feature directions, and rescales each axis by the sample's spread.

**D2. Basis invariance carries the statistics.**
- **What is invariant.** `delta` (path length), `angle` (between two direction vectors in the same coordinates) and Procrustes `shape` (reflections aligned away) are invariant to any orthogonal map `O` applied to every row.
- **Why the rest follows.** The LS-mean fit is linear, so `betas(Y·O) = betas(Y)·O`. RRPP's permutation indices come from `default_rng(seed)` regardless of the outcome's columns.
- **Consequence.** Every orthonormal basis of the span gives the same statistics and p-values, and the two arms of a unit see identical permutation draws. The spec's basis-invariance scenario tests this directly.

**D3. The setting is a top-level evaluation parameter.**
- **Shape.** `SimulationEvaluationParams.latent_coordinates`: `None` (scores), `"scores"` normalized to `None`, or `"orthonormal"`. It joins `grid._SIGNATURE_OMIT_WHEN_NONE`.
- **Why top-level.** This is the `orientation_estimator` pattern: a measurement choice whose later default flip needs the same mechanics (a hash-omitted spelling, explicit pins in committed configurations).
- **The core.** `fit_pls_latent_space` gains a keyword `coordinates` (default `"scores"`). On `"orthonormal"` it replaces the matrix after the final fit and returns a `LinearLatentMap` whose `transform` is `inputs(rows)·Q`, so `population_latent` follows. When set, the metadata records `integration_params.latent_coordinates`. `motco.analysis` passes nothing, so it is unchanged.
- **Alternative: `integration_params["coordinates"]`.** It would also be invisible while absent and axis-able. But integration parameters describe how the space is fitted, and the default flip would need its pins in a nested dictionary.

**D4. Attribution is rejected under orthonormal coordinates.**
- **The conflict.** Attribution's PLS-captured component is the oblique projection onto `span(P)` (Context). The orthonormal measurement corresponds to the orthogonal projector `Q·Qᵀ` onto `span(W)`, and the two differ unless `P` and `W` span the same space. Pairing the two would let the attributed orientation disagree with the tested one.
- **Why reject now.** Redefining the captured component as `Q·Qᵀ` changes attribution's outputs and needs its own validation. Rejecting before fitting keeps every path honest until then.
- **Forced rank.** A forced rank stays allowed, because the coordinate step is independent of how the rank was chosen.

**D5. One fit, two arms per unit.**
- **The unit.** It fits the `pls` integration once in scores and builds the orthonormal integration from it with `pls_orthonormal_from_pls` (no refit, equal to a from-scratch fit; spec). It evaluates both with `evaluate_semisynthetic_trajectory(..., latent=)` at the unit's RRPP seed, serially (`n_jobs = 1`), so the draws match.
- **Why the scores arm reproduces exactly.** It is the paper grade's evaluation and the span-split check's production arm, unchanged.
- **Consistency checks.** These catch an implementation error that reproduction cannot, because reproduction only reads the scores arm: `QᵀQ = I`, and `T = Z·(Qᵀ·R)` to 1e-8.
- **Cost.** The second arm adds one RRPP pass, about 10–15 % of a unit.
- **Alternative: a design-grid column in the study engine.** It would refit PLS per column (twice the compute) and would need the study-configuration key, which is a non-goal.

**D6. The profiles are the stage-span paper grade and the weak design's controls.**
- **Paper grade.** Its planted-scale plants are exact where the sample scaler looks, and its leak check measured the shear on these exact datasets. Any off-target change between arms is therefore attributable to the coordinates. Its anchors and translation give size, and its magnitude cells give `delta` power. The unit builder imports `stage_span_leak_check.column_cells` and `unit_params` unchanged, with `span_scale="planted"` and `span_capture=True`.
- **Weak design.** At SEA-AD's stage signal (δ = 0.11) the fitted span is mostly noise, the rank varies, and the shear is least constrained. Size there is the minimum a production change owes before `motco analyze` can follow. Its units come from `span_split_check.simulation_units` / `latent_space_comparison` unchanged. No power is read there; the comparison found none for any plant.
- **Excluded.** The surgery profile: its plants keep 0.2–0.8 % of their energy in the span, so the coordinates cannot matter to its power, and its anchors are the shared anchors.

**D7. The decision.**
- **O1 and O2** reuse `stage_span_paper_grade_decision.judge_rates` (imported unchanged) on the same families as P1 and P3, so the orthonormal arm is held to exactly the rule the scores passed. The weak family (m = 6) is judged separately, because its R and bound differ.
- **O3 is non-inferiority by an exact paired test.**
  - **Why paired.** The arms share data and draws, so the information is in the discordant units. Each cell is a one-sided binomial test: P(Binomial(b + c, ½) ≥ b), where b counts units only the scores arm rejects and c counts units only the orthonormal arm rejects.
  - **Why exact.** When the arms nearly agree, b + c is small and the normal SE degenerates, to 0 when b = c = 0.
  - **Multiplicity.** Over 45 cells, a per-cell 2-SE rule would false-fail with probability about 1 − 0.977⁴⁵ ≈ 0.65. The count-plus-ceiling form (at most the 95th percentile of Binomial(45, 0.025) cells below 0.025, none below 0.05/45) mirrors the paper grade's rule.
- **Adoption needs no rate gain.**
  - **Why.** The gain is structural and already measured at the population level. On an orthonormal basis the planted bend keeps its path to within 0.11 %, against −0.48 % / −0.59 % in scores, and the rotation's residual is +0.66 % against −3.35 %. The rate-level consequences sit within Monte Carlo noise at R = 500: the bend's `delta` reads 0.016–0.036, and the rotation's `shape` 0.070 / 0.086 sits about one SE above α.
  - **The cost of requiring a gain.** Adoption would hinge on one or two cells with little power to separate.
  - **What remains visible.** The watched cells are reported with their paired test (spec), so a reader sees whether the shear explained them.
  - **Alternative: a predeclared gain criterion on the watched cells.** It was rejected because it is underpowered.
- **Order.** Size before specificity before power: a size failure makes the others moot for production. All criteria are always recorded.

**D8. Execution follows the recent validation drivers.**
- **Driver shape.** One script with `run --profile` / `merge` / `decide`. Shards write JSONL records, and resume skips complete unit keys. `settings.json` is bound by hash and records the parents' sha256.
- **Cluster runs.** An sbatch array (`ORTHO_PROFILE`), run on `ing` from a clean worktree at a pushed revision (`PYTHONPATH` on its `src`), with BLAS pinned to one thread.
- **What gets committed.** Records stay on the cluster. Aggregates, the decision, provenance (revision, tree hash, job ids, `sacct`) and the report are committed under `results/orthonormal-coordinates-<date>/`.

## Risks / Trade-offs

- **[Reproduction fails from environment drift since 2026-11-04]** → The gate stops the reading. Run on the same `ing` environment and record the numpy, scipy and sklearn versions. A mismatch is diagnosed before any decision, never waived.
- **[The 12 four-stage orientation units fail again (planted fixed point)]** → They are recorded as failures in both arms. Reproduction compares `n_available` too, and pairing uses units available in both arms.
- **[The orthonormal arm loses power, because the shear stretched the stage directions]** → O3 catches it, and `power_loss` is a valid outcome. The report then documents the trade-off, and production stays on scores.
- **[The arms are near-identical, so adoption buys a re-issue for an invisible gain]** → Accepted by D7. The structural argument and the watched cells are reported, and the flip remains a separate change the user approves.
- **[Eigengap and `angle` null-width readings are coordinate-dependent]** → Committed eigengap-stratified tables would move under a flip. The descriptive tables quantify the change per cell, so the flip's re-issue can be scoped.
- **[Compute overrun]** → A local smoke measures per-unit time for both profiles before launch. The projected core-hours go to the user for approval with a cap.

## Migration Plan

Nothing migrates. The setting is opt-in and absent everywhere.
- **On `adopt_orthonormal`:** a follow-up change flips the harness default, pins every committed configuration and diagnostic to `"scores"` (the `default-net-displacement-orientation` pattern), and opens the `motco analyze` adoption with its own real-noise check.
- **On any other verdict:** the setting stays as an opt-in diagnostic, and the measurement-visibility page records the shear as a measured limitation.
- **Rollback:** remove the field. Nothing reads it while absent.

## Open Questions

- The core-hour cap for the array is set from the smoke timing, at the approval gate.
