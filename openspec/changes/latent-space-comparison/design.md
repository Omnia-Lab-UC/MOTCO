# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **The stage-span construction (0.15.0).** It edits group B's null population stage means in null-pooled-standardized joint coordinates over `span_layers`. It draws no randomness, so a dataset is the null dataset plus the plant. Its orientation rotates B's configuration by θ = 90°·e in the plane of u, the net displacement, and v, the leading direction of the centred configuration once u is removed. Both directions lie inside the null stage span.
- **Span capture (0.15.0) is PLS-only.** `measure_span_capture` takes the fitted PLS model, maps population means to the model's input coordinates and measures shares against `qr(x_rotations_)`. `population_latent` uses `model.transform`.
- **The diagnostic `pca` exists.** It is a label-free SVD of the pooled-standardized `concat` matrix at a required `n_components`. The orientation probe measured it at rank 2(K − 1), together with `pls_within_group`, whose fixed basis broke size (anchor `angle` 0.48–0.88).
- **The probe driver is the template for multi-arm comparisons.** `scripts/orientation_measurement_probe.py` draws one dataset per unit, integrates each arm and evaluates it with `evaluate_semisynthetic_trajectory(..., latent=)` at one shared RRPP seed. It writes resumable shards bound to a settings hash.
- **Two orders of magnitude separate the stage signals.** The simulator's design-guidance generator (per-omic δ = 2.0) gives a held-out stage AUROC of 1.000. SEA-AD astrocytes `cohort_30` read 0.636. The rank-floor check read 0.52 and 0.58 at δ = 0.05 and 0.1 on the SEA-AD table.
- **Costs.** The stage-span pilot ran 16.9 s per single-arm unit at 199 permutations, including one CV fit. The probe ran ten arms on 1,000 units in 28 core-hours.

## Goals / Non-Goals

**Goals:**
- One construction family in which only the direction of the difference relative to the stage span varies, with the intended geometry held fixed.
- Span capture defined for every linear space, so that capture and power are read on the same footing across arms.
- A pilot that reads every arm's visibility, size and protected cost on identical data, with a predeclared advisory reading.

**Non-Goals:**
- Choosing a production space. The reading proposes a paper grade or documents a trade-off; nothing changes in `pls`, its defaults or `motco analyze`.
- Absolute-scale plants. Stage-span effects are relative to the configuration, so at a weak stage signal the plant is small; the plasmode is the absolute-scale check.
- Making `pls_residual_pca` or `pca` available to `motco analyze`.

## Decisions

**D1. Alignment tilts the rotation plane, holding the intended geometry fixed.**
- **The rotation.** w_φ = cos φ·v + sin φ·o, with o a unit direction outside the null stage span (D3), and B's configuration rotated by θ about its centroid in the (u, w_φ) plane. u, v and o are mutually orthogonal, so w_φ is a unit vector orthogonal to u.
- **Fixed geometry.** The net displacement lies along u, so the intended `angle` is θ at every φ. The rotation is rigid, so the intended `delta` and `shape` are 0. Only where the difference points changes with φ, so differences in capture and power across φ are visibility, not effect size.
- **At φ = 90° the difference is outside the pooled span.** With A's centred stage means a_s, B's difference at stage s minus its mean over stages is ((a_s − ā)·u)·e_d, where e_d = (cos θ − 1)u + sin θ·o. Against the centred pooled means p_s = (a_s − ā) + ((a_s − ā)·u)·e_d / 2:
  - e_d·p_s = ((a_s − ā)·u)·[(cos θ − 1) + ‖e_d‖² / 2];
  - ‖e_d‖² = 2(1 − cos θ), so the bracket is 0.
  
  This is the rhombus, exact for every θ. It is the controlled analogue of the relocation surgery: the same rigid θ, with no cascade or relocation noise.
- **Where it holds.** It needs balanced cells and rank K − 1. The SEA-AD table's unequal cells (the kite) and a CV rank above K − 1 leave it partly visible; capture measures by how much.
- **φ = 0 normalizes to absent,** so the edit, truth, cell id and signature equal 0.15.0's.
- **Alternative: an orientation onto a relocated copy of the program** (u permuted within blocks). It ties the axis to the surgery's support and does not hold θ exactly once projected; rejected.

**D2. Alignment applies to orientation only.**
- **Shape.** At three stages every configuration is a triangle, which is planar. Bending the middle stage out of the plane gives a different triangle than bending it within the plane by the same displacement (its distance to the chord becomes √(h² + Δ²) instead of h + Δ). So the intended `shape` would change with φ and confound visibility with effect. Shape stays at φ = 0 as the protected check of each space's shape test. At four or more stages an out-of-span bend is distinct; that belongs to a later change with four-stage designs.
- **Magnitude and translation.** Scaling keeps B's steps parallel to A's, which is in-span by definition. Translation is invisible to all three statistics. Neither has an out-of-span version.

**D3. o lives on features the null program does not use, from its own stream.**
- **Support.** o's support is the span-layer features whose standardized null stage means are constant across stages (|z_s,f − z̄_f| ≤ 1e-10 for all s). The null span is supported elsewhere, so o is exactly orthogonal to it without a projection, and a sparse o stays exactly sparse. Entries are standard normal on the support, then normalized. `span_support` = k picks k support features uniformly first.
- **Norm.** The plant's norm is set by θ and the configuration, not by o, so dense and sparse cells plant the same energy. That is what the sparse-versus-dense contrast needs, and it tests the visibility page's claim that blindness is a direction, not a feature set.
- **Stream.** The draw uses `np.random.default_rng([seed, <alignment stream tag>])`. It never touches the sampling stream, so the dataset stays the null dataset plus the plant. The unused set depends only on the baseline indicators, which are fixed at a seed, so o is the same for every effect and φ at a replicate and support. The φ curve is therefore a path on identical noise.
- **Truth** records φ, the support size, the span layers and a short hash of o's support and stream seed as the direction identifier.
- **Alternatives.**
  - A Gaussian direction projected off the span: dense everywhere and never exactly sparse.
  - A relocated program: D1's alternative.

**D4. `pls_residual_pca` appends within-stage components to the production scores.**
- **Fit.**
  1. Fit `pls` as in production (CV rank, floor 2).
  2. In the model's input coordinates X_in = (X − x_mean) / x_std, subtract the pooled stage means.
  3. Project the result off span(R), with R the rotations.
  4. Take the top r right singular vectors W.
- **The matrix.** The latent matrix is [`model.transform(X)`, X_in·W].
- **Why within-stage residuals.** Between-stage directions are already in R. Within each stage, the two groups' means differ by B − A, which adds variance along that direction to the within-stage covariance. So the residual components can find an out-of-span difference without the group labels.
- **The metric.** The first block is exactly production's scores, so the arm differs from `pls_cv` only by the appended coordinates, and any change in power or size is theirs.
- **r.** The default r = K − 1 gives a total near the probe's 2(K − 1) at the floor.
- **A shared fit.** A helper builds the arm from a prebuilt `pls` integration, so the driver fits the PLS model once per unit (spec: prebuilt fit).
- **Alternatives.**
  - Appending pooled PCs without stage residualization: they mostly re-find the stage directions.
  - A group × stage PLS: it uses the group labels, so it would need a refit inside every permutation, the `pls_within_group` lesson.

**D5. Validity rests on group-blind fixed bases, checked by size.**
- **The rule.** RRPP permutes the reduced model's residuals (group and stage main effects) and keeps each arm's basis fixed. A basis fitted from the features and stage labels depends on the data through pooled or within-stage covariances and stage means, which the permutations nearly preserve. A basis fitted from the group × stage means depends on the test's own target.
- **The candidates.** `pca`, `pls_residual_pca` and `pls` are all of the first kind. That is necessary, not sufficient, so the reading's size criterion checks every arm at every design (C1), and an arm that fails reads `size_excess`. No per-permutation refit is attempted.

**D6. Each linear integration carries its fitted linear map.**
- **The artifact.** `LatentIntegrationResult.artifacts` gains, for `pls`, `pca`, `pls_residual_pca` and `concat`, a linear-map record:
  - the input transform (centre and scale vectors);
  - an orthonormal basis of the kept directions in input coordinates;
  - the function that maps standardized joint rows to latent coordinates.
- **Span capture and `population_latent`** read that record instead of a PLS model.
- **The PLS record** uses `model._x_mean`, `model._x_std`, `qr(x_rotations_)` and `model.transform`, exactly as 0.15.0 does. PLS records stay byte-identical, and a test pins them against `tests/test_span_capture.py`'s fixtures.
- **`concat`'s** basis is the identity, so its shares are 1 by construction. It is recorded so the tables are uniform.
- **`pls_within_group` and `snf`** record none.

**D7. Arms and ranks.**
- **`pls_cv`** is production: CV rank with floor 2, net displacement and pooled standardization.
- **`pca_matched`** uses the unit's `pls_cv` rank, so that it isolates supervision from dimension. It is the primary label-free arm.
- **`pca_fixed`** at 2(K − 1) = 4 repeats the probe's arm for continuity.
- **PCA is the naive one requested:** unweighted blocks, so methylation (367 features) outweighs expression (131). Block weighting is roadmap item 5.
- **No data-driven rank rule.** A Marchenko–Pastur edge assumes white noise. Parallel analysis counts correlated noise as signal. Any variance threshold is arbitrary. A rule is needed only if a label-free space is proposed for real data, which is the paper grade's question.
- **`concat`** is the reference ceiling and is never eligible. The probe left its p-padding question open.

**D8. A driver, with its own seed family.**
- **Why a driver.** `scripts/latent_space_comparison.py` follows the probe's design: one dataset per unit, one PLS fit shared by `pls_cv`, `pca_matched`'s rank and `pls_residual_pca`, and one RRPP seed for all five arms. A study profile would make each arm a design column: three CV fits per unit, and anchors repeated across sibling profiles.
- **Seeds.**
  - The generator seed is a hash of (base seed 2035, design, replicate), so every cell at a replicate shares its null dataset (spec).
  - The RRPP seed is a hash of (base seed, design, cell, replicate).
  - Calibration draws from a separate `calibration` stream.
- **Why a fresh family, not the design-guidance one.** It needs no coupling to the study machinery's seed derivation. It also gives an independent replicate of the stage-span pilot's in-span cells, including the n = 200 anchor whose `angle` read 0.10 there.
- **The link to 0.15.0** is carried by two tests instead of a records reproduction: the φ-absent byte identity, and `pls_cv` equalling a direct production evaluation.
- **Subcommands and outputs.** `calibrate` | `run` | `merge` | `decide`. Shards are sharded round-robin and resumable, bound to the settings hash. sbatch template: `scripts/latent_space_comparison_array.sbatch`.

**D9. Calibrating the weak design's stage signal.**
- **The ladder.** All three per-omic δ are set together, as in the design-guidance δ arm: {0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5}. Each step runs 40 anchor replicates of the SEA-AD table, `pls_cv` fit only.
- **Adoption.** The first step whose median held-out AUROC lies in [0.60, 0.67] is adopted. Otherwise the closest to 0.636 is adopted, marked `band_missed`.
- **The prior.** The rank-floor check read 0.58 at δ = 0.1 on this table, so a step near 0.15 is expected.
- **Amendment (2026-11-04, after calibration job 933435 and before any run unit).** The original ladder put the band
  between its first two steps (δ = 0.1: median AUROC 0.595; δ = 0.15: 0.710), so the rule would have adopted δ = 0.1
  as `band_missed`. Steps 0.11, 0.12, 0.13 and 0.14 were added and calibration was run again on the amended ladder
  {0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5}, under the same adoption rule. The calibration seeds do
  not depend on δ, so the original steps reproduce. The first run is kept beside the adopted one.
- **What the weak design measures.** At the adopted δ, stage-span plants are as weak as the program they reshape. So this design is descriptive: it shows how each space behaves when the stage axis is barely there, which is SEA-AD's regime. It does not enter C2 or C3.

**D10. Designs and cells.**
- **Designs.** Three stages, methylation + expression, `p_dmp` 0.1, ρ = 0, `joint` magnitude, span layers = layers:
  - `n200`: n = 200 balanced, δ = 2.0;
  - `seaad`: the SEA-AD table `[[11, 10, 28], [9, 10, 12]]`, δ = 2.0;
  - `seaad_weak`: the same table at the calibrated δ.
- **Cells,** all under `construction = stage_span`:

  | cell | mode | e | φ | support |
  |---|---|---|---|---|
  | `anchor` | none | 0 | — | — |
  | `translation` | translation | 1.0 | — | — |
  | `magnitude` | magnitude | 0.05 | — | — |
  | `orient_e{0.25,0.5}_phi{0,30,60,90}` (8 cells) | orientation | 0.25, 0.5 | 0, 30, 60, 90 | dense |
  | `shape_e{0.25,0.5}` (2 cells) | shape | 0.25, 0.5 | — | — |
  | `orient_e0.5_phi90_k20` | orientation | 0.5 | 90 | 20 |

- **Why these effects.** They stay at or below e = 0.5, the cap the stage-span pilot proposed: rotations of at most 45° and bends of at most half a step. Magnitude at 0.05 sits where both constructions had intermediate power (0.72 at n = 200, about 0.2 at the SEA-AD table), so a loss is resolvable. The 20-feature support matches the plasmode's 10 + 10 plant.
- **Total:** 3 × 14 × 100 = 4,200 units of five arms.

**D11. The reading.**
- **Eligibility and bound.** The criteria are in the spec. Eligible arms are `pca_matched`, `pca_fixed` and `pls_residual_pca`. The bound is 0.0936 at R = 100.
- **Where each criterion applies.** C1 (size) applies at all three designs: an invalid null disqualifies wherever it appears. C2 (gain) and C3 (protected) apply at the two δ = 2.0 designs only, where `pls_cv` has power to compare against.
- **The two thresholds.** The protected criterion uses −2 pooled SE, the same margin the probe used. The gain threshold of 2 pooled SE keeps a chance-level gain from qualifying an arm.
- **Advisory.** The paper grade's predeclared decision owns any adoption. The reading only proposes, as in the stage-span pilot.
- **Descriptive output.**
  - `pls_cv`'s own control rates;
  - the sparse-minus-dense contrast;
  - the full φ curves.

**D12. The visibility page's corrections ride with this change.**
- **The rhombus paragraph.** It currently says the cancellation happens "stage by stage". It is exact only when B's steps keep A's lengths and mutual angles and are orthogonal to every one of A's steps, which is D1's φ = 90° case. A rotation within A's own features is partly or fully visible. A three-case numerical check on a 50-feature, four-stage configuration gave in-span shares of 0.000 (relocation), 1.000 (in-span rotation) and 0.010 (random orthogonal rotation).
- **The `concat` sentence.** It gains: a projection cannot see a difference `concat` lacks, but it can raise power by discarding noise directions (probe Phase 6 `shape` 0.19 for `pls_cv` against 0.00 for `concat`).

Both are documentation of the capability this change extends (`latent-space-metric-compatibility`).

## Risks / Trade-offs

- **[R = 100 cannot resolve rates near the bound (SE 0.022 at α)]** → The reading only proposes; a paper grade at 500 × 999 decides. Size failures are reported with their counts.
- **[`pca` and `pls_residual_pca` hold a fixed basis under RRPP]** → C1 at every design (D5). The fixed-basis failure mode is known from `pls_within_group`, and the anchor and translation controls caught it there.
- **[The weak design shows no power anywhere]** → It is descriptive by design (D9), and the report says so. It still answers whether any space rescues a SEA-AD-like stage signal on a relative plant.
- **[Naive PCA follows the larger block]** → Kept naive by request. The report states the block feature counts, and block weighting stays roadmap item 5.
- **[`concat` RRPP in 498 dimensions is the slowest arm]** → It is timed in the smoke, and the core-hour projection is replaced before launch. Dropping it from the weak design is the fallback.
- **[A sparse o puts large per-feature shifts on 20 features, including methylation M-values]** → Truth records the plant norm. `population_latent` reports the realized geometry, and the report compares it with the intended `angle` for the sparse cell.
- **[At the SEA-AD table the kite leaves a φ = 90° rotation partly visible to `pls_cv`]** → That is expected and measured. `n200` is the clean case, and the reading requires the gain at both δ = 2.0 designs.
- **[Few unused features at a high `p_dmp`]** → At `p_dmp` 0.1 most features are unused. The generator raises a named error when k exceeds the support, and the dense case needs at least one unused feature.

## Migration Plan

Opt-in throughout, so there is nothing to migrate.
- **Version:** MOTCO 0.16.0.
- **Rollback:** absent the new generator settings and integration method, every dataset, record and report is byte-identical.

## Open Questions

- Whether the paper-grade follow-up adds a four-stage design, where an out-of-span bend is defined (D2). The answer changes that change's grid, not this one's.
