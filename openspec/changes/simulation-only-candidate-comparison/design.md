# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **What each candidate is.**
  - `predictive_orthonormal` (0.24.0) measures on `X_in·qr(R·V)`, the K − 1 stage-predictive directions of the production fit. `pls_coordinates_from_pls` builds it from a score-coordinate fit without refitting, so it costs one RRPP pass.
  - `weighted` (0.21.0) fits stage-PLS on blocks scaled by 1/√p_b with `scale=False`. It is its own cross-validated fit, and its `delta` is in the weighted metric.
  - `sparse_block` (0.22.0) cross-validates rank × keep fraction with a budget per block. The production selection measured it at seven to nine production fits.
- **What has been read in simulation.**
  - `predictive_orthonormal`: seven two-block designs at 100 × 199 (item 18, base seed 2036); the same designs under within-stage draws (item 20); ADNI's simulated design at 500 / 200 × 999 (item 22).
  - `weighted` and `sparse_block`: three two-block designs at 100 × 199 (item 17, base seed 2035).
  - None has been read on the three-block datasets of base seed 2043, which production, the per-feature arms and the group-dispersion nulls have read.
- **The committed parents.**
  - `results/univariate-baseline-2026-11-23/stage_span/rate_table.csv` holds `delta`, `angle` and `shape` singly at p < 0.05 for `n80`, `n200`, `n500` and `n200_delta0.5`, thirteen cells each, controls at 500 and plants at 200. Its `surgery/` table holds the slice's three-block columns, eleven cells at 500.
  - `results/stage-effect-scaling-2026-11-24/` holds the same arms for `n320_delta1` (thirteen cells) and `band_n320` (eleven cells, δ_c = 0.05 from its committed calibration).
  - `results/adni-gain-specificity-2026-11-17/simulation/rate_table.csv` holds `pls_cv` and `predictive_orthonormal` on ten cells of the ADNI-design construction.
- **The generator after 0.26.0.** `generate_omics(noise_scale=)` stretches one call's samples about their exact cell means; the generator calls it once per group. `null_pooled_scale` and `planted_pooled_scale` take a `within_multiplier`. The group-dispersion change listed dispersion by stage as a non-goal.
- **What item 20 measured about stage dispersion.** On Gaussian latent matrices at the SEA-AD cell table, with variance multipliers by stage of (1.120, 1.278, 0.802), the pattern of the predictive coordinates on SEA-AD, the unrestricted null read `delta` 0.077–0.079, and 0.091–0.102 with the multipliers squared. With (0.805, 0.921, 1.137), the scores' pattern, it was conservative. Those matrices were drawn in the latent space; nothing went through a fit.
- **Costs.** The univariate baseline's stage-span units took a mean of 24.9 s on `ing` with the scan and split arms; the stage-effect scaling run's medians were 34.7 s and 113.8 s; the ADNI simulation took 85 s per unit for two arms.

## Goals / Non-Goals

**Goals:**
- Read the three candidates against production at paper grade on simulated data, on shared datasets and draws.
- Cover conditions no candidate has been read at: three blocks, four stages, unbalance, unequal noise by stage.
- Return a map by design family, so a candidate that helps at some designs and harms at others is described as such.

**Non-Goals:**
- Reproducing the SEA-AD excess. Its regime (125 features per donor) is outside the generator.
- A rule that recommends a candidate from a cohort's design. The map is what such a rule would be written from.
- Power under unequal dispersion.

## Decisions

**D1. One dataset per unit, four arms, two fits shared where they can be.**
- **The unit.** Generate once. Fit `pls` in scores with the production parameters and evaluate `pls_cv`. Build `predictive_orthonormal` from that fit with `pls_coordinates_from_pls` and evaluate it with `latent=`. Fit and evaluate `weighted`, and `sparse_block` where the design is in its subset, each with the production parameters changed only in `integration_params`. Every arm is evaluated serially at the unit's RRPP seed.
- **Why `weighted` and `sparse_block` cannot share the fit.** They change the fit itself. The production selection ran them the same way.
- **Consistency checks**, recorded per unit: the predictive arm's rank and AUROC equal `pls_cv`'s; at a fit rank of at most K − 1 its p-values equal `pls_cv`'s up to the metric (recorded, not required: the orthonormal metric differs from the scores').
- **Alternative: add the candidates to the univariate baseline's unit.** That unit also runs the scan, the split test and a PCA arm, which this question does not need, and a committed driver is not edited (its children bind it by hash).

**D2. Committed designs are rerun at their seeds, and production must reproduce.**
- **Why rerun.** The paired tests need per-unit records for `pls_cv`, and committed aggregates hold counts only. Rerunning costs one production fit, which the predictive arm needs anyway.
- **What it buys.** An exact gate on the generator, the fit and the draws in three profiles, and a map that lines up cell for cell with the univariate baseline's and the stage-effect scaling run's.
- **Why these datasets are admissible for a gain reading.** No candidate has been evaluated on them. The gains to confirm were read at base seeds 2035 and 2036.
- **Which rows are compared.** `delta`, `angle` and `shape` singly on every cell: 39 rows per thirteen-cell design, 33 per eleven-cell one. On `adni_design`, both arms' rows.
- **Alternative: fresh datasets everywhere.** It loses the gate and the alignment and saves nothing.

**D3. New conditions are new design names under base seed 2043.**
- **Designs.** `n200_four_stages` and `n200_delta0.5_four_stages`; `seaad_table` and `seaad_table_delta0.5` (F 11 / 10 / 28, M 9 / 10 / 12, through `group_stage_sizes`); `n200_share0.8`. All three blocks, `joint`, `p_dmp` 0.1, ρ = 0, the stage-span construction, planted scale on the in-span cells.
- **Why these.** Four stages is the pilot's stated limitation: the measured dimension is K − 1 = 3 and the CV rank relates to it differently. The SEA-AD table is the case-study cell table and the one unbalanced design with paper-grade production readings. A share of 0.8 at n = 200 held size at equal dispersion (item 26); n = 80 at that share did not and is under diagnosis (open item 16), so it is left out.
- **Seeds.** The comparison's `generator_seed` and `rrpp_seed` at base seed 2043. A test checks that no new design name occurs in any committed settings file, the group-dispersion run's included once it is on `main`.
- **Out-of-span cells at four stages.** The rotation is defined there; an out-of-span bend is not (roadmap open item 12) and is not a cell.

**D4. Stage dispersion is a generator setting, and its profile is null only.**
- **The field.** `stage_dispersion_scale`: a tuple of `n_stages` positive finite SD multipliers. Each group's call passes its per-cell scales, the stage multiplier times the group multiplier where `group_dispersion_scale` is set. `generate_omics(noise_scale=)` accepts a sequence with one value per cell; a scalar keeps its code path.
- **Pooled scales.** The within-cell multiplier becomes Σ_c n_c·m_c² / Σ_c n_c over the realized cells, with m_c the cell's combined multiplier. With a group scale alone this is 0.26.0's `(n_A + s²·n_B) / (n_A + n_B)`.
- **Patterns**, as SD multipliers by stage in stage order:

  | pattern | multipliers | source |
  |---|---|---|
  | `last_narrow` | (1.058, 1.130, 0.896) | square roots of item 20's `predictive` variance multipliers |
  | `last_wide` | (0.897, 0.960, 1.066) | square roots of its `scores` multipliers |
  | `last_narrow_strong` | (1.120, 1.278, 0.802) | its `predictive_doubled`, as SDs |
  | `last_wide_strong` | (0.805, 0.921, 1.137) | its `scores_doubled`, as SDs |

  On the SEA-AD table the last stage is the largest (40 of 80), so `last_narrow` is the pattern under which item 20's unrestricted null over-rejected.
- **Designs.** All four patterns on `seaad_table` and `seaad_table_delta0.5`; the two strong ones on balanced n = 200, where no stage is the largest and the patterns separate dispersion from unbalance. Their equal-dispersion references are the `conditions` designs and the committed `n200`.
- **Why null only.** Size is what the real data put in question. A plant's power under unequal noise mixes the measurement with a changed signal-to-noise ratio by stage and has no committed reference.
- **Why this is not item 20 again.** Item 20 imposed the pattern on the latent matrix. Here the pattern is in the features, and each arm's fit decides what reaches the statistics. A candidate that concentrates the fit's in-sample directions may turn a feature-level pattern into a different latent one. Each arm's residual dispersion by stage in its measured coordinates is recorded per design, as the diagnosis recorded it.
- **Recorded arms.** `pls_cv` and `predictive_orthonormal` under `permutation_strata="stage"`, on the same latent matrix and seed. Item 20's layer 2 found `predictive_orthonormal` qualifying under those draws at 100 replicates; these are its first readings where the draws are expected to matter.
- **Alternative: reuse item 20's Gaussian matrices.** They cannot carry a candidate, because a candidate is a different map from features to coordinates.

**D5. The surgery and ADNI-design profiles.**
- **`surgery`.** The slice's three-block columns at n = 80 and 200, eleven cells at 500 × 999, through the univariate baseline's builders. The pilot read one surgery cell (Phase 6 orientation at e = 1.0, `angle` 0.66 → 0.48 for the predictive arms); these are the same plants on three blocks at paper grade. n = 500 is left out: both sides saturate there in the baseline's map.
- **`adni_design`.** The ADNI gain run's simulated units through its builders, with the `weighted` arm added by changing `integration_params` alone. `sparse_block` is not run: the CV selects the cap of 19 on 1,501 samples.
- **Why the ADNI design stays in.** It is a simulated design, and leaving out the one simulated family where a candidate is known to fail would choose the answer. The proposal states the consequence for `predictive_orthonormal`.

**D6. The decision is scoped by family.**
- **Families.**

  | family | designs |
  |---|---|
  | `core` | `n80`, `n200`, `n500` |
  | `weak_stage` | `n200_delta0.5`, `n320_delta1`, `band_n320` |
  | `four_stages` | `n200_four_stages`, `n200_delta0.5_four_stages` |
  | `unbalanced` | `seaad_table`, `seaad_table_delta0.5`, `n200_share0.8` |
  | `adni_design` | the simulated ADNI table |

- **Why scoped.** The roadmap's own reopening condition asks for "a rule that scopes admissibility by design and is declared before the reading". One global label would repeat the outcome this change was opened to revisit: a single family's failure hiding where a candidate works.
- **Criteria.** S1, S2, S3 and G are the predictive pilot's (its design D7) with `judge_rates` and the exact paired tests imported unchanged, applied within a family. S1 reads the anchors and translations. S2 and G read the in-span target cells (magnitude → `delta`, orientation at φ = 0 → `angle`, shape → `shape`). S3 reads those cells' two off-target statistics: the bound at δ = 2.0 designs, the paired excess over `pls_cv` below δ = 2.0, where production's own off-target rates exceed the bound (item 24).
- **`weighted`'s `delta`.** It is in the weighted metric. The tests compare rejections, not statistics, so the criteria apply; the report carries the note.
- **D, stage dispersion.** Judged once per arm over the ten designs (m = 60). The fallback to a paired test when production fails is fixed here because the production selection's sparse arms were failed on an anchor production failed too, the reading this change revisits.
- **I, incidental routes.** The surgery target cells and the two out-of-span rotation cells per design. Production reads these through components beyond the stage program, which is not a property the measurement promises (`docs/measurement-visibility.md`, "The incidental routes"). A loss there is stated, with the cells, and does not veto. The labelling rule is the univariate baseline's (0.025 over the cells, margin 0.10).
- **Order and recording.** Size before power before specificity before gain. Every criterion is recorded for every arm and family.
- **`sparse_block`.** Run in parts of `core` and `weak_stage` only, so its best verdict is `scoped_candidate`. A family counts as run for an arm when every design of the family has the arm; otherwise the label is computed on the designs it has and marked `partial`, and `partial` does not count toward `qualifies_everywhere`.
- **Alternative: one global rule, as in the pilots.** Rejected above.
- **Alternative: make I a veto.** It would hold the candidates to what production does by accident. The cells are reported, and a reader who needs those routes can see what is lost.

**D7. What is descriptive.**
- Per design and arm: the selected rank (median, q10, q90), the held-out AUROC, the `pooled_stage` and `group_difference` capture, and the sparse arm's keep fraction and support by block.
- Per cell and arm: the median observed statistic beside the population's through the same map.
- Production's absolute off-target rates below δ = 2.0.
- The within-stage readings and each arm's residual dispersion by stage on the dispersion designs.
- Rates by rank stratum where an arm loses a cell.

**D8. Execution follows the recent validation drivers.**
- One script with `run --profile` / `merge --profile` / `decide`; `settings.json` bound by hash with the parents' sha256; round-robin shards; resume by complete unit key; JSONL records holding seeds and statistics only; an sbatch array with a profile variable, an arms variable and a dry-run variable, BLAS pinned to one thread, a throttled array, run on `ing` from a clean worktree at a pushed revision.
- **Tiers.** (1) `committed`, `conditions`, `dispersion` with `pls_cv`, `predictive_orthonormal` and `weighted`; (2) `surgery` and `adni_design`; (3) `sparse_block` on its three designs, appended to tier 1's records by unit key. Each tier has its own smoke timing, cap and approval. The decision rule is committed before tier 1 and is not amended between tiers.
- Committed: aggregates, the decision, provenance and the report under `results/candidate-simulation-<date>/`. Records stay on the cluster.

## Risks / Trade-offs

- **[The simulation cannot reproduce the regime of the real-noise veto]** → Stated as a non-goal and in the report. The dispersion profile reads the mechanism's simulated analogue; the feature count waits for roadmap open item 15. A candidate that qualifies here is not cleared for a case study at 125 features per donor.
- **[Feature-level stage dispersion does not produce the latent pattern item 20 used]** → Each arm's residual dispersion by stage in its measured coordinates is recorded, so the report says what the pattern became. A profile in which nothing fails is a valid outcome, reported as such.
- **[Production fails size under stage dispersion]** → Predeclared: it is reported as `production_dispersion_excess`, and the candidates are judged against it by pairing.
- **[A scoped verdict is read as a recommendation]** → The consequence text of `scoped_candidates` says the arms stay opt-in and that a rule recommending one from a cohort's design is a separate change.
- **[The planted fixed point fails at four stages]** → It cycled in 12 of 1,000 four-stage orientation units at paper grade. Failed units are recorded and counted per cell; the smoke runs every design; a design whose failures exceed a few units is amended before any array unit.
- **[Compute]** → About three times the largest run so far. Tiers isolate the cost; the levers, in order: `sparse_block` (tier 3), the `surgery` profile's n = 80 column, `n500`, and the mild dispersion patterns.
- **[Reproduction fails from environment drift]** → The gate stops the reading. A mismatch is diagnosed, never waived.
- **[The roadmap's reopening condition asks for real-noise size]** → The change amends that condition in the roadmap, dated, as the user's decision of 2026-11-25, and keeps every real-noise reading in place as a limit.

## Migration Plan

Nothing migrates. The generator setting is opt-in and absent by default.
- **On `propose_default_change`:** a separate change flips the named arm's default in the harness and `motco analyze`, with pins for every committed configuration and a re-issue of the case studies.
- **On `scoped_candidates`:** the arms stay opt-in; the measurement-visibility page carries the map; a rule that recommends an arm from a cohort's design is a separate change.
- **On `keep_production`:** the arms stay opt-in diagnostics, now read at paper grade in simulation.
- **Rollback:** remove the field. Nothing reads it while absent.

## Open Questions

- The core-hour cap of each tier is set from its smoke timing, at its approval gate.
