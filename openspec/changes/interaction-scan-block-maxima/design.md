# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **The joint scan.** `interaction_feature_scan` orders the non-degenerate features by observed F (descending, ties by column), draws RRPP's serial Freedman–Lane sequence, and per draw forms every feature's null F from one W'ᵀ·E product. It counts, at every rank, the draws whose suffix maximum reaches the observed F; the adjusted p-values are the running maximum of those counts. It keeps every draw's maximum and its column (`null_max_f`, `null_max_column`), and the B × p null F only on request. One draw costs O(n · rank(F) · p) for the product and O(p) for the counting.
- **Blocks are already named.** `analyze_trajectories` builds `X` with columns `block__feature` in declared block order, and `_split_feature_names` recovers the blocks; the record's per-block hit counts and null-maximum shares use them. The simulation's feature columns follow the same convention.
- **The scan does not need the latent space.** Its inputs are the analysis's unweighted pooled-standardized features (`_unweighted_features`) and the covariate-adjusted full and reduced models (`_covariate_adjusted_design`). In the parent validation's simulation units they were the PLS fit's input coordinates, `(x − x̄)/s` of the same features, a per-feature affine image that the F ignores. Fitting accounted for 10.5 s of a median 11.5 s simulation unit; the scan for 0.07 s (`results/interaction-scan-2026-11-09/simulation/accounting.json`).
- **The parent run.** The interaction-scan validation (`results/interaction-scan-2026-11-09/`, 8,100 units, 0 failures, both of its parents reproduced) commits per profile its settings, `scan_table.csv`, `block_balance.csv` and `false_hit_table.csv` (simulation) or `recovery_table.csv` (plasmode). Its per-unit records were deleted from the cluster after the merge, so reproduction is against these aggregates.
- **The export lives on `ing`** (`/home1/tgonza/plasmode-data/cohort_30`), where `production-latent-space-selection`'s arrays are running under an approved cap.

## Goals / Non-Goals

**Goals:**
- A block-wise combination whose validity follows from the joint scan's (no new assumption), that gives every declared block the same share of the familywise error, and that costs nothing beyond the joint scan.
- Identity with the joint scan for one block, and byte identity without labels.
- A validation whose units run only the scan, paired arm by arm on the same inputs and draws, with the joint arm reproducing the parent exactly.

**Non-Goals:**
- Choosing the scan's default combination, or a case study's reading of it.
- A global test that is more powerful than the Bonferroni split under dependence between blocks (measured here, not adopted).

## Decisions

**D1. Blocks share the familywise error equally by Bonferroni over within-block step-down max-T.**
- **Definition.** Given a block label per column, K is the number of blocks with at least one non-degenerate feature. Within block k, the features are ordered by observed F (descending, ties by column index), and the joint scan's step-down rule is applied to the block alone: at within-block rank j the raw value is (1 + #{draws whose maximum null F over the block's features ranked j or lower ≥ the j-th observed F}) / (1 + B), and p⁽ᵏ⁾ is the running maximum. A feature's adjusted p-value is min(1, K · p⁽ᵏ⁾). Block k's own global p-value p_k is its top feature's p⁽ᵏ⁾, and the scan's global p-value is min(1, K · min_k p_k), which is also the smallest adjusted p-value.
- **Why it is valid.** Each block's step-down controls the familywise error within the block strongly at any level, under the same subset pivotality the joint scan relies on. Testing each at α/K bounds the familywise error over all features by α (union bound), whatever the dependence between blocks.
- **Why equal weights.** The principle is that every declared block has the same chance of a false hit, whatever its size or noise. Weights by feature share would restore the count imbalance this change removes; caller-chosen weights are a non-goal.
- **Resolution.** A hit needs p⁽ᵏ⁾ ≤ α/K, 0.025 at K = 2: at most 24 exceeding draws out of 999. The combination needs B ≫ K/α (40 at K = 2), not the B ≫ p/α that ruled out per-feature minP (parent design D2).
- **Identity.** With one block, K = 1 and every quantity equals the joint scan's. A block whose features are all degenerate gets no maximum and no p-value and does not count in K.
- **Alternatives.**
  - **minP between blocks.** Each block's maximum is converted to its tail probability under the block's own null maxima, and the smallest is compared with the draws. It exploits the dependence between block maxima: at K = 2 the per-block level is 0.0253 for independent maxima (Bonferroni's 0.025) and up to 0.05 for fully dependent ones. Its global test is cheap. Its per-feature step-down is not: calibrating every step against the full block's null maxima cannot resolve features whose F lies below every null block maximum (the joint step-down counts the suffix maxima of the remaining features there instead), and recomputing each block's suffix-maximum distribution at every step needs the B × p null F (80 MB at 10,000 features, 3.6 GB at 450,000). It is measured descriptively (D6), so the report states what the dependence would have bought.
  - **Holm across blocks.** Passing α on from a block is valid only once every hypothesis in it is rejected (gatekeeping), not after one hit; for sparse differences it gains nothing.
  - **Per-feature minP.** Infeasible at B = 999 (parent design D2).

**D2. The block-wise path costs one pass and keeps the joint quantities.**
- **Order.** The non-degenerate columns are ordered once by (block order, F descending, column index); the blocks keep the order of their first column.
- **Per draw.** The joint scan's W'ᵀ·E product and per-feature F* are unchanged. Within each block segment a reverse cumulative maximum gives the suffix maxima, compared with the ordered observed F. That is O(p) per draw. The draw's maximum per block is kept (B × K), and so are the joint maximum and its column, so the per-block share of null maxima stays computable in both modes.
- **Without labels** the joint code path runs untouched, so its outputs are byte-identical to 0.20.0's.
- **Validation before any draw.** The labels must have one entry per column; a mismatch is an error naming both lengths.

**D3. The result.**
- `InteractionScanResult` gains `combination` (`"joint"` or `"block_bonferroni"`), the block names in order, `n_blocks_tested` (K), and per block the maximum F, its column and its own global p-value p_k, plus the B × K null block maxima.
- In block mode, `p_value_global` and `p_value_adjusted` are the by-block values, so `hits()` returns the by-block hits.
- Without labels the new fields are absent from the serialized result, and every existing field is unchanged.

**D4. The analysis and the CLI.**
- **The setting.** `TrajectoryAnalysisParams.interaction_scan_by_block` (default `False`) requires `interaction_scan`. Without it, the parameters are rejected in `__post_init__`, before any fitting, naming both. The resolved parameters carry it only when true. The labels are the declared block names.
- **The record.** The `interaction_scan` block gains `combination` and `n_blocks_tested`. Its per-block entries gain `max_f` and `max_feature`, and with permutations `p_value` (the block's own global p-value). `p_value_global` and the hit list are the by-block ones. Without the setting the block is unchanged.
- **The table.** `interaction_scan.csv` keeps its columns; its `p_value_adjusted` holds the by-block adjusted p-values.
- **The CLI.** `motco analyze --interaction-scan-by-block` sets it, exits cleanly without `--interaction-scan`, and the printed scan line names the combination and K.
- **The input helper.** A private helper builds the scan's inputs from the dataset and parameters without the latent fit: the unweighted pooled-standardized features, the covariate-adjusted full and reduced models, the LS-mean rows, the contrast and the stage levels. `analyze_trajectories` takes the scan's inputs from it, so the driver's scan-only plasmode unit computes the analysis's own matrices. The identity fixtures pin the analysis unchanged, and a test pins the helper against the analysis's scan bit for bit.

**D5. The validation runs scan-only units in three profiles, two arms each.**
- **Arms.** `joint` (no labels) and `by_block` (labels from `block__feature`), on the same inputs at the unit's seed. Rejection is adjusted p ≤ 0.05 in both, the parent's scan rule.
- **`simulation`.** The settings are rebuilt through `interaction_scan_check.simulation_settings` and must equal the committed `results/interaction-scan-2026-11-09/simulation/settings.json` (differing keys named). A unit generates `lsc.unit_dataset`, takes the harness preprocessing's pooled-standardized features on the measured blocks and the trajectory models from the design, and runs both arms at the unit's RRPP seed without a correction, as the parent did. The population-null mask is the parent's `population_null_mask`.
- **`plasmode`.** The settings come through `interaction_scan_check.plasmode_settings`, equal to the committed ones; the export's file hashes must equal the plasmode's settings. A unit is `plasmode_check.build_unit`, the helper of D4 at the parent's parameters, and both arms with `df` at the analysis's seed. The fidelity gate and the anchor latent check are not repeated: they check the latent space, which the scan does not use, and the export is bound by its hashes.
- **`plasmode_unequal`.** The export as loaded under the contract, with the RNA block reduced to 500 of its 5,000 features. The subset is the first 500 of a permutation of the RNA columns, kept in export order and recorded by the sha256 of the feature list. The permutation's seed is derived as the plasmode derives its seeds, by sha256 of the driver's name, base seed 2034 and the purpose `unequal_subset`, so it is fixed for the profile and shared by every unit. RNA is reduced because its null F has the lighter tail, so the joint scan disadvantages it by both count and tail. Units are the plasmode's (anchor, e = 2.0 and 3.0 × 500, 999 permutations, `paper` stream seeds), with the planted features drawn from the reduced block by the plasmode's own unit builder. It has no parent.
- **Records.** Per arm: the global p-value, hits per block, planted hits and hits outside the planted set per block (plasmode profiles), and hits on population-null features per block (simulation). Per unit: the per-block null-maximum shares (joint), the 20 largest F's planted count (identical in both arms, recorded once), the Spearman correlation of the block maxima across draws, and the minP global p-value of D1, computed from the null block maxima. Counts only: no feature name, donor identifier or label vector.
- **Runs** are sharded, resumable and bound to their settings by hash, as in the parent.

**D6. Reproduction and the decision.**
- **Reproduction.** The `joint` arm must reproduce, exactly: in the simulation, `scan_table.csv` (count, rejections, median hits), `block_balance.csv` (mean hits and null-maximum share) and `false_hit_table.csv` (count, rejections); in the plasmode, `scan_table.csv`, `block_balance.csv` and `recovery_table.csv` (false-hit count and rejections, recall, hit precision, top-20 precision). The sha256 of each parent file is recorded, and a changed file is refused by name in `run`, `merge` and `decide`. A failure is `reproduction_failed`, which gates every criterion.
- **S1, simulation size (m = 6).** The by-block global rate at the anchor and translation of `n200`, `seaad` and `seaad_weak`.
- **S2, real-noise size (m = 2).** The by-block global rate at both exports' anchors.
- **S3, strong familywise control (m = 4).** The by-block rate of units with a hit outside the planted set at e = 2.0 and 3.0 on both exports.
- **S1–S3** use the stage-span paper grade's `judge_rates` by import.
- **P, protected power (m = 4).** At e = 2.0 and 3.0 on both exports, b counts the units where only `joint` rejects globally and c those where only `by_block` does. The exact one-sided loss p-value is P(Binomial(b + c, ½) ≥ b). `power_loss` if any is below 0.05/4.
- **G, allocation gain.** On the unequal export, per unit, d is the by-block arm's planted hits in the reduced block minus the joint arm's. The exact one-sided sign test on the non-zero d must give p ≤ 0.025 at e = 2.0 or at e = 3.0.
- **Precedence:** `reproduction_failed`, `size_excess`, `fwer_excess`, `power_loss`, `no_allocation_gain`, `by_block_validated`.
- **Consequence.** Under `by_block_validated` the option is documented as validated for exports whose blocks differ in size or noise, and the joint scan stays the default. Otherwise the option stays opt-in and is documented with the failing criterion and cell.
- **Descriptive tables, which never change the verdict:**
  - block balance per arm: the joint arm's null-maximum shares, and for both arms the rate of units with a false hit in each block at the anchors (by block, each near α/K);
  - recall and hit precision per block per arm, against the parent's committed recovery;
  - the minP comparison: per cell, the global rate under minP beside Bonferroni's, and the median correlation of block maxima;
  - simulation localization and population-null false hits per arm.

**D7. Documentation.**
- **The measurement-visibility page's scan section** gains when to scan by block (blocks of unequal size, or noise known to differ), what a by-block hit claims (the same familywise error over all features, now split equally across blocks), and the measured transfer of detection between blocks.
- **API pages:** `docs/api/sd.md` documents the labels and the result fields, and `docs/api/analysis.md` the setting, the record fields and the flag.
- **Roadmap and repository guide:** roadmap item 16's thread 3 follow-ups record the verdict, and `AGENTS.md` gains the setting and the driver.

## Risks / Trade-offs

- **[Bonferroni ignores dependence between block maxima]** → At K = 2 the per-block level is 0.025 against at most 0.05 under full dependence. The minP table measures the realized cost. If it is material, a minP global test is a follow-up.
- **[The larger or noisier block loses detection]** → P protects the global rate on both exports, and the per-block recall table shows the transfer. The documentation states the trade.
- **[Many blocks]** → Each of K blocks gets α/K, so with many small blocks the split costs more than it gives. Only K = 2 is measured, and the documentation says so.
- **[The simulation's joint arm reproduces an affine image of the parent's inputs]** → Every F is mathematically identical (tested to 1e-10). Rounding can flip a count only at an exact tie with a null value. The gate stays exact, and a mismatch is investigated per unit, not waived.
- **[The unequal export is a construction]** → It is a random subset of real features with the real noise, so it tests count imbalance with RNA's tail. A small ATAC block is not run.
- **[Scan-only units skip the parent's latent checks]** → The scan never reads the latent space. The export is bound by its hashes, and the joint arm's exact reproduction checks the inputs.
- **[Subset pivotality within a block is approximate under Freedman–Lane with covariates]** → As in the parent; S3 measures strong control on real noise in both exports.
- **[Compute estimate]** → The local smoke times all three profiles, and the projection in `notes.md` replaces the estimate before launch.

## Migration Plan

Opt-in throughout, so there is nothing to migrate.
- **Version:** MOTCO 0.23.0.
- **Rollback:** without labels the primitive equals 0.20.0's, and without the setting every analysis output equals 0.22.0's.
- **Branch:** stacked on `proposal/production-latent-space`; rebased onto `main` after PR #86 merges.

## Open Questions

- **A case study's reading.** Whether `sead_sexdiff` reads the scan by block, and how that enters its multiplicity policy, belongs to its own change and changes nothing here.
