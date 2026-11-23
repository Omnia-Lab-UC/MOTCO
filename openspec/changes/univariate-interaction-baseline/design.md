# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- **What the scan returns.** `interaction_feature_scan` gives, per feature, the interaction F, with the numerator and denominator degrees of freedom, the permutation p-value of the largest F and step-down adjusted p-values. It gives no unadjusted per-feature p-value. Degenerate columns carry no F and enter no maximum.
- **What the scan's design set aside.** It rejected Benjamini–Hochberg on permutation p-values for resolution, and on parametric p-values because nothing validated F tails on real counts (`interaction-feature-scan` design D2). False discovery control was a non-goal of both scan changes.
- **The parents.**
  - The scan validation's simulation profile reran the split-test validation's units (the latent-space comparison's designs and seeds, base seed 2035, two blocks, null scale). It committed, per design and cell, the counts of the scan, the three split parts and the three trajectory statistics (`results/interaction-scan-2026-11-09/simulation/comparator_table.csv`). Its per-unit records are gone from the cluster.
  - The design-guidance paper-grade slice ran the feature-set surgeries through the study engine at n ∈ {80, 200, 500}, three stages, two or three blocks, 500 × 999 (base seed 1028). Its `design_point_operating.csv` tables hold available and rejected counts per column, mode, effect and statistic.
  - `stage_span_leak_check` rebuilds a study-engine column's units (generator and evaluation parameters, replicate seeds) from a committed config.
- **Three blocks.** With `layers` and `span_layers` absent, the harness measures methylation, expression and proteomics and the stage-span plant spans them: 367 + 131 + 160 = 658 features. The stage-span paper grade ran one such column (`n200_three_blocks`).
- **The planted scale cannot be combined with an alignment above 0**, so out-of-span cells run under the null scale.
- **Costs.** A production unit at 999 permutations took a median 12 s in the slice (the cross-validation dominates), and the scan validation's unit, which adds the split test and the scan, 11.5 s on two blocks. The latent-space comparison's five-arm unit took 38 s at 199 permutations.

## Goals / Non-Goals

**Goals:**
- Judge the production trajectory test against the default per-feature analysis on matched units, under a rule fixed before any unit.
- Give each side one test at the same level, so that neither wins by having more tries.
- Read the comparison on the blocks the intended cohorts hold, for the two construction families MOTCO's power claims rest on.
- State separately what only production can claim: the kind of a difference.

**Non-Goals:**
- Ranking the three baseline combinations against one another. They are read as the baseline's alternatives.
- Resolving small differences. A label needs a material difference (D6).
- Four-stage designs and the Phase 5 design point (n = 1,200). The design-guidance columns are the sample sizes the planned analyses have.

## Decisions

**D1. The baseline is the per-feature interaction F, combined three ways.**
- **The statistic.** One F per feature for the group × stage interaction, from the same full and reduced models as the trajectory test. It is what a per-feature linear model reports, and it is invariant to each feature's scale.
- **`scan`: the maximum.** The package's validated test, with its permutation p-value. It pays only for the search and targets a difference on few features.
- **`fdr`: any Benjamini–Hochberg discovery at 0.05 on parametric F p-values.** This is the analysis most papers report. As a test of "any difference" it is Simes' test, which gains over the maximum when many features carry a moderate difference, the regime where production led in the committed cells.
  - **Why parametric.** The scan design's resolution argument holds here too: the first Benjamini–Hochberg threshold at 658 features is 7.6·10⁻⁵, below the 1/1000 resolution of 999 draws.
  - **Why that is acceptable here.** The generator's noise is Gaussian, so the F tails are those of the reference distribution. The scan design's objection was to real counts. The arm's size is read at the controls like any other, and it leaves the comparison if it fails.
- **`total`: the pooled features.** The split test's `total` part, the trace pseudo-F over all features. It pays for every noise dimension and is the weakest arm on a concentrated difference.
- **Alternative: a permutation-calibrated Simes statistic.** Valid under any dependence, but it is not an analysis anyone runs, and it needs every draw's per-feature F. It is the follow-up if `fdr` fails size.
- **Alternative: a moderated F** (roadmap open item 6). It is the form many readers expect, but it is not implemented, and its gain is at small n on real variance heterogeneity, which the generator does not have.
- **Alternative: higher criticism or a top-k sum.** Each adds a tuning choice.

**D2. Production enters as one test at 0.05 / 3.**
- **Why an omnibus.** A reader does not know in advance whether a difference is one of size, orientation or shape. Reading three p-values at 0.05 each gives production three tries at a familywise rate near 0.15 against the baseline's one test at 0.05.
- **Why Bonferroni.** It is what a user can do with `motco analyze`'s output today.
- **Recorded beside it:**
  - each statistic at 0.05, so every committed rate has its counterpart;
  - `production_minp`, the smallest of the three p-values calibrated on the joint permutation draws, which shows what Bonferroni costs when the statistics are correlated. It decides nothing.
- **Alternative: the target statistic alone.** It is an oracle: it uses the plant's kind. It is reported, not judged.

**D3. One unit, one production fit, every arm at one seed.**
- **The unit.** It generates the dataset, fits `pls` once with the cross-validation, and evaluates the trajectory statistics with their null distributions. It then runs the split test and the scan on the split-test validation's inputs (the fit's standardized features and the trajectory design's full and reduced models) at the unit's RRPP seed, uncorrected. Both use RRPP's serial draws, so every arm reads the same permutations.
- **`fdr` from the scan's output.** The upper tail of the F distribution at the scan's degrees of freedom gives each feature's p-value; Benjamini–Hochberg runs over the non-degenerate features. Both come from `scipy.stats`, which the package already depends on. Nothing is reimplemented.
- **The fit is not an input to the baseline.** The scan's F is invariant to the per-feature affine map between the standardized features and the fit's input coordinates. A test checks that the scan on features built with no latent fit (`interaction_scan_blocks_check.simulation_inputs`) gives the same F and global p-value.
- **`concat` and `pca_matched`.** Two further evaluations of the same dataset at the same seed, with the unit's production parameters changed only in the integration method (`concat`; `pca` at the unit's production rank), so they measure the same blocks. The comparison's `arm_params` is not used for them, because it fixes the measured blocks at two. They put the earlier baselines in the same table. They were judged in the comparison (`document_tradeoff`) and are not judged again.
- **Production parameters per profile.** `parent` uses the split-test validation's `production_params`. `stage_span` uses the same cross-validation settings with `layers` absent. `surgery` uses the study cell's evaluation parameters. Each keeps the null distributions, which `production_minp` needs and which change no statistic or p-value.
- **A failure** in the fit or any arm is recorded for the unit without aborting the shard.

**D4. Each arm keeps the rejection rule it was validated under.**
- `scan`: global p ≤ 0.05. `total`, `span`, `complement` and the single trajectory statistics: p < 0.05. `fdr`: at least one adjusted p ≤ 0.05. `production`: smallest p ≤ 0.05 / 3. `production_minp`: p ≤ 0.05.
- **Why not one rule.** The parents' counts must reproduce exactly, and at 999 draws the two forms differ only at a p-value of exactly 0.050.

**D5. The profiles.**
- **`parent` gates and describes.** Its units are the ones the question was read from, so a verdict on them would confirm a pattern on the data that suggested it. It reproduces the seven committed tests exactly, which ties this driver's scan, split parts and production arm to the validated ones, and it adds `fdr` and `production` to the cells the roadmap quotes. The weak design is left out: no test had power there, and it is calibrated to a convenience export.
- **`stage_span` is a fresh family on three blocks.**
  - **Designs.** Balanced n = 80, 200 and 500 are the design-guidance sample sizes, so the two construction families are read at the same designs. The SEA-AD table is not a design here.
  - **The weaker stage effect.** δ = 0.5 is one of the slice's stage-effect values at n = 200. On two blocks the held-out AUROC stayed near 1 there while the CV rank inflated and `delta` and `shape` lost power. The baseline does not use the stage fit, so this is the design where the two are most likely to separate.
  - **Cells.** The in-span ladder is the predictive-subspace pilot's (magnitude at e ∈ {0.05, 0.1, 0.25}; rotation and bend at e ∈ {0.1, 0.25, 0.5}), under the planted scale, where a plant is exact under the sample scaler. The out-of-span cells are the comparison's: a rotation at e = 0.5 and φ = 90, and its 20-feature form, under the null scale.
  - **Grade.** Controls at 500 × 999, where the size rule has a per-cell bound of 0.0695. Signal cells at 200 × 999.
- **`surgery` reuses the slice's three-block units.**
  - **Why reuse.** Production's rates there are committed at 500 × 999, so the `pls_cv` arm reproduces them exactly and the baseline is read on the datasets behind the published power figures.
  - **Why all three sample sizes.** Production is strongest at n = 500 with three blocks (`delta` 0.83–1.00, `angle` 0.66–0.74), so leaving it out would favour the baseline.
  - **Cells.** The anchor, translation at e = 1.0, and every committed effect of magnitude, orientation and shape.
- **Seeds.** `stage_span` uses base seed 2043, which no driver has used, with the comparison's seed derivations. The generator seed depends on the design and replicate, so the cells of a replicate share their null dataset.

**D6. The decision.**
- **Size first, per arm.** Anchors and translations are nulls for every arm: a translation is a constant group offset, which has no interaction. `production` and each baseline arm are judged in two families (`stage_span`, m = 8; `surgery`, m = 6) under `judge_rates`.
  - An inadmissible baseline arm is dropped and reported: a default analysis that does not hold size is a finding.
  - If `production` fails, the verdict is `production_size_excess`. If no baseline arm is admissible, it is `no_admissible_baseline`.
- **Detection by paired exact tests.** The arms share the dataset and the draws, so the information is in the units where exactly one of the two rejects. `paired_p_values` (imported unchanged) gives the one-sided gain and loss p-values of `production` against a baseline arm.
- **Labels.** With m = 71 signal cells and A admissible baseline arms:
  - `production_leads`: the gain p-value is at or below 0.025 / m against every admissible arm, and production's rate exceeds each arm's by at least 0.10;
  - `baseline_leads`: the loss p-value is at or below 0.025 / (m·A) against at least one admissible arm, and that arm's rate exceeds production's by at least 0.10;
  - `unresolved`: otherwise.
- **Why the two thresholds differ.** Production leads only if it beats every baseline alternative, an intersection, so each test may use the full per-cell level. The baseline leads if any one alternative wins, a union, so the level is divided by A. The chance of any false label over the whole table is at most 0.05.
- **Why a margin of 0.10.** At 500 replicates a paired test resolves differences of a few points, such as 0.99 against 0.95. Those are not a reason to prefer either analysis.
- **The verdict is by existence.**
  - `production_leads`: at least one cell is labelled `production_leads` and none `baseline_leads`.
  - `baseline_leads`: the reverse.
  - `complementary`: both labels occur.
  - `unresolved`: neither occurs.
- **Why existence and not a count.** The committed cells already suggest each side leads somewhere. A count would hide that behind the mix of cells chosen. The map by plant family and design carries the content.
- **The out-of-span cells stay in the family.** Production is known to be blind there and the scan was validated on them, so `baseline_leads` labels are expected in those eight cells. Leaving out known losses would tilt the verdict toward production. The consequence is that `production_leads` is reachable only if the scan fails where it was validated: the live outcomes are `complementary` and `baseline_leads`, and the open question is whether production still leads anywhere once it is held to one test and the baseline has a combination for dense differences.
- **Alternative: the unpaired 2-pooled-SE rule** of the scan and split validations. It ignores the pairing and has no familywise control over 71 cells.
- **Order.** Gate, size, detection. Every criterion is recorded whatever the verdict.

**D7. The characterization grade is production's alone.**
- **What it reads.** Among the 36 in-span `stage_span` cells, those that `production` detects in at least 0.80 of units. For each, its two off-target statistics (for a rotation, `delta` and `shape`) are judged under `judge_rates` as one family.
- **Grades.** `names_kind` when the family holds; `detects_only` when it fails, meaning a rejection names a difference and not its kind; `no_detected_cell` when no cell reaches 0.80.
- **Why the planted scale only.** Off target has an exact truth there. The surgeries are not geometrically pure (roadmap, Not yet established), so their off-target rates are reported and not graded.
- **Why it does not change the verdict.** Detection and characterization are different claims. The baseline cannot make the second at all.

**D8. What is descriptive.**
- **The map:** per plant family and design, the count of cells under each label and each arm's rate.
- **Every arm's rates per cell:** the three single trajectory statistics, `production`, `production_minp`, `scan`, `fdr`, `total`, `span`, `complement`, and `concat` and `pca_matched` as single statistics and as the same omnibus.
- **Discoveries:** per cell, the median number of scan hits and of `fdr` discoveries, and the realized false discovery proportion on population-null features where that mask is defined.
- **Fit:** the selected rank, the held-out AUROC and the span-capture shares.
- **`parent`:** `fdr` and `production` beside the seven committed tests.

**D9. Execution follows the recent validation drivers.**
- One script with `run --profile` / `merge --profile` / `decide`; `settings.json` bound by hash with the parents' sha256; round-robin shards; resume by complete unit key; JSONL records holding seeds and statistics only, with no per-feature array; an sbatch array with a profile variable and a dry-run variable, BLAS pinned to one thread, run on `ing` from a clean worktree at a pushed revision.
- Committed: aggregates, the decision, provenance and the report under `results/univariate-baseline-<date>/`. Records stay on the cluster.

## Risks / Trade-offs

- **[`fdr` does not hold size under the features' correlation]** → Its controls are read at 500 × 999 in both families. It leaves the comparison if it fails, the report says so, and the permutation-calibrated form (D1) is the follow-up.
- **[The planted fixed point fails on three blocks at δ = 0.5]** → Failed units are recorded and counted per cell; the smoke runs every design; a design whose failures exceed a few units is amended before any array unit, with the amendment in `notes.md`.
- **[Both sides saturate, so most cells are `unresolved`]** → The effect ladders start below production's committed saturation, and the smoke reports each cell's rates at a few units. An uninformative ladder is amended before the arrays.
- **[The verdict is read as a statement about real data]** → The report states beside every claim that it holds at 658 simulated features with Gaussian noise. The real-noise leg is roadmap open item 2.
- **[The baseline gains from the generator's simplicity]** → Per-feature tests are exact on Gaussian noise with one variance per feature, and a moderated statistic has nothing to add there. The report says so.
- **[Reproduction fails from environment drift]** → The gate stops the reading. A mismatch is diagnosed, never waived.
- **[Compute overrun]** → A local smoke measures per-unit time per profile before launch. The levers are, in order: the `concat` and `pca_matched` arms, the `surgery` column at n = 500, and the signal replicates of `stage_span`. The projection goes to the user with a cap.

## Migration Plan

Nothing migrates. No package code changes and no default moves.
- **On `production_leads`:** the documentation may state that the trajectory test detects differences the default analysis does not, at these designs and this dimension.
- **On `baseline_leads`:** the documentation states that integration gives no detection advantage at these designs. A case study reads the per-feature test first, and the trajectory statistics are presented as a characterization only under `names_kind`.
- **On `complementary`:** the documentation states the map, and a case study reports both analyses under a multiplicity policy fixed in advance.
- **On `unresolved`:** nothing is claimed either way.
- **Rollback:** remove the driver, its test and its report.

## Open Questions

- The core-hour cap for the arrays is set from the smoke timing, at the approval gate.
