## Why

The per-feature interaction scan (MOTCO 0.20.0, `scan_validated`) takes one maximum F over every feature, so its familywise bar is set by whichever block tends to give the largest null F. The validation measured that the blocks do not share the bar evenly (`docs/reports/interaction-scan-2026-11-09.md` §5, `results/interaction-scan-2026-11-09/`):
- **On the SEA-AD export** each block holds half the features, but ATAC gives 0.618 of the null maxima. At the anchor ATAC has 0.040 false hits per unit against RNA's 0.022, and at e = 2.0 it has 0.67 hits per unit against RNA's 0.58. ATAC's null F has the heavier tail, so RNA features face a bar set mostly by ATAC's noise.
- **In the simulation** methylation holds 0.737 of the features but gives 0.50 of the maxima, because its correlated features act as fewer tests.

Validity is unaffected. The imbalance decides which block's features can be named. It matters most when the blocks differ in size. SEA-AD's export forced 5,000 + 5,000 features. An export with, say, 450,000 CpGs beside 2,000 proteins would give the CpGs nearly all of the null maxima, and a protein would have to clear a bar sized for 452,000 tests. The scan's design named per-block maxima as the follow-up if the balance was materially off (`interaction-feature-scan` design, Risks and Open Questions), and on the export it is (0.62 / 0.38).

The scan reads the unweighted pooled-standardized features and uses no latent space, so this change does not depend on `production-latent-space-selection`'s verdict.

## What Changes

- **A block-wise combination in `interaction_feature_scan`.** It is opt-in, through the block label of every column.
  - **Within each block**, the scan runs the step-down max-T it already runs, on the block's own features, with the same draws.
  - **Across blocks**, the familywise error rate is split equally: each of the K blocks that has a non-degenerate feature is tested at α/K (Bonferroni). A feature's adjusted p-value is K times its within-block adjusted p-value, capped at 1. The global p-value is K times the smallest block global p-value, capped at 1.
  - **The result** gains, per block: the maximum F and its column, the block's own global p-value, and the null maximum of every draw.
  - **With one block label** the scan equals the joint scan. Without labels it is byte-identical to 0.20.0's.
- **An opt-in in `motco analyze`.** `TrajectoryAnalysisParams.interaction_scan_by_block` (CLI `--interaction-scan-by-block`) requires `interaction_scan` and takes the blocks from the declared omic blocks.
  - **The `interaction_scan` record block** gains the combination, the number of blocks tested and, per block, the maximum F with its feature and the block's own global p-value.
  - **The scan table's adjusted p-values and the hit list** are the by-block ones.
  - **Without the setting** every output equals 0.22.0's, and the resolved parameters carry the setting only when it is on.
- **A validation driver.** `scripts/interaction_scan_blocks_check.py` (`run --profile simulation|plasmode|plasmode_unequal`, `merge`, `decide`) runs every unit through two arms on the same inputs and draws: `joint` (the 0.20.0 scan) and `by_block`. A unit runs only the scan. It builds the scan's inputs without fitting a latent space, so the run is cheap.
  - **`simulation`:** the interaction-scan validation's 6,600 simulation units (three designs, fourteen cells, base seed 2035, 498 features).
  - **`plasmode`:** its 1,500 SEA-AD plasmode units (anchor, e = 2.0 and 3.0 × 500 at 999 permutations, base seed 2034, `paper` stream, `df`).
  - **`plasmode_unequal`:** the same export with the RNA block reduced to a seeded random 500 of its 5,000 features. It runs the anchor, e = 2.0 and 3.0 × 500 at 999 permutations, with 10 planted features per block. This is the unequal-size case the change is for.
  - **Reproduction:** the `joint` arm must reproduce the parent's committed scan, block-balance, false-hit and recovery tables exactly in counts.
- **A predeclared decision.** `interaction_scan_blocks_decision.json` (rule version 1), with `reproduction_failed` gating every criterion:
  - **S1:** simulation size of the by-block global test at the anchor and translation of every design (m = 6).
  - **S2:** real-noise size at the anchors of both exports (m = 2).
  - **S3:** strong familywise control: units with a hit outside the planted set at e = 2.0 and 3.0 on both exports (m = 4).
  - **P:** protected power: the by-block global rate is not below the joint rate at e = 2.0 and 3.0 on both exports (exact paired loss test, m = 4).
  - **G:** allocation gain: on the unequal export, the by-block arm names more planted small-block features than the joint arm (exact paired sign test, at e = 2.0 or 3.0).

  S1–S3 use the stage-span paper grade's count-plus-ceiling rule. The verdict is `size_excess` | `fwer_excess` | `power_loss` | `no_allocation_gain` | `by_block_validated`. Descriptive tables cover per-block null-maximum shares and false hits per arm, recall and precision per block per arm, and how much a minP combination between blocks would have gained over Bonferroni.
- **Documentation.**
  - **`docs/measurement-visibility.md`:** the scan section gains when to scan by block and how a by-block hit is read.
  - **API pages:** `docs/api/sd.md` and `docs/api/analysis.md`.
  - **Report:** a dated report with aggregates under `results/`.
  - **Roadmap and repository guide:** roadmap item 16's thread 3 follow-ups and `AGENTS.md`.
  - **Version:** MOTCO 0.23.0.

**Non-goals.** Each is a separate change:
- making the by-block combination the default. A default flip would change `sead_sexdiff`'s pinned reading and needs its own change after the verdict;
- minP between blocks as the procedure. It is measured here descriptively, for its global p-value only (design D2);
- weights other than equal, chosen by the caller;
- a moderated per-feature statistic, and FDR control;
- applying the scan to SEA-AD in `sead_sexdiff`;
- the scan in the simulation harness or study configs.

## Capabilities

### New Capabilities
- `interaction-scan-block-validation`: the by-block combination's validation. It covers the three profiles, the scan-only units, reproduction of the interaction-scan validation, the predeclared decision and the dated report.

### Modified Capabilities
- `interaction-feature-scan`: the global p-value and the step-down adjustment are stated for the joint combination, and a block-wise combination is added. It covers equal Bonferroni allocation over blocks, within-block step-down, per-block results, and identity with one block or none.
- `real-data-trajectory-analysis`: `analyze_trajectories` and `motco analyze` gain the opt-in by-block scan and its record fields, byte-identical when absent.
- `latent-space-metric-compatibility`: the measurement-visibility documentation states when to scan by block and how a by-block hit is read.

## Impact

- **Code:**
  - `src/motco/stats/interaction.py` (the block-wise combination and its result fields);
  - `src/motco/analysis.py` (the parameter, the record fields, and a private helper that builds the scan's inputs, which `analyze_trajectories` itself calls) and `src/motco/cli.py` (the flag and the printed line);
  - `src/motco/__init__.py` (0.23.0).
- **Scripts:** `scripts/interaction_scan_blocks_check.py` and `scripts/interaction_scan_blocks_check_array.sbatch`. They import `interaction_scan_check`, `span_split_check`, `latent_space_comparison`, `plasmode_check` and `stage_span_paper_grade_decision` and change none of them.
- **Tests:**
  - by-block adjusted p-values against a brute-force per-block step-down;
  - identity with one block and without labels;
  - per-block null maxima and degenerate blocks;
  - the analysis's byte identity without the setting, and its rejection without `interaction_scan`;
  - the input helper against the analysis's own scan;
  - driver reproduction, resume and settings binding;
  - every decision path.
- **Compute:** about 9,600 scan-only units. The estimate is about 10 core-hours (simulation ≈ 2, plasmode ≈ 5, unequal ≈ 3), against the parent's 69.6, because no unit fits a latent space. A local smoke replaces the estimate before launch. The simulation profile can run on the workstation. The plasmode profiles read the export on `ing` and run after the `production-latent-space-selection` arrays drain, so they do not compete with that run's core-hour cap.
- **Data:** the plasmode profiles read the SEA-AD export in place. The files' sha256 must equal the plasmode's committed settings. The unequal subset is recorded by the sha256 of its feature list. Records hold seeds, statistics and counts only, never a donor identifier, label vector or donor-level value.
- **Compatibility:** every committed configuration, fixture and report is unchanged. Both new settings are off by default.
- **Branch:** this change is stacked on `proposal/production-latent-space` (PR #86) because it touches the same analysis and CLI code. It rebases onto `main` once #86 merges.
