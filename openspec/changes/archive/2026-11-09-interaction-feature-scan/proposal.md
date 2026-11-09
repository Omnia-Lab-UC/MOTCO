## Why

Every test MOTCO has pools the group × stage difference over many dimensions, so none of them sees a difference confined to a few features.
- **`delta`, `angle` and `shape`** see only the part of a difference inside the fitted stage-PLS span. A difference on 20 of SEA-AD's 10,000 features keeps a median 2.2 % of its norm there, the share of a random direction.
- **The split test's complement** sums interaction energy over 9,995 coordinates, and 20 features' energy is lost in that noise.

The SEA-AD plasmode (2026-11-02) planted exactly that kind of difference: a ramp to 3 within-stage SDs on 20 features in one group. At e = 3.0, `delta` / `angle` / `shape` rejected in 0.026 / 0.036 / 0.016 of replicates, and the split test's total / span / complement parts in 0.042 / 0.048 / 0.018. That is the null rate everywhere. Yet the information is in the features: the attribution's residual component ranked the planted features with top-20 precision 0.432, against 0.0023 by chance. The workflow never reaches that ranking, because attribution follows a global rejection that never comes.

Roadmap item 16 named the targeted follow-up, thread 3: a per-feature scan of the interaction. Each feature is tested on its own, and the features are combined by their maximum rather than their sum. A sum pays for every noise dimension (its noise grows like √p). A maximum pays only for the search: the per-feature F bar rises from 3.1 to about 14.7 over 10,000 features. A back-of-envelope calculation at the plasmode's design (Gaussian noise, Bonferroni, not a measurement) gives each planted feature power 0.06 at e = 2.0 and 0.51 at e = 3.0. At least one of the 20 would be detected with probability 0.71 and 1.00.

The method is standard: Westfall–Young step-down max-T over Freedman–Lane per-feature F statistics. What this change adds is consistency with the analysis, so the scan's draws, models, covariates and residual corrections are the ones the trajectory statistics and the split test use. It also measures the scan at the designs and on the real noise where the other tests were measured, with a predeclared decision.

## What Changes

- **A per-feature interaction scan in `motco.stats.interaction`.** `interaction_feature_scan` takes an outcome matrix and the full and reduced model matrices.
  - **Per feature** it returns the interaction F, (SS_int/df_int) / (RSS/df_res), computed on that column alone.
  - **With permutations** it returns:
    - the global p-value of the maximum F;
    - per-feature step-down max-T adjusted p-values, which control the familywise error rate;
    - with `return_null`, the per-draw null maxima and their feature indices.
  - **Draws** are RRPP's serial Freedman–Lane sequence at the seed, and the residual corrections are honoured. Under `df` every F, and so every p-value, is unchanged.
  - **Degenerate features** (zero residual variance) are excluded from the maximum and reported, not given a p-value.
  - **Cost** is O(n · rank(F) · p) per draw, and no n × p permuted matrix is formed.
- **An opt-in in `motco analyze`.** `TrajectoryAnalysisParams.interaction_scan` (CLI `--interaction-scan`) runs the scan on the analysis's standardized features, with the covariate-adjusted full and reduced models and the analysis's seed, permutations and residual correction.
  - **The record** gains an `interaction_scan` block after `measurement` and any `span_split` block. It holds:
    - the feature count, degenerate count and df;
    - the maximum F, its feature and the global p;
    - the features whose adjusted p is ≤ 0.05, with their F;
    - per block, the hit count and the share of null maxima;
    - the seed, the correction and `draws_match_rrpp`.
  - **The artifact** gains `interaction_scan.csv`: one row per feature with its block, F, adjusted p and covariate-adjusted group difference at each stage, in standardized units.
  - **Every other output is unchanged.** Without the setting, every output equals 0.19.0's.
  - **The scan is allowed with latent adjustment and with attribution.** It uses neither the basis nor the bootstrap.
- **A validation driver.** `scripts/interaction_scan_check.py` (`run --profile simulation|plasmode`, `merge`, `decide`) reruns the split-test validation's units (2026-11-05) with the scan added.
  - **`simulation`:** the latent-space comparison's three designs and fourteen cells at base seed 2035, 498 features. The anchor and translation run 500 replicates and the other cells 100, all at 999 permutations.
  - **`plasmode`:** the SEA-AD astrocyte `cohort_30` export under the plasmode's contract, base seed 2034, `paper` stream. The anchor, e = 2.0 and e = 3.0 run 500 replicates × 999 permutations. Each unit is one `analyze_trajectories` with the split test and the scan on.
  - **Reproduction:** each profile must reproduce the split-test run's committed split and trajectory tables exactly.
- **A predeclared decision.** `interaction_scan_decision.json` (rule version 1) reads four criteria, with `reproduction_failed` gating all of them:
  - **S1:** simulation size of the global scan at the anchor and translation of every design (m = 6).
  - **S2:** real-noise size at the plasmode anchor (m = 1).
  - **S3:** strong familywise control under the plasmode's partial null. At e = 2.0 and 3.0, the rate of replicates with any hit outside the planted set is judged (m = 2).
  - **G:** sparse gain. At e = 3.0, the scan's global rate must exceed the best of the six existing tests on the same units by more than 2 pooled SE.

  S1–S3 use the stage-span paper grade's count-plus-ceiling rule. The verdict is `size_excess` | `fwer_excess` | `no_sparse_gain` | `scan_validated`. Descriptive tables cover recovery of the planted set against the attribution's committed precision, the block balance of null maxima, the simulation's localization against the split parts and production statistics, and false hits against population truth.
- **Documentation.**
  - **`docs/measurement-visibility.md`:** a section on the scan's reading. It covers what a hit and a non-rejection mean, that the scan gives no geometry, and how its claim differs from attribution's.
  - **API pages:** `docs/api/sd.md` and `docs/api/analysis.md`.
  - **Report:** a dated report with aggregates under `results/`.
  - **Roadmap and repository guide:** roadmap item 16's thread 3 and `AGENTS.md`.
  - **Version:** MOTCO 0.20.0.

**Non-goals.** Each is a separate change:
- applying the scan to SEA-AD. `sead_sexdiff` does that, pinning this release, only under `scan_validated`, and fixing its multiplicity policy across the scan, the split test and the trajectory statistics first;
- reopening Phase 6 or lifting its attribution leg. The plasmode's verdict (`attribution_leg_failed`) stands, and the scan is a new test with its own claim;
- a moderated (variance-shrunk) per-feature statistic, per-block maxima and FDR control;
- pathway- or set-level aggregation;
- the scan in the simulation harness or study configs;
- thread 2 (a group-aware space refitted per permutation).

## Capabilities

### New Capabilities
- `interaction-feature-scan`: the per-feature interaction scan. It covers the per-feature F, degenerate features, the step-down max-T adjustment and global p-value, Freedman–Lane draws with RRPP's sequence and residual corrections, and the result's contents.
- `interaction-scan-validation`: the scan's validation. It covers the simulation and plasmode profiles, reproduction of the split-test run, the predeclared decision and the dated report.

### Modified Capabilities
- `real-data-trajectory-analysis`: `analyze_trajectories` and `motco analyze` gain the opt-in scan, its record block and its per-feature table, byte-identical when absent.
- `module-contracts`: `motco.stats.interaction` exports the scan beside the split test.
- `latent-space-metric-compatibility`: the measurement-visibility documentation states what the scan tests and how its hits and non-rejections are read.

## Impact

- **Code:**
  - `src/motco/stats/interaction.py` (the scan, sharing the split test's validation and model bases) and the `motco.stats` namespace;
  - `src/motco/analysis.py` (parameter, record block, per-feature table, output file) and `src/motco/cli.py` (flag, printed line);
  - `src/motco/__init__.py` (0.20.0).
- **Scripts:** `scripts/interaction_scan_check.py` and `scripts/interaction_scan_check_array.sbatch`. They import `span_split_check`, `latent_space_comparison`, `plasmode_check` and `stage_span_paper_grade_decision` and change none of them.
- **Tests:**
  - the per-feature F against direct computation;
  - step-down adjusted p-values against a brute-force implementation;
  - draws against `RRPP`, and F's invariance under `df`;
  - degenerate features and input validation;
  - analysis byte identity without the setting, unchanged outputs with it, and the analysis path against the primitive;
  - driver reproduction, resume and settings binding;
  - every decision path.
- **Compute:** about 6,600 simulation units and 1,500 plasmode units, estimated at about 65–70 core-hours on `ing`. The split-test run spent 25.6 and 36.9 core-hours on the same units, and the scan adds seconds per unit. A local smoke replaces the estimate before launch.
- **Data:** the plasmode profile reads the SEA-AD export in place. The files' sha256 must equal the plasmode's committed settings. Records hold seeds, statistics and planted-set counts only, never a donor identifier, label vector or donor-level value.
- **Compatibility:** every committed configuration, fixture and report is unchanged. The analysis option is off by default, and the function is new.
