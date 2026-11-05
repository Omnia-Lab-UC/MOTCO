## Why

The production measurement sees a group difference only through its component inside the fitted stage-PLS span. The latent-space comparison (2026-11-04, `document_tradeoff`) found no space that sees a difference outside that span without losing power on what stage-PLS already sees:
- every low-dimensional space read the null rate at φ = 90° (`pls_cv` `angle` 0.08 / 0.11, capture 0.008 / 0.026);
- the full feature space (`concat`) saw it (1.00 / 0.90), but its `angle` was liberal at the SEA-AD-like stage signal (0.21 / 0.22).

The SEA-AD plasmode reported the same blindness on real blocks: a planted orientation on 20 of 10,000 features went undetected up to 3 within-stage SDs (`angle` 0.036). So the case study's non-rejections (44 analyses) say nothing about a sex × stage difference outside the span. They also cannot say whether such a difference exists at all.

This change does not look for a better space. It asks the out-of-span question with its own test. The group × stage interaction energy in the standardized features is split exactly into the part inside the production span and the part in its orthogonal complement. Each part is tested by Freedman–Lane permutation of a pivotal ratio. Both projections are fixed and never use the group labels, so no basis is refitted.

The test answers three questions:
- **Anything:** is there a group × stage difference anywhere in the measured features?
- **Where:** is it inside or outside the fitted span?
- **A benchmark:** a group-aware space (roadmap thread 2) must beat this test to be worth its per-permutation refit.

## What Changes

- **A split interaction test in `motco.stats`.** A new module, `motco.stats.interaction`, provides `interaction_split_test`. It takes an outcome matrix, the full (group × stage) and reduced (group + stage) model matrices, and an orthonormal basis.
  - It returns, for the parts `total`, `span` and `complement`, the interaction sum of squares, the full model's residual sum of squares and the trace pseudo-F. With permutations, it adds their permutation p-values.
  - `total` = `span` + `complement` holds exactly for both sums of squares.
  - Permutations follow `RRPP`'s Freedman–Lane residual randomization, with the same index sequence as its serial path at the same seed, and they honour `residual_correction`.
  - The pseudo-F is the test statistic. Its ratio makes it robust to the reduced model's residual shrinkage, which a sum of squares alone carries into a liberal null that grows with dimension. The sum-of-squares p-value is returned as a diagnostic.
- **An opt-in in `motco analyze`.** `TrajectoryAnalysisParams.span_split_test` (CLI `--span-split-test`) runs the test on the analysis's own pooled-standardized features, with the fitted PLS span as the basis and the covariate-adjusted models.
  - The record gains a `span_split` block: per part, the dimension, the sums of squares, the F, the p-values, plus the observed in-span share of the interaction energy.
  - Latent scores, rank, statistics, p-values and stage means are unchanged. Without the setting every output equals 0.17.0's.
  - It is rejected with latent adjustment before fitting.
- **A validation driver.** `scripts/span_split_check.py` (`run --profile simulation|plasmode`, `merge`, `decide`) has two profiles:
  - **`simulation`** reuses the latent-space comparison's three designs, fourteen cells, calibrated weak δ and seed family (base seed 2035). The anchor and translation run 500 replicates, every other cell 100 (the comparison's own units), all at 999 permutations. Each unit refits `pls_cv`, evaluates it and runs the split test at the unit's RRPP seed. `pls_cv`'s rates from the first 199 draws must reproduce the comparison's committed rate table exactly.
  - **`plasmode`** reuses the SEA-AD plasmode's export contract, null-making step, plant and seeds (base seed 2034, `paper` stream), at the anchor, e = 2.0 and e = 3.0 × 500 replicates × 999 permutations. Each unit runs `motco analyze` at the case-study configuration with the split test on and attribution off. The trajectory rejection rates must reproduce the plasmode's committed summary exactly.
- **A predeclared decision.** `span_split_decision.json` (rule version 1) reads three criteria:
  - C1: simulation size of all three parts at the anchor and translation at every design;
  - C2: real-noise size at the plasmode anchor;
  - C3: out-of-span gain, the complement's rate for the 45° rotation at φ = 90° against `pls_cv`'s `angle` at both δ = 2.0 designs.

  C1 and C2 use the count-plus-ceiling bound rule. The verdict is `split_validated` | `size_excess` | `no_out_of_span_gain`. `reproduction_failed` gates all three. The decision also writes descriptive tables:
  - localization of each part across φ and modes;
  - the observed in-span share;
  - the comparison's five arms at the shared replicates;
  - sparse-plant power on real blocks;
  - the sum-of-squares diagnostic.
- **Documentation.**
  - **`docs/measurement-visibility.md`:** gains the split test and its reading. A complement rejection means a difference outside the *fitted* span, which at a weak stage signal includes differences along the true program that the fitted span misses.
  - **API pages:** `docs/api/sd.md` and `docs/api/analysis.md`.
  - **Report:** a dated report with aggregates under `results/`.
  - **Roadmap and repository guide:** roadmap item 16 and `AGENTS.md`.
  - **Version:** MOTCO 0.18.0.

**Non-goals.** Each is a separate change:
- applying the test to SEA-AD. `sead_sexdiff` does that in its own change, pinning this release, and only under `split_validated`;
- a group-aware or JIVE-style space refitted per permutation (thread 2);
- a sparse feature-wise interaction scan (thread 3);
- orthonormal span coordinates for the trajectory statistics (item 15);
- the split test in the simulation harness or study configs;
- geometry (`delta`/`angle`/`shape`) of an out-of-span difference;
- four-stage designs.

## Capabilities

### New Capabilities
- `interaction-span-split`: the split interaction test. It covers the interaction energy and pseudo-F per part, the exact decomposition, the Freedman–Lane permutation with the RRPP draw sequence and residual corrections, and the result's contents.
- `span-split-validation`: the validation of the split test. It covers the simulation and plasmode profiles, reproduction of their parent runs, the predeclared decision and the dated report.

### Modified Capabilities
- `real-data-trajectory-analysis`: `analyze_trajectories` and `motco analyze` gain the opt-in split test, its record block and its rejection with latent adjustment, byte-identical when absent.
- `module-contracts`: `motco.stats.interaction` joins the public import surface.
- `latent-space-metric-compatibility`: the visibility documentation states what the split test measures, how its parts are read, and that "in span" means the fitted span.

## Impact

- **Code:**
  - `src/motco/stats/interaction.py` (new) and the `motco.stats` namespace;
  - `src/motco/analysis.py` (parameter, validation, record block) and `src/motco/cli.py` (flag, printed summary);
  - `src/motco/__init__.py` (0.18.0).
- **Scripts:** `scripts/span_split_check.py` and `scripts/span_split_check_array.sbatch`. They import `latent_space_comparison` and `plasmode_check` for their designs, seeds and export contract, and change neither.
- **Tests:**
  - the decomposition and statistic definitions against direct computation;
  - the draw sequence against `RRPP`;
  - residual corrections, including F's invariance under `df`;
  - analysis byte identity without the setting, and unchanged outputs with it;
  - the analysis path against the stats function on the same data;
  - driver reproduction, resume and settings binding;
  - every decision path.
- **Compute:** about 6,600 simulation units and 1,500 plasmode units. The estimate is about 70 core-hours on `ing`: the simulation unit is dominated by the PLS CV fit (about 20 s in the comparison), and the plasmode unit is one analysis without attribution. A local smoke replaces the estimate before launch.
- **Data:** the plasmode profile reads the SEA-AD astrocyte `cohort_30` export in place under the plasmode's contract. The files' sha256 must equal the plasmode's committed settings, and no donor identifier or label vector is recorded.
- **Compatibility:** every committed configuration, fixture and report is unchanged. The analysis option is absent by default, and the stats module is new.
