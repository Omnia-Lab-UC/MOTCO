# Proposal

## Why

RRPP permutes the reduced model's residuals across all samples, which assumes they are exchangeable across stages. The predictive size-excess diagnosis (roadmap item 19, `both_restore`) found that this assumption is what fails in the predictive subspace on the SEA-AD plasmode anchor. Restricting the draws within stage removed the whole `delta` excess (0.130–0.136 → 0.056–0.062), and under the same draws production held size and gained plant power (`delta` 0.51 → 0.63 at e = 0.5). Its verdict opens this change: within-stage draws as an opt-in option, with their own size and power readings before `predictive_orthonormal` is re-read.

Three things are missing.
- **The draws are a driver's loop, not a package option.** `scripts/predictive_size_excess_check.py` carries them; `RRPP`, the harness and `motco analyze` cannot request them.
- **Nothing has read them where the truth is known.** The plasmode permutes its labels within stage, so the restricted draws match its null by construction. No simulated design has been read under them, and `shape` held the anchor's rule by one or two counts of 500.
- **The mechanism was described, not tested.** The diagnosis saw unequal residual dispersion by stage, in opposite directions in the scores and in the predictive coordinates.

An exploratory probe on 2026-11-16 (`notes.md`; not a result) fed Gaussian matrices with those dispersion patterns to both nulls on the SEA-AD cell table, 2,000 replicates at 199 draws:
- **Equal dispersion:** both nulls read 0.04–0.06 on every statistic.
- **The largest stage the narrowest** (the predictive coordinates' pattern): the unrestricted draws read `delta` 0.066–0.070 and `shape` 0.057–0.074, and 0.092–0.097 with the inequality roughly doubled. The within-stage draws read 0.039–0.054.
- **The largest stage the widest** (the scores' pattern): the unrestricted draws are conservative, 0.034–0.043, and `delta` and `shape` fall to 0.014–0.019 with the inequality doubled. The within-stage draws read 0.038–0.050.

So unequal dispersion by stage mis-sizes the unrestricted null in a direction the design sets, and the within-stage null holds in every case. The conservative direction would explain production's gain under within-stage draws. The measured inequality gives the direction of the predictive excess but only a fifth to a quarter of its size, so the anchor still has something the toy lacks. Unequal dispersion across disease stages is plausible in any real cohort, so the option concerns the production measurement and not only the predictive coordinates.

## What Changes

- **Opt-in permutation strata in `RRPP`.** `RRPP(..., strata=)` takes one label per sample. Each draw permutes the reduced model's residuals among the samples of each stratum, visiting the strata in ascending label order from the one seeded generator. One stratum reproduces today's draws exactly, and so does `None`, the default. The reduced fit, the residual correction, the orientation estimator and the full-model estimation are untouched, on the serial and the parallel path. With the stage positions as labels the draws are the diagnosis's, so its within-stage arms are reproducible through the package.
- **The same setting in the harness.** `SimulationEvaluationParams.permutation_strata` (`"stage"` or absent) passes each sample's stage position to `RRPP`. It is validated before any fit, hash-omitted while absent (`grid._SIGNATURE_OMIT_WHEN_NONE`) and recorded in the runtime metadata when set. Only the null moves: the fit, the rank, the coordinates and the observed statistics are unchanged, so it combines with every integration method and every `latent_coordinates` value.
- **The same setting in `motco.analysis` and `motco analyze`.** `TrajectoryAnalysisParams.permutation_strata` and `--permutation-strata stage`, with the covariate-adjusted models and the residual correction. It is omitted from `to_dict()` while unset, and `measurement.permutation_strata` records the value and the stage sizes when set. It is rejected with the split interaction test and with the interaction scan, whose own nulls were validated under unrestricted draws.
- **A paired pilot.** `scripts/within_stage_permutation_check.py` (`run --profile …` | `merge --profile …` | `decide`; sbatch array beside it) fits production once per unit and reads eight arms on the same dataset and RRPP seed: the scores (`pls_cv`) and the three predictive coordinate systems, each under the unrestricted and the within-stage draws, through the package option.
- **Four profiles.**
  - **`simulation`:** a fresh seed family (base seed 2038) on the predictive-subspace pilot's seven stage-span designs, eleven cells and surgery cell, at 100 × 199. 7,800 units. The pilot's own datasets selected `predictive_orthonormal`, so they are spent for judging it.
  - **`weak_size`:** the pilot's 1,000 weak-design control units (500 × 999). The four unrestricted arms must reproduce the pilot's committed counts exactly.
  - **`plasmode`:** the SEA-AD anchor on 1,000 within-stage label draws at 999 permutations, at the committed rank of 5 without cross-validation. On replicates 0–499 the unrestricted arms must reproduce the pilot's committed counts and the within-stage arms the diagnosis's, exactly. Replicates 500–999 are label draws no run has read: the fresh real-noise size reading, `shape` included.
  - **`dispersion`:** Gaussian latent matrices on the SEA-AD cell table with stage-unequal dispersion, where the truth is known and the two nulls are expected to differ. Five patterns (equal; the diagnosis's medians for `predictive_orthonormal` and for the scores; each with its log ratios doubled) at 2 and 5 dimensions, 1,000 × 999, both nulls. 10,000 units, with no fit.
- **A predeclared decision in two layers** (`within_stage_permutation_decision.json`, rule version 1).
  - **The null, on production's scores** (`pls_cv__within_stage` against `pls_cv`): size in four families under the stage-span paper grade's count-plus-ceiling rule (`judge_rates`, imported unchanged): the simulation's anchors and translations (m = 42), `weak_size` (m = 6), the fresh plasmode anchor (m = 3) and the dispersion cells (m = 30). Then protected power on the 54 target cells of the six designs at δ ≥ 0.25, and specificity on their off-target readings. Gains are recorded and not required.
  - **The candidate, under the same draws** (`predictive_orthonormal__within_stage` against `pls_cv__within_stage`), read only when the null is validated: the pilot's criteria S1–S3 and G with the reference changed to the same-draw production arm, so a gain is the coordinates' and not the null's.
  - **Verdicts:** `reproduction_failed` gates everything. The null's first failure gives `size_excess`, `power_loss` or `specificity_excess`. Otherwise `propose_paper_grade` when the candidate qualifies, else `within_stage_validated`, naming the candidate's first failure. Descriptive beside it: the other two predictive systems under the same criteria, the unrestricted arms on the dispersion cells and on the fresh anchor, the pooled 1,000-draw anchor, null widths, and the diagnosis's plant rows.
- **Documentation and records.** A dated report `docs/reports/within-stage-permutation-<date>.md`; aggregates under `results/within-stage-permutation-<date>/`; a section of `docs/measurement-visibility.md` on what each null assumes; `docs/api/`; roadmap open items 1 and 2; `AGENTS.md`; `CHANGELOG.md`; MOTCO 0.25.0.

**Non-goals.**
- **Any default.** The unrestricted draws stay the default of `RRPP`, the harness and `motco analyze`. Whether the within-stage null becomes the case study's or the default is a later change, with its fixture re-issues.
- **The SEA-AD re-issue.** Applying the option to the case study belongs to `sead_sexdiff`, after a verdict.
- **Within-stage draws for the split interaction test and the scan.** Both draw their own nulls and were validated under unrestricted draws. Extending them needs their own readings.
- **A strata-aware residual correction.** `df` and `leverage` keep their definitions under strata (design D4).
- **Strata other than the stage,** and `motco de` or the study configs as users of the option.
- **The paper grade of `predictive_orthonormal`,** its rank rule for a fixed measured dimension, and the ADNI and ROSMAP anchors of roadmap open item 2.
- **Why the anchor's excess is larger than unequal dispersion explains.** The dispersion profile sizes the mechanism; the remainder stays with open item 2.

## Capabilities

### New Capabilities
- `rrpp-stratified-permutation`: the opt-in exchangeability strata of the RRPP primitive. It covers the stratified draws and their defined sequence, identity while absent and with one stratum, what the strata leave unchanged, the parallel path and input validation.
- `within-stage-permutation-validation`: the paired pilot of the within-stage null against the unrestricted one. It covers the eight arms on one fit, the four profiles with their reproductions of the committed parents, the unit record and its consistency checks, the two-layer predeclared decision with its descriptive tables, and the dated report.

### Modified Capabilities
- `simulation-evaluation-harness`: the evaluation parameters gain the opt-in `permutation_strata` setting with its validation, its record and its invisibility while absent.
- `real-data-trajectory-analysis`: the analysis parameters and `motco analyze` gain the opt-in `permutation_strata` setting with its validation, rejections, record and invisibility while absent.
- `latent-space-metric-compatibility`: the documentation states what each permutation null assumes, what the within-stage draws hold fixed, and the pilot's verdict.

## Impact

- **Code:**
  - `src/motco/stats/permutation.py`: `strata` in `RRPP` and `_RRPPWorker`, with one shared draw helper and the validation;
  - `src/motco/simulations/evaluation.py`: the typed field, its validation, the stage positions passed to `RRPP`, the runtime metadata;
  - `src/motco/simulations/grid.py`: `permutation_strata` in `_SIGNATURE_OMIT_WHEN_NONE`;
  - `src/motco/analysis.py`: the typed field, its validation and rejections, the record;
  - `src/motco/cli.py`: `--permutation-strata` and the summary line;
  - `src/motco/__init__.py` and `CHANGELOG.md`: 0.25.0.

  `motco.stats.interaction`, the generator and the study engine are not touched.
- **Scripts:** `scripts/within_stage_permutation_check.py` and `scripts/within_stage_permutation_array.sbatch`. They import unchanged `predictive_subspace_check` (settings, units, arm parameters, consistency checks, decision helpers), `predictive_size_excess_check` (its within-stage loop, as the independent reference), `plasmode_check` and `production_latent_space_check` (the export contract, labels, exact paired tests) and `stage_span_paper_grade_decision` (`judge_rates`).
- **Tests:** the stratified draws (strata respected, one stratum equal to the unrestricted draws, equality with the diagnosis's loop, both corrections, the parallel path, validation); byte identity while absent in the primitive, the harness and the analysis; harness–analysis equivalence under the setting; the rejections; the CLI flag; driver units, seeds, parents, resume, merge, every reproduction outcome and every decision path.
- **Compute:** about 9,800 fitted units and 10,000 unfitted ones. The estimate before the smoke is about 57 core-hours on `ing` (simulation about 46, `weak_size` about 8, plasmode and dispersion under 2 each). A local smoke replaces the estimate, and launch waits on the user's approval with a cap.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged. The setting is absent by default everywhere, and the parents' committed files are read, hashed and refused by name if changed.
