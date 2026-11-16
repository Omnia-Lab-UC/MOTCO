# Tasks

## 1. The primitive: permutation strata in `RRPP`

- [x] 1.1 Add `strata` to `RRPP` and `_RRPPWorker` (design D2, D3): validation of the length and of missing labels beside the other input checks; the row groups in ascending label order; one shared draw helper used by the serial loop and the worker; the literal `rng.permutation(n)` call kept when `strata` is `None`. Document the argument, the draw sequence and the note on small strata in the docstring. Verify in `tests/test_rrpp_strata.py` that:
  - over many draws every index stays in its row's stratum and each draw is a permutation; a one-row stratum never moves;
  - with one stratum the three null distributions equal `RRPP`'s without strata at the same seed, bytewise;
  - reordering rows together with their labels still visits the strata in ascending label order;
  - with the stage indicators in the reduced model and no correction, every draw's stratum means equal the reduced fit's;
  - a wrong length and a missing label raise before any draw, naming the counts or the row.
- [x] 1.2 Pin identity while absent with a frozen copy of the 0.24.0 serial loop inside the test (the `_frozen_…` pattern, not a stored digest): on fixed inputs and seed, `RRPP(...)` and `RRPP(..., strata=None)` return null distributions whose `.tobytes()` equal the frozen loop's, with and without `return_eigengaps` and under both orientation estimators. Verify the test passes under `OMP_NUM_THREADS=1` and under default threads.
- [x] 1.3 Tie the primitive to the diagnosis's construction: on a three-stage design with unequal cells and a covariate, `RRPP(..., strata=<stage positions>)` equals `predictive_size_excess_check.restricted_null` on `within_stage_indices` at the same seed, without a correction and with `df`; and with `df` it equals the stratified call on `ŷ_red + κ·e`. Verify both in `tests/test_rrpp_strata.py`.
- [x] 1.4 Cover the parallel path: two stratified runs with one seed and worker count return identical distributions; on a response whose residuals are constant within each stratum every parallel draw returns the unpermuted statistics. Verify in `tests/test_rrpp_strata.py`, and that `tests/test_residual_correction.py` still passes unchanged.

## 2. The harness setting

- [x] 2.1 Add `PERMUTATION_STRATA` and `SimulationEvaluationParams.permutation_strata` (design D5): validation in `__post_init__`; the stage positions from `design.stage_levels` passed to `RRPP` as `strata`; `runtime_metadata["permutation_strata"]` only when set. Verify in `tests/test_permutation_strata.py` that:
  - with and without the setting, the latent matrix, the latent metadata and the observed statistics are identical, and the p-values and null summary under `stage` equal a direct `RRPP(..., strata=)` call on that latent matrix and seed;
  - it works with `latent=` prebuilt (no refit, patched call count), with `latent_coordinates="predictive_orthonormal"` and with `concat`;
  - `permutation_strata="group"` raises `SimulationEvaluationError` naming the allowed values before any fit;
  - with `permutations=0` the evaluation completes.
- [x] 2.2 Add `permutation_strata` to `grid._SIGNATURE_OMIT_WHEN_NONE`. Verify that two cells differing only in the setting have different ids and signatures, that only the set one's record carries the key, and that `tests/test_committed_config_identity.py`, `tests/test_post_floor_config_identity.py`, `tests/test_orientation_default_pins.py`, `tests/test_rank_floor_default_pins.py` and `tests/test_phase5_magnitude_profiles.py` pass unchanged.

## 3. The analysis setting and the CLI

- [x] 3.1 Add `TrajectoryAnalysisParams.permutation_strata` (design D6): validation and the two rejections in `__post_init__`; omission from `to_dict()` while unset; the working frame's stage positions passed to `RRPP` with the covariate-adjusted models and the residual correction; `measurement.permutation_strata` with the value and the stage sizes in declared order, only when set and permutations are drawn. Verify in `tests/test_analysis_permutation_strata.py` that:
  - only the null moves (rank, AUROC, statistics, scores, loadings and stage means identical; p-values those of the stratified `RRPP` call);
  - the analysis and the harness agree under `stage` on latent matrices, statistics, p-values, rank and null summary;
  - with covariates and `df`, the p-values equal `predictive_size_excess_check.infer(..., indices=within_stage_indices(...))` on the analysis's latent matrix, exactly;
  - with four stages declared in a non-sorted order, and with eleven stages, the strata follow the declared order;
  - the setting with `span_split_test` and with `interaction_scan` raises before any fitting, naming both settings; an unknown value raises naming the allowed values;
  - with `permutations=0` the record has no `measurement.permutation_strata`;
  - it combines with attribution, `covariate_latent_adjustment`, `block_weighting`, `sparsity` and each predictive `latent_coordinates` value.
- [x] 3.2 Confirm invisibility while absent: `tests/test_analysis_identity.py`, `tests/test_attribution_identity.py`, `tests/test_analysis.py` and the other `tests/test_analysis_*.py` pass unchanged, and a default `to_dict()` carries no strata key.
- [x] 3.3 Add `motco analyze --permutation-strata stage` and the summary line (design D7). Verify in `tests/test_analysis_cli.py` that the flag writes `measurement.permutation_strata` and the analysis's p-values, that the outputs without the flag equal the pre-change ones, and that the flag with `--interaction-scan` or `--span-split-test` exits with the message and writes nothing. Regenerate the `analyze --help` block of `docs/api/cli.md` and verify it matches the command's output.
- [x] 3.4 Bump `__version__` to 0.25.0 with its `CHANGELOG.md` entry (package changes only); update `AGENTS.md` (the `permutation.py` bullet, the harness options, the `motco analyze` options table and the `_SIGNATURE_OMIT_WHEN_NONE` list) and `docs/api/` for the new arguments. Verify with `tests/test_release_notes.py`, `uv run mkdocs build --strict` and the pre-commit gate (`ruff`, `mypy`, the fast tests).

## 4. Driver: settings, units and arms

- [x] 4.1 Create `scripts/within_stage_permutation_check.py` with `run`, `merge`, `decide` and `preflight` (design D13): a `settings.json` per profile bound by hash; round-robin shards; resume of complete unit keys; refusal of a resume or merge under other settings, naming the keys; the sha256 of the pilot's three `settings.json`, its `weak_size` and `plasmode` `rate_table.csv` and the diagnosis's `settings.json` and `rate_table.csv`, each refused by name if changed. Verify with tests for resume, a settings mismatch and a changed parent.
- [x] 4.2 Build the fitted profiles' units through the pilot's builders, imported unchanged (design D9, D10): `simulation_settings(base_seed=2038)` and `simulation_units`; `weak_size_settings` / `weak_size_units`; `plasmode_settings(export, replicates=1000)` restricted to the anchor, with cells `anchor` (replicates 0–499) and `anchor_fresh` (500–999). Verify that:
  - the counts are 7,800, 1,000 and 1,000;
  - no simulation dataset seed equals one of the pilot (2036), of the comparison (2035) or of the probe in `notes.md`;
  - the `weak_size` units and the plasmode units at replicates 0–499 carry the pilot's label, dataset and RRPP seeds;
  - no committed `settings.json` that reads the `paper` stream evaluates a replicate at or above 500;
  - a changed export file is refused before any unit.
- [x] 4.3 Implement the fitted unit (design D8): one production fit; the three predictive integrations from it; each of the four matrices evaluated with `permutation_strata` absent and `"stage"` at the unit's RRPP seed, through `evaluate_semisynthetic_trajectory` in simulation and `analyze_trajectories` at the forced rank on the plasmode; a unit that raises recorded as failed without aborting the shard. Verify that:
  - a simulation unit makes exactly one cross-validation and an anchor unit none (patched call counts);
  - on a unit of the pilot, the four unrestricted arm records equal the pilot's `evaluate_harness_unit` and, on the stand-in export, `evaluate_plasmode_unit` records (statistics and p-values);
  - a record holds seeds and statistics only.
- [x] 4.4 Implement the per-unit checks (spec: "Each unit checks that the two nulls read the same measurement"): the pilot's five checks on the unrestricted arms; `same_statistics` with tolerance 0; `construction_identity` on plasmode units against the diagnosis's loop. Verify with tests that a perturbed arm is recorded as inconsistent with the check's name and deviation and leaves every rate.
- [x] 4.5 Implement the `dispersion` profile (design D11): the SEA-AD cell table; the five variance patterns and two dimensions; the centroids with s solved for a stage share of 0.84 at d = 2; seeds from base seed 2038 by the driver's sha256 derivation; `matrix` and `matrix__within_stage` by `RRPP` at one seed, net displacement, no correction, 999 permutations. Verify that:
  - there are 10,000 units and no fit is made (patched call count);
  - the empirical stage share of variance at d = 2 under `equal` is 0.84 to Monte Carlo error on a large draw;
  - each pattern's realized stage variances match its multipliers on a large draw;
  - both arms of a unit share the matrix and the seed, and the observed statistics are equal;
  - the record holds seeds, statistics and p-values only.
- [x] 4.6 Implement `preflight` for the plasmode (design D10): the cross-validated production fit on replicates 0 and 500 against the forced call, written to `preflight.json` bound to the settings hash. Verify with tests for a pass, a wrong rank and a differing p-value.
- [x] 4.7 Add `scripts/within_stage_permutation_array.sbatch` (a profile variable and a dry-run variable that prints the command; absolute paths, BLAS pinned to one thread, a throttled array, `PYTHONPATH` on the worktree's `src`). Verify with `bash -n` and a dry run of its command lines against a temporary directory.

## 5. Merge, reproduction and decision

- [x] 5.1 Implement `merge` (spec: "Merge writes the rate, paired and cell tables"): the rate table; the paired table of each within-stage arm against its counterpart and of each predictive arm against production under the same draws, through `production_latent_space_check`'s exact tests imported unchanged; the cell table; `accounting.json`. Verify the row counts per design, cell and statistic, and that a pair with no discordant unit has three p-values of 1.
- [x] 5.2 Implement the reproduction checks, written to `reproduction.json`: `weak_size`'s four unrestricted arms against the pilot's committed table; plasmode cell `anchor` against the pilot's anchor rows (12) and the diagnosis's within-stage anchor rows (12); every differing row listed with both values. Verify with tests for both outcomes, and that the committed files' columns resolve every compared row.
- [x] 5.3 Implement the null's criteria (design D12, layer 1): size through `stage_span_paper_grade_decision.judge_rates` (imported unchanged) on the four families; protected power on the 54 target cells; specificity with the bound rule at δ = 2.0 and the bound-and-excess rule in the band; the gain p-values recorded. Verify with tests on synthetic tables for each criterion's pass and fail, that a band reading with a significant excess and a rate at the bound does not count, that an unavailable reading fails its criterion, and that the paired rules agree with the pilot's `judge_paired` on the pilot's committed tables.
- [x] 5.4 Implement the candidate's criteria (layer 2) for `predictive_orthonormal__within_stage` against `pls_cv__within_stage`, and the same grading, recorded only, for the other two predictive within-stage arms. Verify with tests for each of `size_excess`, `power_loss`, `specificity_excess`, `no_gain` and `qualifies`, and that a gain over `pls_cv` alone does not satisfy G.
- [x] 5.5 Implement `decide`, which writes `within_stage_permutation_decision.json` (rule version 1) with the verdict and the consequence statement of design D12. Verify with tests for each of the six verdicts, the gate with every criterion still recorded, that the candidate's reading cannot change a verdict when the null fails, and that the consequence statements are those of design D12.
- [x] 5.6 Implement the descriptive tables inside the decision (spec: "The decision reports descriptive comparisons beside the verdict"), the diagnosis's plant rows read from its committed table. Verify that changing any descriptive reading leaves the verdict unchanged, and add the driver row to `AGENTS.md`'s validation drivers table.

## 6. Runs and report

- [x] 6.1 Run a local smoke with BLAS pinned: the pre-flight on the real export; then two replicates of two simulation designs, four `weak_size` units, four plasmode units (two reproduced, two fresh) and twenty units per dispersion cell, through `run`, `merge` and `decide`. Record in `notes.md` the pre-flight outcome, per-arm and per-unit timing, the projected core-hours per profile and host, and the proposed cap. Launch of the full run waits on the user's approval of the cap and the host.
- [ ] 6.2 Run `plasmode` and `dispersion` where approved, then `simulation` and `weak_size` from a clean worktree at a pushed revision, and `merge` each. Verify 0 failed and 0 inconsistent units, both reproductions passed and the core-hours under the cap. Record the package versions and, on `ing`, the job ids and `sacct` output.
- [ ] 6.3 Run `decide`, then commit aggregates only under `results/within-stage-permutation-<date>/`: each profile's `settings.json`, rate, paired and cell tables, `accounting.json` and `reproduction.json`, the `preflight.json`, the decision and `PROVENANCE.txt`. Write `docs/reports/within-stage-permutation-<date>.md` with the verdict, every criterion's reading for the null and for the candidate, the reproductions, the descriptive comparisons, the `df` factor of design D4, the limitations and what follows. Verify that every number in the report can be read from the committed aggregates and that `uv run mkdocs build --strict` passes with the report in `mkdocs.yml`.
- [ ] 6.4 Add the permutation-null section to `docs/measurement-visibility.md` (spec `latent-space-metric-compatibility`): what each null assumes and holds fixed, both nulls' rates on the dispersion cells, the verdict and its consequence, and the candidate's reading. Update roadmap open items 1 and 2, the "Real-noise evidence rests on one export" entry and the Done list with the verdict, and the `--permutation-strata` help text and docstrings if the verdict is `size_excess` (design D7). Verify the roadmap's item numbers and report links resolve in the built site.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- Under `within_stage_validated` or `propose_paper_grade`, `sead_sexdiff` may re-issue the SEA-AD analysis under `--permutation-strata stage` in its own change, pinning 0.25.0.
- Under `propose_paper_grade`, a separate change runs `predictive_orthonormal` under within-stage draws at 500 × 999 with a rank rule for a fixed measured dimension.
- Whether within-stage draws become a default, and whether the split interaction test and the scan take them, are separate changes with their own readings.
