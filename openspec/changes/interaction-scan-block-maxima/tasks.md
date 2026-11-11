# Tasks

## 1. The block-wise combination (`motco.stats.interaction`)

- [x] 1.1 Add the `blocks` keyword to `interaction_feature_scan` (design D1, D2): labels validated before any draw (one per column, an error naming both lengths), K counted over blocks with a non-degenerate feature, the (block, F descending, column) order, per-draw segment suffix maxima and per-block null maxima, the adjusted p-values min(1, K · p⁽ᵏ⁾) and the global p-value min(1, K · min_k p_k). Leave the unlabelled path untouched. Verify in `tests/test_interaction_scan.py` that every by-block adjusted p-value equals a brute-force per-block step-down on the returned null F, that the global p-value is the smallest adjusted p-value, and that a wrong label length raises before any draw.
- [x] 1.2 Extend `InteractionScanResult` (design D3): `combination`, block names, `n_blocks_tested`, per-block maximum F, column and global p-value, and the B × K null block maxima; `hits()` reads the by-block p-values; serialization adds the block fields only with labels. Verify that one label on every column gives the joint scan's global p-value, adjusted p-values and hits exactly, that a wholly degenerate block reports no maximum and leaves K at the other blocks, that every per-block null maximum is attained by a column of its block, and that the existing scan and split-test tests pass unchanged.

## 2. The analysis and the CLI

- [x] 2.1 Factor the scan's input construction out of `analyze_trajectories` into a private helper (design D4): unweighted pooled-standardized features, covariate-adjusted full and reduced models, LS-mean rows, contrast and stage levels, built without the latent fit; `analyze_trajectories` takes the scan's inputs from it. Verify that `tests/test_analysis_identity.py` and `tests/test_analysis_interaction_scan.py` pass unchanged, and that a scan on the helper's inputs equals the analysis's scan bit for bit, with covariates and with block weighting.
- [x] 2.2 Add `TrajectoryAnalysisParams.interaction_scan_by_block` (default `False`): rejected in `__post_init__` without `interaction_scan`, naming both; omitted from `to_dict()` while false; labels from the declared blocks; the record's `combination`, `n_blocks_tested` and per-block `max_f`, `max_feature` and `p_value`. Verify in a new `tests/test_analysis_scan_by_block.py` that the analysis's by-block results equal the primitive's with the declared blocks as labels, that every non-scan output and the scan's F and group differences are identical between the joint and by-block runs, that the table's adjusted p-values are the by-block ones, and that the identity analyses carry no by-block key.
- [x] 2.3 Add `motco analyze --interaction-scan-by-block`: a clean non-zero exit without `--interaction-scan`, naming both flags, and the combination and K in the printed scan line. Verify with CLI tests for the round trip and the refusal.
- [x] 2.4 Bump the version to 0.23.0 and run the pre-commit gate (`ruff`, `mypy`, fast tests).

## 3. The validation driver

- [x] 3.1 Write `scripts/interaction_scan_blocks_check.py` settings and units (design D5): `simulation` and `plasmode` settings rebuilt through `interaction_scan_check` and required equal to the committed ones (differing keys named); the export hash check; the `plasmode_unequal` subset (seed by sha256 of the driver name, 2034 and `unequal_subset`; 500 RNA features in export order; the feature list's sha256); the sha256 of every parent file, refused by name when changed. Verify in `tests/test_interaction_scan_blocks_check.py` that settings rebuild and bind, that a changed parent file is refused by name, and that the subset is fixed and recorded.
- [x] 3.2 Implement the scan-only units: simulation (`lsc.unit_dataset`, harness preprocessing, design models, both arms at the RRPP seed, population-null hits per block), plasmode and unequal (`plasmode_check.build_unit`, the helper of 2.1, both arms with `df`, planted counts per block per arm), plus per unit the joint null-maximum shares, the top-20 planted count, the block-maxima correlation and the minP global p-value; counts only, failures recorded. Verify on synthetic fixtures that both arms share every F, that the simulation inputs give the same F as the PLS input coordinates to 1e-10, that no record holds a feature or donor name, and that a raising unit is recorded as failed.
- [x] 3.3 Implement `run` (sharded, resumable, settings-bound, one job per unit), `merge` (the tables, accounting and `reproduction.json` against the parent's committed tables, exact) and the sbatch array `scripts/interaction_scan_blocks_check_array.sbatch` with a dry-run switch. Verify resume, the settings-mismatch refusal, and that reproduction passes on records rebuilt from the parent's tables and fails naming the cell on a perturbed count.
- [x] 3.4 Implement `decide` (design D6): `reproduction_failed`, S1–S3 through `stage_span_paper_grade_decision.judge_rates`, P's exact paired loss p-values, G's exact sign test, the verdict precedence and its consequence, and the descriptive tables. Verify every verdict path, the precedence, and that changing a descriptive rate leaves the verdict unchanged.

## 4. Smoke and run

- [x] 4.1 Run a local smoke of each profile (a few units per cell, 99 permutations) and write the timing projection to `notes.md`, replacing the proposal's estimate. Verify the projection is recorded with the revision it was timed at.
- [x] 4.2 Run the `simulation` profile on the workstation and the plasmode profiles on `ing` after the `production-latent-space-selection` arrays drain, from a clean worktree at a pushed revision carrying 0.23.0; merge each. Verify 0 unexplained failures and both reproductions passed, and record job ids, `sacct` output and package versions.
- [x] 4.3 Run `decide` and commit aggregates only under `results/interaction-scan-blocks-<date>/`: settings, the merged tables, accounting, reproduction, the decision, `PROVENANCE.txt` and `sacct.txt`. Verify that no per-unit record is committed.

## 5. Report and documentation

- [ ] 5.1 Write `docs/reports/interaction-scan-blocks-<date>.md`: the verdict, every criterion, both reproductions, per-block balance and recall in both arms, the minP comparison and what follows from the verdict. Verify that every number it quotes can be read from the committed aggregates.
- [ ] 5.2 Add the by-block reading to the scan section of `docs/measurement-visibility.md` (design D7), the labels and result fields to `docs/api/sd.md`, the setting, record fields and flag to `docs/api/analysis.md`, the verdict to roadmap item 16's thread 3 follow-ups, and the setting and driver to `AGENTS.md`. Verify with `uv run mkdocs build` and the pre-commit gate.

## Workflow follow-up

- Rebase onto `main` once PR #86 merges.
- Sync the delta specs and archive the change after the report lands.
