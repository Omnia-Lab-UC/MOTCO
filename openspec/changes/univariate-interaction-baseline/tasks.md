# Tasks

## 1. Driver: profiles, units and seeds

- [x] 1.1 Create `scripts/univariate_baseline_check.py` with `run --profile parent|stage_span|surgery`, `merge --profile` and `decide` (design D9): `settings.json` bound by hash; round-robin shards; resume of complete unit keys; refusal of a resume or merge under other settings, naming the keys; the parents' sha256 recorded, with a refusal naming any changed file. Verify in `tests/test_univariate_baseline_check.py` the unit counts (4,400; 12,800; 16,500), resume, the settings-mismatch refusal and a refused parent whose hash differs.
- [x] 1.2 Implement the `parent` builder (design D5) over `interaction_scan_check`, `span_split_check` and `latent_space_comparison`, imported unchanged: the `n200` and `seaad` designs, fourteen cells, the anchor and translation at replicates 0–499 and the other cells at 0–99, 999 permutations, the scan validation's settings rebuilt and required equal to `results/interaction-scan-2026-11-09/simulation/settings.json`. Verify on sampled keys that the generator parameters and both seeds equal the scan validation's for the same design, cell and replicate, and that a changed `settings.json` or `comparator_table.csv` is refused by name.
- [x] 1.3 Implement the `stage_span` builder (design D5) at base seed 2043 with the comparison's `generator_seed` and `rrpp_seed`: the four designs (balanced n = 80, 200, 500 at δ = 2.0; n = 200 with every per-omic δ at 0.5), `layers` and `span_layers` absent, three stages, the thirteen cells with the planted scale at alignment 0 and the null scale at alignment 90, controls at 500 replicates and the other cells at 200. Verify that:
  - each design measures methylation, expression and proteomics, 658 features, with span layers equal to the measured blocks;
  - the designs differ from one another only in sample size and stage effect;
  - no dataset seed equals one of the comparison's, the split-test validation's, the scan validation's or the predictive-subspace pilot's;
  - cells of one design and replicate share their null dataset.
- [x] 1.4 Implement the `surgery` builder (design D5) over `stage_span_leak_check.column_cells` and `unit_params`, imported unchanged, from `design_guidance_core_magnitude_study.json` and `design_guidance_core_geometry_study.json`: the three columns with `layers` absent at n = 80, 200 and 500, the eleven cells, replicates 0–499, 999 permutations, the configs and both `design_point_operating.csv` tables bound by sha256. Verify on sampled keys that the generator parameters, evaluation parameters and replicate seeds equal the study engine's enumeration for the same column, mode, effect and replicate, and that the anchor is read from the geometry profile.

## 2. The unit and its arms

- [x] 2.1 Implement the unit (design D3): generate once, fit `pls` once with the cross-validation, evaluate the trajectory statistics with their null distributions, then the split test and the scan on `span_split_check.split_inputs` at the unit's RRPP seed, uncorrected; a failure recorded without aborting the shard. Verify on a small design that:
  - the `pls_cv` record equals a plain production evaluation of the same unit (rank, statistics, p-values);
  - exactly one cross-validated fit is made (patched call count);
  - the scan's largest F and global p-value equal a direct `interaction_feature_scan` call on the same inputs and seed, and equal the scan on `interaction_scan_blocks_check.simulation_inputs`, built with no latent fit;
  - a patched failure in one arm is recorded and the next unit runs.
- [x] 2.2 Implement the `production` and `production_minp` readings (design D2, D4): the smallest of the three p-values at 0.05 / 3, the three statistics singly at p < 0.05, and the smallest p-value calibrated on the unit's joint null draws. Verify with tests on hand-built p-values and null draws: the omnibus cases of the spec, `production_minp` equal to a direct rank computation, and its p-value never below `1 / (1 + permutations)`.
- [x] 2.3 Implement the `fdr` arm (design D1, D3) from the scan's per-feature F and degrees of freedom with `scipy.stats`: upper-tail F p-values, Benjamini–Hochberg over the non-degenerate features, rejection when the smallest adjusted p is at or below 0.05, and the number of discoveries. Verify against a hand-computed case of eight features, with a degenerate column left out of the family and counted, and with a case where the maximum does not reject and `fdr` does.
- [x] 2.4 Add the recorded arms (design D3, D8): the split test's `span` and `complement`, and `concat` and `pca_matched` evaluated at the same seed with the unit's production parameters changed only in the integration method (`concat`; `pca` at the unit's `pls_cv` rank), each singly and as the omnibus. Record hits and discoveries on population-null features through `interaction_scan_check.population_null_mask` where it is defined. Verify that `pca_matched`'s dimension equals the unit's `pls_cv` rank, that both arms measure the unit's blocks (three in `stage_span` and `surgery`, two in `parent`), that the record carries every field the spec names and no per-feature array, and that each arm's rejection uses the rule of design D4.
- [x] 2.5 Add `scripts/univariate_baseline_array.sbatch` (a profile variable and a dry-run variable that prints the command; absolute paths, BLAS pinned to one thread, a throttled array, `PYTHONPATH` on the worktree's `src`). Verify with `bash -n` and a dry run of its command lines against a temporary directory.

## 3. Merge and reproduction

- [x] 3.1 Implement `merge`: the rate table, the paired table (discordant counts over units available in both arms, `production_latent_space_check.paired_p_values` imported unchanged), the cell table (rank quantiles, AUROC, capture shares, median hits and discoveries, realized false discovery proportion) and `accounting.json`. Verify with tests on synthetic records for the columns, failed units counted by cell and excluded from rates, a paired row with no discordant unit, and duplicate handling.
- [x] 3.2 Implement the reproduction checks, written to `reproduction.json`, each listing every differing row with both counts:
  - `parent`: `scan`, `total`, `span`, `complement`, `delta`, `angle` and `shape` against `results/interaction-scan-2026-11-09/simulation/comparator_table.csv`, exact on count and rejections;
  - `surgery`: `pls_cv` `delta`, `angle` and `shape` against the two `design_point_operating.csv` tables, exact on available and rejected counts.

  Verify with tests for both outcomes of each check, and that the committed files' columns resolve every compared row.

## 4. Decision

- [x] 4.1 Implement `decide`, which writes `univariate_baseline_decision.json` (rule version 1, design D6):
  - the gate on any failed reproduction;
  - size per arm in two families through `stage_span_paper_grade_decision.judge_rates`, imported unchanged (m = 8 and 6), and each arm's admissibility;
  - the label of each of the 71 signal cells from the paired p-values, the levels 0.025 / 71 and 0.025 / (71·A) and the 0.10 margin;
  - the verdict in its order, with the consequence statement of each verdict.

  Verify with tests on synthetic tables for: each verdict; the gate still recording size, labels and grade; an inadmissible baseline arm leaving the labels and changing A; production beating two of three arms; a significant difference below the margin; an unavailable reading; and the levels at A = 1, 2 and 3.
- [x] 4.2 Implement the characterization grade (design D7): the in-span `stage_span` cells with a `production` rate of at least 0.80, their off-target statistics as one `judge_rates` family, and the three grades. Verify with tests on synthetic tables for each grade, the off-target pairs per mode, the exclusion of `surgery` and alignment-90 cells, and that changing the grade's inputs leaves the verdict unchanged.
- [x] 4.3 Implement the descriptive tables (design D8): the map by plant family and design, every arm's rate per cell, discoveries, the fit table and the `parent` table. Write the driver's module docstring (profiles, arms, seeds, parents, outputs, decision) and add the driver's row to the validation-drivers table of `AGENTS.md`. Verify with tests that each table has one row per cell and arm it names and that changing a recorded arm's rates leaves the verdict and grade unchanged, and with `uv run ruff check scripts/univariate_baseline_check.py tests/test_univariate_baseline_check.py` and the pre-commit gate: `uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.

## 5. Runs and report

- [x] 5.1 Run a local smoke of the three profiles at 2 replicates per cell with BLAS pinned in the shell (every design and column, at 999 permutations for a few units of each). Record per-unit timing by profile and arm, any failed planted fixed point by design, and each signal cell's rates. Confirm that `merge`, the reproduction checks and `decide` run, mismatching only because of the smoke's size. Verify that the timing, the core-hour projection, the proposed cap and any amendment are written into `notes.md` and put to the user for approval before any array is submitted.
- [ ] 5.2 Run the three profiles on `ing` through the sbatch array from a clean worktree at a pushed revision, then `merge` each. Verify that both reproductions pass and that failures occur only where the smoke found them. Record the job ids, `sacct` output and package versions.
- [ ] 5.3 Run `decide`, then commit aggregates only under `results/univariate-baseline-<date>/`: settings, rate, paired and cell tables, accounting, reproduction, the decision, `PROVENANCE.txt` and `sacct.txt`. Write `docs/reports/univariate-baseline-<date>.md` with the verdict, the grade, each arm's size reading, both reproductions, the map, the descriptive tables, the limits and what follows from the verdict. Verify that every number the report quotes can be read from the committed aggregates.
- [ ] 5.4 Add the default-analysis section to `docs/measurement-visibility.md` (spec `latent-space-metric-compatibility`): what each side tests, which side leads by plant family and design, what only the trajectory test can claim with the grade, the baseline's three combinations, and the limits. Update roadmap open item 13 and the matching "Not yet established" bullet with the verdict and what it permits, and add the report to the navigation in `mkdocs.yml`. Verify with `uv run mkdocs build --strict` and the pre-commit gate.

## Workflow follow-up

- Sync the delta specs and archive the change once the report is reviewed.
- Under `baseline_leads` or `complementary`, a separate change may add the false-discovery reading to the scan's output and to `motco analyze`, if the `fdr` arm was admissible and led where `scan` did not.
- Under any verdict, the real-noise reading of the same arms belongs to the first multi-omic ROSMAP or ADNI export (roadmap open item 2).
