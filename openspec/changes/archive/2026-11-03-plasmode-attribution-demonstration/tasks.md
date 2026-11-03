# Tasks

## 1. Baseline

- [x] 1.1 Confirm the starting state: `src/motco/__init__.py` reads 0.14.0, and the pre-commit gate (`uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`) passes on a clean `main`. Record the commit in the change's notes.
- [x] 1.2 Confirm the API surface needs no library change (design D7). Do this in a scratch script under the session scratchpad, not in the repo, on a small synthetic dataset. Verify that `analyze_trajectories` with covariates and attribution returns `attribution`, and that `latent_artifacts` holds the estimator, the features and the scales. Verify too that `analyze_orientation_attribution` exposes the per-feature principal-contrast values, `sign_stability` and `top_k_selection_frequency`.

## 2. Input contract and fidelity gate

- [x] 2.1 Create `scripts/plasmode_check.py` with the module docstring stating the input contract (design D2), the constants and `build_settings`. The constants are:
  - base seed 2034, the ladder, B = 200, top-k 20 and |K| = 10 per block;
  - the six covariates, `random_state` 1203 and `df`;
  - the gate's rank 5 and AUROC 0.635556;
  - the thresholds of design D9.

  Settings hash every constant and every input file's sha256. Verify that `uv run python scripts/plasmode_check.py --help` lists `fidelity`, `calibrate`, `run`, `merge` and `decide`.
- [x] 2.2 Implement the export loader. It reads the manifest, the blocks and the metadata, joins the covariates and derives `log10_n_rna` and `log10_n_atac`, and raises a named error for each contract violation in the spec. Verify with a new `tests/test_plasmode_check.py` on a synthetic stand-in export written to `tmp_path`, with the astrocyte cell table and 30 features per block:
  - the derived columns equal `np.log10` of the counts;
  - each violation (donor mismatch, unknown stage, three groups, missing covariate) raises its named error;
  - the recorded hashes equal `hashlib.sha256` of the files.
- [x] 2.3 Implement `fidelity`: one analysis on a permuted assignment with `permutations=0` and no attribution, whose rank and AUROC are compared with the gate's values and written to `fidelity.json`. Verify with tests that monkeypatch the analysis result: a match writes the file, and a mismatch exits non-zero naming both values and writes no calibration or record.

## 3. Null-making step and planted orientation

- [x] 3.1 Implement the seed derivation (base seed, stream, replicate, purpose) and the within-stage label permutation, with the real-labels guard and its redraw counter (design D3, D10). Verify in tests:
  - the cell table is preserved on the stand-in export;
  - the covariates and blocks stay on their donors' rows;
  - the same replicate gives identical labels at every strength;
  - the `paper` and `calibration` streams differ at the same replicate;
  - a monkeypatched draw equal to the real labels is redrawn and counted.
- [x] 3.2 Implement the plant (design D4). It draws K (10 per block, uniform without replacement) and the signs from the `plant` seed, computes `σ_f` as the pooled within-stage SD, and applies `e · c(stage) · σ_f · s_f` to the second group's donors. Verify in tests:
  - a middle-stage planted donor at e = 2.0 differs from the export by exactly `σ_f · s_f` on K and nowhere else;
  - e = 0 returns blocks equal to the export byte for byte (`pd.testing.assert_frame_equal` with `check_exact=True`);
  - K and the signs are identical across strengths of one replicate.

## 4. Unit evaluation and metrics

- [x] 4.1 Implement the analysis call and the consistency checks:
  - `analysis_params`: covariates, attribution at B/top-k with the `bootstrap` seed, the `rrpp` seed, `min_components=2`, net displacement and `residual_correction="df"`;
  - the projection check (1e-8 relative);
  - the anchor-unit latent check against `fidelity.json`.

  A failed check raises an error naming it and writes no record. Verify with tests that monkeypatch the analysis result to break each check.
- [x] 4.2 Implement the metrics (design D8). For every arm, component and block they are:
  - the planted-set precision over min(k, |K|);
  - the mean sign stability and top-k frequency over K;
  - the sign agreement among the recovered planted features.

  Add the descriptive unadjusted arm on `latent_artifacts`. Verify with tests on handcrafted contrast and bootstrap arrays whose metric values are computed by hand, including a K entirely outside the top-k (precision 0, and stability still defined).
- [x] 4.3 Implement `evaluate_unit` and the record line. The line holds:
  - the unit key and seeds, the redraw count, K's names and signs by block, and the realized cosine of the plant with the pooled displacement;
  - the statistics, p-values, rank, AUROC, eigengap flag and valid bootstrap counts;
  - the metrics.

  It holds no donor identifier and no label vector. Verify with one end-to-end unit on the stand-in export at 9 permutations and B = 5: the line parses as JSON, and its keys contain no donor id or label vector.

## 5. Records, resume and merge

- [x] 5.1 Implement `run`: sharded by `--index`/`--count`, resumable, with workers. It refuses a resume when the settings hash or the export hashes differ. Verify in tests: a second `run` on a complete shard appends nothing, and a changed threshold or a changed export file exits non-zero and appends nothing.
- [x] 5.2 Implement `merge`. It writes the per-unit metrics locally and `plasmode_summary.csv`, which gives the mean and SE per cell, arm, component, block and metric, over all units and over the units whose `angle` rejects, plus `merge.json`. Verify in tests: two merges give byte-identical files, and the anchor's planted-set precision appears as the false-recovery reference row.

## 6. Calibration

- [x] 6.1 Implement `calibrate` (design D5). It climbs the ladder with 40 units per step at 99 permutations, with attribution off, from the `calibration` stream. It adopts the first step whose `angle` rate is at least 0.90, falls back to the cap with `cap_reached`, and sets the descriptive cell to the step below or to `e*`/2. Every iteration is appended to `calibration.json`, and the adopted one is named in the settings. Verify in tests with a stubbed `run_once` that implements each calibration scenario of the spec: the first qualifying step, the cap reached and the lowest step qualifying.

## 7. Decision and documentation of the script

- [x] 7.1 Implement the pure `plasmode_decision` and `decide`. They read criteria (i), (ii) and (iii) (design D9), report the descriptive readings, produce the verdict `attribution_leg_lifted` or `attribution_leg_failed` with every failing criterion and cell, and mark the file indicative below paper grade. Verify with a test for each decision scenario in the spec: size holds and size fails on `shape`, detection fails at the cap, recovery holds and stability fails, lifted, failed on two criteria, smoke indicative, and descriptive readings that do not move the verdict.
- [x] 7.2 Add the script to `AGENTS.md`'s architecture, under `analysis.py`, beside the adjusted-attribution check. Describe its subcommands, construction, decision and verdicts, with no result yet. Verify the pre-commit gate passes, and that `uv run pytest tests/test_plasmode_check.py -q` passes.

## 8. Smoke, timing and the approval gate

- [x] 8.1 Run `fidelity` on the real export (`../sead_sexdiff/results/motco_input/cohort_30`) into the scratchpad. Verify it passes with rank 5 and AUROC 0.635556. If it fails, stop and report both values to the user.
- [x] 8.2 Run the smoke locally (R = 20 × 99 permutations × B = 50 at the anchor and at e = 1.0, with BLAS threads pinned to 1), then `merge` and `decide`. Verify that the decision file is marked indicative, that no consistency check failed, and that every anchor unit passed the latent check.
- [x] 8.3 Time one paper-grade unit (999 × B = 200). Extrapolate to calibration (up to 200 units at 99 permutations) and to the paper grade (1,500 units), and choose local workers or `ing` (design D11). Present the compute plan to the user and wait for approval before any calibration or paper-grade unit.

## 9. Calibration and paper-grade run

- [x] 9.1 Run `calibrate` on the real export on the approved machine. Verify that `calibration.json` records every iteration and the adopted `e*` (or `cap_reached`), and that the settings name it.
- [x] 9.2 If the run goes to `ing`, write `scripts/plasmode_array.sbatch` with absolute paths and the node exclusions, copy the export, and verify its sha256 against the settings before submitting. Run the paper grade (3 cells × 500 at 999 × B = 200). Verify that every shard completes (`sacct` `COMPLETED 0:0`, or every local worker exits 0) and that every unit has a record.
- [x] 9.3 Run `merge` and `decide`. Copy the committed files into `results/plasmode-<date>/` (design D11), and write `PROVENANCE.txt` naming the uncommitted records and metrics with their sha256, the MOTCO commit, the `sead_sexdiff` commit and the export hashes. Verify that the decision is not indicative, and that its settings hash equals `settings.json`'s.

## 10. Report and roadmap

- [x] 10.1 Write `docs/reports/plasmode-<date>.md` in the form of the adjusted-attribution report. It covers:
  - the construction and input contract, the null-making step and the plant, stating `e*` in SD units;
  - calibration;
  - each criterion with its value, bound, SE and R;
  - the verdict, read verbatim from `plasmode_decision.json`;
  - the descriptive readings;
  - the limitations of design's Risks: a stylized plant, conditional size, optimism on rejecting units.

  Every number must be read from `results/plasmode-<date>/`. Verify each number against the committed files.
- [x] 10.2 Update `docs/roadmap.md`:
  - item 12, done, with the verdict and the report link;
  - the Phase 6 Work bullet "For significant orientation …", whose withheld status is lifted or reported failed;
  - "Current position".

  Also update `AGENTS.md` with the run's verdict and results paths. Verify the pre-commit gate passes and `uv run mkdocs build --strict` succeeds.

  Done 2026-11-03: the gate passes. `mkdocs build --strict` aborts on 8 warnings that `main` has too (broken `results/` links in older Phase 5 pages); this change adds none (notes.md).

## Workflow follow-up

- Sync the delta spec into `openspec/specs/plasmode-attribution-demonstration/` and archive the change after review.
- Open the PR with the problem statement, a summary, the commands run and the verdict.
