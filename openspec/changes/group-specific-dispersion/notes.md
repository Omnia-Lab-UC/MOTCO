# Notes

## Exploratory probes (2026-11-24; not results)

Run on the workstation from `main` at `f9a1001` (MOTCO 0.25.0) before the proposal was written, with BLAS pinned to one thread. The scripts were scratch files outside the repository and are not committed. Their seeds are spent: no unit of this change may reuse them (the driver derives its seeds from base seed 2043 through `latent_space_comparison.generator_seed` and `rrpp_seed`).

The generator has no dispersion setting yet, so each probe edited a dataset after generation. `s` multiplies group B's noise standard deviation; the share is group A's (`group_ratio`), so at a share of 0.8 the wider group B holds a fifth of the samples.

### 1. Dispersion through the production fit

**Question.** Does a difference in group dispersion alone make the production test reject, and does it depend on balance and on the stage effect?

**Construction.**
- Null datasets from `generate_semisynthetic_trajectory` (mode `none`, n = 200, three stages, `p_dmp` 0.1, every per-omic δ at 2.0 or 0.5).
- Group B's rows moved to `mean + s·(row − mean)` about their exact population cell means (`population_trajectories`), methylation in M-values through the clipped logit and back. This is the construction the change proposes for the generator (design D1).
- Production: stage-PLS on methylation and expression with the study configs' cross-validation (3 inner and 4 outer folds, 5 repeats, at most 20 components, `random_state` 1203), rank floor 2, net displacement, 199 permutations. Rejection at p < 0.05.
- 100 replicates per cell, seeds from a hash of `["dispersion-probe", 990008, δ, share, n, …, replicate]`. The Monte Carlo SE at a rate of 0.05 is 0.022. Twelve workers on eight cores; a unit took 22–27 s (13–14 s with four workers).

**Rejection rates.**

| δ | share | s | `delta` | `angle` | `shape` | rank q25 / median / q75 |
|---|---|---|---|---|---|---|
| 2.0 | 0.5 | 1 | 0.03 | 0.01 | 0.02 | 2 / 2 / 2 |
| 2.0 | 0.5 | 2 | 0.17 | 0.02 | 0.01 | 3 / 3 / 3 |
| 2.0 | 0.8 | 2 | 0.44 | 0.16 | 0.10 | 2 / 3 / 3 |
| 2.0 | 0.2 | 2 | 0.01 | 0.00 | 0.00 | 3 / 3 / 3 |
| 0.5 | 0.5 | 1 | 0.05 | 0.04 | 0.03 | 3 / 4 / 4 |
| 0.5 | 0.5 | 2 | 0.65 | 0.16 | 0.07 | 8 / 9 / 11 |
| 0.5 | 0.8 | 2 | 0.71 | 0.38 | 0.10 | 5 / 5 / 6 |
| 0.5 | 0.2 | 2 | 0.10 | 0.01 | 0.00 | 10 / 12 / 18 |

Every rate is a count out of 100. The held-out stage AUROC was at least 0.992 in every unit.

**What it suggests.**
- With equal groups, `delta` exceeds at δ = 2.0 and rejects in about two thirds of null datasets at δ = 0.5. `angle` and `shape` hold at δ = 2.0.
- With the wider group the smaller, all three statistics exceed.
- With the wider group the larger, the test is conservative at δ = 2.0. At δ = 0.5 `delta` reads 0.10 while `angle` and `shape` stay near zero.
- The selected rank rises with the scale: from 2 to 3 at δ = 2.0, and at δ = 0.5 from 4 to 9 with equal groups and to 12 when the wider group is the larger.

### 2. Within-group draws

**Question.** Do RRPP draws that permute the reduced residuals within each group restore size?

**Construction.**
- The datasets of probe 1's construction on all three blocks (base 990009), the production fit forced at a rank (2 at equal dispersion, 3 under a scale at δ = 2.0, 4 at δ = 0.5) with no cross-validation.
- Both nulls on one fit, one set of observed statistics and one seed at 199 draws: the harness's unrestricted draws, and `RRPP(..., strata=<group labels>)` on the same latent matrix and design, with the harness's extraction and p-value helpers.
- 200 replicates per cell. The Monte Carlo SE at a rate of 0.05 is 0.015.

**Rejection rates** (`delta` / `angle` / `shape`):

| δ | share | s | rank | unrestricted | within group |
|---|---|---|---|---|---|
| 2.0 | 0.5 | 1 | 2 | 0.055 / 0.030 / 0.055 | 0.050 / 0.025 / 0.060 |
| 2.0 | 0.8 | 1 | 2 | 0.075 / 0.055 / 0.065 | 0.065 / 0.060 / 0.055 |
| 2.0 | 0.5 | 2 | 3 | 0.075 / 0.060 / 0.010 | 0.015 / 0.055 / 0.000 |
| 2.0 | 0.8 | 2 | 3 | 0.370 / 0.190 / 0.065 | 0.030 / 0.070 / 0.000 |
| 2.0 | 0.2 | 2 | 3 | 0.005 / 0.015 / 0.005 | 0.045 / 0.035 / 0.090 |
| 0.5 | 0.5 | 2 | 4 | 0.530 / 0.050 / 0.030 | 0.515 / 0.050 / 0.030 |
| 0.5 | 0.8 | 2 | 4 | 0.575 / 0.240 / 0.095 | 0.230 / 0.020 / 0.005 |

**What it suggests.**
- At δ = 2.0 the within-group draws remove the excess, and the conservatism where the wider group is the larger. `shape` reads 0.090 there, two standard errors above the level.
- At δ = 0.5 with equal groups they change nothing for `delta` (0.530 against 0.515), and with the wider group the smaller they leave 0.230.
- So part of the excess is the exchange of residuals between groups, and part is not. A guess, not tested: the pooled stage fit follows the wider group's noise more, so that group's in-sample path is longer in scores that a residual permutation does not refit. Roadmap item 19 found an in-sample fit behind another size excess.
- Equal-group `delta` at δ = 2.0 and s = 2 reads 0.075 here and 0.17 in probe 1. The two differ in blocks (three against two) and in how the rank is chosen (forced against cross-validated).

### 3. Correlation scaling

**Question.** The originally planned grid had a correlation scaling row beside the covariance one. What does shrinking group B's correlations do in the measured space?

**Construction.**
- Null datasets as in probe 1 at equal groups (seeds 770000 + replicate). Group B's rows redrawn from N(mean, c·Σ + (1 − c)·diag Σ) with a fresh generator (880000 + replicate), so its variances are kept and its correlations multiplied by c.
- The production fit forced at rank 2 (δ = 2.0) or 4 (δ = 0.5), no permutation. Read: the within-cell variance of the scores, group B over group A, summed over components and per component. 30 datasets per cell, on two and on three blocks.

**Within-cell variance of the scores, B over A** (median; two blocks / three blocks):

| δ | rank | c = 1, s = 1 | c = 1, s = 2 | c = 0.5 | c = 0 |
|---|---|---|---|---|---|
| 2.0 | 2 | 0.99 / 0.97 | 3.97 / 3.87 | 0.68 / 0.67 | 0.28 / 0.35 |
| 0.5 | 4 | 1.00 / 1.02 | 4.27 / 4.36 | 0.55 / 0.52 | 0.05 / 0.04 |

Per component at c = 0 and δ = 0.5 the ratio runs from 0.01 to 0.31.

**What it suggests.**
- A scale s gives a ratio near s² in every component, as it must.
- Less correlation in a group is less dispersion of its scores, by a factor the fit sets. In total dispersion it is a scale below one on that group, which the proposed cross reads as a scale above one on the other.
- With no correlation the ratio is far from uniform across components, which a scalar scale does not produce (design D3).

### What the probes are not

One sample size, 100 or 200 replicates at 199 permutations, two blocks in probe 1 and forced ranks in probes 2 and 3, a scale of 2 only, and edits made outside the generator. They motivated the setting, the scales, the columns and the recorded arm. The run redraws every cell on its own seeds at 500 × 999 under the predeclared rule.

## Implementation notes (2026-11-24)

- **The balanced designs are balanced to one sample per stage.** The parent's sizing gives each stage `round(n / 3)` samples and each group `round(share × stage)` (Python rounds half to even). At n = 200 the stages hold 67 / 67 / 66, so share 0.5 gives group A 34 / 34 / 33 and group B 33 / 33 / 33: the scaled group is two samples smaller in every `balanced_unequal` design at n = 200. The committed `n200` reference has the same table. The report states it beside the balanced readings; it is far from the 0.65 share of the nearest unbalanced design.
- **Task 2.5's draw check, amended before any unit.** The task asked for a response whose residuals are constant within group. The reduced model holds the group term, so such residuals are zero, and every null, restricted or not, would return the observed statistics. The test instead captures each permuted response the draws evaluate, maps it back to the reduced residual each row received, and checks that every draw is a permutation that never crosses a group, while the unrestricted draws do (`tests/test_group_dispersion_check.py::test_every_draw_keeps_each_residual_inside_its_group`). The arm's own definition is unchanged.
- **What `same_statistics` compares.** The refit at the recorded rank uses `integration_params["forced_components"]`, which skips the cross-validation and fits the same final model, so the observed `delta`, `angle` and `shape` are bit-identical to production's (checked on one unit of each profile in the tests). With one stratum the arm's path reproduces production's p-values exactly.

## Smoke (task 4.1, 2026-11-24)

Both profiles at 2 replicates per cell, every design, 999 permutations, run from the worktree (MOTCO 0.26.0, uncommitted, on `f9a1001`) on the workstation with BLAS pinned in the shell, 10 workers on 16 threads (load about 10, with the fast test suite running beside it). 128 units: 72 core (18 designs × 2 cells × 2) in 3 min 10 s and 56 columns (28 designs × 1 cell × 2) in 2 min 53 s of wall time.

- **No failure.** No unit failed, so no planted fixed point failed at any core design (the translation runs under the planted scale with the dispersion's within-cell multiplier). Every within-group refit had `same_statistics` true (128 of 128).
- **`merge`, the reproduction check and `decide` ran.** The reproduction compared 36 rows (2 cells × 18 parent arms) and differed on 28, only because of the smoke's size: 2 units against the committed 500. The 8 that match are the `concat` arms, which neither run records. `decide` read `reproduction_failed`, as it must here. The families do not hold at 2 replicates, where any rejection exceeds the bound. The smoke also exposed a bug: the map's column order came back alphabetical from the stored settings (sorted keys). It was fixed before any run, with a test (`test_the_map_order_survives_stored_settings`).

**Per-unit time by column** (seconds; the CV fit is most of it, and the within-group arm adds under a second):

| profile | column | units | designs | mean unit | max unit | mean fit | mean within-group | projected core-h at 500 reps |
|---|---|---|---|---|---|---|---|---|
| core | n200 | 72 | 18 | 24.6 | 36.6 | 21.1 | 0.9 | 123 |
| columns | n80 | 14 | 7 | 18.9 | 27.1 | 15.8 | 0.9 | 18 |
| columns | n500 | 14 | 7 | 37.5 | 50.9 | 33.1 | 0.9 | 37 |
| columns | n200_delta0.5 | 14 | 7 | 25.8 | 35.2 | 22.4 | 0.9 | 25 |
| columns | band_n320 | 14 | 7 | 31.1 | 42.8 | 27.6 | 0.8 | 30 |

**Projection.** 233 core-hours at the smoke's per-unit times (the parent's unit 225, the within-group arm 8). The stage-effect scaling run's arrays took about 1.4 times their unpacked projection on `ing`, which gives about 330 core-hours, inside the proposal's 250–350.

**Proposed cap: 400 core-hours** for both arrays together, about 1.2 times the inflated projection. Suggested arrays: `core` 200 shards (90 units each, about 1 h per task at the inflated time) and `columns` 140 shards (100 units each), throttled at `%64`. If the cap binds, the levers stay those of the design (Risks): first the n = 500 column, then the translation cells at s = 1.25, then the within-group arm. No amendment to the designs, cells, grade or rule is proposed.

**Approved 2026-11-24.** The user approved the 400 core-hour cap and the arrays as proposed (`core` 200 shards, `columns` 140, `%64`), with no amendment. The run starts from the pushed branch's revision in a clean worktree on `ing`.

Smoke records: `/tmp/claude-1000/-home-tomas-Documents-Repos-MOTCO/6a5069eb-91d3-47bb-937e-577127cb07df/scratchpad/gd-smoke/` (scratch, not committed).

## Run (task 4.2, submitted 2026-11-24)

- **Revision.** `5aab30e` on `feat/group-specific-dispersion` (pushed), MOTCO 0.26.0, in a clean worktree on `ing`: `/home1/tgonza/MOTCO-group-dispersion` (its `.venv` a symlink to `/home1/tgonza/MOTCO/.venv`; `PYTHONPATH` on the worktree's `src`). Python 3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3. The new and stage-span test files passed on a compute node first (155 passed).
- **Run root.** `/home1/tgonza/MOTCO-group-dispersion-run`, settings pre-created by `run --limit 0`: `core/settings.json` sha256 `9cc9a813…6741751`, `columns/settings.json` sha256 `098763cc…d79b9730`.
- **Arrays** (submitted 20:05 -03): `core` job **958479** (`--array=0-199%64`, N_SHARDS=200) and `columns` job **958544** (`--array=0-139%64`, N_SHARDS=140), partition `512x1024`, n11 excluded, 1 CPU and 4 GB per task.
- **First check** (20:17, 12 min in): 3,935 units completed (2,064 core, 1,637 columns), none failed, every within-group refit consistent, mean unit 17 s (below the smoke's 25 s), no task failed and no error in the logs. At that pace the 32,000 units take about 150 core-hours, well under the cap.
- **To finish** (from the `ing` worktree, after every task of both arrays completes; resubmit failed array ids first, complete units are skipped):

  ```bash
  cd /home1/tgonza/MOTCO-group-dispersion && export PYTHONPATH=$PWD/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  RUN=/home1/tgonza/MOTCO-group-dispersion-run
  sacct -j 958479,958544 -X -o JobID,State,ExitCode,Elapsed,TotalCPU,NodeList > $RUN/sacct.txt
  .venv/bin/python scripts/group_dispersion_check.py merge --profile core --out-dir $RUN
  .venv/bin/python scripts/group_dispersion_check.py merge --profile columns --out-dir $RUN
  .venv/bin/python scripts/group_dispersion_check.py decide --out-dir $RUN
  ```

  Then copy the aggregates (not the records) to `results/group-dispersion-<date>/` in a local worktree of the branch, with `PROVENANCE.txt` and `sacct.txt`, and write the report (task 4.3).
