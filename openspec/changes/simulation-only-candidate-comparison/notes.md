# Notes

## Implementation notes (2026-11-25)

- **Dependency.** `group-specific-dispersion` (0.26.0) merged to `main` as PR #110 (`5e5bcd0`); this branch was rebased onto it, so its proposal commit is now `a43b669`; the pre-rebase remote commit `173e18d` was merged with `-s ours` (`6e2c03e`) so the push is a fast-forward.
- **Generator (0.27.0).** `stage_dispersion_scale` and the per-cell `noise_scale`. Without a stage scale the pooled scales take 0.26.0's group term through the same function as before (`dispersion_within_multiplier` dispatches to `group_dispersion_within_multiplier`), so every group-dispersion dataset is bit-identical; the general term `Σ_c n_c·m_c² / Σ_c n_c` is `cell_dispersion_within_multiplier`. The committed-config identity, orientation, rank-floor and magnitude pins and the group-dispersion tests pass unchanged.
- **Design amendment before any unit: the four-stage designs are `n200_stages4` and `n200_delta0.5_stages4`.** The design named them `n200_four_stages` and `n200_delta0.5_four_stages`, but `n200_four_stages` is a study column of the committed orthonormal-coordinates run (`results/orthonormal-coordinates-2026-11-09/paper_grade/settings.json`). Its datasets differ (another base seed), but design D3 asks that no new design name occur in a committed settings file, so the two designs were renamed; `design.md` D3 and D6 carry the new names and the amendment.
- **Parents rebuilt.** Each parent driver's settings builder, called at its committed grade (`univariate_baseline_check.stage_span_settings` and `surgery_settings`, `stage_effect_scaling_check.saturated_settings` and `band_settings`, `adni_gain_specificity_check.simulation_settings`), rebuilds its committed `settings.json` exactly; the driver refuses to build a profile otherwise, naming the differing keys.
- **`adni_design` arms.** The gain run's unit (`adc.build_unit`, `ags.evaluate_arms`, `ags.common_checks`, `ags.check_construction`), with the gain settings' arms restricted to its two operating arms (`pls_cv`, `predictive_orthonormal`): their within-stage arms are not rerun, and the reproduction compares those two arms' 60 rows. `weighted` is the CV call with `block_weighting = "equal_variance"`. A failing check marks every arm of the unit `inconsistent` and leaves its rates, as the gain run did.
- **Records per (unit, arm).** A later `--arms sparse_block` appends to the profile's records; a resume reads its own shard file, so each arm selection keeps one shard count, and a selection running beside another takes its own. Each arm selection is recorded in `settings.json` outside the hash.
- **Checks before the smoke.** One unit of each kind at 19 permutations: the `pls_cv` p-values, rank and AUROC equal the parent drivers' units (`univariate_baseline_check.evaluate_unit`; `adni_gain_specificity_check.evaluate_simulation_unit`, both arms); the predictive arm makes no CV fit; `weighted` runs its own. Pre-commit gate: 3,645 passed, 3 skipped; the change's tests 139 passed.

## Smoke (task 5.1, 2026-11-25)

Every profile at 2 replicates per cell with every arm (`--arms all`) and 999 permutations, from the committed revision `4cc00ef` (clean tree) on the office workstation, BLAS pinned in the shell, 10 workers on 24 threads with other sessions' test suites running beside it (load 12–30). Records: `~/motco-runs/candidate-simulation-smoke/` (scratch, not committed).

- **No failure anywhere.** `committed` 152 units, `conditions` 130, `dispersion` 40 and `surgery` 44: 0 failed, 0 inconsistent, 0 missing records. No planted fixed point failed, at four stages either (`n200_stages4`, `n200_delta0.5_stages4`: every in-span and out-of-span cell), so no design is amended.
- **`merge` and the reproduction checks ran.** `committed` compared 228 rows and `surgery` 66, all differing only because of the smoke's size (2 units against the committed 200 or 500).
- **Per-arm time** (median seconds on this machine): generation 0.1–0.2; the production fit 7.3 (n = 80) to 20.0 (n = 500), shared by `pls_cv` and `predictive_orthonormal`, whose evaluations take 0.4 each; `weighted` about one production fit (7.6–20.4); `sparse_block` 72 at n = 200 and 113 at `band_n320`, about seven production fits.

**Projection** (the smoke's per-unit times weighted by each cell's full replicate count; `ing` up to 2.1× slower per core):

| tier | profiles and arms | units | core-h here | core-h on `ing` (×2.1) |
|---|---|---|---|---|
| 1 | `committed`, `conditions`, `dispersion` × `pls_cv`, `predictive_orthonormal`, `weighted` | 44,800 | 282 | 593 |
| 2 | `surgery` (above) and `adni_design` (its smoke still running at submission; the 19-permutation check took 81 s for the three calls under load) | 13,600 | 57 + ~60 | 120 + ~125 |
| 3 | `sparse_block` on `n200`, `n200_delta0.5`, `band_n320` | 9,200 | 215 | 452 |

About 615 core-hours here and at most about 1,290 on `ing`, inside the proposal's 900–1,300. Proposed caps: tier 1 700, tier 2 300, tier 3 550 core-hours.

**Approved 2026-11-25.** The user asked for everything to be left synced and running on `ing` before going home, and chose to launch all three tiers at once. The caps above are recorded as proposed; they were not reviewed separately.

## Run (tasks 5.2–5.4, submitted 2026-11-25)

- **Revision.** `4982bc4` on `feat/simulation-only-candidate-comparison` (pushed; the rebased branch was synced by merging the pre-rebase remote commit with `-s ours`, no history rewritten), MOTCO 0.27.0, in a clean detached worktree on `ing`: `/home1/tgonza/MOTCO-candidate-simulation` (`.venv` a symlink to `/home1/tgonza/MOTCO/.venv`; `PYTHONPATH` on the worktree's `src`). Python 3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3. The change's tests and the group-dispersion tests passed on a compute node first (170 passed at `6e2c03e`; `4982bc4` changes only the sbatch's usage text and its test).
- **Run root.** `/home1/tgonza/MOTCO-candidate-simulation-run`, settings pre-created by `run --limit 0` with each profile's arms (`committed` also with `sparse_block`): sha256 `committed` `9307a628…`, `conditions` `9ac8b4c2…`, `dispersion` `5492614f…`, `surgery` `57beb7e4…`, `adni_design` `09012ff1…`. Submissions are logged in `submissions.txt` there; `submit.sh` submits the six arrays two minutes apart.
- **A first submission was cancelled.** Arrays 959017 (`committed`) and 959082 (`conditions`) passed `CSC_ARMS` inside `sbatch --export`, which splits on commas, so the tasks ran `pls_cv` alone. Both were cancelled within minutes; their 49 `pls_cv` records are valid and kept (resume skips them; the resubmitted `committed` array uses the same 400 shards). The sbatch's usage text now passes `CSC_ARMS` through the environment, and its test checks the note. `committed/settings.json` records the selection `["pls_cv"]` among its arm selections.
- **Arrays** (partition `512x1024`, n11 excluded, 1 CPU per task, 4 GB, 8 GB for `adni_design`): tier 1 `committed` **959131** (400 shards, `%64`), `conditions` **959196** (320, `%48`), `dispersion` (200, `%32`); tier 2 `surgery` (220, `%32`), `adni_design` (130, `%16`); tier 3 `committed` × `sparse_block` (300 shards, `%48`; its own shard count, so no task writes a tier-1 file). The job ids of the last four are in `submissions.txt`.
- **To finish** (from the `ing` worktree, after every task completes; resubmit failed array ids first, complete (unit, arm) keys are skipped):

  ```bash
  cd /home1/tgonza/MOTCO-candidate-simulation && export PYTHONPATH=$PWD/src OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
  RUN=/home1/tgonza/MOTCO-candidate-simulation-run
  sacct -j <ids from submissions.txt> -X -o JobID,State,ExitCode,Elapsed,TotalCPU,NodeList > $RUN/sacct.txt
  for P in committed conditions dispersion surgery adni_design; do
    .venv/bin/python scripts/candidate_simulation_check.py merge --profile $P --out-dir $RUN
  done
  .venv/bin/python scripts/candidate_simulation_check.py decide --out-dir $RUN
  ```

  Then copy the aggregates (not the records) to `results/candidate-simulation-<date>/` with `PROVENANCE.txt` and `sacct.txt`, and write the report (task 5.5).
