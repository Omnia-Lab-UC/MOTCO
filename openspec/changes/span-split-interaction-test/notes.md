# Notes

## Baseline (2026-11-05)

- Worktree `.claude/worktrees/span-split-interaction`, branch `feat/span-split-interaction-test`, from `main` at `41237e4`
  (PR #78).
- `src/motco/__init__.py` was `0.17.0`; bumped to `0.18.0` (task 5.3).
- Local runs need the zlib `LD_LIBRARY_PATH` workaround and BLAS pinned to one thread.

## Design notes from implementation

- **Observed quantities and the null share one computation.** Each part's interaction SS and RSS are quadratic forms of
  its reduced-residual Gram matrix against an orthonormal basis of col(F), ordered interaction directions first, then
  the reduced model's: SS = the first df_int forms, RSS = tr(G) − all of them. A draw permutes the weight rows by
  `argsort(idx)` instead of forming `G[idx][:, idx]`, so a draw costs O(n² · rank(F)) per part. The tests check the null
  against a brute-force refit of Ŷ_R + E[π] and against RRPP's captured outcomes at every draw (1e-10).
- **An empty basis is rejected** beside the spec's empty complement: a zero-dimensional span has no defined pseudo-F.
- **Zero energy.** `span_share` is absent only when the total interaction SS is exactly 0. A part with SS = RSS = 0 has
  an undefined pseudo-F (NaN, recorded as `null`) and no F p-value; its SS p-value is still reported.
- **The analysis at 0 permutations** passes no residual correction to the test (none is applied, as in RRPP), so the
  block's `residual_correction` is `null` there, and an undefined `leverage` correction cannot fail an analysis that
  draws nothing.
- **`draws_match_rrpp`** is `n_jobs == 1`: the split test always draws the serial sequence, while RRPP's parallel path
  spawns per-worker seeds.
- **Rejection rules follow the parents.** The artifacts say "α = 0.05" without the inequality. The simulation profile
  counts p < α (the comparison's `rate_table`), the plasmode profile p ≤ α (`plasmode_check._rejects`), for the split
  test and the production statistics alike, so every rate pairs with its parent's on the same units and the
  reproduction checks compare like with like. With 999 draws the two differ only at p = 0.050 exactly.
- **The plasmode's uncorrected SS diagnostic** reruns the split test without the correction on inputs rebuilt from the
  analysis (`analysis._working_frame`, `_covariate_adjusted_design`, `_span_split_inputs`). The unit checks that the
  rebuilt inputs reproduce the analysis block's observed sums exactly, else it fails.
- **Failed units are recorded in both profiles.** The plasmode's own driver raises on a failed consistency check; here a
  failed check (including the anchor latent check) is a failed record counted by cell, as the spec asks of the
  simulation profile, so the merge reports it instead of a silent gap.
- **Parents are checked before any unit.** The simulation settings refuse a calibration whose sha256 is not the
  comparison's `calibration_sha256`, a parent `settings.json` that fails its hash, and designs, cells, seeds or CV knobs
  that differ from the parent's. The plasmode settings refuse an export whose file hashes differ (naming the files), a
  fidelity result other than the parent's and a cell the plasmode did not run. The parents' `rate_table.csv` and
  `plasmode_summary.csv` sha256 are in the settings and re-checked at merge.

## Gate (task 5.3)

`ruff check src/ tests/` and `mypy src/motco/` clean. `mkdocs build`: the same 8 warnings as on `main` (none new).
`MOTCO_TEST_PERMS=99 pytest tests/ -m "not slow"`: 2377 passed, 1 skipped (R `shapes`), 32 deselected (9 min 50 s, BLAS
pinned).

## Smoke, timing and projection (task 6.1, 2026-11-05; local, 8-core laptop, BLAS pinned)

- **Fidelity gate** on the local export (`../sead_sexdiff/results/motco_input/cohort_30`): passed, rank 5, AUROC
  0.635556, 131 s. The export's five file hashes equal `results/plasmode-2026-11-02/settings.json`'s.
- **Simulation smoke** (`--reps 2 --control-reps 2 --perms 19 --workers 7`): 84/84 units, 0 failed, 3 min 26 s wall;
  median unit 15.7 s (fit 15.0 s, evaluation 0.29 s, split 0.011 s). The reproduction check ran and failed on all 126
  rows, each on count only (2 against the comparison's 100), as the smoke's size requires.
- **Plasmode smoke** (`--reps 2 --perms 19 --workers 6`): 6/6 units, 0 failed, every anchor through the latent check;
  median unit 138 s. The reproduction check ran and failed on its 9 rows by n (2 against 500).
- **`decide`** on the smoke: `reproduction_failed`, the gate working as specified.
- **Timing at 999 permutations.** Simulation: 240 units (anchor and translation × 40 × three designs, 7 workers), median
  16.5 s (fit 15.3 s, evaluation 0.68 s, split 0.14 s), q90 20.1 s. Plasmode: 3 units (3 workers), 133–134 s each.
- **Projection.** Simulation 6,600 × 16.5 s = 30 core-hours locally; plasmode 1,500 × 134 s = 56 core-hours locally.
  `ing`'s EPYC ran the comparison's stage-PLS fit at about 20 s against 15 s here, so × 1.3: about 40 + 73 = 113
  core-hours, budget 135 with the usual 20 % for node contention. This replaces the design's 70 (the plasmode unit is
  dominated by the CV fit on 80 × 10,000, not by attribution).

## Pre-check finding: the complement part is liberal when the CV rank exceeds K − 1 (2026-11-05, R = 40)

The 240 timing units are null units (anchor and translation share the null datasets; the translation adds no
interaction energy: the total interaction SS is within 4 % of the anchor's at every replicate, lower only by the
larger `_x_std`). Pseudo-F rejection (p < 0.05) by the unit's CV rank:

| design | rank 2 (= K − 1): span / complement | rank > 2: span / complement | median in-span share (rank 2 / > 2) |
|---|---|---|---|
| `n200` | 1/40 / 0/40 | 0/40 / 4/40 | 0.004 / 0.258 |
| `seaad` | 0/26 / 0/26 | 1/54 / 3/54 | 0.007 / 0.190 |
| `seaad_weak` | 3/45 / 2/45 | 1/35 / 5/35 | 0.255 / 0.264 |

Pooled, the complement rejects 12/129 = 0.093 null units at rank > K − 1 against 2/111 = 0.018 at rank K − 1; the
span and total parts stay at or below 0.067 everywhere. A component beyond K − 1 takes a large share of the null
interaction energy into the span (at `n200` replicate 2, rank 3 puts 632 of 1,107 there against 21 at rank 2), which
is the dependence design D3 names as the empirical risk. At R = 40 this was suggestive, not established; the user chose
a larger local check before any cluster run.

## Larger null check: no size excess (2026-11-05, R = 200 per design and cell)

Scratch driver (session scratchpad, not committed) over the profile's own units: the anchor and the translation, replicates
0–199 at each design, 199 permutations, 1,200 units, BLAS pinned, 14 workers. Each unit runs the split test with three
bases at its RRPP seed: `production` (as designed, `qr(x_rotations_)` at the CV rank), `pls_k1` (the first K − 1 PLS
rotations) and `stage_means` (the centred size-weighted pooled stage means, D3's exactly independent span). Pseudo-F
rejection (p < 0.05); the per-cell bound at R = 200 is 0.081:

| basis | cell | total | span | complement | complement at rank > 2 | median in-span share |
|---|---|---|---|---|---|---|
| `production` | anchor | 20/600 | 18/600 | 18/600 | 8/249 | 0.090 (0.339 at rank > 2) |
| `production` | translation | 18/600 | 16/600 | 21/600 | 17/426 | 0.220 (0.263 at rank > 2) |
| `pls_k1` | anchor | 20/600 | 18/600 | 26/600 | 16/249 | 0.015 |
| `pls_k1` | translation | 18/600 | 15/600 | 18/600 | 14/426 | 0.029 |
| `stage_means` | anchor | 20/600 | 17/600 | 20/600 | 12/249 | 0.013 |
| `stage_means` | translation | 18/600 | 17/600 | 18/600 | 14/426 | 0.019 |

- **No size excess under the basis as designed.** The largest design-level rate of any part is 0.045 (`seaad_weak`
  anchor complement); pooled rates are 0.027–0.035, slightly conservative. At rank > K − 1 the complement reads
  0.032 (anchor) and 0.040 (translation).
- **The pre-check's excess was a chance cluster.** The same units at replicates 0–39 reproduce it at 199 draws (anchor
  4/44 both times; translation 8/85 at 999 draws, 7/85 at 199); replicates 40–199 do not. Anchor and translation share
  their null datasets, so the pre-check's two cells were not independent evidence.
- **Reading, not validity: the null in-span share is large at rank > K − 1.** A fitted component beyond K − 1 carries a
  median 0.26–0.34 of the *null* interaction energy, against 0.01–0.03 for the K − 1 bases (0.08–0.30 at `seaad_weak`,
  where the span is fitted mostly to noise). An observed `span_share` must therefore be read against the share the null
  gives at that rank and design, never against 0. Task 6.5 documents this; the descriptive share-against-capture table
  reports it.
- **The SS diagnostic behaves as D1 predicts.** The uncorrected interaction SS rejects the null complement at 0.060
  (`seaad`, both cells) and 0.065 / 0.105 (`seaad_weak` anchor / translation), against 0.035 or less for the pseudo-F.
- **Decision:** the design stands; no basis revision is needed for size. Paper grade (R = 500, 999 draws) and the real
  noise of the plasmode (rank 5) remain to be measured by C1 and C2.

## Paper-grade run (tasks 6.2–6.5, 2026-11-05)

- **Launch.** Commit `b84d23c` pushed to the cluster clone's side ref `incoming/span-split` (not GitHub); clean detached
  worktree `/home1/tgonza/MOTCO-span-split`, `.venv` symlinked, `PYTHONPATH` on its `src`; driver, module and analysis
  sha256 equal to the workstation's; the export copy's five sha256 equal to the plasmode's settings. Fidelity job 934741
  passed (rank 5, AUROC 0.635556). Arrays 934742 (simulation, 200 shards) and 934743 (plasmode, 150 shards), each
  throttled to 64, n2 and n11 excluded at submission; n2 was re-allowed for pending tasks at 14:16 on the user's
  question (its `completing` state was another user's stuck job, the node healthy), and no task landed on it.
- **Incident.** The 51 simulation tasks first placed on n1 sat about 7 minutes in `D` on `rpc_wait_bit_killable`
  after opening their shard files, while n1 answered small NFS reads in 45 ms; they resumed on their own. Startup
  congestion with 128 tasks starting together; nothing was cancelled.
- **Outcome.** 350/350 tasks `COMPLETED 0:0`; 6,600 + 1,500 units, 0 failed, 0 duplicates; 25.6 + 36.9 = 62.5 recorded
  core-hours (median unit 11.1 s and 80.0 s on EPYC, faster than the laptop's 16.5 s and 134 s, so the projection of
  113 was high). Both reproduction checks passed (126 and 9 rows). Merge and decide: job 935094.
- **Decision: `split_validated`** (C1 0 of 18 above the bound, C2 0 of 3, C3 1.00 vs 0.08 and 0.68 vs 0.10).
  Report `docs/reports/span-split-2026-11-05.md`; aggregates `results/span-split-2026-11-05/`.
- **Descriptive notes the report carries.** The plasmode's complement is conservative on real noise (0.012); no part
  sees the 20-feature plant (0.02–0.05); the null in-span share grows with the fitted rank (0.006 at rank 2, 0.24 at
  rank 5), so the share is read against the null; the production `angle` at the `seaad` design reads 0.072 / 0.070
  at the anchor / translation (R = 500), just above its per-cell bound, which this change does not judge.
