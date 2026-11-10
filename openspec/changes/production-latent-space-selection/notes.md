# Notes

## Amendments before any unit (2026-11-10)

Decided by the user while implementing, folded into the proposal, design, spec and tasks:

- **S3 reads the target cells at e ≤ 0.5 only (m = 24).** The stage-span pilot (2026-11-04) measured production's
  own off-target leak at e = 1.0 under this null-scale construction (the full bend's `delta` 0.17 at n = 200,
  `angle` 0.10–0.12 at n = 80 and on the SEA-AD table), above the bound, so S3 at m = 32 would have failed every
  arm whatever it changes. The e = 1.0 cells stay in S2 and G; their off-target rates are descriptive
  (`descriptive.large_plant_off_target`).
- **S1 and S3 use the stage-span paper grade's count-plus-ceiling rule** (`judge_rates`, imported unchanged), not a
  per-reading bound: at R = 100 a reading at exactly α exceeds α + 2·SE with probability 0.028, so a per-reading
  bound fails an arm at α with probability about 0.45 at m = 21 and 0.50–0.60 at m = 24–32. S1 is two families
  (simulation m = 18, plasmode m = 3), both of which must hold.
- **The plant scenario was corrected.** The span plant adds e times the pooled centred stage configuration to the
  `Male` rows' own; under permuted labels that is (1 + e) × pooled only up to sampling, so the 1e-8 statement is
  about the shift, not the result.

## Smoke and compute projection (task 4.1, 2026-11-10)

Run on the workstation (AMD Ryzen 7 5800X, 8 cores / 16 threads) from the working tree at `d0825f0` plus the
uncommitted driver; Python 3.11.15, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, motco 0.22.0. BLAS pinned to
one thread. Outputs in the session scratchpad (`smoke-2026-11-10/`), not committed.

| Profile | Grade | Units | Workers | Wall | Contention |
|---|---|---|---|---|---|
| `simulation` | 2 reps × 199 perms, full CV | 60 | 7 | 52 min | the pre-commit gate ran beside it (pytest with unpinned BLAS, about 9 cores); workers ran at about 72 % of a CPU |
| `plasmode` | 2 reps × 999 perms, full CV, real export | 8 | 4 | 56 min | none |

**Mechanics.** 68 units × 4 arms, 0 failed records, 0 duplicates. `merge` wrote every table for both profiles.
The simulation reproduction compares 63 rows and every one differs only by the smoke's size (count 2 against the
parent's 100). The plasmode fidelity check passed on both anchor units (`pls_cv` rank 5, AUROC 0.635556). `decide`
read both profiles and returned `reproduction_failed`, as it must at this size, with S1–G recorded for every arm.

**Timing** (seconds per arm, record means):

| Profile | `pls_cv` | `weighted` | `sparse_block` | `weighted_sparse_block` | Unit mean |
|---|---|---|---|---|---|
| `simulation` (contended) | 22.1 | 23.0 | 152.3 | 151.6 | 349.3 |
| `plasmode` | 82.1 | 80.3 | 762.8 | 751.3 | 1,676.5 |

By design (simulation unit mean): n200 384, seaad 330, seaad_weak 335. The SEA-AD-sized designs cost nearly as much
as n = 200: the CV's cost follows the number of fits and the 498 features, not n. A sparse CV costs about 7× a plain
one in simulation and 9× on the plasmode's 10,000 features (the design assumed about 5×): every (rank, fraction)
candidate is refitted in every fold (45,360 power-method calls in one five-repeat CV at n = 200). The sparse arms
are about 85 % (simulation) and 90 % (plasmode) of a unit.

**Selected keep fractions** (every step is used; 0.5 is the most selected sparse step in both profiles):

| Profile | Arm | 1.0 | 0.5 | 0.2 | 0.1 | 0.05 |
|---|---|---|---|---|---|---|
| `simulation` | `sparse_block` | 32 | 18 | 4 | 1 | 5 |
| `simulation` | `weighted_sparse_block` | 42 | 9 | 2 | 6 | 1 |
| `plasmode` | `sparse_block` | 1 | 4 | 0 | 0 | 3 |
| `plasmode` | `weighted_sparse_block` | 0 | 4 | 1 | 0 | 3 |

**Convergence warnings.** The simulation smoke logged 58 `ConvergenceWarning`s from `SparsePLSDA`; the plasmode
logged none. Traced on one unit (n200 magnitude, replicate 1, five repeats): 17 of 45,360 power-method calls reach
`max_iter` (0.04 %), all in the weighted sparse arm at the smallest budget (18 + 7 features) on 100-sample inner
folds; the unweighted arm, the anchor, translation and the weak design had none. A capped call is a slowly drifting
2-cycle on a fixed support (iterates two steps apart differ by 5e-4, one step by 9.5e-3; the two points have cosine
0.9998): with the support and the discarded entry fixed the thresholded power step is linear and non-symmetric, and
its two leading eigenvalues can have near-equal modulus and opposite sign. sklearn's rule (squared change < 1e-6)
then never stops it; mixOmics has the same rule and caps at 100 iterations. The capped calls take 1.7 s of the 89 s
spent in the power method, so they neither drive the cost nor move a selection measurably. No change for this run.

**Weighting on the plasmode.** The export's blocks have 5,000 features each, so `equal_variance` gives them equal
weights: on real data the `weighted` arm is production without sklearn's per-fold re-standardization, not a block
rebalancing. That alone changes the CV rank at the anchor (5 → 7, AUROC 0.6356 → 0.6384). The simulation's blocks
(367 and 131 features) carry the weighting question; the report must say what the plasmode's `weighted` arm
measures.

**Anchor fits repeat.** Every plasmode anchor unit fits the same four label-free spaces (identical rank, fraction
and AUROC at both replicates); only the RRPP draws and labels differ. The anchor cell is 100 × 1,700 s ≈ 47
core-hours of identical fits. Reusing them needs a prebuilt-latent entry in `analyze_trajectories`, a package
change outside this run.

**Projection.** Uncontended workstation cost: simulation 349 s × 0.72 ≈ 250 s per unit × 3,000 ≈ 210 core-hours;
plasmode 1,677 s × 400 ≈ 186 core-hours. An EPYC 7662 core on `ing` is slower than the 5800X: the orthonormal
check's smoke measured 19 s for one n = 200 production CV there against about 15 s here, a factor of about 1.25.
**Projected on `ing`: about 260 (simulation) + 235 (plasmode) ≈ 490 core-hours**, about twice the proposal's
180–260.

**Levers** (not applied; the user decides):

| Lever | Saving | Cost to the check |
|---|---|---|
| Ladder 1.0 / 0.2 / 0.05 | about 40 % of every unit → ≈ 300 core-hours | drops 0.5 (35 of the smoke's 136 sparse selections, the most selected sparse step in both profiles) and 0.1 (7): it changes the arm, not just its price |
| Ladder 1.0 / 0.5 / 0.1 | about 40 % → ≈ 300 core-hours | keeps the most selected step; drops 0.2 and 0.05 (7 and 12 of 136; 0.05 is 6 of the plasmode's 16) |
| Plasmode plant cells at 50 replicates | ≈ 90 core-hours | G's three plasmode cells lose half their discordant units at a 0.05/11 threshold |

**Approved by the user (2026-11-10): run as specified (full ladder, 100 replicates), cap 650 core-hours; the plasmode's weighted arms stay as predeclared and the report states what they measure.** The recommendation was:
run as specified, cap 650 core-hours — about 1.3 × the
projection, leaving room for node packing (once 2.5× on n1, ADNI-design check 2026-11-02). If the recorded
core-hours reach 650 before the arrays drain, stop submitting, merge what exists and report.

**Proposed submission** (from a clean worktree at a pushed revision carrying this driver; after approval):
- `simulation`: 300 shards of 10 units (about 55 min each on EPYC), `--array=0-299%96`.
- `plasmode`: 100 shards of 4 units (about 2.3 h each), `--array=0-99%64`; the export copied to the cluster and its
  hashes checked by the driver.
- Wall time about 5–6 h at 96 + 64 concurrent tasks. Settings pre-created with `run --limit 0`; `logs/` exists in
  the worktree; out-dir outside the worktree; records stay on the cluster.
