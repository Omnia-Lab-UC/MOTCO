# Notes

## Smoke (task 3.1, 2026-11-25)

Run on the workstation (`nixos-office`) in the worktree of `feat/unbalance-equal-dispersion-diagnosis`, on `4f0bdee` with the driver uncommitted: Python 3.11.15, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, MOTCO 0.26.0. Command: `run --reps 2 --out-dir <scratch> --workers 8`, BLAS pinned to one thread in the shell, 999 permutations, every design. The machine was shared: the load average was 48 on 24 cores throughout (another session's test suite), so these times overstate an idle node's.

- **16 units in 52 s wall** (8 workers). None failed, n = 42 included. Both refits were consistent in every unit: the within-group arm and the corrected arms read `same_statistics` true, and no correction failed. The CV rank was 2 in 15 units and 3 in one (`n80_r0.8`, replicate 1).
- **`merge`, the reproduction check and `decide` ran.** The reproduction compared 22 rows (21 rate rows of `n80_r0.8` and its within-group table row) and differed on 18, only because of the smoke's size: 2 units against the committed 500. The 4 that match are the `concat` arms, which neither run records. No row was unresolved. `decide` read `reproduction_failed`, as it must here. Every family was recorded and failed, since every reading is unavailable below 900 units. The re-read was not performed and `unbalance_excess` stands.

**Per-unit time by design** (seconds; the CV fit is most of it):

| design | larger / smaller per stage | mean unit | max unit | mean fit | mean within-group | mean corrected | replicates | projected core-h |
|---|---|---|---|---|---|---|---|---|
| `n80_r0.8` (rerun) | 22 / 5 | 28.3 | 28.7 | 22.7 | 1.0 | 1.6 | 500 | 3.9 |
| `n80_r0.8_b` | 22 / 5 | 26.6 | 27.0 | 21.5 | 1.1 | 1.6 | 1,000 | 7.4 |
| `n200_r0.925` | 62 / 5 | 24.2 | 24.4 | 20.3 | 0.6 | 1.3 | 1,000 | 6.7 |
| `n42_r0.5` | 7 / 7 | 24.2 | 25.4 | 19.2 | 1.0 | 1.4 | 1,000 | 6.7 |
| `n80_r0.74` | 20 / 7 | 18.2 | 19.2 | 14.8 | 0.5 | 1.1 | 1,000 | 5.1 |
| `n200_r0.895` | 60 / 7 | 24.3 | 25.9 | 19.9 | 0.7 | 1.4 | 1,000 | 6.7 |
| `n80_r0.65` | 18 / 9 | 20.3 | 20.5 | 15.2 | 1.0 | 1.4 | 1,000 | 5.6 |
| `n80_r0.5_b` | 14 / 13 | 20.1 | 21.4 | 14.8 | 1.2 | 1.6 | 1,000 | 5.6 |

The corrected arms (a second refit and two sets of 999 draws) take 6% of a unit. The CV fit costs about the same at n = 42 as at n = 200: its cost is the number of fits, not the sample.

**Projection.** 48 core-hours at the smoke's per-unit times, against the proposal's estimate near 40. Two references bracket it. The group-dispersion run recorded 187 core-hours on `ing` against a 233 projection from its own smoke. The stage-effect scaling run took about 1.4 times its projection, which would give 67 here.

**Proposed cap: 80 core-hours**, about 1.2 times the inflated projection. Suggested array: 150 shards (50 units each, about 20 min per task at the smoke's times), throttled at `%48`. No amendment to the designs, grade, rule or re-read is proposed.

**Seeds, checked while writing the tests.** No new design name is one of the group-dispersion run's, the univariate baseline's or the stage-effect scaling run's. The 7,000 new generator seeds are distinct. Two new designs share a committed design's generator: `n80_r0.8_b` shares `n80_r0.8`'s and `n80_r0.5_b` shares the baseline's `n80`. Neither shares a seed with its committed twin. One 31-bit seed coincides across designs whose generators differ: `n42_r0.5` replicate 469 has the seed of `band_n320_r0.2_s2` replicate 495 (n = 320, δ = 0.05, s = 2). Their datasets are not shared, so this is recorded and not amended.

**Approved 2026-11-25.** The user approved the 80 core-hour cap and the array as proposed (150 shards, `%48`), with no amendment, and the commit and push of the driver. The run starts from the pushed branch's revision in a clean worktree on `ing`.

Smoke records: `/tmp/claude-1000/-home-tomas-Documents-Repos-MOTCO/419832ff-6a00-4bff-95f9-c7693419319b/scratchpad/ud-smoke/` (scratch, not committed).
