# Notes

## Label-free probe of the stage fit (2026-11-16; not a result)

Run on the workstation from `main` at `1df3bff` (MOTCO 0.25.0; Python 3.11, the repository's `.venv`) before the proposal was written, serial, with BLAS pinned to one thread. The scripts were scratch files outside the repository and are not committed. The export is `adni_sexdiff/results/motco_input`, read in place; nothing participant-level was printed or written.

**Question.** Before fixing the design points: what does the production stage fit select on subsamples of ADNI's blocks, how strongly does it fit in sample, and what does a unit cost?

**Why it reads no result.** The latent fit is stage-supervised and covariate-blind, so the rank, the held-out stage AUROC and the in-sample stage R² do not depend on the group labels. Every analysis call ran on within-stage-permuted sex labels, and no trajectory statistic or p-value was printed or stored. The decision of this change (design D8–D10) reads rejection rates only.

**Construction.**
- Participants subsampled within sex × stage cells, without replacement, at the proportional tables of design D3: n = 80 (16 subsamples), 200 (12) and 500 (8), and the full cohort once.
- Scratch seeds `default_rng([20261116, n, replicate])`. They are unrelated to the driver's seeds, which derive by sha256 from base seed 2039.
- The case-study configuration through `analyze_trajectories`: blocks `p180` and `nmr` (`identity`), stage order CN / MCI / AD, covariates `age`, `education`, `apoe4`, `phase`, `statin` (the last three categorical), rank floor 2, the harness CV defaults.
- In-sample stage R²: `1 − ‖Y − Ŷ‖² / ‖Y − Ȳ‖²` with `Y` the one-hot stage labels in declared order and `Ŷ` the fitted model's prediction on its own input features, at the CV rank and at forced ranks 5, 8, 13 and 19.

**The full cohort reproduces the case study's fit.** Rank 8 at a held-out stage AUROC of 0.599514, the value in `adni_sexdiff`'s primary record (0.5995142920769246, written by MOTCO 0.10.0). In-sample stage R² 0.086.

**CV rank and held-out AUROC on subsamples.**

| n | subsamples | CV ranks selected | AUROC min / median / max |
|---|---|---|---|
| 80 | 16 | 2 (×10), 3 (×4), 4, 8 | 0.460 / 0.558 / 0.637 |
| 200 | 12 | 2 (×5), 3 (×2), 4 (×2), 5 (×2), 13 | 0.512 / 0.553 / 0.629 |
| 500 | 8 | 2 (×3), 3, 4, 5, 6, 7 | 0.536 / 0.559 / 0.592 |
| 1,501 | 1 | 8 | 0.600 |

**Median in-sample stage R² by fit rank** (range over subsamples in brackets where quoted in the design).

| fit rank | n = 80 | n = 200 | n = 500 | n = 1,501 |
|---|---|---|---|---|
| CV rank | 0.162 | 0.093 | 0.051 | 0.086 |
| 5 | 0.386 | 0.197 | 0.097 | 0.051 |
| 8 | 0.614 (0.566–0.662) | 0.345 (0.292–0.373) | 0.162 (0.154–0.189) | 0.086 |
| 13 | 0.811 | 0.500 | 0.247 | 0.126 |
| 19 | 0.933 (0.902–0.945) | 0.620 (0.573–0.667) | 0.314 (0.285–0.354) | 0.157 |

For comparison, the SEA-AD anchor's fit reads 0.840 at rank 5 on 80 donors and 10,000 features, at a held-out AUROC of 0.636.

**Timing** (one core).

| n | production CV | one forced-rank analysis at 999 permutations |
|---|---|---|
| 80 | 7.9 s | 0.34 s |
| 200 | 9.2 s | 0.35 s |
| 500 | 12.3 s | 0.37 s |
| 1,501 | 22.3 s | 0.49 s |

No analysis raised: every subsample's covariate design had full column rank.

**What it suggests for the design.**
- The stage signal is weak at every n (held-out AUROC about 0.55 to 0.60), and the CV rank is unstable below the full cohort.
- At 80 participants the CV mostly selects rank 2, the number of stages minus one. There the predictive subspace is the scores' own two dimensions in another metric, so an anchor read only at the CV rank would hardly test the guess of roadmap item 19.
- At a fixed rank the in-sample fit falls steeply with n while the held-out AUROC stays where it is. Rank 19, the highest rank the CV can select, brackets SEA-AD's 0.840 from above at n = 80; rank 8, the cohort's own, reaches 0.614.
- About 7 core-hours for 2,000 units (design D13).

**A pitfall the probe hit.** Its first pass sorted the stage labels as strings (AD, CN, MCI), which is not the order of the fitted model's columns (CN, MCI, AD), and returned negative R² values. `predictive_subspace_check.stage_r2` and the size-excess diagnosis's helpers sort their labels the same way; they were written for SEA-AD's `0` / `1` / `2`. The driver hands them each row's stage position (design D6).

**What it is not.** A handful of subsamples per n on scratch seeds, with no rejection rate. It fixed the design points, the two fixed ranks and the compute estimate, and nothing else.
