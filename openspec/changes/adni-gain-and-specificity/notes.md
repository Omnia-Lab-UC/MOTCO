# Notes

## Probes before the proposal (2026-11-17; not results)

Run on the workstation from `main` at `7e35abc` (MOTCO 0.25.0; Python 3.11.15, the repository's `.venv`) before the proposal was written, with BLAS pinned to one thread in the shell and 4 to 11 workers. The scripts were scratch files outside the repository and are not committed. The export is `adni_sexdiff/results/motco_input`, read in place through `adni_anchor_check.load_cohort`; nothing participant-level was printed or written, and the scratch outputs hold aggregates only.

**Seeds.** Scratch base seeds 20261117 (probe 1), 20261118 (probe 2) and 20261119 (probe 3), passed to the anchor driver's and the ADNI-design driver's own functions. They are unrelated to the run's seeds: the plants profile uses the anchor's family (base seed 2039) and the simulation profile base seed 2040.

**Why none is a result.** Probe 1 reads no trajectory statistic: every quantity depends on the stage labels and the covariates only, and each analysis call ran on within-stage-permuted sex labels. Probes 2 and 3 read rejection rates of both arms, on scratch seeds, at 8 to 16 replicates and 199 permutations. They fixed the plant's definition, the design points, the effect ladders and the compute estimate. The decision (design D12, D13) reads the run's units only.

### Probe 1: the reference stage configuration (label-free)

**Question.** If a plant is defined on the participants outside the analysed sample, how reproducible is that configuration, how much of it does the fitted span keep, and how large is a unit plant?

**Construction.**
- Participants subsampled within sex × stage cells at the cohort's proportions, 12 subsamples per size; the reference is the complement.
- The production stage fit through `analyze_trajectories` at the case-study configuration, at the CV rank and forced at ranks 8 and 19.
- The stage configuration in the analysis sample's pooled-standardized coordinates (checked equal to the analysis's own feature matrix to 1e-13): first the stage means, then the stage coefficients of `features ~ stage + covariates` with the six covariate columns centred.
- Reliability: the reference split in two within stage; the correlation of the two centred configurations, and the eigenvalues of their symmetrized cross-product, whose expectation is the true configuration's.
- Capture: the share of a centred configuration's energy inside `qr(x_rotations_)` (scores) and inside the predictive subspace's basis, through `LinearLatentMap`.

**Unadjusted stage means.**

| analysed | reference | reliability | λ₂ / λ₁ | held-out capture, scores, CV rank / rank 8 | own-sample capture, scores, CV rank |
|---|---|---|---|---|---|
| 80 | 1,421 | 0.50 | 0.05 | 0.46 / 0.69 | 0.91 |
| 200 | 1,301 | 0.54 | 0.16 | 0.48 / 0.72 | 0.94 |
| 500 | 1,001 | 0.48 | 0.07 | 0.71 / 0.77 | 0.99 |
| 750 | 751 | 0.39 | 0.00 | 0.73 / 0.80 | 1.00 |
| 1,000 | 501 | 0.24 | 0.00 | 0.74 / 0.76 | 1.00 |

**Covariate-adjusted stage means.**

| analysed | reference | energy, adjusted / unadjusted | cosine | reliability | λ₂ / λ₁ | held-out capture, scores | held-out capture, predictive | own capture, scores | own capture, predictive |
|---|---|---|---|---|---|---|---|---|---|
| 200 | 1,301 | 0.51 | 0.87 | 0.24 | 0.27 | 0.23 / 0.52 | 0.14 / 0.04 | 0.69 / 0.95 | 0.58 / 0.38 |
| 500 | 1,001 | 0.53 | 0.87 | 0.16 | 0.00 | 0.45 / 0.63 | 0.10 / 0.07 | 0.92 / 0.95 | 0.70 / 0.35 |
| 750 | 751 | 0.58 | 0.84 | 0.17 | 0.00 | 0.49 / 0.61 | 0.13 / 0.09 | 0.90 / 0.96 | 0.55 / 0.37 |
| 1,000 | 501 | 0.64 | 0.86 | 0.20 | 0.00 | 0.63 / 0.65 | 0.07 / 0.07 | 0.95 / 0.95 | 0.47 / 0.41 |

Capture pairs are at the CV rank and at rank 8.

**Also read at 750 participants.** CV ranks 2 to 9 across the 12 subsamples, median held-out AUROC 0.574, production CV 16.9 s on one core (20.0 s at 1,000). A unit plant on the adjusted held-out configuration: stretch, shift norm 1.04 / 0.84 / 1.08 pooled SDs at CN / MCI / AD over the 298 features, largest entry 0.22; rotation 1.48 / 1.19 / 1.52, largest entry 0.30; bend 0.46 / 0.19 / 0.47, largest entry 0.10.

**What it says.**
- The covariates account for about half of the stage configuration's energy, and the adjusted configuration is not parallel to the unadjusted one.
- The adjusted configuration is close to one-dimensional as far as the cohort can show: the second eigenvalue of the halves' cross-product is zero from 500 participants up.
- A held-out adjusted configuration has a split-half reliability of 0.16 to 0.24. By Spearman–Brown about 0.3 of the full reference's energy is signal.
- The sample's own configuration lies in the fitted span almost entirely; a held-out one does not.
- The energy share understates the predictive subspace, which keeps 0.07 to 0.14 of a held-out configuration at the CV rank and still read a held-out stretch about as often as the scores in probe 2. Its two directions are the ones along which the fitted stage labels move, not the ones along which the stage means differ.

### Probe 2: plants on ADNI's noise (rates, both arms)

**Question.** At which effects is each kind of plant read, on each reference, and what does a unit cost?

**Construction.** The anchor driver's units (subsample, within-stage label permutation) at 500 and 750 participants and at the full cohort. The plant of design D4 on the second group's rows. `pls_cv` cross-validated and `predictive_orthonormal` forced at its rank, `df`, 199 permutations, rejection at p ≤ 0.05.

**First pass: unadjusted stage means** (16 replicates at 500 and 750, 12 at the cohort). Recorded for the pitfall it found.

| cell | arm | `delta` / `angle` / `shape` rejections |
|---|---|---|
| `cohort`, in-sample stretch 2.0 | `pls_cv` | 12 / 9 / 4 of 12 |
| | `predictive_orthonormal` | 12 / 1 / 1 of 12 |
| `n750`, in-sample stretch 2.0 | `pls_cv` | 16 / 10 / 6 of 16 |
| | `predictive_orthonormal` | 16 / 0 / 4 of 16 |
| `n750`, held-out rotation toward a random direction, 1.0 | `pls_cv` | 16 / 13 / 15 of 16 |
| | `predictive_orthonormal` | 15 / 10 / 15 of 16 |

A stretch of unadjusted stage means is not a stretch of the covariate-adjusted configuration the trajectory model measures, so production's `angle` rejections here are not a leak. The second pass built every plant on the adjusted configuration.

**Second pass: covariate-adjusted configuration** (16 replicates per cell). Rejections of 16 for `delta` / `angle` / `shape`, `pls_cv` then `predictive_orthonormal`:

| cell | `n750`, `sample` | `n750`, `heldout` | `cohort`, `sample` |
|---|---|---|---|
| stretch 0.25 | 3 / 1 / 1, 5 / 1 / 1 | 2 / 1 / 1, 0 / 0 / 2 | 5 / 1 / 1, 9 / 1 / 1 |
| stretch 0.5 | 10 / 2 / 1, 11 / 0 / 1 | 5 / 2 / 1, 2 / 0 / 2 | 14 / 1 / 1, 15 / 1 / 2 |
| stretch 1.0 | 16 / 3 / 2, 16 / 0 / 1 | 13 / 3 / 1, 13 / 2 / 1 | 16 / 1 / 0, 16 / 2 / 4 |
| stretch 2.0 | 16 / 5 / 3, 16 / 0 / 0 | 16 / 7 / 3, 16 / 6 / 3 | 16 / 3 / 3, 16 / 1 / 4 |
| rotation 0.25 | 0 / 1 / 3, 3 / 1 / 1 | 0 / 0 / 1, 0 / 0 / 0 | 0 / 1 / 2, 0 / 6 / 1 |
| rotation 0.5 | 0 / 6 / 3, 2 / 6 / 2 | 1 / 1 / 3, 2 / 2 / 3 | 0 / 2 / 4, 0 / 6 / 2 |
| rotation 1.0 | 2 / 4 / 8, 2 / 1 / 4 | 12 / 6 / 7, 13 / 6 / 5 | 0 / 6 / 6, 0 / 1 / 1 |
| bend 0.5 | 0 / 1 / 3, 2 / 1 / 0 | 1 / 1 / 1, 0 / 0 / 0 | 0 / 0 / 3, 1 / 1 / 6 |
| bend 1.0 | 0 / 2 / 6, 3 / 1 / 5 | 1 / 1 / 1, 1 / 1 / 1 | 0 / 1 / 6, 1 / 1 / 10 |
| out-of-span rotation 0.25 | 3 / 2 / 1, 3 / 3 / 11 | — | 6 / 9 / 1, 5 / 14 / 15 |
| out-of-span rotation 0.5 | 12 / 11 / 7, 15 / 7 / 16 | — | 16 / 16 / 8, 14 / 16 / 16 |
| translation 1.0 | 0 / 1 / 0, 0 / 0 / 1 | — | 0 / 0 / 0, 0 / 1 / 0 |

**The fit under a plant** (median CV rank, median held-out AUROC; the anchor reads rank 8 and 0.60 at the cohort):

| cell | `n750`, `sample` | `cohort`, `sample` |
|---|---|---|
| stretch 1.0 | 13, 0.75 | 18, 0.75 |
| rotation 1.0 | 2, 0.51 | 3, 0.52 |
| out-of-span rotation 0.25 | 19, 0.64 | 19, 0.68 |

**Scale drift** (the largest relative change of a pooled SD, median unit): 0.1 to 1.6 % at effects up to 1.0, and 3.2 to 4.1 % at a stretch of 2.0.

**Timing.** With 10 or 11 workers on 16 threads the production CV took 20 to 34 s at 750 participants and 33 to 57 s at the cohort; a forced call took under a second. The compute estimate of design D17 uses probe 1's 17 s and item 21's 21.5 s for the cohort, both measured at six or seven workers.

**What it says.**
- An in-sample stretch is read from 0.25 to 0.5 and a held-out one from 0.5 to 1.0; a held-out stretch of 0.5 is read about half as often as an in-sample one.
- Rotations and bends are read rarely by production at every effect. A held-out rotation of 90° is read by `delta`, which is not a leak: a held-out plant is mostly the reference's noise (probe 1), so its intended geometry is not the declared one.
- A rotation of 90° in sample cancels half of the pooled stage signal: the AUROC falls to 0.52 and the CV to rank 2 or 3.
- A rotation toward a direction outside the stage span does not stay outside the fitted span on ADNI's noise. The CV selects rank 19, and both arms reject. The candidate reads it as `shape`.
- The translation reads at the null rate for both arms.

### Probe 3: stage-span plants at the ADNI-design simulation

**Question.** Does the stage-span construction compose with item 9's construction, and at which effects is it read?

**Construction.** `adni_design_check.build_settings` at its defaults (stage effect 0.02, `s_z` 0.1, `s_p` 1.75) with `construction: "stage_span"`, `span_layers` methylation and expression, `span_scale: "planted"`; ADNI's table; `pls_cv` and `predictive_orthonormal` forced at its rank through `analyze_trajectories` with the six covariates and `df`; 8 replicates at 199 permutations.

| cell | `pls_cv` | `predictive_orthonormal` | truth's intended `delta` / `angle` / `shape` |
|---|---|---|---|
| anchor | 2 / 1 / 1 | 2 / 1 / 0 | 0 / 0° / 0 |
| stretch 0.5 | 2 / 0 / 0 | 4 / 0 / 0 | 0.35 / 0° / 0 |
| stretch 1.0 | 3 / 0 / 0 | 7 / 0 / 0 | 0.69 / 0° / 0 |
| stretch 2.0 | 8 / 0 / 0 | 8 / 0 / 2 | 1.38 / 0° / 0 |
| rotation 0.5 | 2 / 1 / 1 | 2 / 3 / 1 | 0 / 45° / 0 |
| rotation 1.0 | 2 / 1 / 1 | 3 / 7 / 3 | 0 / 90° / 0 |
| bend 0.5 | 2 / 0 / 1 | 2 / 0 / 0 | 0 / 0° / 0.23 |
| bend 1.0 | 1 / 0 / 1 | 2 / 0 / 1 | 0 / 0° / 0.36 |

Rejections of 8 for `delta` / `angle` / `shape`.

- **Rank and signal.** The CV selected rank 19 in all 64 units; the anchor's median held-out AUROC was 0.61. The population path length is 0.69 in planted-standardized units, so a unit stretch's shift has a norm of about 0.2 at each stage, against about 1.0 for a unit stretch of ADNI's adjusted held-out configuration, most of which is noise.
- **Larger effects build.** A bend of 2.0, a rotation of 2.0 and a stretch of 3.0 generate; a bend of 2.0 has an intended `shape` of 0.51.
- **The anchor is the parent's.** At a seed the stage-span anchor's blocks and metadata equal those of `adni_design_check`'s anchor exactly.
- **The parent reproduces.** Under MOTCO 0.25.0 the parent's own evaluator on its committed settings gives the committed p-values for anchor replicates 0 and 1: 0.29 / 0.153 / 0.195 and 0.187 / 0.746 / 0.387.
- **Timing.** A unit took 35 s on an idle workstation core and about 100 s with twelve processes running. Item 9 measured 90 s on `ing` (174.6 core-hours for 7,000 units), which the estimate uses.

### What the probes suggest for the design

- Build every plant on the covariate-adjusted configuration (design D3).
- Read specificity on in-sample plants, whose truth is exact, and fairness on held-out plants, whose direction no fit was built on (design D5).
- Two design points: 750 participants, the largest sample with a held-out reference of its own size, and the full cohort (design D2).
- Effect ladders that put production between the floor and the ceiling (design D6, D9).
- Record the out-of-span cells and do not judge them (design D7). This is the one choice the probe's rates bear on: the candidate's `shape` readings there would fail a specificity criterion that included them.
- About 40 core-hours on the workstation and 70 on the cluster (design D17).

### Pitfalls the probes hit

- **Unadjusted stage means.** The first pass of probe 2 produced apparent off-target rejections that were the plant's own (above).
- **The energy share as a measure of what the predictive subspace sees.** It is small by construction and does not track the rates (probe 1).
- **A random direction is not outside the fitted span once planted.** On SEA-AD a sparse out-of-span plant kept 2.2 % of its norm in the span (item 12). At ADNI's sample sizes the stage fit moves toward the plant.

## Smokes (task 7.1, 2026-11-17)

Code: the worktree `MOTCO-adni-gain` on `proposal/adni-gain-and-specificity` at `7e35abc` plus this change's uncommitted tree (tree hash of `src`, `scripts`, `tests` `137d8cf55569b837`, the same on both machines at the time of the smokes). Workstation: Python 3.11.15, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3. `ing`: Python 3.11.16 with the same package versions (the main checkout's `.venv`). BLAS pinned to one thread in the shell on both. Outputs stay off the repository: `~/motco-runs/adni-gain-smoke/out` and `ing:/home1/tgonza/MOTCO-adni-gain-smoke/out`.

### Plants profile, on the workstation (the export read in place)

- **Pre-flight passed.** On cohort replicates 0 and 1 the cross-validated fit selected rank 8 at a mean held-out AUROC of 0.599514, and the call forced at rank 8 gave the same statistics and p-values (CV 22.0 s and 21.6 s, forced 0.5 s).
- **Run.** Four replicates per cell at 999 permutations, seven workers: 120 of 120 units completed, 0 failed, 0 inconsistent, no subsample or label redraw; 0.79 recorded core-hours. Every check passed on every unit, `plant` included.
- **The cohort anchor is the anchor run's unit.** The four cohort anchor units equal the anchor driver's own smoke units (`~/motco-runs/adni-anchor-smoke`, same seeds, 999 permutations) in all four arms' statistics and p-values, 16 of 16 arm records. The merge's reproduction against the committed 500-unit counts fails at four units, as it must; the smoke's decision is `gate_failed` for that reason and for the simulation profile it did not see.
- **Timing (median seconds per unit, seven workers).** `n750`: 21.0 for controls, 20.6 for plants (the production CV 19.5; each forced arm about 0.5). `cohort`: 2.3 for the forced anchor, 30.9 for the translation and 31.3 for plants.
- **Projection for the full run** (7,200 units at those medians): `n750` controls 5.8 and plants 18.3 core-hours, `cohort` anchor 0.3, translation 4.3 and plants 17.4: **46.1 core-hours, about 6.6 h at seven workers** (the design estimated 40 and six).
- **Descriptive readings at four units** (medians; none decides anything): CV rank 5.5 on the `n750` anchor and 11.5 under an in-sample stretch of 1.0 (AUROC 0.59 and 0.76); 2.5 at the cohort under a 90° in-sample rotation (AUROC 0.52) and 19 under an out-of-span rotation of 0.25 (AUROC 0.68). The fitted span keeps 0.95 of the sample's configuration and 0.48 of the held-out one at `n750`, the predictive subspace 0.14 of the held-out one; it keeps 0.49 of an out-of-span plant of 0.25 at the cohort. Scale drift 0.2 to 1.8 % at the cells listed.

### The held-out reliability: what the probe's 0.17 was

The smoke's four `n750` anchor units read a split-half correlation of −0.24 to 0.18 and an eigenvalue ratio of 0.22 to 0.54, against the probe's 0.17 and 0.00 (design Context and D5). Two findings, on label-free readings of `n750` subsamples outside the run's replicates (12 for the eigenvalues, replicates 500 to 511; 24 for the correlations, 500 to 523; no trajectory statistic read):

- **The eigenvalues belong to the stage-space cross-product.** The driver first took `(AᵀB + BᵀA)/2` over features, whose eigenvalues carry ± noise pairs of the order of ‖A‖·‖B‖. In stage space, `(A·Bᵀ + B·Aᵀ)/2`, the second eigenvalue is at most zero in 4 of 12 subsamples and the median ratio is 0.000, the probe's reading. The driver now uses the stage-space product, and a test pins it (noise-free halves give correlation 1 and the configuration's squared singular values). The smoke's four records carry the old reading; nothing else in them changes.
- **The probe's reliabilities were Spearman–Brown step-ups.** The mean split-half correlation is 0.086 adjusted and 0.238 unadjusted (medians 0.023 and 0.180); stepped up by 2r/(1 + r) they are 0.16 and 0.38, the probe's 0.17 and 0.39. So 0.17 is already the full held-out reference's reliability, and D5's "about 0.3 by Spearman–Brown" applies the step-up twice: **about 0.16 of a held-out configuration's energy is the cohort's stage program, not a third.** The rules are unaffected (no held-out off-target reading is judged, and G reads the simulation alone); the held-out plants are more noise-dominated than D5 says. Each record now carries both the split-half correlation and its step-up (`full_reference`), and the cell table reports both medians. The report should quote these, and design D5's sentence is the one to correct when the artifacts are next revised.

### Simulation profile, on `ing`

Run from the scratch copy `ing:/home1/tgonza/MOTCO-adni-gain-smoke` (the tree above, rsynced, with the four parent files byte-identical to the committed ones; the main checkout's `.venv` symlinked; nothing else on `ing` touched).

- **The test module on a compute node:** 87 of 88 passed. The one failure built the anchor driver's own settings, which bind the within-stage run's tables, and those were not copied to the scratch directory; no code is involved.
- **Pre-flight passed** (`srun`, 2 min 10 s): anchor replicates 0 and 1 of the ADNI-design run, through its own evaluator under its committed settings, gave its committed p-values exactly (0.29 / 0.153 / 0.195 and 0.187 / 0.746 / 0.387).
- **Timing shard:** two replicates per cell, 20 units over ten array tasks: 20 of 20 completed, 0 failed, 0 inconsistent; 0.26 recorded core-hours. The CV selected rank 19 in 19 units and rank 2 in one 90° rotation (as the probe saw: a full rotation cancels the pooled stage signal).
- **n4 stalled.** The eight tasks first placed on n4 (array 944507) sat in uninterruptible NFS page reads at Python import for 14 minutes (28 CPU-seconds, no unit); they were cancelled and resubmitted with `--exclude=n4,n11` (array 944578), where every task completed in 2.5 to 3.6 minutes. The full run should exclude n4 as well as n11.
- **Unit time:** 43.4 to 49.6 s per unit on n2 and n3 (median by cell 44 to 49 s), against the 90 s the design assumed from item 9.
- **Projection for the full run** (2,600 units): **33 core-hours at these unit times.** Item 9 and the design-guidance runs measured units 2.5 times slower on packed nodes and 5 to 13 minutes of Python startup per task when many start at once; at the packed-node cost the run would take about 83 core-hours. A cap of 90 core-hours covers both.
- **Combined decision.** With the simulation aggregates copied back (settings, pre-flight, tables, accounting; no records), `decide` on the workstation read both profiles: `gate_failed`, for the reproduction alone (four cohort units cannot match 500-unit counts), with every criterion evaluated. The machinery spans both machines as the full run needs.
