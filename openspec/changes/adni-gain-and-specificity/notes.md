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
