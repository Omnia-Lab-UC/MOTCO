# Design

## Context

See `proposal.md` (Why) for the motivation. The facts below shape the approach.

- **The export and the case study's configuration** are those of the ADNI anchor (`openspec/changes/archive/2026-11-17-adni-real-noise-anchor/design.md`, Context): 1,501 participants, the `p180` (133) and `nmr` (165) blocks, Female 237 / 310 / 125 and Male 223 / 448 / 158 at CN / MCI / AD, five covariates in six design columns, rank floor 2, net displacement, the harness CV defaults, rejection at p ≤ α. `scripts/adni_anchor_check.py` reads it through `ADNI_CONTRACT` and builds subsampled units on within-stage-permuted labels.
- **What item 21 read at the full cohort** (rejections of 500, `results/adni-real-noise-anchor-2026-11-17/rate_table.csv`), on one fit forced at rank 8:

  | arm | `delta` | `angle` | `shape` |
  |---|---|---|---|
  | `pls_cv` | 25 | 32 | 28 |
  | `predictive_orthonormal` | 29 | 20 | 18 |
  | `pls_cv__within_stage` | 35 | 38 | 28 |
  | `predictive_orthonormal__within_stage` | 27 | 28 | 21 |

  Its decision gives `predictive_orthonormal` a size scope from 80 participants to the cohort and declares the scoping rule S1–S5: a reading speaks for its export and design points, nothing is interpolated between points, and size alone reopens nothing.
- **What item 9 read in simulation** (`results/adni-design-2026-11-02/run/`): at ADNI's table with six covariate columns, a phase-like nuisance and a stage effect calibrated to ADNI's AUROC, size holds and the CV selects rank 19 in 494 of 500 anchor units. Its plants are feature surgeries: `delta` reaches 0.530 at a doubled stage effect, and `angle` and `shape` stay at the null rate at every effect. Its per-unit p-values are committed.
- **The generator's stage-span edit is a public function.** `motco.simulations.semisynthetic.edit_stage_configuration(z, mode, effect, alignment_degrees=, out_of_span=)` takes a stage configuration (stages × features) and returns the stretched, rotated, bent or translated one. Its companion `out_of_span_direction` needs features whose stage means are exactly constant, which real blocks do not have.
- **What the probes read** (`notes.md`; scratch seeds, 8 to 16 replicates, 199 permutations where a rate is quoted; none is a result):

  | at 750 participants | value |
  |---|---|
  | energy of the covariate-adjusted held-out configuration over the unadjusted one | 0.58 (cosine between them 0.84) |
  | reliability of the adjusted configuration between two halves of the 751 held out | 0.17 |
  | second over first eigenvalue of the halves' cross-product | 0.00 |
  | share of the sample's own adjusted configuration inside the fitted span, CV rank / rank 8 | 0.90 / 0.96 |
  | share of the held-out adjusted configuration inside it, CV rank / rank 8 | 0.49 / 0.61 |
  | the same two shares inside the predictive subspace at the CV rank | 0.55 / 0.13 |
  | norm of a unit stretch's shift at CN / MCI / AD, pooled SDs over 298 features | 1.04 / 0.84 / 1.08 (largest entry 0.22) |
  | production CV, one core | 17 s (22 s at the cohort) |

  Rejections of 16 for `delta` / `angle` / `shape`, production then `predictive_orthonormal`:

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

  In simulation at ADNI's table (rejections of 8; rank 19 in every unit): stretch 0.5 / 1.0 / 2.0 by `delta` 2 / 3 / 8 against 4 / 7 / 8; rotation 0.5 / 1.0 by `angle` 1 / 1 against 3 / 7; bend 0.5 / 1.0 by `shape` 1 / 1 against 0 / 1. The anchor of the stage-span construction is byte-identical to item 9's at a seed, and item 9's committed p-values for anchor replicates 0 and 1 reproduce under MOTCO 0.25.0.
- **The rules.**
  - `judge_rates` at R = 500: a reading *exceeds* above 34 rejections; a family of 6 or 9 readings holds with at most one exceedance and none above 37 or 38. At R = 200 a reading exceeds above 16, and a family of 32 holds with at most two exceedances and none above 20.
  - The pilot's paired rule (`predictive_subspace_check.judge_paired`): on m one-sided exact paired p-values, at most the 95th percentile of Binomial(m, 0.025) below 0.025 (two at m = 16, 30 and 32) and none below 0.05 / m.
  - The gain rule (`production_latent_space_check.judge_gain`): at least one exact paired gain p-value at most 0.05 / m.
- **Data rules.** ADNI is controlled-access. `adni_sexdiff` forbids copying its MOTCO inputs between machines, and MOTCO's repository is public. The simulation uses ADNI's cell table and phase shares only, which `adni_design_check.TABLES` already commits.

## Goals / Non-Goals

**Goals:**
- A reading of gain, protected power and specificity for `predictive_orthonormal` against production at ADNI's design, in simulation and on ADNI's noise, under a rule fixed before any unit runs.
- The first real-noise off-target readings whose truth is exact, and the first real-noise plants no fitted space was built on.
- A statement of what production detects on ADNI's noise that the case study can cite, whatever the verdict.
- Functions a later cohort's driver can import unchanged.

**Non-Goals:**
- Ranking measurement spaces in general. The reading is for one candidate at one cohort's design.
- Explaining a failure. A failing verdict names its readings and stops.
- Estimating ADNI's true stage program. A reference configuration is an input to a plant, not a result.
- Editing a committed driver. Parents are imported and bound, never changed.

## Decisions

### D1. One driver, two profiles, on existing options

`scripts/adni_gain_specificity_check.py` with `--profile plants | simulation`. Every arm is one `analyze_trajectories` call. `motco` stays at 0.25.0.

- **Why one driver.** One decision reads both profiles, as the predictive-subspace and within-stage pilots' did. The profiles share the arms, the records, the tables and the criteria.
- **Why two profiles.** The roadmap's reopening condition names a gain in simulation and specificity in both. The plants profile must run where the export is; the simulation profile may run on the cluster.
- **Alternative: two changes.** The simulated gain alone decides nothing, and the real-noise reading alone cannot meet the condition. Split, one would wait for the other with its rule already written.

### D2. Plants units are the anchor's units

The plants profile imports from `adni_anchor_check`, unchanged: `ADNI_CONTRACT`, `load_cohort`, the design-point tables, `build_unit` (the subsample within cells, the redraw for an absent covariate level, the within-stage label draw and its redraw), `unit_seeds`, the arm parameters and the per-unit checks. Its settings carry the anchor's seed family (base seed 2039).

| point | participants | Female CN / MCI / AD | Male CN / MCI / AD | held out |
|---|---|---|---|---|
| `n750` | 750 | 118 / 155 / 63 | 111 / 224 / 79 | 751 (119 / 155 / 62 and 112 / 224 / 79) |
| `cohort` | 1,501 | 237 / 310 / 125 | 223 / 448 / 158 | none |

- **Why the anchor's seeds.** The cohort anchor's 500 units are then item 21's units, so their counts must reproduce its committed table (D10), and every cohort plant is added to a unit whose null was read there. The driver convention asks later drivers to import earlier ones and their seed functions unchanged.
- **Why `n750`.** It is the largest sample whose held-out reference is as large as itself. It is not on the anchor's ladder, and S1 of item 21 forbids interpolating between 500 and the cohort, so its anchor is read here (D12).
- **Why not 80, 200 and 500.** At 80 and 200 the CV selects rank 2 in 57 % and 46 % of units (item 21), where the candidate measures the scores' own two dimensions, and the held-out AUROC is 0.52 to 0.54. A third point at 500 would add about 20 core-hours for a second held-out reading 250 participants below the first.
- **Own seeds.** Two purposes the anchor lacks derive from this driver's name and base seed 2040, used by no committed driver: `direction` (D7) and `reference_split` (D16).

### D3. The stage configuration is covariate-adjusted

A configuration is the K stage means of the pooled-standardized joint features adjusted for the case study's covariates: the stage coefficients of the least-squares fit `features ~ stage indicators + covariate columns`, with the covariate columns coded as `analyze_trajectories` codes them and centred on the reference's own means. The group label is not in the model.

- **Why adjusted.** The trajectory test compares covariate-adjusted cell means. A plant that edits unadjusted stage means is the claimed edit of a configuration the test does not measure. On ADNI the two differ by half: the adjusted held-out configuration has 0.58 of the unadjusted one's energy at a cosine of 0.84. The first probe built its stretch on unadjusted means, and production's `angle` rejected a stretch of 2.0 in 9 of 12 cohort units; on the adjusted configuration it rejected 3 of 16.
- **Only the centred configuration matters.** A stretch, a rotation and a bend depend on the stage means' deviations from their centroid, and the translation on their net displacement and mean step. Another covariate reference row moves every stage by the same vector.
- **A test pins it.** The adjusted stage means of the planted group's rows before and after a plant differ by the plant's shift exactly, because the shift is constant within a stage.
- **Alternative: unadjusted stage means,** as in SEA-AD's span plant (item 17), which was read for `delta` power only. Rejected for the reason above.

### D4. The plant

With `z` a configuration in the analysis sample's pooled-standardized coordinates (the per-feature mean and SD of the unplanted sample, the scaler the analysis fits), the plant is `edit_stage_configuration(z, mode, effect) − z`, mapped back through the pooled SD and added to every row of the second declared group at each stage, under the unit's permuted labels.

| mode | edit | intended difference |
|---|---|---|
| `magnitude` (stretch) | scale about the centroid by 1 + e | `delta` |
| `orientation` (rotation) | rigid rotation by 90°·e in the plane of the net displacement and the configuration's second direction | `angle` |
| `shape` (bend) | the middle stage moves e mean steps off the chord between its neighbours, at fixed path length | `shape` |
| `translation` | every stage moves e mean steps along the net displacement | none: a group main effect |

- **The planted group is the second one,** as in SEA-AD's span plant. Under permuted labels the choice is immaterial.
- **Covariates are not touched.** They stay on their rows, and the plant is constant within a group × stage cell, so the covariate adjustment leaves it whole.
- **The null scale.** The edit is exact in the unplanted sample's coordinates. The analysis standardizes the planted sample, whose pooled SD the plant moves: by about 1.5 % on the most affected feature at effects up to 1.0 in the probe (the median unit), and by 3 to 4 % at a stretch of 2.0. A stretch is exact under any per-feature rescaling; a rotation and a bend are exact to that order. Each unit records the largest relative change.
- **Alternative: port the generator's planted scale.** Its fixed point is defined on population cell means and a known within-cell variance. A sample version would be a new construction with its own identity checks, for a correction below the sampling noise of 200 units.
- **The plant moves the stage fit.** Half of it enters the pooled stage means. In the probe a stretch of 1.0 raised the cohort's held-out AUROC from 0.60 to 0.75 and its CV rank from 8 to a median of 18, and a 90° rotation lowered them to 0.52 and a median rank of 3. That is what a real difference does to a stage-supervised fit, so every planted unit runs its own cross-validation (D8).

### D5. Two references for a plant

| | `sample` | `heldout` |
|---|---|---|
| configuration of | the unit's analysis sample | the participants outside it, standardized by the sample's scaler |
| points | `n750`, `cohort` | `n750` |
| the planted group's expected configuration | the declared edit of the other group's, exactly | the other group's plus the edit of an estimate |
| direction | in the fitted span by construction (0.90 kept) | independent of every fit (0.49 kept) |
| read for | specificity; protected power; power as an upper bracket | protected power; recorded gain; power as a lower bracket |

- **Why `sample` is exact.** Each unit permutes the group labels within stage, so each group is a random part of each stage and its expected adjusted configuration is the sample's. Adding the edit of that configuration to one group makes its expected configuration the edit of the other's. A stretch then has no intended `angle` or `shape`, and a rotation no intended `delta` or `shape`, up to D4's scale order.
- **Why `sample` is not enough.** The stage-PLS fit is built on the sample's stage means, so the plant lies where every stage-supervised space already looks. The roadmap excludes SEA-AD's in-sample plant from ranking spaces for this reason.
- **Why `heldout` is fair.** No analysed participant contributes to its stage means or its covariate adjustment. A space cannot have been fitted on it.
- **Why `heldout` is not exact.** Between two halves of the 751 held out, the adjusted configuration's reliability is 0.17, so by Spearman–Brown about 0.3 of the full reference's energy is the cohort's stage program and the rest is estimation noise. A held-out stretch therefore carries an intended rotation of its own, and its off-target rates are not specificity readings: in the probe a held-out rotation of 90° was rejected by `delta` in 12 of 16 units by production and 13 by the candidate.
- **The bracket.** A difference along the cohort's true program lies between the two: the sample's configuration is that program plus noise the fit shares, the held-out one is the program plus noise it does not. Each target rate is reported on both references side by side.
- **Alternative: held-out plants only,** the roadmap's first wording. It leaves no exact off-target truth and no plant at the cohort.
- **Alternative: a planted stage program on stage-permuted noise.** Permuting the stage labels removes the real program, so a planted one is the population's and every truth is exact. It also removes each stage's own noise and the covariates' association with stage, the two things items 19 and 20 found the null to depend on. It stays the lead if the two references disagree.
- **Alternative: one fixed half as the reference for every unit.** Every rate would be conditional on one split, and subsamples of the other half would overlap more.

### D6. Cells, effects and replicates

| cells | reference | effects | replicates | points |
|---|---|---|---|---|
| `anchor` | — | — | 500 | both |
| `translation` | `sample` | 1.0 | 500 | both |
| stretch | `sample` | 0.25, 0.5, 1.0 | 200 | both |
| rotation | `sample` | 0.25, 0.5, 1.0 | 200 | both |
| bend | `sample` | 0.5, 1.0 | 200 | both |
| stretch | `heldout` | 0.5, 1.0, 2.0 | 200 | `n750` |
| rotation | `heldout` | 0.5, 1.0 | 200 | `n750` |
| bend | `heldout` | 1.0 | 200 | `n750` |
| out-of-span rotation (D7) | `sample` | 0.25, 0.5 | 200 | both |

4,200 units at `n750` and 3,000 at `cohort`, all at 999 permutations. Cells are named `<reference>_<mode>_e<effect>`. The cells of one design point at one replicate share the sample, the labels and the permutation seed.

- **Why these effects.** They put production's target rate between the floor and the ceiling in the probe: an in-sample stretch reads 3, 10 and 16 of 16 at 0.25, 0.5 and 1.0, a held-out one 5, 13 and 16 at 0.5, 1.0 and 2.0. Rotations and bends stay low for production at every effect (at most 6 and 6 of 16 in sample), which is where a gain would show; a rotation of 90° is kept because it cancels the pooled stage signal and is production's largest off-target reading.
- **Why the translation is a control.** It moves every stage of one group by the same vector, so there is no interaction and all three statistics are off-target, whatever the reference's noise. It is the stage-span construction's hardest negative control and joins the size family (D12).
- **Why 500 and 200.** The controls are judged by counts against a size bound and need the anchor's grade. The plants are judged by exact paired tests on the discordant units. An excess of 0.05 in rate with few reverse discordances is below 0.025 at 200 units (12 units against 2) and not at 100 (6 against 1).
- **Why a held-out bend at one effect.** The probe read 1 of 16 at both effects for both arms. One cell keeps a bend in the protected-power family.

### D7. The out-of-span rotation is recorded, not judged

`o` is a seeded standard normal vector over the joint features (purpose `direction`), projected off the configuration's centred stage means and normalized. The cell is `edit_stage_configuration(z, "orientation", e, alignment_degrees=90, out_of_span=o)`: the planted group's net displacement turns toward `o` by 90°·e. The two effects of a replicate share `o`.

- **Why it is in.** It is the fourth cell open item 2 lists, and the only one where the two measurements can differ in what they keep of a difference outside the stage plane.
- **What the probe read.** On ADNI's noise the fit follows the plant: the CV selects rank 19, and both arms reject. At the cohort at 0.25, production reads 6 / 9 / 1 of 16 and the candidate 5 / 14 / 15. The candidate keeps `angle` and reads the rotation as `shape` as well.
- **Why it is not judged.** The measurement-visibility page states that a difference outside the stage span reaches the statistics by incidental routes with a mixed signature, for production too, and that the split test's complement is the detector there. A specificity criterion on these cells would judge the candidate against a property production is documented not to have. The predictive-subspace pilot treated its out-of-span cell the same way.
- **What is reported.** Both arms' rates for the three statistics, the paired counts, the CV rank and the share of the plant inside each measured space. The standing sentence of every verdict points to them (D13).
- **Alternative: judge them in the specificity criterion.** On the probe's reading the candidate would fail there. It is recorded as the stricter rule; adopting it is a change of rule version before any unit runs.
- **Why not `out_of_span_direction`.** It draws on features whose stage means are exactly constant. Real blocks have none.

### D8. Four arms on one fit

| arm | coordinates | draws | role |
|---|---|---|---|
| `pls_cv` | scores | unrestricted | reference |
| `predictive_orthonormal` | predictive, orthonormal | unrestricted | candidate |
| `pls_cv__within_stage`, `predictive_orthonormal__within_stage` | as named | within stage | recorded |

Every arm runs the case study's configuration with the `df` residual correction (κ = 1.007 at 750 participants and 1.003 at the cohort) at the unit's permutation seed.

- **One cross-validation per unit.** `pls_cv` cross-validates on the unit's data, planted or not, and the other three arms are forced at the rank it selects, so their final fit is the CV's.
- **The cohort anchor is forced.** Its fit does not see the labels, so its 500 units run forced at the pre-flight's rank 8, as in item 21.
- **Why the within-stage arms.** They cost a forced call each. Item 20 read that null `specificity_excess` in simulation and item 21 saw production fail under it at the cohort; its real-noise off-target rates on exact plants are the reading item 20's verdict lacks. They decide nothing here.
- **Alternative: the fixed-rank arms of item 21.** No case study forces a rank, and their grade was `reference_exceeds`.

### D9. The simulation profile

The construction is `adni_design_check`'s, imported unchanged: its settings at the calibrated strengths (stage effect 0.02, `s_z` 0.1, `s_p` 1.75) and ADNI's table, its unit builder and its seed derivation. Three generator settings are added: `construction: "stage_span"`, `span_layers` the two measured blocks, `span_scale: "planted"`.

| cells | effects | replicates |
|---|---|---|
| `anchor`, `translation` (1.0) | — | 500 |
| stretch | 0.5, 1.0, 2.0 | 200 |
| rotation | 0.25, 0.5, 1.0 | 200 |
| bend | 1.0, 2.0 | 200 |

2,600 units at 999 permutations, base seed 2040. The four arms of D8, with the parent's six covariates.

- **Why stage-span plants.** The feature surgeries' rotations and bends keep almost none of their energy in the stage-supervised span (stage-span pilot), so item 9's null-rate `angle` and `shape` measure visibility. The stage-span construction edits the population stage means where the measurement looks, and its truth records the intended `delta`, `angle` and `shape`. It is the simulated counterpart of D4.
- **Why the planted scale.** Under it a rotation and a bend have no intended `delta` in the coordinates the scaler estimates (roadmap item 13 (a)), so the off-target truth is exact.
- **Why these effects.** The simulated program is weak in pooled-SD units (population path length 0.69). In the probe production's `delta` reads a stretch in 2, 3 and 8 of 8 units at 0.5, 1.0 and 2.0. Nothing reads a bend up to 1.0, so the ladder goes to 2.0.
- **The anchor is item 9's null dataset.** The stage-span construction draws no random number, so at a seed its anchor is byte-identical to the parent's. This is what lets the pre-flight bind the two runs (D10).
- **`df` on every arm.** Item 9 ran without it, as the case study's primary was issued. Here every arm of both profiles carries it, so the two profiles read one configuration.
- **Alternative: a second table at 750 participants,** to set the simulation beside the held-out plants at one size. The rule does not need it, and it adds about 35 core-hours. It is a lead if the profiles disagree.
- **Alternative: a band of stage effects at ADNI's table.** Those are other designs. The predictive-subspace pilot read a band at its own tables.
- **Where it runs.** On `ing`, from a clean worktree at a pushed revision, through `scripts/adni_gain_specificity_array.sbatch`. It reads no ADNI file.

### D10. Pre-flights and the reproduction gate

- **Plants pre-flight.** As item 21's: on cohort replicates 0 and 1 the cross-validated fit must select rank 8 at a mean held-out AUROC of 0.599514 to six decimals, with statistics and p-values equal to the call forced at rank 8. `run --profile plants` refuses to start without it.
- **Simulation pre-flight.** Anchor replicates 0 and 1 of item 9, analysed by the parent's own unit evaluator under its committed settings, must give its committed p-values exactly (0.29 / 0.153 / 0.195 and 0.187 / 0.746 / 0.387 for `delta` / `angle` / `shape`). `run --profile simulation` refuses to start without it.
- **Anchor reproduction.** At the merge of the plants profile, the cohort anchor's rejection counts for the four arms must equal the table of the Context above.
- **Parents, bound by sha256** in the settings and refused by name when changed: the anchor run's `rate_table.csv` (the reproduction) and `adni_anchor_decision.json` (the size scope the `reopened` consequence cites), and the ADNI-design run's `settings.json` and `adni_design_pvalues.csv` (the simulation pre-flight). The export's five file hashes are bound as in item 21.
- **Why.** Each check ties the export, the row order, the code and the environment to a committed run before a new number is read. A mismatch is diagnosed, never waived.

### D11. Checks on every unit

- From the anchor driver, unchanged: `rank`, `same_fit`, `measured_dimension`, `same_statistics` (each within-stage arm against its unrestricted counterpart) and `cell_table`.
- `plant` (plants profile): after planting, the second group's raw stage means have moved by the declared shift and the first group's have not, to 1e-8.
- `construction` (simulation profile): the dataset's truth records the stage-span construction at the cell's mode and effect.

A unit failing a check is `inconsistent`, names the check, the arm and the deviation, and leaves every rate.

### D12. The criteria

Rejection is p ≤ 0.05. Every criterion is evaluated whichever fails first.

| criterion | family | readings | holds with |
|---|---|---|---|
| S1 size | plants: the `n750` anchor and the translation at both points | 9 | at most one above 34 of 500, none above 38 |
| | simulation: the anchor and the translation | 6 | at most one above 34, none above 37 |
| S2 protected power | the target statistic of every plant: 8 simulated, 6 held-out, 16 in-sample | 30 | at most two loss p-values below 0.025, none below 0.05 / 30 |
| S3 specificity | the two off-target statistics of the simulated plants | 16 | at most two excess p-values below 0.025, none below 0.05 / 16 |
| | the two off-target statistics of the in-sample plants | 32 | at most two below 0.025, none below 0.05 / 32 |
| G gain | the target statistic of the simulated plants | 8 | at least one gain p-value at most 0.05 / 8 |

S1 is judged for `pls_cv` and for `predictive_orthonormal`. S2, S3 and G compare the candidate with `pls_cv` on the same units by exact tests on the discordant ones, through the pilot's functions imported unchanged.

- **Why the cohort anchor is not judged again.** Its units are item 21's. It enters through the reproduction.
- **Why S3 is a paired excess.** Production itself may leak: in the probe it rejected an in-sample rotation of 90° by `shape` in 6 and 8 of 16 units. An absolute bound would then close the candidate where it leaks less. The pilot judged its band the same way. Production's absolute rates are graded in its own statement (D15).
- **Why S3 reads the in-sample plants and not the held-out ones.** D5: only the first have an exact off-target truth.
- **Why S2 reads all three kinds of plant.** A loss is a loss wherever the plant lies, so one family holds the simulated, the held-out and the in-sample cells.
- **Why G reads the simulation alone.** The roadmap's condition names a gain at the design in simulation, where the plant is the population's. An in-sample plant rewards the space that leans most on the in-sample fit, and a held-out one is mostly reference noise. Gains on both are judged by the same rule and recorded beside G.
- **What the rules tolerate.** At exact size a correctly sized arm fails the plants size family with probability at most 0.07 and the simulated one 0.06. The paired families carry the pilot's levels.
- **Unavailable readings** fail their family, as in `judge_rates`.
- **Below grade.** A decision with fewer than 500 replicates on a control, 200 on a plant or 999 permutations is marked `indicative`.

### D13. The verdict

The first row that applies:

| verdict | when | consequence, fixed with rule version 1 |
|---|---|---|
| `gate_failed` | a pre-flight is missing or failed, the anchor reproduction failed, or a unit is inconsistent | No reading is interpreted. The cause is fixed and the affected units are rerun under the same rule. |
| `production_excess` | `pls_cv` fails S1 | Production's unrestricted null with `df` exceeds size on a control at the named readings. No reading of `predictive_orthonormal` is interpreted. A diagnosis of that control comes before any other reading on ADNI. If the control is a translation on ADNI's blocks, the roadmap states that a group main effect along the stage program reaches the interaction statistics at ADNI's design. |
| `size_excess` | production holds S1; the candidate fails it | `predictive_orthonormal` exceeds size on a control this change adds, where production holds. Item 21's hold on the global null does not extend to that control. Open item 1's second reopening condition is not met for ADNI, and it stays an opt-in diagnostic. |
| `power_loss` | the candidate fails S2 | `predictive_orthonormal` loses target power against production at the named plants. The condition is not met for ADNI, and it stays an opt-in diagnostic. The measurement-visibility page records which differences its two measured dimensions give up at ADNI's design. |
| `specificity_excess` | the candidate fails S3 | `predictive_orthonormal` adds off-target rejections over production at the named cells, where the kind of the planted difference is known exactly. The condition is not met for ADNI, and it stays an opt-in diagnostic. |
| `no_gain` | the candidate fails G | No simulated target cell at ADNI's table gains at the family threshold, so the measured dimension buys nothing at ADNI's design in simulation. The condition is not met, and it stays an opt-in diagnostic. A gain recorded on ADNI's noise without a simulated one is a lead on the simulation of ADNI's design, not a qualification. |
| `reopened` | none of the above | At ADNI's design `predictive_orthonormal` holds size on the added controls, loses no target power, adds no off-target rejection where the truth is exact, and gains over production in simulation. With item 21's size reading, open item 1's second reopening condition is met for ADNI's export at 750 participants and the full cohort. It stays an opt-in diagnostic: a separate change may propose its paper-grade run with a rank rule for a fixed measured dimension, and only that run can change a default. |

- **The standing sentence,** appended to every consequence: "In every verdict the SEA-AD anchor's reading stands for SEA-AD's design, `keep_production` stands, and the out-of-span cells are reported beside the verdict and were not judged."
- **Why production gates the candidate.** As in item 21: a reading of the candidate means something only where the reference holds.
- **Why this order.** It is the pilot's (`size_excess` → `power_loss` → `specificity_excess` → `no_gain`), so the verdict names the most basic property that fails. All failing criteria are listed with it.

### D14. What the verdict licenses

Written into the decision with the rule version, beside item 21's S1–S5.

- **L1. The reading speaks for ADNI's export, the two design points and the simulated table.** Nothing is extended to 80, 200 or 500 participants, where only size was read.
- **L2. At the full cohort, power and specificity were read on in-sample plants and in simulation.** No held-out plant exists there. The held-out readings are those of 750 participants.
- **L3. `reopened` permits a proposal.** It changes no default and adds no option. The later change carries the rank rule open item 1 asks for.
- **L4. Cohorts do not vouch for one another.** The SEA-AD veto stands for SEA-AD's design, and ROSMAP has no reading.
- **L5. The out-of-span cells license nothing.** Outside the stage span the split test's complement is the detector, and its reading on ADNI's noise is not this change's.

### D15. Production's statement

Read on `pls_cv` alone on ADNI's blocks. It is written whatever the verdict, except under `gate_failed`, and never changes it.

- **Power.** Per design point, reference and target statistic: `informative_from`, the smallest effect whose rejection rate reaches 0.80 (item 9's floor), or null with the largest rate and its effect. One sentence per statistic in the form item 9 gave `adni_sexdiff`, naming the reference: a held-out plant bounds the detectable effect from above, an in-sample one from below.
- **Off-target.** Production's 32 off-target readings on the in-sample plants as one family under `judge_rates` at R = 200 (at most two above 16 of 200, none above 20): `specific` when it holds, else `leaks` with the cells and statistics.
- **Why a rule.** A table of production's rates with no rule would be read after the fact. The off-target family's false-alarm rate is at most 0.08 at exact size.
- **What it is for.** The case study reads its non-rejections on a simulation of its design. This is the reading on its own noise, and it is production's, so it holds under every verdict on the candidate.

### D16. Descriptive readings

None changes the verdict or production's statement.

- **The bracket:** each target rate on `sample` beside `heldout` at the same mode and effect, per arm.
- **Recorded gains:** the gain rule on the held-out (6) and on the in-sample (16) target readings.
- **Held-out off-target rates,** beside the reference's reliability.
- **The out-of-span cells** (D7).
- **The within-stage arms:** their control counts under S1's families, and their paired comparisons with the unrestricted arms on every plant.
- **Rates by rank class:** among units whose CV rank is 2 and among those above it, where the two coordinate systems differ.
- **Capture,** per cell: the share of the plant's energy and of its reference configuration inside the fitted span and inside the predictive subspace, through `LinearLatentMap`. The probe shows why the share is reported and not judged: the predictive subspace keeps 0.13 of a held-out configuration while reading its stretch about as often as the scores, because it follows the directions that discriminate stages and not those along which their means differ.
- **Reliability of the held-out configuration,** per `n750` anchor unit: the correlation of the centred adjusted configuration between two halves of the reference, split within stage from the unit's `reference_split` seed, and the two leading eigenvalues of the halves' cross-product.
- **The fit under a plant:** CV rank distribution and median held-out AUROC per cell.
- **In simulation:** the truth's intended geometry per cell beside the observed statistics' medians.

### D17. Execution, compute and data handling

- **Subcommands.** `preflight`, `run` and `merge` take `--profile`; `decide` reads both merged profiles. `--export-dir` is required wherever the plants profile is read and has no default. Each profile writes `settings.json` bound by hash under `<out>/<profile>/`; a resume, merge or decision under other settings exits naming the keys. Round-robin shards, resume by complete unit key, `--workers` with BLAS pinned to one thread in the shell.
- **Records.** One JSONL line per unit: seeds, redraw counts, the cell, ranks, AUROC, per arm the statistics, p-values, `angle` null 95th percentile and seconds, the plant's norm per stage, its capture shares and the scale drift, the checks and timing. No participant identifier, label vector, subsample membership, fold assignment or shift vector. A test searches a stand-in unit's record for every identifier of the stand-in export and for any field with one entry per participant or per feature.
- **Committed.** Under `results/adni-gain-specificity-<date>/`: per profile the settings, the pre-flight record, the rate, paired and cell tables and `accounting.json`; the decision; `PROVENANCE.txt`. Records stay off the repository.
- **Compute, before the smokes.**

  | part | units | per unit | core-hours |
  |---|---|---|---|
  | plants, `n750` | 4,200 | about 19 s (CV 17 s) | 22 |
  | plants, `cohort` anchor | 500 | about 5 s (forced) | 1 |
  | plants, `cohort` other cells | 2,500 | about 24 s (CV 22 s) | 17 |
  | simulation, on `ing` | 2,600 | about 95 s (CV 90 s) | 69 |

  About 40 core-hours on the workstation, six hours at seven workers, and about 70 on the cluster.
- **Where it runs.** The plants profile on the workstation, where `adni_sexdiff` made the export, from a clean worktree at a pushed revision, launched as a `systemd-run --user` unit with its scripts and output outside `/tmp` (item 21's first launch died with a service restart). The simulation profile on `ing` as an sbatch array, with a dry-run variable that prints the command.

## Risks / Trade-offs

- **[An in-sample plant rewards the stage-supervised fit]** → It never decides the gain. Its target rates are reported as the upper bracket beside the held-out ones.
- **[A held-out plant is mostly the reference's estimation noise]** → Its role is a direction no fit was built on, not the cohort's program. Its reliability is reported with every rate, and none of its off-target rates is judged.
- **[The cohort's own design has no held-out plant]** → L2 says so. The cohort is read on in-sample plants and in simulation; the held-out readings are at 750.
- **[Two subsamples of 750 share about half their participants, and every rate is conditional on ADNI's 1,501]** → Stated with every rate. Independent real noise needs ROSMAP.
- **[The simulation is a structural analogue]** → Its CV rank is 19 against ADNI's 8, its blocks are InterSIM's, and its stage program is weaker in pooled-SD units. G follows the roadmap's condition; the recorded real-noise gains say whether the analogue agrees.
- **[A paired excess passes a candidate that leaks wherever production leaks]** → Both arms' absolute off-target rates are reported, and production's are graded in its statement.
- **[The out-of-span cells are not judged, and the probe saw the candidate read them as `shape`]** → They are reported with paired counts in the report and the visibility page, the standing sentence names them, and D7 records the stricter rule.
- **[Probes read rates of both arms before the rule was written]** → They ran on scratch seeds at 16 replicates and 199 permutations, are recorded in `notes.md` as not results, and the run uses the anchor's and a fresh seed family. The rule's families were fixed by what each plant's truth supports (D5, D7), and D7 states the one place where the probe bears on a choice.
- **[The plant changes the stage fit, so a cell's rank is not the anchor's]** → It is part of the measurement. The rank and AUROC are reported per cell.
- **[A rotation and a bend are exact only to the order of the scale drift]** → About 1.5 % of a pooled SD on the most affected feature at effects up to 1.0 in the probe. The drift is recorded per unit, and both arms read the same data.
- **[200 units per plant]** → The paired tests carry the comparison. A decision below grade is marked `indicative`.
- **[The covariate adjustment of a small reference is unstable]** → The smallest reference stage holds 141 participants for six covariate columns. The fit is checked for full column rank, and a failure fails the unit.
- **[A six-hour run on a shared workstation]** → Seven workers with BLAS pinned, a detached unit, resume by complete unit, and a cap approved by the user.
- **[A participant-level value reaches a record or a committed file]** → A plant is a function of participants' data, so no shift vector or configuration is recorded, only norms and shares. The test searches for identifiers and per-participant or per-feature fields, and the commit task checks every committed file.
- **[`reopened` is read as adoption]** → Every consequence says the option stays an opt-in diagnostic, and L3 names what is still missing.

## Migration Plan

Nothing migrates: no package code, configuration, fixture or committed record changes.

- **On `reopened`:** a separate change proposes a rank rule for a fixed measured dimension and the paper-grade run of `predictive_orthonormal` at ADNI's design.
- **On `no_gain`, `specificity_excess`, `power_loss` or `size_excess`:** open item 1 closes for ADNI. Open item 2 continues with ROSMAP for production alone.
- **On `production_excess`:** a diagnosis change for the named control, before any other reading on ADNI.
- **Under any verdict:** `adni_sexdiff` may cite production's statement in its own change, and the out-of-span cells are the starting point of the split test's reading on ADNI's noise.
- **Rollback:** remove the driver, its sbatch array and its test. Nothing imports them.

## Open Questions

- The worker count and core-hour cap of each run are set from its smoke's timing at the approval gates.
