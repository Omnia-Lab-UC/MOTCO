# Design

## Context

See `proposal.md` (Why) for the motivation. The facts below shape the approach.

- **The export** (`adni_sexdiff/results/motco_input`, written by its `export` step):
  - `metadata.csv` (`rid`, `sex`, `stage`) and `covariates.csv` (`rid`, `age`, `education`, `apoe4`, `phase`, `statin`);
  - `p180.csv` (133 metabolites) and `nmr.csv` (165), participants × features with `rid` first, on the scale MOTCO reads as `identity`;
  - `manifest.json`, which records each output's sha256 with a `participant_level` flag and the counts: 1,501 participants, Female 237 / 310 / 125 and Male 223 / 448 / 158 at CN / MCI / AD, 298 joint features.

  Its layout is not the SEA-AD export's, so `plasmode_check.load_export` cannot read it.
- **The case study's configuration** (`adni_sexdiff`, pinned to MOTCO 0.10.0): the five covariates with `apoe4`, `phase` and `statin` categorical (six design columns), rank floor 2, net displacement, the harness CV defaults, rejection at p ≤ α, and no residual correction, which did not exist yet. Its primary selected rank 8 at a held-out stage AUROC of 0.5995.
- **What the SEA-AD anchor read** (rejections of 500, `results/within-stage-permutation-2026-11-16/plasmode/rate_table.csv`), on a rank-5 fit with in-sample stage R² 0.840, held-out AUROC 0.636 and 125 features per donor:

  | arm | `delta` | `angle` | `shape` |
  |---|---|---|---|
  | `pls_cv` | 21 | 24 | 13 |
  | `predictive_orthonormal` | 67 | 9 | 38 |
  | `pls_cv__within_stage` | 19 | 27 | 19 |
  | `predictive_orthonormal__within_stage` | 29 | 17 | 32 |

- **What the probe read on ADNI** (`notes.md`; label-free, a handful of subsamples per n):

  | | n = 80 | n = 200 | n = 500 | n = 1,501 |
  |---|---|---|---|---|
  | features per participant | 3.7 | 1.5 | 0.6 | 0.2 |
  | CV ranks selected | 2 in 10 of 16; up to 8 | 2 to 13 | 2 to 7 | 8 |
  | median held-out AUROC | 0.558 | 0.553 | 0.559 | 0.600 |
  | median in-sample stage R² at rank 8 | 0.614 | 0.345 | 0.162 | 0.086 |
  | median in-sample stage R² at rank 19 | 0.933 | 0.620 | 0.314 | 0.157 |
  | production CV, one core | 7.9 s | 9.2 s | 12.3 s | 22.3 s |

- **The options exist.** `analyze_trajectories` takes `latent_coordinates="predictive_orthonormal"` (0.24.0), `permutation_strata="stage"` (0.25.0), `residual_correction="df"` (0.13.0) and `forced_components`, and they combine. A forced-rank call at 999 permutations takes 0.34 to 0.49 s on these sample sizes.
- **`df` at these designs.** The reduced model has ten columns (intercept, group, two stage indicators, six covariate columns), so κ is 1.069 / 1.026 / 1.010 / 1.003 at n = 80 / 200 / 500 / 1,501.
- **The size rule.** `judge_rates` at R = 500: a reading *exceeds* above 34 rejections. A family of m readings holds with at most one exceedance (two at m = 21) and none above the ceiling: 36 / 37 / 38 / 39 / 40 rejections at m = 3 / 6 / 9 / 12 / 21.
- **Data rules.** ADNI is controlled-access. `adni_sexdiff` commits no participant-level file and forbids copying its MOTCO inputs between machines. MOTCO's repository is public.
- **Helpers written for SEA-AD sort stage labels as strings.** `predictive_subspace_check.stage_r2` and the size-excess diagnosis's out-of-fold helpers build their one-hot columns from sorted labels. SEA-AD's are `0` / `1` / `2`. ADNI's sort as AD, CN, MCI, which is not the fitted model's column order; the probe's first pass returned negative R² values this way.

## Goals / Non-Goals

**Goals:**
- A size reading of the production scores and of `predictive_orthonormal`, under both nulls, on a second cohort's real noise and at design points whose fits are independent, under a rule fixed before any unit runs.
- A reading of item 19's guess that does not depend on which rank the CV happens to select.
- The scoping rule the roadmap requires, declared with the decision.
- Functions a later cohort's driver can import unchanged.

**Non-Goals:**
- Any power or specificity reading. No feature is moved; there is no plant.
- Ranking measurement spaces.
- Explaining a failure. A failing verdict names its lead and stops.
- Editing a committed driver. Parents are imported and bound, never changed.

## Decisions

### D1. A driver on existing options; the package is not touched

Every arm is one `analyze_trajectories` call. `motco` stays at 0.25.0.

- **Why.** The options the arms need were added and pinned by the changes that validated them. Going through the public entry point reads the code a case study would run.
- **Alternative: a shared-fit path in the driver** (fit once, measure many times through the stats primitives). It would save the forced refits, which cost under half a second each, and would add an inference path that needs its own identity check. Rejected.

### D2. The contract is a value, and ADNI's is the first

A `CohortContract` names the identifier, group and stage columns, the two groups, the stage order, the blocks with their files, the covariates with their kinds and levels, and the rank and AUROC the cohort's fit must reproduce (D5). `load_cohort(directory, contract)` returns the blocks and one metadata frame in the metadata file's row order. The subsampling, the unit builder and the decision take the contract as an argument.

- **Checks, in order, each raising a named `plasmode_check.ExportContractError` subclass:**
  - the manifest exists, is the `export` step's, records `motco_transform: identity`, and lists the four participant-level files with their sha256;
  - each of those files exists and hashes to the manifest's value;
  - the metadata has the three columns, no repeated participant, exactly the two declared groups, no stage outside CN / MCI / AD, no empty cell, and the manifest's participant count and sex × stage table;
  - each block lists the metadata's participants in the same order, has the manifest's feature count and holds finite numbers;
  - the covariates list the same participants, with no missing value and no level outside the case study's (`apoe4` 0 / 1 / 2, `phase` ADNI1 / ADNI2GO, `statin` 0 / 1).
- **Binding.** The sha256 of the five files read goes into `settings.json`, so a resume on another export is refused by name. The pre-flight (D5) ties the export to the case study's fit.
- **Why a value.** ROSMAP is the next cohort of open item 2. With the contract as an argument its driver imports this one's functions unchanged, as the driver convention asks, instead of copying them around other constants.
- **Alternative: adapt `plasmode_check.load_export`.** It is a committed parent with SEA-AD's blocks, covariates and manifest layout as module constants. Editing it would break its children's reproduction.

### D3. Four design points, at the cohort's proportions

| point | participants | Female CN / MCI / AD | Male CN / MCI / AD | features per participant |
|---|---|---|---|---|
| `n80` | 80 | 13 / 16 / 7 | 12 / 24 / 8 | 3.7 |
| `n200` | 200 | 31 / 41 / 17 | 30 / 60 / 21 | 1.5 |
| `n500` | 500 | 79 / 103 / 42 | 74 / 149 / 53 | 0.6 |
| `cohort` | 1,501 | 237 / 310 / 125 | 223 / 448 / 158 | 0.2 |

500 replicates at 999 permutations per point; 2,000 units.

- **Subsampling.** Within each sex × stage cell, without replacement, by the cell's own draw from the unit's subsample seed. The table is the cohort's scaled to n and rounded by largest remainder, ties broken in the table's order.
- **Why the cohort's proportions.** The question is ADNI's design at fewer participants. Its largest stage is the middle one, where SEA-AD's is the last, and item 20 showed that the unrestricted null's error follows how the largest stage's dispersion compares with the others'.
- **Alternative: SEA-AD's table on ADNI's noise at n = 80.** It would separate the cell table from the noise source. It is a second question and waits for a reading that asks it.
- **Why the points are not nested.** Independent subsamples make each point its own reading. Nesting would pair the points and correlate the readings one family is judged on.
- **Why the cohort point is in.** It is the case study's design and the only point where the CV keeps eight components. Its 500 units share one fit and differ in the labels, so its rates are conditional on that fit, as SEA-AD's anchor rates are.
- **An absent covariate level.** A subsample with no participant at some level of a categorical covariate would lose a design column. It is refused and redrawn from the next attempt's seed, up to 100 attempts, and the redraw is counted. At n = 80 the probability is about 2 in 10,000 per draw, from `apoe4` = 2; above it is negligible.
- **Sampling fractions.** The smallest cell (Female AD, 125) gives 7, 17 and 42 participants: 6 %, 14 % and 34 %. Subsamples at n = 500 overlap (Risks).

### D4. The null-making step and the seeds

- **Labels.** `plasmode_check.permute_within_stage`, imported unchanged, on the unit's sample: the sex labels are permuted among the participants of each stage, stage by stage in declared order. A draw equal to the sample's real labels is refused and redrawn from the next attempt's seed.
- **Why within stage.** It is the plasmode's step. It keeps the cell table and each stage's own noise, and it makes the group effect and the interaction null while leaving the stage effect in place. Whatever sex difference the cohort holds stays in the data as within-stage noise.
- **Seeds.** `sha256([driver name, 2039, point, replicate, purpose(, attempt)])`, truncated as `plasmode_check._seed` does. Purposes: `subsample`, `labels`, `rrpp`, `folds`. Base seed 2039 is used by no committed driver. Every arm of a unit takes the unit's `rrpp` seed, so arms under one null share their draws.
- **The probe's seeds** were scratch seeds of another form and are unrelated.

### D5. One cross-validation per subsampled unit; the cohort's fit pinned by the pre-flight

- **Subsampled points.** `pls_cv` runs the production cross-validation once. Every other operating arm is forced at the rank it selects, so its final fit is the CV's final fit.
- **The cohort point.** The fit does not see the group labels, so all 500 units fit the same model. `preflight` runs the cross-validated analysis on replicates 0 and 1 and requires:
  - rank 8 at a mean held-out AUROC of 0.599514 to six decimals, the case study's primary record;
  - statistics and p-values equal to the analysis forced at rank 8, exactly.

  Every cohort unit then runs forced at rank 8. `run` refuses to start without a passed `preflight.json` under the same settings.
- **Why pin to the case study.** It ties the export, the row order and the environment to the analysis whose design this is. The probe already reproduced both values under 0.25.0.
- **Alternative: cross-validate every cohort unit.** About 3 core-hours to recompute one label-free integer.

### D6. Thirteen arms on a subsampled unit, nine on a cohort unit

| arm | fit rank | coordinates | draws | correction | family |
|---|---|---|---|---|---|
| `pls_cv` | CV | scores | unrestricted | `df` | operating, reference |
| `predictive_orthonormal` | CV | predictive, orthonormal | unrestricted | `df` | operating, candidate |
| `pls_cv__within_stage`, `predictive_orthonormal__within_stage` | CV | as named | within stage | `df` | operating |
| `pls_rank8`, `pls_rank19` | 8, 19 | scores | unrestricted | `df` | fixed rank, reference |
| `predictive_orthonormal_rank8`, `predictive_orthonormal_rank19` | 8, 19 | predictive, orthonormal | unrestricted | `df` | fixed rank, candidate |
| the four above with `__within_stage` | 8, 19 | as named | within stage | `df` | fixed rank |
| `pls_cv__uncorrected` | CV | scores | unrestricted | none | control |

At the cohort point the CV rank is 8, so the rank-8 arms are the operating arms and are not run twice.

- **Why `df` on every decisive arm.** It is the validated configuration at six covariate columns and 80 samples (item 8), and at 1,501 it changes the residual scale by 0.3 %. One configuration across the ladder keeps the points comparable.
- **Why the control arm.** At the cohort it is the configuration the case study was issued under. Below it, the paired comparison with `pls_cv` reads what `df` does on real noise, which bears on roadmap open item 9. In simulation the uncorrected `delta` read 0.068 and 0.064 at the two anchors of SEA-AD's configuration (item 7 (c)), close to the bound, so the arm is not a test of this anchor's sensitivity.
- **Why fixed ranks.** At a weak stage signal the CV mostly returns the floor, and a fit at the floor is barely fitted in sample (median R² 0.16 at n = 80). SEA-AD's excess arose on a fit at 0.84. The fixed ranks produce such fits on ADNI's noise, with the held-out signal unchanged.
- **Why 8 and 19.** Both are ADNI's and neither is tuned. Rank 8 is what the CV selects on the cohort. Rank 19 is the highest the CV can select and the rank the ADNI-design simulation selected (item 9). Together they span in-sample stage fits from about 0.16 to 0.93 across the points, on both sides of SEA-AD's 0.84.
- **Alternative: rank 13**, whose fit at n = 80 (0.81) is the nearest to SEA-AD's. It would be a rank chosen from the probe to match.
- **Alternative: the other two predictive coordinate systems.** On SEA-AD the three agreed to within three rejections per statistic. They would add eight arms and no decision.
- **Stage labels handed to imported helpers.** Every helper that sorts labels receives each row's stage position in the declared order, as a string. A test pins the in-sample stage R² against a direct computation on labels that do not sort in declared order.

### D7. Checks on every unit

- `rank`: every forced arm reports the rank it was forced at.
- `same_fit`: per fit rank, every arm's fitted model equals the first one fitted at that rank (rotations and response loadings to 1e-8). For the operating arms that model is the CV's.
- `measured_dimension`: each predictive arm measures min(rank, 2) dimensions.
- `same_statistics`: each within-stage arm's observed statistics equal its unrestricted counterpart's exactly, and `pls_cv__uncorrected`'s equal `pls_cv`'s.
- `cell_table`: the realized group × stage table is the point's.

A unit failing a check is `inconsistent`, names the check, the arm and the deviation, and leaves every rate. The pilot's check functions are used where their signatures allow; otherwise the check is restated here and a test pins it against the pilot's function on a stand-in unit.

### D8. The size rule and the size scope

- **Rejection.** p ≤ 0.05, the plasmode's and the case study's rule, exact at 999 permutations.
- **Families.** An arm is judged on a set of design points by `judge_rates` over the set's `delta`, `angle` and `shape` readings as one family:

  | set | readings | holds with |
  |---|---|---|
  | all four points | 12 | at most one above 34 of 500, none above 39 |
  | `n200`, `n500`, `cohort` | 9 | at most one above 34, none above 38 |
  | `n500`, `cohort` | 6 | at most one above 34, none above 37 |
  | `cohort` | 3 | at most one above 34, none above 36 |

- **Size scope.** The first of those sets, in that order, on which the arm holds. It is recorded as the point the set starts at; an arm whose cohort readings fail has no scope.
- **Why upper sets that contain the cohort.** The cohort is the design a case study has. A scope that skipped a point, or stopped short of the cohort, would be a set chosen after the reading.
- **Why one family per set and not one per point.** At exact size a point's three readings fail their own rule with probability 0.038, so an arm would fail at least one of four points with probability 0.14. The twelve readings as one family fail with probability 0.07 (independent readings assumed). Two arms decide the verdict.
- **What the rule tolerates.** One reading of up to 39 rejections (0.078) among twelve. SEA-AD's excess was 67.
- **Unavailable readings.** A reading with no available unit fails its set, as in `judge_rates`.

### D9. The verdict

Read on `pls_cv` and `predictive_orthonormal`, the operating arms under the unrestricted draws. The first row that applies:

| verdict | when | consequence, fixed with rule version 1 |
|---|---|---|
| `gate_failed` | the pre-flight is missing or failed, or a unit is inconsistent | No reading is interpreted. The cause is fixed and the affected units are rerun under the same rule. |
| `production_excess` | `pls_cv`'s scope does not start at `n80` | Production's unrestricted null with `df` exceeds size on ADNI's noise at the named readings. No reading of `predictive_orthonormal` is interpreted at a point where production fails. A diagnosis of those points comes before any other reading on ADNI. If the cohort is among them, the roadmap states that a rejection at ADNI's design lacks a size-controlled null on its own noise; the case study's primary non-rejections stand, since an excess cannot produce them. |
| `excess_persists` | production holds; `predictive_orthonormal` has no scope | The excess appears at ADNI's own design, at 0.2 features per participant, so item 19's prediction for ADNI is refuted. The veto rests on two cohorts and the second reopening condition of open item 1 is closed for ADNI. `predictive_orthonormal` stays an opt-in diagnostic. |
| `excess_scoped` | production holds; the candidate's scope starts at `n200`, `n500` or `cohort` | `predictive_orthonormal` holds size on ADNI's blocks from the named point up to the cohort and exceeds below it. The size part of open item 1's second condition is met for ADNI inside that scope only. It stays an opt-in diagnostic; the gain and specificity readings at ADNI's design are the next change and are confined to the scope. |
| `excess_absent` | both hold on all four points | `predictive_orthonormal` holds size on ADNI's blocks at every design point under RRPP's unrestricted draws: SEA-AD's excess does not transfer to ADNI. The size part of open item 1's second condition is met for ADNI. It stays an opt-in diagnostic until the gain and specificity readings at ADNI's design, the next change. The fixed-rank grade says whether the hold is owed to ADNI's noise or to the ranks its CV selects. |

- **Recorded beside the verdict, never changing it:** the scopes of the two within-stage operating arms. A within-stage arm that does not hold on all four points is named with its exceeding readings, and the documentation of the opt-in null records it.
- **Why production gates the candidate.** A reading of the candidate means something only where the reference holds. Item 17's sparse arms were stopped at a rate production also read.
- **Why the verdict reads the unrestricted draws.** They are the default and the null open item 1's second condition is about. The within-stage draws were read `specificity_excess` in item 20 and can open nothing here.
- **In every verdict** the SEA-AD anchor's reading stands for SEA-AD's design, and `keep_production` stands.
- **Below paper grade.** A decision on fewer than 500 replicates or 999 permutations at any point is marked `indicative`.

### D10. The fixed-rank reading

One family of 21 readings per arm: three statistics at ranks 8 and 19 on `n80`, `n200` and `n500`, and at rank 19 on `cohort`. It holds with at most two readings above 34 of 500 and none above 40. The first row that applies:

| grade | when | statement, fixed with rule version 1 |
|---|---|---|
| `reference_exceeds` | the scores under unrestricted draws fail the family | An over-fitted rank mis-sizes the scores themselves on ADNI's noise. The predictive arms are reported without a grade, and the roadmap records that production's size depends on a rank its CV does not select here. |
| `absent` | `predictive_orthonormal` under unrestricted draws holds it | No excess at either fixed rank. An over-fitted stage fit is not sufficient on ADNI's noise: SEA-AD's excess needs a property of its blocks this ladder lacks. The lead is a feature-subsampling ladder on SEA-AD's own export, then ROSMAP. |
| `reproduced` | it fails, and holds under within-stage draws | ADNI's noise gives the predictive subspace the excess once the fit is over-fitted, and the within-stage draws remove it. Item 19's mechanism is not specific to SEA-AD's blocks. The report tabulates the excess against the in-sample stage fit. |
| `reproduced_unrestored` | it fails under both nulls | The excess appears at an over-fitted rank and the within-stage draws do not remove it at the named readings. Unequal stage dispersion is not the whole mechanism on ADNI's noise; the dispersion table and the out-of-fold shares are the leads. |

- **Why a grade and not a table to interpret.** A reading with no rule would be explained after the fact. The family rule keeps its false-alarm rate at 0.05 at exact size.
- **Why it does not change the verdict.** No case study runs a forced rank. The grade says what the operating reading is owed to.
- **In the family or not.** The cohort's rank-8 readings are the operating arms' and are counted in the verdict only.

### D11. The rule that scopes a size reading

Declared here because the roadmap requires it before the reading. It is written into the decision file with the rule version.

- **S1. A size reading speaks for its export and its design points.** The points are declared before the run. Nothing is interpolated between them or extended beyond them.
- **S2. An option's scope on a cohort is its size scope** (D8). An option with no scope there is not admissible on that cohort.
- **S3. Cohorts do not vouch for one another.** A hold on one export does not lift a failure on another and does not extend to a cohort without its own reading. The SEA-AD anchor's reading stands for SEA-AD's design.
- **S4. Size is one of three readings.** Reopening an option at a design also needs a gain over production at that design in simulation and specificity there, each under a decision declared before its reading. Reopening means that a separate change may propose its paper-grade run; only that run can change a default.
- **S5. Covariates describe; they do not extend.** Features per participant, the selected ranks and the in-sample stage fit are reported beside a scope. They license no reading at a design that was not read.

- **Why scope by export and point, not by a covariate.** Two cohorts differ in omic type, feature count, noise and batch structure at once. A rule of the form "admissible below so many features per participant" would claim a transfer that the roadmap says no reading supports.
- **Alternative: no scoping, one failure anywhere vetoes everywhere.** It is the rule in force, and it stays in force for SEA-AD's design. The roadmap asks for the scoped form because ADNI and ROSMAP are the cohorts MOTCO was meant for.

### D12. Descriptive readings

None changes the verdict or the grade.

- **Rates by rank class.** Each operating arm's rates among units whose CV rank is 2 and among those above it, with their counts. At rank 2 the candidate measures the scores' own subspace.
- **The control arm** against `pls_cv`: rates and the exact paired test.
- **Paired comparisons:** each within-stage arm against its unrestricted counterpart, and each predictive arm against the scores at the same rank and null, through `production_latent_space_check`'s exact tests.
- **The point table:** units by CV rank, the median held-out AUROC, and per fit rank (CV, 8, 19) the median in-sample stage R² beside the median out-of-fold explained share (four stage-stratified folds from the unit's `folds` seed, `predictive_size_excess_check`'s helpers); per coordinate system and fit rank the median residual dispersion of each stage over the pooled one; each arm's median `angle` null 95th percentile; the redraw counts.
- **The cross-cohort table:** SEA-AD's committed anchor counts on its reproduced and its fresh label draws (the within-stage run's `rate_table.csv`, cells `anchor` and `anchor_fresh`, bound by sha256) beside each ADNI point's, with features per participant, fit rank and in-sample stage R².

### D13. Execution, compute and data handling

- **Subcommands.** `preflight`, `run`, `merge`, `decide`, with `--export-dir` required and no default. `settings.json` is bound by hash and carries the export's file hashes and the parent table's; a resume, merge or decision under other settings exits naming the keys. Round-robin shards, resume by complete unit key, `--workers` with BLAS pinned to one thread in the shell.
- **Records.** One JSONL line per unit: seeds, redraw counts, the cell table, ranks, AUROC, the fit readings, per arm the statistics, p-values, `angle` null 95th percentile and seconds, the checks and timing. No participant identifier, label vector, subsample membership or fold assignment. A test searches a stand-in unit's record for every identifier of the stand-in export.
- **Committed.** Under `results/adni-real-noise-anchor-<date>/`: `settings.json`, `preflight.json`, the rate, paired and point tables, `accounting.json`, the decision and `PROVENANCE.txt`. Records stay off the repository.
- **Compute, before the smoke.** A subsampled unit is one CV plus twelve forced calls and twelve small out-of-fold fits: about 13, 15 and 18 s at the three points. A cohort unit is nine forced calls, about 5 s. About 7 core-hours for the 2,000 units, about an hour at seven workers.
- **Where it runs.** On the workstation, where `adni_sexdiff` made the export. The export is not copied, so there is no sbatch array. The run is launched detached, with its scripts and output outside `/tmp`.

## Risks / Trade-offs

- **[Subsamples of one cohort are not independent cohorts]** → Every rate is conditional on ADNI's 1,501 participants, and at `n500` two subsamples share about a third of their participants. The report says so; independent real noise needs ROSMAP.
- **[The cohort's 500 units share one fit]** → Stated with every cohort rate, as for SEA-AD's anchor. `n500` is the largest point whose fits vary.
- **[At small n the CV rank is mostly 2, where the candidate measures the scores' subspace]** → The fixed-rank arms read the over-fitted regime, and the rates by rank class are reported beside each scope.
- **[A correctly sized arm fails the ladder family 7 % of the time, and two arms decide]** → It is the stage-span rule's own level. A verdict resting on one reading within a few counts of the ceiling is reported with its counts; the rule is not relaxed, and a rerun on fresh seeds is a new change.
- **[The within-stage arms hold because the anchor's labels are drawn the way their null assumes]** → They decide nothing here. Their scopes are recorded, and the evidence for that null stays item 20's simulation and dispersion profiles.
- **[The fixed ranks are not what a case study runs]** → They grade a mechanism, not an operating characteristic, and cannot change the verdict.
- **[The ladder's most over-fitted point is not SEA-AD's]** → At most 3.7 features per participant against 125, on metabolites instead of RNA and ATAC. An `absent` grade says an over-fitted fit is not sufficient here; it does not say what SEA-AD's blocks add, and its statement names the next reading.
- **[`df` over-corrects by up to 1.9 % under within-stage draws]** → Conservative, and quoted from the within-stage report.
- **[The p180 run offset]** → The scores follow `phase`, which the trajectory model adjusts for. Under permuted labels it is noise whose mix differs by stage; the dispersion table shows how unequal the stages are in each coordinate system.
- **[An imported helper mis-orders ADNI's stage labels]** → Stage positions are passed everywhere, and a test pins the R² on labels that sort out of order.
- **[The environment drifts from the case study's]** → The pre-flight stops the run. A mismatch is diagnosed, never waived.
- **[A participant-level value reaches a record or a committed file]** → Records hold statistics only, the test searches for identifiers, and the commit task checks every committed file.
- **[A hold is read as reopening the predictive subspace]** → Every consequence statement says it stays an opt-in diagnostic, and S4 names the two readings still missing.

## Migration Plan

Nothing migrates: no package code, configuration, fixture or committed record changes.

- **On `excess_absent` or `excess_scoped`:** a separate change reads the gain and the specificity at ADNI's design, in simulation and on plants in a held-out reference span (open item 2's next slice), confined to the scope.
- **On `excess_persists`:** no change follows for the predictive subspace on ADNI. Open item 2 continues with the plants for production alone.
- **On `production_excess`:** a diagnosis change for the named points, before any other reading on ADNI.
- **Rollback:** remove the driver and its test. Nothing imports them.

## Open Questions

- The worker count and the core-hour cap of the full run are set from the smoke's timing at the approval gate.
