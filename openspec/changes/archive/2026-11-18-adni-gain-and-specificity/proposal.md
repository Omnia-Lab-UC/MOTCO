# Proposal

## Why

The ADNI anchor (roadmap item 21, `excess_absent`) read one of the three things that reopen `predictive_orthonormal` at a design. Its decision fixed the rule: a size reading on a cohort's noise opens nothing without a gain over production at that design in simulation and specificity there, each under a decision declared before its reading. Size is met for ADNI. The rest is unread, and so is most of what production itself does on a second cohort:

- **No real-noise specificity reading exists.** Off-target rates have been judged in simulation only. On real noise the decisions have read the global null and one plant's target statistic.
- **Real-noise power rests on one plant on one export.** SEA-AD's in-span plant is a magnitude edit along the in-sample stage means, which rewards any space fitted on those means (roadmap open item 2). No rotation or bend has been planted on real noise.
- **The gain is unread where production measures the most dimensions.** ADNI's design is that one: rank 8 on the data and 19 in its simulation, against the two the predictive subspace measures. Its simulation (item 9) used feature surgeries that neither space sees.
- **The case study's non-rejections have no real-noise power statement.** `adni_sexdiff` calls them uninformative on the simulation alone.

Open item 2 names this slice next: plants on ADNI's blocks defined against held-out participants, their off-target statistics, and the gain at ADNI's design.

Probes on scratch seeds (2026-11-17, `notes.md`; none is a result) shaped the design:

- **The covariates carry half of ADNI's stage program.** The covariate-adjusted stage configuration has 0.51 to 0.64 of the unadjusted one's energy. A plant built on unadjusted stage means is therefore not the edit it claims under the covariate-adjusted trajectory model: a stretch of 2.0 was rejected by production's `angle` in 9 of 12 cohort units, and in 3 of 16 once the plant was built on the adjusted configuration.
- **A held-out stage configuration is mostly estimation noise.** Between two halves of the held-out participants the adjusted configuration has a reliability of 0.16 to 0.24, so about a third of a held-out plant's energy is the cohort's stage program. Its direction is fair to every space, but its off-target truth is not exact.
- **An in-sample plant has an exact truth and a rewarded direction.** Under within-stage label permutation both groups' expected configuration is the sample's own, so an edit of it is exact. At 750 participants the fitted span keeps 0.90 of the sample's own adjusted configuration and 0.49 of a held-out one.
- **The two kinds of plant are readable at different effects.** At 750 participants production's `delta` rejected a stretch of 0.5 in 10 of 16 units in sample and in 5 of 16 held out.
- **Stage-span plants compose with the ADNI-design simulation.** At ADNI's table the CV selected rank 19 in every unit, and a 90° rotation was rejected by `angle` in 1 of 8 units by production and in 7 of 8 by `predictive_orthonormal`.

## What Changes

- **A driver, and no package change.** `scripts/adni_gain_specificity_check.py` (`preflight` | `run` | `merge` | `decide`, with `--profile plants | simulation`). Every arm is one `analyze_trajectories` call on options that exist. `motco` stays at 0.25.0.
- **Profile `plants`: ADNI's blocks with a planted difference.** It imports the anchor driver's contract, subsampling, label draw and seeds unchanged. Each unit permutes the sex labels within stage and adds a known shift to the rows of the second group, stage by stage. No analysis runs under the real labels.
  - **Two design points:** 750 participants subsampled within cells (`n750`), which leaves 751 held out, and the full cohort of 1,501 (`cohort`).
  - **The plant edits a covariate-adjusted stage configuration** in pooled-standardized coordinates through the generator's own `edit_stage_configuration`: a stretch about the centroid, a rotation in the configuration's plane, or a bend of the middle stage.
  - **Two references.** `sample` edits the analysis sample's own configuration: exact off-target truth, available at both points. `heldout` edits the configuration of the participants outside the sample: a direction no fitted space was built on, available at `n750` only. Together they bracket a plant along the cohort's true stage program.
  - **Controls:** the anchor, and a translation along the net displacement, which plants a group main effect and no interaction.
  - **Recorded cells:** a rotation whose plane leaves the stage span.
- **Profile `simulation`: stage-span plants at ADNI's design.** The item 9 construction, imported unchanged at its calibrated strengths and ADNI's table, with the generator's stage-span construction under the planted scale in place of the feature surgeries. It runs on the cluster from a fresh seed family.
- **Arms.** Production's scores (`pls_cv`) and `predictive_orthonormal` on the same fit at the unit's CV rank, with the `df` correction, under RRPP's unrestricted draws. The same two under within-stage draws are recorded.
- **Replicates.** 500 × 999 on every control and 200 × 999 on every plant: 7,200 units on ADNI's blocks and 2,600 in simulation.
- **Pre-flights and a reproduction gate.** The cohort's fit must select rank 8 at AUROC 0.599514, as in item 21. The cohort anchor's 500 units are item 21's and must reproduce its committed counts. Two of item 9's anchor units must reproduce its committed p-values before any simulated unit runs.
- **One predeclared verdict on `predictive_orthonormal`** in `adni_gain_specificity_decision.json` (rule version 1), the first that applies:
  - `gate_failed`: a pre-flight or the reproduction failed, or a unit is inconsistent.
  - `production_excess`: production exceeds on a control.
  - `size_excess`: the candidate exceeds on a control this change adds.
  - `power_loss`: it loses target power against production on a plant.
  - `specificity_excess`: it adds off-target rejections over production, in simulation or on the in-sample plants.
  - `no_gain`: no simulated target cell at ADNI's table gains.
  - `reopened`: none of the above.

  Each verdict's consequence for roadmap open item 1 is written before any unit runs.
- **A second predeclared reading, of production alone on ADNI's noise.** It states from which effect each statistic detects a held-out and an in-sample plant, and whether production's own off-target rates on the in-sample plants hold the size rule. It never changes the verdict.
- **Descriptive readings beside it:** the share of each reference configuration and of each plant that the fitted span and the predictive subspace keep, the held-out configuration's reliability, the bracket between the two references, rates by rank class, the within-stage arms, and the out-of-span cells.
- **Documentation and records.** A dated report `docs/reports/adni-gain-specificity-<date>.md`; aggregates under `results/adni-gain-specificity-<date>/`; the verdict in roadmap open items 1 and 2, in the entry on real-noise power and specificity, and in the Done list; a section of `docs/measurement-visibility.md`; one row in the validation-drivers table of `AGENTS.md`.

**Non-goals.**

- **Changing a default or adopting `predictive_orthonormal`.** `reopened` lets a later change propose its paper-grade run with a rank rule; nothing here changes what `motco analyze` does.
- **A rank rule for a fixed measured dimension.** Open item 1 asks for it with any reopening; it is that later change's.
- **The split test and the interaction scan on ADNI's noise.** The out-of-span cells are where their real-noise reading belongs.
- **Design points below 750 participants.** At 80 and 200 the CV selects rank 2 in about half the units and the stage fit is near chance. Their size reading stands; nothing else is read there.
- **The other two predictive coordinate systems and the fixed-rank arms.** Only `predictive_orthonormal` at the CV rank has a path to reopening.
- **A plasmode with a planted stage program on stage-permuted noise,** the design that would make every held-out truth exact. It is the lead if the two references disagree.
- **ROSMAP.** It has no export yet.
- **Any change in `adni_sexdiff`.** Its reading rule may cite production's statement in its own change.

## Capabilities

### New Capabilities

- `cohort-design-power-validation`: power and specificity readings of the trajectory tests at a cohort's design, on plants added to the cohort's exported blocks under within-stage-permuted labels and on stage-span plants in a simulation of the design. It covers the plant and its two references, the design points and cells, the paired arms, the pre-flights and the reproduction gate, the per-unit checks, the tables, the criteria, the predeclared verdict on a candidate measurement with production's own statement, and the dated report. This change adds ADNI's reading; further cohorts and candidates extend it.

### Modified Capabilities

- `latent-space-metric-compatibility`: the measurement-visibility documentation reports the real-noise power and specificity readings on ADNI's blocks, the gain at ADNI's design in simulation, what each reference of a plant can and cannot show, and what the verdict licenses.

## Impact

- **Code:** none under `src/motco/`. No version bump and no `CHANGELOG.md` entry.
- **Scripts:** `scripts/adni_gain_specificity_check.py` and `scripts/adni_gain_specificity_array.sbatch` (the simulation profile only). The driver imports unchanged `adni_anchor_check` (contract, loader, design-point tables, subsample and label draws, seeds, arm parameters, per-unit checks), `adni_design_check` (the construction, its seeds and unit builder), `predictive_subspace_check` (arm records, the paired judgements), `production_latent_space_check` (the exact paired tests, the gain rule), `latent_space_comparison` (settings hashing, units) and `stage_span_paper_grade_decision` (`judge_rates`). It binds by sha256 the anchor run's `rate_table.csv` and `adni_anchor_decision.json` and the ADNI-design run's `settings.json` and `adni_design_pvalues.csv`, and refuses a changed file by name.
- **Tests:** `tests/test_adni_gain_specificity_check.py`, on the anchor test's synthetic stand-in export and on a small simulated table: the plant under each mode and reference, its exactness, the cells, seeds and matched units, the arms and their fit counts, the checks, records free of participant data, resume and refusals, the tables, every criterion, every verdict and production's statement.
- **Compute:** about 40 core-hours on the workstation for the plants profile, where the export is (about six hours at seven workers), and about 70 on `ing` for the simulation profile. A smoke of each profile replaces the estimates, and both runs wait on the user's approval.
- **Data:** ADNI is controlled-access. The export is read in place and never copied, so the plants profile has no sbatch array. Records hold seeds and statistics only: no participant identifier, label vector, subsample membership or plant vector. Only aggregates are committed.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged.
