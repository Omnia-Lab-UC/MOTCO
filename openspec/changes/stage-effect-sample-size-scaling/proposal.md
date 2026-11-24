# Proposal

## Why

At a weak stage effect the trajectory test degrades, and every reading of that is at n = 200 or 80:
- **The CV rank inflates.** Its median is 2 at δ = 2.0, 4 at δ = 0.5 and 11 at δ = 0.25 (design-guidance slice §8; n = 200, two blocks).
- **`delta` and `shape` lose power.** At e = 0.1 `delta` falls from 0.94 to 0.31 at δ = 0.5.
- **A rejection stops naming the kind of a difference.** At n = 200 and δ = 0.5 a 45° in-span rotation is read as `shape` in 0.665 of units, which made the univariate baseline's grade `detects_only`.

The case studies are weakly staged (held-out stage AUROC 0.636 and 0.600, against 1.000 at the core stage effect), and ADNI has 1,501 participants. Whether more samples restore the power, the rank and the naming is unread (roadmap open item 14).

The generator gives the question a reference. The stage program's means are δ times a fixed pattern over a fixed noise covariance, so a design's information about the program, and about a plant defined relative to it, grows as δ√n. If a weak stage effect were only a smaller effective sample, a design with half the stage effect and four times the samples would read the same, cell by cell. Nothing committed tests that, because no two committed designs share δ√n:
- **Stretches line up on e·δ·√n** across the baseline's four designs: 0.285 at 0.71, 0.345 at 0.89, 0.815 at 1.41, 0.810 at 1.77, 0.750 at 1.79 and 0.985 at 2.24. That is a description.
- **A timing probe suggests the rank does not follow δ√n** (2026-11-24; two null datasets per design at a seed outside every run; not a reading). It selected rank 3 at δ = 0.5 with n = 1,280, where the matched design at δ = 2.0 with n = 80 selects 2, and rank 19, the cap, at δ = 0.1 with n = 1,280.

The answer decides what a weakly staged cohort is told: collect more samples, or expect a rank rule to be needed (roadmap, On hold). It needs no package change, so it is read before the other simulated conditions (open items 8 and 15).

## What Changes

- **A validation driver.** `scripts/stage_effect_scaling_check.py` (`calibrate` | `run --profile saturated|band` | `merge --profile …` | `decide`; sbatch array `scripts/stage_effect_scaling_array.sbatch`). The unit, its arms and its tables are the univariate baseline's, imported unchanged: one dataset, one production fit, `production` as one test at 0.05 / 3, the per-feature arms `scan`, `fdr` and `total`, and the recorded arms.
- **Profile `saturated`: matched designs where the stage AUROC is at its ceiling.** Three blocks, three stages, equal groups, the baseline's thirteen cells and grade (controls at 500 × 999, other cells at 200 × 999), base seed 2043 with new design names. 16,000 units.
  - **The gate.** The baseline's `n200_delta0.5` design is rerun and must reproduce its committed counts exactly.
  - **Four new designs**, each matched in δ√n to a committed one:

    | new design | matched to | δ√n |
    |---|---|---|
    | δ = 1.0, n = 320 | δ = 2.0, n = 80 | 17.9 |
    | δ = 0.5, n = 1,280 | δ = 2.0, n = 80 | 17.9 |
    | δ = 1.0, n = 800 | δ = 2.0, n = 200 | 28.3 |
    | δ = 0.25, n = 800 | δ = 0.5, n = 200 | 7.1 |

- **Profile `band`: designs below the AUROC ceiling.** A `calibrate` step (fits only, no permutation) finds the stage effect δ_c whose held-out AUROC at n = 320 is nearest the latent-space comparison's target, 0.636, inside its band, 0.60–0.67. Three designs: δ_c at n = 320; δ_c / 2 at n = 1,280, its matched design; and δ_c at n = 1,280. Eleven cells each under the planted scale: the anchor, a translation, stretches at e ∈ {0.5, 1.0, 2.0}, and rotations and bends at e ∈ {0.25, 0.5, 1.0}. 8,400 units.
- **A predeclared decision.** `stage_effect_scaling_decision.json` (rule version 1):
  - **Gate.** `reproduction_failed` if the rerun design differs from its committed counts.
  - **S, size.** `production`, `fdr` and `total` at the anchors and translations of the new designs, in two families (`saturated`, m = 8; `band`, m = 6), under the count-plus-ceiling rule (`judge_rates`, imported unchanged). `production_size_excess` if production fails either.
  - **T, the matched comparison.** For each of the five matched pairs and each signal cell (53 pair-cells), `production`'s rate in the weaker-staged design against its match, by a one-sided exact test in each direction with the familywise error held at 0.05 and a margin of 0.10: `weak_costs` | `weak_gains` | `unresolved`. The `fdr` arm is the control of the pairing: a label counts toward the verdict unless `fdr` carries the same label in that pair-cell.
  - **Verdict.** `weak_stage_costs` | `weak_stage_gains` | `mixed` | `scaling_unrefuted`, by which counted labels occur.
  - **N, naming per design.** The baseline's characterization grade (`names_kind` | `detects_only` | `no_detected_cell`), computed for each design on its own, the four committed ones included. It is reported beside the verdict and does not change it.
  - **M, the map per design.** The baseline's labels of `production` against the admissible per-feature arms on each new design's signal cells, by its paired exact tests.
  - **Descriptive beside it:** per design δ, n, δ√n, the held-out AUROC, the CV rank and the span capture; every arm's rate per cell; each statistic's rate against e·δ·√n; rates by rank stratum.
- **Documentation and records.** A dated report `docs/reports/stage-effect-scaling-<date>.md`; aggregates under `results/stage-effect-scaling-<date>/`; a section of `docs/measurement-visibility.md` on what a weak stage effect costs; roadmap open item 14; a row in the validation-drivers table of `AGENTS.md`.

**Non-goals.**
- **Any package change.** `motco.stats`, the harness, `motco.analysis` and the CLI are untouched and the version does not move.
- **A rank rule.** The run says whether sample size restores what a weak stage effect costs. A rule for the fit rank under a fixed measured dimension stays on hold; `weak_stage_costs` with the costs in the inflated-rank strata is what would reopen it.
- **A default or recommendation flip.** No verdict changes what `motco analyze` runs.
- **Real noise.** That leg waits for multi-omic ROSMAP and ADNI exports (roadmap open item 2).
- **Unbalance, group-specific dispersion and feature count** (roadmap open items 8 and 15), two-block and four-stage designs, and samples beyond n = 1,280.
- **Matching a cohort.** The `band` designs sit at the case studies' held-out AUROC; they do not reproduce either cohort's cell table, covariates or rank. ADNI's design has its own simulation (roadmap item 9).

## Capabilities

### New Capabilities
- `stage-effect-scaling-validation`: the judged reading of how a weaker stage effect trades against sample size for the production trajectory test. It covers the matched designs and their reproduction of the committed parent, the calibrated designs below the AUROC ceiling, the predeclared decision with its control arm, the per-design naming grade and map, the descriptive tables and the dated report.

### Modified Capabilities
- `latent-space-metric-compatibility`: the measurement-visibility documentation states what a weak stage effect costs the trajectory test, whether sample size buys it back, and how a design is placed by its held-out AUROC, rank and δ√n.

## Impact

- **Code:** none under `src/`.
- **Scripts:** `scripts/stage_effect_scaling_check.py` and `scripts/stage_effect_scaling_array.sbatch`. They import unchanged:
  - `univariate_baseline_check` (the stage-span generator and cells, the unit and its arms, the rate, paired and cell tables, the labelling rule and the characterization grade);
  - `latent_space_comparison` (seeds, the calibration unit, its target and band, settings binding);
  - `stage_span_paper_grade_decision` (`judge_rates`).
- **Committed parent bound by sha256:** `results/univariate-baseline-2026-11-23/stage_span/` (`settings.json`, `rate_table.csv`, `cell_table.csv`).
- **Tests:** `tests/test_stage_effect_scaling_check.py`: designs, unit counts and seeds per profile; that each new design differs from its match only in δ and n and has the same δ√n; the rerun design's units equal the baseline's; the calibration's choice and its refusal outside the band; the exact tests against hand-computed cases; resume and settings refusals; a refused parent; both outcomes of the reproduction; and every decision path on synthetic tables.
- **Dependencies:** none added. `scipy` (1.16.3, locked) provides the exact test.
- **Compute:** 24,400 units and 240 calibration fits. The probe measured 25 s per unit at n = 320, 40 s at n = 800 and 75 s at n = 1,280 on the workstation (seven workers), which gives about 180 core-hours for `saturated` and 136 for `band`. The baseline's stage-span array ran at 1.4 times its unpacked projection, so the estimate before the smoke is 320–450 core-hours on `ing`. A local smoke replaces it, and launch waits on the user's approval with a cap.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged. The parent's committed files are read, hashed and refused by name if changed.
