# plasmode-attribution-demonstration

## Why

Phase 6 closed on 2026-11-02 with its attribution leg withheld ([exit review](../../../docs/reports/phase6-exit-review-2026-11-02.md) §1 G1, §6). Neither case study has a global trajectory claim to trace, so `motco analyze --attribution` has never run on real data at the configuration in which a claim would be read. Item 11 (`covariate-adjusted-attribution`, MOTCO 0.14.0, `adjusted_attribution_validated`) removed the capability gap: attribution now runs with declared covariates, on semi-synthetic blocks.

This change is roadmap item 12, the plasmode demonstration. It plants a known group × stage difference in the real exported blocks, covariates and design of SEA-AD astrocyte `cohort_30`, analyses each unit as the case study would, and asks whether attribution traces the detected difference back to the planted features. It either lifts the withheld leg or reports it failed, naming the criterion that failed. The exit review left three things to this proposal, and it fixes them here before any run: the null-making step, the detection rate and the recovery thresholds.

## What Changes

- **Substrate and input contract.** A new script, `scripts/plasmode_check.py`, reads `sead_sexdiff`'s export directory (`results/motco_input/cohort_30`) through a documented contract:
  - `manifest.json`, the `rna` and `atac` blocks (5,000 log2-CPM features each, `identity` transform), `metadata.csv` and `covariates.csv`;
  - the primary's six covariate columns (`age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna`, `log10_n_atac`). The two log10 columns are derived as `sead_sexdiff` derives them.

  The sha256 of every input file is recorded, and a resume on a different export is refused. Nothing donor-level is committed to MOTCO.
- **Label-free fidelity gate.** The latent space is stage-supervised and covariate-blind, so it does not depend on the sex labels. Before any unit, one analysis on permuted labels without permutations must reproduce the primary's CV-selected rank (5) and mean held-out stage AUROC (0.635556). Every unplanted unit is held to the same check. No quantity that depends on the real sex labels is computed or read.
- **Null-making step for criterion (i).** Each replicate permutes the sex labels **within stage**. The permutation keeps every cell size (`[[11, 10, 28], [9, 10, 12]]`) and leaves each donor's blocks and covariates on the donor's own row. A draw that reproduces the real labels is refused and redrawn. The analysis under the real labels is never run.
- **Planted orientation.** In the donors labelled `Male` in the permuted assignment, each feature in a planted set K is shifted by `e · c(stage) · σ_f · s_f`:
  - K is 10 RNA and 10 ATAC features, drawn at random per replicate;
  - `s_f` is a random sign;
  - `σ_f` is the feature's pooled within-stage SD;
  - `c` is 0, 0.5 and 1 at stages 0, 1 and 2.

  The shift is a group-specific linear ramp on features outside the shared trajectory. It turns that group's net displacement and also lengthens it. At e = 0 the blocks equal the export byte for byte. The strength `e` is in within-stage SD units at the last stage.
- **Calibration of the strength.** The ladder is e ∈ {0.5, 1.0, 1.5, 2.0, 3.0}. It is climbed with 40 calibration units per step at 99 permutations, with attribution off, so calibration cannot see criterion (iii). The adopted strength `e*` is the first step at which `angle` rejects in at least 0.90 of the calibration units. If no step reaches that rate, `e*` is the cap of 3.0 and the run proceeds; criterion (ii) then decides. Every iteration is recorded.
- **Units and arms.** Three cells × 500 replicates at 999 permutations, with attribution at B = 200 and top-k = 20:
  - the anchor (e = 0);
  - `e*`;
  - the ladder step below `e*` (`e*`/2 when `e*` is the lowest step), descriptive only.

  The cells share each replicate's permutation, K, signs and seeds (base seed 2034, which no committed check uses). Each unit is analysed **once** by `analyze_trajectories` at the case-study configuration: the six covariates, rank floor 2, net displacement and `residual_correction="df"`, with covariate-adjusted attribution. That analysis is the primary arm. An unadjusted attribution on its `latent_artifacts` is recorded as a descriptive arm, which cannot lift the leg (exit review D6).
- **Predeclared lift rule.** `plasmode_decision.json` returns `attribution_leg_lifted` when all of the following hold, and otherwise `attribution_leg_failed`, naming every failing criterion and cell:
  - **(i) size:** at the anchor, the rejection rate of each of `delta`, `angle` and `shape` is at most α + 2·√(α(1 − α)/R), which is 0.0695 at R = 500;
  - **(ii) detection:** at `e*`, `angle` rejects in at least 0.80 of the units (the roadmap's power floor);
  - **(iii) recovery and stability:** on the units at `e*` whose `angle` rejects, on the observed principal-orientation contrast, each lower bound (mean − 2·SE) must clear its threshold:
    - the planted-set precision of the top-20, |top-20 ∩ K| / 20, at least 0.50;
    - the mean bootstrap sign stability over K, at least 0.90;
    - the mean bootstrap top-k frequency over K, at least 0.50.

  The verdict is indicative below 500 replicates, 999 permutations or B = 200. The following are reported but cannot change the verdict:
  - the unadjusted arm;
  - the PLS-captured and residual components;
  - the per-transition blocks;
  - the readings that are not conditioned on rejection;
  - `delta` and `shape` at `e*`;
  - the descriptive cell.
- **Run, report and roadmap.** The sequence is:
  1. a local smoke run and a one-unit timing, then an approval gate;
  2. calibration, then the paper-grade run (on local workers or on `ing`, chosen from the timing), then merge and decide;
  3. a dated report and `results/plasmode-<date>/` holding settings, calibration, the per-cell summary, the decision and PROVENANCE. Per-unit records stay on the run machine, named by sha256.

  Updates follow to `docs/roadmap.md` (item 12 and the Phase 6 Work bullet) and `AGENTS.md`.
- **Out of scope.**
  - Planting magnitude or shape differences; the orientation kind comes first, and the others are follow-ups.
  - Any analysis or reading of SEA-AD's real sex × stage result.
  - ADNI.
  - Pathway or module enrichment.
  - Any change to `src/motco/`. MOTCO 0.14.0 already does everything the demonstration needs, so no release is made.
  - Any change to `sead_sexdiff`. It supplies the export as it stands.

## Capabilities

### New Capabilities

- `plasmode-attribution-demonstration`: the plasmode check on a real case-study export. It covers:
  - the input contract and the label-free fidelity gate;
  - the within-stage null-making step and the planted orientation construction;
  - strength calibration, units, arms and metrics;
  - the predeclared lift rule and its two verdicts.

### Modified Capabilities

None. The check calls `analyze_trajectories` and the attribution primitive as released in 0.14.0, and no library requirement changes.

## Impact

- **Code:**
  - a new `scripts/plasmode_check.py` with the subcommands `fidelity`, `calibrate`, `run`, `merge` and `decide`;
  - `scripts/plasmode_array.sbatch`, only if the run goes to `ing`;
  - no change under `src/motco/`.
- **Tests:** a new `tests/test_plasmode_check.py` on a synthetic stand-in export in `tmp_path`, with no real data. It covers:
  - contract validation;
  - the permutation's invariants;
  - the plant's formula and its byte identity at e = 0;
  - seed matching and resume refusal;
  - the metrics;
  - the calibration step;
  - every decision scenario.
- **Data:** the SEA-AD export is read in place, and on `ing` from a copy whose sha256 is verified. It is never committed. SEA-AD is open data.
- **Results and docs:**
  - `results/plasmode-<date>/` (aggregates only);
  - `docs/reports/plasmode-<date>.md`;
  - `docs/roadmap.md`;
  - `AGENTS.md`.
- **Downstream:**
  - Under `attribution_leg_lifted`, the Phase 6 G1 leg is recorded as demonstrated on planted data, and the paper cites the report.
  - Under `attribution_leg_failed`, the leg is reported failed, with the criterion named.
  - Neither outcome reopens Phase 6 or changes any case-study reading.
- **Compute:** unknown until the smoke timing. At about 30 s per unit, the 1,500 paper-grade units come to about 12.5 core-hours, plus up to 200 calibration units at 99 permutations.
