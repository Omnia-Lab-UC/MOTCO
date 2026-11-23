# Proposal

## Why

The first comparison a reader makes is against the default analysis with no integration: a group × stage interaction test on each feature, corrected for multiple testing. MOTCO has never been judged against it. The package has that test as the interaction scan (0.20.0), but it was validated as a complement for sparse differences, not as a baseline.

The scan's own validation ran it beside `delta` / `angle` / `shape` on the same units (`docs/reports/interaction-scan-2026-11-09.md` §6; methylation and expression at 498 features, stage-span plants, 100 replicates per signal cell). The cells do not rank the two:
- **Production leads where a difference spreads thinly over many features.** Magnitude at e = 0.05 reads `delta` 0.74 against the scan's 0.09 at n = 200.
- **The scan leads on rotations, the in-span ones included.** At e = 0.25 on the SEA-AD table it reads 1.00 against `angle`'s 0.47.

That is a description, not a reading (roadmap open item 13):
- no decision was declared for the comparison;
- each trajectory statistic was tested at 0.05 on its own, three tests against the scan's one;
- the reading most analyses report, any discovery at a false discovery rate of 0.05, was not an arm;
- proteomics was not measured, and the feature-set surgeries behind the Phase 5 and design-guidance power claims predate the scan.

The roadmap now records that both real datasets read so far are convenience data, and that the intended sample is ROSMAP and ADNI with methylation, expression and proteomics. The simulation on those three blocks is therefore where the baseline is judged first. What MOTCO may claim against the default analysis depends on the answer, so it is read before the measurement refinements of open items 3 to 5.

## What Changes

- **A baseline validation driver.** `scripts/univariate_baseline_check.py` (`run --profile parent|stage_span|surgery` | `merge --profile …` | `decide`; sbatch array `scripts/univariate_baseline_array.sbatch`). Each unit generates one dataset, fits production once and reads every arm at one RRPP seed.
- **Production enters as one test.** `production` rejects when the smallest of the `delta`, `angle` and `shape` p-values is at or below 0.05 / 3. The three statistics are also recorded singly.
- **Three baseline arms**, each one test at 0.05 on the per-feature interaction F of the standardized joint features, under the trajectory test's full and reduced models:
  - `scan`: the maximum over features, by the interaction scan's permutation p-value;
  - `fdr`: Benjamini–Hochberg at 0.05 on the per-feature parametric F p-values, rejecting when any feature is a discovery;
  - `total`: the features pooled, by the split interaction test's `total` part.
- **Recorded arms that decide nothing:** the split test's `span` and `complement`, and the trajectory statistics measured on `concat` and on PCA at the production rank, so one table holds every baseline on the same units.
- **Three profiles.**
  - **`parent`**: the scan validation's simulation units for `n200` and the SEA-AD table (two blocks, base seed 2035, 14 cells; controls at 500 × 999, other cells at 100 × 999). 4,400 units. Every committed count of the scan, the three split parts and the three trajectory statistics must be reproduced exactly. Its new readings are descriptive, because these are the units the question came from.
  - **`stage_span`**: a fresh seed family (base seed 2043) on all three blocks, three stages. Four designs: balanced n = 80, 200 and 500 at δ = 2.0, and n = 200 at δ = 0.5. Thirteen cells each: the anchor, a translation, and three effects each of magnitude, in-span rotation and bend under the planted scale; an out-of-span rotation and its 20-feature form under the null scale. Controls at 500 × 999, other cells at 200 × 999. 12,800 units.
  - **`surgery`**: the three-block columns of the design-guidance paper-grade slice at n = 80, 200 and 500 (its seeds, 500 × 999). Eleven cells each: the anchor, translation at e = 1.0, magnitude at e ∈ {0.02, 0.05, 0.1}, and orientation and shape at e ∈ {0.25, 0.5, 1.0}. 16,500 units. The `pls_cv` arm must reproduce the slice's committed counts exactly.
- **A predeclared decision.** `univariate_baseline_decision.json` (rule version 1):
  - **Gate.** `reproduction_failed` if any `parent` or `surgery` reproduction fails.
  - **S, size.** `production` and each baseline arm at the anchors and translations of `stage_span` (m = 8) and of `surgery` (m = 6), under the stage-span paper grade's count-plus-ceiling rule (`judge_rates`, imported unchanged). An arm is admissible when both families hold. A baseline arm that is not admissible leaves the comparison and is reported.
  - **D, detection.** Over the 71 signal cells of `stage_span` and `surgery`, exact one-sided paired tests of `production` against each admissible baseline arm on the units where exactly one of the two rejects. A cell is labelled `production_leads` when production gains against every admissible arm, `baseline_leads` when it loses against at least one, and `unresolved` otherwise, with the familywise error over all labels held at 0.05.
  - **Verdict.** `production_leads` | `baseline_leads` | `complementary` | `unresolved`, by whether each label occurs in at least one cell. `production_size_excess` and `no_admissible_baseline` come before them.
  - **C, characterization grade.** The scan has no geometry, so naming the kind of a difference is production's claim alone. Over the in-span `stage_span` cells that `production` detects in at least 0.80 of units, the off-target statistics are judged under the count-plus-ceiling rule: `names_kind` | `detects_only` | `no_detected_cell`. The grade is reported beside the verdict and does not change it.
  - **Descriptive beside it:** the map of labels by plant family and design; every arm's rate per cell; discoveries per unit and the realized false discovery proportion; selected rank, held-out AUROC and span capture.
- **Documentation and records.** A dated report `docs/reports/univariate-baseline-<date>.md`; aggregates under `results/univariate-baseline-<date>/`; a section of `docs/measurement-visibility.md` on the default analysis as a baseline; roadmap open item 13; a row in the validation-drivers table of `AGENTS.md`.

**Non-goals.**
- **Any package change.** `motco.stats`, the harness, `motco.analysis` and the CLI are untouched and the version does not move. The `fdr` reading lives in the driver. Adding it to the scan's output is a later change, opened only if the verdict gives it a use.
- **A default or recommendation flip.** No verdict changes what `motco analyze` runs. The verdict fixes what the documentation may claim.
- **Real noise.** That leg waits for multi-omic ROSMAP and ADNI exports (roadmap open item 2). SEA-AD's and ADNI's metabolomics plants are not arms here.
- **Dimension.** The generator has 658 features, fixed by its InterSIM reference. Neither method is read at the dimension real blocks have; the report states that limit beside every claim.
- **A moderated per-feature statistic** (roadmap open item 6) and per-block maxima as further arms.
- **Localization.** Which features differ, and how attribution compares with the scan's hits, is not judged.

## Capabilities

### New Capabilities
- `univariate-baseline-validation`: the judged comparison of the production trajectory test against the per-feature interaction test with no integration. It covers the arms on one unit, the three profiles and their reproduction of the committed parents, the unit record, the predeclared decision with its characterization grade and descriptive tables, and the dated report.

### Modified Capabilities
- `latent-space-metric-compatibility`: the measurement-visibility documentation states how the trajectory test compares with the default per-feature analysis, what each can and cannot claim, and the validation's verdict, grade and limits.

## Impact

- **Code:** none under `src/`. `motco.stats.interaction`, the harness and `motco.analysis` are imported as they are.
- **Scripts:** `scripts/univariate_baseline_check.py` and `scripts/univariate_baseline_array.sbatch`. They import unchanged:
  - `latent_space_comparison` (designs, seeds, `unit_dataset`);
  - `span_split_check` and `interaction_scan_check` (the parent's units, the split and scan inputs, the population-null mask);
  - `stage_span_leak_check` (study-engine columns and unit parameters);
  - `production_latent_space_check` (the exact paired p-values);
  - `stage_span_paper_grade_decision` (`judge_rates`).
- **Committed parents bound by sha256:** `results/interaction-scan-2026-11-09/simulation/` (`settings.json`, `comparator_table.csv`) and `results/design-guidance-slice-2026-10-29/core_{magnitude,geometry}/report/design_point_operating.csv`, with the two study configs.
- **Tests:** `tests/test_univariate_baseline_check.py`: unit counts and seeds per profile, the arms on a small design (the `pls_cv` record equals a plain production evaluation; one cross-validated fit per unit), the `fdr` arm against a hand-computed case, resume and settings refusals, a refused parent, both outcomes of each reproduction, and every decision path on synthetic tables.
- **Dependencies:** none added. `scipy` is already a runtime dependency; the locked version (1.16.3) provides the F distribution and Benjamini–Hochberg.
- **Compute:** about 33,700 units with one production CV each. The design-guidance slice ran 44,000 such units in 172 core-hours (median 12 s), and a unit here adds a split test, a scan and two extra measurement passes. The estimate before the smoke is 130–200 core-hours on `ing`. A local smoke replaces the estimate, and launch waits on the user's approval with a cap.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged. The parents' committed files are read, hashed and refused by name if changed.
