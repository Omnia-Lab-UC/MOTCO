# Proposal

## Why

Three measurement candidates were kept out of production by real-noise readings on convenience data, not by simulation:

- **`predictive_orthonormal`** (roadmap items 18 and 22). It failed only the SEA-AD plasmode's size (`delta` 0.134 against production's 0.042) and, on ADNI's metabolomics export, protected power on an in-sample 90° rotation. In the pilot's simulation it held size, lost no protected power, passed both specificity families and gained in 10 of 18 target cells at δ = 2.0, at 100 replicates.
- **Sparse stage-PLS** (`sparse_block`, item 17). Its `size_excess` is one reading, `delta` 0.11 at 100 replicates on the SEA-AD anchor, where production read 0.11 on the same replicates.
- **Block weighting** (`weighted`, item 17). Its `no_gain` came from a gain rule that counted only the weak design and the SEA-AD plant. At δ = 2.0 it gained (`angle` 0.46 → 0.90–1.00), and that gain was never judged.

Neither export is the data MOTCO was designed for (roadmap, Case studies), each is one dataset at one design, and the SEA-AD anchor's replicates share one fit. On 2026-11-25 the user decided to re-read these candidates in simulation alone, across more conditions, with the real-noise readings kept as recorded limits and not as vetoes.

The simulated evidence is thin where it matters:

- Every favourable reading is at 100 replicates × 199 permutations, three stages, two blocks and one cell table per design.
- No candidate has been read at four stages, on three blocks, under unbalance beyond the SEA-AD table, or under unequal noise.
- The one mechanism the real data taught is reproducible in simulation and has not been read through the generator: when the largest stage has the narrowest noise the unrestricted null over-rejects `delta` up to 0.102 on Gaussian matrices (item 20).

One outcome is known before any unit runs. In the simulation of ADNI's design `predictive_orthonormal` already fails specificity at paper grade: a doubled stretch is read as `shape` in 0.295 of units against production's 0.06 (item 22, S3 on simulated plants: 3 excess readings against 2 allowed). That reading is simulated, so it stays in. A default change for `predictive_orthonormal` is therefore not reachable; what this change can establish is at which simulated designs each candidate is better than production, and at which it is worse.

The other three candidates of the 2026-11-25 summary were rejected by simulation alone and are not re-run: orthonormal span coordinates (item 15), the within-stage null (item 20) and the latent covariate adjustment (item 3).

## What Changes

- **An opt-in stage dispersion scale in the generator.** `SemiSyntheticTrajectoryParams.stage_dispersion_scale` holds one positive multiplier of the noise standard deviation per stage, applied to both groups. The means do not move and no random number is drawn. It composes with `group_dispersion_scale` (0.26.0) as a product per cell. All ones normalizes to `None`, the default, which is hash-omitted, so every committed dataset, truth record, cell id, signature and config is unchanged. Under `stage_span` the pooled scales use each cell's own noise variance. MOTCO 0.27.0.
- **A validation driver.** `scripts/candidate_simulation_check.py` (`run --profile …` | `merge --profile …` | `decide`; an sbatch array beside it). One dataset per unit, read by four arms at the unit's RRPP seed:
  - `pls_cv`: production;
  - `predictive_orthonormal`: built from the production fit, with no second fit;
  - `weighted`: its own cross-validated fit with `block_weighting = "equal_variance"`;
  - `sparse_block`: its own cross-validated sparse fit with the `block` scope and the default ladder, on a named subset of designs.
- **Five profiles, 58,400 units**, at 999 permutations, controls at 500 replicates and plants at 200 unless stated:

  | profile | designs | cells | units |
  |---|---|---|---|
  | `committed` | the univariate baseline's `n80`, `n200`, `n500` and `n200_delta0.5`, and the stage-effect scaling run's `n320_delta1` and `band_n320`, at their committed seeds | their committed thirteen (eleven in the band) | 18,800 |
  | `conditions` | new: four stages at n = 200 (δ = 2.0 and 0.5); the SEA-AD cell table (δ = 2.0 and 0.5); a group share of 0.8 at n = 200 | the baseline's thirteen | 16,000 |
  | `dispersion` | new, null only: four stage-dispersion patterns on the SEA-AD table (δ = 2.0 and 0.5) and the two stronger ones at balanced n = 200 | anchor and translation | 10,000 |
  | `surgery` | the design-guidance slice's three-block columns at n = 80 and 200, at their committed seeds | their eleven, at 500 replicates | 11,000 |
  | `adni_design` | the ADNI gain run's simulated design, at its committed seeds | its ten | 2,600 |

  - **Reproduction gates.** On `committed`, `surgery` and `adni_design` the `pls_cv` counts, and on `adni_design` the `predictive_orthonormal` counts, must equal the committed rate tables exactly. No candidate has been read on the `committed` or `surgery` datasets.
  - **`sparse_block`** runs on `n200`, `n200_delta0.5` and `band_n320` only (9,200 units); elsewhere it is recorded as not run.
  - **Recorded arms** on `dispersion` only: `pls_cv` and `predictive_orthonormal` under within-stage draws, on the same fit and seed.
- **A predeclared decision, scoped by design family.** `candidate_simulation_decision.json` (rule version 1):
  - **Gate.** `reproduction_failed` if any reproduction differs.
  - **Families.** `core` (δ = 2.0, three stages, balanced), `weak_stage`, `four_stages`, `unbalanced` and `adni_design`.
  - **Per arm and family**, the earlier pilots' criteria with their rules unchanged: S1 size of the controls (count-plus-ceiling), S2 protected power on the in-span target cells (exact paired loss), S3 specificity of the off-target statistics (the bound at δ = 2.0, exact paired excess over `pls_cv` below it), G gain (one Bonferroni-protected exact paired gain). Label: `size_excess` | `power_loss` | `specificity_excess` | `no_gain` | `qualifies` | `not_run`.
  - **D, stage dispersion.** Each arm's controls under the count-plus-ceiling rule. If production itself fails, the family becomes a paired excess test against `pls_cv` for every candidate, and production's failure is reported as its own finding.
  - **I, incidental routes** (a rider, not a veto): on the feature-set surgeries and the out-of-span rotations, each candidate is labelled `candidate_loses` | `candidate_gains` | `unresolved` against `pls_cv` per cell.
  - **Per arm verdict.** `qualifies_everywhere` (every family run and qualifying, D held) | `scoped_candidate` (at least one family qualifies) | `keep_production`.
  - **Overall.** `propose_default_change` naming one arm | `scoped_candidates` | `keep_production`.
- **Documentation and records.** A dated report `docs/reports/candidate-simulation-<date>.md`; aggregates under `results/candidate-simulation-<date>/`; a section of `docs/measurement-visibility.md` with the map by family; `docs/api/simulations.md`; the roadmap, including the amended reopening condition for `predictive_orthonormal` and the real-noise readings restated as limits; `AGENTS.md`; `CHANGELOG.md`.

**Non-goals.**
- **Any default.** No verdict changes what the harness or `motco analyze` runs. `propose_default_change` opens a separate change with pins for every committed configuration.
- **Real noise.** No export is read. The SEA-AD and ADNI readings stay in the record as they are.
- **Many features per sample.** The SEA-AD excess was read at 125 features per donor, and no simulated design has reproduced it; the generator holds 658 features (roadmap open item 15). This run cannot say whether a candidate holds size there, and the report says so.
- **`predictive`, `predictive_fitted`, `weighted_sparse_block`, orthonormal span coordinates and the latent adjustment.** Each failed in simulation or has no gain to confirm.
- **The within-stage null as a candidate.** It is recorded on the dispersion profile only.
- **A rank rule, a CLI flag for `latent_coordinates`, attribution or the split test on the predictive subspace.**
- **`motco simulate` as a user of the new setting**, and dispersion by feature.

## Capabilities

### New Capabilities
- `candidate-simulation-validation`: the judged, simulation-only reading of the measurement candidates against the production measurement. It covers the four arms on one dataset, the five profiles and their reproduction of the committed parents, the decision scoped by design family with its dispersion criterion and incidental-route rider, the recorded arms and the dated report.

### Modified Capabilities
- `numpy-omics-generator`: a sampling call accepts one noise scale per cell.
- `semisynthetic-trajectory-generator`: the generator gains the opt-in `stage_dispersion_scale`, with its validation, its truth record, its composition with the group scale, its invisibility while absent and the stage-span scales under it.
- `trajectory-power-study`: study configurations accept `generator.stage_dispersion_scale`, hash-omitted while absent and kept by null cells.
- `latent-space-metric-compatibility`: the measurement-visibility documentation states the simulation-only map of the candidates by design family, and which limits the real-noise readings and the feature count leave on it.

## Impact

- **Depends on** `group-specific-dispersion` (0.26.0, `generate_omics(noise_scale=)` and the within-cell multiplier of the pooled scales), which is on its own branch and not yet merged to `main`. This change is applied after that merge.
- **Code:**
  - `src/motco/simulations/generator.py`: `noise_scale` accepts one value per cell;
  - `src/motco/simulations/semisynthetic.py`: the field, its validation and normalization, the per-cell scales passed to both groups, the truth record and the pooled scales' within-cell term;
  - `src/motco/simulations/grid.py` and `src/motco/simulations/study/config.py`: the field in the omit-when-`None` sets;
  - `src/motco/__init__.py` and `CHANGELOG.md`: 0.27.0.

  `motco.stats`, `motco.analysis`, `motco.cli` and the evaluation harness are not touched.
- **Scripts:** `scripts/candidate_simulation_check.py` and `scripts/candidate_simulation_array.sbatch`. They import unchanged `univariate_baseline_check` (the stage-span and surgery builders and cells), `stage_effect_scaling_check` (its designs and calibration), `adni_gain_specificity_check` (the simulated design and its seeds), `production_latent_space_check` (the weighted and sparse arms' parameters, the paired tests), `latent_space_comparison` (seed functions, settings binding) and `stage_span_paper_grade_decision` (`judge_rates`).
- **Committed parents bound by sha256:** `results/univariate-baseline-2026-11-23/{stage_span,surgery}/`, `results/stage-effect-scaling-2026-11-24/{saturated,band,calibration}/` and `results/adni-gain-specificity-2026-11-17/simulation/`.
- **Tests:** `tests/test_stage_dispersion_scale.py` (the field) and `tests/test_candidate_simulation_check.py` (units, seeds, parents, arms on one dataset, resume, merge, every reproduction outcome and every decision path).
- **Compute:** about 58,400 units. Before the smoke the estimate is 900–1,300 core-hours on `ing`, roughly three times the largest run so far, of which `sparse_block` is about a third and `weighted` about a quarter. The run is launched in three tiers, each with its own cap and approval: (1) `committed`, `conditions` and `dispersion` without `sparse_block`; (2) `surgery` and `adni_design`; (3) `sparse_block`. A tier that is not launched leaves its families `not_run`.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged. The setting is absent by default, and the parents' committed files are read, hashed and refused by name if changed.
