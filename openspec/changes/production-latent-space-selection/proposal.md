# Proposal

## Why

Two candidate changes to the production stage-PLS space exist as opt-in settings without a verdict: block weighting (`block-weighting`, roadmap item 5, which requires a size and power check before the case study may use it) and sparse stage-PLS (`regularized-stage-pls`). Each is a change to the fit's input or regularization, not to the stage objective, so each keeps what the latent-space comparison (2026-11-04) found irreplaceable and asks what it left open: how much noise the production fit carries at the regime where real data sit (a 0.10 share of the stage program at SEA-AD's stage signal; 0.23 → 0.07 of `delta` power lost on the SEA-AD table when two empty directions are added). The question for production is whether either setting, or both, holds size, loses no in-span power where production has it, and gains power where production has none — on matched data, under a predeclared rule, before any default moves.

## What Changes

- **A paired four-arm check on one dataset and one set of draws per unit.** `scripts/production_latent_space_check.py` (`run --profile simulation|plasmode` | `merge --profile …` | `decide`; sbatch array `scripts/production_latent_space_array.sbatch`, `PLS_PROFILE`) evaluates every unit under the weighting × sparsity factorial, each arm its own stage-PLS fit at net displacement, rank floor 2 and span capture, at the unit's RRPP seed, serially:
  - `pls_cv`: production;
  - `weighted`: `block_weighting: equal_variance`;
  - `sparse_block`: `sparsity: {scope: block, default ladder}`;
  - `weighted_sparse_block`: both.
- **Two profiles.**
  - **`simulation`** reuses the latent-space comparison's three designs (n = 200 and the SEA-AD table at δ = 2.0; the SEA-AD table at the committed weak δ = 0.11), base seed 2035 and seed functions, at 100 × 199. Cells: the comparison's seven in-span cells by name (`anchor`, `translation`, `magnitude` e = 0.05, `orient_e0.25_phi0`, `orient_e0.5_phi0`, `shape_e0.25`, `shape_e0.5`), so their datasets and draws are the comparison's and the `pls_cv` arm must reproduce its `rate_table.csv` rows exactly, plus three larger plants at every design for a power gradient the weak design lacks (`magnitude_e0.25`, `orient_e1.0_phi0`, `shape_e1.0`). 10 cells × 3 designs × 100 = 3,000 units. Out-of-span cells (φ > 0, sparse support) are excluded: every arm is stage-supervised and none can see them by construction.
  - **`plasmode`** reads the SEA-AD astrocyte export through the plasmode's input contract (file hashes equal to the committed plasmode settings'), permutes the sex labels within stage as the plasmode does (base seed 2034, the `paper` stream, replicates 0–99) and analyses each unit with `analyze_trajectories` at the case-study configuration (six covariates, floor 2, net displacement, `df`, attribution off) under the four arms. Cells: the `anchor` and a **real-noise span plant** at e ∈ {0.25, 0.5, 1.0}: the donors labelled `Male` at stage *s* are shifted by `e · (m̄_s − m̄)` in the analysis's standardized units — the stage-span magnitude edit with the pooled sample stage configuration standing in for the population — so B's trajectory is A's stretched by 1 + e about the centroid on real noise, in the stage span by construction. Every anchor unit's `pls_cv` arm must reproduce the plasmode's label-free latent space (rank 5, AUROC 0.635556). 4 cells × 100 = 400 units at 999 permutations.
- **A predeclared decision.** `production_latent_space_decision.json` (rule version 1) judges each candidate arm against `pls_cv`:
  - **S1 size:** anchor and translation `delta`/`angle`/`shape` at the three designs (m = 18) and the plasmode anchor (m = 3), each ≤ α + 2·√(α(1 − α)/R), the comparison's rule;
  - **S2 protected power:** at the δ = 2.0 designs, the 16 target cells (magnitude `delta`, orientation `angle`, shape `shape`) by a one-sided exact paired loss test under the count-plus-ceiling rule (at most the 95th percentile of Binomial(16, 0.025) cells below 0.025, none below 0.05/16);
  - **S3 specificity:** the off-target rates of those cells (m = 32) under the S1 bound;
  - **G gain:** the weak design's 8 target cells and the plasmode's 3 span-plant `delta` cells by a one-sided exact paired gain test; G holds when any cell's gain p-value is at most 0.05/11.

  Per arm: `size_excess` | `power_loss` | `specificity_excess` | `no_gain` | `qualifies`, in that order. Overall: `propose_paper_grade` naming the simplest qualifying arm (`weighted`, then `sparse_block`, then `weighted_sparse_block`) | `keep_production`; `reproduction_failed` gates everything. Descriptive beside it: per cell and arm the selected rank and keep fraction (median, q10/q90), the held-out stage AUROC, the `pooled_stage` capture and the population target where available, the paired differences for every cell, the `delta` metric note, and on the plasmode the per-block support sizes.
- **Documentation and records.** A dated report `docs/reports/production-latent-space-<date>.md`; aggregates under `results/production-latent-space-<date>/`; a new section of `docs/measurement-visibility.md` (which stage-supervised space, and what weighting and sparsity change); roadmap item 5 (the check weighting owed) and a new roadmap item; `AGENTS.md`. The package version is unchanged unless package code is touched; the report pins the version it ran on.

**Non-goals.**
- **The paper grade and any default flip.** On `propose_paper_grade`, a separate change runs the named arm at 500 × 999 (the stage-span paper grade's cells and the plasmode's grade), and only its verdict may flip the harness and `motco analyze` defaults, with pins for every committed configuration and the SEA-AD re-issue in `sead_sexdiff`.
- **Joint-scope sparsity and `mfa` weighting as arms.** The factorial names one kind of each; the alternatives are measured only if the named kind fails and a follow-up asks for them.
- **Out-of-span questions.** The split interaction test and the scan own them.
- **A padded or higher-dimensional simulator.** The p ≫ n regime is measured on the real export, not by adding noise features to the generator.
- **Label-aware spaces refitted per permutation** (the comparison's open thread 2).

## Capabilities

### New Capabilities
- `production-latent-space-selection`: the paired check of the production stage-PLS space against its weighted and sparse variants: the arms, both profiles and their reproduction of the committed parents, the real-noise span plant, the unit record, the predeclared decision with its descriptive tables, and the dated report.

### Modified Capabilities
None. The arms are the settings `block-weighting` and `regularized-stage-pls` define; nothing in the package changes behaviour.

## Impact

- **Scripts:** `scripts/production_latent_space_check.py` and `scripts/production_latent_space_array.sbatch`. The script imports `latent_space_comparison` (designs, seeds, `unit_dataset`, `cell_spec`, settings helpers) and `plasmode_check` (`load_export`, label permutation, seeds, the analysis-parameter builder, the fidelity values) unchanged; the real-noise span plant is its own function with tests.
- **Tests:** `tests/test_production_latent_space_check.py`: unit counts (3,000; 400), the shared cells' seeds and datasets equal to the comparison's, the four arms' parameters, the span plant (stage means of the planted group stretched by exactly 1 + e about the centroid in standardized units; the anchor untouched at e = 0; `Female` rows untouched), the fidelity gate, settings binding and parent hashes, resume, merge tables, both reproduction outcomes, every decision path including precedence and the gate.
- **Compute:** about 3,000 simulation units at roughly two production CV fits plus two sparse CV fits (about five production fits each) per unit, and 400 plasmode units at the same four fits on 80 × 10,000; the estimate before the smoke is 180–260 core-hours on `ing`. A local smoke replaces the estimate, the ladder length and the plasmode replicate count are the levers if it overruns, and launch waits on the user's approval with a cap.
- **Dependencies on other changes:** `block-weighting` (0.21.0) and `regularized-stage-pls` (0.22.0) must be merged first; the check runs on the version that carries both.
- **Compatibility:** no package output changes; the parents' committed files are read, hashed and refused by name if changed.
