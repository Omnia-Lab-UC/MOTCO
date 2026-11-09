## Why

The trajectory statistics are measured in PLS score coordinates, `T = X_in·R` with `R = W·(PᵀW)⁻¹`. The weights `W` are orthonormal and `PᵀW` is unit upper triangular, so `T = (X_in·W)·U`, where `U` is unit upper triangular. That is a **shear**: each score axis after the first leans toward the earlier ones. The scores cover the same subspace as `W`, but they do not preserve the lengths and angles of the standardized feature space inside it.

Size is unaffected. Under the null both groups share one population trajectory, every linear map keeps two equal configurations equal, and RRPP permutes in the same coordinates. Specificity is affected: a difference that is pure in the standardized features can register in another statistic. The stage-span paper grade's leak decomposition (2026-11-04) measured this as the B − A path-length difference for a plant whose true difference is exactly 0:
- **The bend at e = 0.5:** −0.48 % of A's path on the SEA-AD table and −0.59 % at four stages in scores, against −0.01 % and +0.06 % on an orthonormal basis of the same span.
- **The 45° rotation on the SEA-AD table:** −3.35 % in scores against +0.66 % on the orthonormal basis.
- **The unexplained cell:** the rotation's off-target `shape` rises steadily with the angle, to 0.070 at n = 200 and 0.086 at four stages at 45°, at the 0.0695 per-cell bound. The plant keeps `shape` exactly in the planted coordinates, so the leak enters through the measurement. The shear is the candidate the decomposition supports, but the run did not isolate it.

This is the one route the construction cannot remove (roadmap item 15). The other linear spaces in the harness (`pca`, `pls_within_group`, `concat`) already measure in orthonormal coordinates, so only the production space is sheared. The fix is to measure on `X_in·Q`, with `Q` an orthonormal basis of the fitted span. `delta`, `angle` and `shape` are invariant to orthogonal maps, so every orthonormal basis gives the same statistics. The measurement then depends only on the fitted subspace and not on how PLS parametrizes it. It also equals the geometry that the stage-span constructions and the documentation already define.

It changes the production measurement, so it must first be measured against the scores on matched data under a predeclared decision. The change should come before the paper's numbers are written: adopting it afterwards would mean re-issuing them.

## What Changes

- **Opt-in orthonormal coordinates in the harness.** `SimulationEvaluationParams.latent_coordinates` takes these values:
  - absent, the default: today's scores;
  - `"scores"`: normalized to absent;
  - `"orthonormal"`.

  Under `"orthonormal"`, the PLS latent matrix is the model's input coordinates `X_in = (X − x_mean)/x_std` projected on `Q = qr(x_rotations_)`. This is the basis `LinearLatentMap.from_pls` already uses, and it equals `W` up to column signs. Nothing else moves:
  - the fitted model, the cross-validation, the selected rank and the forced-rank path are unchanged;
  - the latent matrix, the `pls_latent` checkpoint, the configuration spectrum, the statistics, RRPP and the `population_latent` checkpoint are all computed in the new coordinates;
  - the span-capture shares are unchanged, because they depend only on the subspace.

  It applies to `pls` integration only. With any other method, or with attribution, it is rejected before any fit. While absent, it is left out of cell identifiers, signatures, records and runtime metadata (`grid._SIGNATURE_OMIT_WHEN_NONE`), so every committed configuration, fixture and report is unchanged.
- **A prebuilt-fit path.** `pls_orthonormal_from_pls(...)` builds the orthonormal integration from an existing `pls` integration without refitting, equal to fitting it from scratch. A validation unit therefore fits one PLS model and measures it twice.
- **A validation driver.** `scripts/orthonormal_coordinates_check.py` (`run --profile paper_grade|weak` | `merge` | `decide`, sbatch array) runs two arms per unit, on one dataset, one PLS fit and one RRPP seed: `scores`, the production measurement, and `orthonormal`. Its profiles are:
  - **`paper_grade`** rebuilds the stage-span paper grade's units from its matched seeds. That is 65 cells × 500 replicates × 999 permutations, under the planted scale, at its five columns: the anchor, translation, orientation and shape at e ∈ {0.1, 0.25, 0.5}, and magnitude at {0.02, 0.05, 0.1}. The `scores` arm must reproduce the committed `design_point_operating.csv` rejection counts exactly.
  - **`weak`** rebuilds the latent-space comparison's `seaad_weak` anchor and translation units (δ = 0.11, held-out AUROC 0.617; 2 cells × 500 × 999). This tests size at the stage signal where real data sits. The `scores` arm must reproduce the span-split check's committed `production_table.csv` counts exactly.
- **A predeclared decision.** `orthonormal_coordinates_decision.json` (rule version 1) judges the `orthonormal` arm. It reuses the paper grade's count-plus-ceiling bound rule:
  - **O1, size:** anchors and translation at the five columns (m = 60), and the weak design's two controls (m = 6);
  - **O2, specificity:** the off-target rates of the stage-span orientation, shape and magnitude cells (m = 90);
  - **O3, power non-inferiority:** the target rate of each of the 45 stage-span target cells against `scores`, by an exact paired (McNemar) test under a count-plus-ceiling rule.

  The verdict is `adopt_orthonormal` | `size_excess` | `specificity_excess` | `power_loss`, with `reproduction_failed` gating all of them. Descriptive tables give the paired off-target and target differences per cell, with the rotation's `shape` named. They also give the population-level path-length residual in each arm, and the changes in eigengap and `angle` null width.
- **Documentation and records.**
  - **`docs/measurement-visibility.md`:** the score-coordinate section gains the measured comparison and the verdict.
  - **API page:** `docs/api/simulations.md`.
  - **Report:** a dated report, with aggregates under `results/orthonormal-coordinates-<date>/`.
  - **Roadmap and repository guide:** roadmap item 15 and `AGENTS.md`.
  - **Version:** MOTCO 0.19.0.

**Non-goals.** Each of these is a separate change, and each follows only on `adopt_orthonormal`:
- flipping the harness default, with pins for every committed configuration and diagnostic;
- `motco analyze` and `TrajectoryAnalysisParams`, with a real-noise size check at the case study's configuration and the SEA-AD re-issue;
- attribution under orthonormal coordinates, whose PLS-captured reconstruction is an oblique projection onto the loadings' span and would need an orthogonal counterpart;
- the study-configuration key and design axis.

Also out of scope: the surgery profiles, whose plants keep 0.2–0.8 % of their energy in the span, and `pls_residual_pca`.

## Capabilities

### New Capabilities
- `orthonormal-coordinates-validation`: the validation of orthonormal span coordinates against the production scores. It covers the two profiles and their reproduction of the parent runs, the paired two-arm unit, the predeclared decision and its descriptive tables, and the dated report.

### Modified Capabilities
- `simulation-evaluation-harness`: the PLS measurement gains the opt-in `latent_coordinates` setting, with its normalization, rejection rules, invisibility while absent and the prebuilt-fit path. The single-pass requirement is restated so the measurement matrix is the estimator's training scores or, under orthonormal coordinates, its training inputs on the orthonormal basis of the same estimator.
- `latent-space-metric-compatibility`: the score-coordinate documentation states the measured comparison and the decision, replacing "not validated".

## Impact

- **Code:**
  - `src/motco/simulations/evaluation.py`: the parameter, its validation, the coordinate step after the final fit, the prebuilt-fit helper, and the linear map's `transform` in the measurement's coordinates;
  - `src/motco/simulations/diagnostics.py`: a `LinearLatentMap` constructor for the orthonormal coordinates;
  - `src/motco/simulations/grid.py`: the omit set;
  - `src/motco/__init__.py`: 0.19.0.

  `motco.analysis` and `motco.stats` are not touched.
- **Scripts:** `scripts/orthonormal_coordinates_check.py` and `scripts/orthonormal_coordinates_array.sbatch`. They import `stage_span_leak_check` (column cells, unit parameters), `stage_span_paper_grade_decision` (`judge_rates`) and `span_split_check` / `latent_space_comparison` (the weak design's units) unchanged.
- **Tests:**
  - the orthonormal matrix spans the scores' subspace and equals `X_in·W` up to column signs;
  - statistics are invariant to the choice of orthonormal basis;
  - normalization and rejections;
  - byte identity while absent, including committed cell ids and signatures;
  - prebuilt-fit equality;
  - driver reproduction, resume and settings binding;
  - every decision path.
- **Compute:** about 33,500 units (32,500 paper-grade, 1,000 weak), one PLS fit and two RRPP passes each. The paper grade ran at about 17 s per unit (229 core-hours for 48,000 units), so the estimate is 170–230 core-hours on `ing`. A local smoke replaces the estimate before launch, and launch waits on the user's approval.
- **Compatibility:** every committed configuration, fixture, record and report is unchanged. The setting is absent by default, and nothing outside the harness reads it.
