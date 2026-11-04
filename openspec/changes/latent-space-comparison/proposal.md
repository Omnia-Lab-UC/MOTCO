## Why

The stage-span pilot (2026-11-04) established that the production measurement sees a group difference only through its component inside the stage-supervised PLS span. Given a difference it can see, the tests are sensitive. That blindness is the scope of the space's question ("does each group traverse the shared stage program the same way?"). It raises the next question: does a different space see more, and what does it cost in size and in power on what stage-PLS already sees?

The repository has compared spaces once. The orientation probe ran label-free `pca` at rank 2(K − 1) against one construction, a relocation surgery whose out-of-span share was not controlled. Comparing spaces on any single construction is circular. A plant inside the span favours stage-PLS, and a plant outside it favours a label-free space. Which construction is "true" is not knowable from simulation: how real group differences sit relative to the stage program is the unknown.

This change parametrizes that unknown, the **alignment** of the difference with the stage span, and measures every candidate space across it. It runs at the simulator's stage signal and at one calibrated to the SEA-AD cohort's. The result is a visibility map for each space rather than a ranking that the construction decides.

## What Changes

- **Alignment axis in the stage-span construction.** A new opt-in `span_alignment` (φ, in degrees, 0–90) applies to the stage-span `orientation`.
  - The rotation's plane tilts from the in-span direction v toward a unit direction o on features the null stage program does not use: w_φ = cos φ·v + sin φ·o.
  - The intended population geometry is the same at every φ (`angle` θ, `delta` and `shape` 0). Only where the difference points changes.
  - φ = 0 is the current stage-span rotation, byte-identical. At φ = 90° the difference lies outside the pooled stage span: it is the controlled analogue of the relocation surgery.
  - A companion `span_support` sets o's support: dense over every unused feature (the default), or k of them.
  - o is drawn from a separate seeded stream. A dataset therefore stays the null dataset plus the plant, and o is shared across φ and effect at a replicate.
  - Both fields are hash-omitted while absent.
- **A new diagnostic space, `pls_residual_pca`.** It takes the production stage-PLS scores and appends r principal-component scores of the stage-residualized features, orthogonal to the PLS span (default r = K − 1). Its basis uses the stage labels and the features, never the group labels. It is diagnostic, never production.
- **Span capture for every linear space.** `span_capture` and the `population_latent` checkpoint extend from `pls` to `pca`, `pls_residual_pca` and `concat`. Each space records its own orthonormal basis; for `concat` the basis is the identity and the share is 1. PLS records are unchanged.
- **Comparison driver.** `scripts/latent_space_comparison.py` (`calibrate` | `run` | `merge` | `decide`, with an sbatch array) runs five arms per unit on one dataset with one RRPP seed:
  - `pls_cv`, the production measurement;
  - `pca_matched`, PCA at that unit's `pls_cv` rank;
  - `pca_fixed`, PCA at rank 2(K − 1), the probe's rank;
  - `pls_residual_pca`;
  - `concat`, the reference ceiling.

  The run uses three designs at three stages on methylation + expression: n = 200 at the generator's stage effect δ = 2.0; the SEA-AD table at δ = 2.0; and the SEA-AD table at a δ calibrated so the anchor's held-out stage AUROC matches SEA-AD's (0.636). It runs 14 cells per design (listed in design D10) at pilot grade, 100 replicates × 199 permutations, with base seed 2035.
- **Advisory reading.** `latent_space_comparison_reading.json` reads each eligible arm on three criteria, each against `pls_cv`:
  - size at the anchor and the translation control, at every design;
  - out-of-span gain: orientation `angle` power at φ = 90°;
  - no loss on protected statistics: in-span orientation `angle`, `shape` and magnitude `delta`.

  The verdict is `propose_paper_grade` when at least one arm qualifies, else `document_tradeoff`. Descriptive tables give each arm's capture and power against φ, and a sparse-versus-dense contrast that tests the visibility page's claim that blindness is a direction, not a set of features.
- **Report and documentation.**
  - **Report:** a dated report under `docs/reports/` with aggregates under `results/`.
  - **`docs/measurement-visibility.md`:** gains each space's visibility profile, plus two corrections:
    - The rhombus extends to three or more stages only when B's program moves onto features A's does not use. A rotation within A's own features stays partly or fully visible.
    - A projection cannot see a difference `concat` lacks, but it can raise power by discarding noise directions.
  - **Roadmap:** gains item 14.
  - **Version:** MOTCO 0.16.0.

**Non-goals.** Each is a separate change:

- any change to the production measurement, its defaults or `motco analyze`;
- SNF and graph-native statistics;
- block weighting (roadmap item 5);
- a data-driven PCA rank rule;
- the paper grade, which this pilot's reading proposes;
- a cross-fitted estimate of where a real cohort sits on the alignment axis;
- the alignment axis for `shape` (design D3).

## Capabilities

### New Capabilities
- `latent-space-comparison`: the multi-arm comparison of candidate measurement spaces across alignment. It covers the arms and the designs, the stage-signal calibration, matched seeds, the advisory reading and the dated report.

### Modified Capabilities
- `semisynthetic-trajectory-generator`: the stage-span orientation gains `span_alignment` and `span_support`.
- `simulation-evaluation-harness`: the diagnostic `pls_residual_pca` integration, and span capture extended from PLS to every linear integration method.
- `realized-geometry-diagnostics`: the `population_latent` checkpoint and the `span_capture` shares are defined for any linear basis, not only the PLS rotations.
- `latent-space-metric-compatibility`: the visibility documentation states the rhombus's condition at three or more stages, the information-versus-power reading of `concat`, and each candidate space's visibility profile.

## Impact

- **Code:**
  - `src/motco/simulations/semisynthetic.py`: the two parameters, the out-of-span direction and its stream, and truth;
  - `src/motco/simulations/evaluation.py`: `pls_residual_pca`, a linear basis on every linear integration result, and span capture for those methods;
  - `src/motco/simulations/diagnostics.py`: span capture from a basis;
  - `src/motco/simulations/grid.py`: `_SIGNATURE_OMIT_WHEN_NONE` gains `span_alignment` and `span_support`.
- **Scripts:** `scripts/latent_space_comparison.py` and `scripts/latent_space_comparison_array.sbatch`.
- **Tests:**
  - byte identity for φ absent and φ = 0 and for PLS span capture;
  - the geometry invariance across φ;
  - the out-of-span property at φ = 90°;
  - `pls_residual_pca` validity (group-label invariance, orthogonality);
  - the driver's per-arm identity with direct evaluations;
  - the reading's verdict paths.
- **Compute:** 3 designs × 14 cells × 100 replicates = 4,200 units of five arms on `ing`. That is an estimated 50 core-hours, from the stage-span pilot's 17 s per single-arm unit and the orientation probe's per-arm cost. The smoke timing replaces the estimate before launch.
- **Compatibility:** every committed configuration, fixture and dated report reproduces byte-identically. The new fields are hash-omitted while absent, and `pls` span capture records are unchanged.
- **Docs:** `docs/measurement-visibility.md`, `docs/api/simulations.md`, `docs/roadmap.md`, `AGENTS.md` and the dated report.
