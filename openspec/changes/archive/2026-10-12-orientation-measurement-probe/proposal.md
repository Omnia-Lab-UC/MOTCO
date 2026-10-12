# Proposal

## Why

At the Phase 6 SEA-AD design (n = 80, three stages), orientation/`angle` power is 0.300 at the strongest
construction the generator can build. The Phase 6 records show that the generator and sampling are not what
limits it.

- **The rotation is generated exactly and survives sampling.** The orientation surgery realizes exactly
  arccos(1 − e) in the population: 41.4° / 60° / 75.7° / 90° at e = 0.25 / 0.50 / 0.75 / 1.00. In the observed
  standardized space it is still about 90° at e = 1.
- **The latent space collapses it.** The stage-supervised, group-blind PLS space reduces that 90° to a median
  of 32° when CV keeps rank ≤ 3, against an anchor median of 18°. When CV happens to keep more than three
  components, the median is 76°, and power is 0.506 instead of 0.082. A space fitted to pooled stage means
  spans the symmetric combination of the two groups' stage subspaces, so a group-specific direction is
  averaged away.
- **The statistic has a heavy null tail.** The signed PC1 angle has a 95th percentile of about 137° at the
  anchor, in both spaces. This comes from PC1 instability in near-isotropic three-stage configurations, where
  PC1 can swap with PC2 or flip sign. On an across-replicate oracle threshold, the signed observed-space angle
  separates e = 1 from the anchor in 0.000 of replicates. Its folded counterpart separates 0.996.
- **Folding is not an option.** Orientation is the *direction* of change, not an axis. The candidate
  alternative estimator is the net displacement s_K − s_1. It equals the transition vector at K = 2, as the
  PC1 estimator does, and it stays well-defined when the configuration is isotropic.

A ρ sweep was considered first and is dropped. The Phase 5 design-point pilot already showed that baseline
continuity does not help, and the evidence above locates the loss in the measurement, not the generator.

Which measurement design best recovers trajectory geometry is an open question the paper must answer;
stage supervision is a candidate, not a requirement. This change is a **probe**: an exploratory,
matched-seed comparison that decides which designs earn a paper-grade comparison study (a second, later
change).

## What Changes

- Add an opt-in **net-displacement orientation estimator** to `estimate_difference` and `RRPP`. The default
  stays the PC1 estimator, and default outputs are byte-identical.
- Add two **diagnostic integration methods** to the simulation evaluation harness. Both are probe-only, not
  production:
  - `pca`: label-free principal components of the standardized, concatenated selected blocks, at a declared
    rank.
  - `pls_within_group`: stage-supervised PLS fitted separately within each group. The groups' weight vectors
    are orthonormalized into a common basis spanning S_A ⊕ S_B.
- Add an `orientation_estimator` evaluation parameter. It is omitted from records and signatures while it
  holds the default.
- Record the latent realized-geometry checkpoint for every constructed latent space, not only `pls`.
- Add a probe driver, `scripts/orientation_measurement_probe.py`. It crosses five measurement spaces with two
  orientation estimators on matched seeds at two designs:
  - Spaces: stage-PLS with CV rank (production), stage-PLS at forced rank 2(K − 1), `pls_within_group`,
    `pca` at rank 2(K − 1), and `concat` (ceiling reference only).
  - Estimators: PC1 and net displacement.
  - Designs: Phase 6 SEA-AD (n = 80) and Phase 5 (n = 1200, four stages).
  - Cells: anchor, orientation e = 0.5 and 1.0, magnitude and shape e = 1.0. Each cell runs 100 replicates ×
    199 permutations.
- Add a predeclared, advisory **advancement rule**. It decides which space × estimator arms advance to the
  paper-grade comparison, and writes a machine-readable decision.
- Add an evidence script that reproduces the checkpoint and flip-tail tables from a study `merged.jsonl`.
- Run the probe and commit its results directory and a findings report. Update the roadmap.

## Capabilities

### New Capabilities
- `orientation-measurement-probe`: the matched-seed probe driver, its arms and cells, the advancement rule,
  its outputs and provenance, and the motivating-evidence reproduction.

### Modified Capabilities
- `trajectory-orientation-invariance`: adds the opt-in net-displacement orientation estimator alongside the
  PC1 estimator. The existing PC1 requirements are unchanged and remain the default.
- `simulation-evaluation-harness`: adds the diagnostic `pca` and `pls_within_group` integration methods, the
  `orientation_estimator` evaluation parameter, and the latent geometry checkpoint for non-PLS latent spaces.

## Impact

- **Code:**
  - `src/motco/stats/trajectory.py` and `src/motco/stats/permutation.py`: estimator keyword, passed through
    the RRPP worker.
  - `src/motco/simulations/evaluation.py`: new integration methods, estimator parameter, and latent
    checkpoint wiring.
  - `src/motco/simulations/diagnostics.py`: latent scope for non-PLS spaces, if needed.
  - New `scripts/orientation_measurement_probe.py` and `scripts/orientation_checkpoint_evidence.py`.
- **Compatibility:** every default path is byte-identical.
  - Existing evaluation records, cell ids, signatures and committed configs do not change;
    `tests/test_committed_config_identity.py` must stay green.
  - The new methods are not accepted by the production forced-rank guard or the study engine's acceptance
    logic. Wiring them into the study engine belongs to the follow-up change.
- **Compute:** about 2,500 evaluations per design (5 cells × 100 replicates × 5 spaces; the two estimators
  share each space fit and run separate RRPP passes, about 5,000 RRPP runs per design). Expect a few
  core-hours at n = 80 and a few tens at n = 1200. Runs on `ing`.
- **Docs:** `docs/reports/orientation-measurement-probe-<date>.md`, a roadmap entry, and a CLAUDE.md
  architecture note for the new methods and the estimator.
