# Design

## Context

See proposal.md (Why) for the evidence. The relevant current state:

- `estimate_difference` computes `angle` through `_estimate_orientation`, which returns PC1 signed by net
  displacement. `RRPP` calls `estimate_difference` per permutation through the picklable `_RRPPWorker`.
  `diagnostics.geometry_from_means` calls `estimate_difference` for every realized-geometry checkpoint.
- The harness builds one latent matrix per evaluation (`concat`, `snf`, `pls`) and runs RRPP on it, with no
  refit per permutation. The `pls` latent checkpoint is computed only because `latent_matrix` is passed when
  `integration_method == "pls"` (`evaluation.py:264`).
- Parameter signatures serialize every dataclass field. `grid._SIGNATURE_OMIT_WHEN_NONE` drops listed fields
  while they are `None`, which is how `group_stage_sizes` was added without moving committed cell ids.
- Probe precedent: `scripts/latent_rank_probe.py`. It holds generated data fixed, varies one measurement
  knob, and writes a CSV and Markdown file under `results/<probe>-<date>/` with `PROVENANCE.txt`.

## Goals / Non-Goals

**Goals:**
- Every new path is opt-in. Default behavior, records, signatures and committed configs stay byte-identical.
- The arms are compared paired: same datasets and same permutation draws, so differences between arms are
  measurement differences only.
- The new methods and the estimator are library features (harness parameters), not probe-local code, so the
  follow-up study change can wire them into the study engine without reimplementing them.

**Non-Goals:**
- No study-engine support for the new methods: no config schema, design-grid axes, report columns or
  acceptance targets. That belongs to the follow-up change.
- No per-permutation refit of `pls_within_group`. The probe measures Type I empirically, and refitting is the
  fallback if the anchor inflates.
- No PCA rank rule. The probe fixes the rank to `2(K − 1)`, so the space type is not confounded with a rank
  rule.
- No p-padding (null features) axis, and no changes to production defaults or to the attribution path.

## Decisions

### D1. Net displacement is a selectable estimator inside `estimate_difference`, not a separate statistic

`estimate_difference(..., orientation_estimator="pc1" | "net_displacement")` dispatches between
`_estimate_orientation` and a new `_estimate_net_displacement`. The new function uses the last minus the
first LS-mean row of each group's `levels`, in the order `build_ls_means` produces. `RRPP` takes the same
keyword and stores it on `_RRPPWorker`. `geometry_from_means` takes it too, so the checkpoints measure the
same estimand as the test.

- *Alternative: a fourth statistic next to `angle`.* This would compute both estimators in one RRPP pass.
  It was rejected because it would widen every return tuple and record schema in production code for a probe
  question. The cost is a second RRPP pass per space, which repeats `delta`/`shape` work; that is acceptable
  at 199 permutations.
- **Stage order.** Stages are ordered by `build_ls_means`'s sorted string representation. The harness's
  numeric-aware level ordering already gives "0" < "1" < … for both designs (K ≤ 4). The test suite pins
  this for the estimator.
- **Zero length.** A zero-length displacement returns NaN orientation. The existing availability handling
  then reports `angle` unavailable, as it does for a degenerate PC1.

### D2. The estimator selection is `None`-defaulted and normalized

`SimulationEvaluationParams.orientation_estimator: str | None = None`. `None` means PC1. Validation
normalizes an explicit `"pc1"` to `None`, so both spellings serialize identically. The field is added to
`_SIGNATURE_OMIT_WHEN_NONE`, and the record writes it only when it is not `None`.

- *Alternative: default `"pc1"` plus a separate "omit when equals default" list.* Rejected: it adds a second
  omission mechanism where one already exists.

### D3. `pca` fits on exactly the matrix `concat` produces

`pca` takes the `concat` outcome matrix (M-value methylation, pooled standardization, selected blocks) and
returns the scores of a plain SVD-based PCA at `integration_params["n_components"]` (required, 1 ≤ k ≤
min(n, p)). There is no extra scaling, because the blocks are already standardized. Sharing the `concat`
matrix makes `concat` versus `pca` a pure "reduce or not" contrast.

### D4. `pls_within_group` uses pooled standardization and per-group weight vectors

For each group, a PLS-DA on stage is fitted with `scale=False` to that group's rows of the pooled
standardized matrix, at `integration_params["components_per_group"]` (default `n_stages − 1`). Its
`x_rotations_` are taken, both groups' rotations are stacked, and the stack is orthonormalized by QR with a
rank tolerance. The pooled matrix, centered by its pooled mean, is projected onto the resulting basis.

- **Why `scale=False`.** Per-group rescaling would give the two groups different metrics and make the two
  subspaces incomparable.
- **Why pooled centering.** The projected scores of both groups then share an origin, so translation between
  groups stays a translation.
- **Fixed per-group rank, not CV.** With 9–28 samples per group × stage cell at n = 80, a nested CV inside
  each group is fragile. A fixed `n_stages − 1` matches what the pooled CV selects at both designs.
- **Under the null,** the two groups' subspaces nearly coincide, and the remaining dimensions are noise. That
  over-retention cost is exactly what the anchor and specificity criteria are there to catch.
- *Alternative: PLS supervised on the group × stage cell label.* Rejected for the probe. Its axes directly
  separate the groups, so it would need a per-permutation refit to be valid, which multiplies cost by the
  permutation count.

### D5. The latent checkpoint is passed for the new methods only

The harness passes `latent_matrix` to `calculate_realized_geometry` for `pls`, `pca` and `pls_within_group`.
The scope key stays `pls_latent` for `pls`; the new methods use `pca_latent` and `pls_within_group_latent`.
Records for `snf` and `concat` evaluations do not change.

### D6. The probe driver evaluates every arm per generated dataset

`scripts/orientation_measurement_probe.py` loops over design → cell → replicate. For each unit it:

1. generates the dataset once from the matched-seed derivation (the same replicate index gives the same
   generator seed across cells, as in the study engine);
2. evaluates the five spaces;
3. for each space, runs RRPP under both estimators with one shared RRPP seed.

The forced-rank PLS arm uses the existing `forced_components`. One JSONL line is written per (unit, arm).
Resume keys on (design, cell, replicate). A settings hash, stored in a header file, refuses resumes with
different settings.

- **Parallelism.** Shard by (design, cell, replicate block) with a `--shard i/N` option, run as a SLURM array
  on `ing`. Each shard writes its own JSONL, and a `--merge` step builds the outputs. This avoids a
  workstation run at n = 1200.
- **Cost.** One RRPP pass per estimator per space: 2 designs × 5 cells × 100 replicates × 5 spaces × 2
  estimators = 10,000 evaluations' worth of RRPP at 199 permutations. CV runs only in the `pls` CV arm. A
  smoke run times one unit per design before launch, and the launch budget is recorded in `PROVENANCE.txt`.

### D7. The decision is computed from the rate table, not from records directly

The driver first builds the rate table: rate, SE and count per design × cell × arm × statistic. The
advancement rule (spec) reads only that table and writes `advancement_decision.json` (per arm: criteria,
values, thresholds, pass/fail, and an `advances` flag) plus a CSV mirror. Keeping the rule a pure function of
the table makes it unit-testable on synthetic tables.

### D8. The motivating-evidence script reads committed-format records

`scripts/orientation_checkpoint_evidence.py --merged <merged.jsonl> --column baseline` reproduces the
checkpoint and flip-tail tables (spec) from `realized_geometry.checkpoints` and
`integration_metadata.selected_lv`. Column selection uses `cell_metadata.design_point`, mirroring the
exploration snippet. The script is run on the Phase 6 orientation profile, and its output is quoted in §1 of
the findings report.

## Risks / Trade-offs

- **[The Phase 5 arm under joint magnitude is not the Phase 5 paper-grade run]** → The reference arm's own
  rates are recomputed in the probe, and every comparison is against them. The deviation is stated in the
  outputs and the report.
- **[100 replicates give a pooled SE of about 0.06 near 0.5]** → The 2-SE gain criterion detects roughly
  ≥ 0.17 improvements. Smaller effects are left to the paper-grade follow-up. The probe ranks arms; it does
  not certify them.
- **[`pls_within_group` may inflate Type I without a refit]** → The anchor criterion fails the arm, and the
  report recommends either a refit variant or dropping the arm for the follow-up.
- **[PCA at 2(K − 1) may miss the trajectory, if stage signal is not high variance under the InterSIM
  covariance]** → The latent checkpoint per arm shows whether the rotation reached the space, which separates
  "the space missed it" from "the test could not detect it".
- **[`concat` looks best at p = 498 but would degrade at real p]** → It is reference-only by rule. The
  p-padding axis is deferred to the follow-up explicitly.
- **[Net displacement ignores interior stages, so a shape change at the endpoints would leak into it]** →
  The shape cell's `angle` rate under each estimator is part of the specificity criterion.

## Migration Plan

Nothing to migrate: every change is additive and opt-in. Rollback is a revert. Committed configs and their
cell identities are guarded by `tests/test_committed_config_identity.py`, which must pass unchanged.

## Open Questions

- Whether the follow-up study keeps `concat` as a ceiling or replaces it with a p-padded `concat`. This is
  decided from the probe report, and changes nothing here.
