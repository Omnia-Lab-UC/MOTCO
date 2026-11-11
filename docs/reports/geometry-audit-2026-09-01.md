# Geometry audit: statistics, constructions, and documentation

**Run date:** 2026-09-01
**Code revision:** `cd583a7` (branch `main`, the merge of `diagnose-angle-null-pivotality`)
**Records read:** `results/angle-pivotality-2026-09-01/merged.jsonl` (500 replicates) and
`results/phase4-2026-08-27/merged.jsonl` (1,900 units) — read-only; no new study cells were run
**New computation:** three targeted experiments — a deterministic estimator check, and replicate
regenerations from persisted `generator_seed`s pushed through the identical PLS integration. Every
regenerated statistic was matched against its persisted record (angles to 1e-8, shapes to 4+ decimals)
before any conclusion was drawn.
**Scope:** `stats/` (trajectory, permutation, design, pls, snf, attribution), `simulations/`
(semisynthetic, generator, evaluation, preprocessing, diagnostics, specificity, pivotality, grid), all six
reports under `docs/reports/`, the OpenSpec spec set, and the test suite
**Affects:** [Phase 5 readiness](https://github.com/Omnia-Lab-UC/MOTCO/blob/c440070dcb3c5bcc6ddd2227cf7adbef18573b6a/docs/phase5-readiness.md) items 1–5
**Environment:** Python 3.11.15, numpy 2.3.5, pandas 2.3.3, scikit-learn 1.8.0, scipy 1.16.3; no R at runtime

## Summary

The statistical core — RRPP, the design coding, the Procrustes solver, the sign-anchor fix, the M-value
contract, the matched-seed and signature machinery — is sound, and the 2026-09-01 pivotality diagnosis is
correct as far as it goes. The audit found that it stops one step short, and found four further defects:

- **F1.** The per-replicate `angle` null width that the pivotality report calls *"invisible to every
  integration diagnostic currently persisted"* is in fact predicted (Spearman −0.81 on the 16 extreme
  replicates) by the **relative eigengap of the latent stage-mean configuration** — and that eigengap is
  near zero *by construction*, because the baseline draws independent per-stage indicator programs whose
  stage means form a near-regular simplex.
- **F2.** The orientation and translation effect axes are **silently censored** by the relocation clamp:
  the orientation surgery saturates at e ≈ 0.69, and **80 of 100 replicate pairs at e = 0.75 and e = 1.00
  are byte-identical datasets** with identical statistics. The flat top of the orientation power curve is
  partly guaranteed by construction.
- **F3.** The `shape` reflection contract is **ambient-dimension-dependent**: "reflections retained as
  distinct" is true only where the trajectory configuration spans the full space (the PLS latent
  checkpoint) and silently vacuous at every pre-integration checkpoint, so the Phase-4 localization table
  compares two different statistics across its rows. Measured on regenerated pilot replicates, however,
  reflection is **ruled out** as the orientation→shape mechanism (allowing reflections changes the latent
  shape distance in 0 of 100 replicates) — the latent shape response is genuine projection-induced
  deformation.
- **F4.** The tested `angle` statistic (PC1-axis divergence) and its attribution machinery
  (per-adjacent-transition contrasts) decompose **different orientation quantities** for ≥ 3 stages; the
  planned "signed principal orientations" extension was never implemented.
- **F5 (synthesis).** The load-bearing assumption *surgery mode ⇒ target latent statistic* fails at five
  compounding layers, and `angle` is the only statistic exposed at all five — which is why every fix so far
  (sign anchor, studentization candidates, more replicates) moved it by points rather than tens of points.
  The deepest layer is the independent-indicator baseline, which makes the orientation estimand itself
  nearly degenerate per draw.

Plus a set of documentation and spec defects (D1–D9), including README examples that raise on execution
and a stale GPA description in `CLAUDE.md`.

A remediation plan, expressed as proposed OpenSpec changes in priority order, closes the report.

---

## F1 — The angle-null width has a measurable driver

The pivotality report established that the RRPP null for `angle` tracks its own observed statistic (slope
0.811), that the tracking is load-bearing, and that *"nothing the harness records today identifies which
replicates will be resolvable"* — selected dimensionality and CV AUROC both saturate and neither correlates
with the null width. That last gap is closed here: the missing quantity is the **relative eigengap
(λ₁−λ₂)/Σλ of the centered stage-mean configuration** in the latent space, computable per replicate (and
per permutation) from matrices already in hand.

### Evidence: recomputed from the committed pivotality records

The 8 narrowest-null and 8 widest-null orientation replicates were regenerated from their persisted
`generator_seed`s and pushed through the identical PLS integration (`cv1_splits=3, cv2_splits=4,
n_repeats=5, max_components=20, random_state=1203`). Fidelity first: every recomputed observed angle
matched the persisted `pair_statistics["angle"]` (e.g. 8.56 vs 8.56, 138.46 vs 138.46), so this is the
study pipeline reproduced bit-for-bit, not an approximation.

Per replicate, the relative eigengap of the pooled latent stage-mean configuration against the recorded
null width:

| Group | mean eigengap (pooled config) | mean recorded `angle` q95 |
|---|---:|---:|
| 8 narrowest nulls | 0.097 | 6.6° |
| 8 widest nulls | 0.046 | 174.0° |

- Spearman(log q95, pooled eigengap) = **−0.81**; Pearson on log–log = **−0.88**. (Correlations computed
  within the 16 selected extremes; a full-cell run over all 100 replicates confirms the association — see
  the addendum at the end of this section.)
- A 2× separation in eigengap maps to a **26×** separation in null width.
- The replicates whose *per-group* eigengap is near zero (0.002–0.013) are exactly the ones carrying huge
  observed angles (138°–157°). Observed statistic and null co-vary because both are driven by the same
  anisotropy — that shared dependence **is** the measured 0.811 slope, mechanistically.

### The construction sets the eigengap near zero

`_baseline_methyl` draws **independent Bernoulli(p_dmp) indicators per stage** ("intentionally not forced
continuous", per the module docstring). Independent stage programs of equal expected size put the stage
means at the vertices of a **near-regular simplex**: every pairwise distance is statistically identical, the
configuration is nearly isotropic, and PC1 — the orientation estimator — is decided by the few-percent
random asymmetry of the draw.

Direct measurement (200 indicator draws at k = 4 stages, p_dmp = 0.2): median relative eigengap **0.026**
(q5 0.005, q95 0.050), against **1.0** for a straight trajectory; net displacement is ~⅓ of path length.
PC1's estimator variance scales like noise/eigengap, so per-replicate angle nulls spanning 5°–177° are the
*expected* behavior of this construction. The heavy near-180° null mass (21/700, unchanged by the sign fix)
follows from the same degeneracy: once the axis itself wanders, the sign anchor cannot stabilize the angle.

### Consequences

1. **Readiness items 1 and 3 asked for exactly this quantity.** It is one SVD of a k×d matrix per
   configuration, computable per replicate *and per permutation*, and on real data it doubles as an
   "orientation resolvability" qualifier for reporting.
2. **The item-4 design-point question changes shape.** More samples shrink the noise on stage means, so
   power does rise with n — but the required n scales like the inverse square of the eigengap, and the
   gap's lower tail (draws near perfect isotropy) is what caps the power curve. Sample size buys those
   replicates almost nothing; the baseline's continuity (F5) is the stronger lever.
3. The spec expectation that *"null configurations never approach antiparallel"*
   (`trajectory-orientation-invariance`) is enforced only by a small-noise 2-D test
   (`tests/test_trajectory_orientation.py:109`). The production operating regime violates its spirit
   lawfully, because noise there is not small relative to the anisotropy.

### Addendum: full-cell confirmation

The same computation over **all 100 orientation replicates** (no selection): Spearman(log q95, pooled
eigengap) = **−0.75** (Pearson log–log −0.62), so the association is a property of the cell, not of the
extremes. The eigengap also predicts the *outcome*: replicates that reject `angle` average gap **0.053**,
replicates that fail average **0.035**. The observed angle tracks the smaller per-group gap (Spearman
−0.33), consistent with the shared-anisotropy mechanism. The same run carries the F3 reflection probe over
the whole cell: allowing reflections changes the latent shape distance in **0 of 100** replicates (mean
improper/proper ratio 1.000), and the reflection-allowed observed shape still exceeds the recorded
shape-null q95 in **99 of 100**.

---

## F2 — The orientation and translation effect axes are censored

`_relocate_rows` clamps the relocation count to the destination pool:
`k = min(round(e·|src_pool|), |dst_pool|)`. With p_dmp = 0.2 and four stages, ~59% of CpGs are stage-active
and only ~41% are inactive, so the orientation surgery **saturates at e ≈ 0.69**. The translation mode's
candidate pool (stage-inactive CpGs on untouched genes) is smaller still.

### Evidence: the Phase-4 records' own truth metadata

| Mode | e = 0.25 | 0.50 | 0.75 | 1.00 |
|---|---:|---:|---:|---:|
| orientation, realized relocations (mean) | 54.4 | 108.9 | **146.9** | **149.3** |
| orientation, nominal (uncensored) | 54 | 108 | 163 | 217 |
| translation, realized set size (mean) | 18.0 | 34.3 | **37.7** | **38.2** |
| translation, nominal | 18 | 37 | 55 | 73 |

Because matched seeds share the generator seed across effect sizes and the clamp makes k identical whenever
it binds, **80 of 100 replicate pairs at e = 0.75 and e = 1.00 are byte-identical datasets with identical
angle statistics** (same `generator_seed`, angle equal to 1e-9). The Phase-4 "0.67 vs 0.68" comparison
between those cells differs on 20 replicates of fresh information.

The other modes are clean: shape/relocate realizes ≈ its nominal count at every effect
(18.2 / 36.5 / 54.6 / 72.8 against nominal 18 / 37 / 55 / 73 — its per-stage source pool is smaller than
the destination pool), and magnitude scales δ with no pool at all.

### Consequences

- The flat top of the orientation power curve is partly guaranteed by construction, and the preregistered
  monotonicity gate is near-vacuous between the top cells.
- "Power at effect 1.00" actually characterizes an e ≈ 0.69-equivalent construction.
- A Phase-5 grid inheriting this axis would duplicate ~80% of its compute between those cells while
  reporting them as independent measurements.
- The realized sizes *are* recorded (`truth_metadata.transform`) — they have never been surfaced in a
  report. This is the Phase-2 lesson ("report realized geometry, not requested labels") recurring one level
  further upstream, at requested-effect → realized-surgery.

---

## F3 — The reflection contract is rank-dependent; reflection is *not* the orientation→shape mechanism

The shape spec and the Procrustes audit state that *"reflections are retained as distinct shapes by
default."* That is true precisely when the trajectory configuration spans the full ambient space
(dimension ≤ k−1) and **silently false otherwise**: when the configuration is rank-deficient, the
determinant correction in `_proper_procrustes_distance` flips a null-space singular vector, which costs
nothing in the objective — the "proper rotation" constraint is vacuous and the estimator returns the
reflection-allowed distance.

### Evidence: the production estimator on deterministic configurations

A bent 4-stage planar trajectory vs its mirror image, embedded by zero-padding:

| Ambient dimension | shape distance, config vs mirror |
|---:|---:|
| 2 | 0.431 |
| 3 | **0.000** |
| 10 | **0.000** |
| 660 | **0.000** |

A rank-3 configuration vs its mirror: positive in 3-D, 0.000 in ≥ 4-D. So at every pre-integration
checkpoint (`population_native`, `population_standardized`, `observed_standardized`: hundreds of
dimensions, configurations of rank ≤ 3) shape is reflection-**invariant**, while at `pls_latent` (3 dims,
rank-3 configurations) it is reflection-**sensitive**. The invariance test pins only the 2-D full-rank
case, so the suite passes while the contract fails in the very spaces the realized-geometry checkpoints
use. The Phase-4 localization table therefore compares two different statistics across its rows.

### Evidence: but reflection does not explain orientation→shape

The natural suspicion — that orientation→shape "first becomes material at the PLS checkpoint" because a
mirror-like difference is invisible pre-integration and visible in 3-D — was tested on the regenerated pilot
replicates. Allowing reflections in the latent space changed the shape distance by **0.0% in 100 of 100**
(the optimal alignment is already a proper rotation), and the reflection-allowed distance still exceeded
the replicate's recorded shape-null q95 in 99 of 100. The recomputed proper distances matched the persisted
observed shapes to four decimals. The latent shape response is genuine configuration deformation.

### What that leaves for readiness item 2

With reflection eliminated, the remaining account is the one the Phase-4 report labeled "a property of a
rank-3 stage-supervised space": an orthogonal (rotation-only) difference between two feature-space programs
does not survive a rank-3 projection *as* a rotation — only rotations within the retained subspace do. At
e = 1 the two groups' programs are nearly disjoint, so their relative rotation lies mostly outside any
shared 3-D subspace, and its image under the projection returns as deformation, which the shape statistic
(nearly pivotal, per the pivotality tables) then detects at rate ~0.99. Under this account the response is
an expected property of measuring an out-of-subspace rotation in a rank-limited space — predeclarable for
Phase 5 as a projection artifact, and testable: it should shrink as retained rank grows toward the joint
program rank.

---

## F4 — The angle test and its attribution decompose different orientations

`_estimate_orientation` returns PC1 of the k-stage configuration signed by net displacement: the tested
`angle` is a **principal-axis divergence**. `stats/attribution.py` explains **per-adjacent-transition
normalized contrasts** (`d/‖d‖` per transition). These coincide only for k = 2. For a bent trajectory PC1
is not any transition direction — and under the random-simplex baseline (F1) the divergence between the two
notions is maximal: transition directions are well-defined per step while PC1 is nearly degenerate.

Two documentation consequences: the orientation-invariance spec's purpose line ("a difference in direction
of progression") and its "known angle is recovered" scenario hold only for straight trajectories; and the
roadmap's Phase-3 plan to extend attribution with "signed principal orientations" was not implemented.
Nothing here is a wrong computation — but a significant 4-stage `angle` is currently interpreted through a
decomposition of a related, different quantity.

---

## F5 — Synthesis: why the assumed feature→latent mapping kept failing

The study's load-bearing assumption — *surgery mode ⇒ its target latent statistic* — fails at five distinct
layers:

1. **Requested → realized surgery.** The clamp censors the orientation and translation axes (F2). The
   effect knob is not the effect delivered.
2. **Methylation surgery → multi-omic geometry.** The CpG→gene→protein OR-derivation converts a pure
   methylation relocation into a mixed multi-omic program change (documented in Phase 2 — this layer is
   understood).
3. **Feature geometry → standardized joint geometry.** Pooled per-feature scaling, and magnitude scaling
   only one of three blocks, mix size into direction (documented — the 39.9° population angle of the
   magnitude construction).
4. **Joint geometry → latent geometry.** The stage-supervised space saturates at rank k−1 and retains ~6%
   of the group orientation contrast; rotation-only differences that lie outside the retained subspace
   re-enter as deformation (F3's mechanism for orientation→shape). Angles between projections are not
   angles between originals.
5. **Latent geometry → statistic.** PC1 on a near-isotropic configuration is an ill-determined estimand;
   observed statistic and permutation null co-vary through the shared anisotropy (F1), which *is* the
   measured 0.811 tracking.

`delta` survives all five because path length is a sum of norms — insensitive to axis identity, projection
rotation, and configuration anisotropy. `shape` survives because the constructed bend is large and its null
is nearly pivotal. `angle` is structurally exposed at every layer, which is why every fix so far moved it
by points rather than tens of points.

**The deepest of the five is the baseline.** "Intentionally not forced continuous" was a defensible
generality choice, but it makes the orientation estimand itself nearly meaningless at the population level:
a trajectory with independent stage programs has no direction to differ in. It is also the less
biologically representative case for staged progression, where successive stages share and accumulate
programs and the configuration has a dominant trend — precisely the regime where PC1 is well-conditioned
and angle nulls are tight. A **baseline-continuity axis** (from fully independent to strongly trending
indicator draws) would turn the orientation question from "which replicates got lucky draws" into a
designed, monotone operating characteristic — while keeping the current isotropic baseline as the
stress-test endpoint.

---

## D — Documentation and spec defects

| ID | Location | Issue |
|---|---|---|
| D1 | `CLAUDE.md:73` | Says shape uses "iterative Procrustes GPA" — superseded 2026-08-04 by strict pairwise OPA; contradicts the audit report and the code. |
| D2 | `README.md:52–57` | `motco simulate` example passes `--prop-affected-features`, which the CLI does not accept, and the comment says the command "requires R + InterSIM" while the parser help says "numpy generator, no R" (`cli.py:321`). |
| D3 | `README.md:304–315` | Python example imports `generate_semisynthetic_trajectory_from_intersim` and passes `prop_affected_features` — neither exists in `src/` anymore; the example raises on execution. |
| D4 | `docs/roadmap.md:28` | The resolved shape-invariance item is filed under "Not yet established". |
| D5 | `openspec/specs/semisynthetic-trajectory-generator/spec.md:5` | Purpose line still says datasets are generated "from InterSIM outputs"; the requirements below it correctly say no R at runtime. |
| D6 | `grid.py:328`, `showcase.py:71–74` | Matched-seed docstrings claim family cells generate "the *same* dataset" / that "the only thing that changes is group B's transform". Only the baseline indicators and group assignment match: the transforms consume different RNG amounts, so all sampled values differ across modes (the sole common-random-numbers pair is none↔magnitude, whose transforms draw nothing). The pivotality config note "the comparison across modes is on the same generator draws" overstates likewise. |
| D7 | `CLAUDE.md` / `evaluation.py` vs `docs/roadmap.md` | CLAUDE.md and the evaluation docstring call SNF a "production latent-space method" while rung 2 found SNF leaks independently and the roadmap defers SNF pending graph-native statistics (Phase 7). One message should win. |
| D8 | `evaluation.py:368`, `design.py` | Stage levels sort as strings; ordering silently breaks at n_stages ≥ 10 ("10" < "2"). Latent hazard only, but path length and shape depend on stage order. |
| D9 | `README.md:270`, shape spec | Shape described as removing "orientation differences" with reflections retained — needs the F3 rank qualifier to be true as stated. **Resolved 2026-09-03** by `unify-shape-reflection-policy`: rather than adding the qualifier, the policy was unified — reflections are aligned away at every ambient dimension, and README, spec, and docstrings now say so. |

---

## S — Audited and found sound

- **S1 — RRPP.** The reduced model (additive group+stage) is the right null for the interaction; permuting
  whole residual rows is exchangeable under the generator's homoskedastic cells; keeping the latent space
  fixed under permutation is defensible because integration is group-blind (stage-supervised only) — and
  empirically Type I holds. The add-one right-tail p-value is applied consistently.
- **S2 — Design coding.** `get_model_matrix` / `build_ls_means` / contrast indexing are mutually consistent
  (drop-first dummies, group-major LS-mean rows).
- **S3 — Procrustes solver.** The rotation solution and determinant handling are algebraically correct; the
  independent recomputation reproduced all 16 persisted shape statistics exactly.
- **S4 — Orientation sign anchor.** Net-displacement anchoring is the right convention; its one degeneracy
  (closed trajectories) is documented. The remaining instability is the axis, not the sign (F1).
- **S5 — M-value contract.** Applied uniformly across preprocessing, concat, SNF, and diagnostics; the
  rung-ladder conclusion is honored in the production path.
- **S6 — Study machinery.** Parameter signatures (including `n_jobs`), resumable sharding, seed derivation,
  and record persistence reproduce byte-identically — verified by regenerating replicates from persisted
  seeds and matching statistics to 1e-8.
- **S7 — The pivotality report's logic.** The studentization no-op argument, the fixed-threshold
  catastrophe, and the "tracking is load-bearing" reading are all correct. F1 extends it with the missing
  covariate; it contradicts nothing.
- **S8 — Magnitude path.** The one construction/estimator pair with no unexplained behavior: the estimator
  is more specific than the (deliberately mixed) construction, as reported.

---

## Remediation plan

Proposed as OpenSpec changes, in priority order. P1–P2 are small, additive, and unblock the Phase-5
design-point work; P3 is contract hygiene the paper-grade study's cross-checkpoint claims depend on; P4 is
the substantive scientific decision; P5–P6 are cleanups that can land any time.

### P1 — `record-latent-config-spectrum` (F1 → readiness items 1, 3, 4)

Persist the eigenspectrum of the centered stage-mean configuration — pooled and per group — in every
evaluation record, beside `null_summary`.

- Scope: `evaluation.py` (compute after `estimate_difference`, from `obs_vect` already in hand),
  `grid.py` record schema (version bump so resume signatures cannot mix contracts), study report tables
  (power stratified by eigengap for orientation cells).
- Also record the per-permutation eigengap summary (quantiles over the RRPP draws) — one SVD of a k×d
  matrix per permutation, negligible next to `estimate_difference`.
- Acceptance: the pivotality analysis re-run over new records shows log(null q95) regressing on the
  recorded eigengap with the association found here; the Phase-5 report contract gains the stratified
  power table.
- Non-goal: no change to any statistic or test.

### P2 — `fix-effect-axis-censoring` (F2 → readiness item 4 precondition)

Make the requested-vs-realized surgery relationship explicit and non-degenerate.

- Decide the policy per mode: fail loudly when the clamp would bind (preferred for study configs), or
  renormalize the axis to realized fractions, or raise headroom (e.g. lower p_dmp for the orientation
  cell, or relocate within the active set). The current silent `min()` is the one indefensible option.
- Surface `transform.orientation_relocated` / `translation_set_size` in `summarize_rejection_rates` and the
  study report; annotate or dedupe cells whose realized construction is identical in distribution.
- Acceptance: enumerating the Phase-4 grid under the new policy either errors on the censored cells or
  reports distinct realized effects per cell; no two power cells share > 5% identical replicates.
- Note: changes generator semantics → parameter signatures will (correctly) refuse to resume old shards.

### P3 — `unify-shape-reflection-policy` (F3 → readiness item 2, Phase-5 report contract)

**Landed 2026-09-03**, taking the preferred option below: the determinant correction is gone,
`SHAPE_STATISTIC_VERSION = 2` guards resume compatibility, and the follow-up probe ships as
`scripts/latent_rank_probe.py` (results under `results/latent-rank-probe-2026-09-03/`).

Pick one reflection policy that holds at every checkpoint.

- Preferred: **allow reflections everywhere** (drop the determinant correction). Measured cost on real
  replicates: zero (100/100 identical distances); benefit: one statistic across all checkpoints, matching
  the standard identification of reflections in shape spaces whose dimension exceeds configuration rank.
- Alternative: keep the proper-rotation estimator but document the rank threshold in the spec and annotate
  every cross-checkpoint comparison with its regime.
- Either way: add mirror-pair tests at ambient dimension > rank (the current suite pins only 2-D), and
  update `trajectory-shape-invariance` spec + README + docstrings.
- Include the item-2 follow-up probe: re-measure latent shape at forced larger `n_components`; if the
  orientation→shape response decays toward the population value, predeclare it as a projection artifact.

### P4 — `add-baseline-continuity-axis` (F5 → the 0.80-floor question, readiness item 4)

Add a generator axis controlling stage-program continuity, from independent (current behavior, the
stress-test endpoint) to strongly trending (e.g. Markov persistence of indicators across stages, or
cumulative programs).

- The design decision to make explicitly: orientation power is then reported *as a function of baseline
  continuity*, and the honest scientific claim becomes "MOTCO resolves orientation differences when the
  trajectory has a direction to differ in" — with the eigengap (P1) as the observable that carries that
  statement to real data.
- Realized-geometry diagnostics already in place will show how continuity moves the eigengap distribution.
- This is the change most likely to move the orientation operating characteristics by tens of points;
  n-scaling alone attacks the noise term while the gap's lower tail persists.

### P5 — `align-orientation-attribution` (F4)

Either add a principal-orientation component to attribution (contrast of the PC1 axes reconstructed
through the fitted loadings, beside the existing per-transition contrasts), or narrow the documented
estimand of `angle` to "principal-axis divergence" and state explicitly that attribution explains
transitions. Small; matters for how a significant orientation result is interpreted in Phase 6.

### P6 — `docs-sync-audit` (D1–D9)

One documentation batch: CLAUDE.md GPA line (D1), README simulate example and Python examples (D2, D3),
roadmap bullet placement (D4), semisynthetic spec purpose line (D5), matched-seed docstring wording in
`grid.py`/`showcase.py` and the pivotality config note (D6), one consistent SNF status message (D7), the
stage-label natural-sort guard plus test (D8 — the one code change in the batch), and the reflection
qualifier in shape descriptions (D9, folds into P3 if that lands first).

### Sequencing and cautions

- P1 + P2 before any new power runs; P3 before Phase-5 report drafting; P4 decides the Phase-5 design
  point; P5/P6 anytime.
- **Caution on latent-space sizing (item 3):** any group-aware supervision or sizing puts group information
  into the measurement space and voids the fixed-latent-space RRPP conditioning that currently makes the
  permutation test defensible (S1). That road requires integrating inside the permutation loop or a proof
  of invariance — do not take it casually.
- The zero-effect anchor, `n_jobs`-in-signature, and matched-seed machinery need no changes.

---

## Reproduction

All three experiments read committed study records and re-derive datasets from persisted seeds; none
re-runs RRPP.

**1. Reflection rank-dependence** (pure estimator check):

```python
import numpy as np
from motco.stats.trajectory import _estimate_shape

config = np.array([[0.0, 0.0], [1.0, 0.2], [2.0, -0.1], [3.0, 0.8]])
mirror = config * np.array([1.0, -1.0])
for d in (2, 3, 10, 660):
    A = np.hstack([config, np.zeros((4, d - 2))])
    B = np.hstack([mirror, np.zeros((4, d - 2))])
    V = np.vstack([A, B])
    print(d, _estimate_shape(V, [[0, 1, 2, 3], [4, 5, 6, 7]])[0, 1])
# 2 → 0.431;  3, 10, 660 → 0.000
```

**2. Eigengap vs recorded angle-null width** (regeneration from persisted seeds):

```python
import json
import numpy as np
from motco.simulations.semisynthetic import SemiSyntheticTrajectoryParams, generate_semisynthetic_trajectory
from motco.simulations.evaluation import (
    SimulationEvaluationParams, integrate_semisynthetic_dataset, build_simulation_trajectory_design)
from motco.stats.trajectory import estimate_betas

PLS = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 5, "max_components": 20, "random_state": 1203}
recs = [json.loads(l) for l in open("results/angle-pivotality-2026-09-01/merged.jsonl")]
orient = [r for r in recs if (r.get("cell_metadata") or {}).get("trajectory_mode") == "orientation"
          and r.get("phase") == "power_primary"]

def relgap(cfg):
    X = cfg - cfg.mean(0, keepdims=True)
    e = np.linalg.svd(X, compute_uv=False) ** 2
    return (e[0] - e[1]) / e.sum()

for r in orient:  # or a subset
    ds = generate_semisynthetic_trajectory(SemiSyntheticTrajectoryParams(
        seed=int(r["generator_seed"]), trajectory_mode="orientation",
        n_samples=300, n_stages=4, group_effect_size=1.0, p_dmp=0.2))
    latent = integrate_semisynthetic_dataset(ds, SimulationEvaluationParams(
        integration_method="pls", integration_params=PLS, permutations=0))
    design = build_simulation_trajectory_design(ds.metadata)
    Y = latent.matrix.to_numpy(float)
    stages = ds.metadata["stage"].to_numpy()
    pooled = np.vstack([Y[stages == s].mean(0) for s in sorted(set(stages))])
    # sanity: estimate_difference here reproduces r["pair_statistics"] exactly
    print(r["generator_seed"], relgap(pooled), r["null_summary"]["angle"]["q95"])
```

The reflection probe (F3) is the same regeneration with the group LS-mean configurations
(`build_ls_means @ betas`), `_center_scale_unit` on each, and the Procrustes residual computed with and
without the determinant correction.

**3. Effect-axis censoring** (record scan):

```python
import json
from collections import defaultdict
import numpy as np

rows = defaultdict(list)
for l in open("results/phase4-2026-08-27/merged.jsonl"):
    r = json.loads(l)
    cm = r.get("cell_metadata") or {}
    tr = (r.get("truth_metadata") or {}).get("transform") or {}
    if cm.get("trajectory_mode") == "orientation":
        rows[cm["effect_size"]].append(tr.get("orientation_relocated"))
for e in sorted(rows):
    print(e, np.mean([v for v in rows[e] if v is not None]))
# 0.25 → 54.4, 0.5 → 108.9, 0.75 → 146.9, 1.0 → 149.3 (nominal 54/108/163/217)
```

The 80/100 duplication figure comes from joining orientation `power_primary` records on
`replicate_index` across effects 0.75 and 1.00 and counting pairs with equal `generator_seed` and equal
`pair_statistics["angle"]` (tolerance 1e-9).
