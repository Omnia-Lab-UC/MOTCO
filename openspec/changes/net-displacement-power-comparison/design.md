# Design

## Context

See proposal.md (Why). The current state that shapes the approach, verified against the code and the
committed Type I run:

- **Everything the comparison needs to enumerate and run already exists.** The explicit-points design grid,
  `evaluation.orientation_estimator` as a contracted design axis, matched-seed families with a shared
  zero-effect anchor, design-point reporting (`design_point_operating.csv`, `render_design_point_power`), and
  `bound: "upper"` on Type I / specificity targets all shipped in the last three changes. A config with two
  explicit points over `evaluation.integration_method` × `evaluation.orientation_estimator` and the four
  modes enumerates 36 cells today with no engine change.
- **The two `pls_cv` columns share everything but `angle`.** The estimator enters `estimate_difference` only
  through `_estimate_orientation`; `delta` and `shape` are computed from the same betas and LS-means, and
  the RRPP seed is `evaluation.seed` (constant 3), so for every (cell, replicate) the two columns' `delta`
  and `shape` p-values and observed statistics are identical by construction. The comparison is therefore
  entirely about the `angle` column of every mode: its power under orientation, its leak under magnitude
  and shape, and its size under the anchor and translation.
- **The existing decision rules do not fit.** `rank_decision` is the closest (gain / anchor / protected
  against a reference column) but it keys columns on an integer rank axis, has no mandatory-control
  criterion and only a `gain` target. `anchor_decision` judges only the anchor. The probe's advancement rule
  (`scripts/orientation_measurement_probe.py`) had the right five criteria but lives outside the engine and
  reads the probe's own tables. The helpers a new rule needs are already in `targets.py`: `RateEstimate`,
  `_rate_estimate`, `_pooled_se`, `_top_effect_records`, and `phase4.type_i_inflation_bound`.
- **Cost.** From the Type I run's `PROVENANCE.txt` (recorded per-column medians under contention): `pls_cv`
  18.4–18.8 s per unit at Phase 6 and 73.6–75.3 s at Phase 5; `concat` 74.9 s and 142.1 s. The Phase 5
  paper-grade run was about 150 core-hours; the Type I study 133.
- **Parent records.** The Phase 6 baseline column re-computes `results/phase6-small-n-2026-10-09/` (modes
  orientation/shape/translation at 0.25–1.0, anchor) and `results/phase6-small-n-magnitude-2026-10-09/`
  (magnitude at 0.25, 0.5, 1.0; its axis has no 0.75). The Phase 5 baseline column re-computes
  `results/phase5-2026-09-10/` for orientation/shape/translation and the anchor (those modes ignore
  `magnitude_kind`) and `results/phase5-magnitude-2026-10-08/` for magnitude at 0.25 and 1.0; magnitude at
  0.5 and 0.75 under `joint` is new at Phase 5.

## Goals / Non-Goals

**Goals:**
- A predeclared, config-declared, machine-evaluated rule for "net displacement wins", with the translation,
  magnitude and shape `angle` controls as hard criteria rather than descriptive rows.
- Two profiles that re-measure the parents' cells exactly in the baseline column, so the only new numbers are
  the candidate column's `angle` rates and the comparison is paired on identical datasets.
- Every existing config, record and report byte-identical; the identity fixture gains two entries.

**Non-Goals:**
- No engine support for per-mode effect axes or per-column cell subsets; the single `effect_sizes` list is
  accepted (D3).
- No `concat`, `pca`, forced-rank or `pls_within_group` column (proposal; D3).
- No production-default flip in this change. If the verdict is `adopt_net_displacement`, the flip is a
  separate change that re-pins signatures, as `default-joint-magnitude-construction` did.
- No change to attribution, the gate, the Phase 5 report template, or the Phase 6 interpretability statement
  before the findings report.

## Decisions

### D1. The rule is `acceptance.column_decision`, a general candidate-vs-baseline rule, not an estimator-specific one

`acceptance.column_decision: {candidate: {axis: value, ...}, target: {trajectory_mode, statistic},
criterion: "gain" | "no_loss", gain_se_multiplier?, loss_se_multiplier?, type_i_bound?: {alpha,
se_tolerance}, controls: [{trajectory_mode, statistic}, ...], protected?: [...], name?}` →
`ColumnDecisionRule`. The candidate is matched against `config.design_points()` with `_same_axis_value`
after `_normalize_axis_values` (so `"pc1"` ≡ `null`); it must match exactly one point and that point must not
be the baseline. Partial coordinates are allowed (every named axis must be a declared design axis) so a rule
can name `{"evaluation.orientation_estimator": "net_displacement"}` on a wider grid.

`evaluate_column_decision(records, rule, *, alpha)` in `targets.py` groups completed design-point records by
`record_design_point` (baseline included via its stamped coordinates), finds the candidate group by the same
matching, measures each column with the rank rule's `_measure_rank_column` pattern extended to the control
pairs, and evaluates the four criteria from the spec. The anchor and control bounds use
`type_i_inflation_bound(alpha_b, k, R)` — the same function the gate, the upper-bound targets and the anchor
rule use. Other design columns (none in the committed profiles) get measures and no verdict. Output:
`column_decision.json` (rule payload, verdict, rationale, baseline/candidate/other columns with criteria) and
`column_decision.csv` (one row per column, criterion columns flattened like `RankDecision.to_frame`).
`scripts/motco_study.py report` writes both when the rule is declared, after the anchor decision.

- *Alternative: extend `rank_decision` with a generic axis and a controls list.* Rejected: it would change a
  rule whose output is committed (`phase5-latent-rank-2026-09-08`) and whose semantics (smallest qualifying
  rank, `keep_cv`) do not apply to a two-column comparison.
- *Alternative: a cross-config script only, reading `design_point_operating.csv`.* Rejected for the same
  reason the anchor rule lives in the config: every predeclared rule is in the config and its verdict is in
  the report.
- *Why `criterion` is per config.* Phase 6 is where the gain is claimed (0.29 → 0.67 in the probe) and the
  rule should require it. At Phase 5 PC1 already saturates at 0.84–0.88 and the probe's net-displacement
  rate was 0.89–0.93; requiring a 2-pooled-SE gain there (about 0.045) would fail the candidate for not
  improving a test that is already at its ceiling. `no_loss` is the honest Phase 5 question; the cross-design
  script combines the two.

### D2. Mandatory controls are the `angle` statistic under every non-orientation mode, at the top effect

`controls: [translation/angle, magnitude/angle, shape/angle]`, each held to the one-sided bound α + k·√(α(1 −
α)/R) with the pair's own available R. Translation is the control the roadmap predeclares (0.068 vs 0.0695 at
Phase 6 in the Type I run). Magnitude/`angle` and shape/`angle` are the estimator's two specificity exposures:
net displacement reads only the endpoint stages, so a shape surgery on the endpoints or a size change could
move it. Under PC1 the Phase 5 shape → `angle` rate is 0.89, so the *baseline* fails this control; the rule
judges the candidate only and the report shows both.

- *Why the bound and not "not above the baseline by more than 2 pooled SE" (the probe's form).* The probe's
  form would let a candidate inherit the baseline's leak (0.89 at Phase 5). The bound is absolute and is
  what the paper needs: specificity at level α.
- *Protected pairs.* `magnitude/delta` and `shape/shape` are declared for the rule's generality and because
  the spec requires the criterion, but between two `pls_cv` columns they are identical records and pass with
  loss 0. The findings report states this (D5) rather than presenting it as evidence.

### D3. Two columns, the parents' five-point axis, one `effect_sizes` list for all four modes

Points: `{method: pls, estimator: null}` (baseline) and `{method: pls, estimator: net_displacement}`. Modes
`magnitude`, `orientation`, `shape`, `translation`; effects `[0.0, 0.25, 0.5, 0.75, 1.0]`. Per design: 2
Type I baselines + 2 × (1 anchor + 16) = 36 cells, 18,000 units.

- **Axis.** The engine has one `effect_sizes` list. Phase 6's parents split magnitude (0.02–1.0) from the
  other modes (0.25–1.0) because magnitude/`delta` rises in 0.02–0.25. Here magnitude's `delta` is
  estimator-independent and already measured; what the comparison needs from the magnitude cells is the
  `angle` leak where the size change is large, which 0.25–1.0 covers (Phase 6 `delta` 0.956–1.000 there).
  So both profiles use the orientation/shape axis; Phase 6 magnitude at 0.75 is a new cell (the parent axis
  had none), stated in the profile notes.
- **Why five points at Phase 5 and not the pilots' four.** Every baseline cell then has a parent record
  (except magnitude 0.5/0.75), so the baseline column's reproduction check covers the whole grid and the
  candidate's `angle` curve is on the same abscissa as the paper-grade Phase 5 curve. The cost of 0.75 is
  about 70 core-hours at Phase 5; the go/no-go (Risks) drops it if the timing job projects over budget.
- **No `concat`.** Proposal; the ceiling numbers are cited from the probe and the Type I study.
- **Design axes.** Two axes (`evaluation.integration_method`, `evaluation.orientation_estimator`) rather than
  the Type I profiles' four. Cell ids differ from the Type I profiles' either way (the stamped point enters
  the id), so nothing is lost; the method axis is kept so the baseline point reads as `pls × pc1` in every
  table and so a `concat` point could be added later without changing the axis set.
- **Acceptance block.** `type_i: [{alpha: 0.05, se_tolerance: 2, bound: "upper"}]`; `specificity`: the
  parents' five pairs (translation × 3, magnitude × `angle`/`shape`), all `bound: "upper"`; no `power`
  target (the Phase 6 0.80 floors are known to be missed under PC1 and the column rule is the predeclared
  rule here); `gate.enabled: false`; `attribution.enabled: false` (the driver table is a property of the PLS
  fit, already in the parents); no `report_contract` (nothing for it to govern; the profile notes say not to
  pass `--n-jobs`). `surgery_censoring` default (`error`); every cell is under headroom at these designs (the
  parents enumerate the same cells).
- **Seeds.** Phase 6: `base_seed 800`, family `phase6-small-n`; Phase 5: `base_seed 600`, family
  `phase5-primary`, `magnitude_kind: joint`. `evaluation.seed = 3` is also the RRPP seed, so the baseline
  column reproduces the parents' p-values exactly where a parent cell exists (Context), and the candidate
  column draws the same permutation indices.

### D4. Cross-design selection is a separate script, mirroring `net_displacement_follow_up.py`

`scripts/net_displacement_power_follow_up.py --phase6 <dir>/report/column_decision.json --phase5 ... --out
comparison_decision.json`. Pure function over the two parsed decisions; verdict `adopt_net_displacement`
iff both are `candidate_qualifies`, else `keep_pc1` with per-design verdict and failing criteria/pairs;
rule version, `advisory: true`, inputs with sha256. Unit-tested on synthetic decision dicts for the two spec
scenarios. The script does not decide the production default; it names the follow-up change.

### D5. The findings report adds two consistency checks that the engine does not compute

Both are computed from `merged.jsonl` by a small script committed under `scripts/` (or a documented
one-liner in the report's reproduction section), as the Type I report did for parent reproduction:

1. **Shared statistics are identical.** For every (cell coordinates, replicate) pair across the two columns,
   `delta` and `shape` p-values and observed statistics are equal. Any difference is a bug in the estimator
   plumbing, not a finding.
2. **Baseline reproduces the parents.** Baseline-column records equal the parent runs' records in generator
   seed and p-values wherever a parent cell with the same mode and effect exists (Context lists which).

The report's main table is `angle` rate by mode × effect × estimator per design; the control table is the
rule's `column_decision.csv`; covariates come from `design_point_operating.csv` (eigengap, `angle` null q95
median/IQR, selected rank). A power figure is allowed and expected (unlike the Type I report).

### D6. Run and commit, same contract as every committed study

Timing job on `ing` first: per column per design, the anchor (replicates 0, 1) and the orientation e = 1.0
cell (replicate 0); 4 single-CPU jobs; nodes under memory pressure or NFS stalls excluded as in the Type I
run. Budget written to `PROVENANCE.txt` before launch; BLAS pinned to one thread; `STUDY_N_JOBS` unset; 50
shards per design (100 for Phase 5 if the queue allows); `--error-policy record`. Output under
`results/net-displacement-power-<date>/{phase6,phase5}/` with `merged.jsonl`, `report/`, `PROVENANCE.txt`,
and `comparison_decision.json` at the top level.

## Risks / Trade-offs

- **[Phase 5 is ≈ 370 core-hours, the largest single study so far]** → Timing job go/no-go: if the Phase 5
  projection with a 20 % margin exceeds 450 core-hours, drop e = 0.75 from the Phase 5 profile before
  launch (28 cells, ≈ 290 core-hours), re-pin its identity entry, and say so in the profile notes and the
  report. Phase 6 (≈ 95 core-hours) runs regardless and is the design where the claim lives.
- **[The baseline column re-computes ~9,000 parent units per design]** → Accepted, as in the Type I study:
  the engine requires the baseline column, the comparison must be paired, and exact reproduction is the
  end-to-end determinism check. Reusing parent records by cell-id translation would be a new engine
  mechanism with its own provenance questions; not worth it for one study.
- **[Net displacement may fail a mandatory control at the bound by one or two rejections]** → That is the
  rule working. The report shows the rate, the bound and the margin; the verdict is `keep_pc1` and the
  follow-up is a diagnosis, not a re-run with a looser bound.
- **[The Phase 5 `no_loss` criterion is weaker than `gain`; a reader could see it as moving the goalposts]**
  → It is predeclared here with its reason (PC1 saturates at Phase 5), the rule payload is echoed in the
  decision output, and the report presents the Phase 5 gain/loss with its pooled SE either way.
- **[Phase 6 magnitude at e = 0.75 and Phase 5 magnitude at 0.5/0.75 are new cells with no parent record]**
  → Their `delta` is reported as descriptive and the reproduction check is stated as "wherever a parent
  record exists"; nothing in the rule reads them except their `angle` leak.
- **[A candidate that qualifies at both designs still leaves the Phase 6 orientation power well under 0.80
  (probe: 0.67)]** → The verdict selects the default change; it does not restore the 0.80 claim. The revised
  interpretability statement must quote the measured power, and the report says what a non-rejection then
  argues against.
- **[`column_decision` with a partial candidate on a wider grid could match several points]** → Rejected at
  load time (spec scenario).

## Migration Plan

No migration: the config schema gains one optional acceptance key; every committed config loads, enumerates
and reports unchanged, which `tests/test_committed_config_identity.py` and the byte-identity report fixture
enforce. Rollback is reverting the change; results directories are additive.
