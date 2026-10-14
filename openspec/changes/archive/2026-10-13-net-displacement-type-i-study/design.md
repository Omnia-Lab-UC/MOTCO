# Design

## Context

See proposal.md (Why). The relevant current state, verified against the code and by enumerating the
committed Phase 6 profile with extra axes:

- **The axis mechanism is already generic.** `grid._apply_axis_value` replaces any dataclass field for
  `generator.<field>` / `evaluation.<field>`, and `StudyConfig.axis_baseline_value` reads any field. Adding
  `evaluation.orientation_estimator: [null, "net_displacement"]` to the Phase 6 design grid enumerates 28
  cells today; adding `evaluation.integration_method: ["pls", "concat"]` crosses to 54. The estimator field
  is in `_SIGNATURE_OMIT_WHEN_NONE`, so a `null` column is byte-identical to no axis.
- **The design grid only crosses.** `DesignGrid.axes` → `itertools.product` in `StudyConfig.design_points()`.
  `_design_point_cells` and `_design_anchor_cell` consume `design_points()` and nothing else; the readers
  (`record_design_point`, `resolve_operating_by_design_point`, the decision rules) consume the per-record
  `cell_metadata.design_point` dict and never the config's grid shape.
- **A null-only study is expressible without a new knob.** `trajectory_modes: ["translation"]`,
  `effect_sizes: [0.0, 1.0]` with the shared anchor gives, per column, one anchor plus one translation cell
  (6 cells on the Phase 6 profile with one design axis). The engine also always emits the two
  `type_i_baseline` cells (`none` at 0 and `translation` at the top effect). They are measured in the
  baseline column but carry their own `control:<cell_id>` seed families, so they are independent datasets,
  not copies of the primary anchor and translation cells (corrected during apply).
- **Readers' PLS assumptions.** `resolve_operating_by_design_point` reads `selected_lv` and
  `component_selection` from `integration_metadata`, tolerating absence (`None`); `_component_selection_of`
  returns `None` when no record carries a marker. `assert_production_component_selection` admits a forced
  record when `record_forced_rank(record)` (read from the record's design point) equals its recorded
  `selected_lv`; this is grid-form agnostic. The `pca` method records `pca_latent` checkpoints; `concat`
  records no latent checkpoint. Whether the report figures (`render_design_point_power`,
  `render_rank_ladder`) survive a column with no `selected_lv` is the one thing not verified yet.
- **Decision-rule precedent.** `evaluate_design_point_decision` and `evaluate_rank_decision` are pure
  functions over records, write `<name>.json` + `.csv` through `scripts/motco_study.py report`, and return a
  `no_design_grid` verdict when there is nothing to evaluate.
- **Seeds.** Generator seeds derive from `(matched_seeds.primary_family, replicate_index)`, not from the
  cell, so two configs sharing a family and base seed generate identical datasets per replicate index. The
  RRPP seed is `evaluation.seed` (`run_simulation_replicate` passes it unchanged when set; every committed
  profile sets `3`), so it does *not* derive per cell: two cells with the same evaluation parameters on the
  same dataset draw the same permutations (corrected during apply).

## Goals / Non-Goals

**Goals:**
- One config per design expresses the five columns; everything else (sharding, resume, merge, identity pins,
  design-point reporting) is the existing engine unchanged.
- Every existing config, record and report output stays byte-identical. The identity fixture gains two
  entries and loses none.
- The Type I result and the follow-up selection are machine-readable and reproducible from committed files.

**Non-Goals:**
- No engine knob to suppress the `type_i_baseline` cells (see Risks). Not worth a new enumeration path
  for 1,000 units per design.
- No power cells, no attribution, no gate, no p-padding axis, no `pls_within_group` column, no per-permutation
  refit. No change to any production default or to the Phase 6 interpretability statement.
- No new SLURM template: `scripts/motco_study_array.sbatch` already runs any study config.

## Decisions

### D1. Explicit design points are a second form of `design_grid`, not a new top-level block

`design_grid` accepts either `{"axes": {...}}` or `{"points": [{axis: value, ...}, ...]}`.
`DesignGrid` gains a `points: tuple[Mapping[str, Any], ...]` field; `axes` is derived from the points
(ordered by first appearance in the first point, values in order of first appearance) so that every
existing consumer of `design_grid.axes` — `_validate_design_grid`'s baseline check, `prefer` validation,
`rank_decision` axis validation, `baseline_design_point`, the README — keeps working unchanged.
`design_points()` returns the listed points in listed order when `points` is set and the product otherwise.
Validation (load time, before any compute): both forms declared → error; a point whose key set differs from
the first point's → error naming the point and axis; two points equal under `_same_axis_value` → error; no
point equal to the baseline → error naming the baseline point; every axis passes `_validate_axis_namespace`
and `_normalize_axis_values` (so `null` for a nested key is normalized the same way as in the crossed form).

- *Alternative: allow `axes` plus an `exclude` list of combinations.* Rejected: the study needs 5 of 24
  combinations, so the exclusion list would be longer than the point list and harder to read.
- *Alternative: one config per measurement column (baseline = that space × PC1, one estimator axis).*
  Rejected: eight configs, each re-running its own PC1 baseline column and two `type_i_baseline` cells —
  roughly double the compute — and the pairing across configs would be implicit rather than enumerated.
- **Dump/reload.** `points` serializes as the list given; `axes` is not written when `points` is set, so the
  committed config stays the explicit form.

### D2. The estimator and method axes need tests and reader hardening, not new mechanism

Enumeration already works (Context). The work is: (a) tests pinning the three spec scenarios (estimator
applied per point, baseline byte-identical, forced column admitted through the explicit form); (b) a smoke
run through `scripts/motco_study.py report` on a tiny explicit-points study with `concat` and `pca`
columns, fixing whatever figure or table code assumes `selected_lv` is present. The expected fix sites are
`render_design_point_power` (it annotates the median eigengap, which is available for every method) and
`render_rank_ladder` (it must skip rows whose `component_selection` is `None` rather than treat them as
forced). The ladder figure will be rendered for these profiles because `forced_components` is a declared
axis; that is accepted — it is a correct figure of a one-rung ladder — rather than special-cased.

### D3. The anchor decision is a config-declared rule evaluated by the report, like the other two rules

`acceptance.anchor_decision: {alpha, se_multiplier, statistics, name?}` → `AnchorDecisionRule`.
`evaluate_anchor_decision(records, rule, *, alpha)` groups completed records by `design_point` (baseline
column included via its stamped coordinates), takes the anchor as `trajectory_mode == "none"` at effect 0,
computes rate/SE/R per statistic, bound = α + k·√(α(1 − α)/R) with the *available* R (so a failed replicate
widens the bound slightly rather than silently tightening the test), and classifies. Other cells at the
point with `negative_control` metadata or mode `translation` are reported descriptively under
`other_controls`. Output: `anchor_decision.json` (rule version, bound formula, per-point entries, per-point
verdict) and `anchor_decision.csv`. `scripts/motco_study.py report` writes it when the rule is declared. The
rule never touches `evaluate_phase4_gate` or `evaluate_targets`.

- *Alternative: compute the bound in the follow-up script from `design_point_operating.csv`.* Rejected: the
  rule would then live outside the config, unlike every other predeclared rule, and the report would not
  carry the verdict.
- *Why R per (point, statistic) rather than `n_replicates`.* The probe used the realized count; matching it
  keeps the bound comparable across the two studies.

### D4. The profiles are null-only through the existing negative-control machinery

`trajectory_modes: ["translation"]`, `effect_sizes: [0.0, 1.0]`, `shared_zero_effect_anchor: true`. Per
column: the anchor (`none`, 0.0) and `translation` at 1.0. Translation is the second null for `angle` and
`shape` and the headroom check applies to it at every column (it is evaluation-only variation, so headroom is
identical across columns). `acceptance.specificity` declares the three translation targets and
`acceptance.type_i` the α = 0.05 / 2-SE target, exactly as the Phase 6 profile does; `acceptance.power` is
empty and `gate.enabled` is false. `report_contract` is omitted (no attribution, so no driver report), which
means `--n-jobs` is only warned about; the profile notes say not to pass it.

- **Phase 6 seeds.** `base_seed: 800`, `primary_family: "phase6-small-n"`, generator block copied from
  `phase6_small_n_study.json`. Generator seeds are keyed on `(base_seed, family, replicate_index)`
  (`derive_replicate_seed`), so the baseline anchor and translation cells generate the same datasets as the
  Phase 6 study's at every replicate index. They are *not* the same cells: `make_simulation_cell` hashes the
  metadata into the cell id, and the stamped `design_point` coordinates differ between the two configs, so
  the cell ids and signatures differ. The evaluation parameters, and with them the RRPP seed (`evaluation.seed
  = 3`), are identical, so the baseline column's records must reproduce the parent study's p-values exactly,
  replicate by replicate; the report presents that exact reproduction as the within-family consistency check
  (as the magnitude re-measurement did for the Phase 5 anchor).
- **Phase 5 seeds.** `base_seed: 600`, `primary_family: "phase5-primary"`, `magnitude_kind: "joint"` (the
  probe's Phase 5 design). At zero effect the generator ignores `magnitude_kind`, so the anchor datasets equal
  the paper-grade anchor's; the cell id differs because the signature includes `magnitude_kind`, but the
  evaluation parameters (RRPP seed included) are identical, so the anchor and translation records must
  reproduce the Phase 5 paper-grade records exactly. Task 5.3 checks the first three replicates' datasets by
  hash against a fresh enumeration of `phase5_power_study.json`.
- **Ranks.** `pls_forced` and `pca` use rank 2(K − 1): 4 at Phase 6, 6 at Phase 5, as in the probe.
- **Layers.** Phase 6 sets `integration_params.layers: ["methylation", "expression"]` in the baseline
  evaluation; the design points inherit it because they only override the four named axes.

### D5. Follow-up selection is a separate script over the two decision files

`scripts/net_displacement_follow_up.py --phase6 <dir>/report/anchor_decision.json --phase5 ... --out
follow_up_decision.json`. Pure function over the two parsed decisions; the lead column is identified by its
coordinates (`integration_method: pls`, `forced_components: null`, `n_components: null`,
`orientation_estimator: net_displacement`). Writes verdict, localization per design, rule version, input
paths and their sha256. Unit-tested on synthetic decision dicts for the two spec scenarios.

- *Alternative: a cross-config rule inside the engine.* Rejected: the engine's reports are per config, and a
  two-config rule has no home there.

### D6. Run and commit, same contract as every committed study

Timing job on `ing` (two anchor units and one translation unit per column per design), budget written to
`PROVENANCE.txt` before launch, BLAS pinned to one thread, `STUDY_N_JOBS` unset, 50 shards per design,
memory-pressured nodes avoided. `results/net-displacement-type-i-<date>/{phase6,phase5}/` each with
`merged.jsonl`, `report/`, `PROVENANCE.txt`; `follow_up_decision.json` at the top level. The findings report
follows the Phase 5 template's provenance and reproduction sections and adds the anchor table, the
exact-reproduction check against the parent studies' records, and the verdict.

## Risks / Trade-offs

- **[The `type_i_baseline` cells cost ~1,000 units per design, mostly `pls_cv`]** → Accepted. They are a
  second, independently seeded measurement of the baseline column's anchor and translation rates, reported
  descriptively under the anchor decision's `other_controls`.
- **[The baseline column re-computes the parent studies' records (~1,000 `pls_cv` units per design)]** →
  Accepted: it is the column every other one is read against, the engine requires it, and its exact
  reproduction of the parent records is a free end-to-end determinism check on the cluster.
- **[Phase 5 compute is the dominant cost (≈ 100–130 core-hours) for a null-only answer]** → The alternative
  is to run Phase 6 first and Phase 5 only if Phase 6 holds. Rejected for the primary run: the roadmap and the
  probe predeclare both designs, and a Phase 6-only pass would leave the Phase 5 excess (0.08) unsettled. The
  timing job is the go/no-go; if the projection exceeds 150 core-hours the Phase 5 profile drops `pca` and
  `concat` to 10 cells and the report says so.
- **[A report figure or reader assumes PLS metadata and fails on `concat` / `pca` columns]** → Smoke run
  (D2) before any cluster time; the failure mode is a crash, not a wrong number.
- **[`render_rank_ladder` is produced for a study that is not a rank ladder]** → Accepted; it is correct for
  the one forced rung. The findings report does not use it.
- **[The bound at R = 500 (0.0695) is tighter than at R = 100 (0.0936); a true 0.08 excess will be called]**
  → That is the intent: the question is whether the 0.08–0.10 excess is real, and at R = 500 the SE is
  0.0097, so 0.08 is three SE above α.
- **[Translation at e = 1.0 is pool-limited; its headroom is checked but the realized surgery differs by
  design]** → Already the case in the parent profiles; `realized_surgery.csv` is written as usual.
- **[The explicit form lets a config list a point that is also a baseline-adjacent duplicate under
  `_same_axis_value` normalization (e.g. `"pc1"` vs `null`)]** → Duplicate detection runs after
  normalization, so such a pair is rejected as a duplicate point.

## Migration Plan

No migration: the config schema gains an optional form and an optional acceptance block; every committed
config loads and enumerates unchanged, which `tests/test_committed_config_identity.py` enforces. Rollback is
reverting the change; results directories are additive.

