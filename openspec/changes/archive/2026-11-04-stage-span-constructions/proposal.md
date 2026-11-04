## Why

MOTCO's production measurement projects the features onto a stage-supervised PLS basis fitted to both groups pooled. It sees a group difference only through the part that lies inside that basis's span. At rank K − 1 in a balanced design, that span is the span of the pooled stage means. An equal-length rotation then cancels: for steps `a` and `b` with ‖a‖ = ‖b‖, the pooled direction `a + b` and the difference `a − b` are perpendicular, the diagonals of a rhombus. The simulator's `orientation` surgery is such a rotation, a relocation of stage-changing CpGs. The plasmode plant is sparse and random, so most of it falls outside the span as well. The repository has measured the consequence three times without naming its cause:

- the orientation probe's Phase 5 latent angle was 6.1° against a population 90° (`docs/reports/orientation-measurement-probe-2026-10-12.md` §4);
- the plasmode plant kept 2.2 % of its norm in the rank-5 span (`docs/reports/plasmode-2026-11-02.md`);
- `angle` power tracks the CV-selected rank rather than the effect (the design-guidance pilot and slice).

Our power curves therefore mix two questions:

- how much of the planted difference the measurement can see;
- how sensitive the test is to what it sees.

A scratchpad prototype separated them. It edited group B's latent stage means (size, a 30° rotation, a bend), pulled the edit back to the features through the PLS loadings, and re-ran `analyze_trajectories` from scratch. Each edit was recovered in its own statistic and in no other (p = 0.005 at 199 permutations; off-target p ≥ 0.25), on the stage-changing features only.

The plasmode report names the missing check: "a dense difference aligned with the stage signal was not tested and would be a new check". This change adds that construction to the simulator. It also measures, for every construction, how much of the plant the fitted basis sees.

## What Changes

- **Stage-span construction in the generator.** A new opt-in `construction` field: `surgery` is the default, and `stage_span` is new. Under `stage_span`, `magnitude`, `orientation`, `shape` and `translation` edit group B's exact population stage means inside the span of B's own stage configuration. The edit is made in pooled-standardized joint coordinates over declared `span_layers` and added to B's cell means before sampling.
  - `magnitude` scales the configuration about its centroid by 1 + e.
  - `orientation` rotates it rigidly by θ = 90°·e in the plane of its net displacement and its next principal direction.
  - `shape` moves the middle interior stage off its chord by e mean steps, then restores the path length and the centroid.
  - `translation` offsets every stage by e mean steps along the net displacement.

  The construction consumes no randomness, so a stage-span dataset equals the null dataset at the same seed plus the plant. Truth records the exact intended population geometry. Absent the field, every dataset, cell id, signature and record is byte-identical.
- **Span-capture diagnostic.** A new opt-in evaluation setting, `span_capture`, adds two things to the realized geometry of a PLS evaluation:
  - a `population_latent` checkpoint: the population means projected through the replicate's fitted basis;
  - a `span_capture` block: the shares of the population group × stage difference and of the pooled stage displacement that lie inside the fitted span, with the rank and the random-direction share k/p.

  The study report gains `span_capture.csv` when records carry it.
- **Consistency guard.** A stage-span dataset may only be measured on its span layers. The study loader and the evaluation refuse a mismatch before fitting.
- **Pilot.** Four post-floor profiles at 100 × 199 in the `design-guidance` seed family: {surgery, stage-span} × {magnitude, geometry}. Each has three design columns: n = 200 and n = 80 balanced, and the Phase 6 SEA-AD table `[[11, 10, 28], [9, 10, 12]]`, all at three stages with methylation + expression. The surgery siblings reproduce the 2026-10-28 pilot's records at the shared columns.
- **Retrofit check.** `scripts/latent_retrofit_check.py` promotes the prototype: an edit in a replicate's own fitted basis, pulled back and re-measured with the basis frozen and refitted. It is descriptive, and it compares the fitted-basis footprint with the population stage-span plant.
- **Advisory reading and report.** `scripts/stage_span_pilot_reading.py` writes `stage_span_pilot_reading.json` (`propose_paper_grade` | `revise_construction`), and a dated report records the run.
- **Visibility note.** A concept page, `docs/measurement-visibility.md`, states what the stage-PLS measurement can and cannot see and links the evidence above. The roadmap and `AGENTS.md` are updated.

**Non-goals.** These are follow-up changes:

- the paper-grade profiles and their predeclared decision, which is set from the pilot's reading;
- a latent plasmode on the SEA-AD export, which needs a split-sample basis because real data has no population means;
- any change to `motco analyze`, RRPP or the production measurement.

## Capabilities

### New Capabilities
- `stage-span-construction-validation`: the stage-span pilot profiles, the latent retrofit check, the advisory reading rule and the dated pilot report.

### Modified Capabilities
- `semisynthetic-trajectory-generator`: the trajectory modes gain the opt-in stage-span construction beside the default feature-set surgery.
- `realized-geometry-diagnostics`: an opt-in `population_latent` checkpoint and `span_capture` shares for PLS evaluations.
- `simulation-evaluation-harness`: the `span_capture` setting, and refusal of a stage-span dataset measured on blocks other than its span layers.
- `trajectory-power-study`: configurations accept `generator.construction`, `generator.span_layers` and `evaluation.span_capture`; the anchor of a stage-span profile is the surgery anchor; reports write `span_capture.csv`.
- `latent-space-metric-compatibility`: documents which group differences the stage-supervised PLS measurement can see.

## Impact

- **Code:**
  - `src/motco/simulations/semisynthetic.py`: params, the construction, population trajectories and truth;
  - `src/motco/simulations/generator.py`: an optional per-cell mean shift in `generate_omics` and `omic_population_means`;
  - `src/motco/simulations/diagnostics.py` and `evaluation.py`: the checkpoint, capture shares and guard;
  - `src/motco/simulations/grid.py`: `_SIGNATURE_OMIT_WHEN_NONE` gains `construction`, `span_layers` and `span_capture`;
  - `src/motco/simulations/study/` (`config`, `enumerate`, `report`).
- **Scripts:** `scripts/latent_retrofit_check.py` and `scripts/stage_span_pilot_reading.py`.
- **Configs:** four profiles under `examples/trajectory_power_study/`, registered in `POST_FLOOR_CONFIGS` and pinned by `tests/data/post_floor_config_identity.json`.
- **Compute:** the pilot is about 8,400 units on `ing`, an estimated 60 core-hours from the design-guidance pilot's 26 core-seconds per unit.
- **Compatibility:** every committed configuration, fixture and dated report reproduces byte-identically. The new fields are hash-omitted while unset. The production measurement and `motco analyze` are unchanged. The version moves to MOTCO 0.15.0.
- **Docs:** `docs/measurement-visibility.md` (in the MkDocs nav), `docs/roadmap.md` (new item 13 and the "Not yet established" entry), `AGENTS.md` architecture, and the dated report under `docs/reports/`, with aggregates under `results/`.
