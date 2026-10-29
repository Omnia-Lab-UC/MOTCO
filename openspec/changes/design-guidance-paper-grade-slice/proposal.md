# Proposal

## Why

The design-guidance pilot (`design-guidance-grid`, 2026-10-28) measured every column descriptively at 100 replicates × 199 permutations: SE ≤ 0.05 on every rate, an anchor bound of 0.0936, two anchor excesses it could neither clear nor flag, and seven translation `angle` rates sitting above the paper-grade bound that its SE cannot resolve. The paper's design-guidance statement needs the slice where the rates are neither saturated nor at the null, at the sizes the planned analyses have, at paper grade (500 × 999, SE ≤ 0.022, bound 0.0695), with the selected rank recorded per cell. The pilot report §7 defines that slice; this change runs it.

## What Changes

- **Two paper-grade sibling profiles** under `examples/trajectory_power_study/`, derived from the pilot's core profiles and sharing its seed family (`design-guidance`, `base_seed` 1028), so every design column evaluates the pilot's datasets at replicates 0–99 and 400 new ones:
  - design grid `n_samples` ∈ {80, 200, 500} × three stages × `layers` ∈ {methylation + expression, all three} — six columns, a zero-effect anchor per column;
  - magnitude at e ∈ {0.02, 0.05, 0.1}; orientation, shape and translation at e ∈ {0.25, 0.5, 1.0};
  - 500 replicates × 999 permutations; production measurement named explicitly (PLS with the CV rank floored at 2, `net_displacement`, `joint`);
  - predeclared rules only: `acceptance.anchor_decision` at α 0.05, k = 2 (bound 0.0695), and a one-sided `acceptance.type_i` rule for translation's `angle` at every column. No power targets, no gate: the slice is guidance, not a gate.
  - 26 + 62 cells, 44,000 units (about 380 core-hours on `ing`; a smoke run fixes the shard plan). Translation runs at every geometry effect because a profile has one effect list; the pilot's geometry profile did the same, and the extra two translation effects per column (6,000 units) are descriptive controls.
- **Column-scoped Type I entries.** An `acceptance.type_i` entry may name a negative-control `trajectory_mode` (`none` or `translation`) and a `statistic`; it is then evaluated at every design column, the baseline included, on that control at its largest effect, under the entry's bound form. Today `type_i` entries read only the Type I baseline phase, so a per-column translation rule cannot be predeclared. Entries without a mode are unchanged and every committed acceptance report regenerates byte-identically.
- **Figure script grade label.** `scripts/design_guidance_figure.py` gains `--grade` (default `pilot`, so the committed pilot figures are unchanged) for the figure titles and renders a single stage count as one panel row.
- **Anchor reproduction by design column.** `scripts/anchor_reproduction.py` gains `--by-column`: anchors are joined on (design column, replicate index), one summary line per column, so the sibling-profile check (identical p-values) and the pilot check (same generator seeds and observed statistics at replicates 0–99, p-values not comparable across permutation counts) are committed tools rather than session scripts. Without the flag the output is unchanged.
- **Second priority: the stage-effect arm.** Two sibling profiles over explicit design points on the joint δ ∈ {0.25, 0.5, 1.0} at n = 200, three stages, two blocks (14 + 32 cells, 23,000 units, about 230 core-hours), committed and pinned with the slice, submitted only on the user's decision after the slice's numbers are in. δ = 0.25 is their config baseline, so no column repeats the slice; the slice's n = 200 / two-block column is the δ = 2.0 reference.
- **Paper-grade run and dated findings report** (`results/design-guidance-slice-<date>/`, `docs/reports/design-guidance-slice-<date>.md`): rates ± SE per column, stage AUROC and selected rank per cell, every column's anchor verdict and translation `angle` verdict, the sibling and pilot reproduction checks, the design-guidance figures at paper grade, and the guidance statement the paper will carry. Roadmap item 7(a), `AGENTS.md` and both READMEs updated; version bumped to 0.12.0 for the engine additions.
- Out of scope: the ρ, group-ratio, stage-proportion, `p_dmp`, four-stage and n = 1000 columns (pilot report §7 gives the reasons); per-mode effect lists in the engine (noted as the alternative to running translation at three effects); group-specific dispersion; the case-study-configuration Type I check.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

- `trajectory-power-study`:
  - extends "Study evaluates results against pre-specified targets" with column-scoped Type I entries (negative-control mode and statistic, evaluated at every design column);
  - extends "Design-guidance figure" with the grade label and the single-stage-count layout;
  - adds the design-guidance paper-grade slice profiles, the stage-effect arm profiles, anchor reproduction by design column, and the dated paper-grade findings report.

## Impact

- `src/motco/simulations/study/config.py`, `targets.py`: optional `trajectory_mode` / `statistic` on `TypeIControlTarget`, parsed and validated; column-scoped evaluation in `evaluate_targets`. Byte-identical for entries without a mode (acceptance report CSV/JSON, config echo).
- `scripts/design_guidance_figure.py` (`--grade`), `scripts/anchor_reproduction.py` (`--by-column`).
- `examples/trajectory_power_study/`: four new post-floor configs (`design_guidance_core_{magnitude,geometry}_study.json`, `design_guidance_stage_effect_{magnitude,geometry}_study.json`), appended to `POST_FLOOR_CONFIGS` and `tests/data/post_floor_config_identity.json`; never to the 19-config historical fixture.
- `tests/`: column-scoped target tests (`test_study_targets.py`, `test_study_config.py`), profile structure tests, figure and reproduction-script tests, fixture regeneration.
- Compute: about 380 core-hours on `ing` for the slice (plan for 450 with the pilot's 15–20 % overrun), about 230 more if the stage-effect arm is run. Submitted from a clean worktree of the pushed branch after the user approves the smoke-based plan.
- Docs: `docs/roadmap.md` (item 7(a) → done; Phase 5 sensitivity deliverable at paper grade), `AGENTS.md` study paragraph, `src/motco/simulations/study/README.md` (acceptance fields, figure and reproduction scripts), `examples/trajectory_power_study/README.md` (new section), the dated report, `results/design-guidance-slice-<date>/` (`report/`, `PROVENANCE.txt`, `figures/`; JSONL stays on the cluster).
