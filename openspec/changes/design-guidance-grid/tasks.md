# Tasks

## 1. Baseline on the unmodified tree

- [x] 1.1 Run the fast test suite on the unmodified tree and record that the pin tests pass with 19 historical configs. Verify with `MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`.
- [x] 1.2 Capture the pre-change report outputs:
  - build the small design-grid study used by `tests/test_study_end_to_end.py` (or an equivalent two-axis fixture study) and report it;
  - record the sha256 of every file the report writes in a scratchpad file, for the byte-identity check in 3.2.

  Verify that the hash list exists and that `design_point_operating.csv` is in it.

## 2. Design-guidance profiles and the pin policy

- [x] 2.1 Write the four profiles under `examples/trajectory_power_study/`, following design D1–D4 and D8:
  - `design_guidance_core_magnitude_pilot.json`
  - `design_guidance_core_geometry_pilot.json`
  - `design_guidance_arms_magnitude_pilot.json`
  - `design_guidance_arms_geometry_pilot.json`

  Shared settings:
  - **Baseline:** n 200, 3 stages, `group_ratio` 0.5, `p_dmp` 0.1, ρ 0, every per-omic δ 2.0, `magnitude_kind: "joint"`. Evaluation: `pls`, `orientation_estimator: "net_displacement"`, `min_components: 2`, `layers: ["methylation", "expression"]`, CV knobs as Phase 6, 199 permutations, `n_jobs` 1, `seed` 3.
  - **Run size:** 100 replicates, `base_seed` 1028, family `design-guidance`, attribution disabled.
  - **Acceptance:** a `type_i` entry and `anchor_decision` (α 0.05, k 2, all three statistics).
  - **Metadata:** a block naming the design-guidance purpose.

  Verify that each profile loads and enumerates to 98 / 162 / 68 / 112 cells. If a count differs, record the cause.
- [x] 2.2 Write the fixture capture (`tests/post_floor_fixture_capture.py`, modelled on `tests/phase6_fixture_capture.py`) and capture `tests/data/post_floor_config_identity.json`. It holds `capture_revision`, `note`, and `configs` mapping each cell id to its signature and first three matched seeds. Verify that two captures are identical.
- [x] 2.3 Rework the coverage and floor checks in `tests/test_rank_floor_default_pins.py` and `tests/test_orientation_default_pins.py` (design D7):
  - the historical and post-floor fixtures are disjoint, and their union equals `examples/trajectory_power_study/*.json`;
  - the historical set stays at exactly 19;
  - historical configs keep `min_components == 1`;
  - post-floor configs name a non-null floor and a non-null estimator.

  Verify that both modules pass. Then confirm by hand that a stray `examples/trajectory_power_study/tmp.json` absent from both fixtures fails the coverage test, and delete it.
- [x] 2.4 Add `tests/test_post_floor_config_identity.py`, checking identity against the new fixture. Add `tests/test_design_guidance_profiles.py`, covering the spec scenarios:
  - 16 core columns and 11 arm columns, each with one anchor and one cell per (mode, nonzero effect);
  - mode and effect lists per profile;
  - columns differing only in `layers` draw the same generator seeds at every replicate index;
  - every arm point differs from the baseline in exactly one setting, with the three δ counting as one, and the δ points set all three equal;
  - explicit `net_displacement`, `min_components: 2` and `joint`;
  - enumeration passes the headroom check.

  Verify that both modules pass.
- [x] 2.5 Document the profiles and the pin policy:
  - `examples/trajectory_power_study/README.md`: one entry per profile;
  - `src/motco/simulations/study/README.md`: post-floor configs name their estimator and floor and are pinned in `post_floor_config_identity.json`;
  - the `study/` paragraph of `AGENTS.md`: replace "every committed study config pins `"min_components": 1`… covers all 19 configs" with the historical/post-floor wording.

  Verify with `uv run mkdocs build --strict`.

## 3. Stage-signal table

- [x] 3.1 Add `resolve_stage_signal_by_design_point(records)` to `src/motco/simulations/study/spectrum.py` (design D5) and export it from `motco.simulations.study`:
  - bucket records exactly as `resolve_operating_by_design_point` does;
  - per (design point, mode, effect), report `n_replicates`, `n_auroc_available` and the mean, `q33`, median and `q67` of `integration_metadata["cv_mean_auroc"]`;
  - return an empty frame when there are no `power_design` records.

  Verify with unit tests on synthetic records:
  - every coordinate, the baseline column and every anchor have rows;
  - forced-rank, `concat` and `pca` columns report zero available and empty summaries;
  - a study without a design grid yields an empty frame.
- [x] 3.2 Add the frame to the report frames and write `design_point_stage_signal.csv` only when it is non-empty, in both `report.py` and `scripts/motco_study.py report`. Verify:
  - `tests/test_study_design_point_report.py` finds the file beside `design_point_operating.csv` for a design-grid study, and finds no file without one;
  - reporting the 1.2 study again gives sha256 values equal to the recorded ones for every previously written file, with `design_point_stage_signal.csv` the only new file.
- [x] 3.3 Document the new output: add it to the report-outputs list in `src/motco/simulations/study/README.md` and to the `study/` paragraph of `AGENTS.md`. Verify with `uv run mkdocs build --strict`.

## 4. Design-guidance figure

- [x] 4.1 Write `scripts/design_guidance_figure.py --magnitude <report dir> --geometry <report dir> --out-dir <dir>` (design D6). It writes:
  - `design_guidance_magnitude.png`, `design_guidance_orientation.png` and `design_guidance_shape.png`: panels for stages {3, 4} × blocks {2, 3}, n on a log axis, one line per nonzero effect, ±1 SE bars;
  - `design_guidance_control.png`: the anchor and translation's largest effect per statistic, with the α + 2·√(α(1 − α)/R) bound;
  - `design_guidance_power.csv`, holding exactly the plotted rows.

  It exits non-zero, naming the axis and writing nothing, when a report lacks `generator.n_samples`, `generator.n_stages` or `evaluation.integration_params.layers`. Verify with `tests/test_design_guidance_figure.py` on synthetic operating tables:
  - the files are written;
  - the CSV rows equal the plotted rows;
  - the two-block and three-block columns are never merged;
  - a missing axis exits non-zero without writing.
- [x] 4.2 Document the script's usage in `src/motco/simulations/study/README.md`. Verify that the documented command runs on the 4.1 test fixtures.

## 5. Local smoke run

- [ ] 5.1 Make a scratchpad copy of each profile at 2 replicates × 19 permutations (never committed). Run each through `scripts/run_study_shard.py` (one shard) → `motco_study.py merge` → `motco_study.py report` → `scripts/design_guidance_figure.py`. Verify:
  - every report writes `design_point_operating.csv` and `design_point_stage_signal.csv`;
  - the figure script writes its four figures and CSV;
  - every failed replicate is listed per column, especially n = 80 with four stages.

  Any engine failure found here is fixed in this change, with a test.
- [ ] 5.2 Settle the anchor reading and the compute plan:
  - record whether the magnitude and geometry profiles' duplicated anchors share RRPP draws (identical p-values) or only data (design D2, open question);
  - time units at each n and block count;
  - compute shard counts and a core-hour estimate.

  Present the estimate to the user before any cluster submission. Verify that the timing table and the estimate are in the change's notes.

## 6. Pilot run on ing

- [ ] 6.1 After the user's approval, submit the four arrays with `scripts/motco_study_array.sbatch` and the shard counts from 5.2, splitting n = 1000 into its own array if the estimate calls for it. Merge each profile. Verify:
  - every shard completed;
  - the merged record counts equal 9,800 / 16,200 / 6,800 / 11,200 minus recorded failures;
  - no `(cell_id, replicate_index)` is duplicated.
- [ ] 6.2 Report each profile and run the figure script on the two core reports. Commit `results/design-guidance-pilot-<date>/<profile>/report/`, `PROVENANCE.txt` per profile, and `results/design-guidance-pilot-<date>/figures/`. Verify that each `PROVENANCE.txt` names the code revision, config, software versions, record and failure counts, and reproduction commands.

## 7. Findings report and roadmap

- [ ] 7.1 Write `docs/reports/design-guidance-pilot-<date>.md`. It contains:
  - provenance;
  - per-column operating characteristics with Monte Carlo SE;
  - each column's stage signal (CV AUROC) and selected rank;
  - every column's anchor rates and anchor-decision status;
  - the figures;
  - the translation `angle` rate read against the bound;
  - the planned-grid mapping table (every factor of the original grid with its disposition and reason);
  - the proposed paper-grade slice;
  - limitations: pilot precision, block count confounded with feature count, structural-analogue blocks.

  Verify that the report has a section for each scenario of the "Dated design-guidance pilot report" requirement.
- [ ] 7.2 Update `docs/roadmap.md`:
  - add a Next changes entry for the design-guidance pilot, its outcome and the follow-ups;
  - note it under Phase 5's sensitivity deliverables (sample size, stages, blocks, signal density);
  - update the `study/` paragraph of `AGENTS.md` with the pilot's result line.

  Verify with `uv run mkdocs build --strict`.
- [ ] 7.3 Run the pre-commit gate (`uv run ruff check src/ tests/`, `uv run mypy src/motco/`, `MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`) and verify all three pass.

## Workflow follow-up

- Archive the change after the report is merged.
- Follow-up changes:
  - group-specific dispersion in the generator and its null study (change B);
  - the paper-grade rerun of the slice the report proposes;
  - the Type I check at the case study's configuration (floor 2 plus six covariates).
