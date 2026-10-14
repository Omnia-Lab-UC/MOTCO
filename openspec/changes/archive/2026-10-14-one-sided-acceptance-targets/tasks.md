# Tasks

## 1. Byte-identity baseline

- [x] 1.1 Confirm `tests/test_study_targets.py`, `tests/test_study_config.py` and `tests/test_study_end_to_end.py` pass on unchanged code; this is the guard every later task re-runs.
- [x] 1.2 Capture a pre-change fixture: run the `tests/test_study_end_to_end.py` study (Type I and specificity targets declared, no `bound` key) and store its `report/acceptance_report.csv` and `.json` under the scratchpad; also run `scripts/motco_study.py report --config examples/trajectory_power_study/angle_pivotality_diagnostic.json --out-dir results/latent-spectrum-2026-09-02` and confirm `git diff --stat results/latent-spectrum-2026-09-02/report/` is empty before any code change (the empty-table case).

## 2. Config: the `bound` key (design D1, D3, D5)

- [x] 2.1 Add `bound: str = "two_sided"` to `TypeIControlTarget` and `SpecificityTarget` in `study/config.py`, a module-level `_TARGET_BOUNDS = {"two_sided", "upper"}`, and parsing in `_build_type_i_target` / `_build_specificity_target`. Verify tests: an entry without the key loads as `two_sided`; an entry with `"bound": "upper"` loads as `upper`; a dump/reload round-trips both values.
- [x] 2.2 Reject an unknown value with a `StudyConfigError` naming `acceptance.type_i[<i>].bound` or `acceptance.specificity[<i>].bound` and the offending value. Verify one test per target kind (spec scenario *Unknown bound form is rejected*).

## 3. Evaluator: the upper branch (design D2, D4)

- [x] 3.1 Branch `_evaluate_type_i` and `_evaluate_specificity` in `study/targets.py` on `target.bound`; the `upper` branch calls `phase4.type_i_inflation_bound(alpha, se_tolerance, available_replicates)` and records `bound_kind`, `alpha_se`, signed `deviation`, `bound`, with the rationale `rate = <r> <=|> alpha + <k>·SE_alpha (<bound>)`. The `two_sided` branch is the existing code, untouched. Verify `ruff` and `mypy` pass and the import adds no cycle (`python -c "import motco.simulations.study"`).
- [x] 3.2 Add tests for the spec scenarios: a `0.000` rate over 500 replicates with `upper` is met with bound `0.0695` and deviation `−0.050` (both target kinds); the same cell under `two_sided` is not met; a `0.080` rate over 500 with `upper` is not met and the rationale cites `0.0695`; a `0.024` rate with `upper` is met. Verify `tests/test_study_targets.py` passes.
- [x] 3.3 Add a test that a default-form evaluation's `observations` dict has exactly the pre-change keys (no `bound_kind`, no `alpha_se`) and the pre-change rationale text, and that a report mixing one `upper` and one `two_sided` target writes `obs_bound_kind` / `obs_alpha_se` columns with blanks on the two-sided row. Verify it passes.

## 4. Byte identity (design Risks)

- [x] 4.1 Re-run the task 1.2 end-to-end study on the modified code and verify `acceptance_report.csv` and `.json` are byte-identical to the captured fixture; add this comparison as a test in `tests/test_study_end_to_end.py` (regenerate the report twice in one test, once from a config without the key and once from the same config with `"bound": "two_sided"` spelled out, and assert identical bytes). Re-run `scripts/motco_study.py report` on `results/latent-spectrum-2026-09-02` and verify `git diff --stat` on its `report/` is empty.
- [x] 4.2 Verify `tests/test_committed_config_identity.py` passes unchanged (targets never enter cell ids or signatures).
- [x] 4.3 Run the pre-commit gate (`uv run ruff check src/ tests/ && uv run mypy src/motco/ && MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow" --tb=short`) and confirm it is green.

## 5. Documentation

- [x] 5.1 Document `bound` in `src/motco/simulations/study/README.md`: the acceptance table row, both formulas, that `upper` mirrors the gate's control/inflation bound, which observation keys appear only under `upper`, the signed-`deviation` rule for mixed reports, and that flipping the default would require re-issuing every committed report. Verify the README renders and the key names match the code.
- [x] 5.2 Update `docs/roadmap.md`: remove item 2 from *Next changes*, record the change by name and date beside the magnitude re-measurement addendum's reader disagreement, and add to item 3 that its two configs set `"bound": "upper"` on every control target. Verify the remaining items renumber and the committed-config sentence (no committed config opts in) is present.
- [x] 5.3 Verify `openspec validate one-sided-acceptance-targets --strict` passes and the delta spec's scenarios each map to a test added in tasks 2–4.
