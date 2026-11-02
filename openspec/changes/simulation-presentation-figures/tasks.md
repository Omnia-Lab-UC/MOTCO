# Tasks

## 1. Rank-floor equivalence check at n = 1,200 (design D5)

- [x] 1.1 Add `scripts/rank_floor_equivalence_check.py`, the scratch probe's logic: `--config`, `--replicates`, `--workers` and `--out-dir`. It writes a per-replicate CSV (cell, phase, mode, effect, replicate, rank and mean CV AUROC at floor 1 and floor 2) and a summary JSON (count, mismatches, 95% upper bound). Verify with a unit test of the summary on fabricated rows, including one mismatch.
- [ ] 1.2 Commit the 152-dataset result under `results/rank-floor-equivalence-2026-11-02/`, converted from the scratch run's JSONL, with a PROVENANCE that records the revision (`b7df7c7`), the scratch script's SHA-256, the thread pinning and the wall time. Verify that the CSV's row count and mismatch count equal the scratch output.
- [ ] 1.3 Fill in design D5's pending result and its tagging decision, and add one line for the check to `AGENTS.md`. Verify that `openspec validate simulation-presentation-figures` passes.

## 2. Script skeleton, style and provenance (design D1, D3, D7)

- [x] 2.1 Create `scripts/presentation_figures.py` with the figure registry and the CLI: `--out-dir`, repeatable `--figure`, `--no-footer`, `--titles` and `--seed`. An unknown name exits non-zero, lists the valid names and writes nothing. Verify with a test of the unknown-name exit and of a one-figure selection.
- [x] 2.2 Add the style layer: D3's palette, ink and surface constants, the 16:9 canvas, the type sizes, the font fallback, the footer and the atomic PNG (300 dpi) + PDF writer. Verify with tests that pin the palette hex values to D3's validated set, and that check a failing builder leaves no file.
- [x] 2.3 Write `manifest.json` and `captions.md`. Each source carries its SHA-256, run name and measurement tag. The set records the git revision, a clean/dirty flag, the MOTCO and matplotlib versions, the render time, the seed and the palette's validation record. Verify with a test that a rendered figure's manifest entry has a SHA-256 per source and the dirty flag.
- [x] 2.4 Add the script to `AGENTS.md` (Architecture, beside `scripts/design_guidance_figure.py`). Verify that the line names the script, its inputs (committed aggregates only) and its output directory.

## 3. Source registry and adapters (design D2, D6)

- [x] 3.1 Add the registry and an adapter for each source: slice `core_*` and `stage_effect_*` operating and stage-signal tables, `design_guidance_power.csv`, the net-displacement power and Type I tables (Phase 5, Phase 6), `rank_floor_rates.csv`, the corrected-permutation `size_readings.csv`, the ADNI `size_readings.csv` and the diagnosis `diagnosis_rates.csv`. Each returns D2's tidy frame. Verify with a test that every adapter loads its committed file with no missing reading.
- [x] 3.2 Make each adapter fail on a missing file, a missing column, a missing expected reading or a non-finite rate, naming the file and the item. Verify with parametrized tests on modified copies in `tmp_path`.
- [x] 3.3 Attach measurement tags per D5. Phase 5 sources are tagged current or pre-floor according to the committed check, and Phase 6 net-displacement-power sources are tagged pre-floor. Verify with a test of the tag each Phase 5 and Phase 6 source carries.

## 4. Concept and simulated mode figures

- [x] 4.1 Build `mode_geometry` from constructed coordinates: four panels, B in the mode's colour, A in muted ink, statistic chips, and translation reading "none". Verify with a test that it renders and that each B path's colour equals its target statistic's colour (translation: muted ink).
- [x] 4.2 Build `mode_simulated` from `generate_showcase_datasets` with a fixed seed: the combined file and one file per mode (anchor + four), stage-mean paths over faint samples, labelled as a display projection. Verify with a test that it is deterministic at a reduced n, and that it writes six figures.

## 5. Size atlas (design D6)

- [x] 5.1 Build `size_atlas` as a sectioned dot plot: one row per design, anchor filled and translation hollow, three statistic markers per row, the source's bound tick per row, an α line, and an exceeding reading drawn distinct and labelled. Verify with tests that the plotted readings equal D6's selection (count and values) and that an injected exceedance gets the distinct marker and label.

## 6. Specificity and design-guidance figures (design D4)

- [x] 6.1 Build `specificity_n1200` and `specificity_n80` as annotated heatmaps with a shared 0–1 scale, the diagonal outlined and the effect printed per row. Verify with a test that every cell equals its source rate and that both figures share the colour normalization.
- [x] 6.2 Build `power_vs_n`: an overview (three panels, 2 blocks solid and 3 dashed, ±1 SE, log-n axis) and three per-statistic detail files. Verify with a test that it writes four files and that the overview's plotted values equal `design_guidance_power.csv`'s top-effect diagonals.
- [x] 6.3 Build `power_vs_stage_signal`: power, median CV rank and median stage AUROC in three panels on a shared log-δ axis, with δ 2.0 taken from the slice core. Verify with a test that it has three panels and three separate y-axes (no twin axes), and that the δ = 2.0 values equal the core's n = 200 two-block rows.

## 7. Measurement-decision figures (design D4)

- [x] 7.1 Build `estimator_comparison`: Phase 6 and Phase 5 panels, orientation/`angle` and the shape-mode `angle` leak, PC1 dashed with open markers and net displacement solid with filled markers, an α line, and tags per panel. Verify with a test of the spec's top-effect values (0.300 / 0.648 / 0.840 / 0.924).
- [x] 7.2 Build `covariates_null`: covariate columns 1/3/6 on the x axis, faint individual readings and bold pooled lines for uncorrected and `df`, and the pooled bound 0.0597. Verify with a test of the pooled values (0.0545 / 0.066 / 0.0545).

## 8. Render, review and commit

- [ ] 8.1 Run the pre-commit gate (`uv run ruff check src/ tests/ scripts/presentation_figures.py scripts/rank_floor_equivalence_check.py`, `uv run mypy src/motco/`, `MOTCO_TEST_PERMS=99 uv run pytest tests/ -m "not slow"`) and confirm it passes.
- [ ] 8.2 Render the full set from a clean worktree into `results/presentation-figures-<date>/`. Open every PNG at slide size and check for label collisions, overflow and legibility of the smallest text, splitting the size atlas if needed (design Risks). Verify that the manifest records `worktree_clean: true` and that every figure listed in the spec is present.
- [ ] 8.3 Commit the rendered set, its manifest and its captions. Verify that `captions.md` reads as speaker notes for each figure: message, measurement and caveats.

## Workflow follow-up

- Sync the delta spec into `openspec/specs/simulation-presentation-figures/` and archive the change after review.
