# Design

## Context

See `proposal.md` (Why) for the motivation and `specs/simulation-presentation-figures/spec.md` for the requirements.

**What is in git.** Every run commits its aggregates: per-cell rejection rates, Monte Carlo SEs, bounds, verdicts, eigengap, null width, selected rank and stage AUROC. Most script-driven studies also commit their p-values. Harness replicate records (`shard_*.jsonl`, `merged.jsonl`) are gitignored. Some exist on the workstation only, such as `phase5-2026-09-10` (390 MB) and `phase6-small-n-2026-10-09` (274 MB). Every figure in the set can be drawn from the committed aggregates.

**The measurement changed under the evidence.**

```
            Aug        Sep             Oct                          Nov
orientation PC1 ───────────────────────────┤ net displacement (10-15)
magnitude   'all' ─────────────────────────┤ joint (10-09)
rank rule   CV ≥ 1 ───────────────────────────────────┤ CV ≥ 2 (10-23)
null        plain RRPP ─────────────────────────────────────┤ df opt-in (10-30)
```

The design-guidance pilot and slice, the case-study checks, the corrected permutation, the ADNI design and the adjusted attribution all ran under the current measurement. The n = 1,200 evidence predates the floor. `net_displacement_power_phase5.json` (2026-10-14) differs from the current measurement only in the floor: it ran net displacement and joint magnitude over all four modes at 500 × 999.

**Existing code to reuse.** `simulations.showcase.generate_showcase_datasets` generates one seeded dataset per mode from the numpy generator. `results/design-guidance-slice-2026-10-29/figures/design_guidance_power.csv` is the slice's joined power table, written by `scripts/design_guidance_figure.py`.

## Goals / Non-Goals

**Goals:**
- Each figure carries one message, can stand alone on a slide, and says which measurement it shows.
- The set rebuilds from a clean clone in under a minute, with no cluster, no records and no R.
- One visual system across the set, distinct from the case-study repositories' figures.

**Non-Goals:**
- Paper-grade figures. This set is for talks. A manuscript set would revisit layout and panel composition.
- Changing `study/report.py` or any committed `report/` output.
- Real-data figures, interactive figures and dark mode.
- New simulation runs. The one exception is D5's floor-equivalence check, which fits latent spaces on rebuilt Phase 5 datasets and runs no RRPP.

## Decisions

### D1. One script, not package code

`scripts/presentation_figures.py` holds the style, the source registry, one adapter per source, one builder per figure, the writer and the CLI (`--out-dir`, repeatable `--figure`, `--no-footer`, `--titles`, `--seed`). Tests import it the way `tests/test_adni_design_check.py` imports its script.

*Alternatives.* A module under `src/motco/` would make presentation styling and hard-coded `results/` paths part of the package API. Extending `study/report.py` would couple talk figures to the byte-identity of every committed `report/`. A notebook cannot be tested or rendered from the command line.

### D2. A source registry with validating adapters

Each source is one registry entry: path, run name, measurement tag, required columns, and the row filter that selects what is plotted. An adapter reads the file, checks the columns, applies the filter, checks that every expected reading is present and finite, and returns a tidy frame. That frame has one row per plotted value, with `rate`, `se`, `bound`, `verdict`, `n_replicates` and the measurement tag. Builders see only tidy frames. Schema drift therefore fails in the adapter, naming the file and the missing column or reading, never as a silently thinner plot.

*Alternative:* reading the CSVs ad hoc inside each builder. That is fewer lines, but a filter that matches nothing draws an empty panel without complaint, which is how `design_point_power.png` ended up blank.

### D3. Style

| Role | Value |
|---|---|
| Surface | `#fcfcfb` (near-white) |
| Ink primary / secondary / muted | `#1d1d1b` / `#57564f` / `#8a8980` |
| `delta` (magnitude mode) | `#5a4fcf` violet |
| `angle` (orientation mode) | `#e0712c` orange |
| `shape` (shape mode) | `#1a9b7f` teal |
| Controls (anchor, translation), group A | muted ink |
| Heatmap ramp | one neutral slate hue, light → dark, values printed in every cell |

**Validation.** The statistic triplet was validated on 2026-11-02 with the dataviz skill's `validate_palette.js`, using `--mode light --pairs all` against `#fcfcfb`. All checks passed:
- lightness in band;
- chroma at least 0.1;
- worst colour-vision-deficiency ΔE 10.1, teal against orange, protan;
- worst normal-vision ΔE 24.7;
- contrast at least 3:1 for all three colours.

Three other triplets were rejected: one failed the lightness band, one failed deficiency separation at ΔE 4.3, and one passed only in the 6–8 warn band. A test pins these hex values, so a palette change has to pass through this record.

**Secondary categories** never get a hue (spec, *Modes take their target statistic's colour*):

| Category | Encoding |
|---|---|
| Blocks: 2 / 3 | solid / dashed line |
| Estimator: PC1 / net displacement | dashed with open markers / solid with filled markers |
| Correction: uncorrected / `df` | open / filled markers |
| Effect size | marker shape plus a direct label |
| Control: anchor / translation | filled / hollow marker |

**Canvas and type.**
- Canvas 13.333 × 7.5 in, written as PNG at 300 dpi and as PDF.
- Sizes: ticks 15 pt, axis labels 18 pt, panel titles 18 pt, annotations 14 pt, footer 12 pt.
- Font family Inter, then Source Sans 3, then DejaVu Sans, with missing-font warnings silenced.
- Left and bottom spines only, and a light grid on the value axis.

No in-figure title by default, since the slide carries the message; `--titles` adds it. Writes are atomic (temporary file, then `os.replace`), so a failing figure leaves no partial file.

### D4. Figure contents

| Figure | Form | Source(s) | Selection |
|---|---|---|---|
| `mode_geometry` | four small panels, A in muted ink, B in the mode's colour, arrows per segment, a "reads: δ / θ / shape / none" chip | constructed coordinates | three stages. Scale ×1.6 about stage 1; rotate 35° about stage 1; displace the interior stage; translate by a fixed vector |
| `mode_simulated` | anchor + four modes: group stage-mean paths over faint samples in a display projection (per panel, the first two principal axes of the six group × stage means; translation's second axis is the group offset), "display projection, not the measurement" | `generate_showcase_datasets(seed=0, n_samples=1000, n_stages=3, effect_size=1.0)` | combined file + one file per mode |
| `size_atlas` | dot plot, rows = design × control grouped into sections, three statistic markers per row, a bound tick per row, α line | see D6 | production arm per source |
| `specificity_n1200` | annotated heatmap, rows = anchor + four modes, columns = statistics, diagonal outlined | `net-displacement-power-2026-10-14/phase5/report/design_point_operating.csv` | net-displacement column (`power_design`), e = 1.0; anchor e = 0 |
| `specificity_n80` | same layout and scale | slice `core_magnitude` + `core_geometry` `design_point_operating.csv` | n = 80, layers `(methylation, expression)`; magnitude e = 0.1, geometry e = 1.0 (effect printed per row) |
| `power_vs_n` | overview: three panels (statistic), log-n axis, solid 2 / dashed 3 blocks at the top effect, ±1 SE. Details: per statistic, a 2- and a 3-block panel, one line per effect | `design-guidance-slice-2026-10-29/figures/design_guidance_power.csv` | diagonals, `figure` ∈ {magnitude, orientation, shape} |
| `power_vs_stage_signal` | three panels sharing a log-δ axis: target power per statistic, median CV rank, median stage AUROC | slice `stage_effect_*` (δ 0.25 / 0.5 / 1.0) + slice core n = 200, two blocks (δ 2.0); `design_point_operating.csv` + `design_point_stage_signal.csv` | top effects: magnitude e = 0.1, geometry e = 1.0 |
| `estimator_comparison` | two panels (Phase 6 n = 80, Phase 5 n = 1,200), effect on x: orientation/`angle` in orange and the shape-mode `angle` leak in teal, PC1 dashed / ND solid, α line | `net-displacement-power-2026-10-14/{phase6,phase5}/report/design_point_operating.csv` | PC1 = `power_primary`, ND = `power_design` |
| `covariates_null` | categorical x = covariate columns 1 / 3 / 6; faint readings (2 strengths × 2 controls) and a bold pooled line per correction; pooled bound α + 2·√(α(1 − α)/2000) = 0.0597 | `case-study-diagnosis-2026-10-30/delta{0.1,0.05}/diagnosis_rates.csv` | arms `z_noise{0,2,5}` and their `_corrected` twins, `delta` (columns = 1 + noise count) |

The per-statistic power-vs-n detail files restyle what `design_guidance_figure.py` already draws. The overview and every other figure are new.

*Amended during implementation (2026-11-02).* The display projection of `mode_simulated` was first planned as a 2-component PLS-DA on stage. Rendered, it put translation's group B exactly on group A, and it barely showed the orientation rotation: a stage-supervised projection keeps only the stage directions, and both realized effects lie outside them. The figure now projects onto the first two principal axes of the six group × stage means. Translation's offset is orthogonal to the stage geometry (it is the third axis), so that panel keeps the first axis and uses the group offset, orthogonalized against it, as the second. The axis label says so, and so does the caveat. Nothing is measured in this projection.

### D5. The rank floor at n = 1,200

**Equivalence argument.** The pooled PLS is refit once at the final rank. Equal final ranks under floor 1 and floor 2 therefore mean identical latent scores, statistics and p-values. The records keep only the final rank, which is 2–4 in all 34 cells × 500 replicates of the n = 1,200 run and never 1. A fold-level choice of rank 1 that was outvoted would leave no trace, and under the floor that fold would choose another rank.

**Check.** 152 datasets (19 baseline cells × replicates 0–7) were rebuilt from the config's seeds and fitted at floor 1 and floor 2 (BLAS pinned, 8 workers, 1 h 43 min on 2026-11-02). Result: **every rank agrees.** 0 of 152 selected ranks change: rank 3 in 138 datasets and rank 2 in 14 under both floors. The mean CV AUROC is also equal in all 152, so no outer fold picked rank 1 either. The one-sided 95% upper bound on the mismatch rate is 0.0195. The committed script reproduces the first dataset exactly.

**Tagging decision.** The Phase 5 sources (`net-displacement-power-2026-10-14/phase5`, `net-displacement-type-i-2026-10-13/phase5`) are tagged as the current measurement: "CV rank ≥ 2 by floor check (0/152 ranks change)". The tag is derived from the committed `summary.json` at render time, so a summary with any mismatch, or none at all, reverts the tag to pre-floor.

- If every rank agrees, the Phase 5 sources are tagged as the current measurement, citing the check and its sample. Zero mismatches in 152 gives a 95% upper bound of 3/152 ≈ 2% on the mismatch rate.
- If any rank differs, they are tagged pre-floor.

The check is committed as `scripts/rank_floor_equivalence_check.py`, with its per-replicate CSV and a summary under `results/rank-floor-equivalence-2026-11-02/`, so the tag cites a file in git. The first run was made from a scratch copy of the same code at `b7df7c7`. PROVENANCE records that copy and its SHA-256, so it is not rerun.

**Phase 6.** The Phase 6 (n = 80) sources in `estimator_comparison` stay tagged pre-floor. The comparison is between estimators, and no current-measurement claim rests on it.

### D6. Size atlas rows

| Section | Source | Readings |
|---|---|---|
| Balanced, 3 stages, no covariates | slice `core_geometry/report/anchor_decision.csv` | n ∈ {80, 200, 500} × {2, 3} blocks; anchor + translation e = 1.0 |
| Weaker stage signal (n = 200) | slice `stage_effect_geometry/report/anchor_decision.csv` | δ ∈ {0.25, 0.5, 1.0}; anchor + translation e = 1.0 |
| n = 1,200, 4 stages, 3 blocks | `net-displacement-type-i-2026-10-13/phase5/report/anchor_decision.csv` | the column `pls`, CV rank, `net_displacement`; anchor + translation; tag per D5 |
| SEA-AD astrocyte table, weak stage signal, no covariates | `rank-floor-2026-10-23/delta{0.1,0.05}/rank_floor_rates.csv` | floor-2 arm; anchor + translation |
| SEA-AD case-study configuration (6 covariate columns, `df`) | `corrected-permutation-2026-10-30/size_readings.csv` | `tier = reading`, `arm = df`; four Type I runs + the power run's controls |
| ADNI design (n = 1,501, 6 covariate columns, uncorrected, κ ≈ 1.003) | `adni-design-2026-11-02/run/size_readings.csv` | anchor + translation |

Each reading plots against the bound its source recorded. The corrected-permutation per-reading bound (0.0792) differs from the R = 500 bound (0.0695), and both are shown as recorded.

*Amended during implementation.* The two slice sections read `anchor_decision.csv` rather than `design_point_operating.csv`. The rates are the same, but only `anchor_decision.csv` records the anchor's bound and classification, which the spec asks the atlas to show. That table records neither for translation, so translation's bound is the same α + 2·SE_α at its own R and its verdict is derived against it. The manifest's caveats say so. The n = 1,200 section already read this file.

Excluded as well: the Phase 6 Type I profile, which is pre-floor with no equivalence check. The same SEA-AD table appears under the current measurement in the rank-floor and corrected-permutation rows.

Excluded, with the reason in the manifest's caveats:
- `case-study-type-i` and `covariate-type-i`: superseded by the `df` configuration and pre-floor respectively. Their uncorrected excess is what `covariates_null` shows.
- The latent-adjustment arms: a sensitivity analysis, not production.
- The PC1 columns.

### D7. Manifest and captions

`manifest.json` holds:
- per figure: `id`, `message`, `files`, `sources` (each with `path`, `sha256`, `run`, `measurement`) and `caveats`;
- for the whole set: `git_revision`, `worktree_clean`, `motco_version`, `matplotlib_version`, `rendered_at`, `seed`, and the palette with its validation record.

`captions.md` renders the same content as speaker notes: per figure, a title, the message, the measurement and the caveats.

### D8. Tests

`tests/test_presentation_figures.py` uses the Agg backend and covers:
- every figure renders into `tmp_path` from the committed sources;
- an unknown figure name exits non-zero before any file is written;
- a source with a dropped column, or a missing expected reading, raises with the file named, and leaves no output file;
- the value scenarios in the spec:
  - specificity cells equal their source;
  - the estimator top-effect values (0.300 / 0.648 / 0.840 / 0.924);
  - the pooled covariate rates (0.0545 / 0.066 / 0.0545);
- the palette constants equal D3's validated values;
- `mode_simulated` is deterministic, at a reduced n for speed;
- the manifest carries a SHA-256 per source and the dirty flag.

They run in the fast suite, with no `slow` mark.

## Risks / Trade-offs

- **The atlas has about 36 rows on one slide** → grouped sections and anchor/translation folded into one row per design (filled and hollow markers) bring it to about 18 rows. If it is still unreadable at slide size, split it into a no-covariate figure and a covariate figure. Task 8.2 is the eyeball check.
- **Effect axes differ between designs** (magnitude e = 1.0 at n = 1,200 against 0.1 at n = 80) → each row prints its effect. The message is the off-diagonal mixing, not the magnitude row.
- **The floor check samples 152 of about 9,500 replicates** → the tag states the sample and the bound. A mismatch reverts the tag; it does not change any plotted value.
- **Fonts differ between machines** → fallback chain; the committed set is rendered once, on one machine, from a clean worktree.
- **The figure set and the committed `report/` PNGs show the same numbers differently** → the manifest names the source run of every value, and the reports remain the record.

## Migration Plan

Additive. Nothing existing changes. To roll back, delete the script, its test and the rendered directory.
