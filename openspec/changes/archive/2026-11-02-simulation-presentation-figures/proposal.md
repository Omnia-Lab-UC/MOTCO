# Proposal

## Why

The simulation work needs a set of standalone figures to present as preliminary results, first at a symposium on 2026-11-06. What exists today does not serve that. `study/report.py` writes 114 per-run diagnostic PNGs. Their colours carry no stable meaning, their axes use internal names (`trajectory_mode`, `n_samples`), and some of them are blank placeholders (`net-displacement-power-2026-10-14/phase6/report/design_point_power.png`). Fourteen script-driven studies, among them the covariate, permutation-correction, ADNI and attribution checks, have no figure at all. The measurement also changed four times during the work: PC1 to net displacement, `all` to `joint` magnitude, the rank floor, and the opt-in `df` correction. A figure that does not say which measurement it shows can mislead.

## What Changes

- A new script, `scripts/presentation_figures.py`, renders nine standalone figures, each carrying one message:
  1. `mode_geometry`: a conceptual schematic of what magnitude, orientation, shape and translation do to a two-group trajectory, and which statistic reads each.
  2. `mode_simulated`: the same four modes realized by the numpy generator, generated live from a fixed seed and shown in a display projection. One combined file and one file per mode, for slide build-up.
  3. `size_atlas`: every zero-effect and translation control rate, from n = 80 to n = 1,501, against its own α + 2·SE bound.
  4. `specificity_n1200`: the mode × statistic response at n = 1,200.
  5. `specificity_n80`: the same matrix at n = 80, where the statistics mix.
  6. `power_vs_n`: diagonal power against sample size, two blocks against three, plus one detail file per statistic.
  7. `power_vs_stage_signal`: diagonal power, CV-selected rank and stage AUROC as the joint stage effect δ weakens.
  8. `estimator_comparison`: PC1 against net displacement, for orientation power and the shape → `angle` leak, at both designs.
  9. `covariates_null`: covariate-adjusted `delta` size against the number of covariate columns, uncorrected against `df`-corrected.
- Figures 1, 3–6, 8 and 9 were shortlisted during exploration. Figure 2 splits the shortlisted mode figure into concept and simulation, and figure 7 is added because both case studies have weak stage signals.
- Every figure except figure 2 reads only committed aggregate files under `results/`. It does not read gitignored records, and it does not run new simulations or use the cluster. Figure 2 generates one seeded dataset per mode in seconds.
- Every figure carries a measurement tag (estimator, rank rule, magnitude construction, permutation correction, replicates × permutations) and a source line, both in a removable footer and in a manifest.
- A shared presentation style: a light surface, a fixed and validated colour for each statistic, each mode drawn in the colour of the statistic it targets, neutral controls, non-colour encodings for secondary categories, a 16:9 canvas and slide-legible type. Output is PNG at 300 dpi plus PDF.
- A `manifest.json` and a `captions.md` record, for each figure, its message, files, source paths with SHA-256, measurement tags and caveats, along with the git revision and the MOTCO version.
- The rendered set is committed under `results/presentation-figures-<render date>/`.
- A small committed check, `scripts/rank_floor_equivalence_check.py`, with its result under `results/rank-floor-equivalence-2026-11-02/`. It rebuilds sampled n = 1,200 datasets and compares the selected rank at floor 1 and floor 2, so the n = 1,200 figures can say whether they show the current measurement. It fits latent spaces only, runs no RRPP, and changes no existing result.
- Nothing in `src/motco/` changes. No existing output changes.

## Capabilities

### New Capabilities
- `simulation-presentation-figures`: a curated, provenance-tagged set of standalone figures that present the simulation evidence (operating characteristics, specificity, design guidance and the measurement decisions) built from committed aggregates.

### Modified Capabilities

None. Figure 2 calls the existing `trajectory-mode-showcase` generator (`generate_showcase_datasets`) without changing its contract. The design-guidance figure script stays as specified, and the new figure 6 may import its table loaders.

## Impact

- New: `scripts/presentation_figures.py`, `scripts/rank_floor_equivalence_check.py`, `tests/test_presentation_figures.py`, `results/presentation-figures-<date>/`, `results/rank-floor-equivalence-2026-11-02/`.
- Read-only inputs: committed CSV/JSON under `results/` for the design-guidance slice, the net-displacement Type I and power runs, the rank-floor check, the case-study diagnosis, the corrected permutation and the ADNI design.
- `AGENTS.md` gains one line naming the script.
- Dependencies: matplotlib, already required (`>=3.7`). No new packages.
- The real-data figures are out of scope. They live in `sead_sexdiff` and `adni_sexdiff`.
