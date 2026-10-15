# real-data-trajectory-analysis

## Why

The roadmap's next item is the Phase 6 real-data case study (SEA-AD middle temporal gyrus astrocytes, RNA +
ATAC donor pseudobulks, sex × three ADNC stages, n = 80) in the external
[Omnia-Lab-UC/sead_sexdiff](https://github.com/Omnia-Lab-UC/sead_sexdiff) repository. That repository today
contains only data download, subsetting, pseudobulk and DESeq2 tooling; it does not import MOTCO at all, so
the trajectory analysis has still to be written. MOTCO, however, exposes the **validated** pipeline — pooled
block standardization, stage-supervised PLS-DA with the double-CV rank rule, trajectory statistics with the
net-displacement orientation estimator, RRPP, the selected-rank covariate, null summaries, the stage-mean
spectrum and frozen-model attribution — only through `evaluate_semisynthetic_trajectory`, which is bound to
`SemiSyntheticTrajectoryDataset`: three fixed InterSIM block names, a hard-wired logit on the block called
`methylation`, generator truth, and stage labels that sort correctly because they are integers. None of
that holds for RNA + ATAC pseudobulks with stages named *Not AD + Low / Intermediate / High*.

The alternative surface, `motco plsr` + `motco de`, is a *different* pipeline from the one the operating
studies characterized: it standardizes but never applies the M-value contract, picks the rank by a plain
mode of the CV table (no parsimony tie-break), fits scores through a separate code path, measures `angle`
with PC1 (the estimator the paper will not describe) and records neither the selected rank, the null
summary nor the spectrum. The roadmap's case-study item already warns that primitives callers must pass
`orientation_estimator="net_displacement"` by hand. Without a real-data entry point the case study would
have to re-implement the harness glue from the simulation code and prove its own equivalence, which is
exactly the drift the roadmap's interpretability statement forbids. This change gives the case study the
same measurement the paper validates, and it is the only MOTCO-side blocker left before Next changes item 1
can start.

## What Changes

- **New `motco.analysis` module (real-data entry point).** `analyze_trajectories(dataset, params)` runs
  the production pipeline on caller-supplied data: an ordered set of named omic blocks (samples × features
  DataFrames, any names, any count ≥ 1), each with a declared per-block unit transform (`identity`, or
  `beta_to_mvalue` for methylation beta values), sample metadata with a sample-id column, a two-level group
  column and a stage column, and an **explicit stage order** (first and last stage define net displacement,
  so the order is load-bearing and is never inferred by sorting). It fits pooled per-feature scalers,
  concatenates the blocks in the declared order with block-prefixed feature names, builds the PLS latent
  space with the harness's double-CV rank rule (or a `forced_components` rank for the dimensionality
  sensitivity Phase 6 requires), builds the group × stage design, measures `delta`/`angle`/`shape` with
  net displacement as the default estimator, runs RRPP, and returns the same measurement record the harness
  produces: pair statistics, p-values, null summary, latent stage-mean spectrum, selected rank and
  component-selection mode, CV metadata, estimator, group/stage levels, software versions. Optional
  frozen-model orientation attribution with bootstrap stability reuses the existing attribution routine;
  the truth-recovery block is reported unavailable rather than invented. The result carries the latent
  scores and PLS loadings so the fitted representation can be written down and reused for figures.
- **Writer for a versioned analysis artifact.** `write_analysis_outputs(result, out_dir)` writes
  `result.json` (the JSON-safe record plus the resolved parameters), `latent_scores.csv`,
  `pls_loadings.csv`, `stage_means.csv` (LS-mean vectors in the latent space, for `plot_trajectories`) and,
  when attribution ran, the attribution tables. This is the first instance of the roadmap's cross-cutting
  "preserve fitted preprocessing, PLS model, component count, feature order and loadings" item: everything
  needed to reproduce the figures and to audit the fit, without pickling an sklearn estimator.
- **New `motco analyze` CLI subcommand** over the same function: repeatable `--block name=path.csv`,
  `--transform name=kind`, `--metadata`, `--sample-col`, `--group-col`, `--stage-col`, `--stage-order`,
  `--permutations`, `--seed`, `--forced-components`, `--orientation-estimator` (default
  `net_displacement`), attribution flags, `--out-dir`. `ValueError`s become clean `SystemExit` messages like
  the other subcommands.
- **One implementation for simulation and real data.** The harness's PLS integration and block
  preprocessing are refactored into dataset-free cores that both `evaluate_semisynthetic_trajectory` and
  `analyze_trajectories` call, so the production measurement exists once. A regression test feeds one
  semi-synthetic dataset through both entry points and requires identical latent matrices, statistics,
  p-values, selected rank and null summaries. Every committed byte-identity fixture and every committed
  config's cell ids, signatures and matched seeds are unchanged.
- **`motco de` gains `--orientation-estimator`** (`pc1` default, so existing invocations are
  byte-identical; `net_displacement` selectable). The stats primitives' defaults do not move.
- **Documentation.** `docs/api/analysis.md` (new page in the API nav), CLI reference, README quick start
  pointer, `CLAUDE.md`, and `docs/roadmap.md`: Next changes item 1 names `motco analyze` as the entry point
  the case study must use, and the Phase 6 "Work" list records which bullets this change covers.
- **Version** bumps to `0.7.0` (new public module and subcommand). Tagging the release so the case study
  can pin it is left to the maintainer.

## Capabilities

### New Capabilities

- `real-data-trajectory-analysis`: the real-data entry point — input contract (named blocks, declared unit
  transforms, explicit stage order, two groups), the measurement it performs (pooled standardization,
  stage-PLS double-CV latent space, net-displacement trajectory statistics, RRPP, recorded covariates),
  its equivalence to the simulation harness on identical data, optional attribution, the written analysis
  artifact, and the `motco analyze` subcommand.

### Modified Capabilities

- `trajectory-orientation-invariance`: `motco de` SHALL accept an orientation-estimator selection whose
  default is PC1 (unchanged behaviour) and whose `net_displacement` value measures what the harness
  measures.

## Impact

- **New:** `src/motco/analysis.py`, `docs/api/analysis.md`, `tests/test_analysis.py`,
  `tests/test_analysis_cli.py`.
- **Refactored (behaviour-preserving):** `simulations/preprocessing.py` (block-mapping core under the
  existing dataset-facing functions), `simulations/evaluation.py` (`_pls_integration` core split into a
  function of `X`, `y` and `integration_params`; the design builder and the record helpers reused). The
  byte-identity tests (`test_orientation_probe_byte_identity.py`, `test_phase6_byte_identity.py`), the
  committed-config identity tests and the R-reference regression tests gate the refactor.
- **CLI:** `cli.py` gains `analyze`; `de` gains one flag. `docs/api/cli.md` updated.
- **Docs/config:** `mkdocs.yml` nav, `README.md`, `CLAUDE.md`, `docs/roadmap.md`, `src/motco/__init__.py`
  version.
- **Out of scope:** covariate adjustment in the design (the operating characteristics were measured without
  covariates; adding them is a methodological change with its own validation), count normalization or
  feature filtering of pseudobulks (the case-study repository's raw-input contract), missing-data
  imputation (complete matrices are required and validated), SNF or other integration methods on real data
  (`pls` is the production measurement space; `concat`/`pca`/`snf` stay simulation diagnostics), changing
  the primitives' defaults, and any re-run of committed `results/` or `docs/reports/`.
- No compute is needed. No new dependency.
