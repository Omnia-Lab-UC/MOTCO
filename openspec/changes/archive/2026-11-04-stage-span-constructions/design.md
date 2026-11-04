# Design

## Context

See proposal.md (Why) for the motivation. The facts below shape the approach.

- The generator samples every cell as MVN(μ_cell, Σ_omic) from population means built by `omic_population_means` (methylation in M-values before `rev_logit`; expression and proteomics in native units). It already exposes exact `population_trajectories`. Every trajectory mode is feature-set surgery on methylation indicators, with expression and proteomics re-derived through the incidence maps.
- The production measurement standardizes each feature by its pooled sample SD, concatenates the measured blocks (`integration_params["layers"]`) and fits a pooled stage PLS-DA at the CV rank (floor 2). Realized geometry already reports the `population_native`, `population_standardized`, `observed_standardized` and `pls_latent` checkpoints, but nothing measures the population truth through the fitted basis.
- Matched-seed families give every primary cell the same generator seed at a replicate index. Sibling profiles in one family evaluate the same anchor datasets, which is the design-guidance precedent: the anchors are reproduced, never pooled.
- Cell identity omits unset fields through `grid._SIGNATURE_OMIT_WHEN_NONE`. Committed configurations, fixtures and dated reports must reproduce byte-identically.
- The scratchpad prototype showed two things. Edits in a replicate's own basis are recovered by a full re-analysis. Their loadings footprint correlates 0.95–0.98 with B's own stage-change pattern, so a population-defined plant and a basis-defined plant nearly coincide.

## Goals / Non-Goals

**Goals:**
- A construction whose planted difference lies inside the stage span by definition. Its intended geometry is exact and pure in the null-standardized joint coordinates, and its dataset is the null dataset plus the plant.
- A per-replicate measurement of how much of any construction's population difference the fitted basis sees, in the same units as the plasmode's `latent_capture.py`.
- A pilot that pairs both constructions on the same baseline programs and noise, and reproduces the design-guidance pilot where they overlap.

**Non-Goals:**
- Biological realism for the stage-span construction. It is a measurement-defined sensitivity construction; the surgeries stay the default and remain the biologically grounded constructions.
- Changing what the production measurement fits, how it selects the rank, or RRPP.
- A construction axis inside a design grid, the paper-grade decision and the real-data plasmode (proposal Non-goals).

## Decisions

**D1. The plant basis is the population stage span, not a fitted basis.** The edit is defined on B's exact null population stage means, inside the span of B's own stage configuration. This was decided with the user.
- *Alternative: a reference basis fitted once per design on a large independent draw.* It would include the PLS basis's noise-covariance tilt. But it couples the generator to a PLS fit and a declared rank, while the CV rank varies by replicate, and it makes a fitted object part of the truth. The prototype's 0.95–0.98 footprint correlation says it would buy little.
- *Alternative: each replicate's own basis.* That is circular: the plant aligns with that replicate's noise, so power reads optimistic. It is kept only as the descriptive retrofit check (D12).

**D2. Edits are made in pooled-standardized joint coordinates over declared span layers.**
- The rotation and the bend must be pure in the space the measurement standardizes and concatenates. Per-block edits would be joint-pure only for scaling, the lesson of the `magnitude_kind='all'` diagnostic.
- The scale is the null population's pooled SD, combining the within-cell variance (diagonal of `cov_M`, `cov_expr` and `cov_protein`) with the between-cell variance of the null means at the realized cell proportions. That is the population analogue of the sample scaler, and it does not depend on the plant.
- `span_layers` is a generator setting because the span depends on which blocks are measured. The study loader and the harness enforce `span_layers == layers` (D9) so that a plant is never measured in a space where it is no longer pure.

**D3. The plant enters as a per-cell mean shift before sampling.** `generate_omics` and `omic_population_means` gain optional per-cell shift arrays, and an absent shift leaves the existing path untouched. Methylation is shifted in M-values before `rev_logit`. The construction draws no random numbers, so the stream is the `none` stream: group A, the baseline indicators and B's noise are identical, and B's samples are the null samples plus the plant. This makes the stage-span versus surgery comparisons paired at the noise level for `magnitude`, and paired at the baseline program for every mode.

**D4. Mode definitions keep the off-target statistics exactly at zero.**
- *Magnitude:* scale about the centroid by 1 + e, the same factor as the `joint` surgery, so the two magnitude constructions are directly comparable. They differ only by a translation and by the cascade's per-feature pattern.
- *Orientation:* rotate rigidly by θ = 90°·e in the plane spanned by u (the net displacement) and v (the next principal direction inside the configuration's span). Using u makes the production net-displacement `angle` exactly θ. Using v keeps the rotated path inside the span. e = 1 gives 90°, the surgery's population angle at full relocation. Intermediate effects map to different angles under the two constructions, so the reading compares by effect and reports both constructions' intended population angles from truth.
  - *Alternative: rotate toward an unused direction.* That is what the surgery effectively does, and it is the case this change contrasts against.
- *Shape:* bend the middle interior stage (the stage the surgery bends) by e mean steps, perpendicular to its chord inside the span, then rescale to the original path length and recentre. Path length and net-displacement direction are preserved exactly, so `delta` and `angle` stay at zero. Under PC1 the bend can move the principal axis, but the pilot measures with the production estimator.
  - *Alternative: bend without rescaling.* It co-moves `delta`, as the surgery does.
- *Translation:* offset every stage by e mean steps along u. This is the in-span offset most confounded with progression ("B is ahead at every stage"), and so the hardest negative control. The surgery's translation stays outside the program.

**D5. Degenerate spans fail loudly.** `orientation` and `shape` need a second in-span direction, so fewer than three stages, or a numerically one-dimensional configuration, is an error rather than a silent edit outside the span. This matches the repository's censoring policy of refusing a construction that cannot be realized as declared.

**D6. Identity and hashing.**
- `construction` and `span_layers` (generator) and `span_capture` (evaluation) join `_SIGNATURE_OMIT_WHEN_NONE`.
- `"surgery"` and `false` normalize to `None`, the hash-omitted spelling, so every committed configuration and fixture keeps its cell ids.
- The anchor builder clears `construction` and `span_layers`, so sibling anchors share one identity and one dataset.

**D7. Constructions are compared as sibling profiles, not as a design axis.** A construction axis would give two columns whose anchors are the same dataset. Enumeration would then need to deduplicate anchors across design points, which touches the design-point readers and their byte-identity guarantees. Sibling profiles in one seed family already have that semantics (the design-guidance magnitude and geometry siblings), and `scripts/anchor_reproduction.py` already checks it. The loader rejects `generator.construction` as an axis and points to siblings.

**D8. Span capture is measured in the fitted model's input coordinates.**
- `population_latent` standardizes the population means with the replicate's fitted preprocessor (`transform_population`), concatenates the measured blocks and applies `model.transform`. It is the truth seen through the estimated basis, free of sampling noise in the cell means.
- `span_capture` maps the same population means to the model's input coordinates ((x − x_mean)/x_std), takes an orthonormal basis Q of the rotations' span (QR), and reports the shares of Σ_s‖d_s‖² and of the centred average stage means that lie in span(Q). d_s is B − A with its stage mean removed, because a translation is invisible to all three statistics.
- This is the plasmode's `latent_capture.py` definition, so the shares are comparable across the two studies.
- It is opt-in and computed from the evaluation's own fit (no refit), and statistics are unchanged by construction.

**D9. Guarding the measured blocks.** The loader checks every column's resolved `span_layers` against its resolved `layers` (absent means all three). The harness checks again, from the dataset's truth, before fitting, because datasets can reach the harness outside a study.

**D10. The pilot runs in the design-guidance family.**
- Base seed 1028, 100 × 199, the core profiles' generator (p_dmp 0.1, ρ 0, δ 2.0, `joint`).
- The n = 200 and n = 80 columns at three stages and two blocks are cells the 2026-10-28 pilot already ran. The surgery siblings therefore reproduce those records exactly, which is a free regression check that span capture changes no statistic.
- The SEA-AD table column adds the case-study design and unequal cells, the "kite" route by which surgeries become partly visible.
- The columns are explicit points over `generator.n_samples` and `generator.group_stage_sizes`. The table point keeps the default `n_samples`, as the size-table conflict rule requires.
- Effects follow the slice: magnitude 0.02 / 0.05 / 0.1, geometry 0.25 / 0.5 / 1.0.

**D11. The reading's thresholds.**
- *R1:* the stage-span translation holds the k = 2 bound (0.0936 at R = 100). This is the new construction's only negative control; the anchors are the surgery pilot's datasets, whose size has already been read.
- *R2:* median capture of at least 0.90. The plant lies in the population span by definition, so lower capture means the fitted basis misses that span at the design, and the construction does not deliver its premise.
- *R3:* off-target rates at or below the bound. Purity is exact in the population, so off-target rejection beyond chance would mean the sample basis mixes the statistics.
- Target power is reported, not judged: it is n-conditional, and the paper-grade follow-up sets its decision from this pilot.

**D12. The retrofit check is the prototype, promoted.**
- It runs at the n = 200 column, replicates 0–9, with all four edits, both preimages, frozen and refit re-measurement, and 199 permutations.
- It adds the cosine between each footprint and the population plant for the same mode on the same replicate's B configuration. Both are in pooled-standardized feature units.
- It is descriptive. It documents the equivalence D1 relies on and the circular-basis alternative it rejects.

## Risks / Trade-offs

- [The construction is circular by design: it plants where the measurement looks.] → It answers a different question from the surgeries: test sensitivity given visibility. The report and the visibility page frame it that way. R2 is a premise check, not a finding.
- [Stage-span bypasses the methylation → expression → protein cascade.] → Truth names the construction, and the surgeries remain the default. The construction is documented as measurement-defined, not biological.
- [The plant adds between-cell variance, so the sample scaler differs from the null population scale and realized geometry drifts from the intended geometry.] → `population_latent` reports the realized geometry beside the intended one, and the report quotes both. The plant is at most about 1 SD per feature at the pilot effects.
- [Methylation clipping at extreme M-values distorts the plant on a few CpGs.] → Clipping is shared with every construction. Truth records the plant's norm per stage; the reading needs no exactness beyond R2/R3 tolerances.
- [The surgery siblings fail to reproduce the 2026-10-28 records because of environment drift.] → Run from a clean worktree, record versions in PROVENANCE, and report mismatches by column. A mismatch is diagnosed before the reading is written.
- [The SEA-AD table at the design-guidance stage effect is not the case-study configuration.] → The column is labelled as the table under the design-guidance generator. Case-study-configuration claims stay with the existing case-study checks.
- [Wide arrays stall on `ing`'s NFS.] → Throttle the array as earlier runs did (memory: cold-start limit per node).

## Migration Plan

All settings are opt-in and hash-omitted while unset, so no committed artifact changes. The order is:

1. Ship code and tests, plus fixture checks for every committed config.
2. Commit the four profiles and register them in `POST_FLOOR_CONFIGS`.
3. Run a local smoke.
4. Run on `ing`.
5. Run the reproduction checks and the reading.
6. Write the report.

Rollback means reverting the change. Datasets and records without the new fields are unaffected.

## Open Questions

- Whether `motco simulate` should expose `--construction`. The CLI is unchanged here, and adding it later changes no spec in this change.
- The paper-grade follow-up's columns (four stages, three blocks) and effect grid are taken from the pilot's reading.
- The latent plasmode's donor split for its planting basis is designed in its own change.
