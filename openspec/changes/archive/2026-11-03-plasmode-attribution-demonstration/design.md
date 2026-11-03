# Design

## Context

See `proposal.md` (Why) for the motivation. The demonstration composes three released pieces:
- **`analyze_trajectories` (MOTCO 0.14.0).** It takes declared covariates and covariate-adjusted attribution, and exposes `latent_artifacts` in memory: the fitted estimator, the standardized joint features and the original scales.
- **`sead_sexdiff`'s astrocyte export,** `results/motco_input/cohort_30`. It holds 80 donors, `rna` and `atac` blocks of 5,000 TMM log2-CPM features each, metadata `donor_id`/`sex`/`stage` with stage order `0`/`1`/`2`, and a covariate table.
- **The case-study primary's configuration** (`sead_sexdiff` `ANALYSIS_*` at `2d427b0`):
  - covariates `age_at_death`, `pmi`, `rin`, `years_education`, `log10_n_rna` and `log10_n_atac`, with `missing_covariates="error"`;
  - rank floor 2, net displacement, `residual_correction="df"`, 999 permutations;
  - the harness CV defaults (`random_state` 1203).

  Its CV-selected rank is 5 and its mean held-out stage AUROC is 0.635556 (SEA-AD report §3, §10).

Three constraints shape the approach. The exit review (§6, §9) forbids reading any SEA-AD sex × stage result. MOTCO commits aggregates only. And `sead_sexdiff` is not changed.

The closest precedent is `scripts/adjusted_attribution_check.py` (item 11): a calibrated construction, one shipped analysis per unit, descriptive arms on `latent_artifacts`, sharded resumable records, merge, and a pure decision function. This check follows its structure and conventions. It does not import its construction, which is semi-synthetic.

## Goals / Non-Goals

**Goals:**
- Run the shipped analysis, unchanged, on a real case-study substrate with a planted signal whose drivers are known exactly.
- Keep the real sex labels out of every computation.
- Fix every threshold, the null-making step and the strength rule before any unit, so the verdict cannot be tuned.

**Non-Goals:**
- Biological realism of the planted signal. This is a traceability demonstration, not a power study at a realistic effect.
- Testing confounded-covariate removal, which item 11 validated. Under the within-stage permutation the covariates are balanced across the permuted groups in expectation.
- Any library change, release or new CLI flag.

## Decisions

### D1 — Substrate: astrocyte `cohort_30`

The exit review names it: it is the primary cohort, it is open data, and its size under `df` is measured (corrected-permutation report).

Alternatives:
- `cohort_50`: its size under `df` is only inferred.
- Microglia: not the primary.
- ADNI: under a data-use agreement, and its design analogue sits at the CV cap.

### D2 — Input contract read in place; covariates derived in MOTCO

- **The contract.** The check takes `--export-dir` and reads the manifest's block and metadata entries. It requires exactly the primary's blocks, the two groups `Female`/`Male` and the stage order `0`/`1`/`2`. It derives `log10_n_rna` and `log10_n_atac` from `n_rna` and `n_atac`, as `sead_sexdiff`'s `ANALYSIS_LOG10_COLUMNS` does, and records the sha256 of every file.
- **Alternative: have `sead_sexdiff` export the derived columns.** Rejected: it needs a `sead_sexdiff` change, which is out of scope.
- **Why the derivation is safe.** The covariates do not enter the latent space, so the fidelity gate (D6) cannot check them. The derivation is instead pinned by a unit test and quoted in the report, with its `sead_sexdiff` source and commit.

### D3 — Null-making: permute the group labels within stage

Each replicate permutes the `sex` column within each stage. The rows (blocks and covariates) stay attached to their donors. This keeps the cell table `[[11, 10, 28], [9, 10, 12]]`, the stage structure and the covariates' relation to stage and to the features. Any sex × stage structure is destroyed in expectation.

A permuted vector that equals the real one is refused and redrawn, and the redraw is counted. Its probability is about 10⁻²¹, so the guard is cheap and states the rule in code.

Alternatives:
- **Remove the estimated real sex effect.** Rejected: it computes the forbidden result.
- **Donor bootstrap.** Rejected: it keeps any real sex effect.
- **Permute across stages.** Rejected: it breaks the cell sizes the configuration's size was measured at.

Consequence: size at the anchor is size **conditional on the export**. The 500 units are i.i.d. label draws given one data matrix, so the binomial SE in criterion (i) is valid conditionally. This is the question that matters for the case study's own data.

### D4 — The planted kind: a linear ramp on a known set, in one group

`X[i, f] += e · c(stage_i) · σ_f · s_f` for donors i in the second declared group (`Male`, 31 donors) and features f ∈ K.

- **The set K.** 10 features drawn per replicate from each block (a two-omic module). Its signs `s_f` are drawn uniformly.
- **The stage weight.** `c` = 0, 0.5, 1. The ramp starts at the first stage, so it adds no location offset at baseline.
- **The scale.** `σ_f` is the pooled within-stage SD of the export's feature (ddof = number of stages). It is label-free, so the plant's scale does not depend on the permuted labels.
- **Why it is an orientation plant.** A random ±1 pattern on 20 of 10,000 features is nearly orthogonal to the shared trajectory, so the ramp turns the planted group's net displacement. To first order it adds length only in quadrature.
- **What else it moves.** It lengthens the path, so `delta` can reject too. It bends a curved path slightly, so `shape` can move. Both are recorded at `e*` and are descriptive.
- **Why |K| = k = 20.** If |K| exceeded k and the effects were equal, which planted features entered the top-k would be arbitrary. Top-k frequency over K would then be capped near k/|K| even under perfect separation. With |K| = k, top-k frequency over K is the bootstrap's planted-set precision.

Alternatives:
- **The generator's orientation surgery (feature permutation).** Rejected: on real data it swaps values between features with different distributions, and its driver set is not a fixed feature set.
- **A length-preserving rotation of the planted group's stage means in span(d̂, ŵ).** Rejected: its (cos θ − 1)·d̂ term touches every feature on the shared trajectory, so the drivers are no longer exactly K.
- **Orthogonalizing ŵ against the pooled displacement on K.** Rejected: it makes the per-feature effects unequal for a negligible gain, since |cos| is already small. The realized cosine is recorded instead.

### D5 — Strength: calibrated on detection only, with attribution off

- **The ladder.** e ∈ {0.5, 1.0, 1.5, 2.0, 3.0}, with 40 units per step at 99 permutations and from a separate seed stream.
- **The adoption rule.** `e*` is the first step at which `angle` rejects in at least 0.90 of the units. The SE at 0.90 is about 0.047 at R = 40, so the margin is two SEs above criterion (ii)'s 0.80.
- **Attribution is off in calibration.** The choice of `e*` then cannot see criterion (iii), and calibration is cheaper.
- **The cap.** If no step qualifies, `e*` = 3.0 and the run proceeds. A 3-SD ramp on 20 features is already large. If even that is not detected, criterion (ii) fails, which is a finding about the design and is reported on its own terms.
- **Alternative: stop for the user's decision at the cap.** Rejected: the only options there are to extend the ladder, which moves the goalposts, or to accept the failure, which the predeclared rule already does.
- **The descriptive cell.** The step below `e*` (or `e*`/2) gives a dose-response reading at partial detection.

### D6 — Fidelity gate on label-free quantities

The PLS latent space is stage-supervised, covariate-blind and seeded by `random_state` 1203. It is the same for every label assignment. So `fidelity` runs one analysis on a permuted assignment with `permutations=0` and no attribution, and requires rank 5 and AUROC 0.635556 to six decimals. That one check confirms three things: the export, the MOTCO version's latent fit, and the configuration.

Every anchor unit's latent space must equal the gate's, which costs nothing. Planted units are exempt: the plant changes the data, and their rank is recorded.

### D7 — Arms: the shipped analysis, plus a descriptive unadjusted attribution

- **The primary arm** is the shipped analysis with covariate-adjusted attribution (B = 200, top-k 20).
- **The descriptive arm** is `analyze_orientation_attribution` on `latent_artifacts` without covariates, with the same bootstrap seed. It is cheap and shows how much adjustment moves the drivers on real covariates. It cannot lift the leg (exit review D6).
- **No oracle arm.** There is no injection to remove, and the real covariates' contribution is unknown.
- **The consistency check from item 11 is reused.** The adjusted means must project to the stage means to 1e-8 relative.

### D8 — Metrics, and reading (iii) on claims

The metrics are read on the observed principal-orientation contrast, the decomposition of the tested `angle`:
- the planted-set precision, |top-k ∩ K| / 20;
- the mean bootstrap sign stability over all of K;
- the mean top-k frequency over all of K.

**Why all of K, not K ∩ top-k.** The mean over all of K is always defined, and it is stricter: a planted feature the attribution missed drags it down.

**Why (iii) reads only the units whose `angle` rejects.** The Work bullet attributes significant orientation only. The unconditional readings are reported beside them, because conditioning on rejection selects favourable noise.

**Descriptive readings.** These do not enter the decision:
- sign agreement with the planted sign, under the contrast's group order;
- the PLS-captured and residual components;
- the per-transition blocks;
- the eigengap's degeneracy flag;
- the selected rank.

### D9 — Thresholds, fixed now

| criterion | reading | bound | rationale |
|---|---|---|---|
| (i) | anchor rate of each of `delta`, `angle` and `shape` | ≤ α + 2·√(α(1 − α)/R) | exit review §6, as stated |
| (ii) | `angle` rate at `e*` | ≥ 0.80 | the roadmap's power floor (`adni_design_check.POWER_FLOOR`) |
| (iii-a) | planted-set precision, mean − 2·SE | ≥ 0.50 | at least half of the top-20 are planted; chance is about 0.002, and the anchor's reading is reported as the reference |
| (iii-b) | sign stability over K, mean − 2·SE | ≥ 0.90 | the sign is reproduced in nine of ten resamples. Item 11's top-k read 0.93 under a planted orientation and 0.89 under noise, so the bar is set over K rather than over the arm's own top-k, where selection inflates it |
| (iii-c) | top-k frequency over K, mean − 2·SE | ≥ 0.50 | a planted feature stays in the top-20 in a majority of resamples (the stability-selection convention). Item 11 read 0.48 where the truth outnumbered k about fifteen-fold, a case that D4's |K| = k removes |

These bounds are part of the settings hash, so a decision under other bounds is a different decision.

### D10 — Seeds

- **Base seed 2034.** No committed check uses it: 2031–2033 are taken by the corrected-permutation, ADNI-design and adjusted-attribution checks.
- **Per-unit seeds.** A unit's seeds derive from (base seed, stream ∈ {`paper`, `calibration`}, replicate index, purpose ∈ {`labels`, `plant`, `rrpp`, `bootstrap`}) and never from the strength. So the cells of a replicate are matched.
- **The CV seed** stays at the primary's `random_state` 1203 (D6).
- **The RRPP seed** varies by unit. The primary's own `seed` 20261019 is not reused.

### D11 — Execution and outputs

- **Subcommands.**
  - `fidelity`, then `calibrate`;
  - `run`, sharded by `--index`/`--count`, resumable, and refusing a resume under other settings or another export;
  - `merge`, then `decide`.
- **Order of work.**
  1. Smoke: R = 20 × 99 permutations × B = 50, indicative.
  2. A one-unit timing at 999 × 200.
  3. An approval gate with the compute plan.
  4. Calibration, then the paper grade.
- **Where it runs.** If the timing extrapolates beyond a few hours on the 8-core workstation, with BLAS threads pinned to 1, the run goes to `ing`. There it uses absolute paths and excludes the nodes that hang, and the export is copied to `ing` with its sha256 verified against the settings.
- **What is committed.** `results/plasmode-<date>/` holds `settings.json`, `fidelity.json`, `calibration.json`, `plasmode_summary.csv`, `merge.json`, `plasmode_decision.json`, `PROVENANCE.txt`, and `sacct.txt` if `ing` is used. Records and per-unit metrics stay on the run machine; PROVENANCE names them with their sha256.

### D12 — Ownership split, confirmed

MOTCO owns the null-making step, the plant, the rule and the report. `sead_sexdiff` supplies the export as it stands and changes nothing. No donor-level data, label vector or real-label result enters MOTCO's repository.

## Risks / Trade-offs

- **[The plant is stylized; a 20-feature ramp at several SD is not a typical biological effect]** → The report states `e*` in SD units and frames the result as traceability of a known signal on real noise, covariates and design, not as power at a realistic effect.
- **[Reading (iii) on rejecting units is optimistic]** → The unconditional readings are reported beside it. With `e*` at about 0.9 calibration power, the conditioning drops few units.
- **[The shared trajectory is weak (AUROC 0.64), so the groups' principal axes may be unstable under resampling and (iii-b) may fail even when K is recovered]** → That is a property of attribution at this design, and the verdict reports it as such. The per-transition blocks and the eigengap flag are recorded to diagnose it.
- **[The plant lengthens the path, so `delta` may detect it more often than `angle`]** → (ii) reads `angle` only. `delta` and `shape` are reported.
- **[Planting changes the latent space (rank, captured directions)]** → The rank and AUROC are recorded per planted unit. Only anchor units are held to the gate.
- **[Conditional size is not the size over new datasets]** → This is intended (D3). The report says so, and cites the case-study Type I and corrected-permutation reports for the unconditional size.
- **[The local export differs from the one the re-issue analysed]** → The fidelity gate fails before any unit, and the hashes are recorded.
- **[The compute is unknown, since attribution on 10,000 features at B = 200 has not been timed]** → The one-unit timing comes before the approval gate.

## Migration Plan

Nothing is deployed and no release is made. The change adds a script, tests, results and docs. Rollback means deleting them. The roadmap's item 12 and the Phase 6 Work bullet are updated only after the decision.
