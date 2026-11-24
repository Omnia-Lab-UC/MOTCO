# Univariate interaction baseline — report (2026-11-23)

**Change:** `univariate-interaction-baseline` (MOTCO 0.25.0, no package change). **Reads:**
`results/univariate-baseline-2026-11-23/` (`univariate_baseline_decision.json`; per profile `settings.json`,
`rate_table.csv`, `paired_table.csv`, `cell_table.csv`, `accounting.json`; `parent/` and `surgery/reproduction.json`;
`PROVENANCE.txt`, `sacct.txt`). **Parents:** the [interaction-scan validation](interaction-scan-2026-11-09.md) (its
simulation settings and comparator table, bound by sha256) and the
[design-guidance slice](design-guidance-slice-2026-10-29.md) (its two core configs and operating tables). **Concept
page:** [measurement visibility](../measurement-visibility.md#against-the-default-analysis-the-univariate-baseline).

The first comparison a reader makes is against the default analysis with no integration: a group × stage interaction
test on each feature. The package has that test as the interaction scan, which was validated as a complement for
sparse differences and never judged as a baseline (roadmap open item 13). This change judges it as one. On each unit,
one production fit feeds the trajectory test, which enters as one omnibus test (smallest of `delta`, `angle` and
`shape` at 0.05 / 3). Its standardized features feed the per-feature interaction F, combined three ways: the scan's
max-T maximum (`scan`), Benjamini–Hochberg on the F p-values (`fdr`), and the split test's `total` part (`total`).
One decision, fixed before any unit ran, labels each of 71 signal cells by paired exact tests.

**In one paragraph.** **The verdict is `complementary`.**

- **Each side leads somewhere, and where follows how the difference spreads.** Production leads in 16 cells: every
  one is a small stretch of the stage program or a small bend inside its span. For example, a stretch of e = 0.05 at n = 200 reads
  0.815 against `fdr`'s 0.235 and `total`'s 0.045. The baseline leads in 21 cells: the feature-set surgeries'
  orientations and shapes at n = 80 and 200 (`fdr` 1.000 against production's 0.26–0.86) and one shape at n = 500,
  the out-of-span rotations, and in-span rotations at n = 80 and at the weaker stage effect. The other 34 cells are unresolved, most because both
  sides saturate.
- **The default analysis holds size in two of its three forms.** `fdr` and `total` hold in both families. The max-T
  `scan` is inadmissible: on the three-block design at n = 80 it rejects 38 and 35 of 500 null units, above the
  per-cell bound of 0.0695, where its rule allows one such cell. Production holds everywhere (at most 26 of 500).
- **The grade is `detects_only`.** Production detects 26 of the 36 in-span stage-span cells at 0.80 or more. Its
  off-target statistics fail as a family on one cell: a 45° rotation at the weaker stage effect, where `shape` reads
  0.665 and `delta` 0.155. A rejection there names a difference, not its kind.
- **Both reproductions pass exactly** (196 and 99 rows), and 33,700 of 33,700 units completed.
- **The consequence fixed with the rule.** The documentation states the map of which side leads by plant family and
  design. A case study reports both analyses under a multiplicity policy fixed in advance. The reading holds at 658
  simulated features with Gaussian noise; no real-noise reading exists (roadmap open item 2). No verdict changes what
  `motco analyze` runs.

## 1. Configuration and provenance

- **Profile `parent`: the question's own units, reproduced.** The scan validation's simulation units for the `n200`
  and `seaad` designs: two blocks (498 features), base seed 2035, the anchor and translation at 500 replicates and
  twelve stage-span cells at 100. 4,400 units. It gates the run and describes; it decides nothing, because its units
  are the ones the question was read from.
- **Profile `stage_span`: a fresh family on three blocks.** Methylation, expression and proteomics (658 features),
  balanced n = 80, 200 and 500 at stage effect δ = 2.0, and n = 200 at δ = 0.5, base seed 2043 (unused by any driver).
  The anchor and a translation at 500 replicates. The in-span ladder under the planted scale, each cell at 200:
  stretches at e = 0.05, 0.1 and 0.25, rotations at φ = 0 and bends, each at e = 0.1, 0.25 and 0.5. The out-of-span
  rotation at e = 0.5, φ = 90 and its 20-feature form under the null scale. 12,800 units, 44 signal cells. The cells
  of a replicate share their null dataset.
- **Profile `surgery`: the slice's three-block units.** The design-guidance slice's columns at n = 80, 200 and 500
  with three blocks, from its committed configs: the anchor, translation at e = 1.0, magnitude at e = 0.02, 0.05 and
  0.1, and orientation and shape at e = 0.25, 0.5 and 1, every cell at 500. 16,500 units, 27 signal cells. Every unit's
  permutations are drawn at seed 3, as the slice's were (its configs set `evaluation.seed: 3`).
- **Arms on one fit.** `production` (`pls_cv`, rank floor 2, net displacement, 999 permutations) and its three single
  statistics. `scan`, `fdr` and `total` on the fit's standardized features under the trajectory test's full and
  reduced models, at the unit's RRPP seed, uncorrected. Recorded and never judged: the split test's `span` and
  `complement`, `production_minp` (Westfall–Young min-P over the three statistics' joint null), and `pca_matched`, the
  trajectory statistics on PCA at the unit's CV rank.
- **`concat` was not run.** The local smoke measured it at about 80 % of every unit's time: its `shape` aligns the stage
  means in the 658-dimensional ambient space on each draw. With it the run projected 837–1,673 core-hours, against
  156–311 without. At the approval gate it was dropped, and the cap was set at 350 recorded core-hours
  (`openspec/changes/univariate-interaction-baseline/notes.md`). The
  [latent-space comparison](latent-space-comparison-2026-11-04.md) measured `concat` against production on the
  stage-span construction across alignment; this report cites it for that arm.
- **Code and host.** Revision `a8a27dc`, pushed before any unit ran, from a clean detached worktree on `ing`. Python
  3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3, BLAS pinned to one thread.
- **Where it ran.** Arrays 957651 (`parent`, 40 tasks), 957652 (`stage_span`, 100) and 957653 (`surgery`, 120),
  partition `512x1024` with n11 excluded: 260 of 260 `COMPLETED 0:0`, 18:45–21:09. **236.18 recorded core-hours**
  (14.01 / 88.49 / 133.68) against the cap of 350. The median unit took 10.9 s, 23.2 s and 25.0 s by profile, most of
  it the cross-validated fit.

## 2. Accounting and the reproductions

- **Units.** 4,400, 12,800 and 16,500 complete: no failed unit, no duplicate, no failed planted fixed point.
- **`parent` reproduction: passed.** All 196 rows of the scan validation's comparator table are equal: the scan, the
  split test's parts and `delta` / `angle` / `shape`, per design and cell.
- **`surgery` reproduction: passed.** All 99 rows of the slice's operating tables are equal: `delta` / `angle` /
  `shape` per column, mode and effect, none unresolved.
- **The decision regenerates.** `decide` on the committed aggregates alone, on the workstation, writes a decision
  byte-identical to the cluster's.

## 3. The decision

The order is gate, size, detection. Size is judged per arm in two families: `stage_span` (m = 8) and `surgery`
(m = 6). Each family uses the count-plus-ceiling rule (`judge_rates`), with a per-cell bound of 0.0695 at 500
replicates. Then each signal cell is labelled by paired exact tests of production against every admissible baseline
arm:

- **`production_leads`.** The gain p-value is at most 0.025 / 71 = 3.5 × 10⁻⁴ against every admissible arm, and the
  rate exceeds each arm's by at least 0.10.
- **`baseline_leads`.** The loss p-value is at most 0.025 / (71 · A) against at least one arm, and that arm's rate
  exceeds production's by at least 0.10. A is the number of admissible arms; with A = 2 the level is
  1.76 × 10⁻⁴.
- **`unresolved`.** Otherwise.

The verdict is by existence: `complementary` when both labels occur. The out-of-span cells stay in the family. This
makes `complementary` and `baseline_leads` the live outcomes (design D6).

| criterion | reading |
|---|---|
| gate | both reproductions passed; 0 failed units |
| size | `production`, `fdr` and `total` admissible; `scan` inadmissible (`stage_span`) |
| labels | 16 `production_leads`, 21 `baseline_leads`, 34 `unresolved` of 71 |
| verdict | **`complementary`** |
| grade | **`detects_only`** |

## 4. Size per arm

Rejections of 500 on the null cells (anchor; translation, a constant group offset with no interaction):

| profile | design | production | `scan` | `fdr` | `total` |
|---|---|---|---|---|---|
| `stage_span` | n = 80 | 25; 26 | **38; 35** | 31; 31 | 26; 28 |
| `stage_span` | n = 200 | 24; 17 | 25; 26 | 18; 18 | 20; 17 |
| `stage_span` | n = 500 | 15; 21 | 28; 29 | 23; 23 | 17; 18 |
| `stage_span` | n = 200, δ = 0.5 | 13; 15 | 23; 22 | 20; 20 | 30; 30 |
| `surgery` | n = 80 | 15; 17 | 19; 25 | 16; 21 | 20; 24 |
| `surgery` | n = 200 | 24; 18 | 22; 17 | 16; 14 | 22; 22 |
| `surgery` | n = 500 | 22; 17 | 34; 22 | 27; 17 | 25; 22 |

- **`scan` is inadmissible.** 38 / 500 = 0.076 and 35 / 500 = 0.070 both exceed the bound of 0.0695 in the
  `stage_span` family, which allows one exceedance. Its `surgery` family holds.
- **The two exceedances are nearly one observation.** In `stage_span`, the anchor and the translation of a replicate
  are the same null dataset; the translation adds a group offset, which the interaction F does not see. Only the
  permutation seeds differ by cell. `fdr`, which draws no permutation, reads 31 and 31 there. The rule counts cells as
  declared, and the arm leaves the comparison; the report does not waive it.
- **What it means.** At 80 samples and 658 features, the max-T scan under Freedman–Lane draws runs above its level
  (the scan validation read 0.056 at n = 200 on two blocks). The parametric `fdr` on the same F values holds (at most
  0.062), as does `total` (at most 0.060). Production's largest null rate is 0.052.

## 5. The map

Cells per label (production leads / baseline leads / unresolved) by plant family and design:

| family | n = 80 | n = 200 | n = 500 | n = 200, δ = 0.5 |
|---|---|---|---|---|
| stage-span stretch | 2 / 0 / 1 | 1 / 0 / 2 | 1 / 0 / 2 | 2 / 0 / 1 |
| stage-span rotation (in span) | 0 / 2 / 1 | 0 / 0 / 3 | 0 / 0 / 3 | 0 / 1 / 2 |
| stage-span bend | 1 / 0 / 2 | 1 / 0 / 2 | 0 / 0 / 3 | 2 / 0 / 1 |
| stage-span rotation out of span | 0 / 2 / 0 | 0 / 1 / 1 | 0 / 0 / 2 | 0 / 2 / 0 |
| surgery magnitude | 2 / 0 / 1 | 2 / 0 / 1 | 2 / 0 / 1 | — |
| surgery orientation | 0 / 3 / 0 | 0 / 3 / 0 | 0 / 0 / 3 | — |
| surgery shape | 0 / 3 / 0 | 0 / 3 / 0 | 0 / 1 / 2 | — |

### 5.1 Where production leads (16 cells)

Rates of production, `fdr` and `total`:

| profile | design | cell | production | `fdr` | `total` |
|---|---|---|---|---|---|
| `stage_span` | n = 80 | stretch e = 0.05 | 0.345 | 0.060 | 0.050 |
| `stage_span` | n = 80 | stretch e = 0.1 | 0.750 | 0.310 | 0.045 |
| `stage_span` | n = 80 | bend e = 0.1 | 0.365 | 0.105 | 0.045 |
| `stage_span` | n = 200 | stretch e = 0.05 | 0.815 | 0.235 | 0.045 |
| `stage_span` | n = 200 | bend e = 0.1 | 0.860 | 0.420 | 0.055 |
| `stage_span` | n = 500 | stretch e = 0.05 | 0.985 | 0.710 | 0.035 |
| `stage_span` | n = 200, δ = 0.5 | stretch e = 0.1 | 0.285 | 0.060 | 0.065 |
| `stage_span` | n = 200, δ = 0.5 | stretch e = 0.25 | 0.810 | 0.345 | 0.090 |
| `stage_span` | n = 200, δ = 0.5 | bend e = 0.25 | 0.325 | 0.065 | 0.070 |
| `stage_span` | n = 200, δ = 0.5 | bend e = 0.5 | 0.725 | 0.390 | 0.085 |
| `surgery` | n = 80 | magnitude e = 0.05 | 0.382 | 0.050 | 0.044 |
| `surgery` | n = 80 | magnitude e = 0.1 | 0.752 | 0.314 | 0.044 |
| `surgery` | n = 200 | magnitude e = 0.02 | 0.316 | 0.040 | 0.044 |
| `surgery` | n = 200 | magnitude e = 0.05 | 0.822 | 0.146 | 0.050 |
| `surgery` | n = 500 | magnitude e = 0.02 | 0.770 | 0.074 | 0.052 |
| `surgery` | n = 500 | magnitude e = 0.05 | 0.990 | 0.730 | 0.066 |

A stretch or a small bend moves every feature of the stage program a little, in the direction the fit already
measures. Each feature's F is weak, and the trajectory statistic pools them. At a larger effect, each feature's F
crosses its own threshold and the cell saturates for both sides (for example, a stretch at e = 0.1 at n = 200:
production 0.985, `fdr` 0.965). The split test's `total` part, which pools the features without the stage fit, reads
0.035–0.090 in all 16 cells.

### 5.2 Where the baseline leads (21 cells)

| profile | design | cell | production | `fdr` | `total` |
|---|---|---|---|---|---|
| `stage_span` | n = 80 | rotation e = 0.1, φ = 0 | 0.820 | 0.970 | 0.055 |
| `stage_span` | n = 80 | rotation e = 0.25, φ = 0 | 0.855 | 1.000 | 0.270 |
| `stage_span` | n = 80 | rotation out of span | 0.255 | 1.000 | 0.860 |
| `stage_span` | n = 80 | rotation out of span, 20 features | 0.835 | 1.000 | 0.245 |
| `stage_span` | n = 200 | rotation out of span | 0.720 | 1.000 | 1.000 |
| `stage_span` | n = 200, δ = 0.5 | rotation e = 0.25, φ = 0 | 0.410 | 0.995 | 0.125 |
| `stage_span` | n = 200, δ = 0.5 | rotation out of span | 0.035 | 0.840 | 0.280 |
| `stage_span` | n = 200, δ = 0.5 | rotation out of span, 20 features | 0.115 | 1.000 | 0.180 |
| `surgery` | n = 80 | orientation e = 0.25 / 0.5 / 1 | 0.324 / 0.482 / 0.482 | 1.000 | 1.000 |
| `surgery` | n = 80 | shape e = 0.25 / 0.5 / 1 | 0.256 / 0.374 / 0.520 | 1.000 | 0.742 / 0.998 / 1.000 |
| `surgery` | n = 200 | orientation e = 0.25 / 0.5 / 1 | 0.742 / 0.856 / 0.852 | 1.000 | 1.000 |
| `surgery` | n = 200 | shape e = 0.25 / 0.5 / 1 | 0.524 / 0.720 / 0.864 | 1.000 | 1.000 |
| `surgery` | n = 500 | shape e = 0.25 | 0.802 | 1.000 | 1.000 |

- **The surgeries' orientations and shapes are outside the span.** A per-omic feature permutation moves a few features
  far. The fitted span holds 0.4–1.2 % of the group difference (medians), and the split test's `complement` rejects in
  every unit of 17 of the 18 cells (0.774 in the 18th, a shape at e = 0.25 at n = 80). Production reaches 0.26–0.86 only through incidental routes. At n = 500 it reaches 0.95–0.99 on
  orientation, and those cells are unresolved.
- **The out-of-span rotations are where the scan was validated.** Production's 0.035–0.835 in the five cells the
  baseline leads is the blindness the stage-span construction was built to show. At n = 500 production reaches 0.975
  and 1.000, and both cells are unresolved.
- **In-span rotations lead for the baseline only at small samples or a weak stage effect.** At n = 80 a rotation at
  e = 0.1 reads 0.82 for production against 0.97, from 36 units only the baseline rejects and 6 only production
  rejects. At δ = 0.5 a rotation at e = 0.25 reads 0.41 against 0.995. From n = 200 at δ = 2.0 every in-span rotation
  is unresolved: both read 0.96 or more.

### 5.3 Unresolved (34 cells)

Thirty of the 34 have both production and `fdr` at 0.88 or more, the saturated end of each ladder. The other four
sit near the level on both sides: at δ = 0.5 a stretch at e = 0.05 (production 0.075, `fdr` 0.045), a rotation at
e = 0.1 (0.090, 0.095) and a bend at e = 0.1 (0.075, 0.035), and the surgery's magnitude at e = 0.02 at n = 80
(0.070, 0.034). Every rate is in the rate table and the decision's
`cell_rates`.

## 6. The characterization grade: `detects_only`

The grade reads production alone, on the 36 in-span `stage_span` cells, where off target has an exact truth under the
planted scale. Production detects 26 of them at 0.80 or more. Their 52 off-target statistics (two per cell) form one
family under `judge_rates`: a per-cell bound of 0.081 at 200 replicates, a ceiling of 0.105, and three exceedances
allowed.

| design | cell | statistic | rate | reading |
|---|---|---|---|---|
| n = 200 | bend e = 0.1 | `delta` | 0.085 | above the bound |
| n = 200, δ = 0.5 | rotation e = 0.5, φ = 0 | `delta` | 0.155 | above the ceiling |
| n = 200, δ = 0.5 | rotation e = 0.5, φ = 0 | `shape` | 0.665 | above the ceiling |

All three failing readings come from two cells, and the two above the ceiling come from one. That cell is the 45°
rotation at the weaker stage effect, where the CV rank inflates (median 4, 10–90 % range 3–5, against 2 at δ = 2.0).
This is the rank confound the slice found below δ = 2.0, now as off-target leakage. The other 49 readings hold.
Under `detects_only`, a rejection by production names a group × stage difference; which statistic rejects does not
name its kind at every design. The baseline cannot make that second claim at all (design D7).

The surgeries' off-target rates are reported in the decision (`off_target_reported`, 35 readings) and not graded,
because those plants are not geometrically pure.

## 7. Descriptive readings

### 7.1 The recorded arms

- **`production_minp`.** Its rate is never below the omnibus's and at most 0.055 above (n = 80, stretch e = 0.05:
  0.400 against 0.345), in all 85 cells. The Bonferroni omnibus costs little here.
- **`pca_matched`.** The trajectory statistics on PCA at production's rank rarely reach production's power: in one of
  71 signal cells only (the surgery's orientation at e = 1 at n = 80, 0.840 against 0.482). Unsupervised components at
  the same dimension do not hold the stage program that production's fit holds.
- **`span` and `complement`.** `span` tracks production on the in-span cells. `complement` reads at its level on the
  stretches and bends (0.035–0.075) and on rotations up to e = 0.25 (0.035–0.110), and rises to 0.16–0.31 on the
  rotation at e = 0.5. It reads 1.000 on every surgery orientation and shape but one (0.774), and 0.885–1.000 on the
  out-of-span rotations except the 20-feature form at n = 80 (0.245). That is the split test's statement: these
  differences lie outside the fitted span.

### 7.2 Discoveries

On the null cells the share of units with any `fdr` discovery is the arm's rejection rate (0.028–0.062). On the
signal cells, the `fdr` arm's realized false discovery proportion on population-null features averages 0.005–0.045
per cell, below the 0.05 it controls. The `scan` arm controls the family-wise rate instead, and finds fewer features:
for example, a median of 192.5 hits against 223 discoveries on the surgery's orientation at e = 0.5 at n = 80.

### 7.3 The fit

Production's CV rank has median 2 at every δ = 2.0 design (3 in a few n = 80 cells) and 3–4 at δ = 0.5. Its median
held-out AUROC is at least 0.918 in every cell. The in-span stage-span plants keep 0.88–0.99 of the group difference
in the fitted span (medians), the out-of-span rotations 0.003–0.010 (0.03–0.17 for the 20-feature form).

### 7.4 The `parent` profile

The scan validation's own units, now with `fdr` and production's omnibus beside the committed tests (decision,
`descriptive.parent`). They repeat the question's pattern. At n = 200, a stretch reads production 0.66 against `scan`
0.09 and `fdr` 0.06, and a bend at e = 0.25 reads 0.99 against 0.85 and 0.74. An in-span rotation at e = 0.25 reads
0.73 against 1.00 and 1.00, and the out-of-span rotation at e = 0.25 reads 0.01 against 0.76 and 0.70. On this profile
the omnibus's null rates are 0.024–0.036.

## 8. What follows

- **The documentation states the map** (design, migration plan). [Measurement
  visibility](../measurement-visibility.md#against-the-default-analysis-the-univariate-baseline) now says what each side
  tests and where each leads. The trajectory test leads on diffuse stretches and small bends inside the stage program.
  The per-feature test leads on concentrated differences: feature-set surgeries, out-of-span rotations, and in-span
  rotations at small n or a weak stage effect.
- **A case study reports both analyses** under a multiplicity policy fixed before any real-label reading. The
  admissible per-feature combinations at these designs are `fdr` and `total`. The max-T `scan` ran above its level at
  n = 80 on three blocks.
- **What only the trajectory test claims.** Detection of a diffuse difference that no single feature carries. Under
  `detects_only`, the statistic that rejects is not a reliable name for the kind of difference at every design: the
  weak-stage-effect rotation leaks into `shape`.
- **What does not change.** `motco analyze` runs what it ran; no default moves.
- **Follow-ups the change names.** A separate change may add the false-discovery reading to the scan's output and to
  `motco analyze`, since `fdr` was admissible and leads where the scan cannot be read. The real-noise reading of the
  same arms belongs to the first multi-omic ROSMAP or ADNI export (roadmap open item 2). A moderated per-feature
  statistic (roadmap item 6) would be a fourth baseline arm.

## 9. Limitations

- **Simulated features with Gaussian noise.** Every claim holds at 658 simulated features with one variance per
  feature. Per-feature tests are exact there, and a moderated statistic has nothing to add. Real methylation,
  expression and proteomics have heavy tails, mean–variance dependence and correlated features. No real-noise reading
  exists (roadmap open item 2).
- **The map depends on the cells chosen.** The verdict is by existence for that reason; the map, not a count, carries
  the content. Stretches and bends of the stage program favour production, and permutation surgeries and out-of-span
  rotations favour the baseline. How real differences divide between these forms is not known.
- **`concat` was not run** (approval gate, cost). Its comparison against production is the latent-space comparison's,
  on a related construction and other units.
- **The scan's inadmissibility rests on two correlated readings** of one design. The rule was declared that way and is
  applied as declared. A scan at n = 80 on more blocks, or with a corrected residual, was not read.
- **The surgery profile reuses the slice's permutation seed** (3 for every unit). That reproduces the published
  counts exactly. It also means the profile's units differ only in their generator seeds.
