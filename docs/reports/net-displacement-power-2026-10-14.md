# Net-displacement power comparison — 2026-10-14

> Paper-grade power comparison of the net-displacement orientation estimator against PC1 in the production
> measurement space (`openspec/changes/net-displacement-power-comparison/`). It follows the
> [net-displacement Type I study](net-displacement-type-i-2026-10-13.md), whose verdict
> (`paper_grade_power_comparison`) selected it. It changes no production default; it names the change that would.
> Every rate below is read from `results/net-displacement-power-2026-10-14/{phase6,phase5}/report/`
> (`column_decision.csv`/`.json`, `design_point_operating.csv`, `eigengap_stratified_power.csv`,
> `acceptance_report.csv`), `…/{phase6,phase5}/consistency_checks.json` and
> `results/net-displacement-power-2026-10-14/comparison_decision.json`, except where a number is marked as computed
> from `merged.jsonl` (§9 gives the snippet). Rates are rejection rates at p < 0.05 over 500 replicates × 999
> permutations (Monte Carlo SE about 0.0097 at 0.05).

**Runs:** `results/net-displacement-power-2026-10-14/phase6/` (SLURM array 920303) and `…/phase5/` (array 920304) on
`ing` · **Code:** `db5ab69` (run and reports), `1ccded5` (consistency checks) · **Column decisions:**
`candidate_qualifies` at both designs · **Comparison:** `adopt_net_displacement`.

**In one paragraph.** On identical datasets and permutation draws, replacing PC1 with the net-displacement
orientation estimator in the cross-validated stage-PLS space raises orientation/`angle` power at the construction
maximum from 0.300 to 0.648 at the Phase 6 small-n design and from 0.840 to 0.924 at Phase 5. It also clears every
mandatory control at the one-sided bound 0.0695, which PC1 does not: PC1's `angle` rejects the shape construction at
0.138 (Phase 6) and 0.796 (Phase 5); net displacement rejects it at 0.054 and 0.022. Anchor size holds (0.050 at both
designs). `delta` and `shape` are identical record for record between the two columns, as they must be. The
predeclared rule therefore selects **`adopt_net_displacement`**, and the follow-up is the production-default flip.
Two qualifications carry into the case study (§7): at n = 80 the gain is concentrated at the construction maximum
and, within it, in replicates where CV retains more than three components (0.996 there, 0.280 otherwise), so a
non-rejection stays uninformative; and the Phase 6 translation control sits one rejection under the bound (0.068).

## 1. Design

- **Columns:** two explicit design points (`design_grid.points`) on matched datasets:

  | column | `integration_method` | `orientation_estimator` | role |
  |---|---|---|---|
  | `pls_cv × pc1` | `pls` (CV rank) | `pc1` (`null`) | baseline |
  | `pls_cv × nd` | `pls` (CV rank) | `net_displacement` | candidate |

  No `concat`, `pca`, forced-rank or `pls_within_group` column; the probe's ceiling numbers are cited, not re-run.
- **Designs:**
  - **Phase 6** (`examples/trajectory_power_study/net_displacement_power_phase6.json`): three stages,
    `group_stage_sizes = [[11, 10, 28], [9, 10, 12]]`, methylation + expression, base seed 800, family
    `phase6-small-n`.
  - **Phase 5** (`net_displacement_power_phase5.json`): n = 1200, four stages, `group_ratio` 0.5, all three blocks,
    base seed 600, family `phase5-primary`.
  - Both use `p_dmp = 0.1`, ρ = 0, joint magnitude and `surgery_censoring = error`.
- **Cells (36 per design, 18,000 units):** per column, the shared zero-effect anchor and four modes (`magnitude`,
  `orientation`, `shape`, `translation`) at e ∈ {0.25, 0.5, 0.75, 1.0}, plus the engine's two `type_i_baseline`
  cells (independently seeded, measured in the baseline column). That is 1,000 Type I units and 8,500 units per
  column.
- **Seeds:** both columns evaluate the same dataset at each replicate index. `evaluation.seed` is the constant 3 and
  is also the RRPP seed, so both columns draw the same permutation indices and the baseline column re-computes the
  parent studies' records (§6).
- **Rule:** `acceptance.column_decision` (rule version 1), candidate `pls × net_displacement` against the baseline:
  - **target** orientation/`angle` at the top effect (1.0): Phase 6 criterion `gain` (candidate − baseline >
    2 pooled SE); Phase 5 criterion `no_loss` (baseline − candidate ≤ 2 pooled SE), because PC1 is already near its
    ceiling there;
  - **anchor** `delta`/`angle`/`shape` ≤ α + 2·√(α(1 − α)/R) = 0.0695 at R = 500;
  - **mandatory controls** translation/`angle`, magnitude/`angle`, shape/`angle` at the top effect, each ≤ 0.0695;
    judged on the candidate only;
  - **protected** magnitude/`delta`, shape/`shape`: no loss beyond 2 pooled SE.

  The cross-design selection (`scripts/net_displacement_power_follow_up.py`) returns `adopt_net_displacement` iff
  the candidate qualifies at both designs. Both rules are advisory: they change no default.

## 2. Orientation power: `angle` by mode × effect × estimator

From `design_point_operating.csv` (`statistic = angle`, `power_primary` = PC1, `power_design` = net displacement).
Rate (MC SE). The probe's `concat × net_displacement` ceiling reference reached 1.00 at orientation e = 0.5 and 1.0
at both designs (100 replicates × 199 permutations, its own seeds; [probe report](orientation-measurement-probe-2026-10-12.md)
§3). It is descriptive context only, never a candidate.

**Phase 6 (n = 80, K = 3).**

| mode | e | PC1 | net displacement |
|---|---|---|---|
| anchor | 0 | 0.024 (0.007) | 0.050 (0.010) |
| orientation | 0.25 | 0.034 (0.008) | 0.084 (0.012) |
| orientation | 0.50 | 0.058 (0.010) | 0.126 (0.015) |
| orientation | 0.75 | 0.090 (0.013) | 0.178 (0.017) |
| orientation | 1.00 | 0.300 (0.020) | **0.648** (0.021) |
| shape | 0.25 / 0.50 / 0.75 / 1.00 | 0.042 / 0.076 / 0.112 / 0.138 | 0.056 / 0.068 / 0.042 / 0.054 |
| magnitude | 0.25 / 0.50 / 0.75 / 1.00 | 0.016 / 0.010 / 0.014 / 0.008 | 0.036 / 0.024 / 0.024 / 0.030 |
| translation | 0.25 / 0.50 / 0.75 / 1.00 | 0.018 / 0.018 / 0.030 / 0.028 | 0.068 / 0.058 / 0.066 / 0.068 |

**Phase 5 (n = 1200, K = 4).**

| mode | e | PC1 | net displacement |
|---|---|---|---|
| anchor | 0 | 0.046 (0.009) | 0.050 (0.010) |
| orientation | 0.25 | 0.842 (0.016) | 0.922 (0.012) |
| orientation | 0.50 | 0.878 (0.015) | 0.952 (0.010) |
| orientation | 0.75 | 0.872 (0.015) | 0.936 (0.011) |
| orientation | 1.00 | 0.840 (0.016) | **0.924** (0.012) |
| shape | 0.25 / 0.50 / 0.75 / 1.00 | 0.692 / 0.816 / 0.800 / 0.796 | 0.040 / 0.024 / 0.020 / 0.022 |
| magnitude | 0.25 / 0.50 / 0.75 / 1.00 | 0.006 / 0.014 / 0.018 / 0.024 | 0.004 / 0.000 / 0.000 / 0.000 |
| translation | 0.25 / 0.50 / 0.75 / 1.00 | 0.038 / 0.028 / 0.050 / 0.042 | 0.042 / 0.034 / 0.046 / 0.040 |

The per-design power figure is `report/design_point_power.png` (both columns, every mode and statistic).

**Reading.** Net displacement raises orientation power at every effect at both designs. It does so without leaking
into the other constructions: under size, shape and location changes its `angle` rate is at most 0.068 at every
effect at both designs, while PC1's `angle` responds to the shape construction (Phase 6 up to 0.138; Phase 5
0.69–0.82). At Phase 5 PC1's orientation curve is flat at 0.84–0.88 from e = 0.25; net displacement's is flat at
0.92–0.95. Neither reaches the concatenation ceiling (1.00).

**Paired outcomes** (computed from `merged.jsonl`; same dataset, same permutations; counts of 500):

| design | cell | both reject | PC1 only | nd only | neither |
|---|---|---|---|---|---|
| Phase 6 | orientation 1.0 | 141 | 9 | 183 | 167 |
| Phase 6 | shape 1.0 | 9 | 60 | 18 | 413 |
| Phase 5 | orientation 1.0 | 393 | 27 | 69 | 11 |
| Phase 5 | shape 1.0 | 6 | 392 | 5 | 97 |

The orientation gain is not a reshuffle: nearly every PC1 rejection is also a net-displacement rejection (141 of
150 at Phase 6, 393 of 420 at Phase 5). The shape-construction leak is PC1's, and net displacement does not take it on.

## 3. The rule's criteria

From `column_decision.json` / `.csv`. Bound 0.0695 everywhere (R = 500, no failed replicate).

| criterion | pair | Phase 6 PC1 | Phase 6 nd | Phase 5 PC1 | Phase 5 nd |
|---|---|---|---|---|---|
| target | orientation/`angle`, e = 1.0 | 0.300 | 0.648 | 0.840 | 0.924 |
| anchor | `delta` / `angle` / `shape` | 0.026 / 0.024 / 0.020 | 0.026 / 0.050 / 0.020 | 0.036 / 0.046 / 0.048 | 0.036 / 0.050 / 0.048 |
| control | translation/`angle`, e = 1.0 | 0.028 | 0.068 | 0.042 | 0.040 |
| control | magnitude/`angle`, e = 1.0 | 0.008 | 0.030 | 0.024 | 0.000 |
| control | shape/`angle`, e = 1.0 | **0.138** | 0.054 | **0.796** | 0.022 |
| protected | magnitude/`delta`, e = 1.0 | 1.000 | 1.000 | 1.000 | 1.000 |
| protected | shape/`shape`, e = 1.0 | 0.228 | 0.228 | 0.996 | 0.996 |

- **Phase 6, `gain`:** candidate − baseline = +0.348 against a threshold of 2 × pooled SE = 0.059. Passed.
- **Phase 5, `no_loss`:** baseline − candidate = −0.084 against 0.040. Passed (it is a gain of 0.084, just over
  4 pooled SE, but the rule asked only for no loss).
- **Anchor, controls, protected:** passed at both designs. The baseline's shape/`angle` rates (bold) exceed the
  bound at both designs; the rule judges only the candidate, and the table shows both, as the design required.
- **Protected pairs carry no evidence here.** Between two `pls_cv` columns, `delta` and `shape` are the same
  records' statistics (§6), so the protected criterion passes with loss exactly 0 by construction. It is declared
  for the rule's generality, not as a finding.
- **Translation control at Phase 6:** 0.068 against 0.0695, one rejection under the bound, as in the Type I study
  (0.068). It is the narrowest margin in the comparison (§7).

## 4. Covariates: eigengap, `angle` null width, selected rank

From `design_point_operating.csv` (`angle` rows; median pooled eigengap; median and IQR of the per-replicate `angle`
null q95, in degrees; CV-selected rank median (range)). Both columns share the integration, so the eigengap and rank
are the same records' values.

| design | cell | median eigengap | PC1 null q95 (IQR) | nd null q95 (IQR) | selected rank |
|---|---|---|---|---|---|
| Phase 6 | anchor | 0.111 | 84.0° (115.4) | 12.1° (26.6) | 2 (2–3) |
| Phase 6 | orientation 0.5 | 0.160 | 56.3° (104.9) | 30.6° (8.8) | 3 (2–3) |
| Phase 6 | orientation 1.0 | 0.167 | 62.6° (104.1) | 32.1° (7.0) | 4 (2–9) |
| Phase 6 | shape 1.0 | 0.126 | 50.7° (69.7) | 28.7° (10.4) | 3 (2–3) |
| Phase 5 | anchor | 0.053 | 5.5° (4.9) | 0.7° (0.3) | 3 (2–3) |
| Phase 5 | orientation 1.0 | 0.049 | 11.3° (13.9) | 1.2° (0.4) | 3 (3–4) |
| Phase 5 | shape 1.0 | 0.061 | 5.9° (5.1) | 1.0° (0.4) | 3 (2–4) |

Net displacement narrows the `angle` null and removes most of its across-replicate dispersion, the PC1 sign/flip
tail. Two consequences for power:

**Eigengap no longer gates power.** Orientation/`angle` power within eigengap terciles at e = 1.0
(`eigengap_stratified_power.csv`):

| design | tercile | eigengap range | PC1 | net displacement |
|---|---|---|---|---|
| Phase 6 | low | 0.013–0.127 | 0.114 (0.025) | 0.731 (0.034) |
| Phase 6 | middle | 0.127–0.213 | 0.449 (0.038) | 0.784 (0.032) |
| Phase 6 | high | 0.213–0.441 | 0.337 (0.037) | 0.428 (0.038) |
| Phase 5 | low | 0.002–0.038 | 0.754 (0.033) | 0.946 (0.017) |
| Phase 5 | middle | 0.038–0.064 | 0.886 (0.025) | 0.922 (0.021) |
| Phase 5 | high | 0.064–0.181 | 0.880 (0.025) | 0.904 (0.023) |

PC1's low-tercile collapse (Phase 6 0.114; the Phase 6 report's "the cohort's own eigengap carries an orientation
statement") is gone: net displacement does not need a well-determined principal axis. At e ≤ 0.75 its Phase 6 terciles
are flat (0.08–0.21).

**Rank still gates power at the Phase 6 maximum.** Computed from `merged.jsonl` (`integration_metadata.selected_lv`):
at Phase 6 orientation e = 1.0, 257 of 500 replicates retain more than three components (as in the parent study).

| selected rank | replicates | PC1 | net displacement |
|---|---|---|---|
| ≤ 3 | 243 | 0.082 | 0.280 |
| > 3 | 257 | 0.506 | 0.996 |

At e ≤ 0.75 every replicate retains at most three components and power is 0.084–0.178. The estimator raises power in
both rank strata. The headline 0.648 is the mixture, and it depends on how often CV retains the group-specific
variation the maximal construction creates, which the Phase 6 report found a six-donor change moves from 51 % to 8 %.
At Phase 5 the rank is 2–4 throughout and power does not depend on it.

## 5. Acceptance targets (baseline column, descriptive)

From `acceptance_report.csv`. The profiles declare upper-bound Type I and specificity targets on the baseline
column; the gate is off and the column rule does not read them.

- Phase 6: all eleven targets met.
- Phase 5: ten of eleven met. **translation/`shape` in the `type_i_baseline` cell is 0.070 against 0.0695** (one
  rejection over), failing both `type_i_control` and `specificity[translation, shape]`. That cell is independently
  seeded, measured in the baseline column, and `shape` does not depend on the orientation estimator. The Type I
  study's Phase 5 `type_i_baseline` translation cell recorded the same rate (0.036 / 0.056 / 0.070), and the
  matched-seed translation cells here are 0.030–0.046 on `shape` at every effect. It is a Monte Carlo excursion at
  the bound in a statistic this comparison does not touch, not an estimator effect.

## 6. Consistency checks

From `consistency_checks.json` (`scripts/net_displacement_power_checks.py`, checks version 2).

**(a) `delta` and `shape` are identical between the two columns.** The estimator enters only `_estimate_orientation`,
so the two columns' `delta` and `shape` cannot differ. Over 8,500 paired records per design (17,000 comparisons):
**0 p-value mismatches** at both designs. Observed statistics are bit-identical in 11,668 (Phase 6) and 9,221
(Phase 5) comparisons and differ in the rest by at most 9.0 × 10⁻¹¹ and 5.4 × 10⁻¹⁰ relative: BLAS rounding across
separate processes, inside the checks' 1 × 10⁻⁹ tolerance and never enough to move a p-value.

**(b) The baseline column reproduces the parent studies.** It shares their generator seeds, evaluation parameters
and RRPP seed, so wherever a parent record of the same mode, effect and replicate exists, its seed and p-values must
be equal.

| design | parent | records compared | seed / p-value mismatches | observed bit-identical | max rel. diff. |
|---|---|---|---|---|---|
| Phase 6 | `phase6-small-n-2026-10-09` | 6,500 | 0 / 0 | 5,285 | 3.1 × 10⁻⁸ |
| Phase 6 | `phase6-small-n-magnitude-2026-10-09` | 2,000 | 0 / 0 | 1,167 | 3.1 × 10⁻⁸ |
| Phase 5 | `phase5-2026-09-10` | 6,500 | 0 / 0 | 2,981 | 4.7 × 10⁻¹⁰ |
| Phase 5 | `phase5-magnitude-2026-10-08` | 1,500 | 0 / 0 | 897 | 5.1 × 10⁻¹¹ |

8,000 of 8,500 (Phase 6) and 7,500 of 8,500 (Phase 5) baseline records have a parent. The uncovered cells are the
new ones the design names: Phase 6 magnitude 0.75; Phase 5 magnitude 0.5 and 0.75. The 2,000 magnitude records of
`phase5-2026-09-10` are excluded from the comparison: that run generated magnitude under `magnitude_kind: all`, so
its magnitude datasets differ despite equal seeds. The joint-construction parent `phase5-magnitude-2026-10-08`
covers magnitude 0.25 and 1.0.

The PC1 column's numbers are therefore the parent studies' recorded rates (e.g. Phase 6 orientation/`angle` 0.300,
Phase 6 shape/`angle` 0.138), and the two columns are paired measurements of the same 500 datasets per cell.

**Checks history.** The first pass (checks version 1, `db5ab69`) exited 1 on both designs for two reasons, both in
the checks rather than the records: (a) required bit-identical observed statistics and flagged the rounding above; (b)
compared the Phase 5 joint-magnitude cells against the `magnitude_kind: all` parent (4,000 spurious p-value
mismatches, all in magnitude). `1ccded5` adds the 1 × 10⁻⁹ tolerance and pairs parent records only within one
magnitude construction (`truth_metadata.transform.magnitude_kind`); no record changed.

## 7. Decision and interpretability

`comparison_decision.json` (rule version 1; inputs `phase6/report/column_decision.json` sha256 `c980fafb…`,
`phase5/report/column_decision.json` sha256 `bdc3135e…`):

- Phase 6 `candidate_qualifies` (`gain`); Phase 5 `candidate_qualifies` (`no_loss`);
- **verdict `adopt_net_displacement`**.

The verdict selects the follow-up change: **flip `SimulationEvaluationParams.orientation_estimator`'s production
default to `net_displacement`** in a separate change that re-pins the committed signatures, as
`default-joint-magnitude-construction` did. Until that change lands, the default and every committed result stay
PC1.

**What `angle` means changes.** Under net displacement, `angle` is the angle between the two groups' net
displacement vectors (last stage mean − first stage mean) in the latent space, not between their principal axes. It
reads only the endpoint stages. The paper and the case study must describe it that way.

**Revised `angle` statement for the SEA-AD cohort (n = 80, ≥30 nuclei, RNA + ATAC)**, replacing the `angle` row of
the [Phase 6 report](phase6-small-n-2026-10-09.md) §8 once the default flips:

| statistic | rejection may be reported as evidence? | non-rejection may be reported as evidence? |
|---|---|---|
| `angle` | **Yes, that the groups' net directions of change differ.** Type I is controlled (anchor 0.050; Type I study 0.048–0.058 in every space). Unlike PC1, the size, shape and location constructions do not move it: 0.024–0.068 at every effect, all within the bound, against PC1's up to 0.138 under shape. The translation construction sits at 0.058–0.068, near the bound, so a rejection alongside a location difference warrants care. The cohort's eigengap no longer conditions the statement. | **No.** Power is 0.084–0.178 at e ≤ 0.75. At the construction maximum it is 0.648 overall, but 0.280 when CV retains three or fewer components (243 of 500 replicates) and 0.996 above that. Report the cohort's selected rank beside any `angle` result. A non-rejection is uninformative at n = 80. |

At the Phase 5 design (n = 1200) `angle` has power 0.92–0.95 at every effect and is specific: shape → `angle`
falls from 0.80 to 0.02. The orientation construction is still detected by `shape` (0.998–1.000) and `delta`
(0.80–0.91) at Phase 5, and at Phase 6 by `shape` (0.464) and `delta` (0.354) at e = 1.0. Those are properties of the
construction and of the other two statistics, which this change does not touch. The `delta` and `shape` rows of the Phase 6 statement
are unchanged.

## 8. Provenance

- **Configs:** `examples/trajectory_power_study/net_displacement_power_phase6.json` (sha256 `6f5dc3b2…`) and
  `…_phase5.json` (sha256 `245433ec…`), identity-pinned in `tests/data/committed_config_identity.json`.
- **Code:** `db5ab691e97b9c21e6a8876f740bfb24fc699e67` (branch `feat/net-displacement-power-comparison`) for the run,
  merge, reports and comparison decision; `1ccded5` for the consistency checks (the change touches only
  `scripts/net_displacement_power_checks.py`). Clean tree on the cluster. Python 3.11.16, motco 0.6.0, numpy 2.3.5,
  scikit-learn 1.8.0, scipy 1.16.3, pandas 2.3.3, uv 0.12.10.
- **Host:** `ing`, partition `512x1024` (AMD EPYC 7662). Phase 6: 50 single-CPU shards, 2 GB. Phase 5: 100 shards,
  4 GB. BLAS pinned to one thread, `STUDY_N_JOBS` unset, nodes n4 and n11 excluded.
- **Records:** 18,000/18,000 units per design, all 36 cells present, 0 duplicates, 0 failed, 0 censored, every
  signature equal to the enumerated cell's, `n_jobs` 1 everywhere.
- **Cost:** Phase 6 64.4 recorded core-hours (1 h 53 m wall), against a 74.5 projection. Phase 5 319.2 (6 h 38 m
  wall), against 313.0 projected and 375.7 with margin, under the predeclared 450 go/no-go.
- **Merged records:** `phase6/merged.jsonl` sha256 `f71cbaf1…`, `phase5/merged.jsonl` sha256 `9eacdf37…`
  (gitignored; verified identical on the cluster and the workstation).

The SLURM launch lines, timing jobs and unit timings are in `phase6/PROVENANCE.txt` and `phase5/PROVENANCE.txt`.

## 9. Reproduction

From the shards:

```bash
B=results/net-displacement-power-2026-10-14
for d in phase6 phase5; do
  python scripts/motco_study.py merge --out-dir $B/$d
  python scripts/motco_study.py report \
    --config examples/trajectory_power_study/net_displacement_power_$d.json --out-dir $B/$d
done
python scripts/net_displacement_power_follow_up.py \
  --phase6 $B/phase6/report/column_decision.json \
  --phase5 $B/phase5/report/column_decision.json \
  --out $B/comparison_decision.json
python scripts/net_displacement_power_checks.py --merged $B/phase6/merged.jsonl \
  --parent results/phase6-small-n-2026-10-09/merged.jsonl \
  --parent results/phase6-small-n-magnitude-2026-10-09/merged.jsonl \
  --out $B/phase6/consistency_checks.json
python scripts/net_displacement_power_checks.py --merged $B/phase5/merged.jsonl \
  --parent results/phase5-2026-09-10/merged.jsonl \
  --parent results/phase5-magnitude-2026-10-08/merged.jsonl \
  --out $B/phase5/consistency_checks.json
```

The rank-stratified and paired-outcome numbers (§2, §4) come from `merged.jsonl`:

```python
from pathlib import Path
import pandas as pd
from motco.simulations.grid import read_replicate_results
from motco.simulations.study.spectrum import record_design_point

def angle_frame(merged):
    rows = []
    for r in read_replicate_results(Path(merged)):
        if r.phase == "type_i_baseline":
            continue
        rows.append({
            "est": (record_design_point(r) or {}).get("evaluation.orientation_estimator") or "pc1",
            "mode": r.cell_metadata["trajectory_mode"], "e": r.cell_metadata["effect_size"],
            "rep": r.replicate_index, "rej": r.p_values["angle"] < 0.05,
            "lv": r.integration_metadata["selected_lv"],
        })
    return pd.DataFrame(rows)

d = angle_frame("results/net-displacement-power-2026-10-14/phase6/merged.jsonl")
o = d[(d["mode"] == "orientation") & (d.e == 1.0)]
print(o.assign(high=o.lv > 3).groupby(["high", "est"]).rej.agg(["mean", "size"]))      # §4 rank table
print(pd.crosstab(*o.pivot(index="rep", columns="est", values="rej").T.values))      # §2 paired counts
```
