# Notes

## Implementation (2026-11-23)

- The driver is `scripts/univariate_baseline_check.py`, with `scripts/univariate_baseline_array.sbatch` and `tests/test_univariate_baseline_check.py`. The unit counts are 4,400 for `parent`, 12,800 for `stage_span` and 16,500 for `surgery`, with 71 signal cells (44 + 27) and 36 graded cells.
- **The `surgery` RRPP seed is 3 for every unit.** The slice's configs set `evaluation.seed: 3`, and the study engine uses the cell's own seed when one is set (`grid.run_simulation_replicate`). So every slice unit drew its permutations at seed 3 and only the generator seed varied. The driver resolves the seed the same way, so the slice's counts can reproduce. The record's `rrpp_seed` is that resolved seed.
- **The split test runs once per unit, uncorrected** (`interaction_split_test`), not through `span_split_check.split_record`. That function also reruns the test under `df` for an SS diagnostic this change does not read. The reproduced parts read `p_value_f` from the uncorrected run, as the scan validation's comparator table did.
- **`--recorded-spaces` is bound in the settings** (default: `concat` and `pca_matched`). It is the first compute lever of the design's risk list. A space left out records its arms as unavailable, and nothing judged reads them.

## Per-unit probe before the smoke (one unit each, 999 permutations, BLAS pinned, Python 3.11.15)

| Unit | fit | evaluation | split | scan | `concat` | `pca_matched` | unit |
|---|---|---|---|---|---|---|---|
| parent `n200` magnitude r0 (498 features) | 10.6 s | 0.4 s | 0.1 s | 0.1 s | 36.6 s | 0.4 s | 48.5 s |
| stage_span `n200` orient e0.25 φ0 r0 (658) | 12.7 s | 0.5 s | 0.1 s | 0.1 s | 76.1 s | 0.5 s | 90.3 s |
| surgery `n200_three_blocks` orientation e0.5 r0 (658) | 11.2 s | 0.5 s | 0.1 s | 0.1 s | 75.7 s | 0.5 s | 88.3 s |
| stage_span `n500` orient e0.5 φ90 k20 r0 (658) | 18.7 s | 0.5 s | 0.5 s | 0.3 s | 76.7 s | 0.6 s | 97.6 s |

`concat` costs about the same at every sample size, because its `shape` statistic aligns the stage means in the ambient dimension (a 658 × 658 problem) on each of the 999 draws. Every other arm together costs less than the cross-validation.

## Local smoke (task 5.1, 2026-11-23)

Every design and column at 2 replicates per cell (2 control replicates), 999 permutations, both recorded spaces, 12 workers on the workstation's 8 physical cores (Ryzen 7 5800X, 16 threads), BLAS pinned in the shell: 226 units, 0 failed. `merge`, both reproduction checks and `decide` ran.

- **No planted fixed point failed** in any `stage_span` design (`failures_by_cell` and `failures_by_error` are empty in every profile's accounting).
- **Both reproductions differ only by size.** `parent` compared 196 rows and `surgery` 99 (0 unresolved); every row differs because the smoke counts 2 units against the committed 500. The decision reads `reproduction_failed` for that reason, all 71 signal cells are `unresolved`, and the grade (`names_kind`) is read from 2 units per cell, so none of it is a finding.

Wall time per profile: `parent` 430 s, `stage_span` 1,720 s, `surgery` 1,119 s. Mean seconds per unit from the records:

| Profile | Design | unit | `concat` | fit | rest without `concat` |
|---|---|---|---|---|---|
| `parent` | `n200` / `seaad` | 91.2 / 84.1 | 68.8 / 64.2 | 19.9 / 17.7 | 21.2 (profile mean) |
| `stage_span` | `n80` / `n200` / `n500` / `n200_delta0.5` | 180.9 / 191.5 / 225.5 / 173.4 | 158.2 / 158.3 / 170.0 / 141.9 | 19.8 / 29.7 / 49.7 / 28.5 | 35.7 |
| `surgery` | `n80` / `n200` / `n500` three blocks | 187.4 / 191.8 / 197.0 | 165.3 / 161.4 / 146.4 | 19.4 / 27.3 / 45.5 | 34.4 |

Medians by arm, same units: `evaluation` 1.0–1.2 s, `split` 0.1–0.3 s, `scan` 0.1–0.3 s, `pca_matched` 0.7–1.1 s. The smoke's oversubscribed cores ran each unit about twice as slowly as the single-unit probe above (the `n200` fit 19.9 s against 10.6 s). The probe's cost matches `ing`'s unpacked cores: the design-guidance slice ran its production units there at a median of 12 s.

Signal cells where an arm rejected in both units / in at least one (2 replicates, so only a smoke reading):

| Profile | Design | cells | `production` | `scan` | `fdr` | `total` | `pca_matched` | `concat` |
|---|---|---|---|---|---|---|---|---|
| `stage_span` | `n80` | 11 | 4 / 11 | 8 / 9 | 8 / 9 | 4 / 9 | 1 / 3 | 3 / 7 |
| `stage_span` | `n200` | 11 | 11 / 11 | 9 / 10 | 9 / 10 | 4 / 6 | 4 / 4 | 8 / 10 |
| `stage_span` | `n500` | 11 | 11 / 11 | 11 / 11 | 11 / 11 | 6 / 10 | 6 / 6 | 10 / 10 |
| `stage_span` | `n200_delta0.5` | 11 | 1 / 6 | 4 / 8 | 4 / 7 | 1 / 2 | 2 / 3 | 3 / 6 |
| `surgery` | `n80_three_blocks` | 9 | 0 / 9 | 6 / 7 | 6 / 7 | 6 / 6 | 1 / 4 | 3 / 5 |
| `surgery` | `n200_three_blocks` | 9 | 6 / 9 | 7 / 8 | 7 / 8 | 6 / 6 | 1 / 6 | 5 / 7 |
| `surgery` | `n500_three_blocks` | 9 | 9 / 9 | 8 / 8 | 8 / 8 | 6 / 7 | 1 / 7 | 7 / 8 |

The null cells (anchor and translation) rejected once in `stage_span` (`total`, 1 of 16) and once each for `scan` and `fdr` in `surgery` (1 of 12).

### Projection

Units × the smoke's mean seconds per unit, as an upper bound (packed nodes), and half of it, as the bound at the probe's unpacked cost:

| Profile | units | with `concat` | without `concat` |
|---|---|---|---|
| `parent` | 4,400 | 54–107 h | 13–26 h |
| `stage_span` | 12,800 | 343–686 h | 64–127 h |
| `surgery` | 16,500 | 440–880 h | 79–158 h |
| total | 33,700 | **837–1,673 h** | **156–311 h** |

Without `concat` the run sits at the design's estimate (130–200 core-hours, proposal); with it, 4–8 times above. `concat`'s cost does not fall with n: its `shape` aligns the stage means in the ambient dimension on every draw.

### Amendment and cap (approved by the user, 2026-11-23)

The design's first compute lever is the recorded `concat` and `pca_matched` arms. The smoke shows `pca_matched` costs about 1 s per unit and `concat` about 80 % of every unit. The proposal is to drop `concat` from the run (`--recorded-spaces pca_matched`, bound in the settings, no code change). The spec, design and proposal would drop `concat` from the recorded arms. The report would cite the latent-space comparison (`latent-space-comparison-2026-11-04`), which already measured `concat` against production on the stage-span construction across alignment. No judged arm reads `concat`, so the verdict, labels, grade and size readings are unchanged. Proposed cap: **350 recorded core-hours**, above the packed bound of 311 h.

**Approved** on 2026-11-23: `concat` dropped, cap 350 recorded core-hours, launch on `ing`. The driver's default is now `pca_matched` alone (`DEFAULT_RECORDED_SPACES`), so the run's settings carry `"recorded_spaces": ["pca_matched"]` without a flag; the spec, design (D3, open question) and proposal record the amendment.

## Run on `ing` (task 5.2)

Submitted 2026-11-23 from the clean worktree `/home1/tgonza/MOTCO-univariate-baseline`, detached at `a8a27dce822a4a9ae2281f2ad0e58cc7ea5174cd` (tree `2102d7ce80439b194a667cb14c0d41e6989399b3`), sharing `/home1/tgonza/MOTCO/.venv` with `PYTHONPATH` on the worktree's `src` (Python 3.11.16, numpy 2.3.5, scipy 1.16.3, scikit-learn 1.8.0, pandas 2.3.3). The driver's 63 tests passed there first. Run root `/home1/tgonza/MOTCO-univariate-baseline-run`, each profile's settings pre-created (hashes `1f0357bec8e2…` parent, `fe5a8425d7ac…` stage_span, `df600edfdab0…` surgery; `recorded_spaces` `["pca_matched"]`, 999 permutations).

| Profile | Job | Shards | Units |
|---|---|---|---|
| `parent` | 957651 | 40 | 4,400 |
| `stage_span` | 957652 | 100 | 12,800 |
| `surgery` | 957653 | 120 | 16,500 |

Partition `512x1024`, n11 excluded, one CPU and 4 GB per task, throttle 64 per array.

Outcome: 260 of 260 tasks `COMPLETED 0:0` (18:45:44–21:09:20; 176 on n2, 50 on n1, 34 on n4), 33,700 of 33,700 units complete, 0 failed. **236.18 recorded core-hours** (parent 14.01, stage_span 88.49, surgery 133.68) against the cap of 350; 246.6 allocated CPU-hours in `sacct`, startup included. Python startup stalled on NFS for 5–10 min at launch, as on earlier runs. Both reproductions passed: `parent` 196 of 196 rows and `surgery` 99 of 99. `decide`: **`complementary`** (admissible baseline arms `fdr` and `total`; `scan` inadmissible on `stage_span` size at n = 80; 16 / 21 / 34), grade **`detects_only`**. The decision regenerates byte for byte on the workstation from the committed aggregates. Per-unit records stay under `/home1/tgonza/MOTCO-univariate-baseline-run` (sha256 of each `merged.jsonl` in `PROVENANCE.txt`).

`uv run mkdocs build --strict` aborts on two warnings in `docs/reports/phase5-paper-grade-2026-09-10.md` (links to `results/` outside the docs tree). The same two warnings abort it on `main` (`bb58c8c`), and that report is frozen (`tests/test_phase5_magnitude_profiles.py`). The non-strict build of this branch emits no other warning, and the new report and the concept-page anchor resolve.
