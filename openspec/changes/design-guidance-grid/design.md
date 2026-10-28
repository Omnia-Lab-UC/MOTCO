# Design

## Context

See proposal.md for the motivation. The study engine already does everything this pilot needs:
- crossed design grids over any `generator.<field>`, including `n_samples` and `n_stages`;
- the nested layer axis `evaluation.integration_params.layers`;
- explicit design points.

A scratch draft of the four configs loads and enumerates on `main` (`843f057`) with no engine change:

| Profile | Cells | Runs at 100 replicates |
|---|---:|---:|
| core magnitude | 98 | 9,800 |
| core geometry | 162 | 16,200 |
| arms magnitude | 68 | 6,800 |
| arms geometry | 112 | 11,200 |
| **Total** | | **44,000** |

Three things block or degrade the work today:

1. **The pin tests reject any new config.** `tests/test_rank_floor_default_pins.py` and `tests/test_orientation_default_pins.py` require every file in `examples/trajectory_power_study/` to be in the 19-config pre-flip fixture, and to pin `min_components: 1`.
2. **Stage signal is recorded but never reported per column.** Each PLS record carries `integration_metadata["cv_mean_auroc"]`, which `phase4.py` already reads, but no design-point output summarizes it.
3. **The engine's design-point figure would mislead here.** `render_design_point_power` plots two numeric axes; `layers` is not numeric. On the core grid it would mix the two- and three-block columns in one line.

## Goals / Non-Goals

**Goals:**
- Pilot-grade, per-component design guidance under the production measurement.
- No change to any existing output, fixture or committed config.

**Non-Goals:**
- **Per-mode effect sizes in the engine.** Splitting each grid into a magnitude profile and a geometry profile does the same job.
- **Changing `render_design_point_power`.**
- **Committing `merged.jsonl`.** As for earlier studies, only `report/` and `PROVENANCE.txt` go under `results/`.

## Decisions

### D1. Four profiles: magnitude and geometry, for both core and arms

Each config has one effect list. Magnitude needs small effects: at n = 80 with two blocks, `delta` power is 0.206 at e = 0.05, 0.542 at 0.10 and 0.956 at 0.25, and it saturates sooner at larger n. Orientation and shape sit near the null below e = 0.25 at n = 80. So:
- magnitude runs at e ∈ {0.01, 0.02, 0.05, 0.1, 0.25};
- the geometry modes run at e ∈ {0.25, 0.5, 1.0}. e = 0.75 is dropped to save a cell per mode.

The Phase 6 study split the same way.

Alternative considered: one list {0.05, 0.1, 0.25, 0.5, 1.0} for all modes. It costs about 40 % more per column, and still cannot reach the magnitude effects that matter at n ≥ 500.

### D2. One seed family across all four profiles

All four profiles use `base_seed` 1028 and family `design-guidance`. The magnitude and geometry modes are then matched on the same baseline draws, exactly as they would be within one config.

Consequence: at a given column, the magnitude and geometry profiles' anchors have different cell ids but draw the same generator seed. They evaluate the same datasets. Verified on the draft: the n = 80, four-stage, three-block anchor draws seed 1476637963 at replicate 0 in both. So:
- the report reads the geometry profile's anchor;
- it uses the magnitude profile's anchor only as a reproduction check;
- it never pools the two as independent replicates.

### D3. Baseline: n = 200, three stages, two blocks

- **Three stages** match the case study and the planned AD analyses (CN/MCI/AD).
- **Two blocks** match SEA-AD.
- **n = 200** is the smallest planned size where the geometry modes leave the null without magnitude saturating. That gives the arms room to move power both ways.

Alternatives considered:
- n = 80: the arms would read orientation and shape at the null.
- n = 500: magnitude saturates at every arm.

### D4. Arm values, and why `p_dmp` 0.15 rather than 0.2

At three stages, `p_dmp` = 0.2 fails the enumeration headroom check: orientation saturates at e = 0.89 and translation at 0.52. The config would therefore be refused at load. `p_dmp` = 0.15 fits, with saturating effects of 1.39 for orientation and 1.12 for translation, even under 0.5/0.3/0.2 stage proportions.

Other values:
- **Stage effect.** Spans 0.25–1.0 below InterSIM's 2.0. The rank-floor check showed stage signal vanishing near 0.05–0.1 at n = 80.
- **Group ratio, stage proportions, ρ.** Taken from the planned grid and the Phase 5 pilot.

An SEA-AD-style group × stage table is not repeated. The Phase 6 net-displacement power profile already measured it.

### D5. Stage signal as its own CSV, not new operating-table columns

`design_point_stage_signal.csv` is written beside `design_point_operating.csv` whenever that table is. It has one row per (design point, mode, effect), with:
- `n_auroc_available`;
- the mean and the q33 / median / q67 of `cv_mean_auroc`, the same quantiles the eigengap columns use.

Bucketing reuses `resolve_operating_by_design_point`'s keys, so the rows line up with the operating table.

Alternative considered: append the columns to the operating table, as `component_selection` was. Rejected for two reasons:
- it would change the eight committed `design_point_operating.csv` files on any regeneration;
- AUROC does not vary by statistic, so it would be repeated three times per row.

### D6. Design-guidance figure as a script

`scripts/design_guidance_figure.py --magnitude <report dir> --geometry <report dir> --out-dir <dir>` writes:
- `design_guidance_magnitude.png`, `design_guidance_orientation.png` and `design_guidance_shape.png`. Each has panels for {3, 4} stages × {2, 3} blocks, n on a log axis, and one line per effect with ±1 SE bars;
- `design_guidance_control.png`: anchor, plus translation at its largest effect, per statistic, with the α + 2·SE bound;
- `design_guidance_power.csv`.

In the CSV, the `layers` coordinate is read as its CSV spelling: the two-block list, or empty for all three blocks. It is labelled "2 blocks" or "3 blocks".

The figure joins two studies, which the engine never does. Follow-up scripts are the precedent (`scripts/net_displacement_follow_up.py`).

Alternative considered: a faceted engine figure gated on `n_samples` + `n_stages` axes. Rejected: it adds an engine output that is useful for one study family.

### D7. Pin tests: historical set and post-floor set

`tests/data/orientation_default_identity.json` stays as captured. A new fixture, `tests/data/post_floor_config_identity.json`, has the same shape (`capture_revision`, `note`, `configs` → cell id → signature and first three seeds). It holds every config committed after the rank floor, starting with the four profiles. Later changes append to it.

In both pin tests:
- `test_fixture_covers_every_committed_config` becomes: the two fixtures are disjoint, and their union equals the directory listing. The historical set stays at exactly 19.
- `test_every_committed_config_names_its_rank_floor` keeps `== 1` for historical configs. Post-floor configs must name a non-null floor.

A new `tests/test_post_floor_config_identity.py` checks identity against the new fixture. `tests/test_design_guidance_profiles.py` checks profile structure:
- column counts;
- block columns share seeds;
- one-setting arms;
- joint δ;
- explicit production names.

### D8. Acceptance: advisory anchor decision only

Each profile declares:
- a `type_i` entry;
- `acceptance.anchor_decision` with α 0.05, k 2 and all three statistics. At R = 100 the bound is 0.0936.

There are no power targets, no gate and no report contract. The pilot is descriptive and makes no floor claims.

### D9. Execution

1. Smoke a 2-replicate × 19-permutation copy of each profile locally, in the scratchpad and never committed. It runs shard → merge → report → figure script end to end, and times units at each n.
2. Size the shard counts from those timings and run the four arrays on `ing` with `scripts/motco_study_array.sbatch`.
3. Commit records under `results/design-guidance-pilot-<date>/<profile>/` (`report/` and `PROVENANCE.txt`) and the figures under `results/design-guidance-pilot-<date>/figures/`.

## Risks / Trade-offs

- **Smallest cells.** At n = 80, four stages, cells hold 10 samples, so stage-stratified double CV trains on about 5 per cell. CV can fail or select erratically. → The smoke run includes that column. Failures are recorded as failed replicates and reported per column, not hidden.
- **The engine has never run `n_stages` as a design axis.** → The smoke run exercises report generation end to end. Any engine bug found is fixed in this change, with a test.
- **Compute is dominated by n = 1000** (about half the estimated 330 core-hours). → Use the smoke timings to set shards. If the estimate doubles, submit the n = 1000 shards as their own array, so the rest of the pilot can be reported without waiting on them.
- **Pilot precision is limited.** At 100 replicates the SE is ≤ 0.05 and the anchor bound is loose (0.0936). Translation's `angle` rate sits near the paper-grade bound (0.058–0.068). → The report reads rates with their SE, flags translation, and defers claims to the paper-grade slice.
- **Block count is confounded with feature count** (498 vs 498 + protein). → The report states it. A block-weighted comparison is roadmap item 5.

## Open Questions

- **Whether the two profiles' duplicated anchors also share RRPP draws.** If they do, their p-values must be identical. If not, they are the same data with different permutations. Either outcome fits D2's reading, and the implementation records which.
- **How the funded proposal is named in the report's planned-grid mapping.** The default is "the project's originally planned grid (Table 1)". The wording can change without affecting specs or tasks.
