# Proposal

## Why

The project's originally planned simulation grid (the proposal's Table 1) cannot be run as written. It also no longer matches what the studies have shown:

- crossing its rows is out of reach (3,456 null combinations before any power cell);
- its effect sizes sit far below InterSIM's stage effect;
- its power rows do not map to trajectory components;
- two of its "power" rows are null scenarios.

What the paper still lacks is design guidance: how power for each trajectory component depends on sample size, number of stages and number of blocks, under the **current** production defaults. The only sample-size sweep (the Phase 5 design-point pilot, 2026-09-08) ran orientation and translation only, with PC1, without the rank floor and at 100 replicates.

## What Changes

- **Four pilot study configurations** under `examples/trajectory_power_study/`, all at the production defaults, named explicitly: PLS with CV rank, net displacement, rank floor 2 and `joint` magnitude.
  - **Crossed core**, as a magnitude profile and a geometry profile:
    - design grid n ∈ {80, 200, 500, 1000} × stages ∈ {3, 4} × blocks ∈ {2, 3};
    - magnitude at e ∈ {0.01, 0.02, 0.05, 0.1, 0.25};
    - orientation, shape and translation at e ∈ {0.25, 0.5, 1.0};
    - a zero-effect anchor per column.
  - **One-factor-at-a-time arms**, as a magnitude profile and a geometry profile, off the baseline n = 200, 3 stages, 2 blocks:
    - group ratio {0.65, 0.8};
    - unequal stage proportions (0.5/0.3/0.2);
    - `p_dmp` {0.05, 0.15};
    - joint stage effect δ {0.05, 0.1, 0.25, 0.5, 1.0} (0.05 and 0.1 added at the approval gate: the smoke showed the stage signal saturated down to 0.25);
    - baseline continuity ρ {0.5, 0.8}.
  - All four run at pilot grade: 100 replicates × 199 permutations, about 47,200 runs in total.
- **Stage-signal table.** Study reporting writes `design_point_stage_signal.csv` beside `design_point_operating.csv`. It summarizes each column's recorded held-out CV AUROC for stage. This is the stage-separability number that compares directly with a real analysis (SEA-AD 0.64–0.70). `design_point_operating.csv` itself is unchanged.
- **Design-guidance figure.** A script reads the two core reports and plots each mode's target-statistic power against n. It draws one panel per mode and block count, and one line per stage count. The engine's existing design-point figure would mix the block columns.
- **Config-pin tests distinguish historical from new configs.** The rank-floor and orientation pin tests assume every committed config is historical. In future they hold only the configs in their pre-flip fixtures to the historical pins. Configs written since then must name their estimator and floor explicitly, and their identity is pinned in a fixture of their own.
- **Pilot run and dated findings report**, with a section that maps each row of the originally planned grid to what was run instead, and why.
- Out of scope:
  - group-specific dispersion (generator change; a separate change);
  - the paper-grade rerun of a chosen slice (decided from this pilot);
  - the case-study-configuration Type I check;
  - two-stage designs, added noise and supervised feature filtering, all dropped with reasons in the report.

## Capabilities

### New Capabilities

None.

### Modified Capabilities

- `trajectory-power-study`:
  - adds the design-guidance pilot profiles, the stage-signal table, the design-guidance figure and the dated pilot report;
  - narrows the PC1-estimator and unfloored-rank pin requirements to configurations whose results predate those defaults;
  - requires later configurations to name both settings explicitly and to have their identity pinned.

## Impact

- `src/motco/simulations/study/spectrum.py` and `report.py`: new stage-signal frame and CSV, written only when design-point cells exist. Existing outputs are byte-identical.
- `scripts/`: new design-guidance figure script.
- `examples/trajectory_power_study/`: four new configs. The README lists them.
- `tests/`:
  - `test_rank_floor_default_pins.py` and `test_orientation_default_pins.py` are restricted to their fixtures, plus a check that every config is accounted for;
  - new identity fixture and test for the design-guidance configs;
  - stage-signal report tests.
- Compute: roughly 330 core-hours on `ing`. This is extrapolated from the Phase 5 and Phase 6 timings at 999 permutations, so it is an upper bound; a smoke run fixes it.
- Docs:
  - `docs/roadmap.md`: new entry; Phase 5 sensitivity deliverables;
  - `AGENTS.md`: study bullet;
  - `src/motco/simulations/study/README.md`;
  - a dated report under `docs/reports/`;
  - records under `results/design-guidance-pilot-<date>/`.
