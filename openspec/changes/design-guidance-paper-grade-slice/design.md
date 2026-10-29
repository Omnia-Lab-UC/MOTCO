# Design

## Context

See proposal.md for the motivation. The pilot report §7 fixes the slice; what this design settles is how the engine expresses its predeclared rules and how the run reproduces the pilot. Facts checked on `main` (`d7836d3`):

- **`acceptance.type_i` entries are phase-scoped.** `evaluate_targets` reads them only on summaries whose phase starts with `type_i_`: the baseline column's anchor and translation control in `type_i_baseline`, a separate seed derivation from the design columns. Translation at a design point is a `power_design` cell and no Type I entry sees it. `acceptance.anchor_decision` does evaluate every design point, and lists the other negative controls there descriptively, without a verdict. `acceptance.specificity` reads the top-effect cell of one mode across every column, so it cannot judge per column either.
- **The figure script derives its panels from the data** (`_panels` reads the stage counts and block counts present), so a table with one stage count renders one panel row; only the "(pilot; …)" title text is hard-coded. The crossed `design_grid.axes` form accepts a single-valued axis ("at least one value"), so `generator.n_stages: [3]` keeps the stage column in `design_point_operating.csv` and the script's three-axis contract.
- **Draft profiles enumerate on `main` without any engine change** (scratchpad drafts derived from the pilot's core profiles):

  | Profile | Cells | Phases | Units at 500 replicates |
  |---|---:|---|---:|
  | slice magnitude | 26 | 2 Type I + 4 primary + 20 design | 13,000 |
  | slice geometry | 62 | 2 + 10 + 50 | 31,000 |
  | stage-effect magnitude | 14 | 2 + 4 + 8 | 7,000 |
  | stage-effect geometry | 32 | 2 + 10 + 20 | 16,000 |

  No headroom failure. The geometry profile runs translation at every geometry effect, so a column holds 1 + 3 × 3 cells, not the 1 + 3 + 3 + 1 the pilot report §7 counted (see D2).
- **Seed reproduction is verified on the drafts.** With the pilot's `base_seed` 1028 and family `design-guidance`, every primary and design-point cell of the drafts draws the pilot core cell's generator seeds at replicate indices 0–99 (24 / 24 and 60 / 60 cells). The two Type I baseline cells per profile do not: non-matched phases are keyed on their own cell id, which carries `n_replicates`.
- **Cost basis.** The pilot timed units at 199 permutations on `ing`: 10.7 / 13.2 / 19.4 s at n = 80 / 200 / 500 with three stages and two blocks, about 1.2× with three blocks, the RRPP's 199 draws adding 2–5 s. At 999 permutations a unit is about 20–23 s at n = 80, 26–30 s at n = 200 and 40–48 s at n = 500: 44,000 units ≈ 380 core-hours. The pilot's geometry profile overran its estimate by 23 % (rank-inflated units), so the plan budgets 450.

## Goals / Non-Goals

**Goals:**
- Paper-grade rates (SE ≤ 0.022) with predeclared size verdicts at every column, evaluated by the engine rather than read off by hand.
- The slice reproduces the pilot's datasets at its first 100 replicates and the siblings reproduce each other, with committed tools that check both.
- Every engine addition byte-identical when unused.

**Non-Goals:**
- **Per-mode effect lists in the engine** (see D2).
- **Any new engine figure.** The design-guidance figure script gains a label, nothing else.
- **A report contract.** Attribution is off (as in the pilot), and the contract's driver outputs presuppose attribution records; the `n_jobs` guarantee it gives is replaced by the provenance check that every record carries `n_jobs` 1.
- **Committing `merged.jsonl`.** `report/`, `PROVENANCE.txt` and `figures/` only, as for every study.

## Decisions

### D1. Two sibling profiles derived from the pilot's core profiles

`design_guidance_core_magnitude_study.json` and `design_guidance_core_geometry_study.json` are the pilot's core profiles with the design grid narrowed to `n_samples` {80, 200, 500} × `n_stages` {3} × `layers` {two blocks, null}, the magnitude effects to {0.02, 0.05, 0.1}, 500 replicates, 999 permutations, the paper-grade acceptance block (D4) and new metadata. Everything else — baseline, CV knobs, `evaluation.seed` 3, attribution off, gate off — is unchanged, so a reader can diff a `_study` profile against its `_pilot` sibling.

The `_study` suffix follows `phase6_small_n_study.json` (paper grade) against `_pilot`.

Alternative considered: a single profile with per-mode effects. Rejected in the pilot (D1 there) and again here: the engine has one effect list per configuration, and the split is what makes the two anchors per column a bit-exact reproduction check.

### D2. Translation runs at every geometry effect

A configuration has one effect list, so translation runs at e ∈ {0.25, 0.5, 1.0} like orientation and shape: 6,000 units (about 55 core-hours, 13 % of the slice) beyond the pilot report's count, which assumed one translation cell per column. The pilot's geometry profile ran translation the same way; the control figure and the anchor decision read it at its largest effect, and the other two effects are descriptive controls (translation's `delta`/`angle`/`shape` at smaller location offsets).

Alternatives considered:
- **Per-mode effect lists** (`effect_sizes` keyed by mode): about 50 engine lines across enumeration, the Type I negative control, the loader's known keys and the config echo, byte-identical when absent. It would save 6,000 units here and 6,000 on the stage-effect arm. Not taken: it is a second engine feature in a run change, and the pilot declared it a non-goal. If the user prefers it, it replaces this decision without touching any other cell (cell identity does not depend on the configuration's effect list).
- **A third profile for translation at e = 1.0 only** (14 cells): saves 2,000 units net, because it adds a third anchor per column and two more Type I baseline cells, and it complicates the figure script's two-report contract. Rejected.

### D3. The pilot's seed family, so the slice reproduces the pilot

Both profiles keep `base_seed` 1028 and family `design-guidance`. Consequences, verified on the drafts:
- at every design column the slice evaluates the pilot's datasets at replicates 0–99 and 400 fresh ones; the observed statistics and selected ranks at those replicates must equal the pilot's, and the p-values differ only by the permutation count. The report states this as a reproduction check and never pools pilot and slice: the slice supersedes the pilot's rates at these columns;
- the magnitude and geometry profiles' anchors at a column share seeds and the RRPP seed, so their p-values are identical (the pilot verified 0 mismatches over 3,300 comparisons). The report reads the geometry profile's anchor; the magnitude profile's is the reproduction check;
- the Type I baseline cells draw fresh seeds (keyed on their own ids), giving the baseline column one more independent 500-replicate size reading.

Alternative: a new family. Independent of the pilot, but it discards a free, bit-exact provenance check and gains nothing: the pilot is not pooled either way.

### D4. Predeclared rules: anchor decision plus a column-scoped translation rule

Each slice profile declares:
- `anchor_decision` at α 0.05, `se_multiplier` 2, `delta`/`angle`/`shape` (bound 0.0695 at R = 500), as the pilot;
- a phase-scoped `type_i` entry (α 0.05, tolerance 2), as the pilot;
- a **column-scoped** `type_i` entry: `{"alpha": 0.05, "se_tolerance": 2.0, "bound": "upper", "trajectory_mode": "translation", "statistic": "angle"}`.

No power target, no specificity target, no gate, no report contract: the slice is guidance.

**Column-scoped Type I entries** are the one engine addition. `TypeIControlTarget` gains `trajectory_mode: str | None = None` and `statistic: str | None = None`. The loader accepts `trajectory_mode` only among the negative-control modes (`none`, `translation`), accepts `statistic` only beside a mode and only among `delta`/`angle`/`shape`, and names the entry in every error. `_evaluate_type_i` keeps its phase-scoped path untouched for entries without a mode and routes entries with a mode to a column-scoped evaluator that:
- buckets completed records into design columns the way `evaluate_anchor_decision` does (`record_design_point`; `power_primary` is the baseline column; `type_i_baseline` is never read);
- picks, per column, the named mode's cell at effect 0 (`none`) or at the largest effect present (`translation`);
- emits one `TargetEvaluation` per (column, statistic), using the existing two-sided / `_judge_upper` judgement on that cell's summary, with observations `design_point`, `is_baseline`, `cell_id`, `phase`, `trajectory_mode`, `effect_size`, `statistic`, rate, SE, available replicates and the bound fields; `met: None` with a "no available replicates" rationale when the column has no completed record of the control;
- names the evaluation `type_i_control[<mode>@<effect>, <statistic>, <column>]`.

Byte identity: the new dataclass fields are omitted from `_config_to_dict` and `dump_study_config` while `None`, so the configuration echo of every committed profile is unchanged, and an entry without a mode produces the same observations as before (the acceptance CSV gains the new observation columns only when a column-scoped entry exists).

Alternatives considered:
- **`anchor_decision.controls`** (named controls judged per point against the anchor bound, as `column_decision.controls`). Equivalent in outcome; rejected because the user and the pilot report predeclare it as a `type_i` rule, the acceptance report is where predeclared rules are evaluated, and the anchor decision's "other controls are descriptive" contract stays intact.
- **Reading translation off `anchor_decision.csv` by hand**, as the pilot did. Not a predeclared verdict.

### D5. Figure script: a grade label, nothing else

`scripts/design_guidance_figure.py --grade <label>` (default `pilot`) replaces the hard-coded word in both suptitles; the module docstring stops saying "pilot reports". A single stage count already renders one panel row; a test pins it. The committed pilot figures regenerate unchanged under the default.

### D6. Anchor reproduction by design column

`scripts/anchor_reproduction.py --by-column` pairs anchors on (column, replicate index), where the column is the record's `cell_metadata.design_point` for `power_design`, the baseline for `power_primary`, and `type_i_baseline` anchors are excluded (a different draw). It writes the same per-replicate columns plus the column coordinates and prints one summary line per column, listing a column with zero pairs when the other run lacks it; the verdict covers paired columns only. Without the flag the existing single-anchor path is untouched (it still refuses several anchor cells). This commits the pilot's session-only `anchor_check.py` logic, and its lesson: key by design column, never by `generator_seed`, which every column shares at a replicate within a seed family.

### D7. Stage-effect arm as two profiles with δ = 0.25 as the configuration baseline

`design_guidance_stage_effect_magnitude_study.json` and `..._geometry_study.json`: the slice profiles with the three per-omic δ set to 0.25 in the generator block and explicit design points over `generator.delta_methyl` / `delta_expr` / `delta_protein` ∈ {0.25, 0.5, 1.0}. Every column sits at n = 200, three stages, two blocks.

Why 0.25 is the baseline: the engine always enumerates the configuration baseline as the primary column, so a δ = 2.0 baseline would re-run the slice's n = 200 / two-block column (8,000 units, about 75 core-hours, bit-identical records). The slice column is the δ = 2.0 reference; the arm's columns are matched to it through the shared seed family. The Type I baseline phase (anchor and translation at δ = 0.25) is a fresh draw at the arm's baseline. The arm deviates from the "one shared baseline" convention of the pilot profiles, and its test says so.

Priority: committed and pinned with the slice, so it is reproducible whenever it runs; submitted only on the user's decision after the slice's numbers and timings are in. Expected cost about 230 core-hours (23,000 units; the pilot's δ columns ran with inflated ranks and cost more per unit than the core baseline).

### D8. Results and report

- Records under `results/design-guidance-slice-<date>/{core_magnitude,core_geometry}/` (`report/`, `PROVENANCE.txt`), `figures/` (the four figures and `design_guidance_power.csv`, `--grade paper`), `reproduction/` (the by-column CSVs for siblings and pilot); `stage_effect_{magnitude,geometry}/` if run. JSONL stays on the cluster, gitignored.
- `docs/reports/design-guidance-slice-<date>.md` follows the pilot report's structure: provenance; unit, failure and reproduction accounting; per-column operating characteristics ± SE with the design-guidance figures; stage signal and selected rank per column; size at every column (anchor decision and the translation `angle` rule from `acceptance_report.csv`); the guidance statement; the stage-effect arm (or why not run); limitations.
- `PROVENANCE.txt` per profile as the pilot's, naming the code revision as the version (the cluster venv's distribution metadata lags).

### D9. Execution

1. Scratchpad copies of all four profiles at 2 replicates × 19 permutations (never committed) through `run_study_shard.py` → `motco_study.py merge` → `report` → `design_guidance_figure.py --grade smoke` → `anchor_reproduction.py --by-column` (siblings, and against the pilot's merged records at matching columns where available). Any engine failure found is fixed here, with a test.
2. On `ing`, time one anchor and one orientation e = 1.0 unit at 999 permutations at each (n, blocks) column; fix the shard counts at ≤ 2.5 h per shard; present the plan (units, core-hours, shards, wall) to the user before any submission.
3. Push the branch; add a clean worktree on `ing` beside the main checkout (`git -C /home1/tgonza/MOTCO worktree add ../MOTCO-design-guidance-slice <branch>`, `PYTHONPATH` on its `src`, the main checkout's `.venv`); submit the two arrays with `scripts/motco_study_array.sbatch` (partition `512x1024`, 1 CPU, 4 GB, 6 h, BLAS pinned, `STUDY_N_JOBS` unset, `--error-policy record`).
4. Finish: merge, report, figures, both reproduction checks, `PROVENANCE.txt`; rsync `report/`, `PROVENANCE.txt`, `figures/`, `reproduction/`; write the report; update the roadmap, `AGENTS.md`, both READMEs; bump the version; run the gate.
5. The stage-effect arm repeats 2–4 on the user's go.

## Risks / Trade-offs

- **Translation at three effects costs 6,000 units** (D2) → stated in the proposal; the per-mode alternative is documented and swappable before submission without changing any other cell.
- **The pilot's per-unit overrun** (geometry +23 %, rank-inflated units at the largest columns) → the plan budgets 450 core-hours against the 380 estimate, times units at 999 permutations before sizing shards, and sizes shards at 2.5 h under a 6 h limit.
- **An anchor excess at paper grade.** The pilot's two excesses were at 0.10 against 0.0936; at R = 500 the bound is 0.0695 and chance predicts about 0.4 excesses over 18 anchor tests. → The report reads any excess against that expectation and the baseline column's independent Type I-phase draw, and the slice makes no gate claim; an excess at a column is reported as such beside its rates.
- **Translation `angle` near the bound.** The net-displacement studies measured translation `angle` at 0.058–0.068 against 0.0695; a column may exceed. → The column-scoped rule is predeclared and advisory; the report states the verdict per column and the guidance statement carries it ("a rejection beside a location difference warrants care").
- **Column bucketing drifts from the anchor decision's.** → The column-scoped evaluator reuses `record_design_point` and the same baseline rule; a test asserts the two outputs agree on the anchor at every column of a synthetic design-grid study.
- **Cluster state.** The main checkout on `ing` is a dirty rsync and destructive git operations there are blocked → a fresh worktree per branch, as for the pilot; nothing on the main checkout is touched.

## Open Questions

- **Whether to run the stage-effect arm.** Decided by the user after the slice's numbers and the measured unit cost at 999 permutations are in; it changes no spec or task of the slice.
- **The report's date and results directory name** follow the submission date; placeholders `<date>` in tasks resolve then.
