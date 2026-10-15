# Design

## Context

See proposal.md for motivation. These facts shape the approach:

- `SimulationEvaluationParams.orientation_estimator` is `str | None`, default `None`. `__post_init__`
  rejects unknown names and rewrites an explicit `"pc1"` to `None`; `resolved_orientation_estimator`
  returns `self.orientation_estimator or "pc1"`. `grid._SIGNATURE_OMIT_WHEN_NONE` drops the field from
  every hash payload while it is `None`, so a PC1 evaluation hashes exactly as it did before the field
  existed. `evaluate_semisynthetic_trajectory` writes `runtime_metadata["orientation_estimator"]` only when
  the field is not `None`. Every committed record was produced under PC1 and therefore carries no entry.
- The study loader (`study/config.py`) builds the baseline `SimulationEvaluationParams` from the raw
  `evaluation` block by field name and normalizes axis values through `_AXIS_VALUE_NORMALIZERS`, where
  `evaluation.orientation_estimator` maps `"pc1"` to `None`. That makes `null`, `"pc1"` and "key absent" one
  value — true only while the default is PC1. The nested `evaluation.integration_params.<key>` form already
  defines `null` as "key absent, so the default applies".
- The four net-displacement profiles declare explicit design points with
  `"evaluation.orientation_estimator": null` for their baseline column and omit the key from their
  `evaluation` block. The loader requires one point to equal the baseline evaluation; the matching point
  becomes the primary (`power`) column and the others become `power_design` columns with their own cell ids.
- Committed identity is pinned twice already: `tests/data/committed_config_identity.json` (16 configs:
  signatures and matched seeds, captured before `group_stage_sizes` existed) and
  `tests/data/historical_config_signatures.json` (the ten magnitude-pinned configs). The three Phase 6
  configs are covered only through the net-displacement profile tests' seed-family checks and the
  Phase 6 profile test.
- `SimulationEvaluationParams` is constructed directly, outside the loader, in: `showcase.py` (figure
  only), `realized_geometry_study.py` (committed `results/latent-spectrum-2026-09-02/`, checkpoint
  `angle`), `specificity.py` three sites (committed geometry-specificity, cross-talk and
  magnitude-construction diagnostics, `angle` rejection rates), `scripts/latent_rank_probe.py` (committed
  `results/latent-rank-probe-2026-09-03/`, orientation mode), `scripts/orientation_measurement_probe.py`
  (sets the estimator per arm; unaffected), and the two fixture capture modules
  `tests/orientation_probe_fixture_capture.py` and `tests/phase6_fixture_capture.py`, whose fixtures were
  captured under PC1.
- The stats primitives (`estimate_difference`, `RRPP`, `_RRPPWorker`, `diagnostics.geometry_from_means`,
  `calculate_realized_geometry`) take `orientation_estimator: str = "pc1"` and are pinned by the R-reference
  regression tests and `trajectory-orientation-invariance`. `motco de` wraps them and exposes no estimator
  flag.
- The precedent is `default-joint-magnitude-construction` (2026-10-09): pin the old value in the config
  files, capture a pre-flip fixture before any source edit, flip the default, prove the pins are
  behavior-neutral with a test.

## Goals / Non-Goals

**Goals:**

- Every harness caller that relies on the default measures `angle` with net displacement. Every committed
  config, diagnostic and fixture keeps its exact behavior through an explicit `"pc1"`.
- One test fails if any committed config's enumerated cell set, signatures or matched seeds change, and if
  any committed config stops naming its estimator.
- `null` in a study config means the same thing for every evaluation-namespace key: the evaluator's default.

**Non-Goals:**

- Flipping the primitives' default or adding an estimator flag to `motco de`. See D3.
- Removing the `None` representation of PC1. It is the hash-compatibility mechanism; see D1.
- Re-capturing any existing byte-identity fixture. The fixtures are PC1 by construction and stay valid
  once their capture definitions pin `"pc1"`.

## Decisions

### D1. The default is the string `"net_displacement"`; `None` stays the hash-omitted spelling of PC1

The dataclass default changes from `None` to `"net_displacement"`. Nothing else in the dataclass moves:
`"pc1"` still normalizes to `None`, `None` is still omitted from hash payloads, and
`resolved_orientation_estimator` still returns `"pc1"` for `None`. Consequences:

- A PC1 evaluation (explicit `"pc1"` or `None`) hashes and records exactly as before the field existed.
  That is what preserves every committed signature once the configs pin `"pc1"`.
- A default evaluation hashes with `"orientation_estimator": "net_displacement"` and records the estimator
  in `runtime_metadata`, exactly as an explicit net-displacement column does today. The four committed
  net-displacement columns are unchanged.

*Alternative rejected:* swapping the sentinel so that `None` means net displacement and is omitted, with
`"pc1"` hashed as a string. Every committed PC1 signature would change, and the committed net-displacement
columns' signatures would change too. No pin could save them.

*Alternative rejected:* keeping `None` as the default and making `resolved_orientation_estimator` return
`"net_displacement"`. The committed PC1 records would then be indistinguishable in their metadata from new
default records that measured a different statistic.

### D2. A config `null` resolves to the evaluator's default; `"pc1"` names PC1; every committed config pins `"pc1"`

The loader's normalizer for `evaluation.orientation_estimator` changes from `"pc1" → None` to:
`None` → the dataclass field's default (read from `dataclasses.fields`, so a later default change needs no
loader edit) and `"pc1"` → `None`. The dedup rule for explicit points follows: `null` equals
`"net_displacement"`, and `"pc1"` is its own value. The baseline `evaluation` block needs no special
handling: an absent key takes the dataclass default, and `"pc1"` is normalized by `__post_init__`.

Every committed config gains `"orientation_estimator": "pc1"` in its `evaluation` block, next to
`integration_method`. The four net-displacement profiles additionally rewrite their baseline point's `null`
to `"pc1"`. Both edits are required together: a profile whose baseline block omits the key would, after the
flip, have a net-displacement baseline evaluation, so its `null`/`"pc1"` point would no longer be the
primary column and the net-displacement point would silently become it, with different phases and cell
ids. D4's test catches exactly this.

*Alternative rejected:* keeping `null` ≡ `"pc1"` in the loader. It would leave the four profiles' points
untouched, but `null` would then mean "PC1" for this key and "the default" for every other evaluation key,
and a reader of a config would misread the baseline column as the default. The spec sentence that today
equates `null` with PC1 is rewritten rather than preserved.

*Alternative rejected:* a loader-side default keyed on a config date or version. Hidden behavior; a config
would no longer say what it runs. Same reasoning as the magnitude precedent's D1.

### D3. The primitives keep PC1; the harness is the production pipeline

`estimate_difference`, `RRPP`, `_RRPPWorker`, `geometry_from_means` and `calculate_realized_geometry` keep
`orientation_estimator="pc1"`. The harness passes `resolved_orientation_estimator` explicitly at every call,
so the harness default is independent of the primitive default.

Reasons: the roadmap item and the power comparison report both scope the flip to
`SimulationEvaluationParams`; the validated claim is about the production pipeline (stage-PLS latent space
with RRPP), which is the harness; `trajectory-orientation-invariance` defines the primitive `angle` as
principal-axis divergence and is pinned by the R-reference regression tests, which would otherwise need
`"pc1"` pins of their own; and `motco de` has no estimator flag, so its users would be unable to opt back.

Consequence, stated in the spec delta and the roadmap: real-data code that calls the primitives directly,
which the Phase 6 case study may do, must pass `orientation_estimator="net_displacement"` to measure what
the paper describes. The roadmap's case-study item carries that instruction. A later change may give
`motco de` an estimator flag; it is not needed to land this one.

### D4. Identity is checked against a fixture captured before any edit, covering all 19 configs

**Before** any source or config edit, a one-off capture writes
`tests/data/orientation_default_identity.json`: for each of the 19 committed configs, every enumerated
cell's `cell_id` → `parameter_signature` and the matched seeds of replicates 0–2, with a top-level note
naming the capture revision. After the pins and the flip, `tests/test_orientation_default_pins.py`
loads and enumerates each config and asserts the same cell-id set, the same signature per cell, the same
seeds, and that the config's baseline evaluation resolves to PC1 (`orientation_estimator is None`). A
second test requires every committed config to name `evaluation.orientation_estimator`, and every committed
config that declares the key as a design axis to name `"pc1"` rather than `null` for the PC1 column.

The existing `test_committed_config_identity.py` and `test_magnitude_default_pins.py` keep running and must
also pass; they guard the same property from older captures. The new fixture exists because neither covers
the three Phase 6 configs' signatures and the magnitude fixture covers only ten configs.

*Alternative rejected:* re-capturing `committed_config_identity.json` to add the three Phase 6 configs. Its
note records a specific pre-`group_stage_sizes` revision, and rewriting it would lose that provenance.

### D5. Diagnostics and fixture captures pin `"pc1"`; the showcase follows the default

- `realized_geometry_study.py`, the three `specificity.py` sites and `scripts/latent_rank_probe.py` pass
  `orientation_estimator="pc1"`. Their committed results were measured under PC1; a re-run then reproduces
  them. `specificity.py`'s functions do not grow a new parameter: the pin is at the constructor, as the
  committed diagnostics are historical. A future diagnostic that wants net displacement constructs its own
  parameters.
- `tests/orientation_probe_fixture_capture.py` and `tests/phase6_fixture_capture.py` pin `"pc1"` in their
  `EVALUATIONS`, so the byte-identity tests keep comparing like with like and the fixtures are not
  re-captured. Their docstrings say the fixtures are PC1 records.
- `showcase.py` passes nothing and so renders under the default. No showcase figure is committed.
- `scripts/orientation_measurement_probe.py` already names the estimator per arm and is untouched.

### D6. Reports and committed records are frozen; the roadmap carries the revised statement

No `docs/reports/` file, `results/` record or committed `report/` directory is rewritten. The revised
`angle` interpretability row is already written in the power comparison report's §7 and is marked there as
taking effect "once the default flips". `docs/roadmap.md` is the living document: its Phase 6 statement
replaces the pending-revision note with the §7 row, and the case-study item states that the cohort's
selected rank is the reported covariate and that primitives callers pass the estimator explicitly. Report
readers (`spectrum`, `report`, `summary`, the follow-up scripts) do not read `runtime_metadata`'s estimator
entry and need no change; the follow-up scripts match the candidate column by the normalized value
`"net_displacement"`, which D2 leaves unchanged, and the checks script's baseline coordinate `None` is
still what `"pc1"` normalizes to.

## Risks / Trade-offs

- [A config is missed or a design-grid baseline silently swaps columns] → D4's fixture covers all 19 configs
  and asserts the baseline resolves to PC1; the naming test fails on any config that omits the key or writes
  `null` for a PC1 column.
- [Existing tests assume the `None` default] → expected in `test_diagnostic_integration.py` (the default
  payload now contains the key), `test_study_measurement_axes.py` and `test_study_design_points_explicit.py`
  (baseline records now carry the estimator unless the test config pins `"pc1"`; `null` no longer equals
  `"pc1"`), `test_net_displacement_power_checks.py`, and any test study config that omits the key. Each is
  pinned to `"pc1"` where it tests the PC1 path or updated to the new default where it tests the default. No
  assertion is weakened.
- [Fixture captured on already-edited code] → the capture is task 1.1 and runs on the unedited tree; the
  fixture note records the revision.
- [The case study measures PC1 by calling the primitives] → D3's spec sentence and the roadmap item say to
  pass the estimator explicitly; the roadmap records that the case study waits for this change.
- [Readers take `None` in the dataclass for "no estimator chosen"] → the dataclass comment and `CLAUDE.md`
  say that `None` is the hash-omitted spelling of PC1 and that the default is net displacement.
- [The forced-rank guard or report readers key on the estimator] → verified absent (`grep` over
  `study/` finds the key only in `config.py`); D6.

## Migration Plan

Nothing committed breaks. A caller of `SimulationEvaluationParams` that wants PC1 passes
`orientation_estimator="pc1"`; a study config that wants PC1 writes `"pc1"` (a `null` now means the
default). A caller of the primitives is unaffected. To roll back, revert the dataclass default and the
loader normalizer; the pins are harmless under either default and can stay.
