"""The committed Phase 5 magnitude re-measurement profiles enumerate as declared.

Both profiles re-measure the magnitude mode under the ``joint`` construction at
the Phase 5 paper-grade design point. They derive from ``phase5_power_study.json``
so that the shared zero-effect anchor — and every group A baseline — is
byte-identical to the Phase 5 records at the same replicate index.
"""

from __future__ import annotations

import json
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from motco.simulations.grid import derive_replicate_seed
from motco.simulations.semisynthetic import generate_semisynthetic_trajectory, load_reference
from motco.simulations.study import enumerate_study, load_study_config
from motco.simulations.study.enumerate import SEED_FAMILY_KEY

PHASE5 = Path("examples/trajectory_power_study/phase5_power_study.json")
PILOT = Path("examples/trajectory_power_study/phase5_magnitude_pilot.json")
PAPER = Path("examples/trajectory_power_study/phase5_magnitude_remeasurement.json")
BRACKET = Path("results/magnitude-axis-bracket-2026-10-08/effect_axis_bracket.csv")
PILOT_GRID = [0.0, 0.0025, 0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 1.0]
PAPER_GRID = [0.0, 0.0025, 0.005, 0.01, 0.02, 0.25, 1.0]


def _anchor(grid):
    anchors = [cell for cell in grid.cells if cell.metadata.get("zero_effect_anchor")]
    assert len(anchors) == 1
    return anchors[0]


@pytest.mark.parametrize("path", [PILOT, PAPER])
def test_profiles_derive_from_the_phase5_paper_grade_profile(path: Path) -> None:
    raw = json.loads(path.read_text(encoding="utf-8"))
    phase5 = json.loads(PHASE5.read_text(encoding="utf-8"))

    assert raw["metadata"]["derives_from"] == str(PHASE5)
    generator = dict(raw["generator"])
    assert generator.pop("magnitude_kind") == "joint"
    phase5_generator = dict(phase5["generator"])
    assert phase5_generator.pop("magnitude_kind") == "all"
    assert generator == phase5_generator, "generator must be copied from Phase 5 except magnitude_kind"
    assert raw["evaluation"]["integration_params"] == phase5["evaluation"]["integration_params"]
    assert raw["evaluation"]["integration_method"] == "pls"
    assert raw["base_seed"] == phase5["base_seed"]
    assert raw["matched_seeds"] == phase5["matched_seeds"]
    assert raw["trajectory_modes"] == ["magnitude"]
    assert raw["effect_sizes"] == (PILOT_GRID if path == PILOT else PAPER_GRID)
    assert raw["attribution"] == {"enabled": False}
    assert "design_grid" not in raw and raw["axes"] == {}
    assert "surgery_censoring" not in raw["generator"]
    assert "intersim" not in json.dumps(raw).lower()

    bracket = raw["metadata"]["bracket"]
    assert bracket["source"] == str(BRACKET)
    assert bracket["effect_grid"] == raw["effect_sizes"]
    assert set(bracket["realized_joint_delta_population_standardized"]) == {str(e) for e in raw["effect_sizes"]}

    config = load_study_config(path)
    assert config.generator.magnitude_kind == "joint"
    assert config.generator.surgery_censoring == "error"
    assert config.trajectory_modes == ("magnitude",)
    assert config.matched_seeds.primary_family == "phase5-primary"
    assert config.matched_seeds.shared_zero_effect_anchor
    assert [(t.alpha, t.se_tolerance) for t in config.acceptance.type_i] == [(0.05, 2.0)]
    assert not config.attribution.enabled
    # At least two nonzero effects strictly below 0.25, plus 0 and 1.
    below = [e for e in config.effect_sizes if 0 < e < 0.25]
    assert len(below) >= 2 and 0.0 in config.effect_sizes and 1.0 in config.effect_sizes


def test_bracket_csv_backs_the_metadata_values() -> None:
    import pandas as pd

    frame = pd.read_csv(BRACKET)
    assert set(frame["construction"]) == {"all", "joint"}
    counts = frame.groupby("construction")["effect_size"].nunique()
    assert counts["all"] == counts["joint"] >= 401
    for e in PILOT_GRID:
        assert ((frame["effect_size"] - e).abs() < 1e-9).sum() == 2, e
    joint = frame[frame["construction"] == "joint"]
    anchor = frame[(frame["construction"] == "all") & (frame["effect_size"] == 0.0)].iloc[0]
    assert joint["joint_angle"].abs().max() < 1e-4
    assert joint["joint_shape"].abs().max() < 1e-12
    assert abs(anchor["joint_angle"]) < 1e-4
    recorded = json.loads(PILOT.read_text(encoding="utf-8"))["metadata"]["bracket"]
    for key, value in recorded["realized_joint_delta_population_standardized"].items():
        row = joint[(joint["effect_size"] - float(key)).abs() < 1e-9].iloc[0]
        assert row["joint_delta"] == pytest.approx(value, abs=5e-4)


def test_pilot_enumerates_fifty_by_199_with_the_gate_disabled() -> None:
    config = load_study_config(PILOT)
    assert config.n_replicates == 50
    assert config.evaluation.permutations == 199
    assert not config.acceptance.gate.enabled
    assert config.acceptance.power == () and config.acceptance.specificity == ()
    assert config.report_contract is None

    grid = enumerate_study(config)
    phases = Counter(cell.phase for cell in grid.cells)
    assert phases == {"type_i_baseline": 2, "power_primary": 1 + 8}
    assert sum(cell.n_replicates for cell in grid.cells) == 11 * 50
    anchor = _anchor(grid)
    assert anchor.metadata["resolves_modes"] == ["magnitude"]
    nonzero = sorted(
        cell.metadata["effect_size"]
        for cell in grid.cells
        if cell.phase == "power_primary" and cell.metadata["trajectory_mode"] == "magnitude"
    )
    assert nonzero == PILOT_GRID[1:]
    assert all(cell.evaluation_params.permutations == 199 for cell in grid.cells)
    assert all(cell.n_replicates == 50 for cell in grid.cells)


def test_paper_grade_enumerates_500_by_999_with_the_reduced_gate() -> None:
    config = load_study_config(PAPER)
    assert config.n_replicates == 500
    assert config.evaluation.permutations == 999
    assert config.evaluation.n_jobs == 1

    contract = config.report_contract
    assert contract is not None
    assert contract.driver_component == "observed"
    assert contract.cross_replicate_driver_agreement == "descriptive"
    assert contract.forbids_n_jobs_override
    phase5 = json.loads(PHASE5.read_text(encoding="utf-8"))
    assert json.loads(PAPER.read_text(encoding="utf-8"))["report_contract"] == phase5["report_contract"]

    gate = config.acceptance.gate
    assert gate.enabled
    assert gate.control_modes == ("none",)
    assert gate.min_power_at_top == 0.8
    roles = {(rule.trajectory_mode, rule.statistic): rule.role for rule in gate.rules}
    assert roles == {
        ("magnitude", "delta"): "mandatory_power",
        ("magnitude", "angle"): "mandatory_control",
        ("magnitude", "shape"): "mandatory_control",
    }
    acceptance = config.acceptance
    assert {(t.trajectory_mode, t.statistic, t.min_power_at_top) for t in acceptance.power} == {
        ("magnitude", "delta", 0.8)
    }
    assert {(t.trajectory_mode, t.statistic, t.alpha, t.se_tolerance) for t in acceptance.specificity} == {
        ("magnitude", "angle", 0.05, 2.0),
        ("magnitude", "shape", 0.05, 2.0),
    }
    assert acceptance.design_point is None and acceptance.rank_decision is None

    grid = enumerate_study(config)
    phases = Counter(cell.phase for cell in grid.cells)
    assert phases == {"type_i_baseline": 2, "power_primary": 1 + 6}
    assert sum(cell.n_replicates for cell in grid.cells) == 9 * 500
    assert all(cell.evaluation_params.permutations == 999 for cell in grid.cells)
    assert all(cell.evaluation_params.n_jobs == 1 for cell in grid.cells)
    assert all(not cell.evaluation_params.attribution.enabled for cell in grid.cells)


def test_paper_grade_grid_is_a_pilot_measured_subset() -> None:
    """Every paper-grade effect was measured in the pilot (design D5)."""

    pilot = json.loads(PILOT.read_text(encoding="utf-8"))
    paper = json.loads(PAPER.read_text(encoding="utf-8"))
    assert set(paper["effect_sizes"]) <= set(pilot["effect_sizes"])
    assert paper["metadata"]["pilot"] == str(PILOT)
    assert "pilot_rise" in paper["metadata"]["bracket"]


def test_anchor_shares_phase5_anchor_seeds_at_every_replicate_index() -> None:
    phase5 = _anchor(enumerate_study(load_study_config(PHASE5)))
    for path in (PILOT, PAPER):
        config = load_study_config(path)
        grid = enumerate_study(config)
        anchor = _anchor(grid)
        primary = [cell for cell in grid.cells if cell.phase == "power_primary"]
        assert {cell.metadata[SEED_FAMILY_KEY] for cell in primary} == {"phase5-primary"}
        for index in (0, 1, config.n_replicates - 1):
            assert derive_replicate_seed(anchor, index) == derive_replicate_seed(phase5, index)
            for cell in primary:
                assert derive_replicate_seed(cell, index) == derive_replicate_seed(anchor, index)


@pytest.mark.slow
def test_anchor_dataset_at_replicate_zero_equals_the_phase5_anchor() -> None:
    reference = load_reference()
    phase5 = _anchor(enumerate_study(load_study_config(PHASE5)))
    pilot = _anchor(enumerate_study(load_study_config(PILOT)))
    from dataclasses import replace

    seed = derive_replicate_seed(phase5, 0)
    assert derive_replicate_seed(pilot, 0) == seed
    expected = generate_semisynthetic_trajectory(
        replace(phase5.generator_params, seed=seed), reference=reference
    )
    observed = generate_semisynthetic_trajectory(
        replace(pilot.generator_params, seed=seed), reference=reference
    )
    assert pilot.generator_params.magnitude_kind == "joint"
    for name in ("methylation", "expression", "proteomics", "metadata"):
        a, b = getattr(expected, name), getattr(observed, name)
        assert list(a.index) == list(b.index) and list(a.columns) == list(b.columns)
        np.testing.assert_array_equal(a.to_numpy(), b.to_numpy())
    assert expected.truth["deltas"] == observed.truth["deltas"]
    assert expected.truth["indicator_counts"] == observed.truth["indicator_counts"]


def test_phase5_artefacts_are_untouched() -> None:
    """Adding the profiles edits nothing Phase 5 committed."""

    root = Path(__file__).resolve().parents[1]
    paths = [
        "results/phase5-2026-09-10",
        "examples/trajectory_power_study/phase5_report_template.md",
        "docs/reports/phase5-paper-grade-2026-09-10.md",
        "docs/reports/phase5-exit-review-2026-09-10.md",
    ]
    try:
        diff = subprocess.run(
            ["git", "diff", "--stat", "HEAD", "--", *paths],
            capture_output=True, text=True, check=True, cwd=root,
        )
        committed = subprocess.run(
            ["git", "show", f"HEAD:{PHASE5}"], capture_output=True, text=True, check=True, cwd=root,
        )
    except (OSError, subprocess.CalledProcessError):  # pragma: no cover - no git in the sandbox
        pytest.skip("git is not available")
    assert diff.stdout.strip() == "", diff.stdout

    # The Phase 5 config may differ only by the historical ``magnitude_kind: "all"``
    # pin (``default-joint-magnitude-construction``) and ``orientation_estimator:
    # "pc1"`` pin (``default-net-displacement-orientation``), which
    # ``tests/test_magnitude_default_pins.py`` and
    # ``tests/test_orientation_default_pins.py`` prove are signature-neutral.
    def unpinned(text: str) -> dict:
        raw = json.loads(text)
        assert raw["generator"].pop("magnitude_kind", "all") == "all"
        assert raw["evaluation"].pop("orientation_estimator", "pc1") == "pc1"
        return raw

    assert unpinned((root / PHASE5).read_text(encoding="utf-8")) == unpinned(committed.stdout)
