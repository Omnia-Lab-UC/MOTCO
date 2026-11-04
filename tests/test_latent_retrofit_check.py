"""``scripts/latent_retrofit_check.py`` on a small, cheap copy of a stage-span pilot profile."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from motco.stats.pls import fit_plsda_model

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
try:
    import latent_retrofit_check as retrofit
finally:
    sys.path.pop(0)

CONFIG_DIR = Path(__file__).resolve().parents[1] / "examples" / "trajectory_power_study"
PROFILE = CONFIG_DIR / "stage_span_pilot_surgery_geometry.json"


@pytest.fixture(scope="module")
def small_profile(tmp_path_factory) -> Path:
    raw = json.loads(PROFILE.read_text(encoding="utf-8"))
    raw["n_replicates"] = 2
    raw["evaluation"]["integration_params"].update({"n_repeats": 1, "max_components": 4})
    path = tmp_path_factory.mktemp("retrofit") / "profile.json"
    path.write_text(json.dumps(raw), encoding="utf-8")
    return path


@pytest.fixture(scope="module")
def smoke(small_profile, tmp_path_factory) -> tuple[pd.DataFrame, dict]:
    out_dir = tmp_path_factory.mktemp("retrofit-out")
    code = retrofit.main(
        [
            "--config",
            str(small_profile),
            "--column",
            json.dumps({"generator.n_samples": 80}),
            "--replicates",
            "1",
            "--permutations",
            "9",
            "--out-dir",
            str(out_dir),
        ]
    )
    assert code == 0
    frame = pd.read_csv(out_dir / retrofit.CSV_NAME)
    summary = json.loads((out_dir / retrofit.SUMMARY_NAME).read_text(encoding="utf-8"))
    return frame, summary


def test_outputs_carry_the_documented_columns(smoke) -> None:
    frame, summary = smoke
    assert list(frame.columns) == list(retrofit.CSV_COLUMNS)
    assert len(frame) == len(retrofit.EDIT_MODES) * len(retrofit.PREIMAGES)
    assert set(frame["trajectory_mode"]) == set(retrofit.EDIT_MODES)
    assert set(frame["preimage"]) == set(retrofit.PREIMAGES)
    assert summary["decision"] is None
    assert summary["settings"]["column"] == {"generator.n_samples": 80}
    assert {(cell["trajectory_mode"], cell["preimage"]) for cell in summary["cells"]} == {
        (mode, kind) for mode in retrofit.EDIT_MODES for kind in retrofit.PREIMAGES
    }


def test_frozen_basis_recovers_every_edit(smoke) -> None:
    frame, summary = smoke
    assert frame["frozen_recovery_error"].max() <= 1e-8
    assert summary["max_frozen_recovery_error"] <= 1e-8
    for statistic in ("delta", "angle", "shape"):
        np.testing.assert_allclose(frame[f"frozen_{statistic}"], frame[f"intended_{statistic}"], atol=1e-6)


def test_edits_and_footprints_have_their_shape(smoke) -> None:
    frame, _ = smoke
    rows = frame.set_index(["trajectory_mode", "preimage"])
    orientation = rows.loc[("orientation", "loadings")]
    assert orientation["edit_angle"] == pytest.approx(30.0)
    assert orientation["edit_delta"] == pytest.approx(0.0, abs=1e-8)
    magnitude = rows.loc[("magnitude", "minimum_norm")]
    assert magnitude["edit_angle"] == pytest.approx(0.0, abs=1e-5)
    assert rows.loc[("translation", "loadings"), "footprint_constant_share"] == pytest.approx(1.0)
    for mode in ("magnitude", "orientation", "shape"):
        assert rows.loc[(mode, "loadings"), "footprint_constant_share"] == pytest.approx(0.0, abs=1e-10)
    for column in ("footprint_stage_active_share", "footprint_methylation_share"):
        assert frame[column].between(0.0, 1.0).all()
    assert frame["footprint_effective_features"].between(1.0, frame["n_features"].iloc[0]).all()
    assert frame["cosine_population_plant"].between(-1.0, 1.0).all()


def test_both_preimages_score_to_the_edit() -> None:
    rng = np.random.default_rng(1)
    X = rng.normal(size=(45, 12))
    y = np.repeat(["0", "1", "2"], 15)
    X[y == "1", :3] += 1.5
    X[y == "2", 3:6] += 1.5
    model = fit_plsda_model(X, y, n_components=3)
    delta = rng.normal(size=(3, 3))
    rotations = np.asarray(model.x_rotations_)
    for kind in retrofit.PREIMAGES:
        np.testing.assert_allclose(retrofit.preimage(model, delta, kind) @ rotations, delta, atol=1e-10)
    loadings = retrofit.preimage(model, delta, "loadings")
    minimum = retrofit.preimage(model, delta, "minimum_norm")
    assert np.linalg.norm(minimum) <= np.linalg.norm(loadings) + 1e-12


def test_unknown_column_and_effect_are_refused(small_profile) -> None:
    with pytest.raises(retrofit.RetrofitCheckError, match="selects 0 anchor cells"):
        retrofit.anchor_cell(small_profile, {"generator.n_samples": 999})
    with pytest.raises(retrofit.RetrofitCheckError, match="MODE=VALUE"):
        retrofit.parse_effects(["rotation=1"])
