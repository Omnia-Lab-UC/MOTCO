"""Smoke test for ``scripts/orientation_checkpoint_evidence.py`` on a tiny synthetic ``merged.jsonl``."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "orientation_checkpoint_evidence.py"
_spec = importlib.util.spec_from_file_location("orientation_checkpoint_evidence", _SCRIPT)
assert _spec is not None and _spec.loader is not None
evidence = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(evidence)

BASELINE = {"generator.group_stage_sizes": [[11, 10, 28], [9, 10, 12]]}
OTHER = {"generator.group_stage_sizes": [[10, 9, 25], [9, 9, 12]]}


def _geometry(angle: float) -> dict:
    return {"angle": angle, "delta": 0.0, "shape": 0.0}


def _record(mode: str, effect: float, *, observed: float, latent: float, rank: int, p: float, design=BASELINE):
    population = 0.0 if mode == "none" else float(effect * 90.0)
    return {
        "phase": "power_primary" if design is BASELINE else "power_design",
        "status": "completed",
        "cell_metadata": {"trajectory_mode": mode, "effect_size": effect, "design_point": design},
        "truth_metadata": {"n_stages": 3},
        "integration_metadata": {"selected_lv": rank},
        "p_values": {"angle": p},
        "realized_geometry": {
            "checkpoints": {
                "population_native": {"methylation": _geometry(population)},
                "population_standardized": {"joint": _geometry(population)},
                "observed_standardized": {"joint": _geometry(observed)},
                "pls_latent": {"joint": _geometry(latent)},
            }
        },
    }


@pytest.fixture()
def merged(tmp_path: Path) -> Path:
    records = []
    # Anchor: signed angles with a flip tail (one replicate near 150°).
    for observed in (10.0, 20.0, 30.0, 150.0):
        records.append(_record("none", 0.0, observed=observed, latent=observed / 2, rank=3, p=0.5))
    # Orientation e = 1: observed near 90°, latent collapsed at low rank.
    for observed, latent, rank, p in ((88.0, 30.0, 3, 0.2), (92.0, 32.0, 2, 0.3), (90.0, 80.0, 5, 0.01)):
        records.append(_record("orientation", 1.0, observed=observed, latent=latent, rank=rank, p=p))
    # A non-baseline column that must be ignored.
    records.append(_record("orientation", 1.0, observed=0.0, latent=0.0, rank=3, p=0.9, design=OTHER))
    path = tmp_path / "merged.jsonl"
    path.write_text("".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
    return path


def test_evidence_script_renders_both_tables(merged: Path, tmp_path: Path) -> None:
    out = tmp_path / "evidence.md"
    payload = tmp_path / "evidence.json"
    assert evidence.main(["--merged", str(merged), "--out", str(out), "--json", str(payload)]) == 0
    text = out.read_text(encoding="utf-8")
    assert "## Checkpoint table" in text and "## Flip-tail table" in text
    tables = json.loads(payload.read_text(encoding="utf-8"))
    rows = {(r["mode"], r["rank"]): r for r in tables["checkpoints"]}
    assert rows[("orientation", "all")]["n"] == 3
    assert rows[("orientation", "<=3")]["latent"] == pytest.approx(31.0)
    assert rows[("orientation", ">3")]["angle_rejection_rate"] == pytest.approx(1.0)
    assert rows[("orientation", "all")]["population"] == pytest.approx(90.0)
    separability = tables["flip_tail"]["separability"][0]
    # The signed anchor q95 sits in the flip tail, so no e = 1 replicate clears
    # it; folding removes the tail and every replicate clears.
    assert separability["observed_signed"] == 0.0
    assert separability["observed_folded"] == 1.0
