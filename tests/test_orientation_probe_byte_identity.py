"""Default evaluation paths reproduce the pre-orientation-probe output.

``tests/data/orientation_probe_prechange_fixture.json`` was captured by
``tests/orientation_probe_fixture_capture.py`` before the net-displacement
estimator and the diagnostic ``pca`` / ``pls_within_group`` methods existed.
Integer, string and structural content must match exactly; floats are compared
at a tight relative tolerance so the fixture survives BLAS differences.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.orientation_probe_fixture_capture import EVALUATIONS, evaluation_summary
from tests.test_phase6_byte_identity import _assert_close

FIXTURE = json.loads(
    (Path(__file__).resolve().parent / "data" / "orientation_probe_prechange_fixture.json").read_text(
        encoding="utf-8"
    )
)


@pytest.mark.parametrize("name", sorted(EVALUATIONS))
def test_default_evaluation_matches_pre_change(name: str) -> None:
    _assert_close(json.loads(json.dumps(evaluation_summary(EVALUATIONS[name]))), FIXTURE["evaluations"][name])
