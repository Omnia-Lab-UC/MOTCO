"""Byte identity of orientation attribution without covariates.

``tests/data/attribution_no_covariate_identity.json`` was captured with the
0.13.0 attribution (before covariate-adjusted means existed) by running this
module as a script::

    uv run python -m tests.test_attribution_identity

It holds two payloads:

- ``primitive``: every table :func:`~motco.stats.attribution.attribution_frames`
  returns for a fixed synthetic dataset and fitted model, with original scales,
  feature groups and 25 bootstrap replicates;
- ``analysis``: the attribution block of an ``analyze_trajectories`` record
  (without ``runtime.attribution_seconds``, which varies by run) and the parsed
  attribution CSVs, on the ``tests/test_analysis_identity.py`` dataset.

Missing values are held as ``None``. Structure, strings and integers must match
exactly; floats are compared at a tight relative tolerance so the fixture
survives BLAS differences (the same rule as ``tests/test_phase6_byte_identity.py``).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from motco.analysis import analyze_trajectories, write_analysis_outputs
from motco.simulations.evaluation import AttributionDiagnosticSettings
from motco.stats import analyze_orientation_attribution, attribution_frames, fit_plsda_model
from tests.test_analysis_identity import _identity_dataset, _identity_params
from tests.test_phase6_byte_identity import _assert_close

FIXTURE = Path(__file__).parent / "data" / "attribution_no_covariate_identity.json"
GROUPS = ("A", "B")
STAGES = ("t0", "t1", "t2")
N_FEATURES = 8
BOOTSTRAP: dict[str, Any] = {"bootstrap_replicates": 25, "bootstrap_seed": 11, "top_k": 4}


def _none_for_missing(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _none_for_missing(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_none_for_missing(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _frame_payload(frame: pd.DataFrame) -> dict[str, Any]:
    # A JSON round trip turns numpy scalars into Python ones; NaN survives it and is then held as None.
    return _none_for_missing(json.loads(frame.to_json(orient="split", double_precision=15)))


def primitive_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    """A fixed two-group, three-stage dataset whose groups bend differently."""

    rng = np.random.default_rng(20261102)
    rows, labels = [], []
    for group_index, group in enumerate(GROUPS):
        for stage_index, stage in enumerate(STAGES):
            centre = np.zeros(N_FEATURES)
            centre[0] = float(stage_index)
            centre[1 + group_index] = 0.8 * stage_index**2
            centre[3] = 0.5 * stage_index * (1 - 2 * group_index)
            for _ in range(6):
                rows.append(centre + rng.normal(scale=0.4, size=N_FEATURES))
                labels.append((group, stage))
    index = pd.Index([f"s{i}" for i in range(len(rows))])
    X = pd.DataFrame(np.asarray(rows), index=index, columns=[f"f{j}" for j in range(N_FEATURES)])
    metadata = pd.DataFrame(labels, index=index, columns=["group", "stage"])
    return X, metadata


def primitive_payload() -> dict[str, Any]:
    X, metadata = primitive_inputs()
    model = fit_plsda_model(X, metadata["stage"], n_components=2)
    result = analyze_orientation_attribution(
        X,
        metadata,
        model,
        groups=GROUPS,
        stages=STAGES,
        original_scales=np.linspace(0.5, 2.0, N_FEATURES),
        feature_groups={f"f{j}": f"module_{j % 3}" for j in range(N_FEATURES)},
        **BOOTSTRAP,
    )
    return {name: _frame_payload(frame) for name, frame in attribution_frames(result).items()}


def _csv_payload(path: Path) -> Any:
    # A table without rows (aggregate effects without feature groups) is written as an empty file.
    try:
        return _none_for_missing(pd.read_csv(path).to_dict(orient="split"))
    except pd.errors.EmptyDataError:
        return {"raw": path.read_text(encoding="utf-8")}


def analysis_payload(out_dir: Path) -> dict[str, Any]:
    settings = AttributionDiagnosticSettings(enabled=True, **BOOTSTRAP)
    result = analyze_trajectories(_identity_dataset(), _identity_params(permutations=0, attribution=settings))
    paths = write_analysis_outputs(result, out_dir, overwrite=True)
    block = json.loads(paths["result"].read_text(encoding="utf-8"))["attribution"]
    block["runtime"].pop("attribution_seconds")
    tables = {
        key.split("/", 1)[1]: _csv_payload(path)
        for key, path in sorted(paths.items())
        if key.startswith("attribution/")
    }
    return json.loads(json.dumps({"record": block, "csv": tables}))


def capture_identity(out_dir: Path) -> dict[str, Any]:
    """Both payloads, JSON round-tripped exactly as the fixture stores them."""

    return json.loads(json.dumps({"primitive": primitive_payload(), "analysis": analysis_payload(out_dir)}))


def test_primitive_attribution_without_covariates_is_byte_identical_to_0_13_0() -> None:
    expected = json.loads(FIXTURE.read_text(encoding="utf-8"))["primitive"]
    observed = json.loads(json.dumps(primitive_payload()))
    _assert_close(observed, expected, "primitive")


def test_analysis_attribution_without_covariates_is_byte_identical_to_0_13_0(tmp_path: Path) -> None:
    expected = json.loads(FIXTURE.read_text(encoding="utf-8"))["analysis"]
    observed = analysis_payload(tmp_path)
    assert "adjustment" not in observed["record"]
    _assert_close(observed, expected, "analysis")


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as scratch:
        payload = capture_identity(Path(scratch))
    FIXTURE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {FIXTURE}")
