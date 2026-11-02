"""Covariate-adjusted orientation attribution: the primitive, the compact record and the analysis."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

import motco.stats.attribution as attribution_module
from motco.stats import AttributionError, analyze_orientation_attribution, attribution_frames, fit_plsda_model

GROUPS = ("A", "B")
STAGES = ("t0", "t1", "t2")
CELL_SIZES = ((6, 7, 8), (8, 6, 7))
P = 6


def _cell_structure() -> np.ndarray:
    """Group × stage cell means whose groups bend differently (2 × 3 × P)."""

    mu = np.zeros((2, 3, P))
    for g in range(2):
        for s in range(3):
            mu[g, s, 0] = float(s)
            mu[g, s, 1 + g] = 0.8 * s**2
            mu[g, s, 3] = 0.5 * s * (1 - 2 * g)
    return mu


def _confounded(*, noise: float = 0.3, seed: int = 0) -> dict[str, Any]:
    """Cell structure plus ``z · v`` with ``z``'s mean differing by cell, and a second covariate ``w``."""

    rng = np.random.default_rng(seed)
    mu = _cell_structure()
    v = np.array([0.0, 1.5, -1.0, 0.5, 2.0, 0.0])
    rows, labels, z, w = [], [], [], []
    for g in range(2):
        for s in range(3):
            for _ in range(CELL_SIZES[g][s]):
                z_i = 1.0 + 0.7 * s + 0.9 * g + rng.normal(scale=0.5)
                rows.append(mu[g, s] + z_i * v + rng.normal(scale=noise, size=P))
                labels.append((GROUPS[g], STAGES[s]))
                z.append(z_i)
                w.append(rng.normal())
    index = pd.Index([f"s{i}" for i in range(len(rows))])
    X = pd.DataFrame(np.asarray(rows), index=index, columns=[f"f{j}" for j in range(P)])
    metadata = pd.DataFrame(labels, index=index, columns=["group", "stage"])
    covariates = pd.DataFrame({"z": z, "w": w}, index=index)
    return {"X": X, "metadata": metadata, "covariates": covariates, "mu": mu, "v": v}


def _model(X: pd.DataFrame, metadata: pd.DataFrame):
    return fit_plsda_model(X, metadata["stage"], n_components=2)


def _run(data: dict[str, Any], model=None, **kwargs):
    kwargs.setdefault("groups", GROUPS)
    kwargs.setdefault("stages", STAGES)
    return analyze_orientation_attribution(
        data["X"], data["metadata"], model if model is not None else _model(data["X"], data["metadata"]), **kwargs
    )


def _blocks(result) -> dict[tuple[str, str], Any]:
    blocks: dict[tuple[str, str], Any] = {}
    for record in [*result.transitions, result.principal_orientation]:
        block_id = getattr(record, "transition_id", None) or record.component_id
        for component in ("observed", "pls_captured", "residual"):
            blocks[(block_id, component)] = getattr(record, component)
    return blocks


def _assert_same_contrasts(left, right, *, components=("observed", "pls_captured", "residual")) -> None:
    lhs, rhs = _blocks(left), _blocks(right)
    assert lhs.keys() == rhs.keys()
    for key in lhs:
        if key[1] not in components:
            continue
        a, b = lhs[key], rhs[key]
        np.testing.assert_allclose(a.path_lengths, b.path_lengths, rtol=1e-10, atol=1e-10, err_msg=str(key))
        for x, y in zip(a.group_transitions, b.group_transitions):
            np.testing.assert_allclose(x, y, rtol=1e-10, atol=1e-10, err_msg=str(key))
        assert (a.directional_contrast is None) == (b.directional_contrast is None), key
        if a.directional_contrast is not None:
            np.testing.assert_allclose(
                a.directional_contrast, b.directional_contrast, rtol=1e-10, atol=1e-10, err_msg=str(key)
            )


def _mean_frame(values: np.ndarray) -> pd.DataFrame:
    index = pd.MultiIndex.from_tuples([(g, s) for g in GROUPS for s in STAGES])
    return pd.DataFrame(values.reshape(-1, P), index=index, columns=[f"f{j}" for j in range(P)])


# ── Adjusted means ───────────────────────────────────────────────────────────


def test_known_covariate_effect_is_removed() -> None:
    data = _confounded(noise=0.0)
    model = _model(data["X"], data["metadata"])
    adjusted = _run(data, model, covariates=data["covariates"][["z"]])
    structure = _run(data, model, mean_table=_mean_frame(data["mu"]))
    unadjusted = _run(data, model)

    assert adjusted.config.mean_source == "covariate_adjusted"
    _assert_same_contrasts(adjusted, structure)
    observed = unadjusted.principal_orientation.observed.directional_contrast
    expected = structure.principal_orientation.observed.directional_contrast
    assert not np.allclose(observed, expected, atol=1e-3)
    for record, reference in zip(unadjusted.transitions, structure.transitions):
        assert not np.allclose(record.observed.directional_contrast, reference.observed.directional_contrast, atol=1e-3)


def test_adjusted_means_are_the_models_least_squares_means() -> None:
    data = _confounded()
    result = _run(data, covariates=data["covariates"])
    X = data["X"].to_numpy()
    cells = [(g, s) for g in GROUPS for s in STAGES]
    labels = list(zip(data["metadata"]["group"], data["metadata"]["stage"]))
    indicators = np.array([[float(label == cell) for cell in cells] for label in labels])
    C = data["covariates"].to_numpy()
    beta = np.linalg.lstsq(np.column_stack([indicators, C]), X, rcond=None)[0]
    reference = C.mean(axis=0)
    expected = beta[: len(cells)] + reference @ beta[len(cells) :]
    observed = result.means.set_index(["group", "stage"]).loc[cells].to_numpy()
    np.testing.assert_allclose(observed, expected, rtol=1e-10, atol=1e-10)


def test_contrasts_do_not_depend_on_the_reference_row() -> None:
    data = _confounded()
    model = _model(data["X"], data["metadata"])
    default = _run(data, model, covariates=data["covariates"])
    other = _run(data, model, covariates=data["covariates"], covariate_reference={"z": 3.0, "w": -1.0})
    _assert_same_contrasts(default, other)
    shift = other.means[[f"f{j}" for j in range(P)]].to_numpy() - default.means[[f"f{j}" for j in range(P)]].to_numpy()
    np.testing.assert_allclose(shift, np.tile(shift[0], (shift.shape[0], 1)), rtol=1e-10, atol=1e-10)
    assert np.abs(shift[0]).max() > 1e-3


def test_covariates_collinear_within_cells_are_rejected() -> None:
    data = _confounded()
    cell_value = data["metadata"]["stage"].map({"t0": 0.0, "t1": 1.0, "t2": 3.0}) + (
        data["metadata"]["group"] == "B"
    ).astype(float)
    covariates = data["covariates"].assign(z_cell=cell_value)
    with pytest.raises(AttributionError, match="z_cell"):
        _run(data, covariates=covariates)


def _bad_covariates(data: dict[str, Any], kind: str) -> dict[str, Any]:
    covariates = data["covariates"].copy()
    if kind == "rows":
        return {"covariates": covariates.iloc[:-1]}
    if kind == "index":
        return {"covariates": covariates.set_axis([f"x{i}" for i in range(len(covariates))])}
    if kind in ("nan", "inf"):
        covariates.iloc[3, 0] = np.nan if kind == "nan" else np.inf
        return {"covariates": covariates}
    if kind == "duplicated":
        return {"covariates": pd.concat([covariates[["z"]], covariates[["z"]]], axis=1)}
    if kind == "non-numeric":
        return {"covariates": covariates.assign(w=["a"] * len(covariates))}
    if kind == "empty":
        return {"covariates": covariates.iloc[:, :0]}
    if kind == "reference-columns":
        return {"covariates": covariates, "covariate_reference": {"z": 1.0, "age": 2.0}}
    if kind == "reference-length":
        return {"covariates": covariates, "covariate_reference": [1.0]}
    if kind == "reference-non-finite":
        return {"covariates": covariates, "covariate_reference": {"z": 1.0, "w": np.nan}}
    if kind == "reference-without-covariates":
        return {"covariate_reference": {"z": 1.0, "w": 0.0}}
    raise AssertionError(kind)


@pytest.mark.parametrize(
    ("kind", "message"),
    [
        ("rows", "rows"),
        ("index", "aligned"),
        ("nan", "missing or non-finite"),
        ("inf", "missing or non-finite"),
        ("duplicated", "unique"),
        ("non-numeric", "numeric"),
        ("empty", "at least one"),
        ("reference-columns", "reference"),
        ("reference-length", "reference"),
        ("reference-non-finite", "reference"),
        ("reference-without-covariates", "covariate_reference requires covariates"),
    ],
)
def test_invalid_covariate_input_is_rejected(kind: str, message: str) -> None:
    data = _confounded()
    with pytest.raises(AttributionError, match=message):
        _run(data, **_bad_covariates(data, kind))


def test_precomputed_means_and_covariates_are_exclusive() -> None:
    data = _confounded()
    with pytest.raises(AttributionError, match="only one mean source"):
        _run(data, covariates=data["covariates"], mean_table=_mean_frame(data["mu"]))
    with pytest.raises(AttributionError, match="only one mean source"):
        _run(data, covariates=data["covariates"], precomputed_means=_mean_frame(data["mu"]))


def test_adjusted_configuration_and_interpretation_are_recorded() -> None:
    data = _confounded()
    covariates = data["covariates"].rename(columns={"z": "age", "w": "pmi"})
    result = _run(data, covariates=covariates, bootstrap_replicates=50, bootstrap_seed=4, top_k=3)
    configuration = dict(attribution_frames(result)["configuration"].itertuples(index=False))
    assert configuration["mean_source"] == "covariate_adjusted"
    assert configuration["covariate_columns"] == ["age", "pmi"]
    assert configuration["covariate_reference"] == pytest.approx(covariates.mean().to_dict())
    assert list(configuration["covariate_reference"]) == ["age", "pmi"]
    assert configuration["bootstrap_rank_deficient_replicates"] == 0
    interpretation = dict(attribution_frames(result)["interpretation"].itertuples(index=False))
    assert "additive linear covariate model" in interpretation["adjustment_boundary"]
    assert "reference covariate row" in interpretation["adjustment_boundary"]
    assert "unmeasured confounding" in interpretation["adjustment_boundary"]


def test_no_covariates_adds_no_adjustment_rows() -> None:
    data = _confounded()
    frames = attribution_frames(_run(data, bootstrap_replicates=5))
    assert not {"covariate_columns", "covariate_reference", "bootstrap_rank_deficient_replicates"} & set(
        frames["configuration"]["setting"]
    )
    assert "adjustment_boundary" not in set(frames["interpretation"]["field"])


# ── The adjusted bootstrap ───────────────────────────────────────────────────


class _RecordingGenerator:
    """Delegates to a real generator and records every ``choice`` draw."""

    def __init__(self, generator: np.random.Generator, log: list[np.ndarray]) -> None:
        self._generator = generator
        self._log = log

    def choice(self, *args: Any, **kwargs: Any) -> np.ndarray:
        drawn = self._generator.choice(*args, **kwargs)
        self._log.append(np.asarray(drawn).copy())
        return drawn

    def __getattr__(self, name: str) -> Any:
        return getattr(self._generator, name)


def _record_draws(monkeypatch) -> list[np.ndarray]:
    log: list[np.ndarray] = []
    real = np.random.default_rng
    monkeypatch.setattr(
        attribution_module.np.random, "default_rng", lambda seed=None: _RecordingGenerator(real(seed), log)
    )
    return log


def test_adjusted_and_unadjusted_bootstraps_resample_the_same_rows(monkeypatch) -> None:
    data = _confounded()
    model = _model(data["X"], data["metadata"])
    log = _record_draws(monkeypatch)
    _run(data, model, bootstrap_replicates=7, bootstrap_seed=9)
    unadjusted = [draw.copy() for draw in log]
    log.clear()
    _run(data, model, covariates=data["covariates"], bootstrap_replicates=7, bootstrap_seed=9)
    assert len(unadjusted) == 7 * 6 and len(log) == len(unadjusted)
    for left, right in zip(unadjusted, log):
        np.testing.assert_array_equal(left, right)


def test_adjusted_bootstrap_is_reproducible() -> None:
    data = _confounded()
    model = _model(data["X"], data["metadata"])
    first = _run(data, model, covariates=data["covariates"], bootstrap_replicates=10, bootstrap_seed=3)
    second = _run(data, model, covariates=data["covariates"], bootstrap_replicates=10, bootstrap_seed=3)
    pd.testing.assert_frame_equal(first.bootstrap_summaries, second.bootstrap_summaries)
    pd.testing.assert_frame_equal(first.feature_effects, second.feature_effects)


def test_the_adjustment_is_refitted_in_every_replicate(monkeypatch) -> None:
    data = _confounded()
    model = _model(data["X"], data["metadata"])
    draws = _record_draws(monkeypatch)
    fits: list[np.ndarray] = []
    real_fit = attribution_module.fit_nuisance_coefficients

    def spy(X, design):
        fits.append(np.asarray(X, dtype=float).copy())
        return real_fit(X, design)

    monkeypatch.setattr(attribution_module, "fit_nuisance_coefficients", spy)
    replicates = 4
    result = _run(data, model, covariates=data["covariates"], bootstrap_replicates=replicates, bootstrap_seed=5)
    assert len(fits) == 1 + replicates  # the full sample, then one refit per replicate
    np.testing.assert_array_equal(fits[0], data["X"].to_numpy())
    n_cells = len(GROUPS) * len(STAGES)
    for replicate in range(replicates):
        rows = np.concatenate(draws[replicate * n_cells : (replicate + 1) * n_cells])
        np.testing.assert_array_equal(fits[1 + replicate], data["X"].to_numpy()[rows])
    assert result.config.rank_deficient_replicates == 0


def test_rank_deficient_resamples_are_invalid_and_counted() -> None:
    data = _confounded()
    metadata = data["metadata"]
    flag = np.zeros(len(metadata))
    flag[np.flatnonzero((metadata["group"] == "A") & (metadata["stage"] == "t1"))[2]] = 1.0
    covariates = data["covariates"].assign(**{"site[rare]": flag})
    replicates = 40
    result = _run(data, covariates=covariates, bootstrap_replicates=replicates, bootstrap_seed=1)
    deficient = result.config.rank_deficient_replicates
    assert deficient is not None
    assert 0 < deficient < replicates
    assert set(result.bootstrap_summaries["valid_replicates"]) <= set(range(replicates - deficient + 1))
    principal = result.bootstrap_summaries.query("transition_id == 'principal' and component == 'observed'")
    assert set(principal["valid_replicates"]) == {replicates - deficient}
    configuration = dict(attribution_frames(result)["configuration"].itertuples(index=False))
    assert configuration["bootstrap_rank_deficient_replicates"] == deficient


def test_exact_cell_structure_gives_full_sign_stability() -> None:
    data = _confounded(noise=0.0)
    result = _run(data, covariates=data["covariates"][["z"]], bootstrap_replicates=20, bootstrap_seed=2)
    observed = result.feature_effects.query("component == 'observed'")
    nonzero = observed[observed["effect_standardized"].abs() > 1e-8]
    assert len(nonzero) > 0
    assert (nonzero["sign_stability"] == 1.0).all()
    assert (nonzero["valid_bootstrap_replicates"] == 20).all()
