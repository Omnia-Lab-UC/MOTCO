"""``scripts/adni_anchor_check.py``: contract, design points, units, arms, checks, merge, decision.

Every analysis runs on a synthetic stand-in export written by this module in
the ``adni_sexdiff`` export's layout (manifest, metadata, covariates, ``p180``
and ``nmr``), at small CV knobs and few permutations. No real data are read.
The decision tests run on synthetic tables.
"""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
import math
import shlex
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _script(name: str):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    try:
        return __import__(name)
    finally:
        sys.path.pop(0)


aac = _script("adni_anchor_check")
wsp = _script("within_stage_permutation_check")
pc, psc, pse, plsc, lsc, grade = aac.pc, aac.psc, aac.pse, aac.plsc, aac.lsc, aac.grade

REPO = Path(__file__).resolve().parents[1]
PERMUTATIONS = 9
SMALL_CV = {"cv1_splits": 3, "cv2_splits": 4, "n_repeats": 1, "max_components": 6, "random_state": 1203}
STAND_IN_POINTS = (40, 80)
#: Female / Male × CN / MCI / AD, rows shuffled in the files.
STAND_IN_TABLE = [[30, 36, 18], [28, 44, 20]]
STAND_IN_FEATURES = {"p180": 12, "nmr": 15}
#: Participants at the rare ``apoe4`` level, so subsamples of 40 often miss it.
RARE_APOE4 = 4
ADNI_TABLE = [[237, 310, 125], [223, 448, 158]]


# ── The stand-in export ──────────────────────────────────────────────────────


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _participant_id(index: int) -> str:
    return f"PX{index:04d}Q"


def write_export(directory: Path, *, table=STAND_IN_TABLE, features=STAND_IN_FEATURES, seed: int = 7) -> Path:
    """A stand-in export in the ``adni_sexdiff`` layout, with a stage signal and a phase offset."""

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    groups, stages = ("Female", "Male"), ("CN", "MCI", "AD")
    rows = [(group, stage) for g, group in enumerate(groups) for s, stage in enumerate(stages)
            for _ in range(int(table[g][s]))]
    order = rng.permutation(len(rows))
    rows = [rows[index] for index in order]
    n = len(rows)
    ids = [_participant_id(1000 + index) for index in range(n)]
    metadata = pd.DataFrame({"rid": ids, "sex": [g for g, _ in rows], "stage": [s for _, s in rows]})
    apoe4 = rng.choice([0, 1], size=n, p=[0.6, 0.4])
    apoe4[rng.choice(n, size=RARE_APOE4, replace=False)] = 2
    phase = rng.choice(["ADNI1", "ADNI2GO"], size=n)
    covariates = pd.DataFrame({
        "rid": ids,
        "age": np.round(rng.normal(73.0, 6.0, size=n), 1),
        "education": rng.integers(8, 21, size=n),
        "apoe4": apoe4,
        "phase": phase,
        "statin": rng.integers(0, 2, size=n),
    })
    stage_position = np.array([stages.index(s) for _, s in rows], dtype=float)
    offset = (phase == "ADNI2GO").astype(float)
    blocks = {}
    for name, width in features.items():
        direction = rng.normal(size=width)
        second = rng.normal(size=width)
        values = (rng.normal(size=(n, width)) + 0.45 * np.outer(stage_position, direction)
                  + 0.30 * np.outer(stage_position == 1, second) + 0.5 * np.outer(offset, rng.normal(size=width)))
        blocks[name] = pd.DataFrame(values, index=pd.Index(ids, name="rid"),
                                    columns=[f"{name}_{index}" for index in range(width)])
    metadata.to_csv(directory / "metadata.csv", index=False)
    covariates.to_csv(directory / "covariates.csv", index=False)
    for name, frame in blocks.items():
        frame.to_csv(directory / f"{name}.csv")
    write_manifest(directory, table=table, features=features)
    return directory


def write_manifest(directory: Path, *, table=STAND_IN_TABLE, features=STAND_IN_FEATURES, **overrides) -> None:
    directory = Path(directory)
    outputs = {
        name: {"bytes": (directory / name).stat().st_size, "participant_level": True,
               "sha256": _sha256(directory / name)}
        for name in ("metadata.csv", "covariates.csv", "p180.csv", "nmr.csv")
    }
    manifest = {
        "step": "export",
        "settings": {"motco_transform": "identity", "stage_order": ["CN", "MCI", "AD"]},
        "outputs": outputs,
        "counts": {
            "participants": int(sum(sum(row) for row in table)),
            "sex_stage": {group: {stage: int(table[g][s]) for s, stage in enumerate(("CN", "MCI", "AD"))}
                          for g, group in enumerate(("Female", "Male"))},
            "features": {name: {"kept": int(width)} for name, width in features.items()},
        },
    }
    for key, value in overrides.items():
        manifest[key] = value
    (directory / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")


def _snapshot(directory: Path) -> dict[str, tuple[int, str]]:
    return {path.name: (path.stat().st_mtime_ns, _sha256(path)) for path in sorted(Path(directory).iterdir())}


def _copy_export(source: Path, target: Path) -> Path:
    shutil.copytree(source, target)
    return target


def _rewrite(directory: Path, name: str, frame: pd.DataFrame, *, index: bool = False, update_manifest: bool = True,
             **manifest) -> None:
    """Write one file and, by default, refresh its sha256 in the manifest so only the intended check fires."""

    frame.to_csv(directory / name, index=index)
    if update_manifest:
        stored = json.loads((directory / "manifest.json").read_text())
        stored["outputs"][name]["sha256"] = _sha256(directory / name)
        stored.update(manifest)
        (directory / "manifest.json").write_text(json.dumps(stored, indent=1))


def _stand_in_contract(**changes) -> aac.CohortContract:
    return dataclasses.replace(aac.ADNI_CONTRACT, name="stand_in", fidelity_source="the stand-in's own fit", **changes)


@pytest.fixture(scope="module")
def export_dir(tmp_path_factory) -> Path:
    return write_export(tmp_path_factory.mktemp("export") / "motco_input")


@pytest.fixture(scope="module")
def parent_dir(tmp_path_factory) -> Path:
    directory = tmp_path_factory.mktemp("parent")
    for name in ("rate_table.csv", "cell_table.csv"):
        shutil.copy(aac.PARENT_DIR / name, directory / name)
    return directory


@pytest.fixture(scope="module")
def fidelity(export_dir) -> dict:
    """The stand-in's own cohort fit under a within-stage label draw, label-free (as ADNI's 8 and 0.599514)."""

    cohort = aac.load_cohort(export_dir, _stand_in_contract(fidelity_rank=2, fidelity_auroc=0.5))
    settings = aac.build_settings(cohort, replicates=2, permutations=PERMUTATIONS, point_sizes=STAND_IN_POINTS,
                                  fixed_ranks=(2, 4), cv=SMALL_CV, parent_dir=aac.PARENT_DIR)
    built = aac.build_unit(cohort, settings, aac.COHORT_POINT, 0)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    latent = aac._analyse(dataset, aac.analysis_params(settings, 1, permutations=0)).record["latent"]
    return {"rank": int(latent["selected_lv"]), "auroc": round(float(latent["cv_mean_auroc"]), 6)}


@pytest.fixture(scope="module")
def contract(fidelity) -> aac.CohortContract:
    return _stand_in_contract(fidelity_rank=fidelity["rank"], fidelity_auroc=fidelity["auroc"])


@pytest.fixture(scope="module")
def cohort(export_dir, contract):
    return aac.load_cohort(export_dir, contract)


def _fixed_ranks(contract) -> tuple[int, int]:
    """As ADNI's (8, 19): the cohort's own rank and one above it."""

    return int(contract.fidelity_rank), int(contract.fidelity_rank) + 2


def _settings(cohort, parent_dir, **changes) -> dict:
    kwargs = {"replicates": 2, "permutations": PERMUTATIONS, "point_sizes": STAND_IN_POINTS,
              "fixed_ranks": _fixed_ranks(cohort.contract), "cv": SMALL_CV, "parent_dir": parent_dir}
    kwargs.update(changes)
    return aac.build_settings(cohort, **kwargs)


@pytest.fixture(scope="module")
def settings(cohort, parent_dir) -> dict:
    return _settings(cohort, parent_dir)


@pytest.fixture(scope="module")
def subsampled_record(cohort, settings) -> dict:
    return aac.evaluate_unit(aac.Unit("n40", 0), settings, cohort)


@pytest.fixture(scope="module")
def cohort_record(cohort, settings) -> dict:
    return aac.evaluate_unit(aac.Unit(aac.COHORT_POINT, 0), settings, cohort)


@pytest.fixture(scope="module")
def prepared(tmp_path_factory, cohort, settings) -> Path:
    """An output directory with a passed pre-flight."""

    out_dir = tmp_path_factory.mktemp("prepared")
    result = aac.preflight(out_dir, settings, cohort)
    assert result["passed"], result
    return out_dir


def _counting_cv(monkeypatch) -> list:
    import motco.simulations.evaluation as evaluation_module

    calls: list = []
    original = evaluation_module.plsda_doubleCV

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(evaluation_module, "plsda_doubleCV", counting)
    return calls


# ── 1.2 The contract ─────────────────────────────────────────────────────────


def test_adni_contract_is_the_case_study_configuration() -> None:
    contract = aac.ADNI_CONTRACT
    assert (contract.sample_col, contract.group_col, contract.stage_col) == ("rid", "sex", "stage")
    assert contract.groups == ("Female", "Male") and contract.stage_order == ("CN", "MCI", "AD")
    assert contract.block_names == ["p180", "nmr"]
    assert contract.covariate_names == ["age", "education", "apoe4", "phase", "statin"]
    assert contract.categorical_covariates == ["apoe4", "phase", "statin"]
    assert contract.levels == {"apoe4": ("0", "1", "2"), "phase": ("ADNI1", "ADNI2GO"), "statin": ("0", "1")}
    assert (contract.fidelity_rank, contract.fidelity_auroc, contract.fidelity_decimals) == (8, 0.599514, 6)
    assert contract.files_read == ["manifest.json", "metadata.csv", "covariates.csv", "p180.csv", "nmr.csv"]


def test_a_valid_stand_in_loads_in_the_metadata_order(export_dir, contract) -> None:
    cohort = aac.load_cohort(export_dir, contract)
    metadata = pd.read_csv(export_dir / "metadata.csv", dtype=str)
    assert cohort.metadata["rid"].tolist() == metadata["rid"].tolist()
    assert list(cohort.metadata.columns) == ["rid", "sex", "stage", "age", "education", "apoe4", "phase", "statin"]
    for name, frame in cohort.blocks.items():
        assert frame.index.tolist() == metadata["rid"].tolist()
        assert frame.shape[1] == STAND_IN_FEATURES[name]
    assert cohort.table() == STAND_IN_TABLE
    assert sorted(cohort.files) == sorted(contract.files_read)
    assert cohort.files["p180.csv"] == _sha256(export_dir / "p180.csv")


def test_loading_writes_nothing_under_the_export(export_dir, contract) -> None:
    before = _snapshot(export_dir)
    aac.load_cohort(export_dir, contract)
    assert _snapshot(export_dir) == before


def _metadata(directory: Path) -> pd.DataFrame:
    return pd.read_csv(directory / "metadata.csv", dtype=str)


def _covariates(directory: Path) -> pd.DataFrame:
    return pd.read_csv(directory / "covariates.csv", dtype={"rid": str})


def _block(directory: Path, name: str = "p180") -> pd.DataFrame:
    frame = pd.read_csv(directory / f"{name}.csv", index_col=0)
    frame.index = frame.index.astype(str)
    return frame


def _missing_manifest(directory: Path) -> None:
    (directory / "manifest.json").unlink()


def _missing_file(directory: Path) -> None:
    (directory / "nmr.csv").unlink()


def _hash_differs(directory: Path) -> None:
    frame = _block(directory)
    frame.iloc[0, 0] += 1.0
    _rewrite(directory, "p180.csv", frame, index=True, update_manifest=False)


def _transform(directory: Path) -> None:
    stored = json.loads((directory / "manifest.json").read_text())
    stored["settings"]["motco_transform"] = "beta_to_mvalue"
    (directory / "manifest.json").write_text(json.dumps(stored))


def _repeated(directory: Path) -> None:
    frame = _metadata(directory)
    frame.loc[1, "rid"] = frame.loc[0, "rid"]
    _rewrite(directory, "metadata.csv", frame)


def _third_group(directory: Path) -> None:
    frame = _metadata(directory)
    frame.loc[0, "sex"] = "Other"
    _rewrite(directory, "metadata.csv", frame)


def _unknown_stage(directory: Path) -> None:
    frame = _metadata(directory)
    frame.loc[0, "stage"] = "SMC"
    _rewrite(directory, "metadata.csv", frame)


def _empty_cell(directory: Path) -> None:
    frame = _metadata(directory)
    frame.loc[(frame["sex"] == "Female") & (frame["stage"] == "AD"), "stage"] = "MCI"
    _rewrite(directory, "metadata.csv", frame)


def _participant_count(directory: Path) -> None:
    stored = json.loads((directory / "manifest.json").read_text())
    stored["counts"]["participants"] += 1
    (directory / "manifest.json").write_text(json.dumps(stored))


def _cell_table(directory: Path) -> None:
    stored = json.loads((directory / "manifest.json").read_text())
    stored["counts"]["sex_stage"]["Female"]["CN"] += 1
    stored["counts"]["sex_stage"]["Male"]["CN"] -= 1
    (directory / "manifest.json").write_text(json.dumps(stored))


def _block_participants(directory: Path) -> None:
    frame = _block(directory)
    frame.index = [*frame.index[:-1], "PX9999Q"]
    frame.index.name = "rid"
    _rewrite(directory, "p180.csv", frame, index=True)


def _block_order(directory: Path) -> None:
    frame = _block(directory)
    _rewrite(directory, "p180.csv", frame.iloc[::-1], index=True)


def _feature_count(directory: Path) -> None:
    stored = json.loads((directory / "manifest.json").read_text())
    stored["counts"]["features"]["nmr"]["kept"] += 1
    (directory / "manifest.json").write_text(json.dumps(stored))


def _non_finite(directory: Path) -> None:
    frame = _block(directory, "nmr")
    frame.iloc[3, 2] = np.inf
    _rewrite(directory, "nmr.csv", frame, index=True)


def _missing_covariate(directory: Path) -> None:
    frame = _covariates(directory)
    frame.loc[2, "education"] = np.nan
    _rewrite(directory, "covariates.csv", frame)


def _undeclared_level(directory: Path) -> None:
    frame = _covariates(directory)
    frame.loc[2, "phase"] = "ADNI3"
    _rewrite(directory, "covariates.csv", frame)


VIOLATIONS = [
    ("missing_manifest", _missing_manifest, pc.ManifestError, "manifest.json"),
    ("missing_file", _missing_file, pc.ManifestError, "nmr.csv"),
    ("hash_differs", _hash_differs, pc.ManifestError, "p180.csv"),
    ("transform", _transform, pc.ManifestError, "motco_transform"),
    ("repeated_participant", _repeated, pc.DonorMismatchError, "metadata.csv"),
    ("third_group", _third_group, pc.GroupCountError, "Other"),
    ("unknown_stage", _unknown_stage, pc.UnknownStageError, "SMC"),
    ("empty_cell", _empty_cell, pc.EmptyCellError, "Female/AD"),
    ("participant_count", _participant_count, pc.ManifestError, "participants"),
    ("cell_table", _cell_table, pc.ManifestError, "cell table"),
    ("block_participants", _block_participants, pc.DonorMismatchError, "p180.csv"),
    ("block_order", _block_order, pc.DonorMismatchError, "same order"),
    ("feature_count", _feature_count, pc.ManifestError, "nmr.csv"),
    ("non_finite", _non_finite, pc.BlockValueError, "nmr.csv"),
    ("missing_covariate", _missing_covariate, pc.MissingCovariateError, "education"),
    ("undeclared_level", _undeclared_level, pc.CovariateValueError, "ADNI3"),
]


@pytest.mark.parametrize(("name", "mutate", "error", "named"), VIOLATIONS, ids=[v[0] for v in VIOLATIONS])
def test_each_contract_violation_is_named(tmp_path, export_dir, contract, name, mutate, error, named) -> None:
    directory = _copy_export(export_dir, tmp_path / "export")
    mutate(directory)
    with pytest.raises(error, match=None) as raised:
        aac.load_cohort(directory, contract)
    assert type(raised.value) is error, (name, type(raised.value))
    assert named in str(raised.value)
    assert aac.contract_exit(raised.value).code.startswith(f"export contract violation ({error.violation})")


# ── 1.3 Design points ────────────────────────────────────────────────────────


def test_adni_point_tables() -> None:
    assert aac.point_table(ADNI_TABLE, 80) == [[13, 16, 7], [12, 24, 8]]
    assert aac.point_table(ADNI_TABLE, 200) == [[31, 41, 17], [30, 60, 21]]
    assert aac.point_table(ADNI_TABLE, 500) == [[79, 103, 42], [74, 149, 53]]
    points = aac.design_points(ADNI_TABLE)
    assert [point["name"] for point in points] == ["n80", "n200", "n500", "cohort"]
    assert [sum(map(sum, point["table"])) for point in points] == [80, 200, 500, 1501]
    assert points[-1]["table"] == ADNI_TABLE and not points[-1]["subsampled"]


def test_point_tables_sum_to_n_and_break_ties_in_table_order() -> None:
    for n in range(1, 60):
        assert sum(map(sum, aac.point_table(STAND_IN_TABLE, n))) == n
    assert aac.point_table([[1, 1], [1, 1]], 2) == [[1, 1], [0, 0]]
    assert aac.point_table([[1, 1], [1, 1]], 3) == [[1, 1], [1, 0]]
    with pytest.raises(ValueError):
        aac.design_points(ADNI_TABLE, (200, 80))


# ── 1.4 Units: seeds, subsample, labels ──────────────────────────────────────


def test_a_subsample_holds_each_participant_once_in_its_cell(cohort, settings) -> None:
    for replicate in range(4):
        built = aac.build_unit(cohort, settings, "n40", replicate)
        rows = built["rows"]
        assert len(set(rows.tolist())) == rows.size == 40
        assert np.all(np.diff(rows) > 0)
        real = cohort.metadata.iloc[rows]
        assert aac.cell_table(real["sex"], real["stage"], cohort.contract) == aac.point_spec(settings, "n40")["table"]
        assert built["table"] == aac.point_spec(settings, "n40")["table"]


def test_units_are_reproducible_and_points_share_no_seed(cohort, settings) -> None:
    first, again = aac.build_unit(cohort, settings, "n40", 1), aac.build_unit(cohort, settings, "n40", 1)
    other = aac.build_unit(cohort, settings, "n40", 2)
    assert np.array_equal(first["rows"], again["rows"])
    assert first["metadata"].equals(again["metadata"])
    assert not np.array_equal(first["rows"], other["rows"])
    seeds = {point: aac.unit_seeds(settings, point, 1) for point in aac.point_names(settings)}
    for purpose in ("labels", "rrpp", "folds"):
        values = [seeds[point][purpose] for point in seeds]
        assert len(set(values)) == len(values)
    assert seeds["n40"]["subsample"] != seeds["n80"]["subsample"]
    assert "subsample" not in seeds[aac.COHORT_POINT]
    assert seeds["n40"]["labels"] == aac._seed(2039, "n40", 1, "labels", 0)


def test_the_label_draw_keeps_the_table_and_the_rows(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, "n80", 0)
    real = cohort.metadata.iloc[built["rows"]].reset_index(drop=True)
    permuted = built["metadata"]
    assert not np.array_equal(permuted["sex"].to_numpy(), real["sex"].to_numpy())
    assert built["table"] == aac.point_spec(settings, "n80")["table"]
    for stage in cohort.contract.stage_order:
        rows = real["stage"] == stage
        assert sorted(permuted.loc[rows, "sex"]) == sorted(real.loc[rows, "sex"])
    others = [column for column in real.columns if column != "sex"]
    assert permuted[others].equals(real[others])
    for name, frame in built["blocks"].items():
        assert frame.index.tolist() == permuted["rid"].tolist()
        assert frame.equals(cohort.blocks[name].loc[permuted["rid"]])


def test_the_cohort_point_holds_everyone_and_draws_no_subsample(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, aac.COHORT_POINT, 0)
    assert np.array_equal(built["rows"], np.arange(len(cohort.metadata)))
    assert built["subsample_redraws"] == 0
    assert built["metadata"]["rid"].tolist() == cohort.metadata["rid"].tolist()


def test_an_absent_covariate_level_is_redrawn_and_counted(cohort, settings) -> None:
    redraws = [aac.build_unit(cohort, settings, "n40", replicate)["subsample_redraws"] for replicate in range(12)]
    assert any(redraws), redraws
    for replicate, count in enumerate(redraws):
        built = aac.build_unit(cohort, settings, "n40", replicate)
        assert set(built["metadata"]["apoe4"].astype(str)) == {"0", "1", "2"}
        assert count == built["subsample_redraws"]


def test_a_label_draw_equal_to_the_real_labels_is_redrawn(monkeypatch, cohort, settings) -> None:
    original = pc.permute_within_stage
    calls: list = []

    def first_identity(labels, stages, order, rng):
        calls.append(1)
        if len(calls) == 1:
            return np.array(labels, dtype=object)
        return original(labels, stages, order, rng)

    monkeypatch.setattr(pc, "permute_within_stage", first_identity)
    built = aac.build_unit(cohort, settings, "n40", 3)
    assert built["label_redraws"] == 1 and len(calls) == 2


def test_two_thousand_units_at_the_committed_replicates() -> None:
    settings = {"points": aac.design_points(ADNI_TABLE), "replicates": aac.REPLICATES}
    units = aac.enumerate_units(settings)
    assert len(units) == 2000 and len({aac.unit_key(*unit) for unit in units}) == 2000
    assert units[:4] == [aac.Unit("n80", 0), aac.Unit("n200", 0), aac.Unit("n500", 0), aac.Unit("cohort", 0)]


# ── 2.1, 2.2 Arms ────────────────────────────────────────────────────────────


def test_arm_layout() -> None:
    specs = aac.arm_specs()
    assert list(specs) == [
        "pls_cv", "pls_cv__within_stage", "predictive_orthonormal", "predictive_orthonormal__within_stage",
        "pls_rank8", "pls_rank8__within_stage", "predictive_orthonormal_rank8",
        "predictive_orthonormal_rank8__within_stage", "pls_rank19", "pls_rank19__within_stage",
        "predictive_orthonormal_rank19", "predictive_orthonormal_rank19__within_stage", "pls_cv__uncorrected",
    ]
    settings = {"points": aac.design_points(ADNI_TABLE), "arm_specs": specs, "cohort_rank": 8}
    assert len(aac.point_arms(settings, "n80")) == 13
    cohort_arms = aac.point_arms(settings, "cohort")
    assert len(cohort_arms) == 9 and not any("rank8" in arm for arm in cohort_arms)
    assert specs["pls_cv__uncorrected"]["residual_correction"] is None
    assert all(spec["residual_correction"] == "df" for arm, spec in specs.items() if arm != "pls_cv__uncorrected")


def test_arm_parameters_differ_only_in_the_named_settings(settings) -> None:
    production = dataclasses.asdict(aac.arm_params("pls_cv", settings, 11, subsampled=True, production_rank=3))
    expected = {
        "pls_cv": {},
        "pls_cv__within_stage": {"forced_components": 3, "permutation_strata": "stage"},
        "predictive_orthonormal": {"forced_components": 3, "latent_coordinates": "predictive_orthonormal"},
        "predictive_orthonormal__within_stage": {"forced_components": 3, "latent_coordinates":
                                                 "predictive_orthonormal", "permutation_strata": "stage"},
        "pls_cv__uncorrected": {"forced_components": 3, "residual_correction": None},
    }
    rank = settings["fixed_ranks"][1]
    expected[f"pls_rank{rank}"] = {"forced_components": rank}
    expected[f"predictive_orthonormal_rank{rank}__within_stage"] = {
        "forced_components": rank, "latent_coordinates": "predictive_orthonormal", "permutation_strata": "stage"}
    for arm, changes in expected.items():
        params = dataclasses.asdict(aac.arm_params(arm, settings, 11, subsampled=True, production_rank=3))
        differing = {key: params[key] for key in params if params[key] != production[key]}
        assert differing == changes, arm
    assert production["residual_correction"] == "df" and production["seed"] == 11
    assert production["min_components"] == 2 and production["orientation_estimator"] == "net_displacement"
    cohort = dataclasses.asdict(aac.arm_params("pls_cv", settings, 11, subsampled=False, production_rank=3))
    assert cohort["forced_components"] == 3


def test_a_subsampled_unit_makes_one_cross_validation(monkeypatch, cohort, settings) -> None:
    calls = _counting_cv(monkeypatch)
    record = aac.evaluate_unit(aac.Unit("n40", 1), settings, cohort)
    assert record["status"] == "completed", record.get("inconsistent_checks") or record.get("error_message")
    assert len(calls) == 1
    assert len(record["arms"]) == 13
    rank = record["production"]["rank"]
    for arm, entry in record["arms"].items():
        spec = settings["arm_specs"][arm]
        assert entry["status"] == "completed", arm
        assert entry["rank"] == (rank if spec["fit"] == "cv" else spec["fixed_rank"]), arm
    assert record["arms"]["pls_cv"]["fit"] == "cross_validated"
    assert all(entry["fit"] == "forced_at_fixed_rank" for arm, entry in record["arms"].items() if "rank" in arm)


def test_a_cohort_unit_makes_no_cross_validation(monkeypatch, cohort, settings) -> None:
    calls = _counting_cv(monkeypatch)
    record = aac.evaluate_unit(aac.Unit(aac.COHORT_POINT, 1), settings, cohort)
    assert record["status"] == "completed", record.get("inconsistent_checks") or record.get("error_message")
    assert calls == []
    assert len(record["arms"]) == 9
    rank = int(settings["cohort_rank"])
    assert record["production"]["rank"] == rank
    for arm in aac.OPERATING_ARMS:
        assert record["arms"][arm]["rank"] == rank and record["arms"][arm]["fit"] == "forced_at_cohort_rank"
    assert not any(f"rank{rank}" in arm for arm in record["arms"])


def test_the_control_arm_observes_production_exactly(subsampled_record, cohort_record) -> None:
    for record in (subsampled_record, cohort_record):
        control, production = record["arms"]["pls_cv__uncorrected"], record["arms"]["pls_cv"]
        assert control["pair_statistics"] == production["pair_statistics"]
        assert control["residual_correction"] is None and production["residual_correction"] == "df"


def test_a_failing_arm_leaves_the_others(monkeypatch, cohort, settings) -> None:
    original = aac._analyse

    def failing(dataset, params):
        if params.latent_coordinates == "predictive_orthonormal" and params.permutation_strata == "stage":
            raise RuntimeError("planted failure")
        return original(dataset, params)

    monkeypatch.setattr(aac, "_analyse", failing)
    record = aac.evaluate_unit(aac.Unit("n40", 0), settings, cohort)
    failed = [arm for arm, entry in record["arms"].items() if entry["status"] == "failed"]
    assert failed and all(arm.startswith("predictive_orthonormal") and arm.endswith("__within_stage") for arm in failed)
    assert record["arms"]["predictive_orthonormal__within_stage"]["error_message"] == "planted failure"
    assert sum(entry["status"] == "completed" for entry in record["arms"].values()) == 13 - len(failed)


def test_a_failing_unit_is_recorded_and_the_shard_goes_on(monkeypatch, tmp_path, cohort, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    original = aac.build_unit

    def failing(cohort_, settings_, point, replicate):
        if point == "n40" and replicate == 0:
            raise RuntimeError("planted unit failure")
        return original(cohort_, settings_, point, replicate)

    monkeypatch.setattr(aac, "build_unit", failing)
    assert aac.run_shard(out_dir, settings, cohort, 0, 3, limit=2) == 2
    records = aac.read_records([lsc.shard_path(out_dir, 0, 3)])
    statuses = {aac.record_key(record): record["status"] for record in records}
    assert statuses["n40/0"] == "failed"
    assert any(status != "failed" for key, status in statuses.items() if key != "n40/0")


# ── 2.3 Checks ───────────────────────────────────────────────────────────────


def test_every_check_passes_on_a_stand_in_unit(subsampled_record, cohort_record) -> None:
    for record in (subsampled_record, cohort_record):
        assert record["inconsistent_checks"] == []
        for name in aac.CONSISTENCY_CHECKS:
            check = record["consistency"][name]
            assert check["applicable"] and check["passed"], (name, check)


def test_perturbed_models_statistics_and_tables_are_flagged(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, "n40", 0)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    result = aac._analyse(dataset, aac.analysis_params(settings, 5, forced=3))
    model = result.latent_artifacts.model
    perturbed = copy.deepcopy(model)
    perturbed.y_loadings_ = perturbed.y_loadings_ + 1e-6
    assert aac.check_same_fit({"a": model, "b": copy.deepcopy(model)}, ["a", "b"])["passed"]
    flagged = aac.check_same_fit({"a": model, "b": perturbed}, ["a", "b"])
    assert not flagged["passed"] and set(flagged["failing_arms"]) == {"b"}
    entry = {"status": "completed", "pair_statistics": {"delta": 1.0, "angle": 20.0, "shape": 0.1}}
    arms = {"pls_cv": entry, "pls_cv__within_stage": copy.deepcopy(entry), "pls_cv__uncorrected": copy.deepcopy(entry)}
    assert aac.check_same_statistics(arms)["passed"]
    arms["pls_cv__uncorrected"]["pair_statistics"]["angle"] += 1e-12
    flagged = aac.check_same_statistics(arms)
    assert not flagged["passed"] and set(flagged["failing_arms"]) == {"pls_cv__uncorrected"}
    assert aac.check_cell_table([[1, 2], [3, 4]], {"labels": [[1, 2], [3, 4]]})["passed"]
    assert not aac.check_cell_table([[1, 2], [3, 4]], {"labels": [[1, 2], [4, 3]]})["passed"]
    assert not aac.check_rank({"x": {"status": "completed", "rank": 3}}, {"x": 4})["passed"]
    predictive = {"status": "completed", "rank": 5, "measured_dimension": 3}
    assert not aac.check_measured_dimension({"predictive_orthonormal": predictive}, settings["arm_specs"], 3)["passed"]


def test_restated_checks_agree_with_the_pilots(cohort, settings) -> None:
    built = aac.build_unit(cohort, settings, "n40", 0)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    results, arms = {}, {}
    for arm in aac.OPERATING_ARMS:
        params = aac.arm_params(arm, settings, 5, subsampled=False, production_rank=3)
        results[arm] = aac._analyse(dataset, params)
        arms[arm] = aac.arm_record(arm, results[arm], 0.0, settings["arm_specs"][arm], "forced")
    model = results["pls_cv"].latent_artifacts.model
    pilot = psc.consistency_checks(
        model=model, scores=results["pls_cv"].latent_scores.to_numpy(dtype=float), n_stages=3,
        candidates={"predictive_orthonormal": results["predictive_orthonormal"].latent_scores.to_numpy(dtype=float)},
        arms={"pls_cv": arms["pls_cv"], "predictive_orthonormal": arms["predictive_orthonormal"]}, with_capture=False,
        models={"predictive_orthonormal": results["predictive_orthonormal"].latent_artifacts.model},
    )
    ours = aac.check_same_fit({arm: results[arm].latent_artifacts.model for arm in ("pls_cv",
                                                                                     "predictive_orthonormal")},
                              ["pls_cv", "predictive_orthonormal"])
    assert ours["deviation"] == pilot["same_fit"]["deviation"] and ours["passed"] == pilot["same_fit"]["passed"]
    assert aac.check_same_statistics(arms)["deviation"] == wsp.same_statistics(arms)["deviation"]


def test_an_inconsistent_unit_leaves_every_rate(settings) -> None:
    records = [_synthetic_record(settings, "n40", replicate, p=0.01) for replicate in range(3)]
    records[1]["status"] = "inconsistent"
    rows = {(row["arm"], row["statistic"]): row for row in aac.rate_table(records, settings) if row["point"] == "n40"}
    assert rows[("pls_cv", "delta")]["count"] == 2 and rows[("pls_cv", "delta")]["inconsistent"] == 1


# ── 2.4 Fit readings ─────────────────────────────────────────────────────────


def test_stage_r2_is_read_in_declared_stage_order(cohort, settings) -> None:
    from motco import analysis

    built = aac.build_unit(cohort, settings, "n80", 0)
    dataset = aac.make_dataset(built["blocks"], built["metadata"], cohort.contract)
    params = aac.analysis_params(settings, 5, forced=4)
    result = aac._analyse(dataset, params)
    inputs = analysis._trajectory_inputs(dataset, params)
    readings, dispersion = aac.fit_readings(inputs=inputs, results={"pls_cv": result}, settings=settings,
                                            fits=["cv"], folds_seed=17)
    model, features = result.latent_artifacts.model, result.latent_artifacts.features.to_numpy(dtype=float)
    stages = built["metadata"]["stage"].to_numpy()
    assert sorted(set(stages)) != list(cohort.contract.stage_order)  # the labels sort out of declared order
    onehot = (stages[:, None] == np.asarray(cohort.contract.stage_order)[None, :]).astype(float)
    fitted = model.predict(features)
    direct = 1.0 - np.sum((onehot - fitted) ** 2) / np.sum((onehot - onehot.mean(axis=0)) ** 2)
    assert math.isclose(readings["cv"]["in_sample_r2"], direct, rel_tol=1e-12)
    assert readings["cv"]["rank"] == 4 and readings["cv"]["out_of_fold_share"] is not None
    assert set(dispersion["scores"]["cv"]) == {"stage_0", "stage_1", "stage_2", "max_min_ratio"}


def test_folds_hold_everyone_out_once_and_reproduce(subsampled_record, settings) -> None:
    strata = np.repeat(["0", "1", "2"], [13, 17, 10])
    first = pse.stratified_folds(strata, ["0", "1", "2"], 99, folds=settings["folds"])
    assert np.array_equal(first, pse.stratified_folds(strata, ["0", "1", "2"], 99, folds=settings["folds"]))
    assert sorted(set(first.tolist())) == [0, 1, 2, 3] and first.size == strata.size
    assert subsampled_record["folds"] == {"seed": subsampled_record["seeds"]["folds"], "count": 4}
    fits = subsampled_record["fit_readings"]
    assert set(fits) == {"cv", *(f"rank{rank}" for rank in settings["fixed_ranks"])}
    assert all(entry["out_of_fold_share"] is not None for entry in fits.values())


# ── 2.5 Records ──────────────────────────────────────────────────────────────


def _walk(value, n: int, path: str = "") -> list[str]:
    found = []
    if isinstance(value, dict):
        for key, item in value.items():
            found += _walk(item, n, f"{path}.{key}")
    elif isinstance(value, list):
        if len(value) == n:
            found.append(path)
        for index, item in enumerate(value):
            found += _walk(item, n, f"{path}[{index}]")
    return found


def test_a_record_holds_no_participant_data(subsampled_record, cohort_record, cohort) -> None:
    ids = cohort.metadata["rid"].tolist()
    for record, n in ((subsampled_record, 40), (cohort_record, len(ids))):
        text = json.dumps(record)
        assert not any(identifier in text for identifier in ids)
        assert _walk(record, n) == []
        assert "rows" not in record and "labels" not in record


# ── 3.1 Pre-flight ───────────────────────────────────────────────────────────


def test_the_preflight_passes_on_the_stand_in(prepared, settings) -> None:
    record = json.loads((prepared / aac.PREFLIGHT_NAME).read_text())
    assert record["passed"] and record["settings_hash"] == lsc.settings_hash(settings)
    assert [unit["unit"] for unit in record["units"]] == ["cohort/0", "cohort/1"]
    for unit in record["units"]:
        assert unit["fidelity"]["passed"] and unit["forced_rank"] == settings["cohort_rank"]
        assert unit["statistics_differing"] == [] and unit["p_values_differing"] == []


@pytest.mark.parametrize("change", ["rank", "auroc"])
def test_a_preflight_on_another_fit_fails_and_the_run_refuses(tmp_path, export_dir, contract, parent_dir,
                                                              change) -> None:
    wrong = dataclasses.replace(contract, **(
        {"fidelity_rank": contract.fidelity_rank + 1} if change == "rank"
        else {"fidelity_auroc": round(contract.fidelity_auroc + 0.01, 6)}))
    cohort = aac.load_cohort(export_dir, wrong)
    settings = _settings(cohort, parent_dir)
    result = aac.preflight(tmp_path, settings, cohort)
    assert not result["passed"]
    assert not all(unit["fidelity"]["passed"] for unit in result["units"])
    with pytest.raises(SystemExit, match="pre-flight failed"):
        aac.run_shard(tmp_path, settings, cohort, 0, 1)


def test_a_forced_call_that_differs_fails_the_preflight(monkeypatch, tmp_path, cohort, settings) -> None:
    original = aac._analyse

    def drifting(dataset, params):
        result = original(dataset, params)
        if params.forced_components is not None:
            result.record["p_values"]["shape"] = 0.5
        return result

    monkeypatch.setattr(aac, "_analyse", drifting)
    result = aac.preflight(tmp_path, settings, cohort)
    assert not result["passed"] and result["units"][0]["p_values_differing"] == ["shape"]


def test_the_run_waits_on_the_preflight(tmp_path, cohort, settings) -> None:
    with pytest.raises(SystemExit, match="run `preflight` first"):
        aac.run_shard(tmp_path, settings, cohort, 0, 1)


# ── 1.1 Resume and refusals ──────────────────────────────────────────────────


def test_resume_skips_complete_units(tmp_path, cohort, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    assert aac.run_shard(out_dir, settings, cohort, 1, 3, limit=1) == 1
    assert aac.run_shard(out_dir, settings, cohort, 1, 3) == 1
    assert aac.run_shard(out_dir, settings, cohort, 1, 3) == 0
    lines = lsc.shard_path(out_dir, 1, 3).read_text().splitlines()
    keys = [aac.record_key(json.loads(line)) for line in lines]
    expected = [aac.unit_key(*unit) for position, unit in enumerate(aac.enumerate_units(settings)) if position % 3 == 1]
    assert sorted(keys) == sorted(expected)


def test_other_settings_are_refused_naming_the_keys(tmp_path, cohort, parent_dir, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    other = _settings(cohort, parent_dir, permutations=PERMUTATIONS + 10)
    with pytest.raises(SystemExit, match="differing: permutations"):
        aac.run_shard(out_dir, other, cohort, 0, 1)
    with pytest.raises(SystemExit, match="differing: permutations"):
        aac.preflight(out_dir, other, cohort)


def test_a_changed_parent_is_refused_by_name(tmp_path, cohort, parent_dir, export_dir) -> None:
    parent = tmp_path / "parent"
    shutil.copytree(parent_dir, parent)
    settings = _settings(cohort, parent)
    out_dir = tmp_path / "out"
    lsc.ensure_settings(out_dir, settings)
    with (parent / "rate_table.csv").open("a") as handle:
        handle.write("\n")
    for call in (lambda: aac.run_shard(out_dir, settings, cohort, 0, 1), lambda: aac.merge(out_dir, export_dir),
                 lambda: aac.decide(out_dir, export_dir)):
        with pytest.raises(SystemExit, match="Parent seaad_anchor_rates .* changed"):
            call()


def test_another_export_is_refused_by_file(tmp_path, export_dir, contract, parent_dir, settings, prepared) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    other = _copy_export(export_dir, tmp_path / "other")
    frame = _block(other, "nmr")
    frame.iloc[0, 0] += 0.5
    _rewrite(other, "nmr.csv", frame, index=True)
    moved = aac.load_cohort(other, contract)
    with pytest.raises(SystemExit, match=r"export\.files\.manifest\.json, export\.files\.nmr\.csv"):
        aac.run_shard(out_dir, _settings(moved, parent_dir), moved, 0, 1)
    with pytest.raises(SystemExit, match="differ from the settings' for: manifest.json, nmr.csv"):
        aac.run_shard(out_dir, settings, moved, 0, 1)
    with pytest.raises(SystemExit, match="differ from the settings' for: manifest.json, nmr.csv"):
        aac.merge(out_dir, other)


# ── 4.1 Merge ────────────────────────────────────────────────────────────────


def _synthetic_record(settings, point, replicate, *, p=0.5, rank=2, overrides=None) -> dict:
    """A completed unit with every arm's p-values at ``p`` unless ``overrides[arm][statistic]`` says otherwise."""

    overrides = overrides or {}
    arms = {}
    for arm in aac.point_arms(settings, point):
        spec = settings["arm_specs"][arm]
        values = {statistic: overrides.get(arm, {}).get(statistic, p) for statistic in aac.STATISTICS}
        arms[arm] = {"status": "completed", "rank": rank if spec["fit"] == "cv" else spec["fixed_rank"],
                     "measured_dimension": 2, "pair_statistics": {"delta": 1.0, "angle": 10.0, "shape": 0.1},
                     "p_values": values, "angle_null_q95": 30.0 + replicate, "seconds": 0.1}
    fits = {"cv": {"rank": rank, "in_sample_r2": 0.2, "out_of_fold_share": 0.05}}
    fits.update({aac.fit_label(r): {"rank": r, "in_sample_r2": 0.1 * r, "out_of_fold_share": 0.01}
                 for r in settings["fixed_ranks"]})
    disperse = {"stage_0": 1.0, "stage_1": 0.9, "stage_2": 1.1, "max_min_ratio": 1.1 / 0.9}
    return {
        "point": point, "replicate": replicate, "status": "completed", "inconsistent_checks": [],
        "settings_hash": lsc.settings_hash(settings), "subsample_redraws": 1, "label_redraws": 0,
        "production": {"rank": rank, "cv_mean_auroc": 0.55}, "fit_readings": fits,
        "dispersion": {coordinates: {fit: disperse for fit in fits} for coordinates in aac.COORDINATE_SYSTEMS},
        "arms": arms, "runtime_seconds": {"unit": 2.0},
    }


def test_the_rate_table_counts_rejections_at_p_le_alpha(settings) -> None:
    records = [
        _synthetic_record(settings, "n40", 0, overrides={"pls_cv": {"delta": 0.05}}),
        _synthetic_record(settings, "n40", 1, overrides={"pls_cv": {"delta": 0.051}}),
        {"point": "n40", "replicate": 2, "status": "failed"},
    ]
    rows = {(row["arm"], row["statistic"]): row for row in aac.rate_table(records, settings) if row["point"] == "n40"}
    assert rows[("pls_cv", "delta")]["rejections"] == 1 and rows[("pls_cv", "delta")]["count"] == 2
    assert rows[("pls_cv", "delta")]["failures"] == 1
    assert len([row for row in aac.rate_table(records, settings) if row["point"] == "n40"]) == 13 * 3


def test_the_paired_table_names_every_counterpart(settings) -> None:
    pairs = dict(aac.paired_arms(settings))
    for arm in settings["arms"]:
        if arm.endswith("__within_stage"):
            assert pairs.get(arm) == arm.removesuffix("__within_stage") or (arm, pairs.get(arm))
    assert pairs["pls_cv__uncorrected"] == "pls_cv"
    pairs_list = aac.paired_arms(settings)
    assert ("predictive_orthonormal", "pls_cv") in pairs_list
    assert ("predictive_orthonormal__within_stage", "pls_cv__within_stage") in pairs_list
    assert ("predictive_orthonormal__within_stage", "predictive_orthonormal") in pairs_list
    rank = settings["fixed_ranks"][1]
    assert (f"predictive_orthonormal_rank{rank}", f"pls_rank{rank}") in pairs_list
    assert (f"pls_rank{rank}__within_stage", f"pls_rank{rank}") in pairs_list
    records = [_synthetic_record(settings, "n40", replicate, p=0.01 if replicate % 2 else 0.5) for replicate in
               range(4)]
    rows = aac.paired_table(records, settings)
    assert {(row["arm"], row["reference_arm"]) for row in rows if row["point"] == "n40"} == set(pairs_list)
    assert all(row["b"] == row["c"] == 0 and row["p_two_sided"] == 1.0 for row in rows)
    records[0]["arms"]["predictive_orthonormal"]["p_values"]["delta"] = 0.01
    row = next(row for row in aac.paired_table(records, settings) if row["point"] == "n40" and row["statistic"] ==
               "delta" and row["arm"] == "predictive_orthonormal" and row["reference_arm"] == "pls_cv")
    assert (row["b"], row["c"]) == (0, 1) and row["p_gain"] == 0.5


def test_the_point_table(settings) -> None:
    records = [_synthetic_record(settings, "n40", 0, rank=2), _synthetic_record(settings, "n40", 1, rank=3,
                                                                                p=0.01)]
    rows = aac.point_rows(records, settings, 3, 27)
    value = {(row["point"], row["quantity"], row.get("arm", ""), row.get("fit", ""), row.get("coordinates", ""),
              row.get("statistic", ""), row.get("key", "")): row["value"] for row in rows}
    assert value[("n40", "cv_rank_units", "", "", "", "", "2")] == 1
    assert value[("n40", "cv_rank_share_k_minus_1", "", "", "", "", "")] == 0.5
    assert value[("n40", "features_per_participant", "", "", "", "", "")] == 27 / 40
    assert value[("n40", "in_sample_r2_median", "", "cv", "", "", "")] == 0.2
    assert value[("n40", "dispersion_median", "", "cv", "scores", "", "stage_2")] == 1.1
    assert value[("n40", "angle_null_q95_median", "pls_cv", "", "", "", "")] == 30.5
    assert value[("n40", "rank_class_rejections", "pls_cv", "", "", "delta", "above_k_minus_1")] == 1
    assert value[("n40", "rank_class_count", "pls_cv", "", "", "delta", "k_minus_1")] == 1
    assert value[("n40", "subsample_redraws", "", "", "", "", "")] == 2


# ── 5.1 Size scope ───────────────────────────────────────────────────────────


def _decision_settings() -> dict:
    specs = aac.arm_specs()
    return {
        "points": aac.design_points(ADNI_TABLE), "arm_specs": specs, "arms": list(specs), "fixed_ranks": [8, 19],
        "cohort_rank": 8, "alpha": 0.05, "rejection": "p <= alpha", "replicates": 500, "permutations": 999,
        "contract": {"name": "adni"},
    }


def _rates(settings, counts=None, *, default: int = 25, n: int = 500, missing=()) -> list[dict]:
    counts = counts or {}
    rows = []
    for point in aac.point_names(settings):
        for arm in aac.point_arms(settings, point):
            for statistic in aac.STATISTICS:
                if (point, arm, statistic) in missing:
                    continue
                hits = counts.get((point, arm, statistic), default)
                rows.append({"point": point, "arm": arm, "statistic": statistic, "count": n, "rejections": hits,
                             "rate": hits / n})
    return rows


def _scope(settings, counts, arm="predictive_orthonormal", **kwargs):
    return aac.size_scope(aac.index_rates(_rates(settings, counts, **kwargs)), arm, aac.point_names(settings))


def test_the_rule_thresholds_at_500() -> None:
    assert grade.exceedance_threshold(500) == 34
    for m, ceiling in ((3, 36), (6, 37), (9, 38), (12, 39), (21, 40)):
        assert grade.ceiling_count(500, m) == ceiling
    assert grade.allowed_exceedances(12, grade.exceedance_probability(500)) == 1
    assert grade.allowed_exceedances(21, grade.exceedance_probability(500)) == 2


def test_a_scope_starts_at_each_point_or_nowhere() -> None:
    settings = _decision_settings()
    arm = "predictive_orthonormal"
    assert _scope(settings, {})["scope_start"] == "n80"
    excess = {("n80", arm, "delta"): 50}
    assert _scope(settings, excess)["scope_start"] == "n200"
    excess[("n200", arm, "delta")] = 50
    assert _scope(settings, excess)["scope_start"] == "n500"
    excess[("n500", arm, "delta")] = 50
    assert _scope(settings, excess)["scope_start"] == "cohort"
    excess[("cohort", arm, "delta")] = 50
    judged = _scope(settings, excess)
    assert judged["scope_start"] is None and judged["scope"] is None
    assert [family["holds"] for family in judged["families"]] == [False, False, False, False]


def test_a_scope_does_not_skip_a_point() -> None:
    settings = _decision_settings()
    judged = _scope(settings, {("n200", "predictive_orthonormal", "delta"): 50})
    assert judged["scope_start"] == "n500" and "n80" not in judged["scope"]


@pytest.mark.parametrize(("start", "m", "ceiling"), [("n80", 12, 39), ("n200", 9, 38), ("n500", 6, 37),
                                                       ("cohort", 3, 36)])
def test_the_ceiling_of_each_family_size(start, m, ceiling) -> None:
    settings = _decision_settings()
    arm = "predictive_orthonormal"
    points = aac.point_names(settings)
    index = points.index(start)
    # Fail every set before ``start`` with an over-ceiling reading at the point before it.
    base = {(points[index - 1], arm, "angle"): 60} if index else {}
    at = _scope(settings, {**base, ("cohort", arm, "delta"): ceiling})
    assert at["scope_start"] == start and at["families"][index]["m"] == m
    above = _scope(settings, {**base, ("cohort", arm, "delta"): ceiling + 1})
    assert not above["families"][index]["holds"]
    assert above["families"][index]["n_above_ceiling"] == 1


def test_two_exceedances_fail_the_ladder() -> None:
    settings = _decision_settings()
    arm = "predictive_orthonormal"
    one = _scope(settings, {("n80", arm, "delta"): 35})
    two = _scope(settings, {("n80", arm, "delta"): 35, ("n500", arm, "shape"): 35})
    assert one["scope_start"] == "n80" and two["scope_start"] == "n200"


def test_an_unavailable_reading_fails_its_set() -> None:
    settings = _decision_settings()
    judged = _scope(settings, {}, missing={("n200", "predictive_orthonormal", "shape")})
    assert judged["scope_start"] == "n500"
    assert any(row.get("status") == "unavailable" and row["point"] == "n200"
               for row in judged["families"][0]["failing_readings"])


# ── 5.2 Verdict ──────────────────────────────────────────────────────────────


PASSED = {"passed": True}
SUMMARY = {"inconsistent_units": 0, "by_point": {point: {"usable": 500} for point in ("n80", "n200", "n500",
                                                                                       "cohort")}}


def _decide(counts=None, *, preflight=PASSED, summary=SUMMARY, settings=None, **kwargs) -> dict:
    settings = settings or _decision_settings()
    return aac.decision(settings, _rates(settings, counts), [], [], summary, preflight, **kwargs)


D9 = {
    "gate_failed": "No reading is interpreted. The cause is fixed and the affected units are rerun under the same "
                   "rule.",
    "production_excess": "Production's unrestricted null with df exceeds size on ADNI's noise at the named readings. "
                         "No reading of predictive_orthonormal is interpreted at a point where production fails. A "
                         "diagnosis of those points comes before any other reading on ADNI. If the cohort is among "
                         "them, the roadmap states that a rejection at ADNI's design lacks a size-controlled null on "
                         "its own noise; the case study's primary non-rejections stand, since an excess cannot "
                         "produce them.",
    "excess_persists": "The excess appears at ADNI's own design, at 0.2 features per participant, so item 19's "
                       "prediction for ADNI is refuted. The veto rests on two cohorts and the second reopening "
                       "condition of open item 1 is closed for ADNI. predictive_orthonormal stays an opt-in "
                       "diagnostic.",
    "excess_scoped": "predictive_orthonormal holds size on ADNI's blocks from the named point up to the cohort and "
                     "exceeds below it. The size part of open item 1's second condition is met for ADNI inside that "
                     "scope only. It stays an opt-in diagnostic; the gain and specificity readings at ADNI's design "
                     "are the next change and are confined to the scope.",
    "excess_absent": "predictive_orthonormal holds size on ADNI's blocks at every design point under RRPP's "
                     "unrestricted draws: SEA-AD's excess does not transfer to ADNI. The size part of open item 1's "
                     "second condition is met for ADNI. It stays an opt-in diagnostic until the gain and specificity "
                     "readings at ADNI's design, the next change. The fixed-rank grade says whether the hold is owed "
                     "to ADNI's noise or to the ranks its CV selects.",
}


def test_the_consequences_are_design_d9s() -> None:
    assert aac.CONSEQUENCES == D9 and tuple(D9) == aac.VERDICTS
    assert "SEA-AD anchor's reading stands for SEA-AD's design" in aac.STANDING


def test_the_five_verdicts() -> None:
    absent = _decide()
    assert absent["verdict"] == "excess_absent"
    assert absent["consequence"] == f"{D9['excess_absent']} {aac.STANDING}"
    assert "stays an opt-in diagnostic" in absent["consequence"]
    assert "the SEA-AD anchor's reading stands for SEA-AD's design" in absent["consequence"]
    scoped = _decide({("n80", "predictive_orthonormal", "delta"): 67})
    assert scoped["verdict"] == "excess_scoped" and scoped["scope_start"] == "n200"
    assert scoped["named_readings"][0]["point"] == "n80"
    persists = _decide({("cohort", "predictive_orthonormal", "delta"): 67})
    assert persists["verdict"] == "excess_persists"
    production = _decide({("n80", "pls_cv", "shape"): 60, ("cohort", "predictive_orthonormal", "delta"): 67})
    assert production["verdict"] == "production_excess"
    assert [(row["point"], row["statistic"]) for row in production["named_readings"]] == [("n80", "shape")]
    assert production["operating"]["candidate_interpreted_at"] == ["n200", "n500", "cohort"]
    gated = _decide({("n80", "pls_cv", "shape"): 60}, preflight=None)
    assert gated["verdict"] == "gate_failed"


def test_the_gate_still_records_the_scopes() -> None:
    gated = _decide({("n80", "predictive_orthonormal", "delta"): 67},
                    summary={**SUMMARY, "inconsistent_units": 2})
    assert gated["verdict"] == "gate_failed" and gated["gate"]["failures"] == ["2 inconsistent unit(s)"]
    assert gated["operating"]["candidate"]["scope_start"] == "n200"
    assert gated["operating"]["reference"]["scope_start"] == "n80"
    assert set(gated["within_stage_operating"]) == set(aac.WITHIN_STAGE_OPERATING_ARMS)
    assert _decide(preflight={"passed": False})["verdict"] == "gate_failed"


def test_a_failing_within_stage_arm_leaves_the_verdict() -> None:
    result = _decide({("n80", "predictive_orthonormal__within_stage", "shape"): 70})
    assert result["verdict"] == "excess_absent"
    entry = result["within_stage_operating"]["predictive_orthonormal__within_stage"]
    assert entry["scope_start"] == "n200" and entry["named_readings"][0]["statistic"] == "shape"


def test_a_decision_below_paper_grade_is_indicative() -> None:
    assert not _decide()["indicative"]["indicative"]
    small = {**_decision_settings(), "replicates": 4}
    marked = _decide(settings=small)["indicative"]
    assert marked["indicative"] and "4 replicates per point" in marked["reasons"][0]


# ── 5.3 The fixed-rank grade ─────────────────────────────────────────────────


def test_the_fixed_rank_family_has_21_readings_without_the_cohorts_rank8() -> None:
    settings = _decision_settings()
    readings = aac.fixed_rank_readings(settings, "scores", "unrestricted")
    assert len(readings) == 21
    assert not any(point == "cohort" and "rank8" in arm for point, arm, _ in readings)
    assert sum(point == "cohort" for point, _, _ in readings) == 3


def test_the_four_grades_and_the_verdict_unchanged() -> None:
    def grade_of(counts):
        result = _decide(counts)
        assert result["verdict"] == "excess_absent"
        return result["fixed_rank"]

    assert grade_of({})["grade"] == "absent"
    reproduced = grade_of({("n80", "predictive_orthonormal_rank19", "delta"): 70})
    assert reproduced["grade"] == "reproduced"
    assert reproduced["named_readings"][0]["fit"] == "rank19"
    unrestored = grade_of({("n80", "predictive_orthonormal_rank19", "delta"): 70,
                           ("n80", "predictive_orthonormal_rank19__within_stage", "delta"): 70})
    assert unrestored["grade"] == "reproduced_unrestored"
    reference = grade_of({("n200", "pls_rank8", "shape"): 70, ("n80", "predictive_orthonormal_rank19", "delta"): 70})
    assert reference["grade"] == "reference_exceeds" and reference["graded_arms"] == []
    assert reference["statement"] == aac.GRADE_STATEMENTS["reference_exceeds"]


def test_the_cohorts_operating_readings_are_not_in_the_grade() -> None:
    assert _decide({("cohort", "predictive_orthonormal", "delta"): 70})["fixed_rank"]["grade"] == "absent"


def test_the_grade_names_the_in_sample_fit() -> None:
    settings = _decision_settings()
    points = [{"point": "n80", "quantity": "in_sample_r2_median", "fit": "rank19", "value": "0.93"}]
    result = aac.decision(settings, _rates(settings, {("n80", "predictive_orthonormal_rank19", "delta"): 70}), [],
                          points, SUMMARY, PASSED)
    assert result["fixed_rank"]["named_readings"][0]["in_sample_r2_median"] == 0.93


# ── 5.4 Scoping rule and descriptive readings ────────────────────────────────


def test_the_scoping_rule_is_written_into_the_decision() -> None:
    result = _decide()
    assert list(result["scoping_rule"]) == ["S1", "S2", "S3", "S4", "S5"]
    assert "do not vouch for one another" in result["scoping_rule"]["S3"]
    assert "size part" in result["size_scope_note"] and "license no reading" in result["descriptive_note"]


def test_descriptive_readings_never_change_the_verdict_or_the_grade() -> None:
    settings = _decision_settings()
    counts = {("n80", "predictive_orthonormal", "delta"): 67}
    base = aac.decision(settings, _rates(settings, counts), [], [], SUMMARY, PASSED)
    noisy = {**counts, ("n80", "pls_cv__uncorrected", "delta"): 90}
    paired = [{"point": "n80", "statistic": "delta", "arm": "pls_cv__uncorrected", "reference_arm": "pls_cv",
               "b": 0, "c": 40, "p_gain": 1e-12}]
    points = [{"point": "n80", "quantity": "rank_class_rejections", "arm": "pls_cv", "statistic": "delta",
               "key": "k_minus_1", "value": "90"}]
    changed = aac.decision(settings, _rates(settings, noisy), paired, points, SUMMARY, PASSED)
    assert (changed["verdict"], changed["fixed_rank"]["grade"]) == (base["verdict"], base["fixed_rank"]["grade"])
    assert changed["descriptive"]["control_arm"][0]["c"] == 40


def test_the_cross_cohort_table_quotes_the_committed_seaad_counts() -> None:
    seaad_rates = lsc._read_csv(aac.PARENT_DIR / "rate_table.csv")
    seaad_cells = lsc._read_csv(aac.PARENT_DIR / "cell_table.csv")
    result = _decide(seaad_rates=seaad_rates, seaad_cells=seaad_cells)
    rows = {(row["cohort"], row["point"], row["arm"]): row for row in result["descriptive"]["cross_cohort"]}
    committed = {(row["cell"], row["arm"], row["statistic"]): (int(row["count"]), int(row["rejections"]))
                 for row in seaad_rates}
    for arm in aac.OPERATING_ARMS:
        for cell in aac.SEAAD_CELLS:
            row = rows[("seaad", cell, arm)]
            for statistic in aac.STATISTICS:
                assert (row[statistic]["count"], row[statistic]["rejections"]) == committed[(cell, arm, statistic)]
            assert row["fit_rank"] == 5 and math.isclose(row["in_sample_r2"], 0.83964)
            assert row["features_per_participant"] == 125.0
    assert rows[("seaad", "anchor", "predictive_orthonormal")]["delta"]["rejections"] == 67
    assert rows[("adni", "n80", "pls_cv")]["delta"] == {"count": 500, "rejections": 25}


# ── End to end on the stand-in: merge, decide, the docstring's commands ─────


def test_merge_and_decide_on_the_stand_in(tmp_path, cohort, settings, prepared, export_dir) -> None:
    out_dir = tmp_path / "out"
    shutil.copytree(prepared, out_dir)
    before = _snapshot(export_dir)
    aac.run_shard(out_dir, settings, cohort, 0, 1)
    merged = aac.merge(out_dir, export_dir)
    assert merged["complete_units"] == merged["expected_units"] == 6 and merged["inconsistent_units"] == 0
    result = aac.decide(out_dir, export_dir)
    assert result["verdict"] in aac.VERDICTS and result["indicative"]["indicative"]
    assert result["gate"]["passed"], result["gate"]
    assert _snapshot(export_dir) == before
    for name in ("rate_table.csv", "paired_table.csv", "point_table.csv", "accounting.json", aac.DECISION_NAME):
        text = (out_dir / name).read_text()
        assert not any(identifier in text for identifier in cohort.metadata["rid"])


def test_the_docstrings_commands_run_against_the_stand_in(monkeypatch, tmp_path, export_dir, contract,
                                                          parent_dir) -> None:
    monkeypatch.setitem(aac.CONTRACTS, "adni", contract)
    monkeypatch.setattr(aac, "POINT_SIZES", STAND_IN_POINTS)
    monkeypatch.setattr(aac, "FIXED_RANKS", _fixed_ranks(contract))
    monkeypatch.setattr(aac, "CV", SMALL_CV)
    monkeypatch.setattr(aac, "PARENT_DIR", parent_dir)
    commands = [line.strip() for line in aac.__doc__.splitlines()
                if line.strip().startswith("python scripts/adni_anchor_check.py")]
    assert [shlex.split(command)[2] for command in commands] == ["preflight", "run", "merge", "decide"]
    out_dir = tmp_path / "smoke"
    for command in commands:
        argv = shlex.split(command.replace("$EXPORT", str(export_dir)).replace("$OUT", str(out_dir)))[2:]
        assert aac.main(argv) == 0, command
    decision = json.loads((out_dir / aac.DECISION_NAME).read_text())
    assert decision["complete"]["complete_units"] == 4 * 3
