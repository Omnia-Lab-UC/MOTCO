"""Tests for ``scripts/release_notes.py`` and the changelog it reads."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

import motco

ROOT = Path(__file__).resolve().parents[1]
CHANGELOG = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")

_spec = importlib.util.spec_from_file_location("release_notes", ROOT / "scripts" / "release_notes.py")
assert _spec is not None and _spec.loader is not None
release_notes = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(release_notes)

SAMPLE = """# Changelog

## [Unreleased]

- not yet versioned

## [0.3.0] - 2026-01-03

### Added

- three

## 0.2.1 - 2026-01-02

*Not published; first released in 0.3.0.*

- two point one

## [0.2.0] - 2026-01-01

- two

## [0.1.0] - 2025-12-31

- one

[Unreleased]: https://example.org/compare/v0.3.0...HEAD
[0.3.0]: https://example.org/compare/v0.2.0...v0.3.0
"""


def test_current_version_has_a_changelog_entry():
    """A version bump must add its CHANGELOG.md entry, or the release workflow stops before publishing."""
    sections = release_notes.changelog_sections(CHANGELOG)
    version = release_notes.parse_version(motco.__version__)
    assert version in sections, f"CHANGELOG.md has no entry for {motco.__version__}"
    assert sections[version][1], f"the CHANGELOG.md entry for {motco.__version__} is empty"


def test_sections_skip_unreleased_and_link_references():
    sections = release_notes.changelog_sections(SAMPLE)
    assert sorted(sections) == [(0, 1, 0), (0, 2, 0), (0, 2, 1), (0, 3, 0)]
    assert sections[(0, 1, 0)][1] == "- one"
    assert sections[(0, 2, 1)][0] == "## 0.2.1 - 2026-01-02"


def test_notes_include_untagged_versions_since_the_previous_tag():
    notes = release_notes.release_notes(SAMPLE, (0, 3, 0), (0, 2, 0), "owner/repo")
    assert notes == (
        "### Added\n\n- three\n\n"
        "## 0.2.1 - 2026-01-02\n\n*Not published; first released in 0.3.0.*\n\n- two point one\n\n"
        "**Full diff:** https://github.com/owner/repo/compare/v0.2.0...v0.3.0\n"
    )


def test_bracketed_headings_of_shipped_versions_lose_their_link_brackets():
    notes = release_notes.release_notes(SAMPLE, (0, 3, 0), (0, 1, 0))
    assert "## 0.2.0\n\n- two" not in notes
    assert "## 0.2.0 - 2026-01-01\n\n- two" in notes
    assert "- one" not in notes
    assert "Full diff" not in notes


def test_the_first_release_takes_every_earlier_entry():
    notes = release_notes.release_notes(SAMPLE, (0, 2, 0), None)
    assert notes == "- two\n\n## 0.1.0 - 2025-12-31\n\n- one\n"


@pytest.mark.parametrize("version", [(0, 2, 2), (9, 9, 9)])
def test_a_missing_entry_is_an_error(version):
    with pytest.raises(ValueError, match="no entry"):
        release_notes.release_notes(SAMPLE, version, None)


def test_previous_release_is_the_latest_earlier_tag():
    tags = ["v0.9.0", "v0.10.0", "v0.18.0", "v0.20.0", "v0.23.0", "nightly", "v1.0.0-rc1"]
    assert release_notes.previous_release((0, 23, 0), tags) == (0, 20, 0)
    assert release_notes.previous_release((0, 20, 0), tags) == (0, 18, 0)
    assert release_notes.previous_release((0, 10, 0), tags) == (0, 9, 0)
    assert release_notes.previous_release((0, 9, 0), tags) is None


def test_v0_23_0_notes_cover_the_unpublished_0_21_and_0_22():
    notes = release_notes.release_notes(CHANGELOG, (0, 23, 0), (0, 20, 0), "Omnia-Lab-UC/MOTCO")
    assert "interaction_scan_by_block" in notes
    assert "\n## 0.22.0 - 2026-11-11\n" in notes and "SparsePLSDA" in notes
    assert "\n## 0.21.0 - 2026-11-11\n" in notes and "BLOCK_WEIGHTINGS" in notes
    assert "0.20.0 -" not in notes
    assert notes.endswith("compare/v0.20.0...v0.23.0\n")


def test_main_exits_nonzero_without_an_entry(tmp_path, capsys):
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(SAMPLE, encoding="utf-8")
    assert release_notes.main(["v9.9.9", "--changelog", str(changelog), "--previous", "v0.3.0"]) == 1
    assert "no entry for 9.9.9" in capsys.readouterr().err
    assert release_notes.main(["v0.3.0", "--changelog", str(changelog), "--previous", "v0.2.0"]) == 0
    assert capsys.readouterr().out.startswith("### Added\n\n- three\n")
