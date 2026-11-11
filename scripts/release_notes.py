"""Build a GitHub release's notes from CHANGELOG.md.

The notes for a tag cover every changelog version after the previous release
tag, up to and including the tagged one. A version that was set on ``main`` but
never tagged therefore appears in the notes of the release that first shipped
it. Standard library only, so the publish workflow can run it before installing
anything:

    python3 scripts/release_notes.py v0.23.0 > release-notes.md

Exits with status 1, so nothing is published, when the changelog has no entry
for the tagged version. The previous tag is read from ``git tag``, so a CI
checkout needs the full history (``fetch-depth: 0``).
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPOSITORY = "Omnia-Lab-UC/MOTCO"

Version = tuple[int, int, int]

_VERSION = re.compile(r"v?(\d+)\.(\d+)\.(\d+)")
_HEADING = re.compile(r"^## \[?(\d+)\.(\d+)\.(\d+)\]?(?: - .*)?$")
_LINK_REFERENCE = re.compile(r"^\[[^\]]+\]: ")


def parse_version(text: str) -> Version:
    """``"v0.23.0"`` or ``"0.23.0"`` -> ``(0, 23, 0)``."""
    match = _VERSION.fullmatch(text.strip())
    if match is None:
        raise ValueError(f"not a version: {text!r}")
    return int(match[1]), int(match[2]), int(match[3])


def format_version(version: Version) -> str:
    return ".".join(str(part) for part in version)


def changelog_sections(text: str) -> dict[Version, tuple[str, str]]:
    """Map each versioned ``## `` section to its heading line and stripped body.

    ``## [Unreleased]`` and the link references at the end belong to no version.
    """
    sections: dict[Version, tuple[str, str]] = {}
    current: tuple[Version, str] | None = None
    body: list[str] = []

    def close() -> None:
        if current is not None:
            sections[current[0]] = (current[1], "\n".join(body).strip())

    for line in text.splitlines():
        heading = _HEADING.match(line)
        if heading or line.startswith("## ") or _LINK_REFERENCE.match(line):
            close()
            current = ((int(heading[1]), int(heading[2]), int(heading[3])), line) if heading else None
            body = []
        elif current is not None:
            body.append(line)
    close()
    return sections


def previous_release(version: Version, tags: list[str]) -> Version | None:
    """The latest ``v*`` tag below ``version``, or ``None`` for the first release."""
    earlier = []
    for tag in tags:
        try:
            tagged = parse_version(tag)
        except ValueError:
            continue
        if tagged < version:
            earlier.append(tagged)
    return max(earlier, default=None)


def release_notes(text: str, version: Version, previous: Version | None, repository: str | None = None) -> str:
    """The tagged version's entry, then each untagged version it first ships, newest first."""
    sections = changelog_sections(text)
    if version not in sections or not sections[version][1]:
        raise ValueError(f"CHANGELOG.md has no entry for {format_version(version)}")
    parts = [sections[version][1]]
    shipped = sorted((v for v in sections if v < version and (previous is None or v > previous)), reverse=True)
    for older in shipped:
        heading, body = sections[older]
        parts.append(re.sub(r"^## \[([^\]]+)\]", r"## \1", heading) + ("\n\n" + body if body else ""))
    if previous is not None and repository:
        parts.append(
            f"**Full diff:** https://github.com/{repository}/compare/"
            f"v{format_version(previous)}...v{format_version(version)}"
        )
    return "\n\n".join(parts) + "\n"


def _git_tags() -> list[str]:
    result = subprocess.run(["git", "tag", "--list", "v*"], cwd=ROOT, capture_output=True, text=True, check=True)
    return result.stdout.split()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build a GitHub release's notes from CHANGELOG.md.")
    parser.add_argument("tag", help="the release tag, for example v0.23.0")
    parser.add_argument("--changelog", type=Path, default=ROOT / "CHANGELOG.md")
    parser.add_argument("--previous", help="the previous release tag (default: the latest earlier v* tag)")
    parser.add_argument("--repository", default=os.environ.get("GITHUB_REPOSITORY", DEFAULT_REPOSITORY))
    args = parser.parse_args(argv)
    try:
        version = parse_version(args.tag)
        previous = parse_version(args.previous) if args.previous else previous_release(version, _git_tags())
        notes = release_notes(args.changelog.read_text(encoding="utf-8"), version, previous, args.repository)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    sys.stdout.write(notes)
    return 0


if __name__ == "__main__":
    sys.exit(main())
