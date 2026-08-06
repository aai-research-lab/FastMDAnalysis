"""Which code a run was actually made from.

A manifest recorded the version string and nothing else. That string is
written by setuptools-scm at *install* time, so an editable install carries
whatever the version was when ``pip install -e .`` was last run and drifts
silently from then on. A real study came back stamped 2.3.0 for a run that
used a feature 2.3.0 did not have: the manifest named a version in which the
run could not have happened, and it is the number the report's reproducibility
section prints.

Where the package is imported from a source checkout, the commit says what the
version string cannot. Where it is not -- a conda or PyPI install -- there is
no checkout to ask, and the version string is the whole answer because the
distribution was built from a tag.

Two decisions are worth stating, because the alternatives are defensible:

**The checkout is found from the package, not the working directory.** A run
started inside some other repository would otherwise record that repository's
commit, which is worse than recording nothing: it is a precise claim about the
wrong code.

**A dirty tree is reported, not refused.** Refusing would block the runs
developers make all day, and this project's habit is to say what a result
cannot support rather than to withhold it. A commit beside uncommitted changes
does not describe the code that ran, so the flag is what makes the commit
honest rather than decorative.
"""

from __future__ import annotations

import subprocess
from functools import lru_cache
from pathlib import Path
from typing import Any

__all__ = ["source_checkout", "source_provenance"]

#: Long enough to be unambiguous in any real repository, short enough to read.
_SHORT = 12

#: A repository the size of a source tree answers both questions in
#: milliseconds. The limit is here so a network filesystem or a repository in a
#: strange state cannot hang a run that was only trying to describe itself.
_PATIENCE_SECONDS = 5.0


def source_checkout() -> Path | None:
    """The checkout this package was imported from, or None if installed.

    Walks up from the package rather than from the working directory: a run
    started inside another repository would otherwise report that
    repository's commit, which is a precise claim about the wrong code.
    """
    here = Path(__file__).resolve()
    for directory in here.parents:
        if (directory / ".git").exists():
            return directory
    return None


def _git(root: Path, *arguments: str) -> str | None:
    try:
        finished = subprocess.run(
            ["git", "-C", str(root), *arguments],
            capture_output=True, text=True, timeout=_PATIENCE_SECONDS)
    except (OSError, subprocess.SubprocessError):
        # No git, or it could not run. Nothing is recorded, which is the
        # honest answer rather than a guess.
        return None
    if finished.returncode != 0:
        return None
    return finished.stdout.strip()


@lru_cache(maxsize=1)
def source_provenance() -> dict[str, Any] | None:
    """The commit a run was made from, or None where there is no checkout.

    Returns ``{"commit": ..., "dirty": bool}``, and ``branch`` where the
    checkout is on one. ``dirty`` is the important field: with uncommitted
    changes the commit does not describe the code that ran, and saying so is
    what keeps the commit from being decorative.
    """
    root = source_checkout()
    if root is None:
        return None

    commit = _git(root, "rev-parse", "HEAD")
    if not commit:
        return None

    # An empty answer means nothing is modified. A failure is not an empty
    # answer, so it is reported as unknown rather than clean.
    changed = _git(root, "status", "--porcelain")
    record: dict[str, Any] = {
        "commit": commit[:_SHORT],
        "dirty": None if changed is None else bool(changed),
    }
    branch = _git(root, "rev-parse", "--abbrev-ref", "HEAD")
    if branch and branch != "HEAD":
        record["branch"] = branch
    return record


def described(record: dict[str, Any] | None) -> str | None:
    """One line a person can read, or None where there is nothing to say."""
    if not record:
        return None
    line = record["commit"]
    if record.get("branch"):
        line += f" on {record['branch']}"
    if record["dirty"] is True:
        line += " (with uncommitted changes, so the commit does not describe "
        line += "the code that ran)"
    elif record["dirty"] is None:
        line += " (whether the tree was modified could not be determined)"
    return line
