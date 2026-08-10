"""What the software calls a flag, and what it tells you to type.

`--force` became `--force-overwrite` in the parser and stayed `--force` in
the three messages that tell somebody to use it. A refusal that names a flag
the help does not is worse than either name on its own: it sends the reader
to look for something that is only an undocumented alias.

Pinned because a rename is exactly the change that leaves messages behind --
the parser is one place and the sentences are wherever they happened to be
written.
"""

from __future__ import annotations

import pathlib

import pytest

SOURCE = pathlib.Path(__file__).resolve().parents[1] / "src" / "fastmdxplora"

#: Spellings that are not this flag.
UNRELATED = ("--force-reinstall", "--forcefield")


def _mentions_of_the_old_name() -> list[str]:
    found = []
    for path in sorted(SOURCE.rglob("*.py")):
        for number, line in enumerate(path.read_text(encoding="utf-8")
                                      .splitlines(), 1):
            if "--force" not in line:
                continue
            if any(other in line for other in UNRELATED):
                continue
            if "--force-overwrite" in line:
                continue
            stripped = line.strip()
            # The alias in the parser, and comments recording the history.
            if stripped.startswith("#") or stripped == '"--force",':
                continue
            found.append(f"{path.relative_to(SOURCE)}:{number}: {stripped}")
    return found


class TestEveryMessageNamesTheFlagTheHelpDoes:
    def test_nothing_still_tells_you_to_type_the_old_one(self) -> None:
        stale = _mentions_of_the_old_name()
        assert not stale, (
            "these tell somebody to use --force, which the help no longer "
            "lists:\n  " + "\n  ".join(stale))

    @pytest.mark.parametrize("module,phrase", [
        ("orchestrator.py", "pass --force-overwrite to overwrite it"),
        ("batch/explorer.py", "--force-overwrite to overwrite them"),
        ("cli/main.py", "Use --force-overwrite"),
    ])
    def test_each_refusal_says_the_new_one(self, module, phrase) -> None:
        text = (SOURCE / module).read_text(encoding="utf-8")
        assert phrase in text


class TestTheOldSpellingStillWorks:
    def test_it_remains_an_alias(self) -> None:
        """Renaming a published flag outright breaks whatever scripts pass
        it. It is no longer advertised; it is still accepted."""
        from fastmdxplora.cli.main import _build_parser

        explore = _build_parser()._subparsers._group_actions[0].choices["explore"]
        assert explore.parse_known_args(["--force"])[0].force is True
