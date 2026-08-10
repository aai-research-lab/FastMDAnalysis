"""What the banner tells you to watch while a study runs.

An umbrella study prepares one system that every window shares, and the
preparation announces the directory it writes to. That is accurate and no use:
it is a second of setup, over before anybody could type the command, and the
twelve windows that follow are written somewhere else. A study of twelve
windows told its author to watch a directory that would not change again.
"""

from __future__ import annotations

import pytest

from fastmdxplora.utils.presenter import (
    SHARED_SETUP_DIRECTORY, _worth_watching)


class TestAPreparationPointsAtTheStudy:
    def test_the_shared_setup_resolves_upwards(self) -> None:
        assert _worth_watching("../runs/study/shared_setup") == "../runs/study"

    def test_a_trailing_separator_makes_no_difference(self) -> None:
        assert _worth_watching("../runs/study/shared_setup/") == "../runs/study"

    def test_an_absolute_path_works_too(self) -> None:
        assert _worth_watching("/tmp/study/shared_setup") == "/tmp/study"

    def test_the_name_is_declared_once(self) -> None:
        """So the banner and whatever creates the directory cannot disagree
        about what it is called."""
        assert SHARED_SETUP_DIRECTORY == "shared_setup"


class TestTheSeparatorStyleSurvives:
    """The first version split the path with `pathlib` and rebuilt it, which
    normalises the separators: a Windows caller who wrote `../runs/study` was
    handed back `..\\runs\\study`. A correct path, and not the one they typed
    -- and three tests failed on Windows CI while passing everywhere else.

    The answer is a prefix of what was given, so it is trimmed off the text
    rather than reconstructed.
    """

    def test_forward_slashes_stay_forward(self) -> None:
        assert "\\" not in _worth_watching("../runs/study/shared_setup")

    def test_a_windows_path_is_recognised_anywhere(self) -> None:
        """`os.path.split` on POSIX does not treat a backslash as a
        separator, so a platform-dependent split would leave this untouched
        when the tests run on Linux and handle it on Windows."""
        assert _worth_watching(r"..\runs\study\shared_setup") == r"..\runs\study"

    def test_a_drive_letter_survives(self) -> None:
        assert _worth_watching(r"C:\runs\study\shared_setup") == r"C:\runs\study"

    def test_a_windows_path_keeps_its_separators(self) -> None:
        assert "/" not in _worth_watching(r"..\runs\study\shared_setup")


class TestEverythingElseIsLeftAlone:
    @pytest.mark.parametrize("path", [
        "../runs/study",
        "../runs/study/runs/window-00",
        "runs/metad-smoke",
    ])
    def test_an_ordinary_run_points_at_itself(self, path: str) -> None:
        assert _worth_watching(path) == path

    def test_a_bare_name_has_no_parent_to_offer(self) -> None:
        """Pointing at "." would be worse than pointing at the directory
        that at least exists."""
        assert _worth_watching("shared_setup") == "shared_setup"

    def test_a_directory_merely_containing_the_word_is_untouched(self) -> None:
        assert _worth_watching("../runs/shared_setup_old") == \
            "../runs/shared_setup_old"


class TestTheBannerUsesIt:
    def test_the_watch_line_goes_through_the_helper(self) -> None:
        from pathlib import Path
        import fastmdxplora.utils.presenter as presenter

        source = Path(presenter.__file__).read_text(encoding="utf-8")
        assert "_worth_watching(output)" in source
