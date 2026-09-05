"""What a study prints while it runs, and what it calls overwriting.

A study of nine windows on a laptop printed one heartbeat sentence a minute
for an hour. Sixty near-identical lines bury the ones that matter -- the
windows finishing, and anything that went wrong -- so it draws a bar instead,
the same way the stage runners already do for the steps inside a run.

And `--force` did not say what it would force. It overwrites results.
"""

from __future__ import annotations


from fastmdxplora.batch.explorer import _clear_progress_line, _progress_line


class TestTheHeartbeatIsABar:
    def test_it_fills_as_runs_finish(self) -> None:
        empty = _progress_line(running=3, done=0, queued=6, total=9, seconds=1)
        half = _progress_line(running=3, done=5, queued=1, total=9, seconds=1)
        full = _progress_line(running=0, done=9, queued=0, total=9, seconds=1)
        assert empty.count("━") == 0
        assert 0 < half.count("━") < full.count("━")
        assert full.count("─") == 0

    def test_it_still_says_the_counts(self) -> None:
        """The bar replaces the sentence, not the numbers: which windows are
        running and how many are waiting is what a study is asked."""
        line = _progress_line(running=3, done=6, queued=0, total=9, seconds=90)
        assert "6/9 done" in line
        assert "3 running" in line
        assert "0 queued" in line
        assert "1m30s" in line

    def test_the_percentage_is_of_runs_not_time(self) -> None:
        line = _progress_line(running=3, done=3, queued=3, total=9, seconds=1)
        assert "33.3%" in line

    def test_no_runs_is_not_a_division(self) -> None:
        assert "0.0%" in _progress_line(
            running=0, done=0, queued=0, total=0, seconds=0)

    def test_it_fits_a_narrow_terminal(self) -> None:
        line = _progress_line(running=3, done=6, queued=0, total=9, seconds=3600)
        assert len(line) <= 80


class TestItOnlyRedrawsWhereThatMeansSomething:
    def test_a_file_gets_no_carriage_returns(self, capsys) -> None:
        """Piped to a log, a redrawn line becomes one endless line. A log is
        read after the fact and wants the history."""
        _clear_progress_line()
        assert capsys.readouterr().out == ""

    def test_the_source_checks_before_redrawing(self) -> None:
        from pathlib import Path
        import fastmdxplora.batch.explorer as explorer

        source = Path(explorer.__file__).read_text(encoding="utf-8")
        assert "if _is_a_terminal():" in source
        assert "_clear_progress_line()" in source


class TestOverwritingSaysWhatItOverwrites:
    @staticmethod
    def _explore_parser():
        from fastmdxplora.cli.main import _build_parser

        return _build_parser()._subparsers._group_actions[0].choices["explore"]

    def test_the_flag_names_the_thing(self) -> None:
        options = [s for a in self._explore_parser()._actions
                   for s in a.option_strings]
        assert "--force-overwrite" in options

    def test_the_old_spelling_still_works(self) -> None:
        """Renaming a published flag outright breaks whatever scripts already
        pass it, and the rename is for clarity, not to make a point."""
        parser = self._explore_parser()
        assert parser.parse_known_args(["--force"])[0].force is True
        assert parser.parse_known_args(["--force-overwrite"])[0].force is True

    def test_both_reach_the_same_setting(self) -> None:
        parser = self._explore_parser()
        old = parser.parse_known_args(["--force"])[0]
        new = parser.parse_known_args(["--force-overwrite"])[0]
        assert old.force == new.force
