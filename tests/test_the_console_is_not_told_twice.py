"""The failure is said once on the console, and kept in full in the log.

An analysis that raised printed ``✗ ERROR Analysis 'rdf' failed`` several
analyses before the results table was drawn, and then appeared again as its
own row in that table -- which since `81c488c` carries the reason beneath
it. The early line said only that something had failed, in a place where
nothing could be done about it yet.

It is suppressed on the console rather than demoted, because the two
audiences want different things. `fastmdxplora.log` is an audit record: an
analysis that raised is an error and belongs at ERROR with its traceback,
and a log that recorded it at DEBUG would describe a run as clean when an
analysis failed in it.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.utils.logging import _NotOnTheConsole


class _Raises(Analysis):
    name = "rdf"

    def compute(self, traj):
        raise ValueError("selection_b 'water and name O' matched no atoms")

    def plot(self, result, ax):  # pragma: no cover - never reached
        raise NotImplementedError


class TestTheConsoleIsNotToldTwice:
    def test_the_record_asks_to_be_kept_off_the_console(self, tmp_path, caplog):
        with caplog.at_level(logging.DEBUG, logger="fastmdxplora.analysis.base"):
            result = _Raises(output_dir=tmp_path).run(None)

        assert result.status == "error"
        failures = [r for r in caplog.records if "failed" in r.getMessage()]
        assert len(failures) == 1, "the failure should be logged once"
        assert getattr(failures[0], "to_console", True) is False

    def test_the_log_keeps_it_at_error_with_the_traceback(self, tmp_path, caplog):
        with caplog.at_level(logging.DEBUG, logger="fastmdxplora.analysis.base"):
            _Raises(output_dir=tmp_path).run(None)

        record = next(r for r in caplog.records if "failed" in r.getMessage())
        assert record.levelno == logging.ERROR, (
            "an analysis that raised is an error in the audit record")
        assert record.exc_info is not None, "the traceback is the useful part"

    def test_the_reason_still_reaches_the_result(self, tmp_path):
        """What the table prints comes from here, not from the log."""
        result = _Raises(output_dir=tmp_path).run(None)

        assert "matched no atoms" in result.message


class TestTheFilterItself:
    @staticmethod
    def _record(**extra):
        record = logging.LogRecord(
            "x", logging.ERROR, __file__, 1, "m", None, None)
        for key, value in extra.items():
            setattr(record, key, value)
        return record

    def test_it_drops_only_what_asks_to_be_dropped(self):
        drop = _NotOnTheConsole()

        assert drop.filter(self._record(to_console=False)) is False

    def test_everything_else_is_untouched(self):
        """Every other record in the codebase has no such attribute."""
        drop = _NotOnTheConsole()

        assert drop.filter(self._record()) is True
        assert drop.filter(self._record(to_console=True)) is True
