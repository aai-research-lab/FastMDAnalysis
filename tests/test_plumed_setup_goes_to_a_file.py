"""PLUMED's forty lines of setup belong beside the run, not on the screen.

At the moment the context takes the biasing force, PLUMED prints everything it
resolved: which atoms the collective variable is built from, the hill width,
the deposition pace, the bias factor, the temperature it inferred. Written
from C++ straight to the file descriptor, so it lands in the middle of a
progress bar and cannot be reached by Python logging.

Throwing it away would be the wrong fix. It is the only independent statement
of what PLUMED actually did -- reading `TORSION between atoms 7 9 15 17` in it
is how the psi selection was confirmed to be the right four atoms in the right
order. So it is captured to a file rather than discarded.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from fastmdxplora.utils.native_output import suppress_native_output


class TestTheCaptureCanKeepWhatItCatches:
    def test_native_writes_reach_the_file(self, tmp_path) -> None:
        destination = tmp_path / "plumed.log"
        with suppress_native_output(into=destination):
            os.write(1, b"PLUMED: Action TORSION\n")
        assert "Action TORSION" in destination.read_text(encoding="utf-8")

    def test_nothing_reaches_the_terminal(self, tmp_path, capfd) -> None:
        with suppress_native_output(into=tmp_path / "plumed.log"):
            os.write(1, b"PLUMED: forty lines of this\n")
        assert "forty lines" not in capfd.readouterr().out

    def test_appending_keeps_an_earlier_block(self, tmp_path) -> None:
        """A run that attaches a force more than once should not lose the
        first block to the second."""
        destination = tmp_path / "plumed.log"
        for line in (b"first\n", b"second\n"):
            with suppress_native_output(into=destination):
                os.write(1, line)
        text = destination.read_text(encoding="utf-8")
        assert "first" in text and "second" in text

    def test_discarding_still_works(self, tmp_path, capfd) -> None:
        """The trajectory loader wants the old behaviour: MDTraj's plugin
        chatter is noise, not provenance."""
        with suppress_native_output():
            os.write(1, b"dcdplugin) detected standard 32-bit DCD\n")
        assert "dcdplugin" not in capfd.readouterr().out

    def test_an_unwritable_destination_does_not_stop_the_run(
            self, tmp_path) -> None:
        """Losing a run over where a log line went would be the worse
        failure."""
        with suppress_native_output(into=tmp_path / "no" / "such" / "dir.log"):
            os.write(1, b"still fine\n")


class TestOurOwnMessageSurvives:
    # A log line raised inside the block is captured with PLUMED's, because
    # logging writes through the same descriptor -- verified by hand, and not
    # asserted here: under pytest, stdout is already captured, so a handler
    # bound to `sys.stdout` never touches descriptor 1 and the test would
    # pass or fail on the harness rather than on the code. What is worth
    # pinning is where the message is raised.

    def test_the_attachment_is_announced_outside_the_block(self) -> None:
        import fastmdxplora.simulation.runner as runner

        source = Path(runner.__file__).read_text(encoding="utf-8")
        capture = source.index("suppress_native_output(into=plumed_log)")
        announcement = source.index("PLUMED enabled: biasing force added")
        assert announcement > capture, (
            "the message must be logged after the descriptors are restored")

    def test_the_helper_no_longer_announces_it(self) -> None:
        import fastmdxplora.simulation.plumed as plumed

        source = Path(plumed.__file__).read_text(encoding="utf-8")
        assert 'logger.info(\n        "PLUMED enabled' not in source
