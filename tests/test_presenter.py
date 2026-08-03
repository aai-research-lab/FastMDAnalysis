"""Tests for the SessionPresenter structural output layer."""

from __future__ import annotations

import io
import re

import pytest

from fastmdxplora.utils.presenter import (
    SessionPresenter,
    _ansi_supported,
    _strip_ansi,
    _visual_width,
)


# ---------------------------------------------------------------------------
# Helper: build a presenter that writes to a captured StringIO
# ---------------------------------------------------------------------------
def _presenter(quiet: bool = False, width: int = 80) -> tuple[SessionPresenter, io.StringIO]:
    buf = io.StringIO()
    p = SessionPresenter(stream=buf, quiet=quiet, width=width)
    return p, buf


# ===========================================================================
# Helper functions
# ===========================================================================
class TestHelpers:
    def test_strip_ansi(self):
        s = "\x1b[31mhello\x1b[0m world"
        assert _strip_ansi(s) == "hello world"

    def test_visual_width_ignores_ansi(self):
        s = "\x1b[31mABC\x1b[0m"
        assert _visual_width(s) == 3

    def test_ansi_supported_no_color_env(self, monkeypatch):
        monkeypatch.setenv("NO_COLOR", "1")
        # Doesn't matter if isatty would say True, NO_COLOR wins
        class FakeTTY:
            def isatty(self):
                return True

        assert _ansi_supported(FakeTTY()) is False

    def test_ansi_supported_non_tty(self):
        # StringIO is not a TTY → no ANSI
        assert _ansi_supported(io.StringIO()) is False


# ===========================================================================
# Banner
# ===========================================================================
class TestBanner:
    def test_emits_box_drawing_characters(self):
        p, buf = _presenter()
        p.banner(System="x.pdb", Output="/tmp/out", Version="0.1.0")
        out = buf.getvalue()
        assert "╭" in out
        assert "╮" in out
        assert "╰" in out
        assert "╯" in out
        assert "│" in out

    def test_includes_all_supplied_fields(self):
        p, buf = _presenter()
        p.banner(System="protein.pdb", Output="/tmp/x", Version="0.1.0")
        out = buf.getvalue()
        assert "protein.pdb" in out
        assert "/tmp/x" in out
        assert "0.1.0" in out

    def test_skips_empty_fields(self):
        p, buf = _presenter()
        p.banner(System="x.pdb", Output="", Version="0.1.0")
        out = buf.getvalue()
        assert "x.pdb" in out
        # Empty value not rendered
        assert "Output:" not in _strip_ansi(out)

    def test_title_is_FastMDXplora(self):
        p, buf = _presenter()
        p.banner(System="x.pdb")
        out = _strip_ansi(buf.getvalue())
        assert "FastMDXplora" in out


# ===========================================================================
# Phase headers and status
# ===========================================================================
class TestPhases:
    def test_phase_start_emits_arrow_and_name(self):
        p, buf = _presenter()
        p.phase_start("setup")
        out = _strip_ansi(buf.getvalue())
        assert "▸" in out
        assert "Phase: setup" in out

    def test_phase_end_emits_status_icon(self):
        p, buf = _presenter()
        p.phase_start("setup")
        p.phase_end("setup", status="ok")
        out = _strip_ansi(buf.getvalue())
        assert "✓" in out
        assert "setup complete" in out

    def test_phase_end_failure_icon(self):
        p, buf = _presenter()
        p.phase_start("setup")
        p.phase_end("setup", status="error")
        out = _strip_ansi(buf.getvalue())
        assert "✗" in out

    def test_phase_end_includes_elapsed_time(self):
        p, buf = _presenter()
        p.phase_start("setup")
        p.phase_end("setup", status="ok")
        # Looking for pattern "(N.N s)"
        assert re.search(r"\(\d+\.\d+ s\)", _strip_ansi(buf.getvalue()))

    def test_step_indents_two_spaces(self):
        p, buf = _presenter()
        p.step("Loaded input PDB")
        line = buf.getvalue()
        assert line.startswith("  ")
        assert "Loaded input PDB" in line

    def test_info_uses_plain_indent_no_icon(self):
        p, buf = _presenter()
        p.info("Loading trajectory... 10000 frames")
        out = _strip_ansi(buf.getvalue())
        assert "  Loading trajectory" in out
        # No status icon since this is plain info, not a step
        assert "✓" not in out
        assert "✗" not in out


# ===========================================================================
# Analysis status table
# ===========================================================================
class TestAnalysisTable:
    def test_row_shows_name_status_path_elapsed(self):
        p, buf = _presenter()
        p.analysis_table_row("rmsd", "ok", "analysis/rmsd/", 0.42)
        out = _strip_ansi(buf.getvalue())
        assert "rmsd" in out
        assert "ok" in out
        assert "analysis/rmsd/" in out
        assert re.search(r"\(0\.\d+ s\)", out)

    def test_rows_align_to_name_width(self):
        p, buf = _presenter()
        # Pass name_width to force alignment
        p.analysis_table_row("rmsd", "ok", "analysis/rmsd/", 0.4, name_width=10)
        p.analysis_table_row("dihedrals", "ok", "analysis/dihedrals/", 0.8, name_width=10)
        lines = [_strip_ansi(line) for line in buf.getvalue().splitlines()]
        # Both rows should have the status column at the same column index
        idx1 = lines[0].index("ok")
        idx2 = lines[1].index("ok")
        assert idx1 == idx2

    def test_error_row_uses_error_icon(self):
        p, buf = _presenter()
        p.analysis_table_row("rmsd", "error", "analysis/rmsd/", 0.0)
        out = _strip_ansi(buf.getvalue())
        assert "✗" in out


# ===========================================================================
# Done footer
# ===========================================================================
class TestDone:
    def test_done_emits_total(self):
        p, buf = _presenter()
        p.banner(System="x.pdb")  # starts the session timer
        p.done()
        out = _strip_ansi(buf.getvalue())
        assert "Done" in out
        # Format is "Done in N.N s." or "Done in Nm Ns."
        assert re.search(r"Done in \S+", out)

    def test_done_without_session_start_uses_zero(self):
        """If banner() was never called, done() reports 0.0 s — should not crash."""
        p, buf = _presenter()
        p.done()  # never called banner
        out = _strip_ansi(buf.getvalue())
        assert "Done" in out


# ===========================================================================
# Quiet mode
# ===========================================================================
class TestQuietMode:
    def test_explicit_quiet_silences_everything(self):
        p, buf = _presenter(quiet=True)
        p.banner(System="x.pdb")
        p.phase_start("setup")
        p.step("Did stuff")
        p.phase_end("setup")
        p.analysis_table_row("rmsd", "ok", "analysis/rmsd/", 0.4)
        p.done()
        assert buf.getvalue() == ""

    def test_env_plain_auto_enables_quiet(self, monkeypatch):
        monkeypatch.setenv("FASTMDX_LOG_STYLE", "plain")
        # Use quiet=None so env autodetection kicks in
        buf = io.StringIO()
        p = SessionPresenter(stream=buf, quiet=None, width=80)
        p.banner(System="x.pdb")
        p.phase_start("setup")
        assert buf.getvalue() == ""

    def test_explicit_quiet_false_overrides_env(self, monkeypatch):
        monkeypatch.setenv("FASTMDX_LOG_STYLE", "plain")
        # explicit quiet=False wins over env
        buf = io.StringIO()
        p = SessionPresenter(stream=buf, quiet=False, width=80)
        p.banner(System="x.pdb")
        assert "FastMDXplora" in _strip_ansi(buf.getvalue())


# ===========================================================================
# Color handling
# ===========================================================================
class TestColor:
    def test_no_color_when_stream_is_not_tty(self):
        p, buf = _presenter()
        p.phase_start("setup")
        assert "\x1b[" not in buf.getvalue()

    def test_no_color_when_env_set(self, monkeypatch):
        # Simulate a TTY by stubbing isatty
        monkeypatch.setenv("NO_COLOR", "1")

        class FakeTTY:
            def isatty(self):
                return True

            def write(self, s):
                pass

            def flush(self):
                pass

        p = SessionPresenter(stream=FakeTTY(), width=80)
        # Internal _color is False due to NO_COLOR
        assert p._color is False


# ===========================================================================
# Context manager
# ===========================================================================
class TestPhaseContextManager:
    def test_normal_exit_closes_with_ok(self):
        p, buf = _presenter()
        with p.phase("setup"):
            p.step("did a thing")
        out = _strip_ansi(buf.getvalue())
        assert "setup complete" in out
        # No error icon
        assert "✗" not in out

    def test_exception_closes_with_error_and_reraises(self):
        p, buf = _presenter()
        with pytest.raises(RuntimeError, match="boom"):
            with p.phase("setup"):
                raise RuntimeError("boom")
        out = _strip_ansi(buf.getvalue())
        assert "✗" in out


# ===========================================================================
# Time formatting
# ===========================================================================
class TestTimeFormat:
    def test_seconds(self):
        assert SessionPresenter._fmt_elapsed(0.42) == "0.4 s"
        assert SessionPresenter._fmt_elapsed(12.345) == "12.3 s"

    def test_minutes(self):
        assert SessionPresenter._fmt_elapsed(83.0) == "1m 23s"

    def test_hours(self):
        assert SessionPresenter._fmt_elapsed(3725.0) == "1h 2m"


# ===========================================================================
# Width awareness
# ===========================================================================
def test_narrow_width_clamps_to_minimum():
    """Absurdly narrow terminals should not crash the banner."""
    p = SessionPresenter(stream=io.StringIO(), width=10)
    # Should clamp to at least 40
    assert p.width >= 40


def test_banner_respects_terminal_width():
    """Banner inner width should not exceed terminal width."""
    p, buf = _presenter(width=50)
    p.banner(System="x" * 100, Output="y" * 100)  # very long values
    lines = buf.getvalue().splitlines()
    # No emitted line should exceed terminal width (allowing for some margin)
    for line in lines:
        if line.strip():
            assert _visual_width(line) <= 50 + 20  # generous slack for box chars


def test_orchestrator_reuses_the_process_presenter() -> None:
    """The banner is shown once per presenter, so there must be only one.

    The orchestrator used to construct its own SessionPresenter, which made
    the startup wordmark print twice: once from the CLI and again when a
    study began.
    """
    from fastmdxplora.orchestrator import FastMDXplora
    from fastmdxplora.utils import presenter as presenter_module

    presenter_module._PRESENTER = None
    cli_presenter = presenter_module.get_presenter()
    fmdx = FastMDXplora(system="1L2Y", output_dir="/tmp/does-not-need-to-exist")
    assert fmdx._presenter is cli_presenter


def test_run_summary_lists_the_gui_only_when_one_was_requested(monkeypatch) -> None:
    """`explore` does not start a server, so it must not advertise one.

    The run summary printed the dashboard address unconditionally, sending
    users to a refused connection.
    """
    import io
    import sys

    from fastmdxplora.utils import presenter as presenter_module

    # Other tests activate dashboard telemetry in this process; make the
    # precondition explicit instead of depending on test order.
    monkeypatch.delenv("FASTMDX_DASHBOARD_ACTIVE", raising=False)
    monkeypatch.delenv("FASTMDX_DASHBOARD_URL", raising=False)

    def render(argv: list[str]) -> str:
        presenter_module._PRESENTER = None
        pr = presenter_module.get_presenter()
        buffer = io.StringIO()
        pr.stream = buffer
        pr._color = False
        pr.width = 120
        monkeypatch.setattr(sys, "argv", ["fastmdx", *argv])
        pr.banner(System="1L2Y", Output="/tmp/run")
        return buffer.getvalue()

    without = render(["explore", "--system", "1L2Y"])
    assert "127.0.0.1" not in without

    with_gui = render(["explore", "--system", "1L2Y", "--dashboard"])
    assert "127.0.0.1" in with_gui


def test_banner_is_one_left_aligned_box() -> None:
    """All sections share a single box, aligned with the wordmark.

    Five separate centred boxes drifted to the middle of a wide terminal and
    read as five unrelated blocks.
    """
    import io
    import sys

    from fastmdxplora.utils import presenter as presenter_module

    presenter_module._PRESENTER = None
    pr = presenter_module.get_presenter()
    buffer = io.StringIO()
    pr.stream = buffer
    pr._color = False
    pr.width = 200
    argv = sys.argv
    sys.argv = ["fastmdx", "explore", "--system", "1L2Y"]
    try:
        pr.banner(System="1L2Y", Output="/tmp/run", Version="2.1.0")
    finally:
        sys.argv = argv
    text = buffer.getvalue()

    # Exactly one box: one top-left corner and one bottom-left corner.
    assert text.count(chr(0x256D)) == 1
    assert text.count(chr(0x2570)) == 1

    # Every box line starts at the same two-space indent as the wordmark.
    box_lines = [ln for ln in text.splitlines() if chr(0x2502) in ln or chr(0x256D) in ln]
    assert box_lines
    assert all(ln.startswith("  ") and not ln.startswith("   ") for ln in box_lines)

    for heading in ("MD EXPLORATION", "SETUP", "SIMULATION", "REPORTING & OUTPUTS"):
        assert heading in text


class TestBannerReflectsTheRun:
    """The summary box must show what was asked for, not what is usual.

    The same option is spelled two ways: `fastmdx explore --setup-forcefield`
    and `fastmdx setup --forcefield`. The banner knew only the first, so every
    per-phase run displayed defaults regardless of its arguments, which is
    worse than showing nothing.
    """

    @staticmethod
    def _banner(argv):
        import contextlib
        import io
        import sys

        from fastmdxplora.utils.presenter import SessionPresenter

        original = sys.argv
        sys.argv = argv
        buffer = io.StringIO()
        try:
            with contextlib.redirect_stdout(buffer):
                SessionPresenter().banner(
                    System="4W52", Output="/tmp/x", Version="0.0.0"
                )
        finally:
            sys.argv = original
        return buffer.getvalue()

    def test_the_per_phase_spelling_is_recognised(self) -> None:
        text = self._banner([
            "fastmdx", "setup", "--system", "4W52",
            "--forcefield", "amber-openff", "--ph", "6.5",
        ])
        assert "amber-openff" in text
        assert "6.5" in text

    def test_the_explore_spelling_still_works(self) -> None:
        text = self._banner([
            "fastmdx", "explore", "--system", "4W52",
            "--setup-forcefield", "amber-openff", "--setup-ph", "6.5",
        ])
        assert "amber-openff" in text
        assert "6.5" in text

    def test_the_resolved_force_field_is_shown_not_the_word_auto(self) -> None:
        """"auto" tells the reader nothing about what was simulated."""
        text = self._banner(["fastmdx", "setup", "--system", "4W52"])
        assert "amber-openff" in text
        # ...while still making clear it was chosen rather than named.
        assert "(auto)" in text

    def test_a_named_force_field_is_shown_as_given(self) -> None:
        text = self._banner([
            "fastmdx", "setup", "--system", "4W52", "--forcefield", "charmm36",
        ])
        assert "charmm36" in text
        assert "(auto)" not in text

    def test_simulation_settings_follow_the_same_rule(self) -> None:
        text = self._banner([
            "fastmdx", "simulate", "--nvt-steps", "1234",
            "--production-steps", "9876",
        ])
        assert "1,234" in text or "1234" in text
        assert "9,876" in text or "9876" in text
