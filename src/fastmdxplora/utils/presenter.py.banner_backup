"""Structural output layer for FastMDXplora.

While :mod:`fastmdxplora.utils.logging` handles per-record log formatting
(timestamps, levels, icons), this module handles the *structural* output
that organizes a session into a human-readable narrative: the opening
banner box, per-phase headers, the analysis-phase status table, and the
closing total. The two layers are complementary — the logger keeps doing
its job, and the presenter prints its own clean lines on top.

The presenter is **silent in quiet mode**. Quiet mode is enabled by
setting ``FASTMDX_LOG_STYLE=plain`` in the environment, or by passing
``quiet=True`` to :class:`SessionPresenter` at construction. In quiet
mode the structural output disappears entirely and the user sees only
the underlying per-line log records — the right behavior for CI, log
file redirection, and grep workflows.

Design principles:
  - **Zero external dependencies.** ANSI escape codes only. No ``rich``,
    no ``colorama``.
  - **Terminal-width aware.** The banner box auto-fits ``shutil.get_terminal_size``;
    the analysis status table aligns to its widest analysis name.
  - **TTY-aware coloring.** Colors disable automatically when stdout is
    redirected (e.g., to a file) or when ``NO_COLOR`` is set.
  - **No state coupling with the logger.** The presenter writes directly
    to stdout. The file log never receives presenter output, keeping
    ``fastmdxplora.log`` audit-friendly.
"""

from __future__ import annotations

import os
import shutil
import sys
import time
from contextlib import contextmanager
from typing import IO, Iterator


# ---------------------------------------------------------------------------
# Color palette — matches the FastMDXplora logging module so the structural
# and per-line layers feel like one tool.
# ---------------------------------------------------------------------------
_C = {
    "reset":   "\x1b[0m",
    "dim":     "\x1b[2m",
    "bold":    "\x1b[1m",
    "blue":    "\x1b[38;5;33m",
    "cyan":    "\x1b[38;5;39m",
    "green":   "\x1b[38;5;76m",
    "yellow":  "\x1b[38;5;214m",
    "red":     "\x1b[38;5;196m",
    "muted":   "\x1b[38;5;244m",
}


def _ansi_supported(stream: IO) -> bool:
    """Return True iff we should emit ANSI escapes on this stream."""
    if os.getenv("NO_COLOR"):
        return False
    if not hasattr(stream, "isatty"):
        return False
    return bool(stream.isatty())


def _strip_ansi(s: str) -> str:
    """Remove ANSI escape sequences. Used for width calculations."""
    out: list[str] = []
    i = 0
    while i < len(s):
        if s[i] == "\x1b" and i + 1 < len(s) and s[i + 1] == "[":
            j = s.find("m", i + 2)
            if j == -1:
                break
            i = j + 1
            continue
        out.append(s[i])
        i += 1
    return "".join(out)


def _visual_width(s: str) -> int:
    """Approximate displayed width (strip ANSI; count chars 1-wide).

    A full East-Asian width implementation is overkill for ASCII status
    output; this approximation is correct for the characters we actually
    emit (box-drawing, arrows, status icons, latin-1).
    """
    return len(_strip_ansi(s))


# ---------------------------------------------------------------------------
# SessionPresenter
# ---------------------------------------------------------------------------
class SessionPresenter:
    """Print structured session output to stdout.

    A single presenter instance is created by the project-level orchestrator
    at the start of a session and lives for that session. It tracks the
    current phase (for indentation) and records start times so
    :meth:`phase_end` and :meth:`done` can report elapsed durations.

    Parameters
    ----------
    stream : file-like, optional
        Where to write. Defaults to :data:`sys.stdout`.
    quiet : bool, optional
        If True, every method becomes a no-op. If omitted, quiet mode is
        auto-enabled when ``FASTMDX_LOG_STYLE=plain``. Passing
        ``quiet=False`` forces the presenter on even with that env var.
    width : int, optional
        Override the auto-detected terminal width.
    """

    # Status icons — same palette as the logger
    _STATUS_ICON = {
        "ok": ("✓", "green"),
        "error": ("✗", "red"),
        "skipped": ("·", "muted"),
        "warning": ("⚠", "yellow"),
    }

    def __init__(
        self,
        stream: IO | None = None,
        *,
        quiet: bool | None = None,
        width: int | None = None,
    ) -> None:
        self.stream: IO = stream if stream is not None else sys.stdout

        # Auto-detect quiet mode from env when not explicitly set
        if quiet is None:
            quiet = os.getenv("FASTMDX_LOG_STYLE", "").strip().lower() == "plain"
        self.quiet: bool = bool(quiet)

        # Terminal width: explicit > shutil > 80 fallback
        if width is None:
            try:
                width = shutil.get_terminal_size(fallback=(80, 24)).columns
            except OSError:
                width = 80
        self.width: int = max(40, width)  # clamp absurdly narrow terminals

        self._color: bool = _ansi_supported(self.stream)
        self._session_start: float | None = None
        self._phase_start: float | None = None
        self._current_phase: str | None = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def banner(self, **fields: str) -> None:
        """Print the opening banner box.

        Parameters
        ----------
        **fields
            ``key=value`` pairs displayed inside the box. Common keys are
            ``System``, ``Output``, ``Version``. Empty values are skipped.
        """
        if self.quiet:
            return
        self._session_start = time.monotonic()
                # Show a compact terminal version of the graphical abstract at the
        # start of interactive FastMDXplora runs. The full PNG is bundled
        # with the package and can be opened with `fastmdx splash --open`.
        title = "FastMDXplora"
        rows = [(k, v) for k, v in fields.items() if v]
        if not rows:
            self._write(self._top_rule(title))
            self._write(self._bottom_rule())
            self._write("")
            return

        label_w = max(len(k) for k, _ in rows)
        # Inner content width is bounded by terminal width: 4 chars overhead
        # for the box frame (│ <content> │).
        max_inner_w = max(20, self.width - 4)

        # Compute the wrapped lines first to determine the inner width.
        body_lines: list[str] = []
        for key, val in rows:
            label = self._color_text(f"{key+':':<{label_w + 1}}", "dim")
            full = f"{label} {val}"
            # Truncate body line if it exceeds the terminal-derived inner width.
            if _visual_width(full) > max_inner_w:
                # Strip ANSI for truncation math, then re-truncate the visible part.
                visible = _strip_ansi(full)
                truncated = visible[: max_inner_w - 1] + "…"
                full = truncated  # color is lost when we truncate, but readability wins
            body_lines.append(full)

        inner_w = min(max(_visual_width(b) for b in body_lines), max_inner_w)

        # Top rule, content, bottom rule
        self._write(self._top_rule(title, inner_w))
        for line in body_lines:
            pad = inner_w - _visual_width(line)
            self._write(f"{self._c('│', 'dim')} {line}{' ' * max(0, pad)} {self._c('│', 'dim')}")
        self._write(self._bottom_rule(inner_w))
        self._write("")

    def phase_start(self, name: str) -> None:
        """Print a phase header. Records the start time for :meth:`phase_end`."""
        if self.quiet:
            return
        self._phase_start = time.monotonic()
        self._current_phase = name
        arrow = self._c("▸", "cyan")
        label = self._c(f"Phase: {name}", "bold")
        self._write(f"{arrow} {label}")

    def phase_end(self, name: str, status: str = "ok", *, message: str | None = None) -> None:
        """Print a phase-complete summary line with elapsed time.

        Parameters
        ----------
        name : str
            The phase that just finished. Must match what was passed to
            :meth:`phase_start` for the timing to be correct.
        status : str
            ``"ok"`` | ``"error"`` | ``"skipped"`` | ``"warning"``.
        message : str, optional
            Override the default ``"{name} complete"`` text.
        """
        if self.quiet:
            return
        elapsed = (time.monotonic() - self._phase_start) if self._phase_start else 0.0
        icon, color = self._STATUS_ICON.get(status, ("·", "muted"))
        msg = message or f"{name} complete"
        line = f"  {self._c(icon, color)} {msg} {self._c(f'({self._fmt_elapsed(elapsed)})', 'muted')}"
        self._write(line)
        self._write("")
        self._phase_start = None
        self._current_phase = None

    def step(self, message: str, *, status: str = "ok") -> None:
        """Print a single indented step inside the current phase.

        Used for fine-grained progress messages within a phase, such as
        ``"Loaded input: protein.pdb"`` during setup or ``"Wrote slides.pptx"``
        during report.
        """
        if self.quiet:
            return
        icon, color = self._STATUS_ICON.get(status, ("·", "muted"))
        self._write(f"  {self._c(icon, color)} {message}")

    def info(self, message: str) -> None:
        """Print a plain indented message (no status icon).

        Useful for things like ``"Loading trajectory... 10000 frames"``
        where a status icon would be misleading.
        """
        if self.quiet:
            return
        self._write(f"  {message}")

    def analysis_table_row(
        self,
        name: str,
        status: str,
        path: str,
        elapsed: float,
        *,
        name_width: int | None = None,
    ) -> None:
        """Print one row of the analysis-phase status table.

        Format::

            ▸ rmsd      ✓ ok    →  analysis/rmsd/      (0.4 s)

        The orchestrator should pass ``name_width`` set to the longest
        analysis name in the plan so all rows align.
        """
        if self.quiet:
            return
        nw = name_width if name_width is not None else max(8, len(name))
        icon, color = self._STATUS_ICON.get(status, ("·", "muted"))
        arrow = self._c("▸", "cyan")
        path_arrow = self._c("→", "muted")
        elapsed_str = self._c(f"({self._fmt_elapsed(elapsed)})", "muted")
        line = (
            f"  {arrow} {name:<{nw}}  "
            f"{self._c(icon, color)} {status:<7}  "
            f"{path_arrow}  {path:<22}  {elapsed_str}"
        )
        self._write(line)

    def done(self, *, message: str = "Done") -> None:
        """Print the closing session-total line."""
        if self.quiet:
            return
        elapsed = (time.monotonic() - self._session_start) if self._session_start else 0.0
        line = self._c(f"{message} in {self._fmt_elapsed(elapsed)}.", "bold")
        self._write(line)

    @contextmanager
    def phase(self, name: str, *, status: str = "ok") -> Iterator["SessionPresenter"]:
        """Context manager combining :meth:`phase_start` and :meth:`phase_end`.

        Yields ``self`` so users can call ``presenter.step(...)`` inside.
        On exception, the phase is closed with ``status="error"``.
        """
        self.phase_start(name)
        try:
            yield self
            self.phase_end(name, status=status)
        except Exception:
            self.phase_end(name, status="error")
            raise

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------
    def _write(self, line: str) -> None:
        self.stream.write(line + "\n")
        try:
            self.stream.flush()
        except Exception:  # noqa: BLE001 -- flush best-effort
            pass

    def _c(self, text: str, color: str) -> str:
        """Wrap text in an ANSI color if colors are enabled."""
        if not self._color or color not in _C:
            return text
        return f"{_C[color]}{text}{_C['reset']}"

    def _color_text(self, text: str, color: str) -> str:
        return self._c(text, color)

    def _top_rule(self, title: str, inner_w: int | None = None) -> str:
        """Top of the banner box: ``╭─ Title ──────────────╮``."""
        title_part = f" {title} "
        if inner_w is None:
            inner_w = max(20, len(title) + 4)
        # left segment: ╭─{title}
        left = "╭─"
        # right segment: ──╮ -- enough dashes to fill width
        # box width = 1(╭) + 1(─) + len(title_part) + dashes + 1(╮)
        # inner content width (between │ │) = inner_w; total = inner_w + 4
        dashes_count = inner_w - len(title_part)
        if dashes_count < 2:
            dashes_count = 2
        right_dashes = "─" * dashes_count
        return self._c(f"{left}{title_part}{right_dashes}╮", "dim")

    def _bottom_rule(self, inner_w: int | None = None) -> str:
        if inner_w is None:
            inner_w = 20
        return self._c(f"╰{'─' * (inner_w + 2)}╯", "dim")

    @staticmethod
    def _fmt_elapsed(seconds: float) -> str:
        """Human-friendly elapsed time: '0.4 s', '12.3 s', '1m 23s', '1h 5m'."""
        if seconds < 60:
            return f"{seconds:.1f} s"
        if seconds < 3600:
            m, s = divmod(int(seconds), 60)
            return f"{m}m {s}s"
        h, rem = divmod(int(seconds), 3600)
        m = rem // 60
        return f"{h}h {m}m"


# ---------------------------------------------------------------------------
# Module-level singleton accessor
# ---------------------------------------------------------------------------
_PRESENTER: SessionPresenter | None = None


def get_presenter() -> SessionPresenter:
    """Return the current session presenter, creating one on first access.

    The default presenter auto-detects quiet mode and terminal width from
    the environment. Callers who need different behavior should construct
    a :class:`SessionPresenter` directly.
    """
    global _PRESENTER
    if _PRESENTER is None:
        _PRESENTER = SessionPresenter()
    return _PRESENTER


def reset_presenter() -> None:
    """Reset the singleton. Used by tests to isolate state."""
    global _PRESENTER
    _PRESENTER = None
