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
    "orange":  "\x1b[38;5;208m",
    "purple":  "\x1b[38;5;99m",
    "white":   "\x1b[38;5;255m",
    "navy":    "\x1b[38;5;18m",
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
        import os as _os
        import sys as _sys
        import time as _time

        if self.quiet:
            return

        self._session_start = _time.monotonic()

        argv = list(_sys.argv[1:])

        def arg_value(*names: str, default: str = "") -> str:
            for i, token in enumerate(argv):
                for name in names:
                    if token == name and i + 1 < len(argv):
                        return str(argv[i + 1])
                    prefix = name + "="
                    if token.startswith(prefix):
                        return token[len(prefix):]
            return default

        def arg_list(name: str, default: str = "") -> str:
            values = []
            i = 0
            while i < len(argv):
                if argv[i] == name:
                    j = i + 1
                    while j < len(argv) and not argv[j].startswith("--"):
                        values.append(argv[j])
                        j += 1
                    break
                i += 1
            return ", ".join(values) if values else default

        system = (
            str(fields.get("System") or fields.get("system") or "")
            or arg_value("-s", "-system", "--system", default="ANY_PDB_ID")
        )
        output = (
            str(fields.get("Output") or fields.get("output") or "")
            or arg_value("--output", default="output_folder")
        )

        version = str(fields.get("Version") or fields.get("version") or "")
        if not version:
            try:
                from fastmdxplora import __version__ as _fmdx_version
                version = str(_fmdx_version)
            except Exception:
                version = "unknown"

        platform = arg_value("--simulate-platform", default="auto")
        precision = arg_value("--simulate-precision", default="mixed")

        setup_ph = arg_value("--setup-ph", default="7.4")
        ion_conc = arg_value("--setup-ion-concentration-M", default="0.15")
        forcefield = arg_value("--setup-forcefield", default="charmm36")

        nvt_steps = arg_value("--simulate-nvt-steps", default="default")
        npt_steps = arg_value("--simulate-npt-steps", default="default")
        prod_steps = arg_value("--simulate-production-steps", default="default")
        timestep = arg_value("--simulate-timestep-fs", default="2.0")
        temperature = arg_value("--simulate-temperature-K", default="300")
        friction = arg_value("--simulate-friction-per-ps", default="1.0")
        traj_interval = arg_value("--simulate-trajectory-interval-steps", default="adaptive")
        checkpoint_interval = arg_value("--simulate-checkpoint-interval-steps", default="default")

        report_title = arg_value("--report-title", default="FastMDXplora Run")
        dashboard_link = _os.environ.get("FASTMDX_DASHBOARD_URL", "")
        started = _time.strftime("%Y-%m-%d %H:%M:%S")

        H = chr(0x2500)
        V = chr(0x2502)
        TL = chr(0x256D)
        TR = chr(0x256E)
        BL = chr(0x2570)
        BR = chr(0x256F)
        CHECK = chr(0x2713)

        box_width = min(max(40, self.width - 2), 112)
        content_w = box_width - 4

        ascii_logo = [
            r" ______        _   __  __ _______   __      _                 ",
            r"|  ____|      | | |  \/  |  __ \ \ / /     | |                ",
            r"| |__ __ _ ___| |_| \  / | |  | \ V / _ __ | | ___  _ __ __ _ ",
            r"|  __/ _` / __| __| |\/| | |  | |> < | '_ \| |/ _ \| '__/ _` |",
            r"| | | (_| \__ \ |_| |  | | |__| / . \| |_) | | (_) | | | (_| |",
            r"|_|  \__,_|___/\__|_|  |_|_____/_/ \_\ .__/|_|\___/|_|  \__,_|",
            r"                                     | |                      ",
            r"                                     |_|                      ",
        ]

        def fit(text: str) -> str:
            text = str(text)
            if len(text) <= content_w:
                return text
            return text[: max(0, content_w - 3)] + "..."

        def top(color: str) -> None:
            self._write(self._c(TL + H * (box_width - 2) + TR, color))

        def bottom(color: str) -> None:
            self._write(self._c(BL + H * (box_width - 2) + BR, color))

        def line(text: str = "", border: str = "cyan", text_color: str | None = None) -> None:
            raw = fit(text)
            pad = " " * max(0, content_w - len(raw))
            body = self._c(raw, text_color) if text_color else raw
            self._write(self._c(V, border) + " " + body + pad + " " + self._c(V, border))

        def title(text: str, border: str) -> None:
            line(text, border, "white")
            line(H * len(text), border, border)

        def kv(label: str, value: str, border: str, value_color: str = "white") -> None:
            line(f"{label:<16} {value}", border, value_color)

        self._write("")

        for logo_line in ascii_logo:
            self._write(self._c(logo_line, "purple"))

        self._write(self._c("Fully Automated SysTem for Molecular Dynamics eXploration", "cyan"))
        self._write("")

        top("green")
        title("CURRENT RUN", "green")
        kv("System", system, "green")
        kv("Output", output, "green")
        kv("Version", version, "green")
        kv("Started", started, "green")
        kv("Platform", f"{platform} ({precision} precision)", "green")
        bottom("green")
        self._write("")

        top("cyan")
        title("SETUP", "cyan")
        kv("pH", setup_ph, "cyan")
        kv("Ion Conc.", f"{ion_conc} M", "cyan")
        kv("Force Field", forcefield, "cyan")
        bottom("cyan")
        self._write("")

        top("orange")
        title("SIMULATION", "orange")
        kv("NVT", f"{nvt_steps} steps", "orange")
        kv("NPT", f"{npt_steps} steps", "orange")
        kv("Production", f"{prod_steps} steps", "orange")
        kv("Timestep", f"{timestep} fs", "orange")
        kv("Temperature", f"{temperature} K", "orange")
        kv("Friction", f"{friction} / ps", "orange")
        kv("Trajectory", f"every {traj_interval} steps", "orange")
        kv("Checkpoint", f"every {checkpoint_interval} steps", "orange")
        bottom("orange")
        self._write("")

        top("purple")
        title("ANALYSIS & REPORT", "purple")
        kv("Dashboard", dashboard_link if dashboard_link else " ", "purple")
        kv("Report", report_title, "purple")
        bottom("purple")
        self._write("")

        top("blue")
        title("REPORTING & OUTPUTS", "blue")
        line(f"{CHECK} Markdown report      {CHECK} HTML summary        {CHECK} PDF figures", "blue", "white")
        line(f"{CHECK} PowerPoint slides    {CHECK} PNG/SVG plots       {CHECK} ZIP result bundle", "blue", "white")
        bottom("blue")
        self._write("")

        badges = [
            "Reproducible",
            "Config-driven",
            "OpenMM CPU/GPU",
            "Energy-aware",
            "Publication-ready",
        ]

        def badge_line(items: list[str]) -> str:
            return "  ".join(self._c(CHECK, "green") + f" {item}" for item in items)

        all_badges = badge_line(badges)
        if _visual_width(all_badges) <= self.width + 20:
            self._write(all_badges)
        else:
            self._write(badge_line(badges[:3]))
            self._write(badge_line(badges[3:]))
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
