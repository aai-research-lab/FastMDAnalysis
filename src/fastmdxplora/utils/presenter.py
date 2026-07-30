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
        self._welcome_shown: bool = False

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def welcome(
        self,
        *,
        dashboard_url: str = "http://127.0.0.1:8765",
        dashboard_enabled: bool = False,
    ) -> None:
        """Print the centered FastMDXplora startup identity once.

        The startup screen is deliberately minimal: one responsive wordmark,
        the product name and tagline, and the dashboard address.  It performs
        no simulation work and does not print the full command help.

        ``dashboard_enabled`` is retained for API compatibility with older
        callers, but the startup design intentionally avoids active/inactive
        status wording.
        """
        if self.quiet or self._welcome_shown:
            return

        self._welcome_shown = True
        _ = dashboard_enabled

        # Five-row block glyphs are rendered as one plain-text block before
        # ANSI styling is applied. This keeps width calculations accurate and
        # prevents individual lines from drifting away from the shared center.
        glyphs = {
            "A": ("01110", "10001", "11111", "10001", "10001"),
            "D": ("11110", "10001", "10001", "10001", "11110"),
            "F": ("11111", "10000", "11110", "10000", "10000"),
            "L": ("10000", "10000", "10000", "10000", "11111"),
            "M": ("10001", "11011", "10101", "10001", "10001"),
            "O": ("01110", "10001", "10001", "10001", "01110"),
            "P": ("11110", "10001", "11110", "10000", "10000"),
            "R": ("11110", "10001", "11110", "10100", "10010"),
            "S": ("01111", "10000", "01110", "00001", "11110"),
            "T": ("11111", "00100", "00100", "00100", "00100"),
            "X": ("10001", "01010", "00100", "01010", "10001"),
        }
        word = "FASTMDXPLORA"

        def render_wordmark(cell_width: int, gap: int) -> list[str]:
            on = "█" * cell_width
            off = " " * cell_width
            separator = " " * gap
            rows: list[str] = []
            for row_index in range(5):
                letters = [
                    "".join(
                        on if bit == "1" else off
                        for bit in glyphs[letter][row_index]
                    )
                    for letter in word
                ]
                rows.append(separator.join(letters).rstrip())
            return rows

        usable_width = max(36, self.width - 4)
        target_width = max(1, int(usable_width * 0.88))

        if usable_width >= 131:
            # Two-column pixels produce the largest wordmark. Select the gap
            # dynamically so the logo occupies most of wide terminals.
            gap = max(1, min(4, (target_width - 120) // 11))
            logo = render_wordmark(cell_width=2, gap=gap)
        elif usable_width >= 71:
            # One-column pixels remain crisp on normal terminals while a
            # variable gap prevents the wordmark from looking undersized.
            gap = max(1, min(4, (target_width - 60) // 11))
            logo = render_wordmark(cell_width=1, gap=gap)
        else:
            logo = ["FASTMDXPLORA"]

        def center_plain(text: str, width: int | None = None) -> str:
            available = self.width if width is None else width
            padding = max(0, (available - _visual_width(text)) // 2)
            return (" " * padding) + text

        def purple_bold(text: str) -> str:
            if not self._color:
                return text
            return f"{_C['purple']}{_C['bold']}{text}{_C['reset']}"

        # Normalize every row to one block width and apply one shared left
        # padding value. The block therefore has a single visual center axis.
        logo_width = max(_visual_width(line) for line in logo)
        normalized_logo = [line.ljust(logo_width) for line in logo]
        logo_left = max(0, (self.width - logo_width) // 2)

        self._write("")
        for logo_line in normalized_logo:
            plain_line = (" " * logo_left) + logo_line
            self._write(self._c(plain_line, "purple"))

        product_name = "FastMDXplora"
        tagline = "Fully Automated SysTem for Molecular Dynamics eXploration"

        def wrap_words(text: str, limit: int) -> list[str]:
            words = text.split()
            lines: list[str] = []
            current: list[str] = []
            current_width = 0
            for word in words:
                added = len(word) if not current else len(word) + 1
                if current and current_width + added > limit:
                    lines.append(" ".join(current))
                    current = [word]
                    current_width = len(word)
                else:
                    current.append(word)
                    current_width += added
            if current:
                lines.append(" ".join(current))
            return lines or [""]

        self._write("")
        self._write(purple_bold(center_plain(product_name)))
        for tagline_line in wrap_words(tagline, max(20, self.width - 4)):
            self._write(self._c(center_plain(tagline_line), "purple"))
        self._write("")

        # The box is centered as one unit; each content line is then centered
        # inside the box. All branding uses the same purple palette.
        box_width = min(
            max(64, _visual_width(tagline) + 8),
            max(36, self.width - 4),
        )
        inner_width = max(1, box_width - 2)

        def fit(text: str) -> str:
            if _visual_width(text) <= inner_width:
                return text
            if inner_width <= 3:
                return text[:inner_width]
            return text[: inner_width - 3] + "..."

        def box_line(text: str) -> str:
            fitted = fit(text)
            visible = _visual_width(fitted)
            left = max(0, (inner_width - visible) // 2)
            right = max(0, inner_width - visible - left)
            return "│" + (" " * left) + fitted + (" " * right) + "│"

        box_left = max(0, (self.width - box_width) // 2)
        box_prefix = " " * box_left
        top = "╭" + ("─" * inner_width) + "╮"
        bottom = "╰" + ("─" * inner_width) + "╯"

        # Keep the wordmark purple, but make the dashboard call-to-action
        # visually distinct with a cyan border, white heading/instruction,
        # and a bright cyan URL.
        self._write(self._c(box_prefix + top, "cyan"))
        self._write(
            self._c(
                box_prefix + box_line("DASHBOARD"),
                "white",
            )
        )
        self._write(
            self._c(
                box_prefix + box_line(dashboard_url),
                "cyan",
            )
        )
        self._write(
            self._c(
                box_prefix
                + box_line(
                    "Open the dashboard to configure and launch a simulation."
                ),
                "white",
            )
        )
        self._write(self._c(box_prefix + bottom, "cyan"))
        self._write("")

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

        def field_value(*keys: str, default: object = "") -> str:
            """Return a banner value from explicit fields, falling back to a default."""
            for key in keys:
                value = fields.get(key)
                if value not in (None, ""):
                    return str(value)
            return str(default)

        try:
            from fastmdxplora.setup.pipeline import DEFAULTS as _SETUP_DEFAULTS
        except Exception:
            _SETUP_DEFAULTS = {}

        try:
            from fastmdxplora.simulation.pipeline import DEFAULTS as _SIM_DEFAULTS
        except Exception:
            _SIM_DEFAULTS = {}

        platform = arg_value(
            "--simulate-platform",
            default=field_value("Platform", "platform", default=_SIM_DEFAULTS.get("platform", "auto")),
        )
        precision = arg_value(
            "--simulate-precision",
            default=field_value("Precision", "precision", default=_SIM_DEFAULTS.get("precision", "mixed")),
        )

        setup_ph = arg_value(
            "--setup-ph",
            default=field_value("pH", "ph", "setup_ph", default=_SETUP_DEFAULTS.get("ph", 7.0)),
        )
        ion_conc = arg_value(
            "--setup-ion-concentration-M",
            default=field_value(
                "Ion Conc.",
                "ion_concentration_M",
                "setup_ion_concentration_M",
                default=_SETUP_DEFAULTS.get("ion_concentration_M", 0.15),
            ),
        )
        forcefield = arg_value(
            "--setup-forcefield",
            default=field_value("Force Field", "forcefield", "setup_forcefield", default=_SETUP_DEFAULTS.get("forcefield", "charmm36")),
        )

        timestep = arg_value(
            "--simulate-timestep-fs",
            default=field_value("Timestep", "timestep_fs", "simulate_timestep_fs", default=_SIM_DEFAULTS.get("timestep_fs", 2.0)),
        )
        temperature = arg_value(
            "--simulate-temperature-K",
            default=field_value("Temperature", "temperature_K", "simulate_temperature_K", default=_SIM_DEFAULTS.get("temperature_K", 300.0)),
        )
        friction = arg_value(
            "--simulate-friction-per-ps",
            default=field_value("Friction", "friction_per_ps", "simulate_friction_per_ps", default=_SIM_DEFAULTS.get("friction_per_ps", 1.0)),
        )

        def maybe_int(value: str) -> int | None:
            if value == "":
                return None
            try:
                return int(value)
            except (TypeError, ValueError):
                return None

        def maybe_float(value: str) -> float | None:
            if value == "":
                return None
            try:
                return float(value)
            except (TypeError, ValueError):
                return None

        def format_steps(value: int | None) -> str:
            if value is None:
                return "default"
            return f"{int(value):,}"

        def resolve_stage_display() -> tuple[int | None, int | None, int | None, int | None]:
            """Resolve displayed stage steps from the same defaults as simulation."""
            try:
                from fastmdxplora.simulation.pipeline import DEFAULTS as _SIM_DEFAULTS
                from fastmdxplora.simulation.pipeline import PRESETS as _SIM_PRESETS
                from fastmdxplora.simulation.runner import plan_stages as _plan_stages
            except Exception:
                nvt = maybe_int(arg_value("--simulate-nvt-steps", default=""))
                npt = maybe_int(arg_value("--simulate-npt-steps", default=""))
                prod = maybe_int(arg_value("--simulate-production-steps", default=""))
                total = None if None in (nvt, npt, prod) else int(nvt) + int(npt) + int(prod)
                return nvt, npt, prod, total

            params = dict(_SIM_DEFAULTS)
            preset = arg_value("--simulate-preset", default="")
            if preset:
                params.update(_SIM_PRESETS.get(preset.lower(), {}))
                params["preset"] = preset.lower()

            overrides = {
                "nvt_steps": maybe_int(arg_value("--simulate-nvt-steps", default="")),
                "npt_steps": maybe_int(arg_value("--simulate-npt-steps", default="")),
                "production_steps": maybe_int(arg_value("--simulate-production-steps", default="")),
                "duration_ns": maybe_float(arg_value("--simulate-duration-ns", default="")),
                "nvt_duration_ns": maybe_float(arg_value("--simulate-nvt-duration-ns", default="")),
                "npt_duration_ns": maybe_float(arg_value("--simulate-npt-duration-ns", default="")),
                "timestep_fs": maybe_float(timestep),
            }
            params.update({k: v for k, v in overrides.items() if v is not None})
            plan = _plan_stages(
                duration_ns=params.get("duration_ns"),
                timestep_fs=float(params.get("timestep_fs", 2.0)),
                nvt_steps=params.get("nvt_steps"),
                npt_steps=params.get("npt_steps"),
                production_steps=params.get("production_steps"),
                nvt_duration_ns=params.get("nvt_duration_ns"),
                npt_duration_ns=params.get("npt_duration_ns"),
            )
            nvt = int(plan["nvt_steps"])
            npt = int(plan["npt_steps"])
            prod = int(plan["production_steps"])
            return nvt, npt, prod, nvt + npt + prod

        nvt_steps, npt_steps, prod_steps, total_steps = resolve_stage_display()

        def resolve_trajectory_display() -> str:
            raw = arg_value("--simulate-trajectory-interval-steps", default="")
            explicit = maybe_int(raw)
            if explicit is not None:
                return f"save frame every {explicit:,} production steps"
            try:
                from fastmdxplora.simulation.runner import trajectory_interval_for as _trajectory_interval_for
                interval = _trajectory_interval_for(int(prod_steps or 0))
                return f"save frame every {interval:,} production steps"
            except Exception:
                return "adaptive frame saving during production"

        trajectory_display = resolve_trajectory_display()

        report_title = arg_value("--report-title", default="FastMDXplora Run")
        # Resolve the dashboard URL from an explicitly supplied field,
        # an environment variable, or the dashboard CLI flags.
        dashboard_link = str(
            fields.get("Dashboard")
            or fields.get("dashboard_url")
            or _os.environ.get("FASTMDX_DASHBOARD_URL", "")
        )

        dashboard_enabled = (
            "--dashboard" in argv
            or "--live-dashboard" in argv
            or _os.environ.get("FASTMDX_DASHBOARD_ACTIVE") == "1"
        )

        if not dashboard_link:
            dashboard_host = arg_value(
                "--dashboard-host",
                "--host",
                default="127.0.0.1",
            )
            dashboard_port = arg_value(
                "--dashboard-port",
                "--port",
                default="8765",
            )

            # 0.0.0.0 and :: are server bind addresses, not useful browser URLs.
            display_host = (
                "127.0.0.1"
                if dashboard_host in {"0.0.0.0", "::", "[::]"}
                else dashboard_host
            )

            dashboard_link = f"http://{display_host}:{dashboard_port}"
        started = _time.strftime("%Y-%m-%d %H:%M:%S")

        H = chr(0x2500)
        V = chr(0x2502)
        TL = chr(0x256D)
        TR = chr(0x256E)
        BL = chr(0x2570)
        BR = chr(0x256F)
        CHECK = chr(0x2713)

        # Every run-information box uses one shared width and one shared
        # indentation. The boxes are centered beneath the startup wordmark,
        # while their internal content remains left-aligned for readability.
        box_width = min(112, max(72, self.width - 24))
        box_width = min(box_width, max(38, self.width - 2))
        content_w = box_width - 4
        box_indent = max(0, (self.width - box_width) // 2)
        box_prefix = " " * box_indent

        def fit(text: str, limit: int | None = None) -> str:
            text = str(text)
            available = content_w if limit is None else max(0, limit)
            if _visual_width(text) <= available:
                return text
            if available <= 3:
                return text[:available]
            return text[: max(0, available - 3)] + "..."

        def top(color: str) -> None:
            self._write(
                box_prefix + self._c(TL + H * (box_width - 2) + TR, color)
            )

        def bottom(color: str) -> None:
            self._write(
                box_prefix + self._c(BL + H * (box_width - 2) + BR, color)
            )

        def line(
            text: str = "",
            border: str = "cyan",
            text_color: str | None = None,
        ) -> None:
            raw = fit(text)
            pad = " " * max(0, content_w - _visual_width(raw))
            body = self._c(raw, text_color) if text_color else raw
            self._write(
                box_prefix
                + self._c(V, border)
                + " "
                + body
                + pad
                + " "
                + self._c(V, border)
            )

        def title(text: str, border: str) -> None:
            line(text, border, "white")
            line(H * len(text), border, border)

        def kv(
            label: str,
            value: str,
            border: str,
            *,
            label_color: str | None = None,
            value_color: str = "white",
        ) -> None:
            label_width = min(16, max(8, content_w // 3))
            label_raw = fit(label, label_width).ljust(label_width)
            value_limit = max(0, content_w - label_width - 1)
            value_raw = fit(value, value_limit)
            visible = _visual_width(label_raw) + 1 + _visual_width(value_raw)
            pad = " " * max(0, content_w - visible)
            self._write(
                box_prefix
                + self._c(V, border)
                + " "
                + self._c(label_raw, label_color or border)
                + " "
                + self._c(value_raw, value_color)
                + pad
                + " "
                + self._c(V, border)
            )

        self.welcome(
            dashboard_url=dashboard_link,
            dashboard_enabled=dashboard_enabled,
        )

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
        kv("NVT", f"{format_steps(nvt_steps)} steps", "orange")
        kv("NPT", f"{format_steps(npt_steps)} steps", "orange")
        kv("Production", f"{format_steps(prod_steps)} steps", "orange")
        kv("Total Planned", f"{format_steps(total_steps)} steps", "orange")
        kv("Timestep", f"{timestep} fs", "orange")
        kv("Temperature", f"{temperature} K", "orange")
        kv("Friction", f"{friction} / ps", "orange")
        kv("DCD Frames", trajectory_display, "orange")
        bottom("orange")
        self._write("")

        # Analysis and report information uses a blue frame, cyan labels,
        # and white values. The box itself shares the same centered layout as
        # every other run-information section; its contents remain left-aligned.
        top("blue")
        title("ANALYSIS & REPORT", "blue")
        kv(
            "Dashboard",
            dashboard_link if dashboard_link else "Not available",
            "blue",
            label_color="cyan",
        )
        kv(
            "Report",
            report_title,
            "blue",
            label_color="cyan",
        )
        bottom("blue")
        self._write("")

        # Use green here to balance the cyan, orange, and blue sections above.
        top("green")
        title("REPORTING & OUTPUTS", "green")
        line(
            f"{CHECK} Markdown report      {CHECK} HTML summary        {CHECK} PDF figures",
            "green",
            "white",
        )
        line(
            f"{CHECK} PowerPoint slides    {CHECK} PNG/SVG plots       {CHECK} ZIP result bundle",
            "green",
            "white",
        )
        bottom("green")
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
        if _visual_width(all_badges) <= box_width:
            self._write(box_prefix + all_badges)
        else:
            self._write(box_prefix + badge_line(badges[:3]))
            self._write(box_prefix + badge_line(badges[3:]))
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
