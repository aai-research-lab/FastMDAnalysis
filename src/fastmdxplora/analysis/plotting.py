"""Plotting utilities for analysis figures.

Provides a single point of style configuration so every analysis figure
looks consistent. The defaults are tuned for compact paper-style output:
moderate fonts, standard scientific axes, consistent line widths, and
readable colorbars without slide-scale typography baked into the images.

The module also pins the matplotlib backend to ``Agg`` when imported in a
headless environment (no DISPLAY, common on HPC nodes and CI), so analyses
do not silently hang waiting for a display.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Optional, Sequence, Union

from fastmdxplora.utils.logging import get_logger

logger = get_logger("analysis.plotting")

# Force a non-interactive backend before pyplot is imported. FastMDXplora
# always writes figures to files and never displays them, so an interactive
# backend is never wanted — and MD commonly runs headless (CI, HPC, servers)
# where interactive backends crash (e.g. "Can't find a usable init.tcl").
# This must happen before pyplot is imported. Respect an explicit MPLBACKEND.
if "MPLBACKEND" not in os.environ:
    import matplotlib

    matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  -- backend set above
from matplotlib.axes import Axes  # noqa: E402
from matplotlib.colorbar import Colorbar  # noqa: E402
from matplotlib.ticker import FixedLocator, MaxNLocator  # noqa: E402
import numpy as np  # noqa: E402


NumericSeq = Optional[Union[Sequence[float], np.ndarray]]


# Tableau 10 — colorblind-aware, distinct in print and on screen
PALETTE = (
    "#4E79A7",
    "#F28E2B",
    "#E15759",
    "#76B7B2",
    "#59A14F",
    "#EDC948",
    "#B07AA1",
    "#FF9DA7",
    "#9C755F",
    "#BAB0AC",
)

PAPER_TICK_SIZE = 9.0
PAPER_LABEL_SIZE = 10.0
PAPER_TITLE_SIZE = 11.0
PAPER_FIGSIZE = (6.5, 4.2)

# Backward-compatible names for tests/imports added while restoring styles.
V11_TICK_SIZE = PAPER_TICK_SIZE
V11_LABEL_SIZE = PAPER_LABEL_SIZE
V11_TITLE_SIZE = PAPER_TITLE_SIZE
V11_FIGSIZE = PAPER_FIGSIZE


def apply_style() -> None:
    """Apply the FastMDXplora plotting style globally.

    Called automatically by :func:`new_figure`; safe to call again to
    reset after user customizations.
    """
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "sans-serif"],
            "font.size": PAPER_TICK_SIZE,
            "figure.facecolor": "white",
            "figure.edgecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "black",
            "axes.labelcolor": "black",
            "text.color": "black",
            "xtick.color": "black",
            "ytick.color": "black",
            "axes.labelsize": PAPER_LABEL_SIZE,
            "axes.titlesize": PAPER_TITLE_SIZE,
            "axes.titleweight": "normal",
            "axes.linewidth": 1.0,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "axes.prop_cycle": plt.cycler(color=PALETTE),
            "grid.color": "#DDDDDD",
            "grid.linewidth": 0.6,
            "grid.linestyle": "--",
            "xtick.labelsize": PAPER_TICK_SIZE,
            "ytick.labelsize": PAPER_TICK_SIZE,
            "legend.fontsize": PAPER_TICK_SIZE,
            "legend.frameon": False,
            "figure.dpi": 100,
            "figure.figsize": PAPER_FIGSIZE,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white",
            "savefig.transparent": False,
            "savefig.pad_inches": 0.08,
            "lines.linewidth": 1.4,
            "lines.markersize": 3.0,
        }
    )


def _clean_array(values: NumericSeq) -> np.ndarray | None:
    if values is None:
        return None
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0:
        return None
    arr = arr[np.isfinite(arr)]
    return arr if arr.size else None


def _nice_step(span: float, max_ticks: int, integer: bool) -> float:
    if span <= 0 or not np.isfinite(span):
        return 1.0
    raw = span / max(1, max_ticks - 1)
    if raw <= 0 or not np.isfinite(raw):
        return 1.0
    magnitude = 10 ** np.floor(np.log10(raw))
    residual = raw / magnitude
    if residual <= 1.5:
        nice = 1.0
    elif residual <= 3.0:
        nice = 2.0
    elif residual <= 7.0:
        nice = 5.0
    else:
        nice = 10.0
    step = nice * magnitude
    if integer:
        step = max(1, int(round(step)))
    return float(step)


def auto_ticks(
    values: NumericSeq,
    *,
    max_ticks: int = 8,
    integer: bool = False,
    include_zero: bool = False,
) -> np.ndarray | None:
    """Return readable tick locations for compact paper-style plots."""
    arr = _clean_array(values)
    if arr is None:
        return None
    lo = float(np.min(arr))
    hi = float(np.max(arr))
    if np.isclose(lo, hi):
        delta = max(1.0, abs(lo) * 0.25)
        lo -= delta
        hi += delta
    step = _nice_step(hi - lo, max_ticks, integer)
    if step <= 0:
        return None
    start = np.floor(lo / step) * step
    stop = np.ceil(hi / step) * step
    if include_zero:
        start = min(start, 0.0)
        stop = max(stop, 0.0)
    ticks = np.arange(start, stop + step * 0.5, step, dtype=float)
    if integer:
        ticks = np.unique(np.round(ticks).astype(int)).astype(float)
    return ticks


def _as_tick_array(values: NumericSeq) -> np.ndarray | None:
    arr = _clean_array(values)
    return arr if arr is not None and arr.size else None


def apply_slide_style(
    ax: Axes,
    *,
    x_values: NumericSeq = None,
    y_values: NumericSeq = None,
    x_ticks: NumericSeq = None,
    y_ticks: NumericSeq = None,
    x_max_ticks: int = 8,
    y_max_ticks: int = 6,
    zero_x: bool = False,
    zero_y: bool = False,
    tick_size: int | float | None = None,
    label_size: int | float | None = None,
    title_size: int | float | None = None,
    x_tick_rotation: float | None = None,
) -> dict[str, np.ndarray]:
    """Apply compact paper-style tick/font defaults to an axes."""
    applied: dict[str, np.ndarray] = {}
    tick_size_val = float(tick_size or PAPER_TICK_SIZE)
    label_size_val = float(label_size or PAPER_LABEL_SIZE)
    title_size_val = float(title_size or PAPER_TITLE_SIZE)

    ticks_x = _as_tick_array(x_ticks)
    if ticks_x is None and x_values is not None:
        ticks_x = auto_ticks(x_values, max_ticks=x_max_ticks, include_zero=zero_x)
    if ticks_x is not None:
        ax.xaxis.set_major_locator(FixedLocator(ticks_x))
        applied["x"] = ticks_x
    elif x_values is not None:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=x_max_ticks))

    ticks_y = _as_tick_array(y_ticks)
    if ticks_y is None and y_values is not None:
        ticks_y = auto_ticks(y_values, max_ticks=y_max_ticks, include_zero=zero_y)
    if ticks_y is not None:
        ax.yaxis.set_major_locator(FixedLocator(ticks_y))
        applied["y"] = ticks_y
    elif y_values is not None:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=y_max_ticks))

    ax.tick_params(axis="x", which="major", labelsize=tick_size_val)
    ax.tick_params(axis="y", which="major", labelsize=tick_size_val)
    ax.xaxis.label.set_fontsize(label_size_val)
    ax.yaxis.label.set_fontsize(label_size_val)
    ax.title.set_fontsize(title_size_val)

    ax._fastmdx_tick_size = tick_size_val  # type: ignore[attr-defined]
    ax._fastmdx_label_size = label_size_val  # type: ignore[attr-defined]

    if x_tick_rotation is not None:
        for label in ax.get_xticklabels():
            label.set_rotation(x_tick_rotation)
            if x_tick_rotation:
                label.set_horizontalalignment("right")

    x_data = _clean_array(x_values)
    if x_data is not None and zero_x:
        x_max = float(np.max(x_data))
        x_min = float(np.min(x_data))
        if x_min >= 0:
            span = max(1.0, x_max - x_min)
            pad = max(0.5, span * 0.01)
            ax.set_xlim(-pad, max(ax.get_xlim()[1], x_max + pad))
    return applied


def match_colorbar_font(colorbar: Colorbar, ax: Axes) -> None:
    """Match colorbar tick and label sizes to the compact paper style."""
    tick_size = float(getattr(ax, "_fastmdx_tick_size", PAPER_TICK_SIZE))
    label_size = float(getattr(ax, "_fastmdx_label_size", PAPER_LABEL_SIZE))
    axis_name = "x" if colorbar.orientation == "horizontal" else "y"
    colorbar.ax.tick_params(axis=axis_name, labelsize=tick_size)
    getattr(colorbar.ax, f"{axis_name}axis").label.set_fontsize(label_size)


def _style_all_axes(fig: plt.Figure) -> None:
    for ax in fig.axes:
        apply_slide_style(ax)


def new_figure(
    *,
    figsize: tuple[float, float] | None = None,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    **subplot_kwargs: Any,
) -> tuple[plt.Figure, plt.Axes]:
    """Create a figure pre-styled with FastMDXplora defaults.

    Parameters
    ----------
    figsize : (float, float), optional
        Figure size in inches. Defaults to (6.5, 4.5).
    title, xlabel, ylabel : str, optional
        Set at creation time so analyses can be one-liners.
    **subplot_kwargs
        Additional keyword arguments passed to ``plt.subplots`` (e.g.
        ``nrows=2, sharex=True``).

    Returns
    -------
    (Figure, Axes)
        For multi-subplot figures the Axes object follows matplotlib's
        normal conventions (array of axes).
    """
    apply_style()
    fig, ax = plt.subplots(figsize=figsize or PAPER_FIGSIZE, **subplot_kwargs)
    fig.patch.set_facecolor("white")
    for axis in np.atleast_1d(ax).ravel():
        if hasattr(axis, "set_facecolor"):
            axis.set_facecolor("white")
    # Apply title/labels to a single Axes; for multi-subplot figures the
    # user should set these per-axis after creation.
    if not isinstance(ax, (list, tuple)) and hasattr(ax, "set_title"):
        if title:
            ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
    return fig, ax


def save_figure(
    fig: plt.Figure,
    path: str | Path,
    *,
    dpi: int = 300,
    close: bool = True,
    write_svg: bool = True,
) -> Path:
    """Save a figure to disk and (by default) close it.

    Closing is the safe default: leaving figures open eventually exhausts
    matplotlib's figure manager when many analyses run in sequence.

    Returns the resolved Path that was written.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    _style_all_axes(fig)
    fig.tight_layout()
    fig.patch.set_facecolor("white")
    for ax in fig.axes:
        ax.set_facecolor("white")
    fig.savefig(
        out,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="white",
        transparent=False,
    )

    # Every publication figure also receives a true vector SVG companion.
    # The SVG is written from the original Matplotlib figure (not converted
    # from the PNG), so text, paths, and axes remain editable and scalable.
    if write_svg and out.suffix.lower() != ".svg":
        svg_out = out.with_suffix(".svg")
        try:
            fig.savefig(
                svg_out,
                format="svg",
                bbox_inches="tight",
                facecolor="white",
                edgecolor="white",
                transparent=False,
            )
            logger.debug("Saved SVG companion: %s", svg_out)
        except Exception as exc:  # noqa: BLE001 - PNG remains the primary artifact
            logger.warning("Could not save SVG companion %s: %s", svg_out, exc)

    if close:
        plt.close(fig)
    logger.debug("Saved figure: %s", out)
    return out
