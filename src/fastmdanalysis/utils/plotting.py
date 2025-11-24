# FastMDAnalysis/src/fastmdanalysis/utils/plotting.py
from __future__ import annotations

from typing import Dict, Optional, Sequence, Union

import numpy as np
from matplotlib.axes import Axes
from matplotlib.ticker import FixedLocator, MaxNLocator

NumericSeq = Optional[Union[Sequence[float], np.ndarray]]

__all__ = ["auto_ticks", "apply_slide_style"]


def _as_clean_array(values: NumericSeq) -> Optional[np.ndarray]:
    if values is None:
        return None
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0:
        return None
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    return arr


def _nice_step(span: float, max_ticks: int, integer: bool) -> float:
    if span <= 0 or not np.isfinite(span):
        return 1 if integer else 1.0
    raw = span / max(1, max_ticks - 1)
    if raw == 0 or not np.isfinite(raw):
        return 1 if integer else 1.0

    magnitude = 10 ** np.floor(np.log10(raw))
    residual = raw / magnitude

    if residual <= 1.5:
        nice = 1.0
    elif residual <= 3:
        nice = 2.0
    elif residual <= 7:
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
) -> Optional[np.ndarray]:

    arr = _as_clean_array(values)
    if arr is None:
        return None

    lo = float(np.min(arr))
    hi = float(np.max(arr))
    if np.isclose(lo, hi):
        delta = max(1.0, abs(lo) * 0.25) or 1.0
        lo -= delta
        hi += delta

    span = hi - lo
    step = _nice_step(span, max_ticks, integer)
    if step <= 0:
        return None

    start = np.floor(lo / step) * step
    stop = np.ceil(hi / step) * step

    if include_zero:
        if start > 0:
            start = 0.0
        if stop < 0:
            stop = 0.0
        if not (start <= 0 <= stop):
            if lo > 0:
                start = 0.0
            elif hi < 0:
                stop = 0.0

    ticks = np.arange(start, stop + (step * 0.5), step, dtype=float)
    if integer:
        ticks = np.unique(np.round(ticks).astype(int)).astype(float)
    return ticks


def _adaptive_font_size(
    count: Optional[int],
    *,
    max_size: float = 32.0,
    min_size: float = 16.0,
) -> float:
    if count is None or count <= 0:
        return max_size
    size = max_size / np.sqrt(max(1.0, count / 2.0))
    return float(np.clip(size, min_size, max_size))


def _set_font_sizes(
    ax: Axes,
    *,
    tick_size_x: float,
    tick_size_y: float,
    label_size_x: float,
    label_size_y: float,
    title_size: float,
) -> None:
    ax.tick_params(axis="x", which="major", labelsize=tick_size_x)
    ax.tick_params(axis="y", which="major", labelsize=tick_size_y)
    ax.xaxis.label.set_fontsize(label_size_x)
    ax.yaxis.label.set_fontsize(label_size_y)
    ax.title.set_fontsize(max(title_size, ax.title.get_fontsize()))


def _to_tick_array(data: NumericSeq) -> Optional[np.ndarray]:
    if data is None:
        return None
    arr = np.asarray(list(data), dtype=float)
    if arr.size == 0:
        return None
    return arr


def apply_slide_style(
    ax: Axes,
    *,
    x_values: NumericSeq = None,
    y_values: NumericSeq = None,
    x_ticks: NumericSeq = None,
    y_ticks: NumericSeq = None,
    x_max_ticks: int = 8,
    y_max_ticks: int = 6,
    integer_x: bool = False,
    integer_y: bool = False,
    zero_x: bool = False,
    zero_y: bool = False,
    tick_size: Optional[Union[int, float]] = None,
    label_size: Optional[Union[int, float]] = None,
    title_size: Optional[Union[int, float]] = None,
    x_tick_rotation: Optional[float] = None,
) -> Dict[str, np.ndarray]:
    applied: Dict[str, np.ndarray] = {}

    ticks_x = _to_tick_array(x_ticks)
    if ticks_x is None:
        ticks_x = auto_ticks(
            x_values,
            max_ticks=x_max_ticks,
            integer=integer_x,
            include_zero=zero_x,
        )
    if ticks_x is not None and ticks_x.size:
        ax.xaxis.set_major_locator(FixedLocator(ticks_x))
        applied["x"] = ticks_x
    else:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=x_max_ticks, integer=integer_x))

    ticks_y = _to_tick_array(y_ticks)
    if ticks_y is None:
        ticks_y = auto_ticks(
            y_values,
            max_ticks=y_max_ticks,
            integer=integer_y,
            include_zero=zero_y,
        )
    if ticks_y is not None and ticks_y.size:
        ax.yaxis.set_major_locator(FixedLocator(ticks_y))
        applied["y"] = ticks_y
    else:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=y_max_ticks, integer=integer_y))

    tick_count_x = len(applied["x"]) if "x" in applied else x_max_ticks
    tick_count_y = len(applied["y"]) if "y" in applied else y_max_ticks

    if tick_size is None:
        tick_size_x = _adaptive_font_size(tick_count_x)
        tick_size_y = _adaptive_font_size(tick_count_y)
    else:
        tick_size_x = tick_size_y = float(tick_size)

    if label_size is None:
        label_size_x = max(tick_size_x * 1.25, 20.0)
        label_size_y = max(tick_size_y * 1.25, 20.0)
    else:
        label_size_x = label_size_y = float(label_size)

    if title_size is None:
        base = max(label_size_x, label_size_y)
        title_size_val = float(np.clip(base * 1.12, 24.0, 40.0))
    else:
        title_size_val = float(title_size)

    _set_font_sizes(
        ax,
        tick_size_x=tick_size_x,
        tick_size_y=tick_size_y,
        label_size_x=label_size_x,
        label_size_y=label_size_y,
        title_size=title_size_val,
    )

    if x_tick_rotation is not None:
        for label in ax.get_xticklabels():
            label.set_rotation(x_tick_rotation)
            if x_tick_rotation and abs(x_tick_rotation) > 0:
                label.set_horizontalalignment("right")

    return applied
