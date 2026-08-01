from __future__ import annotations

import numpy as np

from fastmdxplora.analysis.plotting import (
    PAPER_FIGSIZE,
    PAPER_LABEL_SIZE,
    PAPER_TICK_SIZE,
    PAPER_TITLE_SIZE,
    auto_ticks,
    apply_slide_style,
    match_colorbar_font,
    new_figure,
    save_figure,
)


def test_auto_ticks_include_zero_for_positive_data():
    ticks = auto_ticks([5, 15], max_ticks=5, include_zero=True)

    assert ticks is not None
    assert ticks[0] == 0
    assert ticks[-1] >= 15


def test_apply_slide_style_uses_paper_font_sizes():
    fig, ax = new_figure(title="style test", xlabel="Frame", ylabel="Value")
    ax.plot([0, 1, 2], [0.2, 0.5, 0.1])

    applied = apply_slide_style(
        ax,
        x_values=np.arange(3),
        y_values=[0.2, 0.5, 0.1],
        zero_x=True,
    )

    assert "x" in applied
    assert ax.xaxis.label.get_fontsize() == PAPER_LABEL_SIZE
    assert ax.yaxis.label.get_fontsize() == PAPER_LABEL_SIZE
    assert ax.get_xticklabels()[0].get_fontsize() == PAPER_TICK_SIZE
    assert ax.title.get_fontsize() == PAPER_TITLE_SIZE
    assert ax.get_xlim()[0] < 0


def test_new_figure_defaults_to_paper_size():
    fig, _ = new_figure()

    assert tuple(fig.get_size_inches()) == PAPER_FIGSIZE


def test_save_figure_applies_paper_style_to_colorbar(tmp_path):
    fig, ax = new_figure()
    mesh = ax.imshow(np.arange(4).reshape(2, 2))
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label("Scale")
    apply_slide_style(ax)
    match_colorbar_font(cbar, ax)

    out = save_figure(fig, tmp_path / "figure.png")

    assert out.is_file()
    assert cbar.ax.yaxis.label.get_fontsize() == PAPER_LABEL_SIZE


def test_save_figure_also_writes_true_vector_svg(tmp_path):
    fig, ax = new_figure(title="Vector export", xlabel="Time (ns)", ylabel="Value")
    ax.plot([0.0, 0.5, 1.0], [1.0, 1.5, 1.2])

    png = save_figure(fig, tmp_path / "vector_plot.png")
    svg = tmp_path / "vector_plot.svg"

    assert png.is_file()
    assert svg.is_file()
    text = svg.read_text(encoding="utf-8")
    assert "<svg" in text
    assert "Time (ns)" in text
    assert "Value" in text


def test_style_is_colourblind_safe_and_boxed() -> None:
    """Figures use a closed frame and a palette that survives greyscale.

    The previous Tableau palette put a red and a green next to each other,
    which converge for readers with deuteranopia and in greyscale printing.
    """
    import matplotlib.pyplot as plt

    from fastmdxplora.analysis.plotting import PALETTE, apply_style

    # Okabe-Ito, the standard colourblind-safe set for figures.
    assert PALETTE[:3] == ("#0072B2", "#D55E00", "#009E73")

    apply_style()
    assert plt.rcParams["axes.spines.top"] is True
    assert plt.rcParams["axes.spines.right"] is True
    assert plt.rcParams["xtick.direction"] == "in"
    assert plt.rcParams["ytick.direction"] == "in"
    assert plt.rcParams["xtick.minor.visible"] is True
    assert plt.rcParams["savefig.dpi"] == 300


def test_tick_density_adapts_and_categorical_axes_are_left_alone(tmp_path) -> None:
    """Every figure gets an adaptive locator, except where labels are fixed.

    The adaptive helpers existed but nothing called them: apply_slide_style
    only set a locator when given data, and the shared hook passed none, so
    plots fell back to matplotlib defaults and crowded on long trajectories.
    """
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np
    from matplotlib.ticker import AutoMinorLocator, MaxNLocator

    from fastmdxplora.analysis.plotting import new_figure, save_figure

    # A long series must not accumulate an unbounded number of ticks.
    fig, ax = new_figure()
    t = np.linspace(0, 1000, 20000)
    ax.plot(t, np.sin(t / 50))
    save_figure(fig, tmp_path / "series.png", close=False, write_svg=False)
    assert isinstance(ax.xaxis.get_major_locator(), MaxNLocator)
    assert isinstance(ax.xaxis.get_minor_locator(), AutoMinorLocator)
    assert len(ax.get_xticks()) <= 11

    # A categorical axis keeps its labels where they were put.
    fig2, ax2 = new_figure()
    labels = ["helix", "sheet", "coil", "turn"]
    ax2.bar(labels, [4, 3, 2, 1])
    save_figure(fig2, tmp_path / "bars.png", close=False, write_svg=False)
    assert [t.get_text() for t in ax2.get_xticklabels()] == labels


def test_legend_stays_legible_over_data(tmp_path) -> None:
    """A legend on a full axes needs a backing and room above the data."""
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np

    from fastmdxplora.analysis.plotting import new_figure, save_figure

    fig, ax = new_figure()
    x = np.linspace(0, 10, 500)
    ax.plot(x, np.sin(x), label="a")
    ax.plot(x, np.cos(x), label="b")
    ax.legend(loc="best")
    before = ax.get_ylim()
    save_figure(fig, tmp_path / "legend.png", close=False, write_svg=False)

    legend = ax.get_legend()
    assert legend.get_frame_on()
    assert legend.get_frame().get_alpha() == 0.85
    # Headroom was added so the legend has somewhere to sit.
    assert ax.get_ylim()[1] > before[1]
