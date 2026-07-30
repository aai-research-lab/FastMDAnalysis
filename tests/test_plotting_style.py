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
