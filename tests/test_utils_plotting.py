import numpy as np

from fastmdanalysis.utils.plotting import auto_ticks, apply_slide_style


def test_auto_ticks_integer_scaling():
    ticks = auto_ticks(np.arange(0, 251), max_ticks=8, integer=True)
    assert ticks is not None
    assert np.isclose(ticks[0], 0)
    assert np.isclose(ticks[-1], 250)
    diffs = np.diff(ticks)
    assert np.allclose(diffs, diffs[0])
    assert diffs[0] >= 10


def test_auto_ticks_constant_data():
    ticks = auto_ticks([5, 5, 5], max_ticks=6, integer=True)
    assert ticks is not None
    assert (ticks[0] < 5) and (ticks[-1] > 5)
    assert np.any(np.isclose(ticks, 5))


def test_auto_ticks_include_zero_for_positive_data():
    ticks = auto_ticks([5, 15], max_ticks=5, include_zero=True)
    assert ticks is not None
    assert np.isclose(ticks[0], 0.0)
    assert ticks[-1] >= 15


def test_apply_slide_style_sets_fonts(matplotlib):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot([0, 1, 2], [0, 1, 0])
    applied = apply_slide_style(
        ax,
        x_values=np.arange(0, 51),
        y_values=[0, 1],
        integer_x=True,
        tick_size=20,
        label_size=22,
        title_size=24,
        x_tick_rotation=45,
    )

    assert "x" in applied and applied["x"].size > 0
    assert ax.xaxis.get_label().get_fontsize() == 22
    assert ax.yaxis.get_label().get_fontsize() == 22
    assert ax.title.get_fontsize() >= 24
    assert ax.get_xticklabels()[0].get_rotation() == 45

    plt.close(fig)


def test_apply_slide_style_respects_manual_ticks(matplotlib):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    manual_ticks = [0, 5, 10]
    apply_slide_style(ax, x_ticks=manual_ticks, y_values=[0, 1])
    assert np.allclose(ax.get_xticks(), manual_ticks)
    plt.close(fig)


def test_apply_slide_style_adapts_font_size(matplotlib):
    import matplotlib.pyplot as plt

    fig_sparse, ax_sparse = plt.subplots()
    apply_slide_style(ax_sparse, x_ticks=[0, 50], y_ticks=[0, 1])
    fig_sparse.canvas.draw()
    sparse_size = ax_sparse.get_xticklabels()[0].get_fontsize()

    fig_dense, ax_dense = plt.subplots()
    apply_slide_style(ax_dense, x_ticks=np.linspace(0, 100, 20), y_ticks=[0, 1])
    fig_dense.canvas.draw()
    dense_size = ax_dense.get_xticklabels()[0].get_fontsize()

    assert sparse_size > dense_size

    plt.close(fig_sparse)
    plt.close(fig_dense)
