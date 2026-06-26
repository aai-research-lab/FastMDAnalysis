from __future__ import annotations

from fastmdxplora.report.captions import caption_for_figure


def test_known_energy_caption() -> None:
    caption = caption_for_figure("simulation", "energy_trace")

    assert "Potential" in caption
    assert "simulation stability" in caption


def test_known_rmsd_caption() -> None:
    caption = caption_for_figure("rmsd", "rmsd")

    assert "RMSD" in caption
    assert "reference conformation" in caption


def test_fallback_caption() -> None:
    caption = caption_for_figure("custom_analysis", "custom_plot")

    assert "FastMDXplora" in caption
    assert "Custom Analysis" in caption
