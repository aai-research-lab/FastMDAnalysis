# -*- coding: utf-8 -*-
"""FastMDXplora terminal splash and graphical abstract helpers."""

from __future__ import annotations

import os
import sys
import webbrowser
from importlib import resources
from pathlib import Path
from typing import IO

_ASSET_PACKAGE = "fastmdxplora.assets"
_GRAPHICAL_ABSTRACT = "graphical_abstract.png"


def graphical_abstract_path() -> Path:
    """Return a filesystem path to the packaged graphical abstract image."""
    try:
        return Path(str(resources.files(_ASSET_PACKAGE).joinpath(_GRAPHICAL_ABSTRACT)))
    except Exception:
        return Path(__file__).resolve().parents[1] / "assets" / _GRAPHICAL_ABSTRACT


def open_graphical_abstract() -> bool:
    """Open the graphical abstract with the user's default image/browser app."""
    path = graphical_abstract_path()
    if not path.exists():
        return False
    return bool(webbrowser.open(path.resolve().as_uri()))


def _use_color(stream: IO[str]) -> bool:
    if os.getenv("NO_COLOR"):
        return False
    if os.getenv("FASTMDX_COLOR", "1").strip().lower() in {"0", "false", "no", "off"}:
        return False
    return hasattr(stream, "isatty") and stream.isatty()


def print_startup_splash(
    stream: IO[str] | None = None,
    *,
    width: int | None = None,
    show_asset_hint: bool = True,
) -> None:
    """Legacy terminal workflow splash disabled.

    Normal runs use SessionPresenter.banner() from fastmdxplora.utils.presenter.
    This no-op prevents the old INPUTS/SETUP/SIMULATION/ANALYSIS/VISUALIZATION
    graphical workflow banner from printing in the terminal.
    """
    return
