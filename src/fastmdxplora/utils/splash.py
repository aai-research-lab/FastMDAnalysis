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
    """Print the FastMDXplora terminal welcome banner.

    Set FASTMDX_SPLASH=0 to disable the banner.
    Set FASTMDX_COLOR=0 to disable terminal colors.
    """
    if os.getenv("FASTMDX_SPLASH", "1").strip().lower() in {"0", "false", "no", "off"}:
        return

    stream = stream if stream is not None else sys.stdout
    color = _use_color(stream)

    BOLD = "\033[1m" if color else ""
    RESET = "\033[0m" if color else ""
    BLUE = "\033[38;5;33m" if color else ""
    NAVY = "\033[38;5;18m" if color else ""
    CYAN = "\033[38;5;45m" if color else ""
    GREEN = "\033[38;5;40m" if color else ""
    ORANGE = "\033[38;5;208m" if color else ""
    PURPLE = "\033[38;5;99m" if color else ""
    GRAY = "\033[38;5;245m" if color else ""
    WHITE = "\033[38;5;255m" if color else ""

    stream.write("\n")
    stream.write(f"{BOLD}{NAVY}")
    stream.write(r"""
 ______        _   __  __ _______   __      _
|  ____|      | | |  \/  |  __ \ \ / /     | |
| |__ __ _ ___| |_| \  / | |  | \ V / _ __ | | ___  _ __ __ _
|  __/ _` / __| __| |\/| | |  | |> < | '_ \| |/ _ \| '__/ _` |
| | | (_| \__ \ |_| |  | | |__| / . \| |_) | | (_) | | | (_| |
|_|  \__,_|___/\__|_|  |_|_____/_/ \_\ .__/|_|\___/|_|  \__,_|
                                     | |
                                     |_|
""")
    stream.write(f"{RESET}")

    stream.write(
        f"{CYAN}{BOLD}        Fast, modular, reproducible MD analysis from setup to results{RESET}\n\n"
    )

    stream.write(f"{BLUE}┌──────────────────────────────┐{RESET}        {GREEN}┌──────────────────────────────┐{RESET}\n")
    stream.write(f"{BLUE}│{WHITE}{BOLD}  INPUTS                      {BLUE}│{RESET}        {GREEN}│{WHITE}{BOLD}  1. SETUP                    {GREEN}│{RESET}\n")
    stream.write(f"{BLUE}│{RESET}  [PDB/mmCIF] Structure       {BLUE}│{RESET}        {GREEN}│{RESET}  Build and prepare system    {GREEN}│{RESET}\n")
    stream.write(f"{BLUE}│{RESET}  [SDF/MOL2] Ligand/topology  {BLUE}│{RESET}        {GREEN}│{RESET}  Validate files and params   {GREEN}│{RESET}\n")
    stream.write(f"{BLUE}│{RESET}  [YAML/TOML] Configuration   {BLUE}│{RESET}        {GREEN}│{RESET}  Create reproducible run     {GREEN}│{RESET}\n")
    stream.write(f"{BLUE}│{RESET}  [Force fields] Parameters   {BLUE}│{RESET}        {GREEN}└──────────────┬───────────────┘{RESET}\n")
    stream.write(f"{BLUE}└──────────────┬───────────────┘{RESET}                       {GRAY}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                                       {GRAY}v{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                    {ORANGE}┌──────────────────────────────┐{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                    {ORANGE}│{WHITE}{BOLD}  2. SIMULATION               {ORANGE}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                    {ORANGE}│{RESET}  OpenMM engine CPU/GPU       {ORANGE}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                    {ORANGE}│{RESET}  Minimization, NVT, NPT      {ORANGE}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                    {ORANGE}│{RESET}  Production and restart      {ORANGE}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                    {ORANGE}└──────────────┬───────────────┘{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                                   {GRAY}│{RESET}\n")
    stream.write(f"               {GRAY}v{RESET}                                   {GRAY}v{RESET}\n")

    stream.write(f"{CYAN}┌──────────────────────────────┐{RESET}        {PURPLE}┌──────────────────────────────┐{RESET}\n")
    stream.write(f"{CYAN}│{WHITE}{BOLD}  FASTMDXPLORA CORE           {CYAN}│{RESET}        {PURPLE}│{WHITE}{BOLD}  3. ANALYSIS                 {PURPLE}│{RESET}\n")
    stream.write(f"{CYAN}│{RESET}                              {CYAN}│{RESET}        {PURPLE}│{RESET}  RMSD, RMSF, Rg, H-bonds     {PURPLE}│{RESET}\n")
    stream.write(f"{CYAN}│{RESET}       ╭────────────╮         {CYAN}│{RESET}        {PURPLE}│{RESET}  Contacts, PCA, clustering   {PURPLE}│{RESET}\n")
    stream.write(f"{CYAN}│{RESET}       │ MD WORKFLOW│         {CYAN}│{RESET}        {PURPLE}│{RESET}  Free-energy summaries       {PURPLE}│{RESET}\n")
    stream.write(f"{CYAN}│{RESET}       ╰────────────╯         {CYAN}│{RESET}        {PURPLE}└──────────────┬───────────────┘{RESET}\n")
    stream.write(f"{CYAN}│{RESET}  modular  •  scalable        {CYAN}│{RESET}                       {GRAY}│{RESET}\n")
    stream.write(f"{CYAN}│{RESET}  transparent  •  automated   {CYAN}│{RESET}                       {GRAY}v{RESET}\n")
    stream.write(f"{CYAN}└──────────────┬───────────────┘{RESET}        {BLUE}┌──────────────────────────────┐{RESET}\n")

    stream.write(f"               {GRAY}│{RESET}                        {BLUE}│{WHITE}{BOLD}  4. VISUALIZATION            {BLUE}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                        {BLUE}│{RESET}  Clean plots and maps        {BLUE}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                        {BLUE}│{RESET}  Figures ready for papers    {BLUE}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                        {BLUE}│{RESET}  Dashboards when needed      {BLUE}│{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                        {BLUE}└──────────────┬───────────────┘{RESET}\n")
    stream.write(f"               {GRAY}│{RESET}                                       {GRAY}│{RESET}\n")
    stream.write(f"               {GRAY}└──────────────────┬────────────────────┘{RESET}\n")
    stream.write(f"                                  {GRAY}v{RESET}\n")

    stream.write(f"{NAVY}┌───────────────────────────────────────────────────────────────┐{RESET}\n")
    stream.write(f"{NAVY}│{WHITE}{BOLD}  5. REPORTING & OUTPUTS                                       {NAVY}│{RESET}\n")
    stream.write(f"{NAVY}│{RESET}  Markdown reports  |  HTML summaries  |  PDF figures          {NAVY}│{RESET}\n")
    stream.write(f"{NAVY}│{RESET}  PowerPoint slides |  PNG/SVG plots   |  ZIP result bundles   {NAVY}│{RESET}\n")
    stream.write(f"{NAVY}└───────────────────────────────────────────────────────────────┘{RESET}\n\n")

    stream.write(
        f"{GREEN}{BOLD}✓{RESET} Reproducible  "
        f"{GREEN}{BOLD}✓{RESET} Modular  "
        f"{GREEN}{BOLD}✓{RESET} Publication-ready  "
        f"{GREEN}{BOLD}✓{RESET} Open and extensible\n"
    )

    if show_asset_hint:
        stream.write("\nFull graphical abstract: fastmdx splash --open\n")

    try:
        stream.flush()
    except Exception:
        pass
