# FastMDAnalysis/src/fastmdanalysis/utils/__init__.py
from .io import load_trajectory
from .slideshow import gather_figures, slide_show
from .options import OptionsForwarder, forward_options, apply_alias_mapping
from .plotting import auto_ticks, apply_slide_style

__all__ = [
    "load_trajectory",
    "gather_figures",
    "slide_show",
    "OptionsForwarder",
    "forward_options",
    "apply_alias_mapping",
    "auto_ticks",
    "apply_slide_style",
]

