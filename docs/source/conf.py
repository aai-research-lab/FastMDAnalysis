"""Sphinx configuration for FastMDAnalysis documentation."""
from __future__ import annotations

import os
import sys
from datetime import datetime

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

# Add src to path for autodoc
SRC_DIR = os.path.join(PROJECT_ROOT, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

project = "FastMDAnalysis"
copyright = f"{datetime.now():%Y}, A. Aina and D. Kwan"
author = "A. Aina and D. Kwan"

release = "0.1.0"
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary", 
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "special-members": "__init__",
}

autodoc_mock_imports = [
    "mdtraj",
    "numpy", 
    "matplotlib",
    "matplotlib.pyplot",
    "sklearn",
    "scipy",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "collapse_navigation": False,
    "navigation_depth": 4,
}
primary_domain = "py"

# Auto-generate API documentation on build
def setup(app):
    from sphinx.ext import apidoc
    
    def run_apidoc(_):
        apidoc_main = getattr(apidoc, 'main', None)
        if apidoc_main is None:
            # Fallback for older Sphinx versions
            from sphinx.ext.apidoc import main as apidoc_main
        
        apidoc_main([
            '-f',  # Force overwrite
            '-e',  # Put each module/class on its own page
            '-M',  # Put module documentation before submodule
            '-o', 
            os.path.join(os.path.dirname(__file__)),  # Output to source/
            os.path.join(PROJECT_ROOT, 'src', 'fastmdanalysis')  # Package path
        ])
    
    app.connect('builder-inited', run_apidoc)
