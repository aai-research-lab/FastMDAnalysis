"""Every module the package imports must be declared somewhere.

Three dependencies have reached users undeclared: scipy, arriving through
scikit-learn; netCDF4, which mdtraj declares but conda-forge does not install;
and Pillow, arriving through python-pptx. Each worked until it didn't, and
each was found by accident rather than by a check.

This test is that check. It compares what the source actually imports against
what pyproject.toml and the conda recipe declare, so an undeclared dependency
fails here rather than in somebody's environment.
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PACKAGE = REPO / "src" / "fastmdxplora"

#: Import name -> distribution name, where they differ.
DISTRIBUTION = {
    "yaml": "pyyaml",
    "sklearn": "scikit-learn",
    "pptx": "python-pptx",
    "PIL": "pillow",
    "openff": "openff-toolkit",
    "openmmplumed": "openmm-plumed",
    "umap": "umap-learn",
    "netCDF4": "netcdf4",
}

#: Imports guarded by try/except with an actionable message, and reachable
#: only when the user asks for the feature. These need not be installed.
OPTIONAL = {
    "umap-learn",      # one of three dimensionality-reduction methods
    "openmm-plumed",   # enhanced sampling, off unless asked for
}


def _imported_distributions() -> dict[str, set[str]]:
    found: dict[str, set[str]] = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name.split(".")[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                names = [node.module.split(".")[0]]
            else:
                continue
            for name in names:
                if name in sys.stdlib_module_names or name == "fastmdxplora":
                    continue
                dist = DISTRIBUTION.get(name, name)
                found.setdefault(dist, set()).add(path.name)
    return found


def _declared_in_pyproject() -> set[str]:
    text = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    section = text[text.index("dependencies = ["):text.index("[project.urls]")]
    # Entries appear as "name", "name>=1.2", or "name[extra]>=1.2".
    return {
        m.group(1).split("[")[0].lower()
        for m in re.finditer(r'"([A-Za-z0-9_.\-\[\]]+?)(?:[><=!~,\s]|")', section)
    }


def test_every_import_is_declared_in_pyproject() -> None:
    declared = _declared_in_pyproject()
    missing = {
        dist: sorted(files)
        for dist, files in _imported_distributions().items()
        if dist.lower() not in declared and dist not in OPTIONAL
    }
    assert not missing, (
        "these are imported but declared nowhere in pyproject.toml: "
        f"{missing}. Add them to dependencies, or to an extra if the import "
        "is guarded and the feature is optional."
    )


def test_the_conda_recipe_matches_pyproject() -> None:
    """A conda install must be able to do what a pip install can.

    The conda package exists so that one command gives a working four-phase
    pipeline, so anything the source imports unguarded belongs in it.
    """
    recipe_path = REPO / "recipes" / "fastmdxplora" / "recipe.yaml"
    if not recipe_path.is_file():  # pragma: no cover - recipe is optional
        return
    recipe = recipe_path.read_text(encoding="utf-8")
    run = re.search(r"  run:\n(.*?)\n\ntests:", recipe, re.S)
    assert run, "could not read the run requirements from the recipe"
    declared = {
        line.strip().lstrip("- ").split()[0].lower()
        for line in run.group(1).splitlines()
        if line.strip().startswith("- ")
    }
    # conda splits matplotlib to avoid pulling a GUI toolkit.
    aliases = {"matplotlib": "matplotlib-base"}

    missing = sorted(
        dist for dist in _imported_distributions()
        if dist not in OPTIONAL
        and aliases.get(dist, dist).lower() not in declared
    )
    assert not missing, (
        f"imported but absent from the conda recipe: {missing}. A conda "
        "install would fail at the point of use."
    )
