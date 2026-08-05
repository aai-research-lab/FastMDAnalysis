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

import pytest

# The check reads source files, so its verdict is the same on every Python.
# `sys.stdlib_module_names` arrived in 3.10; rather than carry a duplicate
# list of standard-library modules for 3.9, run the check where the
# interpreter can answer the question itself.
pytestmark = pytest.mark.skipif(
    not hasattr(sys, "stdlib_module_names"),
    reason="needs sys.stdlib_module_names (Python 3.10+); result is version-independent",
)

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
}

#: Imported, and deliberately absent from pyproject.toml because pip cannot
#: install them at all. Naming one in an extra does not make it reachable; it
#: makes the whole extra fail, which is what `pip install
#: "fastmdxplora[ligand]"` did. They are declared in the conda recipe, which
#: is where they can be had, and the import raises with that instruction.
CONDA_ONLY = {
    "openff-toolkit": "no PyPI distribution; conda-forge only",
    "openmm-plumed": "no PyPI distribution; conda-forge only",
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
        if dist.lower() not in declared
        and dist not in OPTIONAL
        and dist not in CONDA_ONLY
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


def test_the_recipe_says_where_the_package_is_actually_built() -> None:
    """The recipe in this repository is a reference copy; the conda-forge
    feedstock builds the package.

    A dependency added to one and not the other gives a conda install that
    fails at the point of use -- which is what happened when WeasyPrint was
    added to the PyPI extras and nowhere else. The relationship is written in
    the recipe so the next release does not have to rediscover it.
    """
    from pathlib import Path

    recipe = Path(__file__).resolve().parents[1] / "recipes" / "fastmdxplora" / "recipe.yaml"
    if not recipe.is_file():  # pragma: no cover - the recipe is optional
        return
    text = recipe.read_text(encoding="utf-8")
    assert "feedstock" in text, (
        "the recipe should say that the feedstock is what builds the package"
    )


def test_a_conda_only_dependency_is_still_in_the_recipe() -> None:
    """Something pip cannot install has to be gettable somewhere.

    `pip install "fastmdxplora[ligand]"` named openff-toolkit, which has no
    PyPI distribution, so the command failed outright rather than installing
    what pip could. It is out of the extra and in the recipe, and this checks
    it did not simply vanish.
    """
    recipe_path = REPO / "recipes" / "fastmdxplora" / "recipe.yaml"
    if not recipe_path.is_file():  # pragma: no cover - the recipe is optional
        return
    recipe = recipe_path.read_text(encoding="utf-8")
    for name in CONDA_ONLY:
        assert name in recipe, (
            f"{name} cannot be installed by pip and is not in the recipe "
            "either, so there is nowhere to get it"
        )


def test_the_reason_is_recorded_for_each() -> None:
    """So the list cannot become a place to hide an undeclared import."""
    assert CONDA_ONLY
    for name, reason in CONDA_ONLY.items():
        assert reason and len(reason) > 10, name


def test_the_alias_package_carries_the_released_version() -> None:
    """`fastmdx` is a thin alias for `fastmdxplora`, and its version is
    written in a file where the main package takes its own from the git tag.

    A hand-written version beside a derived one drifts, and this one did: the
    release workflow refused a 2.3.0 tag because the alias still said 2.2.0.
    That check fires after tagging, which is too late to be comfortable, so
    this one fires in the ordinary test run -- against the changelog, which is
    the thing a release is cut from.
    """
    import re
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    shim = repo / "shim-package" / "pyproject.toml"
    if not shim.is_file():  # pragma: no cover - the alias is optional
        return

    written = re.search(r'^version = "([^"]+)"', shim.read_text(encoding="utf-8"),
                        re.M)
    assert written, "the alias declares no version"

    changelog = (repo / "CHANGELOG.md").read_text(encoding="utf-8")
    released = re.search(r"^## \[(\d+\.\d+\.\d+)\]", changelog, re.M)
    assert released, "the changelog names no released version"

    assert written.group(1) == released.group(1), (
        f"the alias says {written.group(1)} and the changelog's latest "
        f"release is {released.group(1)}. Bump shim-package/pyproject.toml "
        "before tagging, or the release workflow will refuse the tag."
    )


def test_info_reports_every_backend_the_software_reaches_for() -> None:
    """`fastmdx info` exists to say what will work on this machine.

    It listed two backends where the software reaches for eight, so a PyPI
    install reported both present and said nothing about the OpenFF toolkit --
    which a protein-ligand setup needs and which pip cannot install. Somebody
    reading it would conclude their install was complete and find out
    otherwise three phases later.
    """
    from fastmdxplora.cli.main import _BACKENDS

    reported = {
        import_name
        for _group, backends in _BACKENDS
        for _display, import_name, _hint in backends
    }
    for needed in ("openmm", "pdbfixer", "openff.toolkit", "openmmforcefields",
                   "rdkit", "propka", "weasyprint", "markdown"):
        assert needed in reported, f"info does not mention {needed}"


def test_each_backend_says_where_to_get_it() -> None:
    """And says something that works: naming a pip extra for the OpenFF
    toolkit sent people to a command that could not succeed."""
    from fastmdxplora.cli.main import _BACKENDS

    for _group, backends in _BACKENDS:
        for display, import_name, hint in backends:
            assert hint, display
            if import_name == "openff.toolkit":
                assert "conda" in hint and "pip" not in hint, (
                    "the toolkit is not on PyPI, so a pip hint cannot work"
                )


def test_backends_are_grouped_by_what_they_are_for() -> None:
    """A trajectory analysis does not care that OpenMM is absent, and saying
    so unqualified reads as a broken install."""
    from fastmdxplora.cli.main import _BACKENDS

    groups = {group for group, _backends in _BACKENDS}
    assert len(groups) >= 3
    assert all(group and group[0].islower() for group in groups), (
        "the groups read as a phrase after the heading"
    )


def test_a_backend_is_probed_in_a_process_of_its_own() -> None:
    """Importing is the only honest test of whether a backend will load, and
    the failure is not always quiet or catchable.

    WeasyPrint without Pango raises from the dynamic loader and prints five
    lines about its installation guide on the way. Redirecting Python's stderr
    did not stop them and neither did redirecting the file descriptor -- and
    the test written for that fix passed only because it faked the message
    with a print, which checks the assumption rather than the behaviour.

    A subprocess ends the argument: its output is captured whatever writes it,
    and a backend that fails hard cannot take this command with it.
    """
    import inspect

    # `from fastmdxplora.cli import main` gives the entry-point function,
    # which the package re-exports; the module is the longer path.
    from fastmdxplora.cli.main import _probe_backends

    source = inspect.getsource(_probe_backends)
    assert "subprocess.run" in source
    assert "capture_output=True" in source


def test_it_reports_what_is_there_and_what_is_not() -> None:
    from fastmdxplora.cli.main import _probe_backends

    found = _probe_backends(("json", "no_such_module_at_all"))
    assert found["json"][0] == "installed"
    assert found["no_such_module_at_all"][0] == "missing"


def test_a_backend_that_raises_on_import_is_reported_broken(tmp_path) -> None:
    """The case a real machine hits: installed, and will not load."""
    import os
    import sys

    from fastmdxplora.cli.main import _probe_backends

    (tmp_path / "pretend_backend.py").write_text(
        "import sys\n"
        "sys.stderr.write('a library complaining loudly\\n')\n"
        "raise OSError(\"cannot load library 'libgobject-2.0-0'\")\n",
        encoding="utf-8")

    original = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = os.pathsep.join(
        [str(tmp_path)] + ([original] if original else []))
    try:
        state, detail = _probe_backends(("pretend_backend",))["pretend_backend"]
    finally:
        if original:
            os.environ["PYTHONPATH"] = original
        else:
            os.environ.pop("PYTHONPATH", None)

    assert state == "broken"
    assert "libgobject" in detail


def test_the_probe_cannot_be_interrupted_by_what_it_probes(capfd) -> None:
    """The point of the subprocess: whatever the backend writes, and however
    far down it writes it, stays out of this command's output."""
    import os

    from fastmdxplora.cli.main import _probe_backends

    _probe_backends(("json",))
    captured = capfd.readouterr()
    assert "complaining" not in captured.err



def test_a_probe_that_cannot_run_still_answers(capsys) -> None:
    """A table of "unknown" is worse than a noisy answer.

    The subprocess probe failed on a real machine and reported nothing about
    every backend, because the fallback swallowed the reason -- which is the
    defect this software exists to avoid, in the code written to avoid it. It
    now says why the subprocess could not answer and checks in this process
    instead.
    """
    import importlib
    import sys

    module = importlib.import_module("fastmdxplora.cli.main")
    real = sys.executable
    sys.executable = "/nonexistent/python"
    try:
        found = module._probe_backends(("json", "no_such_module_at_all"))
    finally:
        sys.executable = real

    assert found["json"][0] == "installed"
    assert found["no_such_module_at_all"][0] == "missing"
    assert "checked in this process" in capsys.readouterr().out


def test_the_reason_the_subprocess_failed_is_named(capsys) -> None:
    """Reporting that something could not be checked, without saying why,
    leaves nothing to act on."""
    import importlib
    import sys

    module = importlib.import_module("fastmdxplora.cli.main")
    real = sys.executable
    sys.executable = "/nonexistent/python"
    try:
        module._probe_backends(("json",))
    finally:
        sys.executable = real

    printed = capsys.readouterr().out
    assert "FileNotFoundError" in printed or "No such file" in printed


def test_a_backend_writing_to_stdout_cannot_corrupt_the_answer(tmp_path) -> None:
    """WeasyPrint writes its five-line complaint to stdout, not stderr.

    That took four attempts to notice. Redirecting Python's stderr did not
    quiet it, redirecting the file descriptor did not either, and moving the
    probe to a subprocess then failed to parse its own output -- because the
    answer was coming back on the channel the message was using. It travels by
    file now, so nothing a backend prints, on any stream, can reach it.
    """
    import importlib
    import os

    module = importlib.import_module("fastmdxplora.cli.main")

    (tmp_path / "loud_backend.py").write_text(
        "import sys\n"
        "sys.stdout.write('\\n-----\\nsomething loud on stdout\\n-----\\n')\n"
        "sys.stdout.flush()\n"
        "raise OSError(\"cannot load library 'libgobject-2.0-0'\")\n",
        encoding="utf-8")

    original = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = os.pathsep.join(
        [str(tmp_path)] + ([original] if original else []))
    try:
        found = module._probe_backends(("loud_backend", "json"))
    finally:
        if original:
            os.environ["PYTHONPATH"] = original
        else:
            os.environ.pop("PYTHONPATH", None)

    assert found["loud_backend"][0] == "broken"
    assert "libgobject" in found["loud_backend"][1]
    assert found["json"][0] == "installed", (
        "the noise must not take the other answers with it"
    )


def test_the_answer_does_not_come_back_on_stdout() -> None:
    """Which is the channel a failing backend is most likely to use."""
    import importlib
    import inspect

    module = importlib.import_module("fastmdxplora.cli.main")
    source = inspect.getsource(module._probe_backends)
    assert "sys.stdout.write(json.dumps" not in source
    assert "open(sys.argv[2]" in source
