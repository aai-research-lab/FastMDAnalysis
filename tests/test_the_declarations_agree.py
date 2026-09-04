"""Three files declare this package's dependencies, and they drift.

`pyproject.toml`, `environment.yml` and `recipes/fastmdxplora/recipe.yaml`
each state a floor for the same packages, and each has at some point been
corrected on its own. The OpenMM floor is the recorded case: pyproject
explains at length why `>=8.0` is below the real one and corrects it to
`>=8.2`, and the other two kept saying 8.0 for four releases.

CI cannot catch a floor. All nine legs resolve the newest version of
everything, so no leg has ever run at a version this package claims to
support. That makes these declarations the only statement of the range,
and a wrong one is not visible until somebody outside the lab installs it.

Same species as the `cuda-version` lesson in the TODO: a constraint stated
in one place and enforced nowhere travels to the next environment as a
comment, not as a rule.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PYPROJECT = ROOT / "pyproject.toml"
ENVIRONMENT = ROOT / "environment.yml"
RECIPE = ROOT / "recipes" / "fastmdxplora" / "recipe.yaml"


def _pyproject() -> dict:
    return tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))


def _floor(text: str, package: str) -> str | None:
    """The version floor declared for `package` anywhere in a text file."""
    # Covers a TOML list entry (`"numpy>=2.0",`), a conda list entry
    # (`  - numpy>=2.0`) and a recipe entry (`    - numpy >=2.0`).
    match = re.search(
        rf"""^\s*(?:-\s*)?["']?{re.escape(package)}\s*>=\s*([0-9][0-9.]*)""",
        text, re.MULTILINE)
    return match.group(1) if match else None


class TestTheExtrasTableIsClosed:
    """`[all]` may only name extras that exist."""

    def test_all_names_only_real_extras(self) -> None:
        extras = _pyproject()["project"]["optional-dependencies"]
        named: set[str] = set()
        for requirement in extras["all"]:
            for group in re.findall(r"\[([^\]]*)\]", requirement):
                named.update(part.strip() for part in group.split(","))

        undefined = sorted(named - set(extras))

        # pip prints "does not provide the extra 'plumed'" and installs
        # anyway, so this failed as a warning for as long as it existed.
        assert not undefined, (
            f"[all] names extras that do not exist: {undefined}. pip warns "
            f"and continues, so the install silently lacks them."
        )

    def test_plumed_is_not_offered_as_an_extra(self) -> None:
        """Stated as its own assertion because the comment above the table
        says so and the table did not."""
        extras = _pyproject()["project"]["optional-dependencies"]
        assert "plumed" not in extras
        assert not any("plumed" in req for req in extras["all"])


class TestNumpyIsDeclaredHighEnoughToRun:

    def test_the_floor_admits_np_trapezoid(self) -> None:
        """`np.trapezoid` is the NumPy 2.0 rename of `np.trapz`.

        `binding_free_energy()` calls it, so any NumPy the declarations admit
        must have it. The floor said 1.22 in all three files.
        """
        source = (ROOT / "src" / "fastmdxplora" / "simulation"
                  / "binding.py").read_text(encoding="utf-8")
        if "np.trapezoid" not in source:
            pytest.skip("binding.py no longer calls np.trapezoid")

        floor = _floor(PYPROJECT.read_text(encoding="utf-8"), "numpy")
        assert floor is not None
        assert int(floor.split(".")[0]) >= 2, (
            f"binding.py calls np.trapezoid, which needs NumPy 2.0, but the "
            f"declared floor is {floor}."
        )


class TestTheThreeFilesAgree:

    @pytest.mark.parametrize("package", ["numpy", "openmm"])
    def test_the_floors_match(self, package: str) -> None:
        floors = {
            "pyproject.toml": _floor(
                PYPROJECT.read_text(encoding="utf-8"), package),
            "environment.yml": _floor(
                ENVIRONMENT.read_text(encoding="utf-8"), package),
            "recipe.yaml": _floor(RECIPE.read_text(encoding="utf-8"), package),
        }
        stated = {k: v for k, v in floors.items() if v is not None}
        if len(stated) < 2:
            pytest.skip(f"{package} is declared in fewer than two files")

        assert len(set(stated.values())) == 1, (
            f"{package} floors disagree: {stated}. A floor corrected in one "
            f"file and not the others is how openmm stayed at 8.0 in two of "
            f"three for four releases."
        )
