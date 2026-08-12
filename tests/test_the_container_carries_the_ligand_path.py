"""The image exists because wheels cannot reach the ligand path.

The OpenFF toolkit and openmmforcefields are distributed through conda-forge
only. They are not on PyPI, so a machine that cannot reach conda-forge cannot
have them however the install is arranged -- and FastMDXplora on such a
machine prepares proteins perfectly well and refuses the moment a ligand needs
parameterising. That is a correct refusal and a real limit.

A container carries the whole stack in one file, which is the thing that can
cross an airgap. What this file checks is that the definition will actually
produce such an image, and that it fails rather than shipping one that will
not.
"""

from __future__ import annotations

from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DEFINITION = ROOT / "container" / "fastmdx.def"
WORKFLOW = ROOT / ".github" / "workflows" / "container.yml"


@pytest.fixture(scope="module")
def definition() -> str:
    if not DEFINITION.is_file():
        pytest.skip("no container definition")
    return DEFINITION.read_text(encoding="utf-8")


class TestItCarriesWhatWheelsCannot:
    @pytest.mark.parametrize("package", [
        "openff-toolkit", "openmmforcefields",
    ])
    def test_the_conda_only_packages_are_installed(
            self, definition: str, package: str) -> None:
        """The whole reason the image exists."""
        assert package in definition

    @pytest.mark.parametrize("package", [
        "openmm", "openmm-plumed", "mdtraj", "pdbfixer",
    ])
    def test_the_rest_of_the_stack_is_there_too(
            self, definition: str, package: str) -> None:
        assert package in definition

    def test_they_are_solved_together(self, definition: str) -> None:
        """Separate installs are how you get an OpenMM that does not match
        the toolkit built against it."""
        install = definition[definition.index("micromamba install"):]
        block = install[:install.index("micromamba clean")]
        for package in ("openmm", "openff-toolkit", "openmmforcefields"):
            assert package in block


class TestItRefusesToShipAnImageThatWillNotWork:
    def test_the_build_asserts_the_ligand_path_imports(
            self, definition: str) -> None:
        assert "import openff.toolkit" in definition

    def test_the_build_asserts_a_cuda_platform(self, definition: str) -> None:
        """An image that silently lacks CUDA runs everywhere it is taken, on
        the CPU, and nothing says why it is slow."""
        assert 'assert "CUDA" in platforms' in definition

    def test_a_failure_says_what_was_found(self, definition: str) -> None:
        assert "getPluginLoadFailures" in definition


class TestItBringsItsOwnPlugins:
    def test_the_plugin_directory_points_inside(self, definition: str) -> None:
        """A host with a site OpenMM will have OPENMM_PLUGIN_DIR set to a
        path that does not exist in the image."""
        assert "OPENMM_PLUGIN_DIR=/opt/conda/lib/plugins" in definition


@pytest.fixture(scope="module")
def workflow() -> str:
    if not WORKFLOW.is_file():
        pytest.skip("no container workflow")
    return WORKFLOW.read_text(encoding="utf-8")


class TestItIsBuiltWhereItCanBe:
    def test_it_builds_on_a_tag(self, workflow: str) -> None:
        assert 'tags:' in workflow and '"v*"' in workflow

    def test_it_waits_for_the_release_to_reach_conda_forge(
            self, workflow: str) -> None:
        """The image installs from conda-forge, so it is the feedstock that
        has to have caught up -- not PyPI. An image built too early installs
        the previous version and does not say so."""
        assert "api.anaconda.org/package/conda-forge/fastmdxplora" in workflow

    def test_it_installs_from_conda_forge_not_pip(self, definition: str) -> None:
        """One solver. Installing the stack with conda and FastMDXplora with
        pip on top re-resolves its dependencies against what conda already
        put there, and pip does not read conda's metadata properly: it can
        replace a conda-installed numpy with a wheel."""
        assert "pip install" not in definition
        assert "fastmdxplora=" in definition

    def test_it_can_be_run_by_hand(self, workflow: str) -> None:
        """So an image can be rebuilt without cutting a release."""
        assert "workflow_dispatch" in workflow


class TestTheDocumentationSaysWhyItExists:
    def test_it_names_the_reason(self) -> None:
        page = (ROOT / "docs" / "remote.md").read_text(encoding="utf-8")
        assert "conda-forge only" in page
        assert "not on PyPI" in page

    def test_it_says_to_pass_the_gpu_through(self) -> None:
        page = (ROOT / "docs" / "remote.md").read_text(encoding="utf-8")
        assert "--nv" in page
