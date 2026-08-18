"""The manifest has to say which stack produced the numbers.

Written after a discrepancy that could be demonstrated and not explained:
the 3PTB benchmark's interaction table held 292 pair rows from the release
container and 73 when the same trajectory was re-analysed elsewhere, with
identical settings and identical ligand chemistry, and with the
container's own commit giving the second answer. Everything needed to
localise it was recorded except the one thing that decided it.
"""

from __future__ import annotations

import pytest

from fastmdxplora.provenance import (
    RESULT_BEARING_PACKAGES,
    environment_record,
)


class TestWhatIsRecorded:
    def test_the_packages_that_decide_a_number_are_covered(self):
        """Not every dependency: the ones that change a result.

        mdtraj measures the contacts, RDKit perceives the chemistry,
        OpenMM integrates. A version change in any of them can move a
        reported number without moving a line of this project's code.
        """
        for package in ("mdtraj", "openmm", "rdkit", "numpy"):
            assert package in RESULT_BEARING_PACKAGES

    def test_loaded_packages_report_their_version(self):
        import mdtraj  # noqa: F401
        import numpy  # noqa: F401

        record = environment_record()
        assert record["mdtraj"] and record["mdtraj"] != "loaded"
        assert record["numpy"]

    def test_the_interpreter_and_platform_are_there(self):
        record = environment_record()
        assert record["python"].count(".") == 2
        assert record["platform"]

    def test_an_absent_package_is_none_rather_than_missing(self):
        """"This run did not use RDKit" and "this manifest forgot to say"
        are different statements, and a reader cannot tell them apart if
        the key is simply absent."""
        record = environment_record()
        for package in RESULT_BEARING_PACKAGES:
            assert package in record, f"{package} is not reported at all"

    def test_it_reads_what_ran_not_what_could_run(self, monkeypatch):
        """A subprocess probe answers a different question.

        What could be imported is a fact about the machine; what was
        imported is a fact about the run, and only the second explains a
        number.
        """
        import sys

        monkeypatch.setitem(sys.modules, "propka", None)
        assert environment_record()["propka"] is None


class TestItReachesTheManifest:
    def test_the_orchestrator_writes_it(self):
        import inspect
        from pathlib import Path

        from fastmdxplora import orchestrator

        source = Path(inspect.getfile(orchestrator)).read_text(
            encoding="utf-8")
        writer = source[source.index("def _write_manifest"):]
        writer = writer[:writer.index("\n    def ", 1)]

        assert "environment_record()" in writer, (
            "the manifest records the software version and must record the "
            "stack underneath it")
        assert '"environment"' in writer


class TestTheMetadataFallback:
    """Some packages carry no version attribute at all.

    PDBFixer is one, and it decides every protonation state in a run. The
    installed distribution knows even where the module is silent, so the
    record asks it rather than settling for "loaded".
    """

    def test_a_module_without_an_attribute_is_named_anyway(self, monkeypatch):
        import sys
        import types

        # A module that is loaded, carries no `__version__`, and whose
        # distribution is installed: exactly PDBFixer's shape.
        stand_in = types.ModuleType("scipy")
        monkeypatch.setitem(sys.modules, "scipy", stand_in)

        record = environment_record()
        assert record["scipy"] not in (None, "loaded"), (
            "the distribution's metadata should have supplied a version "
            "where the module did not")

    def test_loaded_but_unknown_is_still_distinguishable(self, monkeypatch):
        """A module with neither attribute nor installed distribution.

        "Loaded, version unknown" and "not loaded" say different things
        about a run, so they must not collapse onto each other.
        """
        import sys
        import types

        from fastmdxplora import provenance

        # A name with no attribute and no installed distribution. Chosen
        # rather than borrowed: the first draft used propka, which is
        # installed here, so the fallback found 3.5.1 and the test was
        # asserting the opposite of what it meant.
        monkeypatch.setattr(
            provenance, "RESULT_BEARING_PACKAGES",
            ("fastmdx_not_a_real_package", "openmmplumed"))
        monkeypatch.setitem(
            sys.modules, "fastmdx_not_a_real_package",
            types.ModuleType("fastmdx_not_a_real_package"))
        record = provenance.environment_record()

        assert record["fastmdx_not_a_real_package"] == "loaded"
        assert record["openmmplumed"] is None
