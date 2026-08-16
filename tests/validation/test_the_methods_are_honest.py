"""Checks on the enhanced-sampling cases themselves, needing no simulation.

A case that claims to exercise a method and does not is worse than none: it
reads as evidence. These run on every ordinary `pytest`.
"""

from __future__ import annotations

import pytest

from .methods import COMMON_SIMULATION, INAPPLICABLE, METHODS, TRIPEPTIDE, names


class TestTheCasesAreWellFormed:
    @pytest.mark.parametrize("method", METHODS, ids=lambda m: m.name)
    def test_each_names_a_biasing_block(self, method) -> None:
        simulation = method.config["simulation"]
        biased = {"umbrella", "metadynamics", "steered"} & set(simulation)
        assert len(biased) == 1, (
            f"{method.name} configures {biased or 'no'} biasing")

    @pytest.mark.parametrize("method", METHODS, ids=lambda m: m.name)
    def test_each_promises_something(self, method) -> None:
        assert method.produces, f"{method.name} asserts nothing"

    @pytest.mark.parametrize("method", METHODS, ids=lambda m: m.name)
    def test_the_settings_are_accepted(self, method) -> None:
        """A misspelt key is accepted and ignored by anything that does not
        check, which is how `minimum_ovelap: 0.15` stitched a study at the
        three per cent default while its author believed otherwise."""
        simulation = method.config["simulation"]
        if "umbrella" in simulation:
            from fastmdxplora.simulation.umbrella import check_umbrella_keys

            check_umbrella_keys(simulation["umbrella"])
        if "metadynamics" in simulation:
            from fastmdxplora.simulation.metadynamics import METADYNAMICS_KEYS

            unknown = set(simulation["metadynamics"]) - set(METADYNAMICS_KEYS)
            assert not unknown, unknown

    @pytest.mark.parametrize("method", METHODS, ids=lambda m: m.name)
    def test_the_selections_resolve_on_the_peptide(self, method) -> None:
        """`resid 3` is the natural thing to write for a tripeptide numbered
        1 to 3 in its own file, and MDTraj counts from zero -- it matched no
        atoms and failed a window after the shared system had been built."""
        import tempfile
        from pathlib import Path

        import mdtraj as md

        path = Path(tempfile.mkdtemp()) / "t.pdb"
        path.write_text(TRIPEPTIDE, encoding="utf-8")
        topology = md.load(str(path)).topology

        simulation = method.config["simulation"]
        block = (simulation.get("umbrella") or simulation.get("metadynamics")
                 or simulation.get("steered"))
        for key in ("selection", "selection_a", "selection_b"):
            expression = block.get(key)
            if expression:
                assert len(topology.select(expression)), (
                    f"{method.name}: {key} {expression!r} matches no atoms")


class TestTheRunsAreShortEnoughToBeRun:
    """A suite nobody runs because it takes an hour catches nothing."""

    @pytest.mark.parametrize("method", METHODS, ids=lambda m: m.name)
    def test_no_case_is_long(self, method) -> None:
        steps = method.config["simulation"].get("production_steps", 0)
        total = steps * method.runs
        assert total <= 40_000, (
            f"{method.name} runs {total:,} steps in total, which on a "
            "fourteen-atom peptide is minutes rather than seconds")

    def test_pressure_coupling_is_off(self) -> None:
        """A barostat needs enough water for its fluctuations to average
        out. A tripeptide in a 1.45 nm box has a couple of hundred molecules,
        and the box collapsed below twice the cutoff."""
        assert COMMON_SIMULATION["npt_steps"] == 0


class TestTheGatingIsExercised:
    def test_something_is_expected_not_to_apply(self) -> None:
        """On a real protein every analysis applies, so a corpus of them
        cannot test the gates at all."""
        assert INAPPLICABLE

    @pytest.mark.parametrize("analysis", sorted(INAPPLICABLE))
    def test_the_reason_is_written_down(self, analysis) -> None:
        assert len(INAPPLICABLE[analysis]) > 60

    def test_the_peptide_is_too_short_to_fold(self) -> None:
        """Three residues, and Q needs a pair four apart in sequence."""
        assert TRIPEPTIDE.count("CA ") == 3


class TestTheMethodsAreDistinct:
    def test_each_is_named_once(self) -> None:
        assert len(set(names())) == len(METHODS)

    def test_all_three_methods_are_covered(self) -> None:
        assert set(names()) == {"umbrella", "metadynamics", "steered"}


class TestTheDocumentedDefaultsAreTheRealOnes:
    """A force field named in a docstring reaches a methods section.

    `setup.prepare` said the default was charmm36 for some time, and it
    never was: `auto` resolves to amber-openff, which is the ligand-capable
    choice. Nothing checked, because the claim was prose and the value was
    code, and the two only met in a reader's head.
    """

    def test_auto_resolves_to_what_the_docstring_says(self) -> None:
        from fastmdxplora.setup import prepare
        from fastmdxplora.setup.forcefields import (
            AUTO_FORCEFIELD,
            resolve_forcefield,
        )

        resolved = resolve_forcefield("auto")
        assert resolved.name == AUTO_FORCEFIELD

        stated = prepare.__doc__
        assert AUTO_FORCEFIELD in stated, (
            "the module docstring must name the force field `auto` gives")
        for xml in resolved.xmls:
            assert xml in stated, (
                f"{xml} is what `auto` loads and the docstring does not "
                "mention it")

    def test_no_docstring_claims_an_unregistered_default(self) -> None:
        """charmm36 is a choice, not the default, and saying otherwise sends
        someone to report parameters they did not use."""
        from fastmdxplora.setup import prepare
        from fastmdxplora.setup.forcefields import resolve_forcefield

        text = prepare.__doc__ or ""
        default_xmls = set(resolve_forcefield("auto").xmls)
        other = set(resolve_forcefield("charmm36").xmls) - default_xmls

        for xml in other:
            if xml in text:
                assert "never was" in text or "available by name" in text, (
                    f"{xml} appears in the docstring without saying it is "
                    "not the default")
