"""What went wrong when a simulation failed, read from the state it failed in.

A run ending with "the coordinates are not finite" has told you almost
nothing, and the advice that follows -- lower the timestep, lower the
temperature, raise the friction -- is a list of things that sometimes help,
offered without knowing which applies. Sometimes none of them do, because the
problem is a ligand whose parameters are wrong.
"""

from __future__ import annotations

import numpy as np
import pytest


def _system(residues, atoms_per_residue=4):
    import mdtraj as md

    top = md.Topology()
    chain = top.add_chain()
    for name in residues:
        residue = top.add_residue(name, chain, resSeq=1)
        for index in range(atoms_per_residue):
            top.add_atom(f"C{index}", md.element.carbon, residue)
    return top


def _blow_up(topology, names):
    coordinates = np.zeros((topology.n_atoms, 3))
    for atom in topology.atoms:
        if atom.residue.name in names:
            coordinates[atom.index] = np.nan
    return coordinates


class TestItSaysWhichCauseTheEvidencePointsAt:
    """The four cases need four different remedies, and the same list for all
    of them is the same as no diagnosis."""

    def test_a_ligand_alone_points_at_the_ligand(self) -> None:
        from fastmdxplora.simulation.diagnose import diagnose_failure

        topology = _system(["ALA"] * 20 + ["BNZ"])
        found = diagnose_failure(topology, _blow_up(topology, {"BNZ"}),
                                 stage="NVT")
        assert "ligand" in found.likely_cause
        assert any("SDF" in item or "chemistry" in item for item in found.advice)

    def test_and_says_a_smaller_timestep_will_not_help(self) -> None:
        """Because it will not: no timestep is small enough to fix a wrong
        parameter."""
        from fastmdxplora.simulation.diagnose import diagnose_failure

        topology = _system(["ALA"] * 20 + ["BNZ"])
        found = diagnose_failure(topology, _blow_up(topology, {"BNZ"}),
                                 stage="NVT")
        assert "No timestep is small enough" in found.likely_cause
        assert not any("timestep" in item for item in found.advice)

    def test_lipids_point_at_the_packing(self) -> None:
        from fastmdxplora.simulation.diagnose import diagnose_failure

        topology = _system(["ALA"] * 20 + ["POP"] * 10)
        found = diagnose_failure(topology, _blow_up(topology, {"POP", "ALA"}),
                                 stage="NPT")
        assert "bilayer" in found.likely_cause
        assert any("restrain" in item.lower() for item in found.advice)

    def test_the_whole_system_points_at_the_integration(self) -> None:
        """Which is the case the usual advice is actually for."""
        from fastmdxplora.simulation.diagnose import diagnose_failure

        topology = _system(["ALA"] * 20 + ["HOH"] * 20)
        found = diagnose_failure(topology, _blow_up(topology, {"ALA", "HOH"}),
                                 stage="NVT")
        assert "integration failure" in found.likely_cause
        assert any("timestep" in item for item in found.advice)

    def test_one_residue_kind_points_somewhere_local(self) -> None:
        from fastmdxplora.simulation.diagnose import diagnose_failure

        topology = _system(["ALA"] * 20 + ["TRP"] * 2)
        found = diagnose_failure(topology, _blow_up(topology, {"TRP"}),
                                 stage="minimization")
        assert "One residue type" in found.likely_cause
        assert "TRP" in found.likely_cause

    def test_where_it_points_nowhere_it_says_so(self) -> None:
        """More useful than a confident list that happens not to apply."""
        from fastmdxplora.simulation.diagnose import diagnose_failure

        topology = _system(["ALA"] * 40 + ["SER"] * 10 + ["HOH"] * 5)
        coordinates = np.zeros((topology.n_atoms, 3))
        for atom in topology.atoms:
            if atom.index % 17 == 0:
                coordinates[atom.index] = np.nan
        found = diagnose_failure(topology, coordinates, stage="NVT")
        assert "do not fall into a pattern" in found.likely_cause
        assert "may not apply" in found.likely_cause


class TestWhatItReports:
    def test_it_names_what_went_wrong_and_where(self) -> None:
        from fastmdxplora.simulation.diagnose import diagnose_failure

        topology = _system(["ALA"] * 20 + ["BNZ"])
        text = diagnose_failure(topology, _blow_up(topology, {"BNZ"}),
                                stage="NVT equilibration").as_text()
        assert "NVT equilibration" in text
        assert "BNZ" in text
        assert "What to try:" in text

    def test_a_finite_state_is_reported_as_an_energy_problem(self) -> None:
        """The failure was in the energy rather than the positions -- a force
        that is large but not yet infinite."""
        from fastmdxplora.simulation.diagnose import diagnose_failure

        topology = _system(["ALA"] * 10)
        found = diagnose_failure(topology, np.zeros((topology.n_atoms, 3)),
                                 stage="minimization")
        assert found.n_bad_atoms == 0
        assert "energy rather than the positions" in found.likely_cause


class TestItDoesNotRetry:
    """A run that exploded because its ligand is wrong will explode again more
    slowly at half the timestep, and a rescue that produces a trajectory from
    a broken system is worse than a failure: the failure is visible.
    """

    def test_nothing_here_reruns_anything(self) -> None:
        import inspect

        from fastmdxplora.simulation import diagnose

        source = inspect.getsource(diagnose)
        for attempt in ("simulation.step", "retry", "rerun", "minimizeEnergy"):
            assert attempt not in source, (
                f"a diagnosis should report, not {attempt}"
            )

    def test_and_the_module_says_why(self) -> None:
        import inspect

        from fastmdxplora.simulation import diagnose

        assert "Nothing here retries" in inspect.getdoc(diagnose)


class TestTheRunnerUsesIt:
    def test_a_failed_state_is_diagnosed(self) -> None:
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner._validate_state_finite)
        assert "topology=simulation.topology" in source
        assert "positions=positions" in source

    def test_and_a_diagnosis_that_fails_does_not_hide_the_failure(self) -> None:
        """The failure worth reporting is the simulation's, not the
        diagnosis's."""
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner._validation_error)
        assert "fall through to the general message" in source
