"""Putting a protein in a lipid bilayer, and the things that go wrong quietly.

A membrane protein simulated in water is not the protein: the hydrophobic belt
that sits in the bilayer is exposed to solvent and the helices splay. OpenMM
builds the bilayer itself, which removes the usual dependency on an external
packing tool -- but it does not check the two things that make the result
meaningful, and both fail silently.
"""

from __future__ import annotations

import numpy as np
import pytest


class TestTheProteinHasToBeOrientedFirst:
    """`addMembrane` places the bilayer in the xy plane and assumes the
    protein is already lying along z. A PDB entry usually is not:
    crystallographic axes have no relation to a membrane normal. Embedding one
    sideways packs lipids around it wrong, the run completes, and every number
    describes a structure nobody would recognise.
    """

    @staticmethod
    def _rod(axis: int):
        """A protein-shaped object extended along one axis."""
        import mdtraj as md

        top = md.Topology()
        chain = top.add_chain()
        coordinates = []
        for index in range(40):
            residue = top.add_residue("ALA", chain, resSeq=index + 1)
            for name in ("N", "CA", "C"):
                top.add_atom(name, md.element.carbon, residue)
            base = np.zeros(3)
            base[axis] = index * 0.35
            for offset in (0.0, 0.12, 0.24):
                point = base.copy()
                point[(axis + 1) % 3] += offset
                coordinates.append(point)
        return top.to_openmm(), np.array(coordinates)

    def test_a_protein_along_z_is_accepted(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        from fastmdxplora.setup.membrane import check_orientation

        topology, positions = self._rod(axis=2)
        assert check_orientation(topology, positions) is None

    def test_a_protein_lying_sideways_is_refused(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        from fastmdxplora.setup.membrane import check_orientation

        topology, positions = self._rod(axis=0)
        problem = check_orientation(topology, positions)
        assert problem is not None
        assert "membrane normal" in problem

    def test_the_refusal_says_where_oriented_structures_come_from(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        from fastmdxplora.setup.membrane import check_orientation

        topology, positions = self._rod(axis=0)
        problem = check_orientation(topology, positions)
        assert "opm.phar.umich.edu" in problem
        assert "membrane_orientation_checked" in problem, (
            "and how to proceed anyway"
        )

    def test_something_too_small_to_judge_is_not_refused(self) -> None:
        """The test is the shape of the protein, and a handful of atoms has
        no shape worth reading."""
        pytest.importorskip("openmm", reason="requires the [md] extra")

        import mdtraj as md

        from fastmdxplora.setup.membrane import check_orientation

        top = md.Topology()
        residue = top.add_residue("ALA", top.add_chain(), resSeq=1)
        top.add_atom("CA", md.element.carbon, residue)
        assert check_orientation(top.to_openmm(), np.zeros((1, 3))) is None


class TestTheForceFieldNeedsLipidParameters:
    """Without them the run fails at system creation with a message about a
    residue template for POPC, which names the symptom rather than the missing
    file."""

    def test_a_bundle_that_already_has_them_is_left_alone(self) -> None:
        """The first version decided by looking for "lipid" in the filename.
        `amber14-all.xml` is a bundle carrying lipid17 inside it, so adding
        the file again raised about a duplicate template for a residue nobody
        mentioned. The question is what the force field can build, not what
        its files are called."""
        pytest.importorskip("openmm", reason="requires the [md] extra")

        from fastmdxplora.setup.membrane import membrane_forcefield_files

        given = ["amber14-all.xml", "amber14/tip3p.xml"]
        assert membrane_forcefield_files(given) == given

    def test_one_that_lacks_them_gets_them(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        from fastmdxplora.setup.membrane import membrane_forcefield_files

        found = membrane_forcefield_files(
            ["amber14/protein.ff14SB.xml", "amber14/tip3p.xml"])
        assert any("lipid" in f for f in found)

    def test_a_force_field_with_no_lipids_available_says_so(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        from fastmdxplora.setup.membrane import membrane_forcefield_files

        with pytest.raises(ValueError, match="lipid parameters"):
            membrane_forcefield_files(["some-other-field.xml"])


class TestTheBarostatKnowsItIsAMembrane:
    """An ordinary barostat scales x, y and z together, which squeezes a
    bilayer that should be free to change thickness independently of its area.

    Area per lipid is the number membrane simulations are validated against,
    so getting this wrong invalidates the run without stopping it.
    """

    @staticmethod
    def _topology(names):
        import mdtraj as md

        top = md.Topology()
        chain = top.add_chain()
        for name in names:
            residue = top.add_residue(name, chain, resSeq=1)
            top.add_atom("C", md.element.carbon, residue)
        return top.to_openmm()

    def test_lipids_are_detected_from_the_topology(self) -> None:
        """Rather than taken as a setting: the barostat has to be right
        whether or not anybody remembered to say so."""
        pytest.importorskip("openmm", reason="requires the [md] extra")

        from fastmdxplora.simulation.runner import is_membrane_system

        assert not is_membrane_system(self._topology(["ALA", "HOH"]))
        assert is_membrane_system(self._topology(["ALA", "POP", "HOH"]))

    def test_a_membrane_gets_the_membrane_barostat(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        import openmm as omm
        import openmm.unit as unit

        from fastmdxplora.simulation.runner import _add_barostat

        system = omm.System()
        system.addParticle(1.0)
        _add_barostat({"openmm": omm, "unit": unit}, system,
                      temperature_K=300.0, pressure_bar=1.0, frequency=25,
                      membrane=True)
        assert isinstance(system.getForce(0), omm.MonteCarloMembraneBarostat)

    def test_and_water_gets_the_ordinary_one(self) -> None:
        pytest.importorskip("openmm", reason="requires the [md] extra")

        import openmm as omm
        import openmm.unit as unit

        from fastmdxplora.simulation.runner import _add_barostat

        system = omm.System()
        system.addParticle(1.0)
        _add_barostat({"openmm": omm, "unit": unit}, system,
                      temperature_K=300.0, pressure_bar=1.0, frequency=25)
        assert isinstance(system.getForce(0), omm.MonteCarloBarostat)
        assert not isinstance(system.getForce(0),
                              omm.MonteCarloMembraneBarostat)

    def test_the_membrane_plane_is_coupled_and_the_normal_is_free(self) -> None:
        """Which is what makes it the right barostat: the bilayer keeps its
        area while its thickness moves."""
        pytest.importorskip("openmm", reason="requires the [md] extra")

        import openmm as omm
        import openmm.unit as unit

        from fastmdxplora.simulation.runner import _add_barostat

        system = omm.System()
        system.addParticle(1.0)
        _add_barostat({"openmm": omm, "unit": unit}, system,
                      temperature_K=300.0, pressure_bar=1.0, frequency=25,
                      membrane=True)
        barostat = system.getForce(0)
        assert barostat.getXYMode() == omm.MonteCarloMembraneBarostat.XYIsotropic
        assert barostat.getZMode() == omm.MonteCarloMembraneBarostat.ZFree
        assert barostat.getDefaultSurfaceTension().value_in_unit(
            unit.bar * unit.nanometer) == 0.0, (
            "a bilayer at equilibrium has no surface tension, and imposing "
            "one is a decision to make deliberately"
        )


class TestTheSettingsReachEveryInterface:
    def test_they_are_declared(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        membrane = PHASE_SCHEMAS["setup"].get("membrane")
        assert membrane is not None
        assert membrane.choices and "POPC" in membrane.choices
        assert PHASE_SCHEMAS["setup"].get("membrane_orientation_checked") is not None

    def test_they_reach_the_command_line(self) -> None:
        from fastmdxplora.cli.main import _PHASE_SPEC

        table, _prefix = _PHASE_SPEC["setup"]
        offered = {dest for _flag, dest, _kw in table}
        assert {"membrane", "membrane_orientation_checked"} <= offered

    def test_the_pipeline_passes_them(self) -> None:
        import inspect

        from fastmdxplora.setup import pipeline

        source = inspect.getsource(pipeline)
        assert "membrane=params.get" in source
