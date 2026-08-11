"""A metal in a protein site, and a force field that will not keep it there.

Standard force fields treat a metal ion as a point charge with Lennard-Jones
terms and nothing else. A carboxylate cage is deep enough to hold a divalent
ion that way. A histidine-ligated zinc is not, so it leaves -- and nothing
else in a run reports it: the fold is intact, the RMSD is flat, the radius of
gyration does not move, and the site the study was about has rearranged.

Measured on thermolysin over 20 ns. All four calciums held their carboxylate
sites to within a tenth of an Angstrom from first frame to last. The catalytic
zinc lost His142 in the first production frame and never recovered it, ending
on a different set of residues from the one in the crystal structure.

Not a refusal: non-bonded metals are the right choice for many questions,
including anything away from the site. It is a refusal to be silent.
"""

from __future__ import annotations

import numpy as np
import pytest

app = pytest.importorskip("openmm.app")
from openmm import unit  # noqa: E402

from fastmdxplora.setup.prepare import (  # noqa: E402
    STRUCTURAL_METALS, _warn_metals_are_not_held)

E = app.element


def _system(entries):
    """entries: (resname, atomname, element, x, y, z) with x,y,z in nm."""
    topology = app.Topology()
    chain = topology.addChain()
    coordinates = []
    for resname, atomname, element, x, y, z in entries:
        residue = topology.addResidue(resname, chain)
        topology.addAtom(atomname, element, residue)
        coordinates.append([x, y, z])
    return topology, np.array(coordinates) * unit.nanometer


HIS = ("HIS", "NE2", E.nitrogen, 0.0, 0.0, 0.0)
ASP = ("ASP", "OD1", E.oxygen, 1.0, 0.0, 0.0)


class TestItFiresOnAMetalInASite:
    def test_a_zinc_beside_a_histidine(self) -> None:
        said = _warn_metals_are_not_held(
            *_system([HIS, ("ZN", "ZN", E.zinc, 0.2, 0.0, 0.0)]))
        assert said is not None and "ZN" in said

    def test_several_are_named_with_their_counts(self) -> None:
        said = _warn_metals_are_not_held(*_system([
            HIS, ("ZN", "ZN", E.zinc, 0.2, 0.0, 0.0),
            ASP, ("CA", "CA", E.calcium, 1.24, 0.0, 0.0)]))
        assert "ZN" in said and "CA" in said

    def test_it_says_what_to_do(self) -> None:
        """Naming the problem without the remedies would leave somebody
        knowing their run is suspect and not what to do about it."""
        said = _warn_metals_are_not_held(
            *_system([HIS, ("ZN", "ZN", E.zinc, 0.2, 0.0, 0.0)]))
        assert "restrain" in said
        assert "check the coordination at the end" in said


class TestItStaysQuietWhereThereIsNothingToSay:
    def test_salt_is_not_a_structural_metal(self) -> None:
        """Sodium and chloride are meant to move. A warning about them would
        fire on every run and teach people to ignore it."""
        said = _warn_metals_are_not_held(*_system([
            HIS, ("NA", "NA", E.sodium, 2.0, 0.0, 0.0),
            ("CL", "CL", E.chlorine, 3.0, 0.0, 0.0)]))
        assert said is None

    def test_a_metal_in_the_solvent_is_not_in_a_site(self) -> None:
        said = _warn_metals_are_not_held(
            *_system([HIS, ("MG", "MG", E.magnesium, 1.0, 0.0, 0.0)]))
        assert said is None

    def test_a_protein_with_no_metals_says_nothing(self) -> None:
        assert _warn_metals_are_not_held(*_system([HIS, ASP])) is None

    def test_the_usual_structural_metals_are_listed(self) -> None:
        assert {"ZN", "MG", "CA", "MN", "FE", "CU"} <= STRUCTURAL_METALS

    def test_the_salt_ions_are_not(self) -> None:
        assert not ({"NA", "K", "CL", "BR", "F"} & STRUCTURAL_METALS)


class TestItNeverStopsARunItCannotJudge:
    def test_positions_carrying_units_are_read(self) -> None:
        """OpenMM returns a Quantity of Vec3, so the unit comes off the whole
        array rather than each component. Taking `float()` of a Vec3 raises,
        and a broad except turned that into permanent silence -- the warning
        was there and said nothing, which is the fault it exists to prevent.
        """
        said = _warn_metals_are_not_held(
            *_system([HIS, ("ZN", "ZN", E.zinc, 0.2, 0.0, 0.0)]))
        assert said is not None, "unit-carrying positions were not read"

    def test_plain_numbers_are_read_too(self) -> None:
        topology, positions = _system(
            [HIS, ("ZN", "ZN", E.zinc, 0.2, 0.0, 0.0)])
        bare = np.asarray(positions.value_in_unit(unit.nanometer))
        assert _warn_metals_are_not_held(topology, bare) is not None

    def test_a_topology_that_will_not_walk_is_not_an_error(self) -> None:
        class _Broken:
            def residues(self):
                raise RuntimeError("no")

        assert _warn_metals_are_not_held(_Broken(), np.zeros((1, 3))) is None
