"""Putting a protein in a lipid bilayer, and the things that go wrong quietly.

A membrane protein simulated in water is not the protein: the hydrophobic belt
that sits in the bilayer is exposed to solvent, the helices splay, and the
trajectory describes a molecule that does not exist. So a membrane system has
to be built as one.

OpenMM builds the bilayer itself, which removes the usual dependency on an
external packing tool. What it does not do is check the two things that make
the result meaningful, and both fail silently:

**The protein has to be oriented before it is embedded.** ``addMembrane``
places the bilayer in the xy plane with its normal along z, and assumes the
protein is already lying that way. A structure taken straight from the PDB
usually is not: crystallographic frames have no relation to a membrane normal.
Embed an unoriented protein and the lipids pack around it sideways, the run
proceeds, and every number that follows is about a structure nobody would
recognise. The OPM database publishes structures already oriented, and this
module checks the orientation rather than trusting it.

**The pressure coupling has to be anisotropic.** An ordinary barostat scales
x, y and z together, which squeezes a bilayer that should be free to change
thickness independently of its area. The result is a wrong area per lipid,
which is the number membrane simulations are validated against. OpenMM has
``MonteCarloMembraneBarostat`` for this, and using the ordinary one is a
mistake that runs to completion.
"""

from __future__ import annotations

from typing import Any

__all__ = [
    "LIPIDS",
    "membrane_forcefield_files",
    "check_orientation",
    "OrientationWarning",
]


#: Lipids OpenMM can build a bilayer from, with what each is usually for.
#: Not a list this software maintains -- it is what ``Modeller.addMembrane``
#: accepts, and offering one it does not would fail at the point of use.
LIPIDS: dict[str, str] = {
    "POPC": "the default for most membrane work; a fluid phosphatidylcholine",
    "POPE": "phosphatidylethanolamine, common in bacterial inner membranes",
    "DLPC": "a short-tailed phosphatidylcholine, thinner bilayer",
    "DLPE": "the phosphatidylethanolamine equivalent",
    "DMPC": "dimyristoyl; often used where a thinner membrane is wanted",
    "DOPC": "dioleoyl; two unsaturated tails, more fluid",
    "DPPC": "dipalmitoyl; gel phase at room temperature",
}

#: Force field files carrying lipid parameters, by the protein force field
#: they go with. A membrane system built against a force field with no lipid
#: parameters fails when the system is created, with a message about a
#: residue template rather than about the missing file.
_LIPID_PARAMETERS = {
    "amber14": "amber14/lipid17.xml",
    "charmm36": "charmm36/waters.xml",
}


class OrientationWarning(UserWarning):
    """The protein does not look like it is oriented for a membrane."""


def membrane_forcefield_files(existing: list[str], lipid: str = "POPC") -> list[str]:
    """Add lipid parameters to a force field that lacks them.

    The test is whether the force field can already build the lipid, not what
    its files are called: ``amber14-all.xml`` is a bundle that carries lipid17
    inside it, and adding the file again raises about a duplicate template for
    a residue nobody mentioned. Reading the name would have missed that, and
    did.

    Without the parameters the run fails at system creation with a message
    about a residue template for POPC, which names the symptom rather than the
    missing file.
    """
    files = list(existing)

    try:
        from openmm.app import ForceField

        if any(name.startswith(lipid) for name in ForceField(*files)._templates):
            return files
    except Exception:  # noqa: BLE001 - if it will not load, say so below
        pass

    for family, parameters in _LIPID_PARAMETERS.items():
        if any(family in f.lower() for f in files):
            candidate = files + [parameters]
            try:
                from openmm.app import ForceField

                ForceField(*candidate)
            except Exception as exc:  # noqa: BLE001
                raise ValueError(
                    f"Adding {parameters} for the {lipid} bilayer did not "
                    f"work: {exc}"
                ) from exc
            return candidate

    raise ValueError(
        f"A {lipid} bilayer needs lipid parameters, and the force field files "
        f"given carry none: {', '.join(files)}. The amber14 and charmm36 "
        "families both have them; if you are using another, add its lipid "
        "parameter file to `force_field` yourself."
    )


def check_orientation(topology: Any, positions: Any) -> str | None:
    """Whether the protein looks oriented with the membrane normal along z.

    Returns ``None`` where it does, and a description of the problem where it
    does not.

    The test is the shape of the protein: a membrane protein is longer along
    the normal than across it, because it spans a bilayer. A structure lying
    on its side has its longest principal axis in the membrane plane instead.
    That is a weak test and it is meant to be -- it catches the common case of
    a PDB entry used as it came, which is the one that produces a silently
    meaningless run. It cannot tell a correctly oriented protein from a
    subtly misoriented one, and it does not claim to.
    """
    import numpy as np

    import mdtraj as md

    mdtop = md.Topology.from_openmm(topology)
    protein = mdtop.select("protein")
    if len(protein) < 20:
        return None  # too small to say anything about

    # Positions may arrive as an OpenMM Quantity of Vec3 or as a plain array,
    # depending on where they came from. Stripped to nanometres either way:
    # the shape of a structure does not depend on how it was handed over.
    try:
        from openmm import unit as openmm_unit

        if openmm_unit.is_quantity(positions):
            positions = positions.value_in_unit(openmm_unit.nanometer)
    except ImportError:  # pragma: no cover - only reached without OpenMM
        pass

    coordinates = np.asarray(
        [[float(positions[int(i)][axis]) for axis in range(3)]
         for i in protein], dtype=float)
    centred = coordinates - coordinates.mean(axis=0)
    # Principal axes: the direction the structure is most extended in.
    _values, vectors = np.linalg.eigh(np.cov(centred.T))
    longest = vectors[:, -1]

    # How much of the longest axis lies along z. A membrane protein spanning
    # a bilayer is extended along the normal.
    along_z = abs(float(longest[2]))
    if along_z >= 0.7:
        return None

    return (
        f"The protein's longest axis lies {along_z:.0%} along z, where a "
        "membrane-spanning protein would be most extended along the membrane "
        "normal. `addMembrane` places the bilayer in the xy plane and assumes "
        "the protein is already oriented for it, so a structure taken "
        "straight from the PDB is usually in the wrong frame — "
        "crystallographic axes have no relation to a membrane normal. "
        "Embedding it anyway produces a run that completes and describes a "
        "structure nobody would recognise.\n\n"
        "The OPM database (https://opm.phar.umich.edu) publishes structures "
        "already oriented in a membrane; its coordinates can be used "
        "directly. Pass `membrane_orientation_checked: true` to proceed with "
        "the structure as it is."
    )
