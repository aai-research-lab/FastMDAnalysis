"""What holds a ligand in place, one interaction type at a time.

Counting contacts says how much of the protein the ligand touches. It does not
say what is holding it: a salt bridge that a charge change would destroy, or
hydrophobic packing that tolerates one. That difference is what a medicinal
chemist is asking, and it is the difference between counting and describing.

Each rule here is implemented against a published criterion, named in its own
docstring, and the thresholds are settings rather than constants -- because the
published values disagree, and a reader should be able to see which was used.
Where a widely-used tool departs from the literature, the departure is recorded
rather than quietly adopted.

Geometry only. Which atoms can take part is chemistry, and that is settled
before this runs: see ``ligand_chemistry``, which also records how confidently.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "Contact",
    "hydrogen_bonds",
    "hydrophobic_contacts",
    "donors_and_acceptors",
    "hydrophobic_atoms",
]


#: Elements that donate and accept hydrogen bonds. Carbon is excluded: C-H
#: donors exist and are real, but they are weak enough that including them
#: swamps the count, which is the same reason every tool surveyed leaves them
#: out.
_POLAR = frozenset({"N", "O", "S"})


@dataclass(frozen=True)
class Contact:
    """One interaction, in one frame."""

    kind: str
    frame: int
    #: Atom indices into the trajectory, ligand side first.
    ligand_atom: int
    protein_atom: int
    distance_nm: float
    #: Set where the rule has one; None where it does not.
    angle_deg: float | None = None


def donors_and_acceptors(
    topology: Any, atom_indices: Any
) -> tuple[list[tuple[int, int]], list[int]]:
    """Which atoms can donate a hydrogen bond, and which can accept one.

    A donor is a nitrogen, oxygen or sulphur with a hydrogen bonded to it; an
    acceptor is any of those three. This is the criterion Baker and Hubbard
    used (Prog Biophys Mol Biol 44:97, 1984) and it is what every tool
    surveyed uses, because with explicit hydrogens present it needs no
    perception: the hydrogen is either bonded to the nitrogen or it is not.

    Returns donors as ``(heavy, hydrogen)`` pairs, because the angle is
    measured at the hydrogen and a nitrogen with two hydrogens can donate
    twice.
    """
    wanted = set(int(i) for i in atom_indices)
    atoms = list(topology.atoms)

    attached: dict[int, list[int]] = {}
    for bond in topology.bonds:
        first, second = bond[0], bond[1]
        for heavy, light in ((first, second), (second, first)):
            if light.element is not None and light.element.symbol == "H":
                attached.setdefault(heavy.index, []).append(light.index)

    # A hydrogen bonded to nothing is the signal that connectivity is missing.
    # A hydrogen bonded to a carbon is not: that is an ordinary non-polar
    # hydrogen, and a selection made only of those legitimately cannot donate.
    bonded_atoms = {a.index for bond in topology.bonds for a in (bond[0], bond[1])}

    donors: list[tuple[int, int]] = []
    acceptors: list[int] = []
    orphan_hydrogens = 0
    for atom in atoms:
        if atom.index not in wanted:
            continue
        element = atom.element.symbol if atom.element is not None else ""
        if element == "H":
            if atom.index not in bonded_atoms:
                orphan_hydrogens += 1
            continue
        if element not in _POLAR:
            continue
        acceptors.append(atom.index)
        for hydrogen in attached.get(atom.index, ()):
            donors.append((atom.index, hydrogen))

    # A selection whose hydrogens are bonded to nothing cannot be seen to
    # donate, only to accept, and the count would report one direction as
    # though it were both. This is not hypothetical: renaming a residue in a
    # PDB is enough to lose its standard bonds, and a ligand deposited without
    # CONECT records arrives the same way. Saying so beats returning zero.
    if orphan_hydrogens:
        raise ValueError(
            f"{orphan_hydrogens} hydrogen(s) in this selection are bonded to "
            "nothing, so it cannot be seen to donate a hydrogen bond -- only "
            "to accept one. "
            "Counting under those conditions would report one direction as "
            "though it were both. Give a topology that carries the "
            "connectivity: a PDB with CONECT records, or the topology the "
            "setup phase writes."
        )
    return donors, acceptors


def hydrophobic_atoms(topology: Any, atom_indices: Any) -> list[int]:
    """Carbons whose only neighbours are carbon or hydrogen.

    The criterion PLIP states and ProLIF's SMARTS encodes. A carbon bonded to
    nitrogen, oxygen or fluorine is polarised enough that treating it as
    hydrophobic would count a polar contact twice -- once here and once as a
    hydrogen bond.
    """
    wanted = set(int(i) for i in atom_indices)
    neighbours: dict[int, set[str]] = {}
    for bond in topology.bonds:
        first, second = bond[0], bond[1]
        for atom, other in ((first, second), (second, first)):
            symbol = other.element.symbol if other.element is not None else ""
            neighbours.setdefault(atom.index, set()).add(symbol)

    found = []
    for atom in topology.atoms:
        if atom.index not in wanted:
            continue
        if atom.element is None or atom.element.symbol != "C":
            continue
        around = neighbours.get(atom.index, set())
        if around and around <= {"C", "H"}:
            found.append(atom.index)
    return found


def _angles(traj: Any, triples: Any, periodic: bool) -> np.ndarray:
    import mdtraj as md

    if len(triples) == 0:
        return np.zeros((traj.n_frames, 0))
    return np.rad2deg(md.compute_angles(traj, triples, periodic=periodic))


def _distances(traj: Any, pairs: Any, periodic: bool) -> np.ndarray:
    import mdtraj as md

    if len(pairs) == 0:
        return np.zeros((traj.n_frames, 0))
    return md.compute_distances(traj, pairs, periodic=periodic)


def hydrogen_bonds(
    traj: Any,
    ligand_indices: Any,
    protein_indices: Any,
    *,
    distance_nm: float = 0.35,
    angle_deg: float = 120.0,
    periodic: bool = True,
) -> list[Contact]:
    """Hydrogen bonds between ligand and protein, in every frame.

    The criterion is the literature standard: donor to acceptor within 3.5 A,
    and the donor-hydrogen-acceptor angle above 120 degrees (Baker & Hubbard
    1984; McDonald & Thornton, J Mol Biol 238:777, 1994).

    PLIP uses 4.1 A and 100 degrees, which it says is deliberate -- refined
    against low-resolution crystal structures where hydrogen positions are not
    known. That reasoning does not apply here: these frames have hydrogens
    placed by the force field, so the angle is measured rather than inferred,
    and the stricter criterion is the one the measurement supports. PLIP's
    values remain reachable by setting them.

    Both directions are found: the ligand donating to the protein and the
    protein donating to the ligand are different interactions, and a ligand
    that can only accept is a fact about the ligand worth seeing.
    """
    ligand_donors, ligand_acceptors = donors_and_acceptors(
        traj.topology, ligand_indices)
    protein_donors, protein_acceptors = donors_and_acceptors(
        traj.topology, protein_indices)

    triples: list[tuple[int, int, int]] = []
    pairs: list[tuple[int, int]] = []
    sides: list[tuple[int, int]] = []      # (ligand atom, protein atom)
    for heavy, hydrogen in ligand_donors:
        for acceptor in protein_acceptors:
            triples.append((heavy, hydrogen, acceptor))
            pairs.append((heavy, acceptor))
            sides.append((heavy, acceptor))
    for heavy, hydrogen in protein_donors:
        for acceptor in ligand_acceptors:
            triples.append((heavy, hydrogen, acceptor))
            pairs.append((heavy, acceptor))
            sides.append((acceptor, heavy))

    if not triples:
        return []

    separations = _distances(traj, np.array(pairs), periodic)
    openings = _angles(traj, np.array(triples), periodic)
    close_enough = separations < distance_nm
    straight_enough = openings > angle_deg

    found: list[Contact] = []
    for frame, column in zip(*np.where(close_enough & straight_enough)):
        ligand_atom, protein_atom = sides[column]
        found.append(Contact(
            kind="hydrogen_bond",
            frame=int(frame),
            ligand_atom=int(ligand_atom),
            protein_atom=int(protein_atom),
            distance_nm=float(separations[frame, column]),
            angle_deg=float(openings[frame, column]),
        ))
    return found


def hydrophobic_contacts(
    traj: Any,
    ligand_indices: Any,
    protein_indices: Any,
    *,
    distance_nm: float = 0.40,
    periodic: bool = True,
) -> list[Contact]:
    """Hydrophobic contacts between ligand and protein, in every frame.

    A pair of hydrophobic carbons within 4.0 A, which is PLIP's threshold;
    ProLIF uses 4.5 A. There is no angle: hydrophobic association is entropic
    rather than directional, so there is no geometry to require, which is why
    every tool uses a distance alone and a generous one.

    The count is large by nature -- PLIP notes it can exceed every other type
    combined -- so it is reported per atom pair and left to the caller to
    reduce. Reducing here would bake in one view of which contact represents a
    residue, and the reduction PLIP applies is one choice among several.
    """
    ligand_carbons = hydrophobic_atoms(traj.topology, ligand_indices)
    protein_carbons = hydrophobic_atoms(traj.topology, protein_indices)
    if not ligand_carbons or not protein_carbons:
        return []

    pairs = [(a, b) for a in ligand_carbons for b in protein_carbons]
    separations = _distances(traj, np.array(pairs), periodic)

    found: list[Contact] = []
    for frame, column in zip(*np.where(separations < distance_nm)):
        ligand_atom, protein_atom = pairs[column]
        found.append(Contact(
            kind="hydrophobic",
            frame=int(frame),
            ligand_atom=int(ligand_atom),
            protein_atom=int(protein_atom),
            distance_nm=float(separations[frame, column]),
        ))
    return found
