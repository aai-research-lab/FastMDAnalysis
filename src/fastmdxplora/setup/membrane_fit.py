"""Find the bilayer a structure belongs in, rather than assuming one.

Orienting by principal axes asks which way the protein is longest. That is
the right question for a lone transmembrane bundle and the wrong one for
anything carrying a soluble domain: 6B73's longest axis runs through its G
protein, and rotating onto it put one receptor through the bilayer upside
down. Checking the result afterwards did not help either, because every test
of a supplied orientation compares populations across the whole molecule, so
a T4L fusion or a beta-barrel's polar lumen swamps the belt that matters.

The approach here is the one PPM takes, simplified: search for the slab that
a lipid bilayer would occupy, by maximising the hydrophobic surface inside it
and the polar surface outside. Where a real bilayer fits, the optimum is
deep, narrow, and near the thickness of a real membrane. Where none does, the
best slab is shallow and its position arbitrary -- which is the difference
between a membrane protein and a soluble one, measured rather than assumed.

Two properties make this tractable. Accessible surface area does not change
when a structure is rotated, so it is computed once. And for a given
direction, every candidate slab centre and thickness can be scored together
from one sorted cumulative sum rather than one pass each.

    Lomize et al., Positioning of proteins in membranes: a computational
    approach, Protein Sci 2006.
    Lomize et al., OPM database and PPM web server, Nucleic Acids Res 2012.
    Lomize et al., Spatial arrangement of proteins in planar and curved
    membranes by PPM 3.0, Protein Sci 2022.

This is not PPM. PPM optimises a transfer free energy against an anisotropic
solvent model with dielectric and hydrogen-bonding profiles, and reports an
energy in kcal/mol. What follows keeps the shape of the idea -- fit the slab,
score buried hydrophobic area -- with a single hydrophobic term, because the
published result is that hydrophobic interactions alone usually determine the
position. Where the answer matters, an oriented structure from OPM is still
the better input.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

#: Half-thicknesses to consider, in nanometres. A DOPC bilayer's hydrophobic
#: core is about 3 nm across, and OPM's fitted thicknesses for transmembrane
#: proteins run roughly 2 to 4 nm; outside that a "membrane" is not one.
MINIMUM_HALF_THICKNESS_NM = 1.0
MAXIMUM_HALF_THICKNESS_NM = 2.0

#: How finely to step the slab's centre and half-thickness, in nanometres.
SLAB_STEP_NM = 0.1

#: How many directions to try. Spread over a hemisphere by the Fibonacci
#: construction, because a slab and its mirror are the same slab. About a
#: degree and a half between neighbours at this count, which is finer than
#: the tilt a bilayer imposes.
DIRECTION_COUNT = 1500


@dataclass(frozen=True)
class SlabFit:
    """Where a bilayer would sit, and how well it fits."""

    normal: Any
    """Unit vector along the membrane normal, in the input's frame."""

    centre_nm: float
    """Where the slab's midplane cuts that axis."""

    half_thickness_nm: float
    """Half the hydrophobic thickness."""

    buried_hydrophobic_nm2: float
    """Hydrophobic area the slab covers."""

    buried_polar_nm2: float
    """Polar area it covers, which a real bilayer keeps small."""

    score_nm2: float
    """Hydrophobic minus polar area inside the slab. The quantity maximised."""

    @property
    def thickness_nm(self) -> float:
        return 2.0 * self.half_thickness_nm

    @property
    def hydrophobic_fraction(self) -> float:
        """Of the area the slab covers, how much is hydrophobic.

        A bilayer-spanning protein presents a belt that is mostly hydrophobic;
        a soluble protein has no such band, so its best slab still covers a
        surface of ordinary composition.
        """
        total = self.buried_hydrophobic_nm2 + self.buried_polar_nm2
        if total <= 0:
            return 0.0
        return self.buried_hydrophobic_nm2 / total


def hemisphere_directions(count: int = DIRECTION_COUNT) -> Any:
    """Unit vectors spread evenly over a hemisphere.

    A hemisphere rather than a sphere: a slab is unchanged by flipping its
    normal, so searching both halves would do the same work twice.
    """
    import numpy as np

    index = np.arange(count, dtype=float) + 0.5
    # z from 0 to 1 gives the upper hemisphere; the golden angle spreads the
    # azimuth so no two directions bunch together.
    z = index / count
    radius = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    golden = np.pi * (1.0 + 5.0 ** 0.5)
    angle = golden * index
    return np.column_stack(
        [np.cos(angle) * radius, np.sin(angle) * radius, z])


def _best_slab_along(
    projections: Any,
    hydrophobic: Any,
    polar: Any,
    *,
    step: float,
    minimum_half: float,
    maximum_half: float,
) -> tuple[float, float, float, float]:
    """The best slab along one axis, from cumulative sums.

    Sorting once and accumulating means every centre and thickness is scored
    by two lookups rather than a pass over the atoms, which is what makes a
    search over a thousand directions affordable.
    """
    import numpy as np

    order = np.argsort(projections)
    ordered = projections[order]
    hydrophobic_sum = np.concatenate([[0.0], np.cumsum(hydrophobic[order])])
    polar_sum = np.concatenate([[0.0], np.cumsum(polar[order])])

    low, high = float(ordered[0]), float(ordered[-1])
    if high - low < 2.0 * minimum_half:
        return 0.0, 0.0, 0.0, 0.0

    centres = np.arange(low + minimum_half, high - minimum_half + step, step)
    halves = np.arange(minimum_half, maximum_half + step, step)
    if centres.size == 0 or halves.size == 0:
        return 0.0, 0.0, 0.0, 0.0

    # Edges of every (centre, half-thickness) pair at once.
    lower = centres[:, None] - halves[None, :]
    upper = centres[:, None] + halves[None, :]
    start = np.searchsorted(ordered, lower.ravel(), side="left")
    stop = np.searchsorted(ordered, upper.ravel(), side="right")

    inside_hydrophobic = hydrophobic_sum[stop] - hydrophobic_sum[start]
    inside_polar = polar_sum[stop] - polar_sum[start]
    scores = inside_hydrophobic - inside_polar

    best = int(np.argmax(scores))
    return (
        float(scores[best]),
        float(centres[best // halves.size]),
        float(halves[best % halves.size]),
        float(inside_hydrophobic[best]),
    ), float(inside_polar[best])


def fit_membrane_slab(
    coordinates: Any,
    hydrophobic_area: Any,
    polar_area: Any,
    *,
    directions: Any = None,
    step: float = SLAB_STEP_NM,
    minimum_half: float = MINIMUM_HALF_THICKNESS_NM,
    maximum_half: float = MAXIMUM_HALF_THICKNESS_NM,
) -> SlabFit | None:
    """Fit the bilayer slab that buries the most hydrophobic surface.

    ``coordinates`` are atom positions in nanometres. The two area arrays
    carry each atom's accessible surface split by character -- an atom
    contributes to one or the other, not both -- so the search never
    recomputes a surface that rotation cannot change.
    """
    import numpy as np

    points = np.asarray(coordinates, dtype=float)
    hydrophobic = np.asarray(hydrophobic_area, dtype=float)
    polar = np.asarray(polar_area, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        return None
    if len(points) != len(hydrophobic) or len(points) != len(polar):
        return None
    if len(points) < 20 or (hydrophobic.sum() + polar.sum()) <= 0:
        return None

    if directions is None:
        directions = hemisphere_directions()
    directions = np.asarray(directions, dtype=float)

    centred = points - points.mean(axis=0)

    best_score = -np.inf
    best: SlabFit | None = None
    for normal in directions:
        projections = centred @ normal
        (score, centre, half, buried_hydrophobic), buried_polar = _best_slab_along(
            projections, hydrophobic, polar,
            step=step, minimum_half=minimum_half, maximum_half=maximum_half,
        )
        if score > best_score:
            best_score = score
            best = SlabFit(
                normal=normal.copy(),
                centre_nm=centre,
                half_thickness_nm=half,
                buried_hydrophobic_nm2=buried_hydrophobic,
                buried_polar_nm2=buried_polar,
                score_nm2=score,
            )
    return best


def rotation_onto_z(normal: Any) -> Any:
    """A rotation taking ``normal`` to the z axis.

    Built from a single rotation about the axis perpendicular to both, so
    nothing else about the structure's orientation is disturbed -- the fit
    decides which way is up and no more.
    """
    import numpy as np

    source = np.asarray(normal, dtype=float)
    source = source / np.linalg.norm(source)
    target = np.array([0.0, 0.0, 1.0])

    axis = np.cross(source, target)
    sine = float(np.linalg.norm(axis))
    cosine = float(np.dot(source, target))
    if sine < 1e-12:
        # Already along z, or exactly opposed.
        return np.eye(3) if cosine > 0 else np.diag([1.0, -1.0, -1.0])

    axis = axis / sine
    cross = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return np.eye(3) + cross * sine + cross @ cross * (1.0 - cosine)


#: Elements whose exposed surface a bilayer is happy to touch. Hydrogen is
#: assigned by what it is bonded to: the surface of a methyl is apolar, the
#: surface of a hydroxyl's hydrogen is not.
_APOLAR_ELEMENTS = frozenset({"C", "S", "SE"})


def areas_by_character(topology: Any, positions: Any) -> tuple[Any, Any, Any]:
    """Per-atom accessible surface, split into apolar and polar.

    Computed once and reused for every candidate orientation, because
    rotating a structure does not change how much of it is exposed. Returns
    the coordinates in nanometres alongside the two area arrays, so callers
    do not have to unwrap units twice.
    """
    import mdtraj as md
    import numpy as np

    try:
        from openmm import unit as openmm_unit

        if openmm_unit.is_quantity(positions):
            positions = positions.value_in_unit(openmm_unit.nanometer)
    except ImportError:  # pragma: no cover - only without OpenMM
        pass

    coordinates = np.asarray(
        [[float(p[0]), float(p[1]), float(p[2])] for p in positions], dtype=float)

    mdtop = md.Topology.from_openmm(topology) if not isinstance(
        topology, md.Topology) else topology
    trajectory = md.Trajectory(coordinates[None, :, :], mdtop)
    per_atom = md.shrake_rupley(trajectory, mode="atom")[0]

    heavy_partner = {}
    for first, second in mdtop.bonds:
        for atom, other in ((first, second), (second, first)):
            symbol = getattr(atom.element, "symbol", "") or ""
            if symbol.upper() == "H":
                heavy_partner[atom.index] = other

    apolar = np.zeros(len(per_atom))
    polar = np.zeros(len(per_atom))
    for atom in mdtop.atoms:
        symbol = (getattr(atom.element, "symbol", "") or "").upper()
        if symbol == "H":
            partner = heavy_partner.get(atom.index)
            symbol = (getattr(getattr(partner, "element", None), "symbol", "")
                      or "").upper()
        if symbol in _APOLAR_ELEMENTS:
            apolar[atom.index] = per_atom[atom.index]
        else:
            polar[atom.index] = per_atom[atom.index]

    return coordinates, apolar, polar


def fit_membrane(topology: Any, positions: Any, **kwargs: Any) -> SlabFit | None:
    """Fit a bilayer to a structure, from its topology and coordinates."""
    coordinates, apolar, polar = areas_by_character(topology, positions)
    return fit_membrane_slab(coordinates, apolar, polar, **kwargs)
