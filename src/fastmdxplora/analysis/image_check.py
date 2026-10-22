"""How close the solute comes to its own periodic images, frame by frame.

A periodic box holds one copy of the solute and is surrounded by its images.
Where any part of the solute comes within the nonbonded cutoff of any part of
an image, the solute is interacting with its own copy, and every quantity
measured of those frames describes the cell as much as the molecule. Setup
sizes the box from the starting structure, so a chain that unfolds, or a
complex that turns, can reach its image though the box was sound at the
start. A 10-residue chain in a chignolin study came to 0.261 nm of its own
copy while the check of its two ends alone passed it.

So the whole solute is checked, as ``gmx mindist -pi`` checks it: in each
frame, the closest approach of any solute heavy atom to any heavy atom of the
solute's periodic images, kept as a series, with the frames under the cutoff
counted, and those under 0.5 nm.

**Made whole first.** Each frame's solute is made whole along its bonds
(between heavy atoms; where the topology carries no bonds, along the chain's
atom order, consecutive atoms of a chain being neighbours): every atom is
placed at the minimum image of the atom it is bonded to, which is exact for a
bond shorter than half the cell's narrowest width. Where the solute is
several molecules (chains, a ligand), each keeps the place the loader gave it
beside the others (:func:`fastmdxplora.analysis.loading._made_whole`).

**The images searched.** Every lattice vector ``L = i a + j b + k c`` with
``i, j, k`` from -2 to 2 but not all zero (two cells either way, in any cell
shape, triclinic included), and for each the nearest pair between the solute
and its image ``L`` away, found with a k-d tree. ``L`` and ``-L`` give the
same distance, so half are searched; and an image whose shift is longer than
the closest approach found so far plus the solute's diameter cannot be
closer, so the search stops there. Frames where it did not stop before the
two cells ran out are counted: an image further away than two cells could
not be ruled out in them.

The check's finding (``findings.solute_image``) is the one finding every
other analysis of the same frames reads (the orchestrator hands it to each):
where the solute reaches its image, each mean is given without its error
bar and with this as the reason, through the same path the end-to-end
distance used for its two ends (:meth:`Analysis._withhold_the_error`), and
the end-to-end distance reads it rather than checking the ends alone.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

from fastmdxplora.analysis.base import Analysis, narrowest_width
from fastmdxplora.analysis.orchestrator import register_analysis
from fastmdxplora.analysis.plotting import colour
from fastmdxplora.refusals import StudyError

#: The nonbonded cutoff a solute's image is checked against: the cutoff setup
#: uses unless a force field or the person sets another.
CUTOFF_NM = 1.0

#: Closer than this, the solute is in contact with its own copy, not merely
#: within the cutoff of it.
CLOSE_NM = 0.5

#: How many cells either way the images are searched.
SEARCHED_CELLS = 2

#: Frames made whole and searched at a time, so memory stays bounded.
FRAMES_AT_A_TIME = 1000


def _heavy(topology, atoms: np.ndarray) -> np.ndarray:
    """The heavy atoms among ``atoms``: those with an element that is not
    hydrogen (a virtual site has none)."""
    kept = []
    for index in atoms:
        element = topology.atom(int(index)).element
        if element is not None and element.symbol not in ("H", "D") and element.mass > 0:
            kept.append(int(index))
    return np.asarray(kept, dtype=np.int64)


def walk(topology, atoms: np.ndarray) -> tuple[list[tuple[int, int]], str]:
    """The steps that make ``atoms`` whole, as pairs of positions in
    ``atoms`` (the one placed, the one placed beside), in the order they are
    taken; and how they were found."""
    where = {int(a): k for k, a in enumerate(atoms)}
    neighbours: list[list[int]] = [[] for _ in atoms]
    if topology.n_bonds:
        how = "along its bonds"
        for one, other in topology.bonds:
            i, j = where.get(one.index), where.get(other.index)
            if i is not None and j is not None:
                neighbours[i].append(j)
                neighbours[j].append(i)
    else:
        how = ("along the chain's atom order, the topology carrying no bonds")
        for k in range(1, len(atoms)):
            if (topology.atom(int(atoms[k])).residue.chain.index
                    == topology.atom(int(atoms[k - 1])).residue.chain.index):
                neighbours[k].append(k - 1)
                neighbours[k - 1].append(k)
    steps: list[tuple[int, int]] = []
    seen = np.zeros(len(atoms), dtype=bool)
    for root in range(len(atoms)):
        if seen[root]:
            continue
        seen[root] = True
        queue = [root]
        while queue:
            parent = queue.pop(0)
            for child in neighbours[parent]:
                if not seen[child]:
                    seen[child] = True
                    steps.append((parent, child))
                    queue.append(child)
    return steps, how


def made_whole(xyz: np.ndarray, boxes: np.ndarray,
               steps: list[tuple[int, int]]) -> np.ndarray:
    """``xyz`` (frames, atoms, 3) with each atom moved to the minimum image of
    the one it steps from, in place; ``boxes`` (frames, 3, 3), rows the cell
    vectors."""
    inverse = np.linalg.inv(boxes)
    for parent, child in steps:
        step = xyz[:, child] - xyz[:, parent]
        cells = np.rint(np.einsum("fi,fij->fj", step, inverse))
        step -= np.einsum("fj,fji->fi", cells, boxes)
        xyz[:, child] = xyz[:, parent] + step
    return xyz


def _shifts(cells: int = SEARCHED_CELLS) -> np.ndarray:
    """The lattice shifts searched, one of each pair L and -L."""
    span = np.arange(-cells, cells + 1)
    every = np.array(np.meshgrid(span, span, span, indexing="ij")).reshape(3, -1).T
    first = np.array([row[np.nonzero(row)[0][0]] if row.any() else 0 for row in every])
    return every[first > 0].astype(np.float64)


def nearest_image(points: np.ndarray, box: np.ndarray,
                  shifts: np.ndarray) -> tuple[float, int, int, bool]:
    """The closest approach of ``points`` to their images ``shifts @ box``
    away: the distance, the atom (a position in ``points``), the atom whose
    image it is, and whether the search stopped before the shifts ran out
    (an image further away could not then be closer)."""
    from scipy.spatial import cKDTree

    centre = points.mean(axis=0)
    radius = float(np.sqrt(((points - centre) ** 2).sum(axis=1).max()))
    lattice = shifts @ box
    lengths = np.sqrt((lattice ** 2).sum(axis=1))
    order = np.argsort(lengths, kind="stable")
    tree = cKDTree(points)
    best, near_atom, image_atom = np.inf, -1, -1
    for k in order:
        # |x_j + L - x_i| >= |L| - |x_j - c| - |x_i - c| >= |L| - 2r.
        if lengths[k] - 2.0 * radius >= best:
            return best, near_atom, image_atom, True
        image = points + lattice[k]
        # An image atom within `best` of any atom is within r + best of the centre.
        candidates = np.nonzero(((image - centre) ** 2).sum(axis=1)
                                < (radius + best) ** 2)[0]
        if candidates.size == 0:
            continue
        found, nearest = tree.query(image[candidates], k=1, distance_upper_bound=best)
        which = int(np.argmin(found))
        if found[which] < best:
            best = float(found[which])
            near_atom, image_atom = int(nearest[which]), int(candidates[which])
    return best, near_atom, image_atom, False


def _atom_said(topology, index: int) -> str:
    atom = topology.atom(index)
    residue = atom.residue
    chain = getattr(residue.chain, "chain_id", None)
    label = f"{residue.name}{residue.resSeq} {atom.name}"
    if topology.n_chains > 1 and chain:
        label = f"{label} (chain {chain})"
    return label


def solute_image_record(traj: md.Trajectory, atoms: np.ndarray, *,
                        cutoff_nm: float = CUTOFF_NM) -> tuple[np.ndarray, dict[str, Any]]:
    """Per frame, the closest approach of the heavy atoms among ``atoms`` to
    the heavy atoms of their periodic images, in nm, and the finding.

    The finding (``reaches_its_image``, the frames under ``cutoff_nm`` and
    under :data:`CLOSE_NM`, the closest approach with its frame and atoms,
    and ``said``, the sentence every reader is given) is what every analysis
    of these frames is qualified by.
    """
    if traj.unitcell_vectors is None:
        raise StudyError(
            "This trajectory carries no unit cell, so the solute has no "
            "periodic images to come near.",
            code="analysis.system.inapplicable", analysis="image_check",
            reason="no unit cell")
    heavy = _heavy(traj.topology, np.asarray(atoms, dtype=np.int64))
    if heavy.size == 0:
        raise StudyError(
            "The solute selection holds no heavy atoms, so there is nothing to "
            "check against its images.",
            code="analysis.selection.empty", expression=None, role="solute")
    steps, how = walk(traj.topology, heavy)
    shifts = _shifts()
    n = traj.n_frames
    distances = np.empty(n, dtype=np.float64)
    pairs = np.empty((n, 2), dtype=np.int64)
    open_ended = 0
    for low in range(0, n, FRAMES_AT_A_TIME):
        high = min(low + FRAMES_AT_A_TIME, n)
        boxes = np.asarray(traj.unitcell_vectors[low:high], dtype=np.float64)
        xyz = made_whole(np.asarray(traj.xyz[low:high][:, heavy], dtype=np.float64),
                         boxes, steps)
        for f in range(high - low):
            d, i, j, closed = nearest_image(xyz[f], boxes[f], shifts)
            distances[low + f], pairs[low + f] = d, (i, j)
            open_ended += not closed

    frame = int(np.argmin(distances))
    nearest = float(distances[frame])
    under = int(np.count_nonzero(distances < cutoff_nm))
    close = int(np.count_nonzero(distances < CLOSE_NM))
    times = np.asarray(traj.time, dtype=float) if traj.time is not None else None
    timed = times is not None and times.size == n and bool(np.all(np.isfinite(times)))
    widths = narrowest_width(traj)
    record: dict[str, Any] = {
        "cutoff_nm": float(cutoff_nm),
        "close_nm": CLOSE_NM,
        "frames": n,
        "nearest_nm": nearest,
        "nearest_frame": frame,
        "nearest_time_ns": float(times[frame]) / 1000.0 if timed else None,
        "nearest_atoms": [_atom_said(traj.topology, int(heavy[pairs[frame, 0]])),
                          _atom_said(traj.topology, int(heavy[pairs[frame, 1]]))],
        "frames_under_cutoff": under,
        "fraction_under_cutoff": under / n,
        "frames_under_close": close,
        "fraction_under_close": close / n,
        "heavy_atoms": int(heavy.size),
        "made_whole": how,
        "searched_cells": SEARCHED_CELLS,
        "frames_beyond_not_ruled_out": int(open_ended),
        "narrowest_width_nm": [float(np.min(widths)), float(np.max(widths))],
        "reaches_its_image": bool(under),
    }
    when = (f"frame {frame}, {record['nearest_time_ns']:.4g} ns" if timed
            else f"frame {frame}")
    atoms_said = (f"{record['nearest_atoms'][0]} and the image of "
                  f"{record['nearest_atoms'][1]}")
    if under:
        record["said"] = (
            f"The solute comes within {nearest:.3f} nm of its own periodic image "
            f"({when}: {atoms_said}), inside the {cutoff_nm:g} nm cutoff in "
            f"{under:,} of {n:,} frames ({100.0 * under / n:.2f}%) and under "
            f"{CLOSE_NM:g} nm in {close:,} ({100.0 * close / n:.2f}%). In those "
            "frames it is interacting with its own copy, and what is measured of "
            "it describes the cell as much as the molecule. A larger box, in a "
            "new study (a box cannot grow mid-run), is what would change this.")
    else:
        record["said"] = (
            f"The solute comes no closer than {nearest:.3f} nm to its own periodic "
            f"image ({when}: {atoms_said}), outside the {cutoff_nm:g} nm cutoff in "
            f"every one of {n:,} frames.")
    if open_ended:
        record["said"] += (
            f" In {open_ended:,} frames the solute is long enough against the cell "
            f"that an image further than {SEARCHED_CELLS} cells away could not be "
            "ruled out.")
    return distances, record


class SoluteImage(Analysis):
    """The solute's closest approach to its own periodic images, per frame.

    Parameters
    ----------
    cutoff_nm : float, default 1.0
        The nonbonded cutoff the closest approach is checked against, in nm:
        the one setup uses unless a force field or the person set another.
    selection : str, optional
        MDTraj atom selection for the solute; its heavy atoms are checked.
        The analysis scope's solute by default (the protein and any ligand),
        the protein where the scope is ``all``.

    Output
    ------
    ``image_check.dat``: the closest approach in nm, one value per frame;
    the finding in ``options.json`` under ``findings.solute_image``.
    """

    name = "image_check"
    description = "Closest approach to its own periodic image"
    time_series = True
    #: Images exist only in a periodic box.
    requires_periodic_box = True
    #: And there must be a solute to have them.
    requires_solute = True

    def __init__(self, *, cutoff_nm: float = CUTOFF_NM, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        cutoff_nm = float(cutoff_nm)
        if not cutoff_nm > 0:
            raise StudyError(
                f"`cutoff_nm` is a distance, more than 0 nm; {cutoff_nm:g} was given.",
                code="analysis.option.out_of_range", option="cutoff_nm",
                given=cutoff_nm, frames=None)
        self.cutoff_nm: float = cutoff_nm
        self.options.update(cutoff_nm=self.cutoff_nm)

    _traj_for_plot: md.Trajectory | None = None

    def compute(self, traj: md.Trajectory) -> np.ndarray:
        atoms = (traj.topology.select(self.selection) if self.selection
                 else traj.topology.select("protein"))
        if len(atoms) == 0:
            raise StudyError(
                f"Selection {self.selection or 'protein'!r} matched no atoms, so "
                "there is no solute to check against its images.",
                code="analysis.selection.empty",
                expression=self.selection or "protein", role="solute")
        distances, record = solute_image_record(traj, atoms, cutoff_nm=self.cutoff_nm)
        self.findings["solute_image"] = record
        # Read by this analysis as by every other.
        self.solute_image = record
        return distances

    def run(self, traj: md.Trajectory):
        self._traj_for_plot = traj
        return super().run(traj)

    def plot(self, result: np.ndarray, ax: plt.Axes) -> None:
        x = (self.frame_axis(self._traj_for_plot)[0]
             if self._traj_for_plot is not None else np.arange(result.shape[0]))
        ax.plot(x, result, linewidth=1.2)
        ax.axhline(self.cutoff_nm, color=colour("GUIDE"), linestyle="--", linewidth=0.9,
                   label=f"cutoff, {self.cutoff_nm:g} nm")
        ax.axhline(CLOSE_NM, color=colour("GUIDE"), linestyle=":", linewidth=0.9,
                   label=f"{CLOSE_NM:g} nm")
        ax.set_ylim(bottom=0)
        ax.legend(loc="best", fontsize=7.5)

    def default_xlabel(self) -> str | None:
        if self._traj_for_plot is None:
            return "Frame"
        return self.frame_axis(self._traj_for_plot)[1]

    def default_ylabel(self) -> str | None:
        return "Closest approach to its image (nm)"


register_analysis(SoluteImage.name, SoluteImage)
