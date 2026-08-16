"""The radial distribution function between two selections.

How the density of one group varies with distance from the other, relative
to what it would be if the two ignored each other. A value of one means no
structure at that separation; the first peak is the first solvation shell,
and where it sits and how tall it stands are quantities with published
values for every common water model, which makes this one of the few things
a simulation reports that can be checked against a number somebody else
measured.

**A radial distribution function needs a box.** The normalisation divides
by the bulk density, and the bulk density is the pair count divided by the
volume the system occupies. A trajectory carrying no unit cell has no such
volume, and the curve that comes out of assuming one is a histogram wearing
the units of a g(r).

**And it stops at half the box.** Beyond that separation the minimum-image
convention no longer supplies a complete shell: part of every sphere lies
outside the periodic cell and is counted from the wrong image or not at
all, so g(r) sags towards zero for reasons that have nothing to do with the
liquid. The curve stays smooth and plausible while it does this, which is
why the range is capped here rather than left to the person reading the
plot.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.analysis.orchestrator import register_analysis

#: Pairs beyond this many are subsampled before the histogram. A protein
#: against every water oxygen is tens of millions of pairs per frame, which
#: is minutes of arithmetic for a curve that a random tenth of them
#: resolves to within the width of its own line.
MAX_PAIRS = 400_000


class RadialDistribution(Analysis):
    """g(r) between two atom selections.

    Parameters
    ----------
    selection_a : str, default "protein"
        The first group. Distances are measured from these atoms.
    selection_b : str, default "water and name O"
        The second group. Water oxygens by default, which is the pairing
        a solvation shell is usually read from.
    r_max : float, optional
        Where to stop, in nm. Capped at half the smallest box dimension,
        and defaults to it.
    bin_width : float, default 0.005
        Histogram bin width in nm.
    **kwargs
        Standard base-class options.

    Output
    ------
    ``rdf.dat`` -- two columns, separation in nm and g(r).
    """

    name = "rdf"
    description = "Radial distribution function between two selections"
    #: The two selections are the analysis's own, and a scope selection
    #: would name a third thing that is neither of them.
    honours_selection = False
    default_selection = None
    time_series = False
    #: Without a unit cell there is no bulk density to normalise against.
    requires_periodic_box = True

    def __init__(
        self,
        *,
        selection_a: str = "protein",
        selection_b: str = "water and name O",
        r_max: float | None = None,
        bin_width: float = 0.005,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        self.selection_a = str(selection_a)
        self.selection_b = str(selection_b)
        self.r_max = None if r_max is None else float(r_max)
        self.bin_width = float(bin_width)
        self.options.update(
            selection_a=self.selection_a,
            selection_b=self.selection_b,
            r_max=self.r_max,
            bin_width=self.bin_width,
        )

    def compute(self, traj: md.Trajectory) -> np.ndarray:
        if traj.unitcell_lengths is None:
            raise ValueError(
                "This trajectory carries no unit cell, so there is no volume "
                "to take a bulk density from and no g(r) to normalise "
                "against it. A trajectory stripped of its box, or one built "
                "in vacuum, does this."
            )

        a = traj.topology.select(self.selection_a)
        b = traj.topology.select(self.selection_b)
        for name, selection, found in (
                ("selection_a", self.selection_a, a),
                ("selection_b", self.selection_b, b)):
            if len(found) == 0:
                raise ValueError(
                    f"{name} {selection!r} matched no atoms, so there is no "
                    "distribution between the two groups to report."
                )

        half_box = float(np.min(traj.unitcell_lengths)) / 2.0
        requested = self.r_max if self.r_max is not None else half_box
        r_max = min(requested, half_box)
        capped = self.r_max is not None and self.r_max > half_box

        rng = np.random.default_rng(0)
        pairs = np.array(
            [(int(i), int(j)) for i in a for j in b if int(i) != int(j)],
            dtype=int)
        subsampled = False
        if len(pairs) > MAX_PAIRS:
            keep = rng.choice(len(pairs), MAX_PAIRS, replace=False)
            pairs = pairs[np.sort(keep)]
            subsampled = True

        radii, g = md.compute_rdf(
            traj, pairs, r_range=(0.0, r_max), bin_width=self.bin_width)

        record: dict[str, Any] = {
            "pairs": int(len(pairs)),
            "r_max_nm": r_max,
            "half_box_nm": half_box,
            "subsampled": subsampled,
        }
        if capped:
            record["capped"] = (
                f"r_max was asked for at {self.r_max:g} nm and stopped at "
                f"{half_box:.3g} nm, half the smallest box dimension. Past "
                "that the minimum-image convention supplies only part of "
                "each shell, so the curve falls away for a reason belonging "
                "to the box rather than to the liquid, and it does so "
                "smoothly enough to be mistaken for structure."
            )
        if subsampled:
            record["subsampled_note"] = (
                f"{MAX_PAIRS} pairs of the {len(a) * len(b)} available, "
                "drawn once and held fixed across frames so the curve is "
                "not a different estimator in every frame."
            )

        # The first peak, which is what a published value is quoted for.
        inside = radii > 0.15
        if inside.any() and np.isfinite(g[inside]).any():
            peak = int(np.argmax(np.where(inside, g, -np.inf)))
            record["first_peak_nm"] = float(radii[peak])
            record["first_peak_height"] = float(g[peak])

        self.findings["rdf"] = record
        return np.column_stack([radii, g])

    def plot(self, result: np.ndarray, ax: plt.Axes) -> None:
        ax.plot(result[:, 0], result[:, 1], linewidth=1.2)
        ax.axhline(1.0, color="#888888", linestyle=":", linewidth=0.8)

    def default_xlabel(self) -> str | None:
        return "r (nm)"

    def default_ylabel(self) -> str | None:
        return "g(r)"


register_analysis(RadialDistribution.name, RadialDistribution)
