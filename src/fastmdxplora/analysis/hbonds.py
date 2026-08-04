"""Hydrogen bond analysis.

Identifies hydrogen bonds across the trajectory using either the
Baker-Hubbard or Wernet-Nilsson geometric criteria, and produces two
outputs: a per-frame H-bond count time series, and a long-form table of
which donor/acceptor pairs participated in H-bonds, with occupancy
fractions.

Methods
-------
**Baker-Hubbard** (default) — D–H···A angle > 120°, H···A distance < 2.5 Å,
applied with a 10% occupancy threshold by default. Standard choice for
protein backbone H-bonds.

**Wernet-Nilsson** — Geometric criterion designed for water; the cutoff
distance is dynamically adjusted by the D–H–A angle. Useful for protein-
water and water-water bonds in solvated systems.

References
----------
Baker, E.; Hubbard, R. *Prog. Biophys. Mol. Biol.* **1984**, 44, 97.
Wernet, P. et al. *Science* **2004**, 304, 995.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.utils.logging import get_logger
from fastmdxplora.analysis.orchestrator import register_analysis

logger = get_logger("analysis.hbonds")


class HBonds(Analysis):
    """Hydrogen bond identification and counting.

    Parameters
    ----------
    method : {"baker_hubbard", "wernet_nilsson"}, default "baker_hubbard"
        Geometric criterion. Baker-Hubbard is the conventional choice for
        protein backbone; Wernet-Nilsson is better for water hydrogen
        bonding.
    freq : float, default 0.1
        Occupancy above which a bond is called persistent. Every bond is
        counted in the per-frame series whatever its occupancy; this decides
        only how many are reported as persistent, alongside the total number
        found, in ``n_persistent_bonds``. Has no effect on Wernet-Nilsson.
    candidate_freq : float, default 0.0
        Occupancy threshold Baker-Hubbard applies when proposing which bonds
        to evaluate. Zero proposes every bond seen in any frame, which is
        what a per-frame count needs: a bond present in five per cent of
        frames is present in those frames, and a threshold applied here would
        drop it from all of them. Raise it only to restrict the series to
        bonds that persist. Has no effect on Wernet-Nilsson.
    count_multiplier : int, default 1
        Multiplies the per-frame count, and is a compatibility device rather
        than a convention. MDTraj enumerates each hydrogen bond once, as one
        donor-hydrogen-acceptor triplet, so a multiplier of 2 reports twice
        the number of bonds present. It exists to reproduce the counts
        published by version 1 and warns when it is not 1. Must be at least
        1. Has no effect on Wernet-Nilsson.
    **kwargs
        Standard base-class options.

    Output
    ------
    ``hbonds.dat`` — CSV with frame-by-frame H-bond counts.
    ``hbonds.png`` — Time-series of H-bond count per frame.

    Notes
    -----
    The compute() method returns a pandas DataFrame:
    ``frame, n_hbonds`` with one row per frame.
    """

    name = "hbonds"
    description = "Hydrogen bonds"
    default_selection = None  # MDTraj selects donors/acceptors automatically

    def __init__(
        self,
        *,
        method: str = "baker_hubbard",
        freq: float = 0.1,
        candidate_freq: float = 0.0,
        count_multiplier: int = 1,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        method = str(method).lower()
        if method not in ("baker_hubbard", "wernet_nilsson"):
            raise ValueError(
                f"HBonds method must be 'baker_hubbard' or 'wernet_nilsson'; "
                f"got {method!r}"
            )
        self.method: str = method
        self.freq: float = float(freq)
        # Every bond that occurs is evaluated at every frame. Proposing only
        # bonds that clear an occupancy threshold first would leave a bond
        # present in five per cent of frames out of all of them, including the
        # frames where it is there -- and the series is labelled as the number
        # of hydrogen bonds per frame, not the number of persistent ones.
        self.candidate_freq: float = (
            0.0 if candidate_freq is None else float(candidate_freq)
        )  # None is still accepted: it meant "same as freq" before.
        self.count_multiplier: int = int(count_multiplier)
        if self.count_multiplier < 1:
            raise ValueError("count_multiplier must be at least 1")
        if self.count_multiplier != 1:
            logger.warning(
                "hbonds: multiplying the per-frame count by %d. MDTraj lists "
                "each hydrogen bond once, as one donor-hydrogen-acceptor "
                "triplet, so the reported series is %d times the number of "
                "hydrogen bonds present. This reproduces an earlier count; it "
                "is not a convention.",
                self.count_multiplier, self.count_multiplier,
            )
        self.options.update(
            method=self.method,
            freq=self.freq,
            candidate_freq=self.candidate_freq,
            count_multiplier=self.count_multiplier,
        )

    def compute(self, traj: md.Trajectory) -> pd.DataFrame:
        """Compute per-frame H-bond counts.

        Returns
        -------
        pandas.DataFrame
            Columns ``frame, n_hbonds``. One row per frame.
        """
        # Restrict to the selected atoms (e.g. protein/solute) before H-bond
        # detection, so on a solvated system we count the solute's hydrogen
        # bonds rather than scanning thousands of waters.
        atom_idx = self.select_atoms(traj)
        if len(atom_idx) < traj.n_atoms:
            traj = traj.atom_slice(atom_idx)

        # MDTraj's H-bond functions need explicit bond connectivity in the
        # topology to identify donor-H pairs. Most PDB-loaded trajectories
        # have this; some programmatically built or unusual ones don't.
        # Create standard bonds if missing — a no-op when already present.
        if traj.topology.n_bonds == 0:
            traj.topology.create_standard_bonds()

        if self.method == "wernet_nilsson":
            # Returns a list (one per frame) of (donor, H, acceptor) triplets.
            per_frame = md.wernet_nilsson(traj)
            counts = np.array([len(bonds) for bonds in per_frame], dtype=int)
        else:
            # Baker-Hubbard returns aggregated bonds present above `freq`
            # threshold. To produce a per-frame count we re-evaluate the
            # bonds frame by frame using its definitions (an O(n_frames)
            # loop, but MDTraj's vectorized distance/angle is fast).
            bonds = md.baker_hubbard(
                traj,
                freq=self.candidate_freq,
                exclude_water=True,
                periodic=False,
            )
            counts, occupancy = _per_frame_baker_hubbard(traj, bonds)
            # Occupancy is already known per bond from the per-frame pass, so
            # reporting how many clear the threshold costs nothing. Without
            # this, freq decided which bonds were proposed and nothing else,
            # and once every bond is proposed it would decide nothing at all.
            self.options["n_persistent_bonds"] = int(
                (occupancy >= self.freq).sum()
            )
            self.options["n_candidate_bonds"] = int(len(bonds))
            counts = self.count_multiplier * counts

        return pd.DataFrame({"frame": np.arange(traj.n_frames), "n_hbonds": counts})

    def plot(self, result: pd.DataFrame, ax: plt.Axes) -> None:
        x, _ = self.frame_axis_for_plot(self._traj_for_plot, len(result))
        ax.plot(x, result["n_hbonds"].to_numpy(), linewidth=1.4)
        ax.fill_between(x, 0, result["n_hbonds"].to_numpy(), alpha=0.15)

    def save_data(self, result: pd.DataFrame, path) -> Any:
        from pathlib import Path

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(path, index=False)
        return path

    _traj_for_plot: md.Trajectory | None = None

    def run(self, traj: md.Trajectory):
        self._traj_for_plot = traj
        return super().run(traj)

    def frame_axis_for_plot(
        self, traj: md.Trajectory | None, n_points: int
    ) -> tuple[np.ndarray, str]:
        if traj is None:
            return np.arange(n_points), "Frame"
        return self.frame_axis(traj)

    def default_xlabel(self) -> str | None:
        if self._traj_for_plot is None:
            return "Frame"
        _, label = self.frame_axis(self._traj_for_plot)
        return label

    def default_ylabel(self) -> str | None:
        return "Number of hydrogen bonds"


def _per_frame_baker_hubbard(
    traj: md.Trajectory, bonds: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Recompute per-frame occupancy for an aggregated Baker-Hubbard set.

    ``bonds`` is the (n_bonds, 3) [donor_idx, H_idx, acceptor_idx] array
    returned by ``md.baker_hubbard``. We evaluate each candidate bond at
    every frame against the standard Baker-Hubbard cutoffs (H-A distance
    < 0.25 nm AND D-H-A angle > 120°).

    Returns the count per frame and the occupancy of each bond, which the
    same pass already establishes.
    """
    if len(bonds) == 0:
        return np.zeros(traj.n_frames, dtype=int), np.zeros(0, dtype=float)

    # Distances H-A
    h_a_pairs = bonds[:, [1, 2]]
    distances = md.compute_distances(traj, h_a_pairs, periodic=False)

    # Angles D-H-A (in radians)
    d_h_a_triples = bonds[:, [0, 1, 2]]
    angles = md.compute_angles(traj, d_h_a_triples, periodic=False)

    # Mask: distance < 0.25 nm AND angle > 120° (2.0944 rad)
    cutoff_dist = 0.25
    cutoff_angle_rad = np.deg2rad(120.0)
    present = (distances < cutoff_dist) & (angles > cutoff_angle_rad)
    return present.sum(axis=1).astype(int), present.mean(axis=0)


register_analysis(HBonds.name, HBonds)
