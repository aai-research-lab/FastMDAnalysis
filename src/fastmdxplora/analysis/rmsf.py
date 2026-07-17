"""Root-Mean-Square Fluctuation (RMSF).

Per-residue (default) or per-atom RMSF over the trajectory. The
trajectory is optionally superposed onto a reference (frame 0 or a user-chosen
reference) using the selected atom subset to remove rigid-body motion;
the per-atom RMSF is then the standard deviation of each atom's position
around its mean. The per-residue RMSF reduces per-atom RMSF to one value
per residue by averaging over the atoms that belong to each residue.

This is the standard "flexibility profile" plot used in nearly every MD
publication — peaks indicate flexible loops/termini, troughs indicate
rigid secondary structure.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.analysis.orchestrator import register_analysis


class RMSF(Analysis):
    """Per-residue root-mean-square fluctuation.

    Parameters
    ----------
    ref : int, default 0
        Reference frame for the alignment (superposition) step. The choice
        affects only the bookkeeping; fluctuations are measured relative
        to each atom's mean position over the full trajectory, which is
        invariant under rigid-body alignment.
    per_residue : bool, default True
        If True, collapse the per-atom RMSF down to one value per residue
        by averaging over the residue's atoms. If False, return the
        per-atom array (one value per selected atom).
    align : bool, default True
        If True, superpose frames before calculating fluctuations. Set False
        to calculate directly around the mean positions, preserving
        rigid-body motion for legacy workflows.
    selection : str, optional
        MDTraj atom selection. Defaults to ``"name CA"`` (alpha carbons)
        for protein analysis.
    **kwargs
        Standard base-class options.

    Output
    ------
    A two-column ``rmsf.dat`` (residue_index, rmsf_nm) when ``per_residue``
    is True, or (atom_index, rmsf_nm) when False.

    Examples
    --------
    Standard per-residue plot for proteins::

        rmsf = RMSF()
        rmsf.run(trajectory)

    Per-atom on backbone heavy atoms::

        rmsf = RMSF(per_residue=False, selection="backbone and not element H")
    """

    name = "rmsf"
    description = "Root-mean-square fluctuation"
    default_selection = "name CA"

    def __init__(
        self,
        *,
        ref: int = 0,
        per_residue: bool = True,
        align: bool = True,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        self.ref: int = int(ref)
        self.per_residue: bool = bool(per_residue)
        self.align: bool = bool(align)
        self.options.update(
            ref=self.ref, per_residue=self.per_residue, align=self.align
        )

    def compute(self, traj: md.Trajectory) -> np.ndarray:
        """Compute the RMSF.

        Returns
        -------
        np.ndarray, shape (N, 2)
            Two columns: index (residue or atom number) and RMSF in nm.
        """
        atom_idx = self.select_atoms(traj)

        # Alignment is the normal default: it removes rigid-body translation
        # and rotation so residual variance is conformational. Some legacy
        # workflows intentionally retain this motion around the mean.
        n = traj.n_frames
        ref = self.ref if self.ref >= 0 else n + self.ref
        if not (0 <= ref < n):
            raise ValueError(
                f"Reference frame {self.ref} is out of range for trajectory "
                f"with {n} frames."
            )

        # MDTraj's superpose mutates its target trajectory. Work on a copy so
        # running RMSF cannot alter coordinates seen by later analyses.
        aligned = (
            traj.slice(range(traj.n_frames), copy=True).superpose(
                traj, frame=ref, atom_indices=atom_idx
            )
            if self.align else traj
        )

        # Per-atom RMSF on the selected atoms only:
        # rmsf[i] = sqrt(mean over frames of ||r_i(t) - <r_i>||^2)
        xyz = aligned.xyz[:, atom_idx, :]
        mean_xyz = xyz.mean(axis=0)
        disp = xyz - mean_xyz
        per_atom = np.sqrt(np.mean(np.sum(disp * disp, axis=2), axis=0))

        if not self.per_residue:
            # Return atom_serial, rmsf
            atoms = [traj.topology.atom(int(i)) for i in atom_idx]
            atom_serials = np.array([a.serial for a in atoms])
            return np.column_stack([atom_serials, per_atom]).astype(np.float64)

        # Collapse to per-residue by averaging over each residue's atoms in
        # the selection. Preserve residue order as encountered.
        atoms = [traj.topology.atom(int(i)) for i in atom_idx]
        residues: dict[int, list[float]] = {}
        for atom, val in zip(atoms, per_atom):
            residues.setdefault(atom.residue.index, []).append(float(val))

        # Use residue.resSeq (PDB numbering) for the x axis when available;
        # fall back to topology index otherwise.
        rows: list[tuple[int, float]] = []
        for ridx in sorted(residues):
            res = traj.topology.residue(ridx)
            try:
                label = int(res.resSeq)
            except (AttributeError, TypeError):
                label = ridx
            rows.append((label, float(np.mean(residues[ridx]))))

        return np.array(rows, dtype=np.float64)

    def plot(self, result: np.ndarray, ax: plt.Axes) -> None:
        x = result[:, 0]
        y = result[:, 1]
        ax.plot(x, y, linewidth=1.4, marker="o", markersize=3, markeredgewidth=0)
        ax.fill_between(x, 0, y, alpha=0.12)

    def default_xlabel(self) -> str | None:
        return "Residue" if self.per_residue else "Atom serial"

    def default_ylabel(self) -> str | None:
        return "RMSF (nm)"


register_analysis(RMSF.name, RMSF)
