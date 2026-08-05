"""Ligand pose RMSD (RMSD of the ligand after protein alignment).

This is the headline protein-ligand stability metric: it measures whether the
ligand stays in its binding pose over the trajectory. Each frame is rigidly
aligned onto a reference using the **protein** atoms (so protein tumbling is
removed), and then the RMSD is computed on the **ligand** atoms of the
already-aligned coordinates. A low, flat profile means the ligand holds its
pose; a rising profile means it is drifting or unbinding.

This differs from the standard :class:`~fastmdxplora.analysis.rmsd.RMSD`,
which aligns and measures on the same atom set. Here alignment (protein) and
measurement (ligand) use different selections, which is the correct way to ask
"how much has the ligand moved *relative to the protein*".

Output is a single-column ``ligand_rmsd.dat`` of RMSD values in nanometers,
and a time-series figure.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.analysis.orchestrator import register_analysis


class LigandRMSD(Analysis):
    """Per-frame RMSD of the ligand after aligning on the protein.

    Parameters
    ----------
    ligand_resname : str
        Residue name of the ligand (e.g. ``"LIG"``). Required — this analysis
        only makes sense for a protein-ligand complex. The orchestrator
        supplies it automatically from the setup manifest.
    align_selection : str, default "protein and name CA"
        Atom selection used for the rigid-body alignment (the receptor frame).
        Cα atoms are the standard, robust choice.
    ref : int, default 0
        Reference frame. Negative indices count from the end.
    **kwargs
        Standard base-class options.

    Notes
    -----
    The ``selection`` attribute is not used for the measurement here (the
    measured atoms are always the ligand); alignment is controlled by
    ``align_selection``.
    """

    name = "ligand_rmsd"
    description = "Ligand pose RMSD (after protein alignment)"
    requires_ligand = True
    # Measurement atoms are the ligand, resolved from ligand_resname; this
    # analysis is ligand-only by nature, so it does not use scope.
    default_selection = None
    #: This works out its own atoms, so a general selection has nothing to
    #: apply to.
    honours_selection = False

    def __init__(
        self,
        *,
        ligand_resname: str | None = None,
        align_selection: str = "protein and name CA",
        ref: int = 0,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        if not ligand_resname:
            raise ValueError(
                "LigandRMSD requires `ligand_resname` (the ligand residue "
                "name, e.g. 'LIG'). This analysis applies only to "
                "protein-ligand complexes."
            )
        self.ligand_resname: str = str(ligand_resname)
        self.align_selection: str = str(align_selection)
        self.ref: int = int(ref)
        self.options.update(
            ligand_resname=self.ligand_resname,
            align_selection=self.align_selection,
            ref=self.ref,
        )

    def compute(self, traj: md.Trajectory) -> np.ndarray:
        """Compute per-frame ligand RMSD after protein alignment.

        Returns
        -------
        np.ndarray of shape (n_frames,)
            Ligand RMSD in nanometers.
        """
        ligand_idx = traj.topology.select(f"resname {self.ligand_resname}")
        if len(ligand_idx) == 0:
            raise ValueError(
                f"No atoms matched ligand resname "
                f"{self.ligand_resname!r}; cannot compute ligand RMSD."
            )
        align_idx = traj.topology.select(self.align_selection)
        if len(align_idx) == 0:
            raise ValueError(
                f"Alignment selection {self.align_selection!r} matched zero "
                f"atoms; cannot align on the protein."
            )

        n = traj.n_frames
        ref = self.ref if self.ref >= 0 else n + self.ref
        if not (0 <= ref < n):
            raise ValueError(
                f"Reference frame {self.ref} is out of range for trajectory "
                f"with {n} frames."
            )

        # Align every frame onto the reference using the PROTEIN atoms. This
        # transforms all coordinates (including the ligand) by the same
        # rigid-body operation, so the residual ligand motion is motion
        # relative to the protein frame.
        aligned = traj.superpose(traj, frame=ref, atom_indices=align_idx)

        # RMSD of the LIGAND atoms on the aligned coordinates, vs the
        # reference frame's ligand coordinates. No further alignment.
        ref_xyz = aligned.xyz[ref, ligand_idx, :]
        disps = aligned.xyz[:, ligand_idx, :] - ref_xyz
        rmsd_nm = np.sqrt(np.mean(np.sum(disps * disps, axis=2), axis=1))

        self._resolved_ref = ref
        return rmsd_nm.astype(np.float64)

    def plot(self, result: np.ndarray, ax: plt.Axes) -> None:
        x, _ = self.frame_axis_for_plot(result, self._traj_for_plot)
        ax.plot(x, result, linewidth=1.4, color="#b5651d")
        ref_x = x[self._resolved_ref]
        ax.axvline(
            ref_x, color="#888888", linestyle=":", linewidth=1.0,
            label=f"reference (frame {self._resolved_ref})",
        )
        ax.legend(loc="best")

    # Plot plumbing mirrors RMSD.
    _traj_for_plot: md.Trajectory | None = None
    _resolved_ref: int = 0

    def run(self, traj: md.Trajectory):
        self._traj_for_plot = traj
        return super().run(traj)

    def frame_axis_for_plot(
        self, result: np.ndarray, traj: md.Trajectory | None
    ) -> tuple[np.ndarray, str]:
        if traj is None:
            return np.arange(len(result)), "Frame"
        return self.frame_axis(traj)

    def default_xlabel(self) -> str | None:
        if self._traj_for_plot is None:
            return "Frame"
        _, label = self.frame_axis(self._traj_for_plot)
        return label

    def default_ylabel(self) -> str | None:
        return "Ligand RMSD (nm)"


register_analysis(LigandRMSD.name, LigandRMSD)
