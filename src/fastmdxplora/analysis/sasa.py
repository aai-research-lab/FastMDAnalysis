"""Solvent-Accessible Surface Area (SASA).

Per-frame SASA computed with the Shrake-Rupley rolling-sphere algorithm
(MDTraj's :func:`mdtraj.shrake_rupley`), for the whole molecule, for each
residue per frame, or as each residue's mean over the run -- which is the
summary that says which residues are buried. Outputs the total SASA time
series and, optionally, a per-residue heatmap that shows which residues
become exposed/buried over the simulation.

SASA is a sensitive probe of conformational changes that involve burial
or exposure of hydrophobic surfaces — it can detect folding/unfolding
events, partial unfolding of loops, and binding/unbinding transitions
that don't necessarily show up in RMSD.

References
----------
Shrake, A.; Rupley, J. A. *J. Mol. Biol.* **1973**, 79, 351.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.analysis.orchestrator import register_analysis


#: What a SASA run reports. ``total`` is the whole molecule per frame,
#: ``residue`` every residue per frame, and ``average_residue`` each residue's
#: mean over the run -- which is the summary somebody reads to find out what
#: is buried.
VALID_MODES = ("total", "residue", "average_residue")


class SASA(Analysis):
    """Solvent-accessible surface area.

    Parameters
    ----------
    mode : {"total", "residue", "average_residue"}, default "total"
        ``"total"`` returns one value per frame (sum over all atoms).
        ``"residue"`` returns a per-residue SASA matrix (n_frames × n_residues).
    probe_radius : float, default 0.14
        Probe (solvent) radius in nm. The default 0.14 nm is the water
        radius and is the standard choice for biomolecular SASA.
    n_sphere_points : int, default 960
        Number of points on the unit sphere for the Shrake-Rupley rolling
        ball. Higher is more accurate but slower. 960 is MDTraj's default
        and provides ~1% precision.
    **kwargs
        Standard base-class options.

    Output
    ------
    ``sasa.dat`` — CSV. Either ``frame, sasa_nm2`` (total) or
    ``frame, residue, sasa_nm2`` (per residue, long format).
    ``sasa.png`` — Time series (total) or heatmap (residue).
    """

    name = "sasa"
    description = "Solvent-accessible surface area"
    default_selection = None

    def __init__(
        self,
        *,
        mode: str = "total",
        probe_radius: float = 0.14,
        n_sphere_points: int = 960,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        mode = str(mode).lower()
        if mode not in VALID_MODES:
            raise ValueError(
                f"SASA mode must be one of {VALID_MODES}; got {mode!r}"
            )
        self.mode: str = mode
        self.probe_radius: float = float(probe_radius)
        self.n_sphere_points: int = int(n_sphere_points)
        self.options.update(
            mode=self.mode,
            probe_radius=self.probe_radius,
            n_sphere_points=self.n_sphere_points,
        )

    def compute(self, traj: md.Trajectory) -> pd.DataFrame:
        """Compute SASA per frame.

        Returns
        -------
        pandas.DataFrame
            ``mode="total"``: columns ``frame, sasa_nm2``.
            ``mode="residue"``: columns ``frame, residue, sasa_nm2`` (long form).
        """
        # Restrict to the selected atoms (e.g. protein/solute) before
        # computing SASA — solvent should not contribute to the solute's
        # accessible surface area.
        atom_idx = self.select_atoms(traj)
        if len(atom_idx) < traj.n_atoms:
            traj = traj.atom_slice(atom_idx)

        mode_arg = "atom" if self.mode == "total" else "residue"
        sasa = md.shrake_rupley(
            traj,
            probe_radius=self.probe_radius,
            n_sphere_points=self.n_sphere_points,
            mode=mode_arg,
        )

        if self.mode == "total":
            total = sasa.sum(axis=1)
            return pd.DataFrame(
                {"frame": np.arange(traj.n_frames), "sasa_nm2": total}
            )

        if self.mode == "average_residue":
            # The mean exposure of each residue over the whole run, which is
            # the summary somebody actually reads: which residues are buried
            # and which are on the surface. The per-frame matrix contains it,
            # but reading it off a heatmap by eye is not the same as having
            # it. Version 1 wrote all three from one run.
            residues = list(traj.topology.residues)
            try:
                labels = np.array([int(r.resSeq) for r in residues])
            except (AttributeError, TypeError):
                labels = np.array([r.index for r in residues])
            return pd.DataFrame({
                "residue": labels,
                "mean_sasa_nm2": sasa.mean(axis=0),
                # The spread matters: a residue at 1.0 every frame and one
                # alternating between 0 and 2 have the same mean and are not
                # the same thing.
                "std_sasa_nm2": sasa.std(axis=0),
            })

        # Per-residue: build a long-form table. Residue labels = resSeq
        # (PDB numbering) when available.
        residues = list(traj.topology.residues)
        try:
            labels = np.array([int(r.resSeq) for r in residues])
        except (AttributeError, TypeError):
            labels = np.array([r.index for r in residues])

        n_frames, n_res = sasa.shape
        frames = np.repeat(np.arange(n_frames), n_res)
        residue_col = np.tile(labels, n_frames)
        return pd.DataFrame(
            {
                "frame": frames,
                "residue": residue_col,
                "sasa_nm2": sasa.flatten(),
            }
        )

    def plot(self, result: pd.DataFrame, ax: plt.Axes) -> None:
        if self.mode == "total":
            x, _ = self.frame_axis_for_plot(self._traj_for_plot, len(result))
            ax.plot(x, result["sasa_nm2"].to_numpy(), linewidth=1.4)
            ax.fill_between(x, 0, result["sasa_nm2"].to_numpy(), alpha=0.15)
        elif self.mode == "average_residue":
            # A bar per residue, with the spread over the run drawn on it. A
            # residue at 1.0 every frame and one alternating between 0 and 2
            # have the same mean, and a bar chart without the spread says they
            # are the same.
            residues = result["residue"].to_numpy()
            means = result["mean_sasa_nm2"].to_numpy()
            ax.bar(residues, means,
                   yerr=result["std_sasa_nm2"].to_numpy(),
                   color="#3a7ca5", error_kw={"ecolor": "#1f4257",
                                              "elinewidth": 0.8, "capsize": 2})
            ax.set_ylim(bottom=0)
        else:
            # Per-residue heatmap: pivot long-form -> (residue × frame)
            grid = result.pivot(
                index="residue", columns="frame", values="sasa_nm2"
            ).to_numpy()
            im = ax.imshow(
                grid,
                aspect="auto",
                origin="lower",
                cmap="viridis",
                interpolation="nearest",
            )
            ax.figure.colorbar(im, ax=ax, label="SASA (nm²)", shrink=0.85)

    def save_data(self, result: pd.DataFrame, path) -> Any:
        """Write the table, and for a per-residue run the average beside it.

        That run already contains every number the average needs, so computing
        the surface a second time to get it would cost minutes for arithmetic.
        Version 1 wrote all three outputs from one run; this writes two, and
        the third mode exists for anyone who wants only the summary.
        """
        from pathlib import Path

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(path, index=False)

        if self.mode == "residue":
            summary = (
                result.groupby("residue")["sasa_nm2"]
                .agg(mean_sasa_nm2="mean", std_sasa_nm2="std")
                .reset_index()
            )
            summary.to_csv(
                path.parent / f"{self.name}_average_per_residue.csv",
                index=False,
            )
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
        if self.mode == "average_residue":
            return "Residue"
        if self.mode == "residue":
            return "Frame"
        if self._traj_for_plot is None:
            return "Frame"
        _, label = self.frame_axis(self._traj_for_plot)
        return label

    def default_ylabel(self) -> str | None:
        if self.mode == "average_residue":
            return "Mean SASA (nm²)"
        if self.mode == "residue":
            return "Residue (index in topology)"
        return "SASA (nm²)"


register_analysis(SASA.name, SASA)
