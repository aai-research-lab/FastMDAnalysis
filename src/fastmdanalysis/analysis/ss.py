# FastMDAnalysis/src/fastmdanalysis/analysis/ss.py
"""
SS Analysis Module

Computes secondary structure assignments for each frame using DSSP.
Saves:
  - ss_letters.dat : per-frame secondary-structure letters (comma-separated)
  - ss_numeric.dat : per-frame numeric codes (0..7) suitable for analysis
  - ss.png         : discrete heatmap of SS vs frame
  - ss_README.md   : legend for DSSP codes
  - ss_letter_codes.png : quick reference image of the legend table

The heatmap uses a discrete, high-contrast colormap so each SS letter is
easily distinguished. The colorbar tick labels display the SS letter codes.
The residue index axis is labeled with whole numbers starting at 1.

Usage:
    from fastmdanalysis import SSAnalysis
    analysis = SSAnalysis(trajectory, atoms="protein")
    analysis.run()   # computes SS, writes files, and plots by default
    analysis.plot()  # re-plot if needed with custom options
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple
from pathlib import Path

import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

from .base import BaseAnalysis, AnalysisError
from ..utils.options import OptionsForwarder
from ..utils.plotting import apply_slide_style, auto_ticks


# Order matters: numeric codes 0..7 map to these labels/colors
SS_MAP: Dict[str, int] = {
    "C": 0,  # Coil / Loop (also space maps to coil)
    " ": 0,
    "H": 1,  # Alpha helix
    "B": 2,  # Isolated beta-bridge
    "E": 3,  # Extended strand (beta sheet)
    "G": 4,  # 3-10 helix
    "I": 5,  # Pi helix
    "T": 6,  # Turn
    "S": 7,  # Bend
}
# For the colorbar (index -> label)
SS_TICKS = np.arange(0, 8, dtype=int)
SS_TICK_LABELS = ["C", "H", "B", "E", "G", "I", "T", "S"]

# Distinct, high-contrast colors for the eight states (index order as above)
SS_COLORS = [
    "#9e9e9e",  # C (gray)
    "#e41a1c",  # H (red)
    "#ff7f00",  # B (orange)
    "#377eb8",  # E (blue)
    "#4daf4a",  # G (green)
    "#984ea3",  # I (purple)
    "#f781bf",  # T (pink)
    "#a65628",  # S (brown)
]
SS_CMAP = ListedColormap(SS_COLORS)
SS_NORM = BoundaryNorm(np.arange(-0.5, 8.5, 1.0), SS_CMAP.N)


SS_CODE_ROWS = [
    ("H", "Alpha helix"),
    ("B", "Isolated beta-bridge"),
    ("E", "Extended strand (beta sheet)"),
    ("G", "3-10 helix"),
    ("I", "Pi helix"),
    ("T", "Turn"),
    ("S", "Bend"),
    ("C (or space)", "Coil / Loop (no regular secondary structure)"),
]


class SSAnalysis(BaseAnalysis):
    _ALIASES = {
        "atom_indices": "atoms",
        "selection": "atoms",
    }
    
    def __init__(
        self, 
        trajectory, 
        atoms: Optional[str] = None,
        algorithm: str = "dssp",
        mkdssp_path: Optional[str] = None,
        strict: bool = False,
        **kwargs
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        atoms : str, optional
            MDTraj selection string. If provided, SS will be computed on that subset.
            If None, all atoms are used.
            Aliases: atom_indices, selection
        algorithm : str
            Secondary structure algorithm (default: "dssp"). Currently only DSSP is supported.
        mkdssp_path : str, optional
            Path to mkdssp executable if needed (MDTraj typically auto-detects).
        strict : bool
            If True, raise errors for unknown options. If False, log warnings.
        kwargs : dict
            Passed to BaseAnalysis.
        """
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "atoms": atoms,
            "algorithm": algorithm,
            "mkdssp_path": mkdssp_path,
            "strict": strict,
        }
        analysis_opts.update(kwargs)
        
        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {"atoms", "algorithm", "mkdssp_path", "strict", "output"},
            context="ss",
            warn=warn_unknown,
        )

        atoms = resolved.get("atoms", None)
        algorithm = resolved.get("algorithm", "dssp")
        mkdssp_path = resolved.get("mkdssp_path", None)
        base_kwargs = {k: v for k, v in resolved.items() 
                      if k not in ("atoms", "algorithm", "mkdssp_path", "strict")}

        super().__init__(trajectory, **base_kwargs)
        self.atoms = atoms
        self.algorithm = algorithm.lower()
        self.mkdssp_path = mkdssp_path
        self.strict = strict
        self.data: Optional[np.ndarray] = None  # will hold letter codes (n_frames, n_residues)

    # -------------------------- helpers --------------------------

    def _subset_trajectory(self):
        """Return the trajectory (subselected by atoms if provided)."""
        if self.atoms:
            idx = self.traj.topology.select(self.atoms)
            if idx is None or len(idx) == 0:
                raise AnalysisError(f"No atoms selected using selection: '{self.atoms}'")
            return self.traj.atom_slice(idx)
        return self.traj

    @staticmethod
    def _letters_to_numeric(dssp_letters: np.ndarray) -> np.ndarray:
        """
        Map (n_frames, n_residues) array of DSSP letters to numeric codes 0..7.
        Unknown letters fall back to 0 (coil).
        """
        numeric = np.zeros_like(dssp_letters, dtype=int)
        # Vectorized map via dict: slower for truly huge arrays, but fine here
        # We iterate per frame to avoid building a giant object array
        for i in range(dssp_letters.shape[0]):
            # Convert row of single-char strings -> ints
            numeric[i, :] = [SS_MAP.get(ch, 0) for ch in dssp_letters[i, :]]
        return numeric

    def _write_letters_dat(self, dssp_letters: np.ndarray) -> Path:
        """
        Write per-frame letters to ss_letters.dat (comma-separated).
        """
        out = self.outdir / "ss_letters.dat"
        with open(out, "w") as f:
            for frame_idx in range(dssp_letters.shape[0]):
                row = ",".join(dssp_letters[frame_idx, :].tolist())
                f.write(row + "\n")
        return out

    def _write_numeric_dat(self, numeric: np.ndarray) -> Path:
        """
        Write per-frame numeric codes to ss_numeric.dat as a dense table.
        """
        # Header like: res1 res2 ... resN
        n_res = numeric.shape[1]
        header = " ".join([f"res{i+1}" for i in range(n_res)])
        return self._save_data(numeric, "ss_numeric", header=header, fmt="%d")

    def _generate_readme(self) -> Path:
        """
        Generate ss_README.md with the legend of DSSP codes.
        """
        lines = [
            "# Secondary Structure (SS) Letter Codes",
            "",
            "This document explains the DSSP secondary-structure codes used in "
            "the FastMDAnalysis SS heatmap and data tables.",
            "",
            "| Code | Description |",
            "|------|-------------|",
        ]
        for code, desc in SS_CODE_ROWS:
            lines.append(f"| {code} | {desc} |")
        text = "\n".join(lines) + "\n"
        path = self.outdir / "ss_README.md"
        with open(path, "w") as f:
            f.write(text)
        return path

    def _generate_reference_image(self) -> Path:
        """
        Create a PNG table summarizing DSSP codes.
        """
        fig, ax = plt.subplots(figsize=(10, 4.5))
        ax.axis("off")
        fig.suptitle("Secondary Structure (SS) Letter Codes", fontsize=14, fontweight="bold", y=0.95)
        fig.text(
            0.5,
            0.88,
            "DSSP codes used in the FastMDAnalysis SS heatmap.",
            ha="center",
            va="center",
            fontsize=10,
        )
        table = ax.table(
            cellText=SS_CODE_ROWS,
            colLabels=["Code", "Description"],
            colWidths=[0.22, 0.78],
            loc="center",
            cellLoc="left",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1.0, 1.3)
        for cell in table.get_celld().values():
            cell.get_text().set_ha("left")
            cell.get_text().set_va("center")
            cell.get_text().set_wrap(True)

        out = self._save_plot(fig, "ss_letter_codes")
        plt.close(fig)
        return out

    # ---------------------------- run/plot ----------------------------

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute SS assignments using DSSP; save letters and numeric tables; make plots and legend.

        Returns
        -------
        dict
            {"ss_data": <letters ndarray (n_frames, n_residues)>}
        """
        try:
            subtraj = self._subset_trajectory()

            # Some topologies need inferred bonds for proper geometry handling
            try:
                subtraj.topology.create_standard_bonds()
            except Exception:
                pass

            # md.compute_dssp returns array of shape (n_frames, n_residues) with single-char strings
            dssp_letters = md.compute_dssp(subtraj)
            if not isinstance(dssp_letters, np.ndarray) or dssp_letters.ndim != 2:
                raise AnalysisError("Unexpected DSSP output shape; expected (n_frames, n_residues).")

            self.data = dssp_letters
            self.results = {"ss_data": self.data}

            # Persist data tables
            letters_path = self._write_letters_dat(dssp_letters)
            numeric = self._letters_to_numeric(dssp_letters)
            numeric_path = self._write_numeric_dat(numeric)
            self.results.update({"ss_letters_file": letters_path, "ss_numeric_file": numeric_path})

            # Legend artifacts
            readme_path = self._generate_readme()
            codes_png_path = self._generate_reference_image()
            self.results.update({"ss_readme": readme_path, "ss_codes_plot": codes_png_path})

            # Default heatmap
            self.plot()

            return self.results
        except AnalysisError:
            raise
        except Exception as e:
            raise AnalysisError(f"SS analysis failed: {e}")

    def plot(
        self,
        data: Optional[np.ndarray] = None,
        *,
        title: str = "SS Heatmap",
        xlabel: str = "Frame",
        ylabel: str = "Residue Index",
        filename: str = "ss",
        cmap: Optional[ListedColormap] = None,
    ) -> Path:
        """
        Generate heatmap of SS vs. frame with discrete colorbar.

        Parameters
        ----------
        data
            Letter codes (n_frames, n_residues). If None, uses self.data.
        title, xlabel, ylabel
            Plot labeling.
        filename
            Base filename for the plot (without extension). Default "ss".
        cmap
            Optional colormap override; by default uses SS_CMAP.

        Returns
        -------
        Path
            Path to saved heatmap image.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No SS data available to plot. Please run analysis first.")

        if data.dtype.kind in ("U", "S", "O"):
            # Convert letters to numeric for plotting
            numeric = self._letters_to_numeric(data)
        else:
            # Assume already numeric
            numeric = np.asarray(data, dtype=int)

        # Transpose so rows -> residues, columns -> frames
        Z = numeric.T  # shape: (n_residues, n_frames)
        n_residues = Z.shape[0]

        fig, ax = plt.subplots(figsize=(12, 8))
        im = ax.imshow(
            Z,
            aspect="auto",
            interpolation="none",
            cmap=(cmap if cmap is not None else SS_CMAP),
            norm=SS_NORM,
        )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        cbar = fig.colorbar(im, ax=ax, ticks=SS_TICKS)
        cbar.set_ticklabels(SS_TICK_LABELS)
        cbar.set_label("SS Code")

        frames = np.arange(Z.shape[1], dtype=int)
        residues = np.arange(n_residues, dtype=int)
        if n_residues <= 60:
            res_ticks = residues
        else:
            auto = auto_ticks(residues, max_ticks=30, integer=True)
            res_ticks = auto.astype(int) if auto is not None and auto.size > 0 else residues[:: max(1, n_residues // 30)]

        apply_slide_style(
            ax,
            x_values=frames,
            y_ticks=res_ticks,
            integer_x=True,
            integer_y=True,
        )
        tick_font = ax.get_yticklabels()[0].get_fontsize() if ax.get_yticklabels() else None
        ax.set_yticks(res_ticks)
        ax.set_yticklabels(
            [str(int(v) + 1) for v in res_ticks],
            rotation_mode="anchor",
            fontsize=tick_font,
        )

        fig.tight_layout()

        # Save (BaseAnalysis._save_plot supports filename=)
        try:
            out = self._save_plot(fig, "ss", filename=f"{filename}.png")  # type: ignore[arg-type]
        except TypeError:
            # Backward-compat if _save_plot doesn't accept filename
            out = self._save_plot(fig, filename)

        plt.close(fig)
        return Path(out)

