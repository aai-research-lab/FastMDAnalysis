# FastMDAnalysis/src/fastmdanalysis/analysis/base.py
"""
BaseAnalysis Module

Provides common functionality for all analysis types in FastMDAnalysis:
- consistent output directory handling,
- saving data tables (.dat),
- saving plots (.png),
- a shared AnalysisError exception.

All analysis modules should subclass BaseAnalysis and implement run() and plot().
"""

from __future__ import annotations

from typing import Any, Optional
import logging
from pathlib import Path

import numpy as np

__all__ = ["BaseAnalysis", "AnalysisError"]

logger = logging.getLogger(__name__)


class AnalysisError(Exception):
    """Custom exception class for analysis errors."""
    pass


class BaseAnalysis:
    """
    BaseAnalysis with common operations for analysis modules.

    Attributes
    ----------
    traj : mdtraj.Trajectory
        The MDTraj Trajectory object to analyze.
    output : str
        Name of the output directory for this analysis.
    outdir : pathlib.Path
        Directory where data and plots are saved.
    results : dict
        Container for computed results/artifacts.
    data : Any
        Primary data product (shape/type depends on analysis).
    """

    def __init__(
        self,
        trajectory,
        output: Optional[str] = None,
        *,
        save_data: bool = True,
        store_results: bool = True,
        **kwargs,
    ):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MDTraj Trajectory to analyze.
        output : str, optional
            Output directory name. Defaults to "<analysis>_output".
        kwargs : dict
            Reserved for future extensibility.
        """
        self.traj = trajectory
        default_name = self.__class__.__name__.replace("Analysis", "").lower() + "_output"
        self.output = output or default_name
        self.outdir = Path(self.output)
        self.outdir.mkdir(parents=True, exist_ok=True)

        self.save_data = bool(save_data)
        self.store_results = bool(store_results)

        self.results: dict = {}
        self.data: Any = None

        logger.debug("Initialized %s(outdir=%s)", self.__class__.__name__, str(self.outdir))

    # --------------------------------------------------------------------- I/O

    def _save_plot(self, fig, key: str, *, filename: Optional[str] = None, dpi: int = 300) -> Path:
        """
        Save a matplotlib figure (.png) into the analysis output directory.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure to save.
        key : str
            Logical base name (used if filename is not given), e.g., "rmsd".
        filename : str, optional
            Exact filename (e.g., "my_rmsd.png"). If not provided, use f"{key}.png".
        dpi : int
            Figure DPI for rasterization (default: 300).

        Returns
        -------
        pathlib.Path
            Path to the saved image.
        """
        if filename is None:
            filename = f"{key}.png"
        elif not filename.endswith(".png"):
            # Enforce .png for consistency
            filename = f"{filename}.png"

        path = self.outdir / filename
        fig.savefig(path, bbox_inches="tight", dpi=dpi)
        logger.info("Plot saved to %s", path)
        return path

    def _save_data(
        self,
        data: Any,
        filename: str,
        *,
        header: Optional[str] = None,
        fmt: Optional[str] = None,
        delimiter: str = " "
    ) -> Path:
        """
        Save data as a .dat table in the analysis output directory.

        Parameters
        ----------
        data : array-like or object
            Data to save. If array-like, will be written via numpy.savetxt.
            Non-array data will be stringified and written as-is.
        filename : str
            Base name without extension (".dat" is appended).
        header : str, optional
            Header line (without comment symbol). If None, generated heuristically.
        fmt : str, optional
            Numpy savetxt format (e.g., "%.6f", "%d"). If None, inferred from dtype when possible.
        delimiter : str
            Column delimiter (default: space).

        Returns
        -------
        pathlib.Path
            Path to the saved .dat file.
        """
        path = self.outdir / f"{filename}.dat"

        if not self.save_data:
            return None

        # Try numpy path first
        try:
            arr = np.asarray(data)
            if arr.ndim == 1:
                arr = arr.reshape(-1, 1)

            # Header inference
            if header is None:
                if arr.ndim == 2:
                    header = " ".join([f"col{i+1}" for i in range(arr.shape[1])])
                else:
                    header = "data"

            # Format inference
            if fmt is None:
                if np.issubdtype(arr.dtype, np.integer):
                    fmt = "%d"
                elif np.issubdtype(arr.dtype, np.floating):
                    fmt = "%.6f"
                else:
                    # Fallback to generic if mixed/object dtype
                    fmt = "%s"

            np.savetxt(
                path,
                arr,
                fmt=fmt,
                delimiter=delimiter,
                header=header,
                comments="",  # no leading '#'
            )
            logger.info("Data saved to %s", path)
            return path

        except Exception:
            # Fallback: write generic string representation
            with open(path, "w") as f:
                if header:
                    f.write(f"{header}\n")
                f.write(str(data))
            logger.info("Non-array data saved to %s", path)
            return path

    # ----------------------------------------------------------------- Abstracts

    def run(self):
        """Perform the analysis. Must be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement the run() method.")

    def plot(self, *args, **kwargs):
        """Generate plots for the analysis. Should be implemented by subclasses."""
        raise NotImplementedError("Subclasses must implement the plot() method.")

