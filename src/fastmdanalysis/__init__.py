# FastMDAnalysis/src/fastmdanalysis/__init__.py
"""
FastMDAnalysis – automated MD trajectory analysis.

Documentation: https://fastmdanalysis.readthedocs.io/en/latest/

FastMDAnalysis Package Initialization

This version of FastMDAnalysis allows you to instantiate a single object by providing
the trajectory and topology file paths, along with optional parameters for frame and atom selection.
Frame selection is specified as an iterable with three elements (start, stop, stride); negative indices are supported—
e.g. using frames=(-10, -1, 1) or frames=[-10, -1, 1] will select frames from (n_frames - 10) to the last frame.
An atom selection string may also be provided to use a specific subset of atoms.
All subsequent analyses (rmsd, rmsf, rg, hbonds, cluster, ss, sasa, dimred) use the pre-loaded
trajectory and default atom selection unless overridden.
"""

from __future__ import annotations

from typing import Optional, Tuple, Union, Sequence
import logging

# Optional dependency import to ensure availability at import time (not used directly here).
import mdtraj as md  # noqa: F401

from .analysis import rmsd, rmsf, rg, hbonds, cluster, ss, dimred, sasa
from .utils import load_trajectory  # Extended utility supporting multiple files.
from .utils.logging import setup_library_logging, log_run_header  # convenient re-exports

# -----------------------------------------------------------------------------
# Package version
# -----------------------------------------------------------------------------
try:
    from importlib.metadata import version as _pkg_version, PackageNotFoundError  # type: ignore
except Exception:  # pragma: no cover
    try:  # Python <3.8 backport
        from importlib_metadata import version as _pkg_version, PackageNotFoundError  # type: ignore
    except Exception:  # pragma: no cover
        _pkg_version = None  # type: ignore
        PackageNotFoundError = Exception  # type: ignore

def _resolve_version() -> str:
    # Prefer distribution name in normalized (PEP 503) form; fall back to alt casing.
    for dist_name in ("fastmdanalysis", "FastMDAnalysis"):
        try:
            if _pkg_version:
                return _pkg_version(dist_name)
        except PackageNotFoundError:
            continue
        except Exception:
            continue
    # Last resort if metadata is unavailable (e.g., source checkout)
    return "0+unknown"

__version__ = _resolve_version()

# -----------------------------------------------------------------------------
# Package logging: install a NullHandler so library users don't get warnings.
# The CLI configures handlers/levels; library users can configure logging as well.
# -----------------------------------------------------------------------------
_pkg_logger = logging.getLogger("fastmdanalysis")
if not _pkg_logger.handlers:
    _pkg_logger.addHandler(logging.NullHandler())

# -----------------------------------------------------------------------------
# Expose analysis classes.
# -----------------------------------------------------------------------------
RMSDAnalysis = rmsd.RMSDAnalysis
RMSFAnalysis = rmsf.RMSFAnalysis
RGAnalysis = rg.RGAnalysis
HBondsAnalysis = hbonds.HBondsAnalysis
ClusterAnalysis = cluster.ClusterAnalysis
SSAnalysis = ss.SSAnalysis
DimRedAnalysis = dimred.DimRedAnalysis
SASAAnalysis = sasa.SASAAnalysis

__all__ = [
    "__version__",
    "FastMDAnalysis",
    "RMSDAnalysis",
    "RMSFAnalysis",
    "RGAnalysis",
    "HBondsAnalysis",
    "ClusterAnalysis",
    "SSAnalysis",
    "DimRedAnalysis",
    "SASAAnalysis",
    "load_trajectory",
    "setup_library_logging",
    "log_run_header",
]


def _normalize_frames(
    frames: Optional[Union[Sequence[Union[int, None]], Tuple[Optional[int], Optional[int], Optional[int]]]]
) -> Optional[Tuple[Optional[int], Optional[int], int]]:
    """
    Normalize (start, stop, stride) for slicing.

    - Accepts None or a 3-tuple/list.
    - Converts elements to int or None (for start/stop).
    - Ensures stride is a positive integer (defaults to 1 if None/0).

    We intentionally avoid comparing start/stop to integers; Python slicing
    fully supports None and negative indices.
    """
    if frames is None:
        return None
    if not isinstance(frames, (list, tuple)) or len(frames) != 3:
        raise TypeError("frames must be None or a 3-tuple/list: (start, stop, stride)")

    start, stop, stride = frames

    def _int_or_none(x):
        if x is None:
            return None
        try:
            return int(x)
        except Exception as e:
            raise TypeError("frames elements must be int or None") from e

    start_i = _int_or_none(start)
    stop_i = _int_or_none(stop)
    stride_i = _int_or_none(stride)

    if stride_i is None or stride_i == 0:
        stride_i = 1
    if stride_i < 0:
        stride_i = -stride_i

    return (start_i, stop_i, stride_i)


class FastMDAnalysis:
    """
    Main API class for MD trajectory analysis.

    This class loads an MD trajectory from file paths and optionally subsets the trajectory (frames)
    and the set of atoms used (atom selection). These default selections are then applied to all analyses,
    although each analysis method can override them if desired.

    Parameters
    ----------
    traj_file : str
        Path to the trajectory file (e.g. "trajectory.dcd").
    top_file : str
        Path to the topology file (e.g. "topology.pdb").
    frames : iterable of three (int or None), optional
        (start, stop, stride) to subset frames. Negative indices and None are allowed.
        For example, (-10, -1, 1) selects frames from (n_frames - 10) through the last frame.
        If None, the entire trajectory is used.
    atoms : str or None, optional
        An MDTraj atom selection string (e.g., "protein" or "protein and name CA") specifying which atoms to use.
        If None, all atoms are used.

    Examples
    --------
    >>> from fastmdanalysis import FastMDAnalysis
    >>> fastmda = FastMDAnalysis("trajectory.dcd", "topology.pdb", frames=(0, None, 10), atoms="protein")
    >>> rmsd_analysis = fastmda.rmsd(ref=0)
    """

    def __init__(self, traj_file: str, top_file: str, frames=None, atoms: Optional[str] = None, keep_full_traj: bool = True):
        norm_frames = _normalize_frames(frames)

        self._default_atoms_preapplied = False
        self.default_atoms = atoms

        if not keep_full_traj:
            base_traj = load_trajectory(traj_file, top_file, frames=norm_frames, atoms=atoms)
            if atoms:
                self._default_atoms_preapplied = True
            self.full_traj = None
            self.traj = base_traj
        else:
            full_traj = load_trajectory(traj_file, top_file)
            self.full_traj = full_traj

            if norm_frames is not None:
                start, stop, stride = norm_frames
                base_traj = full_traj[start:stop:stride]
            else:
                base_traj = full_traj

            self.traj = base_traj

        # Optional: common output/figure dirs other utilities may look for
        # (these attributes are probed by the slides utility in analyze.py)
        self.figdir = getattr(self, "figdir", "figures")
        self.outdir = getattr(self, "outdir", "results")

    def _get_atoms(self, specific_atoms: Optional[str]) -> Optional[str]:
        """
        Determine the atom selection string to use; prefer per-call override,
        else fall back to the default selection provided at initialization.
        """
        if specific_atoms is not None:
            return specific_atoms
        if self._default_atoms_preapplied:
            return None
        return self.default_atoms

    # ----------------------------- Analyses -----------------------------------

    def rmsd(
        self,
        reference_frame: Optional[int] = None,
        ref: Optional[int] = None,
        atoms: Optional[str] = None,
        save_data: bool = True,
        store_results: bool = True,
        **kwargs,
    ):
        """
        Run RMSD analysis on the stored trajectory.

        Parameters
        ----------
        reference_frame / ref : int, optional
            Reference frame index for RMSD calculations (default: 0).
            Both names are accepted; `ref` overrides if both provided.
        atoms : str, optional
            Atom selection string for this analysis. If not provided, uses the default atom selection.
        kwargs : dict
            Additional keyword arguments to pass to RMSDAnalysis.

        Returns
        -------
        RMSDAnalysis
            An RMSDAnalysis instance containing the computed results.
        """
        a = self._get_atoms(atoms)
        rf = ref if ref is not None else (reference_frame if reference_frame is not None else 0)
        analysis = RMSDAnalysis(
            self.traj,
            reference_frame=rf,
            atoms=a,
            save_data=save_data,
            store_results=store_results,
            **kwargs,
        )
        analysis.run()
        return analysis

    def rmsf(
        self,
        atoms: Optional[str] = None,
        save_data: bool = True,
        store_results: bool = True,
        **kwargs,
    ):
        """Run RMSF analysis on the stored trajectory."""
        a = self._get_atoms(atoms)
        analysis = RMSFAnalysis(
            self.traj,
            atoms=a,
            save_data=save_data,
            store_results=store_results,
            **kwargs,
        )
        analysis.run()
        return analysis

    def rg(
        self,
        atoms: Optional[str] = None,
        save_data: bool = True,
        store_results: bool = True,
        **kwargs,
    ):
        """Run Radius of Gyration (RG) analysis."""
        a = self._get_atoms(atoms)
        analysis = RGAnalysis(
            self.traj,
            atoms=a,
            save_data=save_data,
            store_results=store_results,
            **kwargs,
        )
        analysis.run()
        return analysis

    def hbonds(self, atoms: Optional[str] = None, **kwargs):
        """Run Hydrogen Bonds (HBonds) analysis."""
        a = self._get_atoms(atoms)
        analysis = HBondsAnalysis(self.traj, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def cluster(
        self,
        methods="dbscan",
        eps: float = 0.5,
        min_samples: int = 5,
        n_clusters: Optional[int] = None,
        atoms: Optional[str] = None,
        plot_style: str = "full",
        combined_plot_name: str = "cluster_combined",
        save_data: bool = True,
        store_results: bool = True,
        feature_mode: str = "coordinates",
        **kwargs
    ):
        """Run clustering analysis on the stored trajectory."""
        a = self._get_atoms(atoms)
        analysis = ClusterAnalysis(
            self.traj,
            methods=methods,
            eps=eps,
            min_samples=min_samples,
            n_clusters=n_clusters,
            atoms=a,
            plot_style=plot_style,
            combined_plot_name=combined_plot_name,
            save_data=save_data,
            store_results=store_results,
            feature_mode=feature_mode,
            **kwargs
        )
        analysis.run()
        return analysis

    def ss(self, atoms: Optional[str] = None, **kwargs):
        """Run Secondary Structure (SS) analysis on the stored trajectory."""
        a = self._get_atoms(atoms)
        analysis = SSAnalysis(self.traj, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def sasa(self, probe_radius: float = 0.14, atoms: Optional[str] = None, **kwargs):
        """Run Solvent Accessible Surface Area (SASA) analysis on the stored trajectory."""
        a = self._get_atoms(atoms)
        analysis = SASAAnalysis(self.traj, probe_radius=probe_radius, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def dimred(self, methods="all", atoms: Optional[str] = None, **kwargs):
        """Run dimensionality reduction analysis on the stored trajectory."""
        a = self._get_atoms(atoms)
        analysis = DimRedAnalysis(self.traj, methods=methods, atoms=a, **kwargs)
        analysis.run()
        return analysis


# Bind analyze method to FastMDAnalysis
from .analysis.analyze import analyze as _analyze  # noqa: E402
FastMDAnalysis.analyze = _analyze  # adds the bound method
