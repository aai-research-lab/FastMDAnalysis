"""
FastMDAnalysis Package Initialization

Exports key analysis modules for API usage.
"""

from .analysis import (
    rmsd,
    rmsf,
    rg,
    hbonds,
    cluster,
    ss,         # renamed secondary_structure to ss
    dimred,
    sasa
)
from .utils import load_trajectory, create_dummy_trajectory

RMSDAnalysis = rmsd.RMSDAnalysis
RMSFAnalysis = rmsf.RMSFAnalysis
RGAnalysis = rg.RGAnalysis
HBondsAnalysis = hbonds.HBondsAnalysis
ClusterAnalysis = cluster.ClusterAnalysis
SSAnalysis = ss.SSAnalysis         # renamed class: SecondaryStructureAnalysis -> SSAnalysis
DimRedAnalysis = dimred.DimRedAnalysis
SASAAnalysis = sasa.SASAAnalysis

class FastMDAnalysis:
    """
    Main API class for MD trajectory analysis.
    Provides wrapper methods to run various analyses.
    """

    def __init__(self):
        # You can initialize common settings here if needed.
        pass

    def _load(self, traj, top):
        from .utils import load_trajectory
        return load_trajectory(traj, top)

    def rmsd(self, traj_path, top, output=None, **kwargs):
        import mdtraj as md
        traj = self._load(traj_path, top)
        analysis = RMSDAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

    def rmsf(self, traj_path, top, output=None, **kwargs):
        traj = self._load(traj_path, top)
        analysis = RMSFAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

    def rg(self, traj_path, top, output=None, **kwargs):
        import mdtraj as md
        traj = self._load(traj_path, top)
        analysis = RGAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

    def hbonds(self, traj_path, top, output=None, **kwargs):
        traj = self._load(traj_path, top)
        analysis = HBondsAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

    def cluster(self, traj_path, top, output=None, **kwargs):
        import mdtraj as md
        traj = self._load(traj_path, top)
        analysis = ClusterAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

    def ss(self, traj_path, top, output=None, **kwargs):
        """
        Run secondary structure analysis (renamed to ss).
        """
        traj = self._load(traj_path, top)
        analysis = SSAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

    def sasa(self, traj_path, top, output=None, **kwargs):
        """
        Run solvent accessible surface area analysis.
        """
        traj = self._load(traj_path, top)
        analysis = SASAAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

    def dimred(self, traj_path, top, output=None, **kwargs):
        """
        Run dimensionality reduction analysis.
        """
        import mdtraj as md
        traj = self._load(traj_path, top)
        from .analysis.dimred import DimRedAnalysis
        analysis = DimRedAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

