"""
FastMDAnalysis Package Initialization

Exports key analysis modules for API usage.
Users can import analysis classes with a single import.
Note: The CLI module is not imported here to avoid circular dependencies.
"""

from .analysis import (
    rmsd,
    rmsf,
    rg,
    hbonds,
    cluster,
    secondary_structure,
    dimred  # This now includes our new module
)
from .utils import load_trajectory, create_dummy_trajectory

# Expose analysis classes for API usage.
RMSDAnalysis = rmsd.RMSDAnalysis
RMSFAnalysis = rmsf.RMSFAnalysis
RGAnalysis = rg.RGAnalysis
HBondsAnalysis = hbonds.HBondsAnalysis
ClusterAnalysis = cluster.ClusterAnalysis
SecondaryStructureAnalysis = secondary_structure.SecondaryStructureAnalysis

class FastMDAnalysis:
    """
    Main API class for MD trajectory analysis.
    Instantiates each analysis type as a method.
    """

    def __init__(self):
        # Initialize common settings if needed
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
        import mdtraj as md
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
        import mdtraj as md
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

    def secondary_structure(self, traj_path, top, output=None, **kwargs):
        import mdtraj as md
        traj = self._load(traj_path, top)
        analysis = SecondaryStructureAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

    def dimred(self, traj_path, top, output=None, **kwargs):
        """
        Run dimensionality reduction analysis.
        """
        import mdtraj as md
        traj = self._load(traj_path, top)
        analysis = DimRedAnalysis(traj, output=output, **kwargs)
        analysis.run()
        return analysis

