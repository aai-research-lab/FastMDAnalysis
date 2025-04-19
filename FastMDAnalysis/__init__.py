"""
FastMDAnalysis Package Initialization

This version of FastMDAnalysis allows you to instantiate a single object by providing
the trajectory and topology file paths, along with optional parameters for frame and atom selection.
Frame selection is specified as an iterable with three elements (start, stop, stride); negative indices are supported—
e.g. using frames=(-10, -1, 1) or frames=[-10, -1, 1] will select frames from (n_frames - 10) to the last frame.
An atom selection string may also be provided to use a specific subset of atoms.
All subsequent analyses (rmsd, rmsf, rg, hbonds, cluster, ss, sasa, dimred) use the pre-loaded
trajectory and default atom selection unless overridden.
"""

import mdtraj as md
from .analysis import rmsd, rmsf, rg, hbonds, cluster, ss, dimred, sasa
from .utils import load_trajectory, create_dummy_trajectory

# Expose analysis classes.
RMSDAnalysis = rmsd.RMSDAnalysis
RMSFAnalysis = rmsf.RMSFAnalysis
RGAnalysis = rg.RGAnalysis
HBondsAnalysis = hbonds.HBondsAnalysis
ClusterAnalysis = cluster.ClusterAnalysis
SSAnalysis = ss.SSAnalysis
DimRedAnalysis = dimred.DimRedAnalysis
SASAAnalysis = sasa.SASAAnalysis

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
    frames : iterable of three int, optional
        An iterable containing three integers (start, stop, stride) to subset frames.
        Negative indices are allowed. For example, (-10, -1, 1) selects frames from (n_frames - 10)
        through the last frame.
        If None, the entire trajectory is used.
    atoms : str or None, optional
        An MDTraj atom selection string (e.g., "protein" or "protein and name CA") specifying which atoms to use.
        If None, all atoms are used.

    Examples
    --------
    >>> from FastMDAnalysis import FastMDAnalysis
    >>> # Load a trajectory, selecting every 10th frame from the start to the end, and using only protein atoms.
    >>> fastmda = FastMDAnalysis("trajectory.dcd", "topology.pdb", frames=(-10, -1, 1), atoms="protein")
    >>> rmsd_analysis = fastmda.rmsd(ref=0)
    """

    def __init__(self, traj_file, top_file, frames=None, atoms=None):
        # Load the full trajectory.
        #self.full_traj = md.load(traj_file, top=top_file)
        self.full_traj = load_trajectory(traj_file, top_file)


        # Subset frames if a frames iterable is provided.
        if frames is not None:
            try:
                # Convert the iterable to a tuple.
                frames = tuple(frames)
            except Exception as e:
                raise ValueError(f"Error converting frames to tuple: {e}")
            if len(frames) != 3:
                raise ValueError("The frames parameter must be an iterable of three integers: (start, stop, stride).")
            
            start, stop, stride = frames
            n_frames = self.full_traj.n_frames

            # Convert negative start index.
            if start < 0:
                start = n_frames + start
            # Convert negative stop index so that -1 includes the last frame.
            if stop < 0:
                stop = n_frames + stop + 1

            self.traj = self.full_traj[start:stop:stride]
        else:
            self.traj = self.full_traj

        # Store the default atom selection.
        self.default_atoms = atoms

    def _get_atoms(self, specific_atoms):
        """
        Determine the atom selection string to use.

        Parameters
        ----------
        specific_atoms : str or None
            The atom selection string specified for a particular analysis call.
            If None, the default atom selection provided during initialization is used.

        Returns
        -------
        str or None
            The atom selection string to be used.
        """
        return specific_atoms if specific_atoms is not None else self.default_atoms

    def rmsd(self, ref=0, atoms=None, **kwargs):
        """
        Run RMSD analysis on the stored trajectory.

        Parameters
        ----------
        ref : int, optional
            Reference frame index for RMSD calculations (default is 0).
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
        analysis = RMSDAnalysis(self.traj, ref_frame=ref, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def rmsf(self, atoms=None, **kwargs):
        """
        Run RMSF analysis on the stored trajectory.

        Parameters
        ----------
        atoms : str, optional
            Atom selection for RMSF analysis (if not provided, uses the default).
        kwargs : dict
            Additional keyword arguments for RMSFAnalysis.

        Returns
        -------
        RMSFAnalysis
            An RMSFAnalysis instance with the computed results.
        """
        a = self._get_atoms(atoms)
        analysis = RMSFAnalysis(self.traj, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def rg(self, atoms=None, **kwargs):
        """
        Run Radius of Gyration (RG) analysis.

        Parameters
        ----------
        atoms : str, optional
            Atom selection for RG analysis; if not provided, uses the default.
        kwargs : dict
            Additional keyword arguments for RGAnalysis.

        Returns
        -------
        RGAnalysis
            An RGAnalysis instance containing the results.
        """
        a = self._get_atoms(atoms)
        analysis = RGAnalysis(self.traj, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def hbonds(self, atoms=None, **kwargs):
        """
        Run Hydrogen Bonds (HBonds) analysis.

        Parameters
        ----------
        atoms : str, optional
            Atom selection for HBonds analysis; if not provided, uses the default.
        kwargs : dict
            Additional keyword arguments for HBondsAnalysis.

        Returns
        -------
        HBondsAnalysis
            An HBondsAnalysis instance with the computed results.
        """
        a = self._get_atoms(atoms)
        analysis = HBondsAnalysis(self.traj, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def cluster(self, methods="dbscan", eps=0.5, min_samples=5, n_clusters=None, atoms=None, **kwargs):
        """
        Run clustering analysis on the stored trajectory.

        Parameters
        ----------
        methods : str or list
            Clustering method(s) to use (e.g., "dbscan" or "kmeans"). If not a list, it is converted to one.
        eps : float, optional
            Maximum distance for DBSCAN clustering (default: 0.5).
        min_samples : int, optional
            Minimum number of samples for DBSCAN clustering (default: 5).
        n_clusters : int, optional
            Number of clusters for KMeans clustering (required if using KMeans).
        atoms : str, optional
            Atom selection for clustering. If not provided, the default is used.
        kwargs : dict
            Additional keyword arguments for ClusterAnalysis.

        Returns
        -------
        ClusterAnalysis
            A ClusterAnalysis instance containing clustering results.
        """
        a = self._get_atoms(atoms)
        from .analysis import cluster  # Import here for clarity.
        analysis = cluster.ClusterAnalysis(self.traj, methods=methods,
                                           eps=eps, min_samples=min_samples,
                                           n_clusters=n_clusters, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def ss(self, atoms=None, **kwargs):
        """
        Run Secondary Structure (SS) analysis on the stored trajectory.

        Parameters
        ----------
        atoms : str, optional
            Atom selection for SS analysis; if not provided, uses the default.
        kwargs : dict
            Additional keyword arguments for SSAnalysis.

        Returns
        -------
        SSAnalysis
            An SSAnalysis instance containing the computed results.
        """
        a = self._get_atoms(atoms)
        from .analysis import ss
        analysis = ss.SSAnalysis(self.traj, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def sasa(self, probe_radius=0.14, atoms=None, **kwargs):
        """
        Run Solvent Accessible Surface Area (SASA) analysis on the stored trajectory.

        Parameters
        ----------
        probe_radius : float, optional
            Probe radius (in nm) used for SASA calculation (default: 0.14).
        atoms : str, optional
            Atom selection for SASA analysis; if not provided, the default is used.
        kwargs : dict
            Additional keyword arguments for SASAAnalysis.

        Returns
        -------
        SASAAnalysis
            A SASAAnalysis instance with the computed results.
        """
        a = self._get_atoms(atoms)
        from .analysis import sasa
        analysis = sasa.SASAAnalysis(self.traj, probe_radius=probe_radius, atoms=a, **kwargs)
        analysis.run()
        return analysis

    def dimred(self, methods="all", atom_selection=None, **kwargs):
        """
        Run dimensionality reduction analysis on the stored trajectory.

        Parameters
        ----------
        methods : str or list
            Dimensionality reduction methods to use (options: "pca", "mds", "tsne").
            If "all" or if "all" is in the provided list, all available methods are applied.
        atom_selection : str, optional
            Atom selection for constructing the feature matrix. If not provided, uses the default atom selection.
        kwargs : dict
            Additional keyword arguments for DimRedAnalysis.

        Returns
        -------
        DimRedAnalysis
            A DimRedAnalysis instance containing the computed 2D embeddings.
        """
        if atom_selection is None:
            atom_selection = self.default_atoms
        from .analysis import dimred
        analysis = dimred.DimRedAnalysis(self.traj, methods=methods,
                                         atom_selection=atom_selection, **kwargs)
        analysis.run()
        return analysis

