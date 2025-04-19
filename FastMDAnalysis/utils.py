"""
Utility Functions for FastMDAnalysis

Provides functions to load trajectories and create dummy trajectories for testing.
This version of load_trajectory has been extended to accept multiple trajectory file inputs –
for example, a list, tuple, or a comma-separated string, as well as a glob pattern.
Any input acceptable by mdtraj.load is supported.
"""

import mdtraj as md
import glob
from pathlib import Path
import numpy as np

def load_trajectory(traj_input, top):
    """
    Load an MD trajectory using MDTraj.

    This function accepts various types of inputs for the trajectory:
      - A single file path (string or pathlib.Path).
      - A list or tuple of file paths.
      - A comma-separated string of file paths.
      - A glob pattern (string with wildcards such as "*" or "?" or "[").

    Parameters
    ----------
    traj_input : str, list, or tuple
        A trajectory file path, a list/tuple of file paths, a comma-separated string,
        or a glob pattern.
    top : str or pathlib.Path
        Path to the topology file.

    Returns
    -------
    mdtraj.Trajectory
        The loaded trajectory.

    Raises
    ------
    ValueError
        If no files are found when using a glob pattern.
    TypeError
        If traj_input is not one of the supported types.
    """
    # If the input is a list or tuple, use it directly.
    if isinstance(traj_input, (list, tuple)):
        files = [str(Path(f).resolve()) for f in traj_input]
        return md.load(files, top=str(Path(top).resolve()))
    # If the input is a string (or pathlib.Path)
    elif isinstance(traj_input, (str, Path)):
        traj_str = str(traj_input)
        # Check if the string contains a comma (i.e. multiple file paths are comma-separated)
        if ',' in traj_str:
            files = [s.strip() for s in traj_str.split(',')]
            return md.load(files, top=str(Path(top).resolve()))
        # Check if the string contains glob wildcards.
        elif any(char in traj_str for char in ['*', '?', '[']):
            files = sorted(glob.glob(traj_str))
            if not files:
                raise ValueError(f"No files found matching the glob pattern: {traj_str}")
            return md.load(files, top=str(Path(top).resolve()))
        else:
            return md.load(traj_str, top=str(Path(top).resolve()))
    else:
        raise TypeError("traj_input must be a string, list, or tuple")

def create_dummy_trajectory(n_frames: int = 5, n_atoms: int = 10) -> md.Trajectory:
    """
    Create a dummy MDTraj Trajectory for testing purposes.

    This function builds a simple topology with one chain (each residue with one CA atom)
    and creates random coordinates for the specified number of frames and atoms.

    Parameters
    ----------
    n_frames : int, optional
        Number of frames in the dummy trajectory (default: 5).
    n_atoms : int, optional
        Number of atoms (residues) per frame (default: 10).

    Returns
    -------
    md.Trajectory
        A dummy trajectory object for testing.
    """
    # Create random coordinates: shape (n_frames, n_atoms, 3)
    xyz = np.random.rand(n_frames, n_atoms, 3)
    
    # Create a simple topology: one chain with one atom per residue.
    from mdtraj.core.topology import Topology
    top = Topology()
    chain = top.add_chain()
    for i in range(n_atoms):
        residue = top.add_residue("GLY", chain)
        top.add_atom("CA", residue)
    
    traj = md.Trajectory(xyz, top)
    return traj

