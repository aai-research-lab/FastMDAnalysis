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
from mdtraj.core.element import get_by_symbol  # Use get_by_symbol to obtain an Element instance

def load_trajectory(traj_input, top):
    """
    Load an MD trajectory using MDTraj.

    Accepts:
      - A single file path (string or pathlib.Path).
      - A list or tuple of file paths.
      - A comma-separated string of file paths.
      - A glob pattern (wildcards such as "*" or "?" or "[").
    
    Parameters
    ----------
    traj_input : str, list, or tuple
        A trajectory file path, a list/tuple of file paths, a comma-separated string, or a glob pattern.
    top : str or pathlib.Path
        Path to the topology file.

    Returns
    -------
    mdtraj.Trajectory
        The loaded trajectory.
    """
    if isinstance(traj_input, (list, tuple)):
        files = [str(Path(f).resolve()) for f in traj_input]
        return md.load(files, top=str(Path(top).resolve()))
    elif isinstance(traj_input, (str, Path)):
        traj_str = str(traj_input)
        if ',' in traj_str:
            files = [s.strip() for s in traj_str.split(',')]
            return md.load(files, top=str(Path(top).resolve()))
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

    Builds a topology with one chain; each residue contains one CA atom.
    Proper element information is provided via mdtraj.core.element.get_by_symbol.
    
    Parameters
    ----------
    n_frames : int, optional
        Number of frames in the dummy trajectory (default: 5).
    n_atoms : int, optional
        Number of atoms (residues) per frame (default: 10).

    Returns
    -------
    md.Trajectory
        A dummy trajectory object.
    """
    # Create random coordinates: shape (n_frames, n_atoms, 3)
    xyz = np.random.rand(n_frames, n_atoms, 3)
    
    # Create a simple topology.
    from mdtraj.core.topology import Topology
    top = Topology()
    chain = top.add_chain()
    element_C = get_by_symbol("C")  # Obtain element instance for Carbon.
    for i in range(n_atoms):
        residue = top.add_residue("GLY", chain)
        top.add_atom("CA", element_C, residue)
    
    traj = md.Trajectory(xyz, top)
    return traj

