"""
Utility Functions for FastMDAnalysis

Provides functions to load trajectories and create dummy trajectories for testing.
"""

import mdtraj as md
from pathlib import Path
import numpy as np

def load_trajectory(traj_path: str, top: str) -> md.Trajectory:
    """
    Load an MD trajectory using MDTraj from the provided file paths.

    Parameters
    ----------
    traj_path : str
        Path to the trajectory file (e.g., DCD, XTC).
    top : str
        Path to the topology file (e.g., PDB).

    Returns
    -------
    md.Trajectory
        The loaded trajectory.

    Raises
    ------
    Exception
        If the trajectory cannot be loaded.
    """
    try:
        traj = md.load(traj_path, top=top)
        return traj
    except Exception as e:
        raise Exception(f"Error loading trajectory: {e}")

def create_dummy_trajectory(n_frames: int = 5, n_atoms: int = 10) -> md.Trajectory:
    """
    Create a dummy MDTraj Trajectory for testing purposes.

    This function builds a topology with a single chain where each residue contains one CA atom,
    and creates random coordinates for the specified number of frames and atoms.

    Parameters
    ----------
    n_frames : int, optional
        Number of frames in the dummy trajectory. Default is 5.
    n_atoms : int, optional
        Number of atoms (residues) in each frame. Default is 10.

    Returns
    -------
    md.Trajectory
        A dummy trajectory object.
    """
    from mdtraj.core.topology import Topology

    # Create dummy coordinates: shape (n_frames, n_atoms, 3)
    xyz = np.random.rand(n_frames, n_atoms, 3)
    
    # Create a simple topology with one chain; each residue has one CA atom.
    top = Topology()
    chain = top.add_chain()
    for i in range(n_atoms):
        residue = top.add_residue("GLY", chain)
        top.add_atom("CA", residue)
    
    traj = md.Trajectory(xyz, top)
    return traj

