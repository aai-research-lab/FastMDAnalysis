"""
Utility functions for FastMDAnalysis.
"""

import mdtraj as md
from pathlib import Path
import numpy as np

def load_trajectory(traj_path: str, top: str):
    """
    Load a trajectory using MDTraj.
    
    Args:
        traj_path (str): Path to the trajectory file.
        top (str): Path to the topology file.
    
    Returns:
        md.Trajectory: The loaded trajectory.
    """
    try:
        traj = md.load(traj_path, top=top)
        return traj
    except Exception as e:
        raise Exception(f"Error loading trajectory: {e}")


def create_dummy_trajectory(n_frames: int = 5, n_atoms: int = 10):
    """
    Create a dummy MDTraj Trajectory for testing purposes.
    The topology is built with a single chain with each residue containing one CA atom.

    Args:
        n_frames (int): Number of frames.
        n_atoms (int): Number of atoms.
    
    Returns:
        md.Trajectory: A dummy trajectory object.
    """
    import mdtraj as md
    import numpy as np
    from mdtraj.core.topology import Topology

    # Create dummy coordinates: shape (n_frames, n_atoms, 3)
    xyz = np.random.rand(n_frames, n_atoms, 3)
    
    # Create a simple topology.
    top = Topology()
    chain = top.add_chain()
    for i in range(n_atoms):
        residue = top.add_residue("GLY", chain)
        top.add_atom("CA", residue)
    
    traj = md.Trajectory(xyz, top)
    return traj

