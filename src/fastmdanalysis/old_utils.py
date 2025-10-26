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


    return traj

