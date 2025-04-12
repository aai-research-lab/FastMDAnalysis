#!/bin/bash
# Example command-line usage for FastMDAnalysis.

# RMSD analysis example:
fastmda rmsd -traj path/to/trajectory.dcd -top path/to/topology.pdb -o rmsd_output --ref 0

# Cluster analysis example with additional arguments:
fastmda cluster -traj path/to/trajectory.dcd -top path/to/topology.pdb -o clusters --eps 0.4 --min_samples 10

