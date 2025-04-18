#!/bin/bash
# Example CLI usage for FastMDAnalysis
# This script demonstrates how to run various MD analyses using the fastmda CLI command.
# Global options --frames and --atoms allow you to specify a frame selection and default atom selection for all analyses.
# Make sure that the fastmda command is available in your PATH (via module load or proper installation).

# Example frame selection: "0,-1,10" selects frames from the first frame to the last frame with a stride of 10.
# Example atom selection: "protein" uses all protein atoms.

###############################################
# RMSD Analysis:
# Computes RMSD vs. frame using a reference frame (default index 0).
###############################################
echo "Running RMSD analysis..."
fastmda rmsd -traj path/to/trajectory.dcd -top path/to/topology.pdb -o rmsd_output --frames 0,-1,10 --atoms "protein" --ref 0
echo "RMSD analysis completed. Output in rmsd_output/"

###############################################
# RMSF Analysis:
# Computes per-atom RMSF; uses the default atom selection or can override with --selection.
###############################################
echo "Running RMSF analysis..."
fastmda rmsf -traj path/to/trajectory.dcd -top path/to/topology.pdb -o rmsf_output --frames 0,-1,10 --atoms "protein"
echo "RMSF analysis completed. Output in rmsf_output/"

###############################################
# Radius of Gyration (RG) Analysis:
###############################################
echo "Running Radius of Gyration analysis..."
fastmda rg -traj path/to/trajectory.dcd -top path/to/topology.pdb -o rg_output --frames 0,-1,10 --atoms "protein"
echo "RG analysis completed. Output in rg_output/"

###############################################
# Hydrogen Bonds (HBonds) Analysis:
###############################################
echo "Running Hydrogen Bonds analysis..."
fastmda hbonds -traj path/to/trajectory.dcd -top path/to/topology.pdb -o hbonds_output --frames 0,-1,10
echo "HBonds analysis completed. Output in hbonds_output/"

###############################################
# Cluster Analysis (DBSCAN):
###############################################
echo "Running Cluster analysis with DBSCAN..."
fastmda cluster -traj path/to/trajectory.dcd -top path/to/topology.pdb -o cluster_output --frames 0,-1,10 --atoms "protein" --methods dbscan --eps 0.5 --min_samples 5
echo "DBSCAN clustering completed. Check cluster_output/"

###############################################
# Cluster Analysis (KMeans):
###############################################
echo "Running Cluster analysis with KMeans..."
fastmda cluster -traj path/to/trajectory.dcd -top path/to/topology.pdb -o cluster_output --frames 0,-1,10 --atoms "protein" --methods kmeans --n_clusters 5
echo "KMeans clustering completed. Check cluster_output/"

###############################################
# Secondary Structure (SS) Analysis:
###############################################
echo "Running Secondary Structure analysis..."
fastmda ss -traj path/to/trajectory.dcd -top path/to/topology.pdb -o ss_output --frames 0,-1,10 --atoms "protein"
echo "SS analysis completed. Output in ss_output/"

###############################################
# SASA Analysis:
# Computes Total SASA vs. frame, per-residue SASA heatmap, and average per-residue SASA.
###############################################
echo "Running SASA analysis..."
fastmda sasa -traj path/to/trajectory.dcd -top path/to/topology.pdb -o sasa_output --frames 0,-1,10 --atoms "protein" --probe_radius 0.14
echo "SASA analysis completed. Output in sasa_output/"

###############################################
# Dimensionality Reduction Analysis:
# Runs specified methods (e.g., PCA and t-SNE) using an atom selection for the feature matrix.
###############################################
echo "Running Dimensionality Reduction analysis..."
fastmda dimred -traj path/to/trajectory.dcd -top path/to/topology.pdb -o dimred_output --frames 0,-1,10 --atoms "protein" --methods pca tsne --atom_selection "protein and name CA"
echo "Dimensionality Reduction analysis completed. Output in dimred_output/"

echo "All analyses have been executed."

