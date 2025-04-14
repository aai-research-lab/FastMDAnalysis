#!/bin/bash
# Example CLI usage for FastMDAnalysis

# Ensure that the FastMDAnalysis CLI command ("fastmda") is available in your PATH.
# You may need to add the package's bin directory to your PATH.

###############################################
# RMSD Analysis:
# Computes the RMSD of each frame relative to the reference frame (default frame index: 0)
###############################################
echo "Running RMSD analysis..."
fastmda rmsd -traj path/to/trajectory.dcd -top path/to/topology.pdb -o rmsd_output --ref 0
echo "RMSD analysis completed. Output in rmsd_output/"

###############################################
# RMSF Analysis:
# Computes the per-atom RMSF. The example uses "c-alpha" atoms.
###############################################
echo "Running RMSF analysis..."
fastmda rmsf -traj path/to/trajectory.dcd -top path/to/topology.pdb -o rmsf_output --selection "c-alpha"
echo "RMSF analysis completed. Output in rmsf_output/"

###############################################
# Radius of Gyration (rg) Analysis:
###############################################
echo "Running Radius of Gyration analysis..."
fastmda rg -traj path/to/trajectory.dcd -top path/to/topology.pdb -o rg_output
echo "Rg analysis completed. Output in rg_output/"

###############################################
# Hydrogen Bonds (hbonds) Analysis:
###############################################
echo "Running Hydrogen Bonds analysis..."
fastmda hbonds -traj path/to/trajectory.dcd -top path/to/topology.pdb -o hbonds_output
echo "H-Bonds analysis completed. Output in hbonds_output/"

###############################################
# Cluster Analysis:
# Example with DBSCAN:
###############################################
echo "Running Cluster analysis with DBSCAN..."
fastmda cluster -traj path/to/trajectory.dcd -top path/to/topology.pdb -o cluster_output --methods dbscan --eps 0.5 --min_samples 5
echo "DBSCAN clustering completed. Check cluster_output/"

###############################################
# Cluster Analysis:
# Example with KMeans:
###############################################
echo "Running Cluster analysis with KMeans..."
fastmda cluster -traj path/to/trajectory.dcd -top path/to/topology.pdb -o cluster_output --methods kmeans --n_clusters 5
echo "KMeans clustering completed. Check cluster_output/"

###############################################
# Cluster Analysis:
# Example with both DBSCAN and KMeans:
###############################################
echo "Running Cluster analysis with both DBSCAN and KMeans..."
fastmda cluster -traj path/to/trajectory.dcd -top path/to/topology.pdb -o cluster_output --methods dbscan kmeans --n_clusters 5
echo "Combined clustering completed. Check cluster_output/"

###############################################
# Secondary Structure (ss) Analysis:
###############################################
echo "Running Secondary Structure (ss) analysis..."
fastmda ss -traj path/to/trajectory.dcd -top path/to/topology.pdb -o ss_output
echo "SS analysis completed. Output in ss_output/"

###############################################
# SASA Analysis:
# Computes Total SASA per frame, per-residue SASA, and average per-residue SASA.
###############################################
echo "Running SASA analysis..."
fastmda sasa -traj path/to/trajectory.dcd -top path/to/topology.pdb -o sasa_output --probe_radius 0.14
echo "SASA analysis completed. Output in sasa_output/"

###############################################
# Dimensionality Reduction (dimred) Analysis:
# Runs all available methods (PCA, MDS, t-SNE) using default atom selection.
###############################################
echo "Running Dimensionality Reduction analysis..."
fastmda dimred -traj path/to/trajectory.dcd -top path/to/topology.pdb -o dimred_output --methods all --atom_selection "protein and name CA"
echo "Dimensionality reduction analysis completed. Output in dimred_output/"

echo "All analyses have been executed."

