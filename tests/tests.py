"""
Unit tests for FastMDAnalysis using real dataset files.

This module tests the primary analyses (rmsd, rmsf, rg, hbonds, cluster, ss, sasa, dimred)
using a real dataset (e.g., ubiquitin) provided by the datasets module.
"""

import unittest
import numpy as np
import mdtraj as md
from pathlib import Path

from FastMDAnalysis import FastMDAnalysis
from FastMDAnalysis.datasets import trp_cage
from FastMDAnalysis.utils import load_trajectory

class TestFastMDAnalysis(unittest.TestCase):
    def setUp(self):
        # Retrieve the ubiquitin dataset file paths.
        traj_path, top_path = trp_cage.traj, trp_cage.top
        # Load the trajectory using the load_trajectory utility function.
        self.traj = load_trajectory(traj_path, top_path)
        # Create standard bonds in case they are missing (needed for hbonds analysis).
        self.traj.topology.create_standard_bonds()
        # Initialize FastMDAnalysis with the ubiquitin dataset.
        self.fastmda = FastMDAnalysis(traj_path, top_path, frames=None, atoms="protein")

    def test_rmsd(self):
        analysis = self.fastmda.rmsd(ref=0)
        self.assertTrue(hasattr(analysis, "data"), "RMSD analysis missing data attribute.")
        self.assertIsInstance(analysis.data, np.ndarray)

    def test_rmsf(self):
        analysis = self.fastmda.rmsf()
        self.assertTrue(hasattr(analysis, "data"), "RMSF analysis missing data attribute.")
        self.assertIsInstance(analysis.data, np.ndarray)

    def test_rg(self):
        analysis = self.fastmda.rg()
        self.assertTrue(hasattr(analysis, "data"), "Radius of gyration analysis missing data attribute.")
        self.assertIsInstance(analysis.data, np.ndarray)

    def test_hbonds(self):
        analysis = self.fastmda.hbonds()
        self.assertTrue(hasattr(analysis, "data"), "Hydrogen bonds analysis missing data attribute.")

    def test_cluster_dbscan(self):
        analysis = self.fastmda.cluster(methods="dbscan")
        self.assertIn("dbscan", analysis.results, "DBSCAN results missing in cluster analysis.")
        self.assertIn("labels", analysis.results["dbscan"], "DBSCAN labels missing in cluster results.")
        self.assertIsInstance(analysis.results["dbscan"]["labels"], np.ndarray)

    def test_cluster_kmeans(self):
        analysis = self.fastmda.cluster(methods="kmeans", n_clusters=3)
        self.assertIn("kmeans", analysis.results, "KMeans results missing in cluster analysis.")
        self.assertIn("labels", analysis.results["kmeans"], "KMeans labels missing in cluster results.")

    def test_cluster_hierarchical(self):
        analysis = self.fastmda.cluster(methods="hierarchical", n_clusters=3)
        self.assertIn("hierarchical", analysis.results, "Hierarchical results missing in cluster analysis.")
        self.assertIn("labels", analysis.results["hierarchical"], "Hierarchical labels missing in cluster results.")

    def test_ss(self):
        analysis = self.fastmda.ss()
        self.assertTrue(hasattr(analysis, "data"), "Secondary structure analysis missing data attribute.")

    def test_sasa(self):
        analysis = self.fastmda.sasa(probe_radius=0.14)
        self.assertTrue(hasattr(analysis, "data"), "SASA analysis missing data attribute.")

    def test_dimred(self):
        analysis = self.fastmda.dimred(methods=["pca", "mds", "tsne"])
        self.assertTrue(hasattr(analysis, "data"), "Dimensionality reduction analysis missing data attribute.")

    def test_plotting(self):
        # Test that the plot() method returns output file paths.
        analysis = self.fastmda.rmsd(ref=0)
        plot_result = analysis.plot()
        self.assertTrue(isinstance(plot_result, dict) or isinstance(plot_result, Path),
                        "Plot method must return a dict or a Path.")
        if isinstance(plot_result, dict):
            self.assertGreater(len(plot_result), 0, "Plot method returned an empty dict.")
            for key, file_path in plot_result.items():
                self.assertTrue(Path(file_path).exists(), f"Plot for {key} not found at {file_path}")
        else:
            self.assertTrue(Path(plot_result).exists(), "Plot file does not exist.")

if __name__ == "__main__":
    unittest.main()

