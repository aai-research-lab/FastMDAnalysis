"""
Cluster Analysis Module

Performs clustering on an MD trajectory using selected atoms (default: protein C-alpha atoms).
Supports one or more clustering methods including 'dbscan' and 'kmeans'.

For DBSCAN, parameters 'eps' and 'min_samples' are used.
For KMeans, parameter 'n_clusters' must be provided.

The module computes a pairwise RMSD distance matrix (using the selected atoms),
forces matrix symmetry, applies the specified clustering algorithm(s),
saves the clustering results and distance (or coordinate) data,
and generates a bar plot showing cluster populations.
"""

from pathlib import Path
import numpy as np
import mdtraj as md
from sklearn.cluster import DBSCAN, KMeans

# Import plotting libraries.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from .base import BaseAnalysis, AnalysisError


class ClusterAnalysis(BaseAnalysis):
    def __init__(self, trajectory, methods='dbscan', eps: float = 0.5, min_samples: int = 5, n_clusters: int = None, **kwargs):
        """
        Initialize the clustering analysis.

        Args:
            trajectory: MD trajectory (an mdtraj.Trajectory object).
            methods (str or list): Clustering method(s) to use. Options: 'dbscan', 'kmeans'.
                                   Default is 'dbscan'. You may also pass a list to run more than one.
            eps (float): Maximum distance between two samples for DBSCAN.
            min_samples (int): Minimum number of samples required in a neighborhood for DBSCAN.
            n_clusters (int): Number of clusters for KMeans (required if 'kmeans' is used).
            kwargs: Additional base analysis arguments.
        """
        super().__init__(trajectory, **kwargs)
        if isinstance(methods, str):
            self.methods = [methods.lower()]
        elif isinstance(methods, list):
            self.methods = [m.lower() for m in methods]
        else:
            raise AnalysisError("Parameter 'methods' must be a string or list of strings.")
        self.eps = eps
        self.min_samples = min_samples
        self.n_clusters = n_clusters
        self.atom_indices = self._validate_indices()
        self.results = {}  # To store the results for each clustering method

    def _validate_indices(self) -> np.ndarray:
        """Validate and return atom indices for protein C-alpha atoms."""
        indices = self.traj.topology.select("protein and name CA")
        if indices is None or len(indices) == 0:
            raise AnalysisError("No valid C-alpha atoms found for clustering.")
        return indices

    def _calculate_rmsd_matrix(self) -> np.ndarray:
        """
        Calculate the pairwise RMSD matrix between frames based on the selected atoms.
        
        For each frame, the RMSD to all other frames is computed.
        The matrix is forced to be symmetric by averaging it with its transpose.
        """
        n_frames = self.traj.n_frames
        distances = np.zeros((n_frames, n_frames))
        for i in range(n_frames):
            ref_frame = self.traj[i]
            distances[i] = md.rmsd(self.traj, ref_frame, atom_indices=self.atom_indices)
        distances = (distances + distances.T) / 2.0
        return distances

    def run(self) -> dict:
        """
        Run the clustering analysis using the specified method(s).

        Returns:
            A dictionary containing clustering results for each method.
            For each method, results include cluster labels, associated data, and the population plot file.
        """
        distances = self._calculate_rmsd_matrix()
        if distances.shape != (self.traj.n_frames, self.traj.n_frames):
            raise AnalysisError("Invalid distance matrix shape.")
        if np.any(np.isnan(distances)):
            raise AnalysisError("Distance matrix contains NaN values.")

        # Run clustering for each requested method.
        for method in self.methods:
            if method == 'dbscan':
                dbscan = DBSCAN(eps=self.eps, min_samples=self.min_samples, metric='precomputed')
                labels = dbscan.fit_predict(distances)
                self.results['dbscan'] = {
                    'labels': labels,
                    'distance_matrix': distances
                }
                self._save_data(labels.reshape(-1, 1), "cluster_labels_dbscan")
                self._save_data(distances, "cluster_distance_matrix_dbscan")
                plot_path = self._plot_population(labels, "dbscan_clusters_population")
                self.results['dbscan']['pop_plot'] = plot_path
            elif method == 'kmeans':
                if self.n_clusters is None:
                    raise AnalysisError("n_clusters must be provided for KMeans clustering.")
                # For KMeans, use flattened coordinates of the selected atoms.
                X = self.traj.xyz[:, self.atom_indices, :]  # shape: (n_frames, n_sel, 3)
                X_flat = X.reshape(X.shape[0], -1)
                kmeans = KMeans(n_clusters=self.n_clusters, random_state=42)
                labels = kmeans.fit_predict(X_flat)
                self.results['kmeans'] = {
                    'labels': labels,
                    'coordinates': X_flat
                }
                self._save_data(labels.reshape(-1, 1), "cluster_labels_kmeans")
                self._save_data(X_flat, "cluster_coordinates_kmeans")
                plot_path = self._plot_population(labels, "kmeans_clusters_population")
                self.results['kmeans']['pop_plot'] = plot_path
            else:
                raise AnalysisError(f"Unknown clustering method: {method}")
        return self.results

    def _plot_population(self, labels, filename, **kwargs):
        """
        Generate a bar plot showing cluster populations from the given labels.
        """
        unique, counts = np.unique(labels, return_counts=True)
        title = kwargs.get('title', "Cluster Populations")
        xlabel = kwargs.get('xlabel', "Cluster ID")
        ylabel = kwargs.get('ylabel', "Number of Frames")
        color = kwargs.get('color', None)
        fig = plt.figure(figsize=(10, 6))
        plt.bar(unique, counts, width=0.8, color=color)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xticks(unique)
        plt.grid(alpha=0.3)
        plot_path = self._save_plot(fig, filename)
        plt.close(fig)
        return plot_path

    def plot(self, method=None, **kwargs):
        """
        Generic plot method for clustering analysis.
        This method overrides the BaseAnalysis plot() method.
        
        If 'method' is not specified, the default is the first clustering method used.
        It returns the cluster population bar plot for that method.
        
        Args:
            method (str): The clustering method for which to generate the plot (e.g. 'dbscan' or 'kmeans').
            kwargs: Additional keyword arguments to customize the population plot.
        
        Returns:
            Path: The file path to the saved population plot.
        """
        if method is None:
            method = self.methods[0]
        if method not in self.results:
            raise AnalysisError(f"No clustering result for method: {method}")
        labels = self.results[method]['labels']
        return self._plot_population(labels, f"{method}_clusters_population", **kwargs)

