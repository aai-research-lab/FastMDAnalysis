"""
Cluster Analysis Module

Performs clustering on an MD trajectory using a specified set of atoms.
By default, the atom selection is "protein and name CA" (unless overridden).
Supports three clustering methods:
  - dbscan: Uses a precomputed RMSD distance matrix.
  - kmeans: Uses the flattened coordinates of selected atoms.
  - hierarchical: Uses hierarchical clustering (Ward linkage) to generate a dendrogram and assign clusters.

For DBSCAN, parameters `eps` and `min_samples` are used.
For KMeans and hierarchical clustering, parameter `n_clusters` must be provided.

This module computes a pairwise RMSD distance matrix (if needed) using the selected atoms,
forces the matrix to be symmetric, and applies the specified clustering algorithm.
It then generates:
  - A bar plot of cluster populations with distinct colors.
  - Two trajectory projection plots:
      * A histogram-style plot where each frame is represented as a vertical bar colored by its cluster.
      * A scatter plot where each frame is plotted at y = 0 and colored by its cluster.
  - For DBSCAN, a heatmap plot of the RMSD distance matrix with an accompanying colorbar.
  - For hierarchical clustering, a dendrogram plot with branches and x-tick labels colored 
    according to the final clusters (branches that are not homogeneous are colored gray).

All computed data and plots are saved, and their file paths are stored in the results dictionary.
"""

import logging
from pathlib import Path
import numpy as np
import mdtraj as md
from sklearn.cluster import DBSCAN, KMeans

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, to_hex
from matplotlib.cm import ScalarMappable

from scipy.cluster.hierarchy import dendrogram, fcluster, linkage

from .base import BaseAnalysis, AnalysisError

# Use the module-specific logger (its configuration is set by the CLI)
logger = logging.getLogger(__name__)

def adjust_labels(labels):
    """Convert cluster labels to 1-based numbering."""
    labels = np.array(labels)
    if labels.min() < 1:
        shift = 1 - labels.min()
        logger.debug("Adjusting labels with shift: %d", shift)
        return labels + shift
    return labels

def get_cluster_cmap(n_clusters: int):
    """
    Return a categorical colormap for clustering.
    
    Uses a predefined set of 12 visually distinct colors for n_clusters ≤ 12; otherwise falls back to nipy_spectral.
    """
    predefined_colors = [
        '#1f77b4',  # Blue
        '#ff7f0e',  # Orange
        '#2ca02c',  # Green
        '#d62728',  # Red
        '#9467bd',  # Purple
        '#8c564b',  # Brown
        '#e377c2',  # Pink
        '#bcbd22',  # Olive
        '#17becf',  # Cyan
        '#e41a1c',  # Distinct Red
        '#377eb8',  # Different Blue
        '#f781bf'   # Magenta-ish
    ]
    if n_clusters <= len(predefined_colors):
        logger.debug("Using predefined colormap for %d clusters", n_clusters)
        return ListedColormap(predefined_colors[:n_clusters])
    else:
        logger.debug("Using fallback colormap for %d clusters", n_clusters)
        return plt.cm.get_cmap("nipy_spectral", n_clusters)

def get_discrete_norm(unique_labels):
    """Return a BoundaryNorm for the discrete cluster labels."""
    boundaries = np.arange(unique_labels[0] - 0.5, unique_labels[-1] + 0.5 + 1, 1)
    logger.debug("Discrete boundaries: %s", boundaries)
    return BoundaryNorm(boundaries, len(boundaries) - 1)

def get_leaves(linkage_matrix, idx, N):
    """Recursively retrieve the leaves (original frame indices) from the linkage matrix."""
    if idx < N:
        return [idx]
    if idx >= 2 * N - 1:
        logger.error("Index %d exceeds maximum allowed internal index %d", idx, 2 * N - 1)
        return []
    try:
        left = int(linkage_matrix[idx - N, 0])
        right = int(linkage_matrix[idx - N, 1])
        return get_leaves(linkage_matrix, left, N) + get_leaves(linkage_matrix, right, N)
    except IndexError:
        logger.error("Index error in get_leaves: idx=%d, N=%d, linkage_matrix.shape=%s", idx, N, linkage_matrix.shape)
        return []

def dendrogram_link_color_func_factory(linkage_matrix, final_labels):
    """
    Return a link_color_func that assigns a hex color if all leaves share the same cluster label,
    or gray if they are heterogeneous.
    """
    N = len(final_labels)
    def link_color_func(i):
        leaves = get_leaves(linkage_matrix, i, N)
        if not leaves:
            logger.error("No leaves found for internal node %d", i)
            return "#808080"
        branch_labels = final_labels[leaves]
        if np.all(branch_labels == branch_labels[0]):
            unique = np.sort(np.unique(final_labels))
            cmap_local = get_cluster_cmap(len(unique))
            norm_local = get_discrete_norm(unique)
            color_hex = to_hex(cmap_local(norm_local(branch_labels[0])))
            logger.debug("Internal node %d: uniform cluster %d, color %s", i, branch_labels[0], color_hex)
            return color_hex
        else:
            logger.debug("Internal node %d: heterogeneous clusters %s", i, branch_labels)
            return "#808080"
    return link_color_func

class ClusterAnalysis(BaseAnalysis):
    def __init__(self, trajectory, methods='dbscan', eps: float = 0.5,
                 min_samples: int = 5, n_clusters: int = None, atoms: str = None, **kwargs):
        """
        Initialize clustering analysis.
        
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            Trajectory to analyze.
        methods : str or list
            Clustering method(s) ("dbscan", "kmeans", "hierarchical").
        eps : float, optional
            DBSCAN epsilon (default: 0.5).
        min_samples : int, optional
            DBSCAN minimum samples (default: 5).
        n_clusters : int, optional
            Number of clusters (required for kmeans and hierarchical).
        atoms : str, optional
            MDTraj atom selection string.
        kwargs : dict
            Additional arguments.
        """
        super().__init__(trajectory, **kwargs)
        if isinstance(methods, str):
            self.methods = [methods.lower()]
        elif isinstance(methods, list):
            self.methods = [m.lower() for m in methods]
        else:
            raise AnalysisError("Parameter 'methods' must be a string or a list of strings.")
        self.eps = eps
        self.min_samples = min_samples
        self.n_clusters = n_clusters
        self.atoms = atoms
        self.atom_indices = self.traj.topology.select(self.atoms) if self.atoms is not None else None
        if self.atoms and (self.atom_indices is None or len(self.atom_indices) == 0):
            raise AnalysisError(f"No atoms found with the selection: '{self.atoms}'")
        self.results = {}

    def _calculate_rmsd_matrix(self) -> np.ndarray:
        """Calculate a symmetric pairwise RMSD matrix using the selected atoms."""
        logger.info("Calculating RMSD matrix...")
        n_frames = self.traj.n_frames
        distances = np.zeros((n_frames, n_frames))
        for i in range(n_frames):
            ref_frame = self.traj[i]
            distances[i] = md.rmsd(self.traj, ref_frame, atom_indices=self.atom_indices) if self.atom_indices is not None else md.rmsd(self.traj, ref_frame)
        logger.debug("RMSD matrix shape: %s", distances.shape)
        return (distances + distances.T) / 2.0

    def _plot_population(self, labels, filename, **kwargs):
        """Generate and save a population bar plot."""
        logger.info("Plotting population bar plot...")
        unique = np.sort(np.unique(labels))
        counts = np.array([np.sum(labels == u) for u in unique])
        cmap = get_cluster_cmap(len(unique))
        norm = get_discrete_norm(unique)
        fig = plt.figure(figsize=(10, 6))
        plt.bar(unique, counts, width=0.8, color=[cmap(norm(u)) for u in unique])
        plt.title(kwargs.get("title", "Cluster Populations"))
        plt.xlabel(kwargs.get("xlabel", "Cluster ID"))
        plt.ylabel(kwargs.get("ylabel", "Number of Frames"))
        plt.xticks(unique)
        plt.grid(alpha=0.3)
        return self._save_plot(fig, filename)

    def _plot_cluster_trajectory_histogram(self, labels, filename, **kwargs):
        """Generate and save a cluster trajectory histogram plot."""
        logger.info("Plotting trajectory histogram...")
        unique = np.sort(np.unique(labels))
        cmap = get_cluster_cmap(len(unique))
        norm = get_discrete_norm(unique)
        image_data = np.array(labels).reshape(1, -1)
        fig, ax = plt.subplots(figsize=(12, 4))
        im = ax.imshow(image_data, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
        ax.set_title(kwargs.get("title", "Cluster Trajectory Histogram"))
        ax.set_xlabel(kwargs.get("xlabel", "Frame"))
        ax.set_yticks([])
        cbar = fig.colorbar(im, ax=ax, orientation="vertical", ticks=unique)
        cbar.ax.set_yticklabels([str(u) for u in unique])
        cbar.set_label("Cluster")
        return self._save_plot(fig, filename)

    def _plot_cluster_trajectory_scatter(self, labels, filename, **kwargs):
        """Generate and save a cluster trajectory scatter plot."""
        logger.info("Plotting trajectory scatter...")
        frames = np.arange(len(labels))
        cmap = get_cluster_cmap(len(np.sort(np.unique(labels))))
        norm = get_discrete_norm(np.sort(np.unique(labels)))
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.scatter(frames, np.zeros_like(frames), c=labels, s=100, cmap=cmap, norm=norm, marker="o")
        ax.set_title(kwargs.get("title", "Cluster Trajectory Scatter Plot"))
        ax.set_xlabel(kwargs.get("xlabel", "Frame"))
        ax.set_yticks([])
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation="vertical", ticks=np.sort(np.unique(labels)))
        cbar.ax.set_yticklabels([str(u) for u in np.sort(np.unique(labels))])
        cbar.set_label("Cluster")
        return self._save_plot(fig, filename)

    def _plot_distance_matrix(self, distances, filename, **kwargs):
        """Generate and save an RMSD distance matrix heatmap."""
        logger.info("Plotting distance matrix heatmap...")
        fig = plt.figure(figsize=(10, 8))
        im = plt.imshow(distances, aspect="auto", interpolation="none", cmap=kwargs.get("cmap", "viridis"))
        plt.title(kwargs.get("title", "RMSD Distance Matrix"))
        plt.xlabel(kwargs.get("xlabel", "Frame"))
        plt.ylabel(kwargs.get("ylabel", "Frame"))
        cbar = plt.colorbar(im, ax=plt.gca())
        cbar.set_label("RMSD (nm)")
        return self._save_plot(fig, filename)

    def _plot_dendrogram(self, linkage_matrix, labels, filename, **kwargs):
        """Generate and save a dendrogram for hierarchical clustering."""
        logger.info("Plotting dendrogram...")
        N = len(labels)
        explicit_labels = np.arange(N)
        def color_func(i):
            leaves = get_leaves(linkage_matrix, i, N)
            if not leaves:
                logger.error("No leaves found for internal node %d", i)
                return "#808080"
            branch_labels = labels[leaves]
            if np.all(branch_labels == branch_labels[0]):
                unique = np.sort(np.unique(labels))
                cmap_local = get_cluster_cmap(len(unique))
                norm_local = get_discrete_norm(unique)
                color_hex = to_hex(cmap_local(norm_local(branch_labels[0])))
                logger.debug("Internal node %d: uniform cluster %d, color %s", i, branch_labels[0], color_hex)
                return color_hex
            else:
                logger.debug("Internal node %d: heterogeneous clusters %s", i, branch_labels)
                return "#808080"
        try:
            fig, ax = plt.subplots(figsize=(12, 6))
            dendro = dendrogram(linkage_matrix, ax=ax, labels=explicit_labels, link_color_func=color_func)
            new_labels = [str(labels[i]) if i < len(labels) else "NA" for i in dendro["leaves"]]
            ax.set_xticklabels(new_labels, rotation=90)
            unique = np.sort(np.unique(labels))
            cmap_local = get_cluster_cmap(len(unique))
            norm_local = get_discrete_norm(unique)
            for tick, i in zip(ax.get_xticklabels(), dendro["leaves"]):
                if i < len(labels):
                    tick.set_color(cmap_local(norm_local(labels[i])))
            ax.set_title(kwargs.get("title", "Hierarchical Clustering Dendrogram"))
            ax.set_xlabel(kwargs.get("xlabel", "Frame (Cluster Assignment)"))
            ax.set_ylabel(kwargs.get("ylabel", "Distance"))
            return self._save_plot(fig, filename)
        except Exception as e:
            logger.exception("Error during dendrogram plotting:")
            raise

    def _save_plot(self, fig, name: str):
        """Save the figure as a PNG file in the output directory and log its path."""
        plot_path = self.outdir / f"{name}.png"
        fig.savefig(plot_path, bbox_inches="tight")
        logger.info("Plot saved to %s", plot_path)
        return plot_path

    def run(self) -> dict:
        """
        Run the clustering analysis for the selected methods.
        """
        if self.results:
            logger.info("Results already computed; returning existing results.")
            return self.results

        try:
            logger.info("Starting clustering analysis...")
            results = {}
            distances = None
            if "dbscan" in self.methods:
                logger.info("Computing RMSD matrix for DBSCAN...")
                distances = self._calculate_rmsd_matrix()

            X_flat = None
            if any(method in self.methods for method in ["kmeans", "hierarchical"]):
                logger.info("Computing feature matrix for KMeans/Hierarchical...")
                X = self.traj.xyz[:, self.atom_indices, :] if self.atom_indices is not None else self.traj.xyz
                X_flat = X.reshape(self.traj.n_frames, -1)
                logger.debug("Feature matrix shape: %s", X_flat.shape)

            for method in self.methods:
                logger.info("Running method: %s", method)
                if method == "dbscan":
                    dbscan = DBSCAN(eps=self.eps, min_samples=self.min_samples, metric="precomputed")
                    labels = dbscan.fit_predict(distances)
                    labels = adjust_labels(labels)
                    logger.info("DBSCAN produced %d labels.", len(labels))
                    method_res = {"labels": labels, "distance_matrix": distances}
                    method_res["pop_plot"] = self._plot_population(labels, "dbscan_pop")
                    method_res["trajectory_histogram"] = self._plot_cluster_trajectory_histogram(labels, "dbscan_traj_hist")
                    method_res["trajectory_scatter"] = self._plot_cluster_trajectory_scatter(labels, "dbscan_traj_scatter")
                    method_res["distance_matrix_plot"] = self._plot_distance_matrix(distances, "dbscan_distance_matrix")
                    results["dbscan"] = method_res
                elif method == "kmeans":
                    if self.n_clusters is None:
                        raise AnalysisError("For KMeans clustering, n_clusters must be provided.")
                    kmeans = KMeans(n_clusters=self.n_clusters, random_state=42)
                    labels = kmeans.fit_predict(X_flat)
                    labels = adjust_labels(labels)
                    logger.info("KMeans produced %d labels.", len(labels))
                    method_res = {"labels": labels, "coordinates": X_flat}
                    method_res["pop_plot"] = self._plot_population(labels, "kmeans_pop")
                    method_res["trajectory_histogram"] = self._plot_cluster_trajectory_histogram(labels, "kmeans_traj_hist")
                    method_res["trajectory_scatter"] = self._plot_cluster_trajectory_scatter(labels, "kmeans_traj_scatter")
                    results["kmeans"] = method_res
                elif method == "hierarchical":
                    if self.n_clusters is None:
                        raise AnalysisError("For hierarchical clustering, n_clusters must be provided.")
                    logger.info("Computing linkage matrix for hierarchical clustering...")
                    linkage_matrix = linkage(X_flat, method="ward")
                    from scipy.cluster.hierarchy import fcluster
                    labels = fcluster(linkage_matrix, t=self.n_clusters, criterion="maxclust")
                    labels = adjust_labels(labels)
                    logger.info("Hierarchical clustering produced %d labels.", len(labels))
                    if len(labels) != self.traj.n_frames:
                        logger.warning("Mismatch: number of labels (%d) != number of frames (%d)", len(labels), self.traj.n_frames)
                    method_res = {"labels": labels, "linkage": linkage_matrix}
                    method_res["pop_plot"] = self._plot_population(labels, "hierarchical_pop")
                    method_res["trajectory_histogram"] = self._plot_cluster_trajectory_histogram(labels, "hierarchical_traj_hist")
                    method_res["trajectory_scatter"] = self._plot_cluster_trajectory_scatter(labels, "hierarchical_traj_scatter")
                    method_res["dendrogram_plot"] = self._plot_dendrogram(linkage_matrix, labels, "hierarchical_dendrogram")
                    results["hierarchical"] = method_res
                else:
                    raise AnalysisError(f"Unknown clustering method: {method}")

            self.results = results
            logger.info("Clustering analysis complete.")
            return results
        except Exception as e:
            logger.exception("Clustering failed:")
            raise AnalysisError(f"Clustering failed: {str(e)}")

    def plot(self, **kwargs):
        if not self.results:
            raise AnalysisError("No clustering results available. Run the analysis first.")
        return self.results

