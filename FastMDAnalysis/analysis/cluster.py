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

# Setup logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def adjust_labels(labels):
    """
    Adjust clustering labels so they are 1-based.
    
    Parameters
    ----------
    labels : array-like
        Original cluster labels.
        
    Returns
    -------
    np.ndarray
        1-based cluster labels.
    """
    labels = np.array(labels)
    if labels.min() < 1:
        shift = 1 - labels.min()
        logger.debug("Adjusting labels with shift: %d", shift)
        return labels + shift
    return labels

def get_cluster_cmap(n_clusters: int):
    """
    Return a categorical colormap for clustering.
    
    For n_clusters ≤ 12, uses a predefined set of 12 visually distinct colors
    (avoiding similar reds/blues and without gray). Otherwise, falls back to "nipy_spectral".
    
    Parameters
    ----------
    n_clusters : int
        Number of clusters.
        
    Returns
    -------
    ListedColormap or Colormap
        A colormap instance.
    """
    predefined_colors = [
        '#1f77b4',  # Blue
        '#ff7f0e',  # Orange
        '#2ca02c',  # Green
        '#d62728',  # Bright Red
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
    """
    Create a BoundaryNorm for discrete 1-based cluster labels.
    
    Parameters
    ----------
    unique_labels : array-like
        Sorted unique 1-based cluster labels.
        
    Returns
    -------
    BoundaryNorm
        A norm object.
    """
    boundaries = np.arange(unique_labels[0] - 0.5, unique_labels[-1] + 0.5 + 1, 1)
    logger.debug("Discrete boundaries: %s", boundaries)
    return BoundaryNorm(boundaries, len(boundaries) - 1)

def get_leaves(linkage_matrix, idx, N):
    """
    Recursively obtain the leaves (original observation indices) for a given index.
    
    Parameters
    ----------
    linkage_matrix : ndarray
        The linkage matrix of shape (N-1, 4).
    idx : int
        An index which may be an original observation (if idx < N) or an internal node (if idx ≥ N).
    N : int
        The number of original observations.
        
    Returns
    -------
    list of int
        List of original observation indices.
    """
    if idx < N:
        return [idx]
    # Valid internal nodes are from N to 2*N-2.
    if idx >= 2 * N - 1:
        logger.error("Index %d exceeds maximum allowed internal index %d", idx, 2 * N - 1)
        return []
    try:
        left = int(linkage_matrix[idx - N, 0])
        right = int(linkage_matrix[idx - N, 1])
        return get_leaves(linkage_matrix, left, N) + get_leaves(linkage_matrix, right, N)
    except IndexError as exc:
        logger.error("Index error in get_leaves: idx=%d, N=%d, linkage_matrix.shape=%s", idx, N, linkage_matrix.shape)
        return []

def dendrogram_link_color_func_factory(linkage_matrix, final_labels):
    """
    Create a link_color_func for the dendrogram.
    
    For a given branch, if all leaves share the same final (1-based) cluster label, return that color;
    otherwise, return gray.
    
    Parameters
    ----------
    linkage_matrix : ndarray
        The linkage matrix.
    final_labels : ndarray
        Final 1-based cluster labels.
        
    Returns
    -------
    function
        A function mapping a linkage index to a color.
    """
    N = len(final_labels)
    def link_color_func(i):
        leaves = get_leaves(linkage_matrix, i, N)
        if not leaves:
            logger.error("No leaves found for branch index %d", i)
            return "#808080"
        branch_labels = final_labels[leaves]
        if np.all(branch_labels == branch_labels[0]):
            unique = np.sort(np.unique(final_labels))
            cmap_local = get_cluster_cmap(len(unique))
            norm_local = get_discrete_norm(unique)
            # Convert RGBA tuple to hex string.
            color_val = cmap_local(norm_local(branch_labels[0]))
            color_hex = to_hex(color_val)
            logger.debug("Internal node %d: uniform cluster %d, color %s", i, branch_labels[0], color_hex)
            return color_hex
        else:
            logger.debug("Internal node %d: heterogeneous clusters %s", i, branch_labels)
            return "#808080"
    return link_color_func

class ClusterAnalysis(BaseAnalysis):
    def __init__(self, trajectory, methods='dbscan', eps: float = 0.5, min_samples: int = 5,
                 n_clusters: int = None, atoms: str = None, **kwargs):
        """
        Initialize clustering analysis.
        
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            Trajectory to analyze.
        methods : str or list
            Clustering methods to use ("dbscan", "kmeans", "hierarchical").
        eps : float, optional
            DBSCAN epsilon (default: 0.5).
        min_samples : int, optional
            DBSCAN minimum samples (default: 5).
        n_clusters : int, optional
            Number of clusters (required for kmeans and hierarchical).
        atoms : str, optional
            MDTraj atom selection string.
        kwargs : dict
            Additional arguments passed to BaseAnalysis.
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
        self.atom_indices = (self.traj.topology.select(self.atoms)
                             if self.atoms is not None else None)
        if self.atoms is not None and (self.atom_indices is None or len(self.atom_indices) == 0):
            raise AnalysisError(f"No atoms found with the selection: '{self.atoms}'")
        self.results = {}
        logger.info("ClusterAnalysis initialized with methods: %s", self.methods)

    def _calculate_rmsd_matrix(self) -> np.ndarray:
        """Calculate the pairwise RMSD matrix between frames."""
        logger.info("Calculating RMSD matrix...")
        n_frames = self.traj.n_frames
        distances = np.zeros((n_frames, n_frames))
        for i in range(n_frames):
            ref_frame = self.traj[i]
            if self.atom_indices is not None:
                distances[i] = md.rmsd(self.traj, ref_frame, atom_indices=self.atom_indices)
            else:
                distances[i] = md.rmsd(self.traj, ref_frame)
        logger.debug("RMSD matrix shape: %s", distances.shape)
        return (distances + distances.T) / 2.0

    def _plot_population(self, labels, filename, **kwargs):
        """Plot cluster populations as a bar plot with distinct colors."""
        logger.info("Plotting population bar plot...")
        unique = np.sort(np.unique(labels))
        counts = np.array([np.sum(labels == u) for u in unique])
        cmap = get_cluster_cmap(len(unique))
        norm = get_discrete_norm(unique)
        colors = [cmap(norm(u)) for u in unique]
        fig = plt.figure(figsize=(10, 6))
        plt.bar(unique, counts, width=0.8, color=colors)
        plt.title(kwargs.get("title", "Cluster Populations"))
        plt.xlabel(kwargs.get("xlabel", "Cluster ID"))
        plt.ylabel(kwargs.get("ylabel", "Number of Frames"))
        plt.xticks(unique)
        plt.grid(alpha=0.3)
        return self._save_plot(fig, filename)

    def _plot_cluster_trajectory_histogram(self, labels, filename, **kwargs):
        """Plot a histogram-style view of cluster assignments across frames."""
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
        """Plot a scatter view of cluster assignments for each frame."""
        logger.info("Plotting trajectory scatter...")
        frames = np.arange(len(labels))
        y_values = np.zeros_like(frames)
        unique = np.sort(np.unique(labels))
        cmap = get_cluster_cmap(len(unique))
        norm = get_discrete_norm(unique)
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.scatter(frames, y_values, c=labels, s=100, cmap=cmap, norm=norm, marker="o")
        ax.set_title(kwargs.get("title", "Cluster Trajectory Scatter Plot"))
        ax.set_xlabel(kwargs.get("xlabel", "Frame"))
        ax.set_yticks([])
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation="vertical", ticks=unique)
        cbar.ax.set_yticklabels([str(u) for u in unique])
        cbar.set_label("Cluster")
        return self._save_plot(fig, filename)

    def _plot_distance_matrix(self, distances, filename, **kwargs):
        """Plot the RMSD distance matrix as a heatmap with a colorbar."""
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
        """
        Plot a dendrogram for hierarchical clustering.
        
        Uses N = len(labels) as the number of original observations.
        Passes explicit labels (0, 1, ..., N-1) to dendrogram to ensure valid leaf order.
        Remaps these indices to the final 1-based cluster labels and colors the x-tick labels accordingly.
        """
        logger.info("Plotting dendrogram...")
        N = len(labels)
        logger.debug("N (number of observations) = %d", N)
        explicit_labels = np.arange(N)
        def color_func(i):
            # i should be a valid internal node index
            leaves = get_leaves(linkage_matrix, i, N)
            if not leaves:
                logger.error("No leaves found for internal node %d", i)
                return "#808080"
            branch_labels = labels[leaves]
            if np.all(branch_labels == branch_labels[0]):
                unique = np.sort(np.unique(labels))
                cmap_local = get_cluster_cmap(len(unique))
                norm_local = get_discrete_norm(unique)
                color_val = cmap_local(norm_local(branch_labels[0]))
                # Convert RGBA tuple to hex string
                color_hex = to_hex(color_val)
                logger.debug("Internal node %d: uniform cluster %d, color %s", i, branch_labels[0], color_hex)
                return color_hex
            else:
                logger.debug("Internal node %d: heterogeneous clusters %s", i, branch_labels)
                return "#808080"
        try:
            fig, ax = plt.subplots(figsize=(12, 6))
            dendro = dendrogram(linkage_matrix, ax=ax, labels=explicit_labels, link_color_func=color_func)
            leaf_order = dendro["leaves"]
            logger.debug("Dendrogram leaf order: %s", leaf_order)
            new_labels = []
            for i in leaf_order:
                if i < len(labels):
                    new_labels.append(str(labels[i]))
                else:
                    logger.error("Leaf index %d out of bounds (len(labels)=%d)", i, len(labels))
                    new_labels.append("NA")
            ax.set_xticklabels(new_labels, rotation=90)
            unique = np.sort(np.unique(labels))
            cmap_local = get_cluster_cmap(len(unique))
            norm_local = get_discrete_norm(unique)
            for tick, i in zip(ax.get_xticklabels(), leaf_order):
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
        """Save the matplotlib figure to a PNG file in the output directory."""
        plot_path = self.outdir / f"{name}.png"
        fig.savefig(plot_path, bbox_inches="tight")
        logger.info("Plot saved to %s", plot_path)
        return plot_path

    def run(self) -> dict:
        """
        Perform clustering analysis using the specified methods.
        
        Returns a dictionary with keys for each method ("dbscan", "kmeans", "hierarchical") and values containing:
          - "labels": 1-based cluster labels.
          - "pop_plot": Path to the population bar plot.
          - "trajectory_histogram": Path to the histogram-style trajectory plot.
          - "trajectory_scatter": Path to the scatter plot trajectory.
          - "distance_matrix_plot": (DBSCAN only) Path to the distance matrix heatmap.
          - "dendrogram_plot": (Hierarchical only) Path to the dendrogram plot.
        """
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
                if self.atom_indices is not None:
                    X = self.traj.xyz[:, self.atom_indices, :]
                else:
                    X = self.traj.xyz
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
                    linkage_matrix = linkage(X_flat, method='ward')
                    from scipy.cluster.hierarchy import fcluster
                    labels = fcluster(linkage_matrix, t=self.n_clusters, criterion='maxclust')
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

