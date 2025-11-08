# FastMDAnalysis/src/fastmdanalysis/analysis/cluster.py
"""
Cluster Analysis Module

Performs clustering on an MD trajectory using a specified set of atoms.
By default, all methods are run:
  - dbscan: Uses a precomputed RMSD distance matrix (nm).
  - kmeans: Uses flattened coordinates (optionally PCA upstream in future).
  - hierarchical: Ward linkage over flattened coordinates.

Key behaviors:
- DBSCAN labels are relabeled to compact positive integers 1..K; if noise exists (-1 in raw),
  it is mapped to K+1 so external consumers always see 1,2,3,... (no negatives).
- Raw sklearn labels (-1 for noise) are also saved for provenance.
- Extensive logging aids debugging and parameter tuning.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional, Tuple, Sequence, Union

import numpy as np
import mdtraj as md

from sklearn.cluster import DBSCAN, KMeans

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, to_hex
from matplotlib.cm import ScalarMappable

from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import squareform

from .base import BaseAnalysis, AnalysisError

# Module logger (configured by CLI / caller)
logger = logging.getLogger(__name__)

# ----------------------------- Label helpers ----------------------------------

def relabel_compact_positive(
    labels_raw: np.ndarray,
    start: int = 1,
    noise_as_last: bool = True
) -> Tuple[np.ndarray, Dict[int, int], Optional[int]]:
    """
    Map labels to contiguous positive integers 1..K; optionally map noise (-1) to K+1.

    Parameters
    ----------
    labels_raw : array-like of int
        sklearn-style labels (DBSCAN uses -1 for noise).
    start : int
        First compact label (default 1).
    noise_as_last : bool
        If True and noise present, it is mapped to K+1.

    Returns
    -------
    labels_compact : np.ndarray
        Positive labels 1..K (and K+1 for noise if present).
    mapping : dict
        Mapping from original non-negative labels to compact labels.
    noise_label : Optional[int]
        The compact label used for noise (None if no noise).
    """
    raw = np.asarray(labels_raw, dtype=int)
    uniq_nonneg = sorted([u for u in np.unique(raw) if u >= 0])
    mapping = {u: (i + start) for i, u in enumerate(uniq_nonneg)}

    labels = np.empty_like(raw)
    for u, m in mapping.items():
        labels[raw == u] = m

    noise_label = None
    if np.any(raw == -1):
        noise_label = start + len(uniq_nonneg)
        labels[raw == -1] = noise_label
    return labels, mapping, noise_label

# ----------------------------- Colormaps/Norms --------------------------------

def get_cluster_cmap(n_clusters: int):
    predefined_colors = ['#e41a1c','#377eb8','#00ff00','#ffd700',
               '#9932cc','#ffa500','#00bfff','#a52a2a',
               '#808080','#000000','#006400','#000080'
               ]
    if n_clusters <= len(predefined_colors):
        logger.debug("Using predefined colormap for %d clusters", n_clusters)
        return ListedColormap(predefined_colors[:n_clusters])
    logger.debug("Using fallback colormap for %d clusters", n_clusters)
    return plt.cm.get_cmap("nipy_spectral", n_clusters)

def get_discrete_norm(unique_labels):
    unique_labels = np.asarray(unique_labels, dtype=int)
    unique_labels.sort()
    boundaries = np.arange(unique_labels[0] - 0.5, unique_labels[-1] + 1.5, 1)
    logger.debug("Discrete boundaries: %s", boundaries)
    return BoundaryNorm(boundaries, len(boundaries) - 1)

# ----------------------------- Dendrogram helpers -----------------------------

def get_leaves(linkage_matrix, idx, N):
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
        logger.debug("Internal node %d: heterogeneous clusters %s", i, branch_labels)
        return "#808080"
    return link_color_func

# ------------------------------- Core class -----------------------------------

class ClusterAnalysis(BaseAnalysis):
    def __init__(
        self,
        trajectory,
        methods: Union[str, Sequence[str]] = "all",
        eps: float = 0.2,
        min_samples: int = 5,
        n_clusters: Optional[int] = None,
        atoms: Optional[str] = None,
        plot_style: str = "full",
        combined_plot_name: str = "cluster_combined",
        save_data: bool = True,
        feature_mode: str = "coordinates",
        **kwargs
    ):
        """
        Parameters
        ----------
        methods : {"all"} | str | list[str]
            Which clustering methods to run. Default "all" expands to
            ["dbscan", "kmeans", "hierarchical"].
        eps : float
            DBSCAN epsilon in **nm** (MDTraj RMSD units). 0.2 nm ≈ 2 Å is a good starting point.
        min_samples : int
            DBSCAN minimum samples (default: 5).
        n_clusters : int, optional
            If None and methods include kmeans/hierarchical, defaults to 3.
        atoms : str, optional
            MDTraj atom selection string used to slice the trajectory for all methods.
        plot_style : {"full", "minimal", "none"}
            Controls which plots are generated. "full" retains the legacy behavior,
            "minimal" emits a single combined scatter summary, and "none" skips plots.
        combined_plot_name : str
            Base name for the combined scatter figure when plot_style="minimal".
        feature_mode : {"coordinates", "distance"}
            Select the feature space for kmeans/hierarchical clustering. "coordinates" flattens
            Cartesian coordinates (legacy behaviour), while "distance" reuses the RMSD distance
            matrix to align with lightweight benchmarking scenarios.
        """
        super().__init__(trajectory, save_data=save_data, **kwargs)

        # Normalize methods
        if isinstance(methods, str):
            methods_norm = [methods.lower()]
        else:
            methods_norm = [m.lower() for m in methods]

        if "all" in methods_norm:
            self.methods = ["dbscan", "kmeans", "hierarchical"]
        else:
            self.methods = methods_norm

        self.eps = float(eps)
        self.min_samples = int(min_samples)
        self.n_clusters = int(n_clusters) if (n_clusters is not None and int(n_clusters) > 0) else None
        self.atoms = atoms
        if plot_style not in {"full", "minimal", "none"}:
            raise AnalysisError("plot_style must be 'full', 'minimal', or 'none'.")
        self.plot_style = plot_style
        self.combined_plot_name = combined_plot_name
        if feature_mode not in {"coordinates", "distance"}:
            raise AnalysisError("feature_mode must be 'coordinates' or 'distance'.")
        self.feature_mode = feature_mode

        self.atom_indices = self.traj.topology.select(self.atoms) if self.atoms is not None else None
        if self.atoms and (self.atom_indices is None or len(self.atom_indices) == 0):
            raise AnalysisError(f"No atoms found with the selection: '{self.atoms}'")

        self.results: Dict[str, Dict] = {}
        self._labels_for_combined: Dict[str, np.ndarray] = {}

    # ----------------------------- Distances/Features --------------------------

    def _calculate_rmsd_matrix(self) -> np.ndarray:
        """Compute a symmetric pairwise RMSD matrix (nm) over frames with low memory overhead."""
        logger.info("Calculating RMSD matrix via mdtraj._rmsd (float32)...")

        from mdtraj import _rmsd

        coords = self.traj.xyz.astype(np.float32, copy=True, order="C")
        if self.atom_indices is not None:
            coords = np.ascontiguousarray(coords[:, self.atom_indices, :], dtype=np.float32)

        T, n_atoms, _ = coords.shape
        traces = _rmsd._center_inplace_atom_major(coords)

        result = np.empty((T, T), dtype=np.float32)

        for ref_idx in range(T):
            column = _rmsd.getMultipleRMSDs_atom_major(
                coords,
                coords,
                traces,
                traces,
                ref_idx,
            )
            result[:, ref_idx] = column.astype(np.float32, copy=False)

        # Ensure exact symmetry and zero diagonal
        result = 0.5 * (result + result.T)
        np.fill_diagonal(result, 0.0)

        logger.debug(
            "RMSD matrix computed: shape=%s min=%.4f max=%.4f nm",
            result.shape,
            float(result.min()),
            float(result.max()),
        )
        return result

    def _distance_diagnostics(self, D: np.ndarray) -> Dict[str, float]:
        tri = D[np.triu_indices_from(D, k=1)]
        if tri.size == 0:
            return {"p25": 0.0, "p50": 0.0, "p75": 0.0, "p90": 0.0, "max": 0.0}
        diags = {
            "p25": float(np.percentile(tri, 25)),
            "p50": float(np.percentile(tri, 50)),
            "p75": float(np.percentile(tri, 75)),
            "p90": float(np.percentile(tri, 90)),
            "max": float(np.max(tri)),
        }
        logger.info("RMSD percentiles (nm): p25=%.3f p50=%.3f p75=%.3f p90=%.3f max=%.3f",
                    diags["p25"], diags["p50"], diags["p75"], diags["p90"], diags["max"])
        return diags

    # ----------------------------- Plotting helpers ----------------------------

    def _plot_population(self, labels, filename, **kwargs):
        logger.info("Plotting population bar plot...")
        unique = np.sort(np.unique(labels))
        counts = np.array([np.sum(labels == u) for u in unique])
        fig = plt.figure(figsize=(10, 6))
        cmap = get_cluster_cmap(len(unique))
        norm = get_discrete_norm(unique)
        plt.bar(unique, counts, width=0.8, color=[cmap(norm(u)) for u in unique])
        plt.title(kwargs.get("title", "Cluster Populations"))
        plt.xlabel(kwargs.get("xlabel", "Cluster ID (compact)"))
        plt.ylabel(kwargs.get("ylabel", "Number of Frames"))
        plt.xticks(unique)
        plt.grid(alpha=0.3)
        return self._save_plot(fig, filename)

    def _plot_cluster_trajectory_histogram(self, labels, filename, **kwargs):
        logger.info("Plotting trajectory histogram...")
        unique = np.sort(np.unique(labels))
        image_data = np.array(labels).reshape(1, -1)
        cmap = get_cluster_cmap(len(unique))
        norm = get_discrete_norm(unique)
        fig, ax = plt.subplots(figsize=(12, 4))
        im = ax.imshow(image_data, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
        ax.set_title(kwargs.get("title", "Cluster Trajectory Histogram"))
        ax.set_xlabel(kwargs.get("xlabel", "Frame"))
        ax.set_yticks([])
        cbar = fig.colorbar(im, ax=ax, orientation="vertical", ticks=unique)
        cbar.ax.set_yticklabels([str(u) for u in unique])
        cbar.set_label("Cluster (compact)")
        return self._save_plot(fig, filename)

    def _plot_cluster_trajectory_scatter(self, labels, filename, **kwargs):
        logger.info("Plotting trajectory scatter...")
        frames = np.arange(len(labels))
        fig, ax = plt.subplots(figsize=(10, 4))
        unique = np.sort(np.unique(labels))
        cmap = get_cluster_cmap(len(unique))
        norm = get_discrete_norm(unique)
        ax.scatter(frames, np.zeros_like(frames), c=labels, s=100, cmap=cmap, norm=norm, marker="o")
        ax.set_title(kwargs.get("title", "Cluster Trajectory Scatter Plot"))
        ax.set_xlabel(kwargs.get("xlabel", "Frame"))
        ax.set_yticks([])
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, orientation="vertical", ticks=unique)
        cbar.ax.set_yticklabels([str(u) for u in unique])
        cbar.set_label("Cluster (compact)")
        return self._save_plot(fig, filename)

    def _plot_distance_matrix(self, distances, filename, **kwargs):
        logger.info("Plotting distance matrix heatmap...")
        fig = plt.figure(figsize=(10, 8))
        im = plt.imshow(distances, aspect="auto", interpolation="none", cmap=kwargs.get("cmap", "viridis"))
        plt.title(kwargs.get("title", "RMSD Distance Matrix (nm)"))
        plt.xlabel(kwargs.get("xlabel", "Frame"))
        plt.ylabel(kwargs.get("ylabel", "Frame"))
        cbar = plt.colorbar(im, ax=plt.gca())
        cbar.set_label("RMSD (nm)")
        return self._save_plot(fig, filename)

    def _plot_dendrogram(self, linkage_matrix, labels, filename, **kwargs):
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
            logger.debug("Internal node %d: heterogeneous clusters %s", i, branch_labels)
            return "#808080"
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

    def _plot_combined_clusters(self, labels_map: Dict[str, np.ndarray]) -> Optional[Path]:
        """Draw a single scatter summary for all requested methods when in minimal mode."""
        order = [m for m in ["kmeans", "dbscan", "hierarchical"] if m in labels_map]
        if not order:
            return None

        fig, axes = plt.subplots(1, len(order), figsize=(5 * len(order), 4))
        if len(order) == 1:
            axes = [axes]

        for ax, method in zip(axes, order):
            labels = labels_map[method]
            frames = np.arange(len(labels))
            unique = np.sort(np.unique(labels))
            cmap = get_cluster_cmap(len(unique))
            norm = get_discrete_norm(unique)
            ax.scatter(frames, labels, c=labels, cmap=cmap, norm=norm, s=60)
            ax.set_title(method.capitalize())
            ax.set_xlabel("Frame")
            ax.set_ylabel("Cluster")
            ax.grid(alpha=0.3)

        fig.suptitle("Clustering Results", fontweight="bold")
        fig.tight_layout()
        out = self._save_plot(fig, self.combined_plot_name)
        plt.close(fig)
        return out

    def _save_plot(self, fig, name: str):
        """Save the figure as a PNG file in the output directory and log its path."""
        plot_path = self.outdir / f"{name}.png"
        fig.savefig(plot_path, bbox_inches="tight")
        logger.info("Plot saved to %s", plot_path)
        return plot_path

    # --------------------------------- Run ------------------------------------

    def run(self) -> dict:
        """
        Run the clustering analysis for the selected methods.
        """
        if self.results:
            logger.info("Results already computed; returning existing results.")
            return self.results

        try:
            logger.info("Starting clustering analysis...")
            results: Dict[str, Dict] = {}
            self._labels_for_combined = {}

            D = None
            diags = {}
            X_flat = None

            need_distance_matrix = "dbscan" in self.methods or self.feature_mode == "distance"
            if need_distance_matrix:
                D = self._calculate_rmsd_matrix()
                diags = self._distance_diagnostics(D)
                if "dbscan" not in self.methods:
                    logger.info("RMSD matrix computed for distance-based feature mode.")
                if "dbscan" in self.methods and diags:
                    if self.eps >= diags.get("p90", float("inf")):
                        logger.warning(
                            "DBSCAN eps=%.3f nm >= p90=%.3f nm — likely to merge clusters.",
                            self.eps,
                            diags["p90"],
                        )
                    elif self.eps <= diags.get("p25", 0.0):
                        logger.warning(
                            "DBSCAN eps=%.3f nm <= p25=%.3f nm — likely to mark many frames as noise.",
                            self.eps,
                            diags["p25"],
                        )

            if self.feature_mode == "coordinates":
                coords = (
                    self.traj.xyz[:, self.atom_indices, :]
                    if self.atom_indices is not None
                    else self.traj.xyz
                )
                X_flat = coords.reshape(self.traj.n_frames, -1)
                logger.debug("Feature matrix shape: %s", X_flat.shape)
            else:
                logger.debug("Using RMSD distance matrix as feature space for clustering.")

            if any(m in self.methods for m in ("kmeans", "hierarchical")):
                if self.n_clusters is None or self.n_clusters < 1:
                    self.n_clusters = 3
                    logger.info(
                        "n_clusters not provided; defaulting to %d for kmeans/hierarchical.",
                        self.n_clusters,
                    )

            for method in self.methods:
                key = method.lower()
                logger.info("Running method: %s", key)

                if key == "dbscan":
                    if D is None:
                        D = self._calculate_rmsd_matrix()
                        diags = self._distance_diagnostics(D)

                    db = DBSCAN(eps=self.eps, min_samples=self.min_samples, metric="precomputed")
                    labels_raw = db.fit_predict(D).astype(int, copy=False)   # sklearn: -1=noise
                    labels_compact, mapping, noise_label = relabel_compact_positive(labels_raw, start=1, noise_as_last=True)
                    n_clusters = int(len(set(labels_compact)) - (1 if noise_label is not None else 0))

                    logger.info("DBSCAN clusters (excluding noise): %d", n_clusters)
                    if noise_label is not None:
                        logger.info("DBSCAN noise mapped to compact label %d.", noise_label)

                    frame_idx = np.arange(labels_compact.size, dtype=int)
                    dbscan_entry = {
                        "labels_raw": labels_raw,              # may contain -1
                        "labels": labels_compact,              # 1..K (noise -> K+1 if present)
                        "n_clusters": n_clusters,
                        "eps_nm": float(self.eps),
                        "min_samples": int(self.min_samples),
                        "distance_percentiles_nm": diags or {},
                    }
                    if self.save_data:
                        dbscan_entry["distance_matrix"] = D
                        dbscan_entry["labels_file_compact"] = self._save_data(
                            np.column_stack((frame_idx, labels_compact)),
                            "dbscan_labels_compact",
                            header="frame label(1..K; noise=K+1)", fmt="%d",
                        )
                        dbscan_entry["labels_file_raw"] = self._save_data(
                            np.column_stack((frame_idx, labels_raw)),
                            "dbscan_labels_raw",
                            header="frame label_raw(-1=noise)", fmt="%d",
                        )
                    results["dbscan"] = dbscan_entry
                    # Plots use compact labels so colorbars show 1..K(+1)
                    if self.plot_style == "full":
                        results["dbscan"]["pop_plot"] = self._plot_population(labels_compact, "dbscan_pop")
                        results["dbscan"]["trajectory_histogram"] = self._plot_cluster_trajectory_histogram(labels_compact, "dbscan_traj_hist")
                        results["dbscan"]["trajectory_scatter"] = self._plot_cluster_trajectory_scatter(labels_compact, "dbscan_traj_scatter")
                        results["dbscan"]["distance_matrix_plot"] = self._plot_distance_matrix(D, "dbscan_distance_matrix")
                    if self.plot_style in {"minimal", "full"}:
                        self._labels_for_combined["dbscan"] = labels_compact

                elif key == "kmeans":
                    if self.n_clusters is None or self.n_clusters < 1:
                        raise AnalysisError("For KMeans clustering, n_clusters must be provided and >=1.")
                    km = KMeans(n_clusters=int(self.n_clusters), random_state=42, n_init=10)
                    if self.feature_mode == "coordinates":
                        fit_matrix = X_flat
                    else:
                        if D is None:
                            D = self._calculate_rmsd_matrix()
                        fit_matrix = D
                    labels0 = km.fit_predict(fit_matrix).astype(int, copy=False)  # 0..K-1
                    labels = labels0 + 1                                      # 1..K
                    frame_idx = np.arange(labels.size, dtype=int)
                    kmeans_entry = {
                        "labels": labels,
                        "n_clusters": int(self.n_clusters),
                        "inertia_": float(km.inertia_),
                    }
                    if self.save_data:
                        kmeans_entry["labels_file"] = self._save_data(
                            np.column_stack((frame_idx, labels)),
                            "kmeans_labels",
                            header="frame label(1..K)",
                            fmt="%d",
                        )
                        if self.feature_mode == "coordinates" and X_flat is not None:
                            kmeans_entry["coordinates_file"] = self._save_data(
                                X_flat,
                                "kmeans_coordinates",
                                header="Flattened coordinates",
                                fmt="%.6f",
                            )
                    results["kmeans"] = kmeans_entry
                    if self.plot_style == "full":
                        results["kmeans"]["pop_plot"] = self._plot_population(labels, "kmeans_pop")
                        results["kmeans"]["trajectory_histogram"] = self._plot_cluster_trajectory_histogram(labels, "kmeans_traj_hist")
                        results["kmeans"]["trajectory_scatter"] = self._plot_cluster_trajectory_scatter(labels, "kmeans_traj_scatter")
                    if self.plot_style in {"minimal", "full"}:
                        self._labels_for_combined["kmeans"] = labels

                elif key == "hierarchical":
                    if self.n_clusters is None or self.n_clusters < 1:
                        raise AnalysisError("For hierarchical clustering, n_clusters must be provided and >=1.")
                    logger.info("Computing Ward linkage for hierarchical clustering...")
                    if self.feature_mode == "coordinates":
                        Z = linkage(X_flat, method="ward")
                    else:
                        if D is None:
                            D = self._calculate_rmsd_matrix()
                        condensed_dist = squareform(D)
                        Z = linkage(condensed_dist, method="ward")
                    labels = fcluster(Z, t=int(self.n_clusters), criterion="maxclust").astype(int, copy=False)  # 1..K
                    frame_idx = np.arange(labels.size, dtype=int)
                    hierarchical_entry = {
                        "labels": labels,
                        "n_clusters": int(self.n_clusters),
                        "linkage": Z,
                    }
                    if self.save_data:
                        hierarchical_entry["labels_file"] = self._save_data(
                            np.column_stack((frame_idx, labels)), "hierarchical_labels",
                            header="frame label(1..K)", fmt="%d",
                        )
                        hierarchical_entry["linkage_file"] = self._save_data(
                            Z, "hierarchical_linkage",
                            header="cluster1 cluster2 distance sample_count", fmt="%.6f",
                        )
                    results["hierarchical"] = hierarchical_entry
                    if self.plot_style == "full":
                        results["hierarchical"]["pop_plot"] = self._plot_population(labels, "hierarchical_pop")
                        results["hierarchical"]["trajectory_histogram"] = self._plot_cluster_trajectory_histogram(labels, "hierarchical_traj_hist")
                        results["hierarchical"]["trajectory_scatter"] = self._plot_cluster_trajectory_scatter(labels, "hierarchical_traj_scatter")
                        results["hierarchical"]["dendrogram_plot"] = self._plot_dendrogram(Z, labels, "hierarchical_dendrogram")
                    if self.plot_style in {"minimal", "full"}:
                        self._labels_for_combined["hierarchical"] = labels

                else:
                    raise AnalysisError(f"Unknown clustering method: {method}")

            if self.plot_style == "minimal":
                combined = self._plot_combined_clusters(self._labels_for_combined)
                if combined is not None:
                    results["combined_plot"] = combined

            # Release large intermediates when they are no longer needed
            if not self.save_data:
                D = None
                X_flat = None

            self._labels_for_combined = {}

            if self.store_results:
                self.results = results
                logger.info("Clustering analysis complete.")
                return results

            self.results = {}
            logger.info("Clustering analysis complete.")
            return results

        except Exception as e:
            logger.exception("Clustering failed:")
            raise AnalysisError(f"Clustering failed: {str(e)}")

    def plot(self, **kwargs):
        if not self.results:
            raise AnalysisError("No clustering results available. Run the analysis first.")
        return self.results

