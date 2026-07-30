"""Conformational clustering.

Clusters frames of the trajectory by structural similarity, producing
one cluster-membership labeling per requested method. Three methods are
supported:

  - **k-means** (default, fast, requires choosing ``n_clusters``)
  - **hierarchical** (agglomerative, Ward linkage by default)
  - **dbscan** (density-based, no pre-specified cluster count)

All methods operate on the pairwise RMSD distance matrix (computed via
MDTraj's QCP algorithm), which is the standard featurization for
conformational clustering. Outputs per method:

  - ``cluster_<method>.dat`` — per-frame integer cluster labels.
  - ``cluster_<method>.png`` — cluster labels as a function of time.
  - ``cluster_<method>_counts.png`` — cluster population bar chart.
  - ``cluster_hierarchical_dendrogram.png`` — hierarchical dendrogram
    when hierarchical clustering is requested and SciPy is available.
  - ``hierarchical_distance_matrix.npy`` and ``hierarchical_linkage.npy`` —
    reproducibility data for dashboard/report-native dendrogram rendering.

Because this analysis produces multiple files per run, it overrides the
base class's :meth:`save_data` and :meth:`_do_plot` methods.

References
----------
Daura, X. et al. *J. Mol. Graph. Model.* **1999**, 18, 122 (RMSD-based MD clustering).
Lloyd, S. *IEEE Trans. Inf. Theory* **1982**, 28, 129 (k-means).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN, AgglomerativeClustering, KMeans

from fastmdxplora.analysis.base import Analysis, AnalysisResult
from fastmdxplora.analysis.orchestrator import register_analysis
from fastmdxplora.analysis.plotting import new_figure, save_figure


VALID_METHODS = ("kmeans", "hierarchical", "dbscan")


class Cluster(Analysis):
    """Conformational clustering by pairwise RMSD.

    Parameters
    ----------
    methods : list of str, default ``["kmeans"]``
        Which clustering algorithms to run. Each method produces its own
        output files. Valid values: ``"kmeans"``, ``"hierarchical"``,
        ``"dbscan"``.
    n_clusters : int, default 5
        Number of clusters (used by k-means and hierarchical).
    eps : float, default 0.2
        DBSCAN distance threshold in nm.
    min_samples : int, default 5
        DBSCAN minimum samples per cluster.
    linkage : {"ward", "complete", "average", "single"}, default "average"
        Hierarchical linkage method. ``"ward"`` requires Euclidean
        distance metric, which is not available when using a precomputed
        RMSD matrix — ``"average"`` is the safe default for RMSD-based
        clustering.
    selection : str, optional
        MDTraj atom selection for the RMSD calculation. Defaults to
        ``"name CA"`` (CA-only is fast and capture the global fold well).
    **kwargs
        Standard base-class options.

    Output
    ------
    Per method, in ``<output_dir>/cluster/``:
      - ``cluster_<method>.dat`` — CSV with ``frame, cluster`` columns.
      - ``cluster_<method>.png`` — Cluster timeline figure.
    """

    name = "cluster"
    description = "Conformational clustering"
    default_selection = "name CA"

    def __init__(
        self,
        *,
        methods: list[str] | None = None,
        n_clusters: int = 5,
        eps: float = 0.2,
        min_samples: int = 5,
        linkage: str = "average",
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        methods = list(methods) if methods else ["kmeans", "hierarchical"]
        methods = [m.lower() for m in methods]
        unknown = [m for m in methods if m not in VALID_METHODS]
        if unknown:
            raise ValueError(
                f"Unknown clustering method(s): {unknown}. Valid: {VALID_METHODS}"
            )
        self.methods: list[str] = methods
        self.n_clusters: int = int(n_clusters)
        self.eps: float = float(eps)
        self.min_samples: int = int(min_samples)
        self.linkage: str = str(linkage)
        self.options.update(
            methods=self.methods,
            n_clusters=self.n_clusters,
            eps=self.eps,
            min_samples=self.min_samples,
            linkage=self.linkage,
        )

    def compute(self, traj: md.Trajectory) -> dict[str, np.ndarray]:
        """Run all requested clustering methods.

        Returns
        -------
        dict
            Maps method name → array of per-frame cluster labels (int).
            DBSCAN's ``-1`` label indicates "noise" (unclustered frames).
        """
        # Build the pairwise RMSD distance matrix once and share across methods.
        atom_idx = self.select_atoms(traj)
        distances = _pairwise_rmsd(traj, atom_idx)
        self._distances = distances  # cached for the plot

        # k-means and hierarchical clustering need at least as many frames
        # (samples) as requested clusters. Fail with an actionable message
        # rather than letting scikit-learn raise an opaque internals error.
        n_frames = traj.n_frames
        partitioning = [m for m in self.methods if m in ("kmeans", "hierarchical")]
        if partitioning and n_frames < self.n_clusters:
            raise ValueError(
                f"Clustering needs at least n_clusters={self.n_clusters} "
                f"frames, but the trajectory has only {n_frames}. Use a longer "
                f"trajectory, or set a smaller n_clusters (e.g. n_clusters="
                f"{max(2, n_frames)} or fewer)."
            )

        results: dict[str, np.ndarray] = {}
        for method in self.methods:
            if method == "kmeans":
                results[method] = _cluster_kmeans(distances, self.n_clusters)
            elif method == "hierarchical":
                results[method] = _cluster_hierarchical(
                    distances, self.n_clusters, self.linkage
                )
            elif method == "dbscan":
                results[method] = _cluster_dbscan(
                    distances, self.eps, self.min_samples
                )

        return results

    # --------------------------------------------------------------
    # Override run() to handle multi-output (one file pair per method)
    # --------------------------------------------------------------
    def run(self, traj: md.Trajectory) -> AnalysisResult:
        from datetime import datetime, timezone

        started = datetime.now(timezone.utc).isoformat()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        options_path = self._write_options_manifest()

        try:
            self.result = self.compute(traj)
            artifacts: list[Path] = []
            for method, labels in self.result.items():
                # Data file: frame, cluster
                data_path = self.output_dir / f"cluster_{method}.dat"
                df = pd.DataFrame(
                    {"frame": np.arange(len(labels)), "cluster": labels}
                )
                df.to_csv(data_path, index=False)
                artifacts.append(data_path)

                # Figure
                fig_path = self.output_dir / f"cluster_{method}.png"
                fig, ax = new_figure(
                    title=f"{self.figure_title()} ({method})",
                    figsize=self._user_figsize,
                )
                _plot_cluster_timeline(ax, labels, method)
                xlabel = self._user_xlabel or "Frame"
                ylabel = self._user_ylabel or "Cluster"
                ax.set_xlabel(xlabel)
                ax.set_ylabel(ylabel)
                save_figure(fig, fig_path)
                artifacts.append(fig_path)
                svg_path = fig_path.with_suffix(".svg")
                if svg_path.is_file():
                    artifacts.append(svg_path)

                counts_path = self.output_dir / f"cluster_{method}_counts.png"
                fig_counts, ax_counts = new_figure(
                    title=f"Cluster populations ({method})",
                    figsize=(5.8, 3.6),
                )
                _plot_cluster_counts(ax_counts, labels)
                save_figure(fig_counts, counts_path)
                artifacts.append(counts_path)
                counts_svg_path = counts_path.with_suffix(".svg")
                if counts_svg_path.is_file():
                    artifacts.append(counts_svg_path)

                if method == "hierarchical":
                    dendro_path = self.output_dir / "cluster_hierarchical_dendrogram.png"
                    skip_path = self.output_dir / "cluster_hierarchical_dendrogram_skipped.json"
                    distance_path = self.output_dir / "hierarchical_distance_matrix.npy"
                    np.save(distance_path, self._distances)
                    artifacts.append(distance_path)
                    try:
                        linkage_matrix = _hierarchical_linkage_matrix(
                            self._distances,
                            self.linkage,
                        )
                        linkage_path = self.output_dir / "hierarchical_linkage.npy"
                        np.save(linkage_path, linkage_matrix)
                        artifacts.append(linkage_path)
                        fig_den, ax_den = new_figure(
                            title="Hierarchical clustering dendrogram",
                            figsize=(6.5, 3.8),
                        )
                        _plot_hierarchical_dendrogram_from_linkage(ax_den, linkage_matrix)
                        save_figure(fig_den, dendro_path)
                        artifacts.append(dendro_path)
                        dendro_svg_path = dendro_path.with_suffix(".svg")
                        if dendro_svg_path.is_file():
                            artifacts.append(dendro_svg_path)
                        if skip_path.exists():
                            skip_path.unlink()
                    except Exception as dendro_exc:  # noqa: BLE001
                        skip_path.write_text(
                            json.dumps(
                                {
                                    "artifact": dendro_path.name,
                                    "status": "skipped",
                                    "reason": str(dendro_exc),
                                },
                                indent=2,
                            )
                            + "\n",
                            encoding="utf-8",
                        )
                        artifacts.append(skip_path)

            finished = datetime.now(timezone.utc).isoformat()

            # Use the first method's outputs as the "primary" data/figure
            # paths in the AnalysisResult; all are listed in the artifacts
            # field of the analysis manifest.
            primary_method = self.methods[0]
            return AnalysisResult(
                name=self.name,
                status="ok",
                data=self.result,
                output_dir=self.output_dir,
                figure_path=self.output_dir / f"cluster_{primary_method}.png",
                data_path=self.output_dir / f"cluster_{primary_method}.dat",
                options_path=options_path,
                artifacts=[options_path, *artifacts],
                message=f"{self.name}: ok ({', '.join(self.methods)})",
                started_at=started,
                finished_at=finished,
            )
        except Exception as exc:  # noqa: BLE001
            finished = datetime.now(timezone.utc).isoformat()
            return AnalysisResult(
                name=self.name,
                status="error",
                output_dir=self.output_dir,
                options_path=options_path,
                message=f"{self.name}: {exc}",
                started_at=started,
                finished_at=finished,
            )

    # Required by the ABC but not used (we override run())
    def plot(self, result: dict[str, np.ndarray], ax: plt.Axes) -> None:
        # Plot the first method on the supplied axes — used only by
        # external callers who instantiate the figure themselves.
        primary = next(iter(result))
        _plot_cluster_timeline(ax, result[primary], primary)


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------
def _pairwise_rmsd(traj: md.Trajectory, atom_idx: np.ndarray) -> np.ndarray:
    """Compute the n_frames × n_frames pairwise RMSD matrix.

    Uses MDTraj's QCP-based ``md.rmsd`` row by row, which aligns each
    frame to the reference before measuring. Returns a symmetric matrix
    in nm.
    """
    n = traj.n_frames
    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        # Aligning frame i against the full trajectory once gives row i
        dist[i] = md.rmsd(traj, traj, frame=i, atom_indices=atom_idx)
    # Symmetrize (MDTraj's RMSD is symmetric up to float precision)
    dist = 0.5 * (dist + dist.T)
    np.fill_diagonal(dist, 0.0)
    return dist


def _cluster_kmeans(distances: np.ndarray, n_clusters: int) -> np.ndarray:
    """K-means on multidimensional scaling embedding of the RMSD matrix."""
    # K-means needs a feature vector representation; project the distance
    # matrix into a Euclidean space via classical MDS (cmdscale).
    embedding = _classical_mds(distances, n_components=min(10, distances.shape[0] - 1))
    model = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
    return model.fit_predict(embedding).astype(int)


def _cluster_hierarchical(
    distances: np.ndarray, n_clusters: int, linkage: str
) -> np.ndarray:
    """Agglomerative hierarchical clustering on the precomputed distance matrix."""
    # `ward` linkage requires a Euclidean (not precomputed) distance.
    # When the user requests ward, embed via MDS first.
    if linkage == "ward":
        embedding = _classical_mds(
            distances, n_components=min(10, distances.shape[0] - 1)
        )
        model = AgglomerativeClustering(
            n_clusters=n_clusters, linkage="ward"
        )
        return model.fit_predict(embedding).astype(int)

    model = AgglomerativeClustering(
        n_clusters=n_clusters,
        metric="precomputed",
        linkage=linkage,
    )
    return model.fit_predict(distances).astype(int)


def _cluster_dbscan(
    distances: np.ndarray, eps: float, min_samples: int
) -> np.ndarray:
    """DBSCAN on the precomputed RMSD matrix.

    Frames classified as noise receive label -1. The number of clusters
    is determined by the algorithm based on ``eps`` and ``min_samples``.
    """
    model = DBSCAN(eps=eps, min_samples=min_samples, metric="precomputed")
    return model.fit_predict(distances).astype(int)


def _classical_mds(distances: np.ndarray, n_components: int) -> np.ndarray:
    """Classical multidimensional scaling on a distance matrix.

    Returns an (n_samples, n_components) embedding in which Euclidean
    distance approximates the input distance. Standard cmdscale formula:
        B = -0.5 * H * D^2 * H,  where H = I - 1/n
    eigendecompose B; coordinates = U * sqrt(lambda).
    """
    n = distances.shape[0]
    d2 = distances ** 2
    h = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * h @ d2 @ h
    eigvals, eigvecs = np.linalg.eigh(b)
    # Take the top-k positive eigenvalues
    idx = np.argsort(eigvals)[::-1][:n_components]
    keep_vals = np.maximum(eigvals[idx], 0)
    return eigvecs[:, idx] * np.sqrt(keep_vals)


def _plot_cluster_timeline(ax: plt.Axes, labels: np.ndarray, method: str) -> None:
    """Plot per-frame cluster labels as a scatter / step plot."""
    frames = np.arange(len(labels))
    unique = sorted(set(labels))
    # Compact colormap; -1 (noise) gets gray
    palette = plt.cm.tab10.colors
    for k in unique:
        mask = labels == k
        color = "#BBBBBB" if k == -1 else palette[k % len(palette)]
        ax.scatter(
            frames[mask],
            labels[mask],
            s=8,
            c=[color],
            edgecolor="none",
            label=f"cluster {k}" if k != -1 else "noise",
        )
    ax.set_yticks(unique)
    if len(unique) <= 10:
        ax.legend(loc="best", fontsize=8, ncol=2)


def _plot_cluster_counts(ax: plt.Axes, labels: np.ndarray) -> None:
    unique, counts = np.unique(labels, return_counts=True)
    palette = plt.cm.tab10.colors
    colors = ["#BBBBBB" if k == -1 else palette[int(k) % len(palette)] for k in unique]
    ax.bar(unique, counts, color=colors, width=0.75)
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Frames")
    ax.set_xticks(unique)


def _plot_hierarchical_dendrogram(
    ax: plt.Axes,
    distances: np.ndarray,
    linkage_method: str,
) -> None:
    z = _hierarchical_linkage_matrix(distances, linkage_method)
    _plot_hierarchical_dendrogram_from_linkage(ax, z)


def _hierarchical_linkage_matrix(
    distances: np.ndarray,
    linkage_method: str,
) -> np.ndarray:
    try:
        from scipy.cluster.hierarchy import linkage
        from scipy.spatial.distance import squareform
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("SciPy is required for dendrogram generation") from exc

    if distances.shape[0] < 2:
        raise ValueError("at least two frames are required for a dendrogram")
    condensed = squareform(distances, checks=False)
    method = "average" if linkage_method == "ward" else linkage_method
    return linkage(condensed, method=method)


def _plot_hierarchical_dendrogram_from_linkage(
    ax: plt.Axes,
    linkage_matrix: np.ndarray,
) -> None:
    try:
        from scipy.cluster.hierarchy import dendrogram
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("SciPy is required for dendrogram generation") from exc

    dendrogram(linkage_matrix, ax=ax, no_labels=True, color_threshold=None)
    ax.set_xlabel("Frame")
    ax.set_ylabel("Distance")
    ax.set_ylim(bottom=0)


register_analysis(Cluster.name, Cluster)
