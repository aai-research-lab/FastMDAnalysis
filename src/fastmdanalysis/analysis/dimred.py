# FastMDAnalysis/src/fastmdanalysis/analysis/dimred.py
"""
Dimensionality Reduction Analysis Module

Computes 2D embeddings of an MD trajectory using one or more methods:
  - PCA
  - MDS
  - t-SNE

Feature construction:
  - Uses the provided atom selection (if any) to slice the trajectory, then flattens
    the 3D coordinates per frame → (n_frames, n_selected_atoms * 3).

For each selected method, a 2D embedding is computed, saved to <outdir>/dimred_<method>.dat,
and a scatter plot colored by frame index is saved to <outdir>/dimred_<method>.png.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Union
import logging
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE

from .base import BaseAnalysis, AnalysisError

logger = logging.getLogger(__name__)


def _normalize_methods(methods: Union[str, Sequence[str]]) -> Sequence[str]:
    """Return a normalized list of methods. 'all' expands to ['pca', 'mds', 'tsne']."""
    if isinstance(methods, str):
        mlist = [methods.lower()]
    else:
        mlist = [m.lower() for m in methods]
    return ["pca", "mds", "tsne"] if ("all" in mlist) else mlist


class DimRedAnalysis(BaseAnalysis):
    def __init__(self, trajectory, methods: Union[str, Sequence[str]] = "all",
                 atoms: Optional[str] = None, **kwargs):
        """
        Parameters
        ----------
        trajectory : mdtraj.Trajectory
            The MD trajectory to analyze.
        methods : {'all', 'pca', 'mds', 'tsne'} or list of these
            Which dimensionality reduction method(s) to run. Default: "all".
        atoms : str or None
            MDTraj selection string used to build the feature matrix. If None, all atoms are used.
        kwargs : dict
            Passed to BaseAnalysis (e.g., output directory).
        """
        super().__init__(trajectory, **kwargs)
        self.methods = list(_normalize_methods(methods))
        self.atoms = atoms

        # Prepare the sub-trajectory used to build features
        if self.atoms:
            sel = self.traj.topology.select(self.atoms)
            if sel is None or len(sel) == 0:
                raise AnalysisError(f"No atoms found for selection: '{self.atoms}'")
            self._feature_traj = self.traj.atom_slice(sel)
        else:
            self._feature_traj = self.traj

        self.results: Dict[str, np.ndarray] = {}
        self.data: Optional[Dict[str, np.ndarray]] = None

    # --------------------------------------------------------------------- run

    def run(self) -> Dict[str, np.ndarray]:
        """
        Compute 2D embeddings using the selected methods.

        Returns
        -------
        dict
            Keys: "pca", "mds", "tsne" (subset depending on methods), values: (n_frames, 2) arrays.
        """
        logger.info(
            "DimRed: starting (methods=%s, atoms=%s, n_frames=%d, n_atoms=%d)",
            ",".join(self.methods),
            self.atoms if self.atoms else "ALL",
            self._feature_traj.n_frames,
            self._feature_traj.n_atoms,
        )

        # Build feature matrix
        X = self._feature_traj.xyz  # (T, N, 3)
        X_flat = X.reshape(self._feature_traj.n_frames, -1)  # (T, N*3)
        T = X_flat.shape[0]

        # PCA
        if "pca" in self.methods:
            pca = PCA(n_components=2, random_state=42)
            emb = pca.fit_transform(X_flat)
            self.results["pca"] = emb
            self._save_data(emb, "dimred_pca", header="x y", fmt="%.6f")

        # MDS (metric, euclidean)
        if "mds" in self.methods:
            # Modern sklearn: prefer explicit params to avoid deprecation warnings.
            mds = MDS(
                n_components=2,
                metric=True,
                n_init=4,
                max_iter=300,
                random_state=42,
                dissimilarity="euclidean",
                normalized_stress="auto"  # handled if available in installed sklearn
            )
            emb = mds.fit_transform(X_flat)
            self.results["mds"] = emb
            self._save_data(emb, "dimred_mds", header="x y", fmt="%.6f")

        # t-SNE
        if "tsne" in self.methods:
            # Robust perplexity for small T:
            # Must be < T and roughly less than T/3; keep in [5, 30].
            if T <= 5:
                perplexity = max(2, T - 1)  # minimal viable
            else:
                perplexity = int(np.clip(T // 10, 5, 30))
            tsne = TSNE(
                n_components=2,
                perplexity=perplexity,
                n_iter=500,             # modest iterations for CI speed
                learning_rate="auto",
                init="pca",
                random_state=42,
                metric="euclidean",
                verbose=0,
            )
            emb = tsne.fit_transform(X_flat)
            self.results["tsne"] = emb
            self._save_data(emb, "dimred_tsne", header="x y", fmt="%.6f")
            logger.info("t-SNE used perplexity=%d (T=%d).", perplexity, T)

        self.data = self.results

        # Generate default plots
        self.plot()

        logger.info("DimRed: done (methods=%s).", ",".join(self.results.keys()))
        return self.results

    # -------------------------------------------------------------------- plot

    def plot(self, data: Optional[Dict[str, np.ndarray]] = None,
             method: Optional[str] = None, **kwargs):
        """
        Generate scatter plots for the 2D embeddings.

        Parameters
        ----------
        data : dict, optional
            If provided, mapping method -> (n_frames, 2) embedding arrays. Defaults to self.data.
        method : {'pca','mds','tsne'} or None
            If specified, plot only that method; else plot all available embeddings.
        kwargs : dict
            title, xlabel, ylabel, marker, cmap (matplotlib options).

        Returns
        -------
        dict | str
            Paths of saved plot(s).
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No dimensionality reduction data available to plot. Run the analysis first.")

        def _plot_one(embedding: np.ndarray, method_name: str):
            title = kwargs.get("title", f"{method_name.upper()} Projection")
            xlabel = kwargs.get("xlabel", "Component 1")
            ylabel = kwargs.get("ylabel", "Component 2")
            marker = kwargs.get("marker", "o")
            cmap = kwargs.get("cmap", "viridis")

            fig, ax = plt.subplots(figsize=(10, 8))
            colors = np.arange(self.traj.n_frames)
            sc = ax.scatter(embedding[:, 0], embedding[:, 1], c=colors, cmap=cmap, marker=marker)
            ax.set_title(title)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.grid(alpha=0.3)
            cbar = fig.colorbar(sc, ax=ax)
            cbar.set_label("Frame")
            fig.tight_layout()
            out = self._save_plot(fig, f"dimred_{method_name}")
            plt.close(fig)
            return out

        if method is not None:
            m = method.lower()
            if m not in data:
                raise AnalysisError(f"Dimensionality reduction method '{m}' not found in results.")
            return _plot_one(data[m], m)

        paths = {}
        for m, emb in data.items():
            paths[m] = _plot_one(emb, m)
        return paths
