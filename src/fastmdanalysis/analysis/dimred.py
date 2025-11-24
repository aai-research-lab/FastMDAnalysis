from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE

from .base import BaseAnalysis
from ..utils.options import OptionsForwarder
from ..utils.plotting import apply_slide_style


PathLike = Union[str, Path]


def _as_list(methods: Union[str, Sequence[str]]) -> List[str]:
    """
    Normalize methods to a lowercase list in canonical order.
    Supports: "all", comma-separated string, or an iterable.
    """
    if isinstance(methods, str):
        if methods.lower() == "all":
            return ["pca", "mds", "tsne"]
        items = [m.strip().lower() for m in methods.split(",")]
    else:
        items = [str(m).strip().lower() for m in methods]

    if "all" in items:
        return ["pca", "mds", "tsne"]

    order = ["pca", "mds", "tsne"]
    items_set = set(items)
    return [m for m in order if m in items_set]


def _auto_tsne_perplexity(n_frames: int, user_value: Optional[int] = None) -> int:
    """
    Choose a sensible t-SNE perplexity. Clamp to [5, 30] and < n_frames.
    """
    if user_value is not None:
        p = int(user_value)
    else:
        # simple heuristic: ~ min(30, max(5, n/10))
        p = max(5, min(30, n_frames // 10 if n_frames > 0 else 30))
    # t-SNE requires perplexity < n_samples
    p = min(p, max(1, n_frames - 1))
    return max(5, p)


class DimRedAnalysis(BaseAnalysis):
    """
    Dimensionality reduction (PCA, MDS, t-SNE) on Cartesian coordinates.

    Parameters
    ----------
    trajectory : mdtraj.Trajectory
        The trajectory to analyze.
    methods : {"all", "pca", "mds", "tsne"} or list of them
        Which embeddings to compute. "all" runs pca, mds, tsne.
        Alias: method (singular)
    atoms : str, optional
        MDTraj atom selection string.
    outdir : str, optional
        Output directory (default: "dimred_output").
    n_components : int, optional
        Number of components for dimensionality reduction (default 2).
    random_state : int, optional
        Random seed for reproducibility (default 42).
    perplexity : int, optional
        t-SNE perplexity (overrides tsne_perplexity).
    n_iter : int, optional
        t-SNE iterations (deprecated, use max_iter).
    max_iter : int, optional
        t-SNE max iterations (default 500).
    metric : str, optional
        Metric for MDS (default "euclidean"). Use "precomputed" for dissimilarity matrix.
    dissimilarity : str, optional
        Alias for metric in MDS.
    strict : bool
        If True, raise errors for unknown options. If False, log warnings.
    """

    _ALIASES = {
        "method": "methods",
        "atom_indices": "atoms",
        "selection": "atoms",
        "tsne_perplexity": "perplexity",
        "tsne_max_iter": "max_iter",
        "n_iter": "max_iter",  # deprecated but map it
        "dissimilarity": "metric",
    }

    def __init__(
        self,
        trajectory,
        methods: Union[str, Sequence[str]] = "all",
        atoms: Optional[str] = None,
        outdir: Optional[PathLike] = None,
        n_components: int = 2,
        random_state: int = 42,
        perplexity: Optional[int] = None,
        max_iter: int = 500,
        metric: str = "euclidean",
        strict: bool = False,
        **kwargs
    ):
        warn_unknown = kwargs.pop("_warn_unknown", False)

        analysis_opts = {
            "methods": methods,
            "atoms": atoms,
            "output": outdir,
            "n_components": n_components,
            "random_state": random_state,
            "perplexity": perplexity,
            "max_iter": max_iter,
            "metric": metric,
            "strict": strict,
        }
        analysis_opts.update(kwargs)
        
        forwarder = OptionsForwarder(aliases=self._ALIASES, strict=strict)
        resolved = forwarder.apply_aliases(analysis_opts)
        resolved = forwarder.filter_known(
            resolved,
            {
                "methods",
                "atoms",
                "output",
                "n_components",
                "random_state",
                "perplexity",
                "max_iter",
                "metric",
                "strict",
            },
            context="dimred",
            warn=warn_unknown,
        )

        methods = resolved.get("methods", "all")
        atoms = resolved.get("atoms", None)
        outdir = resolved.get("output", None)
        n_components = resolved.get("n_components", 2)
        random_state = resolved.get("random_state", 42)
        perplexity = resolved.get("perplexity", None)
        max_iter = resolved.get("max_iter", 500)
        metric = resolved.get("metric", "euclidean")
        
        # Initialize the base class with the trajectory
        super().__init__(trajectory, output=outdir)
        
        self.atoms = atoms
        self.methods = _as_list(methods)
        self.n_components = int(n_components)
        self.random_state = int(random_state)
        self.tsne_perplexity = perplexity
        self.tsne_max_iter = int(max_iter)
        self.mds_metric = metric
        self.strict = strict

        # Results: method -> ndarray of shape (n_frames, 2)
        self.results: Dict[str, np.ndarray] = {}

    # ------------------------------- helpers ---------------------------------

    def _flatten_xyz(self) -> np.ndarray:
        """
        Return (n_frames, n_atoms*3) float32 array for (optionally) selected atoms.
        """
        t = self.traj
        if self.atoms:
            idx = t.topology.select(self.atoms)
            if idx.size == 0:
                raise ValueError(f"Atom selection returned 0 atoms: {self.atoms}")
            t = t.atom_slice(idx, inplace=False)

        # (n_frames, n_atoms, 3) -> (n_frames, n_atoms * 3)
        X = t.xyz.astype(np.float32).reshape((t.n_frames, -1), order="C")
        return X

    def _save_array(self, name: str, arr: np.ndarray) -> Path:
        """Save array to file using BaseAnalysis method."""
        return self._save_data(arr, f"dimred_{name}")

    def _plot_one(self, name: str, emb: np.ndarray) -> Path:
        """
        Scatter colored by frame index; returns the saved PNG path.
        """
        fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

        # Color by frame index with a colorbar
        c = np.arange(emb.shape[0], dtype=np.int32)
        sc = ax.scatter(emb[:, 0], emb[:, 1], s=20, c=c, cmap="viridis", alpha=0.7)
        cb = fig.colorbar(sc, ax=ax)
        cb.set_label("Frame Index")

        # Set titles and labels based on method
        title_map = {
            "pca": "PCA Projection",
            "mds": "MDS Projection", 
            "tsne": "t-SNE Projection"
        }
        title = title_map.get(name, f"{name.upper()} Projection")
        
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel("Component 1", fontsize=12)
        ax.set_ylabel("Component 2", fontsize=12)
        
        # Add grid for better readability
        ax.grid(True, alpha=0.3)
        
        # Remove top and right spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        apply_slide_style(
            ax,
            x_values=emb[:, 0],
            y_values=emb[:, 1],
            x_max_ticks=6,
            y_max_ticks=6,
        )

        # Save plot using BaseAnalysis method
        out = self._save_plot(fig, f"dimred_{name}")
        plt.close(fig)
        return out

    # --------------------------------- API -----------------------------------

    def run(self) -> "DimRedAnalysis":
        """
        Compute the requested embeddings and save numeric outputs immediately.
        Also logs brief progress via BaseAnalysis logger (if configured).
        """
        X_flat = self._flatten_xyz()
        n_frames = X_flat.shape[0]

        # PCA
        if "pca" in self.methods:
            pca = PCA(n_components=self.n_components, random_state=self.random_state)
            emb = pca.fit_transform(X_flat)
            self.results["pca"] = emb.astype(np.float32)
            self._save_array("pca", self.results["pca"])

        # MDS (silence FutureWarning by setting n_init explicitly)
        if "mds" in self.methods:
            mds = MDS(
                n_components=self.n_components,
                n_init=4,  # default value today; avoids FutureWarning("...will change...")
                random_state=self.random_state,
                normalized_stress="auto",
                dissimilarity=self.mds_metric,
            )
            emb = mds.fit_transform(X_flat)
            self.results["mds"] = emb.astype(np.float32)
            self._save_array("mds", self.results["mds"])

        # t-SNE
        if "tsne" in self.methods:
            perplexity = _auto_tsne_perplexity(n_frames, self.tsne_perplexity)
            tsne = TSNE(
                n_components=self.n_components,
                perplexity=perplexity,
                max_iter=self.tsne_max_iter,  # use max_iter (not deprecated n_iter)
                random_state=self.random_state,
                init="pca",
                learning_rate="auto",
            )
            emb = tsne.fit_transform(X_flat)
            self.results["tsne"] = emb.astype(np.float32)
            self._save_array("tsne", self.results["tsne"])

        # Generate plots automatically after computation
        self.plot()
        
        return self

    def plot(self) -> Dict[str, Path]:
        """
        Generate and save plots for each computed embedding, returning a {method: Path} map.
        """
        out: Dict[str, Path] = {}
        for name, emb in self.results.items():
            out[name] = self._plot_one(name, emb)
        return out
