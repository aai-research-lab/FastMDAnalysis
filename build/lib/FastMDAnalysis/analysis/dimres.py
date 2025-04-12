"""
Dimensionality Reduction Analysis Module

This module computes 2D embeddings of an MD trajectory using one or more methods:
  - PCA
  - MDS
  - t-SNE

The analysis uses a default atom selection ("protein and name CA") to build a feature matrix by
flattening the 3D coordinates of the selected atoms. For each chosen method, a 2D embedding is computed,
saved to a file, and a default scatter plot is generated.

Usage Example (API):
    from FastMDAnalysis import FastMDAnalysis
    
    fastmda = FastMDAnalysis()
    
    # Run all dimensionality reduction methods (PCA, MDS, t-SNE) using the default atom selection.
    analysis = fastmda.dimred("protein_traj.dcd", "protein.pdb", methods=["pca", "mds", "tsne"])
    data = analysis.data
    
    # Replot the embeddings with custom options.
    analysis.plot(data, title="Custom DimRed Plot", xlabel="X Component", ylabel="Y Component",
                  marker="s", cmap="plasma")
"""

import numpy as np
import mdtraj as md

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE

from pathlib import Path

from .base import BaseAnalysis, AnalysisError

class DimRedAnalysis(BaseAnalysis):
    def __init__(self, trajectory, methods="all", atom_selection="protein and name CA", **kwargs):
        """
        Initialize the Dimensionality Reduction analysis.

        Args:
            trajectory (mdtraj.Trajectory): The MD trajectory.
            methods (str or list): Which reduction methods to use.
                Options: 'pca', 'mds', 'tsne'. If "all" (default), all three are used.
            atom_selection (str): MDTraj selection string to choose atoms.
                Default is "protein and name CA".
            kwargs: Additional arguments passed to BaseAnalysis.
        """
        super().__init__(trajectory, **kwargs)
        # Determine which methods to run.
        if isinstance(methods, str):
            if methods.lower() == "all":
                self.methods = ["pca", "mds", "tsne"]
            else:
                self.methods = [methods.lower()]
        elif isinstance(methods, list):
            self.methods = [m.lower() for m in methods]
        else:
            raise AnalysisError("Parameter 'methods' must be a string or list of strings.")
        
        self.atom_selection = atom_selection
        self.atom_indices = self.traj.topology.select(self.atom_selection)
        if self.atom_indices is None or len(self.atom_indices) == 0:
            raise AnalysisError("No atoms found using the specified atom selection for dimensionality reduction.")
        self.results = {}
        self.data = None  # Will hold a dictionary with embeddings.

    def run(self) -> dict:
        """
        Compute 2D embeddings using the selected methods.

        The coordinates of the selected atoms are flattened into a feature matrix. Then for each method:
          - PCA, MDS, and t-SNE are applied separately.
          - The 2D embedding is saved to the output directory.
          - The embedding is stored in the results dictionary.
          - A default scatter plot is automatically generated for each method.
        
        Returns:
            dict: A dictionary with keys 'pca', 'mds', and/or 'tsne' mapping to the corresponding 2D arrays.
        """
        # Extract coordinates of selected atoms.
        X = self.traj.xyz[:, self.atom_indices, :]  # shape: (n_frames, n_atoms, 3)
        X_flat = X.reshape(self.traj.n_frames, -1)    # shape: (n_frames, n_atoms*3)
        
        for method in self.methods:
            if method == "pca":
                pca = PCA(n_components=2)
                embedding = pca.fit_transform(X_flat)
                self.results["pca"] = embedding
                self._save_data(embedding, "dimred_pca")
            elif method == "mds":
                mds = MDS(n_components=2, random_state=42, dissimilarity="euclidean")
                embedding = mds.fit_transform(X_flat)
                self.results["mds"] = embedding
                self._save_data(embedding, "dimred_mds")
            elif method in ["tsne", "t-sne"]:
                tsne = TSNE(n_components=2, random_state=42, metric="euclidean")
                embedding = tsne.fit_transform(X_flat)
                self.results["tsne"] = embedding
                self._save_data(embedding, "dimred_tsne")
            else:
                raise AnalysisError(f"Unknown dimensionality reduction method: {method}")
        
        self.data = self.results
        # Generate default scatter plots for all computed methods.
        self.plot()
        return self.results

    def plot(self, data=None, method=None, **kwargs):
        """
        Generate scatter plots for the computed 2D embeddings.

        Args:
            data (dict, optional): A dictionary with embeddings. If not provided, uses self.data.
            method (str, optional): If specified, only replot the embedding for that method ('pca', 'mds', or 'tsne').
            kwargs: Custom plot options (e.g., title, xlabel, ylabel, marker, cmap).

        Returns:
            If method is specified, returns the file path for that plot.
            Otherwise, returns a dictionary mapping each method to its plot file path.
        """
        if data is None:
            data = self.data
        if data is None:
            raise AnalysisError("No dimensionality reduction data available. Please run analysis first.")

        def _plot_embedding(embedding, method_name):
            title = kwargs.get("title", f"{method_name.upper()} Projection")
            xlabel = kwargs.get("xlabel", "Component 1")
            ylabel = kwargs.get("ylabel", "Component 2")
            marker = kwargs.get("marker", "o")
            cmap = kwargs.get("cmap", "viridis")
            
            fig = plt.figure(figsize=(10, 8))
            # Color points using frame index.
            colors = np.arange(self.traj.n_frames)
            sc = plt.scatter(embedding[:, 0], embedding[:, 1], c=colors, cmap=cmap, marker=marker)
            plt.title(title)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.grid(alpha=0.3)
            cbar = plt.colorbar(sc)
            cbar.set_label("Frame Index")
            plot_path = self._save_plot(fig, f"dimred_{method_name}")
            plt.close(fig)
            return plot_path

        plot_paths = {}
        if method:
            m = method.lower()
            if m not in data:
                raise AnalysisError(f"Dimensionality reduction method '{m}' not found in results.")
            return _plot_embedding(data[m], m)
        else:
            for m, emb in data.items():
                plot_paths[m] = _plot_embedding(emb, m)
            return plot_paths

