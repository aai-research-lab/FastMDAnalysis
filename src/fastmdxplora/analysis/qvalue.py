"""Q-value: fraction of native contacts retained per frame.

The Q-value (or Q-fraction) is the canonical folding-state metric in
protein dynamics. For each frame, it reports the fraction of the
reference structure's residue-residue contacts that are still present.
Q ≈ 1 means the fold is intact; Q ≈ 0 means it is fully unfolded.

Both the contacts and the measure follow Best, Hummer, & Eaton (PNAS 2013).
Two residues are natively in contact if any pair of their heavy atoms is
within a cutoff (default 0.45 nm), among residue pairs at least
``min_seq_separation`` apart in sequence (default 4, which excludes
covalent and local contacts that do not probe the tertiary fold).

Each native contact then contributes through a switching function rather
than a threshold::

    Q(X) = 1/|S| * sum over (i,j) in S of
           1 / (1 + exp[beta * (r_ij(X) - lambda * r0_ij)])

with beta = 50 per nm and lambda = 1.8. What that buys is a contact judged
against the distance it had natively, not against one number shared by every
pair: a contact formed at 0.20 nm and one formed at 0.44 nm are each allowed
to stretch by the same factor before they stop counting, and each stops
counting gradually rather than at a step. A single threshold applied to every
pair reports a fold coming apart while the paper's measure still reads it as
intact.

References
----------
Best, R. B.; Hummer, G.; Eaton, W. A. Native contacts determine protein
folding mechanisms in atomistic simulations. *PNAS* **2013**, 110, 17874.
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
from scipy.special import expit

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.analysis.orchestrator import register_analysis


class QValue(Analysis):
    """Fraction of native contacts retained per frame.

    Parameters
    ----------
    ref : int, default 0
        Reference frame defining the "native" state. The contacts present
        in this frame become the denominator of the Q calculation.
    cutoff : float, default 0.45
        Heavy-atom-contact cutoff in nm, used to decide which residue pairs
        are natively in contact in the reference frame. It does not decide
        whether a contact is present later: that is what the switching
        function is for.
    beta : float, default 50.0
        Steepness of the switching function, in 1/nm. The paper's value.
        Larger values approach a hard threshold.
    lambda_factor : float, default 1.8
        How far a contact may stretch relative to its own native distance
        before it stops counting. The paper's value.
    min_seq_separation : int, default 4
        Minimum |i - j| in sequence for a pair to be considered. The
        default 4 excludes local contacts (i±1, i±2, i±3) that don't
        probe the global fold.
    **kwargs
        Standard base-class options.

    Output
    ------
    ``qvalue.dat`` — single-column file of Q per frame (range [0, 1]).
    ``qvalue.png`` — Time series of Q vs. frame/time.
    """

    #: The shortest chain this can say anything about.
    #:
    #: Q measures how much of a fold is intact, from contacts between
    #: residues far enough apart in sequence to probe tertiary structure. A
    #: chain shorter than `min_seq_separation + 1` residues has no such pair,
    #: and no fold for them to describe. Raising there recorded a failed
    #: analysis for a peptide that simply has no tertiary structure -- the
    #: same category error `requires_water` exists to avoid, where "there is
    #: no water here" was a failed phase rather than a question that did not
    #: apply.
    requires_tertiary_structure = True

    name = "qvalue"
    time_series = True
    description = "Native contact fraction (Q)"
    default_selection = None

    def __init__(
        self,
        *,
        ref: int = 0,
        cutoff: float = 0.45,
        beta: float = 50.0,
        lambda_factor: float = 1.8,
        min_seq_separation: int = 4,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        self.ref: int = int(ref)
        self.cutoff: float = float(cutoff)
        self.beta: float = float(beta)
        self.lambda_factor: float = float(lambda_factor)
        self.min_seq_separation: int = int(min_seq_separation)
        self.options.update(
            ref=self.ref,
            cutoff=self.cutoff,
            beta=self.beta,
            lambda_factor=self.lambda_factor,
            min_seq_separation=self.min_seq_separation,
        )

    def compute(self, traj: md.Trajectory) -> np.ndarray:
        """Compute Q per frame.

        Returns
        -------
        np.ndarray of shape (n_frames,)
            Q values in [0, 1]. NaN if the reference has zero contacts
            (the calculation is undefined — e.g. an unfolded reference).
        """
        # Restrict to the selected atoms (e.g. protein/solute) BEFORE
        # enumerating residue pairs. On a solvated system the full trajectory
        # has thousands of water residues; building all-vs-all pairs over them
        # is both meaningless and catastrophically slow (n_res^2). Slicing
        # first drops n_res from ~thousands to the protein's ~hundreds.
        atom_idx = self.select_atoms(traj)
        if len(atom_idx) < traj.n_atoms:
            traj = traj.atom_slice(atom_idx)

        # Build candidate residue pairs satisfying |i-j| >= min_seq_separation
        n_res = traj.n_residues
        ii, jj = np.triu_indices(n_res, k=self.min_seq_separation)
        if len(ii) == 0:
            raise ValueError(
                f"No residue pairs satisfy min_seq_separation={self.min_seq_separation} "
                f"in a trajectory with {n_res} residues."
            )

        pairs = np.column_stack([ii, jj])

        # Resolve negative reference indices
        n_frames = traj.n_frames
        ref = self.ref if self.ref >= 0 else n_frames + self.ref
        if not (0 <= ref < n_frames):
            raise ValueError(
                f"Reference frame {self.ref} is out of range for trajectory "
                f"with {n_frames} frames."
            )

        # Closest-heavy-atom distance between each residue pair, per frame.
        # md.compute_contacts returns (n_frames, n_pairs) distances and the
        # pair indices it actually computed (some may be skipped if no
        # heavy atoms in a residue).
        distances, computed_pairs = md.compute_contacts(
            traj, contacts=pairs, scheme="closest-heavy"
        )

        # Native contact mask: pairs within cutoff in the reference
        ref_dist = distances[ref]
        native_mask = ref_dist < self.cutoff
        n_native = int(native_mask.sum())
        if n_native == 0:
            return np.full(n_frames, np.nan)

        # Each native contact is measured against the distance it had in the
        # reference, through the paper's switching function. Counting instead
        # how many are still inside a single cutoff judges a contact formed at
        # 0.20 nm by the same standard as one formed at 0.44 nm, and reports a
        # fold coming apart well before it has.
        native_distances = distances[:, native_mask]
        r0 = ref_dist[native_mask]
        # expit(-x) is 1/(1+exp(x)) evaluated without overflowing. Written the
        # direct way, a contact stretched much past its native distance sends
        # the exponent over what a float can hold: the answer is still zero,
        # but the warning arrives exactly when someone is watching something
        # come apart, which is when they are least in the mood for it.
        q = np.mean(
            expit(-self.beta * (native_distances - self.lambda_factor * r0)),
            axis=1,
        )
        self._n_native = n_native
        return q.astype(np.float64)

    def plot(self, result: np.ndarray, ax: plt.Axes) -> None:
        x, _ = self.frame_axis_for_plot(self._traj_for_plot, len(result))
        ax.plot(x, result, linewidth=1.4)
        ax.axhline(1.0, color="#888888", linestyle=":", linewidth=0.8)
        ax.set_ylim(-0.02, 1.05)
        # Annotate the number of native contacts in the legend
        n = getattr(self, "_n_native", None)
        if n is not None:
            ax.legend(
                [f"Q (n_native = {n})"],
                loc="best",
            )

    _traj_for_plot: md.Trajectory | None = None

    def run(self, traj: md.Trajectory):
        self._traj_for_plot = traj
        return super().run(traj)

    def frame_axis_for_plot(
        self, traj: md.Trajectory | None, n_points: int
    ) -> tuple[np.ndarray, str]:
        if traj is None:
            return np.arange(n_points), "Frame"
        return self.frame_axis(traj)

    def default_xlabel(self) -> str | None:
        if self._traj_for_plot is None:
            return "Frame"
        _, label = self.frame_axis(self._traj_for_plot)
        return label

    def default_ylabel(self) -> str | None:
        return "Q (fraction native contacts)"


register_analysis(QValue.name, QValue)
