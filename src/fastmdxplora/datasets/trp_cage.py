"""The Trp-cage miniprotein, as a reference system to point studies at.

A 20-residue miniprotein (PDB 1L2Y) widely used to exercise molecular
dynamics pipelines: small enough to simulate quickly and folded enough for
RMSD, radius of gyration and secondary structure to mean something.

No trajectory is bundled. This module described itself as "a placeholder
class with the expected API surface" whose "actual trajectory data files are
bundled in a future release"; that was written for v0.1.0, two major releases
ago, and ``TrpCage.traj`` meanwhile returned the path to a file that had never
existed. Reading it gave a plausible string, and passing it anywhere gave a
file-not-found from inside a trajectory reader.

Bundling a trajectory in a Python distribution to serve as a reference would
be several megabytes of wheel for something the software can produce: it
fetches 1L2Y and simulates it. So the metadata that is true is kept, and the
paths that were not say so when asked.
"""

from __future__ import annotations

from typing import Any

__all__ = ["TrpCage"]


class _NotBundled:
    """An attribute that explains itself instead of handing out a dead path."""

    def __init__(self, what: str) -> None:
        self._what = what

    def __set_name__(self, owner: type, name: str) -> None:
        self._name = name

    def __get__(self, instance: Any, owner: type | None = None) -> str:
        raise FileNotFoundError(
            f"No {self._what} is bundled for Trp-cage, so there is no "
            f"`TrpCage.{self._name}` to read. FastMDXplora produces one:\n\n"
            "    fastmdx explore --system 1L2Y\n\n"
            "and the trajectory is then at <output>/simulation/production.dcd "
            "with its topology beside it."
        )


class TrpCage:
    """The Trp-cage miniprotein as a reference system.

    Carries what is known about the system. It does not carry a trajectory;
    see the module docstring for why, and for the one command that makes one.
    """

    pdb_id: str = "1L2Y"
    n_residues: int = 20
    description: str = "Trp-cage miniprotein (1L2Y), 20 residues"

    traj = _NotBundled("trajectory")
    top = _NotBundled("topology")
