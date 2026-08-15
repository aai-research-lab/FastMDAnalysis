"""Typed protein-ligand interactions, frame by frame.

Counting contacts says how much of the protein a ligand touches. This says
what is holding it: a salt bridge a charge change would destroy, a hydrophobic
packing that tolerates one, a hydrogen bond to a backbone that a substitution
cannot reach. That difference is what a medicinal chemist is asking.

Eight interaction types, each implemented against a published criterion named
in the rule's own docstring, in ``interactions``. What can take part is
chemistry, and that is settled first: see ``ligand_chemistry``, which records
how confidently it knew.

The occupancy of each interaction is reported with the observation behind it.
A contact present in three frames of five hundred and one present in four
hundred and fifty are both "present"; only one means anything, and a fraction
alone does not say which is which.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import mdtraj as md
import numpy as np
import pandas as pd

from fastmdxplora.analysis.base import Analysis
from fastmdxplora.analysis.orchestrator import register_analysis

__all__ = ["ProteinLigandInteractions", "ALL_KINDS"]

#: Every rule, in the order a reader would want them: what holds the ligand
#: most often first. At module level so it can be the declared default of the
#: ``kinds`` option rather than something filled in from None afterwards.
ALL_KINDS = (
    "hydrophobic", "hydrogen_bond", "salt_bridge", "pi_stacking",
    "pi_cation", "halogen_bond", "metal_coordination", "water_bridge",
)


class ProteinLigandInteractions(Analysis):
    """What holds the ligand in place, and how well each contact is observed.

    Parameters
    ----------
    ligand_resname : str
        Ligand residue name. Supplied by the orchestrator.
    protein_selection : str, default "protein"
        The other side of the interaction. Usually the whole protein, and
        worth changing in three cases: a complex where only one chain matters
        (``"protein and chainid 0"``), a domain rather than the whole fold
        (``"protein and resid 40 to 180"``), or a receptor that is not a
        protein at all -- ``"nucleic"`` for a ligand bound to DNA or RNA,
        which is a real and common case that the default name does not
        suggest.
    ligand_chemistry : path, optional
        An SDF stating the ligand's chemistry. Where a trajectory came from
        elsewhere and its residue name is not one the Chemical Component
        Dictionary knows, the bond orders otherwise have to be inferred from
        the coordinates -- which is a guess, and a wrong bond order moves a
        hydrogen and invents or destroys a hydrogen bond.
    ligand_net_charge : int, optional
        The ligand's net charge, where you know it. It matters for the
        interactions that are claims about charge: perceiving it from
        coordinates is ambiguous more often than not, and guanidinium comes out
        as an anion if the first balancing charge is taken.
    kinds : sequence of str, default all eight
        Which interactions to look for. Declared as the tuple rather than
        filled in from None, so anything reading the signature -- a form
        drawing a control, the help, a config template -- can say what it
        would do.
    minimum_occupancy : float, default 0.1
        How often an interaction must appear before it counts towards a binding
        mode. Below this, a single fleeting contact splits one arrangement into
        two modes that are really one.
    periodic : bool, default True
        Measure across the periodic boundary where the trajectory carries a
        unit cell.
    **kwargs
        Standard base-class options.

    Outputs
    -------
    ``pl_interactions.dat``
        One row per interaction, with its occupancy, the number of frames it
        was present, and the number of separate times it formed.
    ``pl_interactions.png``
        Occupancy per interaction, drawn so that thinly observed contacts are
        visibly thinly observed.
    """

    name = "pl_interactions"
    description = "Protein-ligand interactions by type"
    requires_ligand = True
    default_selection = None
    #: This works out its own atoms, so a general selection has nothing to
    #: apply to.
    honours_selection = False

    ALL_KINDS = ALL_KINDS

    def __init__(
        self,
        *,
        ligand_resname: str | None = None,
        protein_selection: str = "protein",
        ligand_chemistry: str | None = None,
        ligand_net_charge: int | None = None,
        kinds: Any = ALL_KINDS,
        minimum_occupancy: float = 0.1,
        periodic: bool = True,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        if not ligand_resname:
            raise ValueError(
                "ProteinLigandInteractions requires `ligand_resname`; it "
                "applies only to protein-ligand complexes."
            )
        self.ligand_resname = str(ligand_resname)
        self.protein_selection = str(protein_selection)
        self.ligand_chemistry = ligand_chemistry
        self.ligand_net_charge = (None if ligand_net_charge is None
                                  else int(ligand_net_charge))
        chosen = tuple(kinds) if kinds else ALL_KINDS
        unknown = [k for k in chosen if k not in self.ALL_KINDS]
        if unknown:
            raise ValueError(
                f"Unknown interaction type(s) {unknown}. "
                f"Valid: {', '.join(self.ALL_KINDS)}"
            )
        self.kinds = list(chosen)
        self.minimum_occupancy = float(minimum_occupancy)
        self.periodic = bool(periodic)
        self.options.update(
            ligand_resname=self.ligand_resname,
            protein_selection=self.protein_selection,
            kinds=self.kinds,
            minimum_occupancy=self.minimum_occupancy,
            periodic=self.periodic,
        )

    def compute(self, traj: md.Trajectory) -> pd.DataFrame:
        from fastmdxplora.analysis import interaction_summary as summary
        from fastmdxplora.analysis import interactions as rules
        from fastmdxplora.analysis.ligand_chemistry import resolve_ligand_chemistry

        ligand = traj.topology.select(f"resname {self.ligand_resname}")
        protein = traj.topology.select(self.protein_selection)
        if len(ligand) == 0:
            raise ValueError(
                f"No atoms match resname {self.ligand_resname!r}."
            )
        if len(protein) == 0:
            raise ValueError(
                f"No atoms match {self.protein_selection!r}."
            )
        water = traj.topology.select("water")

        chemistry = resolve_ligand_chemistry(
            traj, self.ligand_resname, ligand,
            supplied=self.ligand_chemistry,
            run_dir=self._run_directory(),
            net_charge=self.ligand_net_charge,
        )
        # Recorded beside the results, because an interaction computed from
        # perceived bond orders is a different claim from one computed from
        # chemistry that was resolved, and a reader deserves to know which.
        self.findings["ligand_chemistry"] = chemistry.as_record()

        found: list[Any] = []
        refused: dict[str, str] = {}

        def attempt(kind: str, call) -> None:
            if kind not in self.kinds:
                return
            try:
                found.extend(call())
            except ValueError as refusal:
                # A rule that cannot answer says so, and the run continues:
                # a ligand whose charge is unknown still has hydrogen bonds.
                refused[kind] = str(refusal)

        attempt("hydrogen_bond", lambda: rules.hydrogen_bonds(
            traj, ligand, protein, periodic=self.periodic))
        attempt("hydrophobic", lambda: rules.hydrophobic_contacts(
            traj, ligand, protein, periodic=self.periodic))
        attempt("salt_bridge", lambda: rules.salt_bridges(
            traj, chemistry, ligand, protein))
        attempt("pi_stacking", lambda: rules.pi_stacking(
            traj, chemistry, ligand, protein))
        attempt("pi_cation", lambda: rules.pi_cation(
            traj, chemistry, ligand, protein))
        attempt("halogen_bond", lambda: rules.halogen_bonds(
            traj, chemistry, ligand, protein))
        attempt("metal_coordination", lambda: rules.metal_coordination(
            traj, ligand, protein))
        attempt("water_bridge", lambda: rules.water_bridges(
            traj, ligand, protein, water))

        if refused:
            self.findings["not_measured"] = refused

        # Residues the charge and ring tables do not know. Reported rather
        # than refused: one modified residue in a large protein is a footnote,
        # and a selection made entirely of nucleotides is not, and only the
        # person who made the selection can tell which this is.
        unknown = rules.residues_not_covered(traj.topology, protein)
        if unknown:
            self.options["residues_not_examined_for_charge_or_rings"] = unknown

        occupancies = summary.occupancies(found, traj.n_frames)
        modes = summary.binding_modes(
            found, traj.n_frames, minimum_occupancy=self.minimum_occupancy)
        transitions = summary.mode_transitions(modes["per_frame"])

        self.findings["binding_modes"] = modes["modes"][:10]
        self.findings["mode_transitions"] = {
            k: v for k, v in transitions.items() if k != "counts"
        }
        self._occupancies = occupancies

        # The exact residue-level table, taken while the per-frame masks are
        # still in hand. Written beside the pair table rather than replacing
        # it: the pairs say which atoms touch, and this says how often the
        # residue does, and neither can be derived from the other.
        self._by_residue = summary.residue_occupancies(
            found, traj.n_frames,
            lambda index: str(traj.topology.atom(index).residue))

        if not occupancies:
            return pd.DataFrame(columns=[
                "kind", "ligand_atom", "protein_atom", "residue",
                "frames_present", "frames_total", "occupancy", "episodes",
                "standard_error", "well_sampled",
            ])

        rows = []
        for entry in occupancies:
            record = entry.as_record()
            record["residue"] = str(traj.topology.atom(entry.protein_atom).residue)
            rows.append(record)
        frame = pd.DataFrame(rows)
        return frame.rename(columns={"fraction": "occupancy"})

    def save_data(self, result: "pd.DataFrame", path: Path) -> Path:
        """The pair table, and the exact residue table beside it.

        A residue's occupancy is the union of its pairs' frames, which the
        pair table cannot express: read from it, a residue touched through
        eight atoms is somewhere between the largest single pair and the
        sum of them all. The union is knowable at the point the contacts
        are counted, so it is written rather than left to be estimated.
        """
        written = super().save_data(result, path)

        rows = getattr(self, "_by_residue", None)
        if rows:
            import pandas as _pd

            beside = Path(path).with_name("pl_interactions_by_residue.dat")
            _pd.DataFrame(rows).to_csv(beside, index=False)
        return written

    def _run_directory(self) -> Path | None:
        """The run this analysis belongs to, if it is one of ours.

        Its setup phase may already have resolved the ligand's chemistry, which
        beats inferring it from coordinates.
        """
        output = getattr(self, "output_dir", None)
        if not output:
            return None
        here = Path(output)
        for candidate in (here, *here.parents):
            if (candidate / "setup" / "ligands").is_dir():
                return candidate
        return None

    def plot(self, result: pd.DataFrame, ax) -> None:
        if result.empty:
            ax.text(0.5, 0.5, "No interactions found", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_axis_off()
            return

        top = result.head(20).iloc[::-1]
        labels = [f"{row.kind}  {row.residue}" for row in top.itertuples()]
        positions = np.arange(len(top))

        # Thinly observed contacts are drawn hollow. An occupancy resting on
        # one observation should not look like one resting on four hundred,
        # and a footnote is easier to miss than a difference in the bar.
        for position, row in zip(positions, top.itertuples()):
            well = bool(row.well_sampled)
            ax.barh(position, row.occupancy,
                    color="#3fb0ac" if well else "none",
                    edgecolor="#3fb0ac", hatch=None if well else "///")
            error = row.standard_error
            if error is not None and not (isinstance(error, float) and np.isnan(error)):
                ax.errorbar(row.occupancy, position, xerr=error,
                            color="#20504f", capsize=3, fmt="none")

        ax.set_yticks(positions)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlim(0, 1)
        thin = int((~top["well_sampled"].astype(bool)).sum())
        if thin:
            ax.text(0.98, 0.02,
                    f"{thin} hatched: fewer than five separate observations",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8, color="#666666")

    def default_xlabel(self) -> str | None:
        return "Fraction of frames present"

    def default_ylabel(self) -> str | None:
        return None


register_analysis(ProteinLigandInteractions.name, ProteinLigandInteractions)
