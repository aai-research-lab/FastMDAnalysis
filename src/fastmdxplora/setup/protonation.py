"""Decide a ligand's protonation from its pKa *in the complex*.

The usual approach protonates a ligand in isolation: a rule table says a
carboxylic acid is deprotonated above pH 4, so at pH 7 it is a carboxylate.
That is a statement about the molecule in water, and it is often wrong about
the molecule in a pocket. A carboxylate buried against an aspartate, or in a
hydrophobic cavity with no compensating charge, can easily remain neutral;
shifts of two or three pKa units on binding are ordinary rather than exotic.

So the question "what charge state does this ligand adopt?" cannot be answered
by looking at the ligand. It is a property of the complex. PROPKA computes
exactly that, for proteins since 3.0 and for protein-ligand complexes since
3.1, so the honest thing is to ask it about the bound state.

Three outcomes, and only one of them proceeds:

* PROPKA returns a pKa comfortably away from the simulation pH. The state is
  determined; adopt it and say so.
* PROPKA returns a pKa close to the pH. The group is genuinely poised, and
  both states are populated. Refuse.
* PROPKA returns nothing usable for the ligand. Refuse, rather than falling
  back to a rule table that answers a different question.

The third case is not a failure of the design. Empirical pKa prediction on an
arbitrary ligand is not always possible, and a tool that quietly substitutes a
worse method when the better one declines has learned nothing.
"""

from __future__ import annotations

import contextlib
import io
from dataclasses import dataclass
from pathlib import Path

from fastmdxplora.utils.logging import get_logger

logger = get_logger("setup.protonation")

#: A group whose pKa is within this many units of the simulation pH is
#: appreciably populated in both states: one unit corresponds to roughly a
#: 9:1 ratio. Choosing one state would misrepresent the ensemble.
POISED_MARGIN = 1.0


class ProtonationError(RuntimeError):
    """The ligand's protonation at the requested pH could not be settled."""


@dataclass(frozen=True)
class GroupPka:
    """One ionizable group and the pKa it has in this complex."""

    resname: str
    chain: str
    resseq: int
    group_type: str
    pka: float
    model_pka: float

    @property
    def shift(self) -> float:
        """How far the environment moved this group from its solution value."""
        return self.pka - self.model_pka

    def __str__(self) -> str:
        return (
            f"{self.resname} {self.chain}{self.resseq} ({self.group_type}): "
            f"pKa {self.pka:.2f}, shifted {self.shift:+.2f} from {self.model_pka:.2f}"
        )


@dataclass(frozen=True)
class ProtonationState:
    """A settled protonation, with the evidence that settled it."""

    resname: str
    protonated: bool
    groups: tuple[GroupPka, ...]
    reason: str


def _propka():
    try:
        import propka.run

        return propka.run
    except ImportError as exc:  # pragma: no cover - exercised by absence
        raise ProtonationError(
            "PROPKA is needed to determine a ligand's protonation in the "
            "complex and is not installed. Install it (conda install -c "
            "conda-forge propka), or supply the ligand with the protonation "
            "you intend already assigned."
        ) from exc


def ligand_pka(
    complex_pdb: str | Path,
    resname: str,
    *,
    chain: str,
    resseq: int,
) -> list[GroupPka]:
    """Ionizable groups of one ligand copy, with their pKa in this complex.

    Returns an empty list when PROPKA finds no ionizable group for the
    component, which for a ligand like benzene is the correct answer.
    """
    run = _propka()
    resname = resname.upper()

    # PROPKA writes progress to stdout; keep it out of the phase output.
    buffer = io.StringIO()
    try:
        with contextlib.redirect_stdout(buffer):
            molecule = run.single(str(complex_pdb), optargs=("--quiet",), write_pka=False)
    except Exception as exc:  # noqa: BLE001 - any failure means "cannot settle"
        raise ProtonationError(
            f"PROPKA could not analyse the complex containing {resname}: {exc}"
        ) from exc

    if not molecule.conformation_names:
        raise ProtonationError("PROPKA returned no conformation to read.")
    conformation = molecule.conformations[molecule.conformation_names[0]]

    found: list[GroupPka] = []
    for group in conformation.groups:
        atom = group.atom
        if atom.res_name.strip().upper() != resname:
            continue
        if str(atom.chain_id).strip() != str(chain).strip():
            continue
        if int(atom.res_num) != int(resseq):
            continue
        # Backbone pseudo-groups carry a zero model pKa and are not titratable
        # in the sense meant here.
        if group.model_pka == 0.0 and group.pka_value == 0.0:
            continue
        found.append(
            GroupPka(
                resname=resname,
                chain=str(atom.chain_id).strip(),
                resseq=int(atom.res_num),
                group_type=str(group.type),
                pka=float(group.pka_value),
                model_pka=float(group.model_pka),
            )
        )
    return found


def decide(
    resname: str,
    groups: list[GroupPka],
    ph: float,
    *,
    expected_ionizable: bool,
    margin: float = POISED_MARGIN,
) -> ProtonationState:
    """Settle a ligand's protonation, or refuse.

    ``expected_ionizable`` comes from inspecting the ligand's chemistry: if it
    carries an ionizable group but PROPKA reported none, the two disagree and
    the disagreement is itself a reason not to proceed.
    """
    if not groups:
        if expected_ionizable:
            raise ProtonationError(
                f"{resname} carries an ionizable group, but PROPKA returned no "
                "pKa for it in this complex, so its charge state here is "
                "undetermined. Empirical prediction does not succeed on every "
                "ligand.\nSupply the ligand with the protonation you intend, as "
                "an SDF or MOL2 file, and set its net charge."
            )
        return ProtonationState(
            resname=resname,
            protonated=False,
            groups=(),
            reason="no ionizable group; protonation is unambiguous",
        )

    poised = [g for g in groups if abs(g.pka - ph) < margin]
    if poised:
        detail = "\n".join(f"  {g}" for g in poised)
        raise ProtonationError(
            f"{resname} has ionizable groups whose pKa sits within {margin:g} "
            f"unit of pH {ph:g}, so both charge states are appreciably "
            f"populated and neither represents the ensemble:\n{detail}\n"
            "Choose the state you intend and supply the ligand explicitly, or "
            "run at a pH where the group is not poised."
        )

    # Every group is decisively on one side of the pH.
    protonated_groups = [g for g in groups if g.pka > ph]
    detail = "; ".join(str(g) for g in groups)
    state = ProtonationState(
        resname=resname,
        protonated=bool(protonated_groups),
        groups=tuple(groups),
        reason=(
            f"pKa determined in the complex at pH {ph:g}: {detail}"
        ),
    )
    logger.info(
        "%s protonation settled from the complex: %d group(s), %s",
        resname, len(groups),
        "protonated" if state.protonated else "deprotonated",
    )
    return state


def settle(
    complex_pdb: str | Path,
    resname: str,
    *,
    chain: str,
    resseq: int,
    ph: float,
    expected_ionizable: bool,
    margin: float = POISED_MARGIN,
) -> ProtonationState:
    """Determine a ligand's protonation in the complex, or raise."""
    groups = ligand_pka(complex_pdb, resname, chain=chain, resseq=resseq)
    return decide(resname, groups, ph, expected_ionizable=expected_ionizable,
                  margin=margin)
