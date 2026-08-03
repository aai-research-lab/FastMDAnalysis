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
import logging
from dataclasses import dataclass
from pathlib import Path

from fastmdxplora.utils.logging import get_logger

logger = get_logger("setup.protonation")

#: A group whose pKa is within this many units of the simulation pH is
#: appreciably populated in both states: one unit corresponds to roughly a
#: 9:1 ratio. Choosing one state would misrepresent the ensemble.
POISED_MARGIN = 1.0



class _CapturingHandler(logging.Handler):
    """Collects PROPKA's own log records instead of letting them print."""

    def __init__(self) -> None:
        super().__init__()
        self.messages: list[str] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.messages.append(record.getMessage())


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


#: Whether each titratable pattern, when it matches, matched the protonated
#: form. The SMARTS already carry this: a carboxylic acid pattern matches only
#: the neutral acid and a carboxylate only the anion, and the amine patterns
#: exclude ``[N+]``, so they match only the neutral base. The phosphate and
#: sulfonate patterns accept either form and so cannot report which they found;
#: they are absent here and no claim is made about them.
WRITTEN_PROTONATED: dict[str, bool] = {
    "carboxylic acid": True,
    "carboxylate": False,
    "amidine": False,
    "primary or secondary amine": False,
    "tertiary amine": False,
    "imidazole": False,
    "guanidine": False,
    "thiol": True,
    "tetrazole": True,
}


#: Where the proton goes for each titratable group, as (SMARTS, index of the
#: atom within the match). The site is not always the atom the classifying
#: pattern is named for: an amidine and a guanidine protonate on the sp2
#: nitrogen, because the cation is delocalised. Building it on the sp3 nitrogen
#: gives a different molecule which also keeps the C=N double bond, and with it
#: an undefined stereocentre the small-molecule toolkits refuse.
#:
#: Groups absent here have no unambiguous site. The phosphate and sulfonate
#: patterns match either form and name no particular oxygen, so a settled state
#: cannot be placed and the ligand is refused instead.
PROTONATION_SITE: dict[str, tuple[str, int]] = {
    "amidine": ("[NX3][CX3]=[NX2]", 2),
    "guanidine": ("[NX3][CX3](=[NX2])[NX3]", 2),
    "primary or secondary amine": ("[NX3;H2,H1;!$(NC=O);!$(N[a])]", 0),
    "tertiary amine": ("[NX3;H0;!$(NC=O);!$(N[a]);!$([N+])]", 0),
    "imidazole": ("c1cnc[nH]1", 2),
    "carboxylic acid": ("[CX3](=O)[OX2H1]", 2),
    "thiol": ("[SX2H1]", 0),
    "tetrazole": ("c1nnn[nH]1", 1),
}


def _rdkit():
    try:
        from rdkit import Chem

        return Chem
    except ImportError as exc:  # pragma: no cover - exercised by absence
        raise ProtonationError(
            "RDKit is needed to put a ligand into the protonation state settled "
            "in the complex, and is not installed. Install it (conda install -c "
            "conda-forge rdkit), or supply the ligand with --setup-ligand "
            "already in the state you intend."
        ) from exc


def _add_proton(Chem, mol, index):
    """Attach a proton, giving the new hydrogen a position."""
    editable = Chem.RWMol(mol)
    atom = editable.GetAtomWithIdx(index)
    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
    atom.SetNoImplicit(False)
    built = editable.GetMol()
    Chem.SanitizeMol(built)
    return Chem.AddHs(built, addCoords=True, onlyOnAtoms=[index])


def _remove_proton(Chem, mol, index):
    """Strip a proton, preferring an explicit hydrogen if one is present."""
    editable = Chem.RWMol(mol)
    atom = editable.GetAtomWithIdx(index)
    atom.SetFormalCharge(atom.GetFormalCharge() - 1)
    hydrogen = next(
        (n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 1), None
    )
    if hydrogen is not None:
        editable.RemoveAtom(hydrogen)
    else:
        atom.SetNumExplicitHs(max(0, atom.GetNumExplicitHs() - 1))
    built = editable.GetMol()
    Chem.SanitizeMol(built)
    return built


def apply_settled_state(sdf_text: str, chemistry, state) -> tuple[str, int]:
    """Put the reference chemistry into the state settled in the complex.

    Returns the rewritten SDF and its net formal charge.

    The pKa is determined in the pocket, but the file handed to the
    small-molecule force field is the reference chemistry, which carries
    whatever protonation the dictionary happens to hold. Parameterizing that
    file simulates the ligand in a charge state the pocket contradicts, and the
    charge is usually what binds it: the amidinium of benzamidine, the
    carboxylate of retinoic acid, the ammonium of a tamoxifen.

    Raises rather than guessing where the settled state cannot be placed
    unambiguously -- an unrecognised group, or more than one titratable group,
    which a single protonated/deprotonated answer cannot describe.
    """
    groups = tuple(chemistry.titratable_groups)
    if not groups:
        return sdf_text, chemistry.formal_charge

    if len(groups) > 1:
        raise ProtonationError(
            f"{chemistry.resname} carries more than one titratable group "
            f"({', '.join(groups)}), and a single settled answer cannot say "
            "which is protonated and which is not: a molecule with both a basic "
            "amine and an acid is a zwitterion, not uniformly one or the other."
            "\nSupply the ligand with --setup-ligand as an SDF or MOL2 already "
            "in the state you intend, with explicit hydrogens and its net "
            "charge set."
        )

    label = groups[0]
    site = PROTONATION_SITE.get(label)
    if site is None:
        raise ProtonationError(
            f"{chemistry.resname} carries a {label}, whose protonation site is "
            "not one this pipeline knows how to place. Supply the ligand with "
            "--setup-ligand already in the state you intend."
        )

    smarts, offset = site
    if WRITTEN_PROTONATED.get(label) == state.protonated:
        # The reference chemistry is already in the settled state.
        return sdf_text, chemistry.formal_charge

    Chem = _rdkit()
    mol = Chem.MolFromMolBlock(sdf_text, removeHs=False, sanitize=True)
    if mol is None:
        raise ProtonationError(
            f"The reference chemistry for {chemistry.resname} could not be read "
            "as a molecule, so its protonation cannot be adjusted."
        )

    matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    if len(matches) != 1:
        raise ProtonationError(
            f"{chemistry.resname} has {len(matches)} matches for its {label}, "
            "so the settled protonation cannot be placed unambiguously. Supply "
            "the ligand with --setup-ligand already in the state you intend."
        )

    index = matches[0][offset]
    change = _add_proton if state.protonated else _remove_proton
    adjusted = change(Chem, mol, index)

    charge = Chem.GetFormalCharge(adjusted)
    expected = chemistry.formal_charge + (1 if state.protonated else -1)
    if charge != expected:
        raise ProtonationError(
            f"Adjusting {chemistry.resname}'s {label} produced a net charge of "
            f"{charge:+d} where {expected:+d} was intended, so the reference "
            "chemistry is not what this pipeline assumed. Supply the ligand "
            "with --setup-ligand already in the state you intend."
        )
    logger.info(
        "%s: %s its %s to match the state settled in the complex; net charge "
        "%+d.",
        chemistry.resname,
        "protonated" if state.protonated else "deprotonated",
        label, charge,
    )
    return Chem.MolToMolBlock(adjusted), charge


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

    # PROPKA reports progress on stdout and structural complaints through
    # logging. Both are captured: an unmodelled surface side chain is worth
    # knowing about, but it belongs in the debug log rather than interleaved
    # with the phase output.
    buffer = io.StringIO()
    propka_logger = logging.getLogger("propka")
    captured = _CapturingHandler()
    previous_level = propka_logger.level
    propka_logger.addHandler(captured)
    propka_logger.setLevel(logging.WARNING)
    try:
        with contextlib.redirect_stdout(buffer):
            molecule = run.single(str(complex_pdb), optargs=("--quiet",), write_pka=False)
    except Exception as exc:  # noqa: BLE001 - any failure means "cannot settle"
        raise ProtonationError(
            f"PROPKA could not analyse the complex containing {resname}: {exc}"
        ) from exc
    finally:
        propka_logger.removeHandler(captured)
        propka_logger.setLevel(previous_level)

    incomplete = [m for m in captured.messages if "Unexpected number" in m]
    if incomplete:
        # Residues the depositors could not model. Surface side chains with
        # poor density are ordinary; near the binding site they would matter.
        logger.debug("PROPKA noted %d incomplete residues:\n  %s",
                     len(incomplete), "\n  ".join(incomplete))
        logger.info(
            "%d residues in this structure are missing side-chain atoms; "
            "pKa values are computed from what was modelled.", len(incomplete),
        )

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
    """Determine a ligand's protonation in the complex, or raise.

    A ligand with no ionizable group has one protonation state whatever its
    surroundings, so the calculation is skipped: it could not change the
    answer, and running it only invites complaints about a structure whose
    imperfections are irrelevant here.
    """
    if not expected_ionizable:
        return ProtonationState(
            resname=resname,
            protonated=False,
            groups=(),
            reason="no ionizable group; protonation is unambiguous",
        )
    groups = ligand_pka(complex_pdb, resname, chain=chain, resseq=resseq)
    return decide(resname, groups, ph, expected_ionizable=expected_ionizable,
                  margin=margin)
