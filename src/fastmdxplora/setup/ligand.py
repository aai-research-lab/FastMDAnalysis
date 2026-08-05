"""Ligand / cofactor loading for protein-ligand systems.

This module validates and loads small-molecule ligands for parameterization
with the OpenFF small-molecule force fields (via ``openmmforcefields``'
``SystemGenerator``). It deliberately keeps the *loading/validation* concern
separate from the *system build* concern (which lives in
:mod:`fastmdxplora.setup.prepare`): here we turn a ligand file into a
validated OpenFF ``Molecule`` with a known net charge; the prepare step feeds
that molecule to the ``SystemGenerator``.

Supported input formats are SDF and MOL2 — the formats OpenFF reads cleanly
with full bond/charge information. (PDB ``HETATM`` extraction is intentionally
not supported yet; a bare PDB ligand lacks the bond orders OpenFF needs.)

The OpenFF toolkit is an optional dependency. :func:`load_ligand` raises a
clear, actionable :class:`LigandError` if it is not installed, rather than an
opaque ImportError, so the setup phase can degrade gracefully.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

from fastmdxplora.utils.logging import get_logger

if TYPE_CHECKING:  # pragma: no cover - typing only
    pass

logger = get_logger("setup.ligand")

#: Ligand file formats OpenFF reads with full chemical information.
SUPPORTED_LIGAND_FORMATS = ("sdf", "mol2")


class LigandError(Exception):
    """Raised for ligand input problems (format, missing file, charge, deps)."""


def detect_ligand_format(ligand_file: str | Path) -> str:
    """Return the ligand format (``sdf`` or ``mol2``) from the file suffix.

    Raises
    ------
    LigandError
        If the suffix is not a supported ligand format.
    """
    ext = Path(ligand_file).suffix.lower().lstrip(".")
    if ext not in SUPPORTED_LIGAND_FORMATS:
        supported = ", ".join(SUPPORTED_LIGAND_FORMATS)
        raise LigandError(
            f"Unsupported ligand format {ext!r} for {ligand_file!r}. "
            f"Use one of: {supported}. (A ligand embedded in a PDB lacks the "
            f"bond/charge information OpenFF needs; export it to SDF/MOL2.)"
        )
    return ext


def _import_openff() -> Any:
    """Import the OpenFF ``Molecule`` class, or raise a clear LigandError.

    Mirrors :func:`fastmdxplora.setup.prepare._import_openmm` so a missing
    optional dependency degrades with an actionable install message instead
    of an opaque ImportError.
    """
    try:
        from openff.toolkit import Molecule
    except ImportError as exc:  # pragma: no cover - exercised on hosts w/o openff
        raise LigandError(
            "Ligand parameterization needs the OpenFF toolkit, which is not "
            "installed and is not on PyPI: it is distributed through "
            "conda-forge only, so pip cannot fetch it whatever extra is "
            "named.\n\n"
            "    conda install -c conda-forge openff-toolkit "
            "openmmforcefields\n\n"
            "Or install FastMDXplora itself from conda-forge, which brings "
            "the whole ligand path with it."
        ) from exc
    return Molecule


def load_ligand(
    ligand_file: str | Path,
    *,
    name: str = "LIG",
    net_charge: int | None = None,
) -> Any:
    """Load and validate a ligand into an OpenFF ``Molecule``.

    Parameters
    ----------
    ligand_file : path
        Path to an SDF or MOL2 file.
    name : str, default "LIG"
        Residue/molecule name assigned to the ligand.
    net_charge : int, optional
        Formal net charge. If ``None`` (default), the charge is inferred from
        the molecule's formal charges (typical for a correctly prepared SDF).
        Supply this explicitly when the file is ambiguous or you need to
        override the inferred value.

    Returns
    -------
    openff.toolkit.Molecule
        The loaded molecule, with ``.name`` set.

    Raises
    ------
    LigandError
        On a missing file, unsupported format, missing OpenFF toolkit, or an
        unreadable/ambiguous ligand.
    """
    path = Path(ligand_file).expanduser().resolve()
    if not path.exists():
        raise LigandError(f"Ligand file not found: {path}")
    detect_ligand_format(path)

    Molecule = _import_openff()
    try:
        molecule = Molecule.from_file(str(path))
    except Exception as exc:  # noqa: BLE001 - normalize to LigandError
        raise LigandError(
            f"Could not read ligand {path.name!r} as a valid molecule: {exc}. "
            f"Ensure the SDF/MOL2 has explicit hydrogens and bond orders."
        ) from exc

    # Molecule.from_file may return a list when the file holds multiple
    # molecules; we parameterize a single ligand for now (the config is
    # list-shaped so multi-ligand support can layer on later).
    if isinstance(molecule, list):
        if len(molecule) != 1:
            raise LigandError(
                f"Ligand file {path.name!r} contains {len(molecule)} "
                f"molecules; provide a single-molecule SDF/MOL2 (multi-ligand "
                f"support is not yet implemented)."
            )
        molecule = molecule[0]

    molecule.name = name

    # Net charge: infer from formal charges unless the user overrides.
    inferred = _infer_net_charge(molecule)
    if net_charge is not None and inferred is not None and net_charge != inferred:
        logger.warning(
            "Ligand %s: user net_charge=%d overrides inferred charge=%d.",
            name, net_charge, inferred,
        )
    resolved_charge = net_charge if net_charge is not None else inferred
    logger.info(
        "Loaded ligand %s from %s (net charge=%s).",
        name, path.name,
        resolved_charge if resolved_charge is not None else "unknown",
    )
    return molecule


def _infer_net_charge(molecule: Any) -> int | None:
    """Best-effort integer net charge from an OpenFF molecule's total charge.

    Returns ``None`` if the charge cannot be determined as a clean integer.
    """
    try:
        total = molecule.total_charge
        # openff total_charge is a pint/openff Quantity in elementary charge;
        # coerce to a plain number robustly.
        magnitude = getattr(total, "magnitude", total)
        value = float(magnitude)
        rounded = round(value)
        if abs(value - rounded) > 1e-6:
            return None
        return int(rounded)
    except Exception:  # noqa: BLE001 - inference is best-effort
        return None
