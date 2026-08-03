"""PDBFixer wrapper.

PDBFixer wrapper providing
a single function that takes a raw PDB and produces a clean,
fully-protonated PDB suitable for solvation and parameterization.

The implementation calls PDBFixer's standard sequence:

  1. ``removeHeterogens`` (optional; removes everything that is not a
     standard residue, with an option to retain crystallographic waters)
  2. ``findMissingResidues`` (chain-break / loop detection)
  3. ``findMissingAtoms`` (heavy-atom completion)
  4. ``addMissingAtoms``
  5. ``addMissingHydrogens(pH)``  — places hydrogens at the specified pH

The function is strict: it raises rather than returning an error code.
This makes phase-level error
handling easy (the orchestrator's per-phase try/except records the
failure cleanly).

Requires :mod:`pdbfixer` and :mod:`openmm`, both conda-forge packages
in the optional ``[setup]`` extras group.
"""

from __future__ import annotations

from pathlib import Path

from fastmdxplora.utils.logging import get_logger

logger = get_logger("setup.pdbfix")



# Water and monatomic ions are removed as a matter of course and are not worth
# reporting; anything else might be the ligand.
_UNREMARKABLE_RESIDUES = frozenset({
    "HOH", "WAT", "H2O", "TIP", "TIP3", "SOL", "DOD",
    "NA", "K", "CL", "MG", "CA", "ZN", "MN", "FE", "CU", "BR", "IOD", "SO4",
})


def _heterogen_residue_counts(topology) -> dict[str, int]:
    """Count non-water, non-ion heterogen residues about to be discarded."""
    counts: dict[str, int] = {}
    standard = _standard_residue_names()
    for residue in topology.residues():
        name = (residue.name or "").strip().upper()
        if not name or name in standard or name in _UNREMARKABLE_RESIDUES:
            continue
        counts[name] = counts.get(name, 0) + 1
    return counts


def _standard_residue_names() -> frozenset[str]:
    """Amino acids and nucleotides, which are not heterogens."""
    amino = (
        "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
        "TRP TYR VAL HID HIE HIP HSD HSE HSP CYX ASH GLH LYN MSE"
    ).split()
    nucleic = (
        "A C G U T DA DC DG DT DU RA RC RG RU "
        "A3 A5 C3 C5 G3 G5 U3 U5 DA3 DA5 DC3 DC5 DG3 DG5 DT3 DT5"
    ).split()
    return frozenset(amino + nucleic)


def fix_pdb_with_pdbfixer(
    input_pdb: str,
    output_pdb: str,
    *,
    ph: float = 7.0,
    keep_heterogens: bool = False,
    keep_water: bool = False,
    reinstated: tuple[str, ...] = (),
    explained: tuple[str, ...] = (),
) -> None:
    """Strict PDBFixer wrapper: raises on failure.

    Parameters
    ----------
    input_pdb : path-like
        Input PDB file path.
    output_pdb : path-like
        Where to write the fixed PDB. Parent directories are created.
    ph : float, default 7.0
        pH for hydrogen placement. Determines protonation state of
        titratable residues (Asp, Glu, His, Lys, etc.) via PDBFixer's
        residue-template library.
    keep_heterogens : bool, default False
        If True, retain non-standard residues (ligands, cofactors,
        ions). Default removes them.
    keep_water : bool, default False
        If True (and ``keep_heterogens=False``), retain crystallographic
        waters during heterogen removal. Has no effect when
        ``keep_heterogens=True``.

    Raises
    ------
    ImportError
        If ``pdbfixer`` or ``openmm`` is not installed. Install the
        ``[md]`` extras: ``pip install fastmdxplora[md]``, or
        better via conda: ``conda install -c conda-forge pdbfixer openmm``.
    FileNotFoundError
        If ``input_pdb`` doesn't exist.
    Exception
        Re-raises any error from PDBFixer (residue identification
        failures, malformed input, etc.).

    Notes
    -----
    Provides a clean PDBFixer wrapper with the
    ``fix_pdb_with_pdbfixer`` exactly so users moving between the two
    tools see identical results.
    """
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise ImportError(
            "fix_pdb_with_pdbfixer requires pdbfixer and openmm. Install "
            "via conda (recommended): conda install -c conda-forge "
            "pdbfixer openmm — or via pip with the optional [setup] "
            "extras: pip install fastmdxplora[md]."
        ) from exc

    inp = Path(input_pdb)
    out = Path(output_pdb)

    if not inp.exists():
        raise FileNotFoundError(f"Input PDB not found: {inp}")

    logger.info("Fixing PDB with PDBFixer: %s (pH=%s)", inp, ph)

    fixer = PDBFixer(filename=str(inp))
    if not keep_heterogens:
        removed = _heterogen_residue_counts(fixer.topology)
        fixer.removeHeterogens(keepWater=keep_water)
        kept = {name.upper() for name in reinstated}
        # Components the caller has already reasoned about and reported on.
        # Warning that a buffer molecule "was removed, and might be the ligand
        # you meant to simulate" contradicts a decision just explained to the
        # user, and invites them to second-guess a correct one.
        accounted = kept | {name.upper() for name in explained}
        discarded = {n: c for n, c in removed.items() if n.upper() not in accounted}

        if kept:
            # Under the auto policy the components worth simulating have
            # already been parameterized and are added back through the
            # small-molecule path. Warning that they were "removed" would
            # describe a loss that did not happen.
            logger.info(
                "Stripped %s from the structure; they are re-added with "
                "small-molecule parameters.",
                ", ".join(sorted(kept)),
            )

        if discarded:
            summary = ", ".join(
                f"{name} ({count})" if count > 1 else name
                for name, count in sorted(discarded.items())
            )
            # Removing crystallization additives is usually right, but a bound
            # ligand looks identical to a buffer molecule at this stage. Say
            # what went, so a silently apo run cannot be mistaken for a holo
            # one.
            logger.warning(
                "Removed heterogens: %s. Pass --setup-keep-heterogens (or "
                "setup.keep_heterogens: true) to retain them, for example when "
                "one of these is the ligand you intend to simulate.",
                summary,
            )
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=float(ph))

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    logger.info(" - wrote fixed PDB to %s", out)
