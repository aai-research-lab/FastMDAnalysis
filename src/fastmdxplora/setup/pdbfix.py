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


def _drop_untemplated_gaps(fixer) -> None:
    """Cancel rebuilds of residues PDBFixer has no chemical template for.

    ``fixer.missingResidues`` maps an insertion point to the residue names to
    build there. A name with no template cannot be built; leaving it scheduled
    raises an ``AttributeError`` from deep inside ``addMissingAtoms`` that names
    nothing useful. Entries are dropped in place, so the surviving insertion
    indices stay consistent with the topology they were counted against.
    """
    missing = getattr(fixer, "missingResidues", None)
    if not missing:
        return

    # _getTemplate is PDBFixer's own lookup, and covers definitions downloaded
    # for non-standard components as well as the built-in set. Fall back to the
    # public dict if a future version drops the private helper.
    lookup = getattr(fixer, "_getTemplate", None)
    if lookup is None:
        lookup = getattr(fixer, "templates", {}).get

    pruned: dict = {}
    dropped: dict[str, int] = {}
    for insertion_point, names in missing.items():
        buildable = []
        for name in names:
            if lookup(name) is None:
                dropped[name] = dropped.get(name, 0) + 1
            else:
                buildable.append(name)
        if buildable:
            pruned[insertion_point] = buildable

    if not dropped:
        return

    summary = ", ".join(
        f"{name} ({count})" if count > 1 else name
        for name, count in sorted(dropped.items())
    )
    logger.info(
        "Not rebuilding %s: the deposited sequence lists them but no chemical "
        "template is available, which is expected for components removed as "
        "heterogens. The chain is prepared with those positions left as a gap.",
        summary,
    )
    fixer.missingResidues = pruned


def _drop_terminal_extensions(fixer, *, build_termini: bool = False) -> None:
    """Build gaps between resolved residues; do not extend a chain past its end.

    ``findMissingResidues`` compares SEQRES against what was modelled, so every
    unresolved residue is scheduled -- including the ones before the first and
    after the last that anyone could see. Those two cases are not the same
    thing.

    A gap between two resolved residues is a loop: both ends are pinned, and
    what is built has to reach from one to the other. A run past the last
    resolved residue is a disordered terminus, anchored at one end and free at
    the other, and what gets built walks wherever the builder puts it.

    6B73 is the case. SEQRES declares about 429 residues per receptor copy and
    281 and 269 were modelled, so 137 and 149 were scheduled -- mostly
    terminal. The built chains reached 54 and 59 nm from the structure, giving
    a "protein" 51 nm across for 1,104 residues, and `addMembrane` then spent
    ten minutes packing lipids across that face before failing with a NaN that
    named none of this. Its two nanobody chains, whose gaps are internal, were
    built correctly and stayed within 7 nm.

    Residues nobody observed are not restored by guessing where they went. The
    ones dropped here are named, because a methods section should say the
    termini were not modelled.
    """
    missing = getattr(fixer, "missingResidues", None)
    if not missing or build_termini:
        return

    try:
        lengths = {
            index: sum(1 for _ in chain.residues())
            for index, chain in enumerate(fixer.topology.chains())
        }
    except (AttributeError, TypeError):
        return

    kept: dict = {}
    dropped: dict[int, int] = {}
    for insertion_point, names in missing.items():
        try:
            chain_index, offset = insertion_point
        except (TypeError, ValueError):
            kept[insertion_point] = names
            continue
        length = lengths.get(chain_index)
        # Offset 0 inserts before the first resolved residue; offset == length
        # appends after the last. Everything between is a gap with an anchor
        # at both ends.
        if length is not None and 0 < offset < length:
            kept[insertion_point] = names
        else:
            dropped[chain_index] = dropped.get(chain_index, 0) + len(names)

    if not dropped:
        return

    summary = ", ".join(
        f"chain {index}: {count}" for index, count in sorted(dropped.items())
    )
    logger.info(
        "Not building %d unresolved residue(s) at chain termini (%s). They "
        "were never observed, and a terminus has an anchor at one end only, "
        "so what is built there is placed rather than determined. Internal "
        "gaps, which are pinned at both ends, are built as usual. Pass "
        "setup.build_missing_termini: true to build them anyway.",
        sum(dropped.values()), summary,
    )
    fixer.missingResidues = kept


def fix_pdb_with_pdbfixer(
    input_pdb: str,
    output_pdb: str,
    *,
    ph: float = 7.0,
    build_missing_termini: bool = False,
    keep_heterogens: bool = False,
    keep_water: bool = False,
    reinstated: tuple[str, ...] = (),
    explained: tuple[str, ...] = (),
    replace_nonstandard: bool = True,
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

    # removeHeterogens() deletes components from the topology but not from
    # SEQRES, so findMissingResidues() reads the hole they leave as unresolved
    # polymer and schedules it for rebuilding. PDBFixer holds no template for a
    # component it has just discarded, and dies inside addMissingAtoms with
    # "'NoneType' object has no attribute 'topology'" -- a message that names
    # neither the residue nor the cause. 6LU7 reaches it, because its inhibitor
    # chain interleaves standard residues with 010, 02J and PJE.
    #
    # The scheduled rebuilds are pruned rather than the calls reordered.
    # Calling findMissingResidues() before removeHeterogens() also clears the
    # crash, but its keys are (chain index, residue index) counted on the
    # topology that still holds the heterogens; removing them afterwards shifts
    # every later index down by one, and a genuinely unresolved loop is then
    # rebuilt at the wrong place -- appended to the chain terminus instead of
    # inserted into its gap. That failure is silent, and a quietly wrong
    # structure is worse than a loud crash.
    _drop_untemplated_gaps(fixer)
    _drop_terminal_extensions(fixer, build_termini=build_missing_termini)

    # Modified residues are part of the polymer, not ligands: a selenomethionine
    # or an oxidised cysteine belongs in the chain. Left in place they reach the
    # heterogen classifier, which can only refuse them, and a structure with a
    # single modified cysteine becomes unusable. PDBFixer substitutes the
    # standard equivalent, which is the ordinary preparation choice and what
    # every comparable tool does.
    if replace_nonstandard:
        fixer.findNonstandardResidues()
        substitutions = list(getattr(fixer, "nonstandardResidues", []) or [])
        if substitutions:
            described = ", ".join(
                f"{residue.name}{getattr(residue, 'id', '')}->{standard}"
                for residue, standard in substitutions
            )
            logger.info(
                "Replaced %d modified residue(s) with their standard "
                "equivalents: %s", len(substitutions), described,
            )
            fixer.replaceNonstandardResidues()

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=float(ph))

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    logger.info(" - wrote fixed PDB to %s", out)
