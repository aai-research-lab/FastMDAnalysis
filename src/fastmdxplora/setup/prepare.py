"""System preparation: solvate, ionize, parameterize.

Takes a fully-protonated PDB (typically the output of PDBFixer) and
produces a simulation-ready OpenMM ``System`` together with a starting
``State`` and a topology PDB. The outputs are serialized to disk so the
simulation phase can load them without re-doing any of this work.

Pipeline
--------
  1. Load the topology + positions from the prepared PDB.
  2. Build an OpenMM ``ForceField`` from the force-field XMLs.
  3. Use ``Modeller`` to add a water box (with configurable padding) and
     ions for charge neutralization + a target ionic concentration.
  4. Call ``ForceField.createSystem(...)`` to apply parameters and obtain
     an OpenMM ``System`` with all the standard constraints/options.
  5. Build an OpenMM ``Context`` to capture the initial ``State``
     (positions, velocities, box vectors).
  6. Serialize ``System`` and ``State`` to XML and write a topology PDB
     of the solvated system.

Default force fields
--------------------
``charmm36.xml`` + ``charmm36/water.xml`` for protein-only systems —
uses sensible general-purpose defaults giving consistent
parameterization out of the box. Override via the ``force_field`` and
``water_model`` options.

Outputs
-------
  - ``solvated.pdb`` -- topology + positions after solvation
  - ``topology.pdb`` -- alias of ``solvated.pdb`` (the canonical
    "topology" the analysis phase consumes)
  - ``system.xml`` -- serialized OpenMM ``System`` (force field applied)
  - ``state.xml`` -- serialized initial ``State``
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from fastmdxplora.setup.forcefields import (
    available_forcefields,
    require_forcefield_xmls,
    resolve_forcefield,
)
from fastmdxplora.setup.ligand import load_ligand
from fastmdxplora.utils.logging import get_logger

logger = get_logger("setup.prepare")


# ---------------------------------------------------------------------------
# Defaults — sensible general-purpose biomolecular MD settings
# (The default force field now lives in setup/forcefields.py, the single
# source of truth for the named force-field registry.)
# ---------------------------------------------------------------------------
DEFAULT_PROBE_RADIUS_NM = 0.14
DEFAULT_PADDING_NM = 1.0
DEFAULT_IONIC_STRENGTH_M = 0.15
DEFAULT_PH = 7.0


def _import_openmm():
    """Lazy import of OpenMM. Raises a clean ImportError otherwise."""
    try:
        import openmm
        from openmm import unit
        from openmm.app import (
            CutoffNonPeriodic,
            CutoffPeriodic,
            Ewald,
            ForceField,
            HBonds,
            Modeller,
            NoCutoff,
            PDBFile,
            PME,
        )

        return {
            "openmm": openmm,
            "unit": unit,
            "ForceField": ForceField,
            "HBonds": HBonds,
            "Modeller": Modeller,
            "PDBFile": PDBFile,
            "PME": PME,
            "NoCutoff": NoCutoff,
            "CutoffNonPeriodic": CutoffNonPeriodic,
            "CutoffPeriodic": CutoffPeriodic,
            "Ewald": Ewald,
        }
    except ImportError as exc:
        raise ImportError(
            "Setup-phase system parameterization requires OpenMM. Install "
            "via conda (recommended): conda install -c conda-forge openmm "
            "pdbfixer — or via pip with the optional [md] extras: "
            "pip install fastmdxplora[md]."
        ) from exc


def prepare_system(
    prepared_pdb: str | Path,
    output_dir: str | Path,
    *,
    forcefield: str | None = None,
    force_field: list[str] | None = None,
    water_model: str | None = None,
    ligand: str | Path | list[str | Path] | None = None,
    ligand_forcefield: str | None = None,
    ligand_name: str = "LIG",
    ligand_net_charge: int | None = None,
    check_ligand_clashes: bool = True,
    ligand_clash_threshold_nm: float = 0.15,
    solvent_padding_nm: float = DEFAULT_PADDING_NM,
    box_shape: str = "cube",
    ion_positive: str = "Na+",
    ion_negative: str = "Cl-",
    ion_concentration_M: float = DEFAULT_IONIC_STRENGTH_M,
    neutralize: bool = True,
    nonbonded_method: str = "PME",
    nonbonded_cutoff_nm: float = 1.0,
    ewald_error_tolerance: float = 0.0005,
    use_switching_function: bool = True,
    switch_distance_nm: float | None = None,
    dispersion_correction: bool = True,
    remove_cm_motion: bool = False,
    constraints: str = "HBonds",
    rigid_water: bool = True,
    hydrogen_mass_amu: float | None = None,
    temperature_K: float = 300.0,
) -> dict[str, Path]:
    """Solvate, ionize, parameterize, and serialize an OpenMM system.

    Parameters
    ----------
    prepared_pdb : path
        Input PDB (typically the output of :func:`fix_pdb_with_pdbfixer`).
    output_dir : path
        Where to write ``solvated.pdb``, ``topology.pdb``, ``system.xml``,
        and ``state.xml``. Parent directories are created.
    force_field : list of str, optional
        Force-field XML file names recognized by OpenMM. Default is
        ``["charmm36.xml", "charmm36/water.xml"]``.
        Pass the protein force field plus its accompanying water model.
    water_model : str, optional
        Water model name (``"tip3p"``, ``"tip4pew"``, ``"spce"``, etc.)
        used by Modeller. When ``None`` (the default), Modeller picks the
        model that matches the supplied water-model XML.
    solvent_padding_nm : float, default 1.0
        Minimum distance in nm between any solute atom and the periodic
        box wall.
    box_shape : {"cube", "dodecahedron", "octahedron"}, default "cube"
        Periodic box geometry.
    ion_positive, ion_negative : str
        Counter-ions. Defaults to NaCl.
    ion_concentration_M : float, default 0.15
        Target ionic concentration in M (physiological for biomolecules).
    neutralize : bool, default True
        Add ions to neutralize the net solute charge before reaching the
        target concentration.
    nonbonded_cutoff_nm : float, default 1.0
        PME real-space cutoff in nm.
    constraints : {"None", "HBonds", "AllBonds", "HAngles"}, default "HBonds"
        Bond constraints. ``HBonds`` is the standard choice for 2 fs
        timesteps with water.
    rigid_water : bool, default True
        Constrain water bond lengths and angles.
    hydrogen_mass_amu : float, optional
        If set, repartition heavy-atom mass onto hydrogens to enable
        longer timesteps (the standard "HMR" technique, 4 amu allows
        ~4 fs steps with constraints).
    temperature_K : float, default 300.0
        Used to set initial velocities on the State (Maxwell-Boltzmann).

    Returns
    -------
    dict
        Mapping artifact-name -> ``Path``: ``solvated_pdb``,
        ``topology_pdb``, ``system_xml``, ``state_xml``.
    """
    omm = _import_openmm()
    unit = omm["unit"]

    prepared_path = Path(prepared_pdb)
    if not prepared_path.exists():
        raise FileNotFoundError(f"Prepared PDB not found: {prepared_path}")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve the force field. Precedence: an explicit raw `force_field`
    # XML list (power-user escape hatch) wins; otherwise resolve the named
    # `forcefield` selector to its XML set and default water model. A named
    # choice and a raw list together is rejected upstream (setup pipeline).
    ff_choice = None
    if force_field:
        force_field = list(force_field)
    else:
        ff_choice = resolve_forcefield(forcefield)
        require_forcefield_xmls(ff_choice)
        force_field = list(ff_choice.xmls)
        if water_model is None:
            water_model = ff_choice.water_model

    # Normalize the ligand argument to a list (the config is list-shaped from
    # day one; only single-ligand is implemented for now).
    ligands = _normalize_ligands(ligand)
    if ligands:
        # A ligand requires a ligand-capable force field. Raw XML lists can't
        # be introspected for ligand support, so require the named selector.
        if ff_choice is None:
            raise ValueError(
                "A ligand was supplied with a raw `force_field` XML list. "
                "Protein-ligand parameterization needs the named "
                "`forcefield` selector (use forcefield='amber-openff'), "
                "which wires the OpenFF small-molecule generator."
            )
        if not ff_choice.supports_ligand:
            valid = ", ".join(
                n for n in available_forcefields()
                if resolve_forcefield(n).supports_ligand
            )
            raise ValueError(
                f"Force field {ff_choice.name!r} does not support ligands. "
                f"For protein-ligand systems use a ligand-capable force "
                f"field: {valid}."
            )
        if len(ligands) > 1:
            raise ValueError(
                f"{len(ligands)} ligands supplied; only single-ligand "
                f"parameterization is implemented currently."
            )

    constraints_obj = _resolve_constraints(omm, constraints)

    # ----- 1. Load topology + positions -----
    logger.info("Loading prepared PDB: %s", prepared_path)
    pdb = omm["PDBFile"](str(prepared_path))

    # ----- 2. Build force field (+ ligand) -----
    modeller = omm["Modeller"](pdb.topology, pdb.positions)
    system_generator = None
    if ligands:
        # Protein-ligand: load the ligand as an OpenFF Molecule, build a
        # SystemGenerator that combines the protein/water force field with
        # the OpenFF small-molecule force field, and add the ligand into the
        # topology so it is solvated together with the protein.
        sm_ff = ligand_forcefield or ff_choice.small_molecule_forcefield
        logger.info(
            "Building protein-ligand ForceField: %s + small-molecule %s",
            force_field, sm_ff,
        )
        ligand_mol = load_ligand(
            ligands[0], name=ligand_name, net_charge=ligand_net_charge,
        )
        ff, system_generator = _build_ligand_forcefield(
            force_field, sm_ff, ligand_mol,
        )
        n_protein_atoms = modeller.topology.getNumAtoms()
        _add_ligand_to_modeller(omm, modeller, ligand_mol, ligand_name=ligand_name)
        if check_ligand_clashes:
            _check_ligand_clashes(
                modeller, n_protein_atoms, unit,
                threshold_nm=ligand_clash_threshold_nm,
                ligand_name=ligand_name,
            )
    else:
        logger.info("Building ForceField: %s", force_field)
        ff = omm["ForceField"](*force_field)

    # ----- 3. Solvate + ionize with Modeller -----
    logger.info(
        "Solvating (padding=%.2f nm, box=%s, ions=%s/%s @ %.3f M)",
        solvent_padding_nm,
        box_shape,
        ion_positive,
        ion_negative,
        ion_concentration_M,
    )

    add_solvent_kwargs: dict[str, Any] = {
        "padding": solvent_padding_nm * unit.nanometer,
        "boxShape": box_shape,
        "positiveIon": ion_positive,
        "negativeIon": ion_negative,
        "ionicStrength": ion_concentration_M * unit.molar,
        "neutralize": neutralize,
    }
    if water_model is not None:
        add_solvent_kwargs["model"] = water_model

    try:
        modeller.addSolvent(ff, **add_solvent_kwargs)
    except TypeError:
        # Older OpenMM versions (<7.7) don't support boxShape; fall back.
        add_solvent_kwargs.pop("boxShape", None)
        modeller.addSolvent(ff, **add_solvent_kwargs)

    n_atoms_solvated = modeller.topology.getNumAtoms()
    logger.info("Solvated system: %d atoms", n_atoms_solvated)

    # ----- 4. Parameterize: build the OpenMM System -----
    method_map = {
        "nocutoff": "NoCutoff",
        "cutoffnonperiodic": "CutoffNonPeriodic",
        "cutoffperiodic": "CutoffPeriodic",
        "pme": "PME",
        "ewald": "Ewald",
    }
    method_key = method_map.get(str(nonbonded_method).lower())
    if method_key is None:
        raise ValueError(
            f"Unknown nonbonded_method {nonbonded_method!r}. Valid: "
            f"NoCutoff, CutoffNonPeriodic, CutoffPeriodic, PME, Ewald."
        )
    nonbonded_method_obj = omm[method_key]
    is_cutoff_method = method_key in (
        "CutoffNonPeriodic", "CutoffPeriodic", "PME", "Ewald"
    )

    logger.info(
        "Creating OpenMM System (%s, cutoff=%.2f nm)",
        method_key, nonbonded_cutoff_nm,
    )
    create_system_kwargs: dict[str, Any] = {
        "nonbondedMethod": nonbonded_method_obj,
        "constraints": constraints_obj,
        "rigidWater": rigid_water,
        "removeCMMotion": bool(remove_cm_motion),
    }
    # Cutoff/PME/Ewald methods take a cutoff distance + dispersion correction.
    if is_cutoff_method:
        create_system_kwargs["nonbondedCutoff"] = nonbonded_cutoff_nm * unit.nanometer
        create_system_kwargs["useDispersionCorrection"] = bool(dispersion_correction)
    # PME / Ewald take an Ewald error tolerance.
    if method_key in ("PME", "Ewald"):
        create_system_kwargs["ewaldErrorTolerance"] = float(ewald_error_tolerance)
    # Switching function only applies to cutoff-based methods.
    if is_cutoff_method and use_switching_function:
        create_system_kwargs["switchDistance"] = (
            (switch_distance_nm if switch_distance_nm is not None
             else 0.9 * nonbonded_cutoff_nm) * unit.nanometer
        )
    if hydrogen_mass_amu is not None:
        create_system_kwargs["hydrogenMass"] = hydrogen_mass_amu * unit.amu

    # Guard: for periodic cutoff methods, OpenMM requires the nonbonded
    # cutoff to be <= half the smallest box dimension, otherwise it raises a
    # cryptic NonbondedForce error at Context construction. Catch it here
    # with an actionable message. (This normally can't happen with the
    # default 1.0 nm padding + 1.0 nm cutoff, but a user with a small
    # padding or large cutoff could trip it.)
    if method_key in ("CutoffPeriodic", "PME", "Ewald"):
        box_vectors = modeller.topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            try:
                # Smallest box edge length (nm). Take the min diagonal
                # component; works for cubic and triclinic boxes.
                edges_nm = [
                    float(box_vectors[i][i].value_in_unit(unit.nanometer))
                    for i in range(3)
                ]
                min_edge_nm = min(edges_nm)
            except (TypeError, ValueError, AttributeError):
                # Box vectors aren't real numeric quantities (e.g. mocked
                # in tests, or an unexpected type) — skip the guard and let
                # OpenMM handle validation.
                min_edge_nm = None
            if min_edge_nm is not None and nonbonded_cutoff_nm > 0.5 * min_edge_nm:
                raise ValueError(
                    f"Nonbonded cutoff ({nonbonded_cutoff_nm:.2f} nm) exceeds "
                    f"half the smallest periodic box dimension "
                    f"({0.5 * min_edge_nm:.2f} nm; box edge {min_edge_nm:.2f} "
                    f"nm). Increase solvent_padding_nm (currently "
                    f"{solvent_padding_nm:.2f} nm) or decrease "
                    f"nonbonded_cutoff_nm so that the cutoff is at most half "
                    f"the box."
                )

    system = ff.createSystem(modeller.topology, **create_system_kwargs)

    # ----- 5. Capture initial State -----
    # Use a no-op integrator just to obtain a Context for State serialization.
    integrator = omm["openmm"].VerletIntegrator(0.001 * unit.picoseconds)
    context = omm["openmm"].Context(system, integrator)
    context.setPositions(modeller.positions)
    context.setVelocitiesToTemperature(temperature_K * unit.kelvin)
    state = context.getState(
        getPositions=True, getVelocities=True, enforcePeriodicBox=True
    )

    # ----- 6. Serialize to disk -----
    solvated_pdb = out_dir / "solvated.pdb"
    topology_pdb = out_dir / "topology.pdb"  # alias for the analysis phase
    system_xml = out_dir / "system.xml"
    state_xml = out_dir / "state.xml"

    with solvated_pdb.open("w", encoding="utf-8") as fh:
        omm["PDBFile"].writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
    # Symlink-style copy: keep both files (cheap, makes the analysis phase
    # contract simple — it always looks for topology.pdb).
    import shutil as _shutil
    _shutil.copy2(solvated_pdb, topology_pdb)

    with system_xml.open("w", encoding="utf-8") as fh:
        fh.write(omm["openmm"].XmlSerializer.serialize(system))
    with state_xml.open("w", encoding="utf-8") as fh:
        fh.write(omm["openmm"].XmlSerializer.serialize(state))

    logger.info("Wrote: %s, %s, %s, %s", solvated_pdb, topology_pdb, system_xml, state_xml)

    return {
        "solvated_pdb": solvated_pdb,
        "topology_pdb": topology_pdb,
        "system_xml": system_xml,
        "state_xml": state_xml,
    }


def _normalize_ligands(
    ligand: str | Path | list[str | Path] | None,
) -> list[str | Path]:
    """Normalize the ligand argument to a list of paths (possibly empty)."""
    if ligand is None:
        return []
    if isinstance(ligand, (str, Path)):
        return [ligand]
    return list(ligand)


def _build_ligand_forcefield(force_field, small_molecule_ff, ligand_mol):
    """Build an OpenMM ForceField wired for a small-molecule ligand.

    Uses ``openmmforcefields``' ``SystemGenerator`` to combine the protein/
    water force field XMLs with the OpenFF small-molecule force field and the
    ligand molecule. Returns ``(forcefield, system_generator)`` — the
    ``forcefield`` is a standard OpenMM ``ForceField`` (with the ligand
    template generator registered) usable for ``addSolvent`` and
    ``createSystem``.
    """
    try:
        from openmmforcefields.generators import SystemGenerator
    except ImportError as exc:
        from fastmdxplora.setup.ligand import LigandError

        raise LigandError(
            "Ligand parameterization needs openmmforcefields, which is not "
            "installed. Install the ligand extra (pip install "
            "'fastmdxplora[ligand]') or via conda-forge "
            "(conda install -c conda-forge openmmforcefields)."
        ) from exc

    system_generator = SystemGenerator(
        forcefields=list(force_field),
        small_molecule_forcefield=small_molecule_ff,
        molecules=[ligand_mol],
    )
    return system_generator.forcefield, system_generator


def _add_ligand_to_modeller(omm, modeller, ligand_mol, ligand_name="LIG") -> None:
    """Add an OpenFF ligand molecule into the Modeller topology + positions.

    The ligand residue is renamed to ``ligand_name`` so the written topology
    (topology.pdb) and the recorded manifest agree, and so downstream
    ``resname`` selections (e.g. ligand-aware analyses) can find it. OpenFF's
    ``to_openmm()`` otherwise names the residue ``UNK``.
    """
    from openff.units.openmm import to_openmm as _to_openmm

    off_topology = ligand_mol.to_topology()
    omm_topology = off_topology.to_openmm()
    # OpenFF names the ligand residue 'UNK'; set it to the configured name so
    # topology.pdb and setup_parameters.json are consistent and resname-based
    # selection works.
    for residue in omm_topology.residues():
        residue.name = ligand_name
    # Conformer positions -> OpenMM Quantity. OpenFF guarantees at least the
    # conformer loaded from the SDF/MOL2 file.
    positions = _to_openmm(ligand_mol.conformers[0])

    # Record the residue count before adding so we can re-assert the ligand
    # residue name on the MERGED topology — modeller.add() does not reliably
    # preserve the input topology's residue names across all OpenMM versions.
    n_residues_before = modeller.topology.getNumResidues()
    modeller.add(omm_topology, positions)
    for i, residue in enumerate(modeller.topology.residues()):
        if i >= n_residues_before:
            residue.name = ligand_name


def _check_ligand_clashes(
    modeller, n_protein_atoms: int, unit, *, threshold_nm: float, ligand_name: str,
) -> None:
    """Fail at setup if the ligand pose severely overlaps the protein.

    FastMDXplora simulates the protein-ligand complex as provided: the ligand
    coordinates in the SDF/MOL2 must already be a feasible bound pose (e.g.
    from a co-crystal structure or docking). If the supplied pose places
    ligand atoms on top of protein atoms, energy minimization cannot relieve
    the overlap and the simulation diverges to NaN several steps later. We
    detect that here and stop with an actionable message, rather than letting
    it surface as an opaque integration failure downstream.

    Parameters
    ----------
    n_protein_atoms : int
        Number of atoms in the topology before the ligand was added; atoms at
        index >= this are the ligand.
    threshold_nm : float
        Minimum allowed ligand-protein interatomic distance. Pairs closer than
        this count as a clash.
    """
    import math

    try:
        positions = modeller.positions
        n_total = modeller.topology.getNumAtoms()
        # Coordinates in nm as plain floats.
        coords = [
            (
                p.x if hasattr(p, "x") else p[0],
                p.y if hasattr(p, "y") else p[1],
                p.z if hasattr(p, "z") else p[2],
            )
            for p in positions.value_in_unit(unit.nanometer)
        ]
    except (AttributeError, TypeError):
        # Positions aren't real numeric quantities (e.g. mocked in tests) —
        # skip the geometric check.
        return

    protein = coords[:n_protein_atoms]
    ligand = coords[n_protein_atoms:n_total]
    if not ligand or not protein:
        return

    thr_sq = threshold_nm * threshold_nm
    n_clashes = 0
    min_dist_sq = math.inf
    for lx, ly, lz in ligand:
        for px, py, pz in protein:
            d2 = (lx - px) ** 2 + (ly - py) ** 2 + (lz - pz) ** 2
            if d2 < min_dist_sq:
                min_dist_sq = d2
            if d2 < thr_sq:
                n_clashes += 1

    if n_clashes:
        min_dist_nm = math.sqrt(min_dist_sq)
        raise ValueError(
            f"Ligand {ligand_name!r} clashes with the protein: {n_clashes} "
            f"ligand-protein atom pair(s) are closer than "
            f"{threshold_nm:.2f} nm (closest {min_dist_nm:.3f} nm). "
            f"FastMDXplora simulates the complex as provided — the ligand "
            f"coordinates must be a feasible bound pose (from a co-crystal "
            f"structure or docking), not an arbitrary position. Provide a "
            f"properly placed ligand, or if the contact is acceptable, lower "
            f"`ligand_clash_threshold_nm` or set `check_ligand_clashes=False`."
        )
    logger.info(
        "Ligand-protein clash check passed (closest contact %.3f nm).",
        math.sqrt(min_dist_sq) if min_dist_sq != math.inf else 0.0,
    )


def _resolve_constraints(omm: dict, constraints: str):
    """Map a string ``constraints`` argument to an OpenMM enum value."""
    if constraints is None or str(constraints).lower() == "none":
        return None
    mapping = {
        "hbonds": omm["HBonds"],
    }
    # Add more constraints if OpenMM exposes them in this install
    try:
        from openmm.app import AllBonds, HAngles

        mapping["allbonds"] = AllBonds
        mapping["hangles"] = HAngles
    except ImportError:
        pass

    key = str(constraints).lower()
    if key not in mapping:
        raise ValueError(
            f"Unknown constraints option {constraints!r}. Valid: "
            f"None, HBonds, AllBonds, HAngles."
        )
    return mapping[key]
