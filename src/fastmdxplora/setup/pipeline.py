"""Setup pipeline: prepare a simulation-ready system.

This module's public surface is the :func:`run` function called by the
FastMDXplora orchestrator. Starting in v0.2.0 it performs the full
chemistry pipeline:

  1. Resolve input form: file path, PDB ID (4 chars), or sequence.
  2. For a PDB ID, fetch the structure from RCSB.
  3. Run PDBFixer to repair missing residues/atoms and add hydrogens
     at the requested pH (see :func:`fastmdxplora.setup.pdbfix.fix_pdb_with_pdbfixer`).
  4. Solvate, ionize, parameterize with OpenMM and serialize the
     resulting ``System`` + ``State`` XMLs plus a topology PDB
     (see :func:`fastmdxplora.setup.prepare.prepare_system`).

Defaults are sensible general-purpose settings so users
between the two tools see consistent parameterization.

Graceful degradation
--------------------
PDBFixer and OpenMM are conda-forge-only packages in the optional
``[setup]`` extras. When they are not installed, this module still
classifies the input, writes the parameters manifest, and reserves the
canonical artifact paths. It records the ImportError in the manifest
and emits a presenter warning so users see exactly what is missing.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any

from fastmdxplora.dependencies import dependency_error_message, missing_dependencies
from fastmdxplora.utils.logging import get_logger

if TYPE_CHECKING:
    from fastmdxplora.orchestrator import FastMDXplora

logger = get_logger("setup")


# Default parameters. The CHARMM36 + CHARMM36 water choice matches
# sensible protein-only defaults so users
# between the two tools see identical out-of-the-box parameterization.
DEFAULTS: dict[str, Any] = {
    # PDBFixer options
    "ph": 7.4,
    "replace_nonstandard_residues": True,
    "heterogens": "auto",
    "protonation_margin": 1.0,
    "keep_heterogens": False,
    "keep_water": False,
    "fixed_pdb": None,             # skip PDBFixer; use this already-fixed PDB
    # System preparation options
    "forcefield": "auto",      # named selector (resolved to XMLs + water)
    "force_field": None,           # raw XML-list override (power users)
    "water_model": None,           # default derived from the named forcefield
    # Ligand / cofactor (protein-ligand systems; list-shaped, single impl)
    "ligand": None,                # path or list of SDF/MOL2 ligand files
    "ligand_forcefield": None,     # OpenFF small-molecule FF (default per FF)
    "ligand_name": "LIG",          # residue/molecule name
    "ligand_net_charge": None,     # inferred from SDF unless set
    "check_ligand_clashes": True,  # fail setup on a clashing ligand pose
    "ligand_clash_threshold_nm": 0.15,  # min ligand-protein contact (nm)
    "solvent_padding_nm": 1.0,
    "box_shape": "cube",
    "ion_positive": "Na+",
    "ion_negative": "Cl-",
    "ion_concentration_M": 0.15,
    "neutralize": True,
    # createSystem pass-throughs
    "nonbonded_method": "PME",
    "nonbonded_cutoff_nm": 1.0,
    "ewald_error_tolerance": 0.0005,
    "use_switching_function": True,
    "switch_distance_nm": None,    # default: 0.9 * cutoff
    "dispersion_correction": True,
    "remove_cm_motion": True,   # OpenMM's own default
    "constraints": "HBonds",
    "rigid_water": True,
    "hydrogen_mass_amu": None,
    "temperature_K": 300.0,
}


def _classify_input(system: str | None) -> str:
    """Return one of ``{"pdb_file", "pdb_id", "sequence"}``.

    Heuristics:
      - existing path with .pdb / .cif extension -> pdb_file
      - 4-character alphanumeric string -> pdb_id
      - longer alphabetic-only string -> sequence
    """
    if system is None:
        raise ValueError("setup phase requires a system input")

    p = Path(system)
    if p.exists() and p.suffix.lower() in {".pdb", ".cif", ".pdbx"}:
        return "pdb_file"
    if len(system) == 4 and system.isalnum():
        return "pdb_id"
    if system.isalpha():
        return "sequence"
    raise ValueError(
        f"Could not classify system input {system!r}. Expected a PDB file path, "
        "a 4-character PDB ID, or a one-letter amino-acid sequence."
    )


def _fetch_pdb_from_rcsb(pdb_id: str, dest: Path) -> Path:
    """Fetch ``{pdb_id}.pdb`` from RCSB to ``dest``. Raises on HTTP error."""
    import urllib.request

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    logger.info("Fetching PDB from RCSB: %s", url)
    dest.parent.mkdir(parents=True, exist_ok=True)
    urllib.request.urlretrieve(url, dest)  # noqa: S310 -- trusted URL
    return dest


def _resolve_input(
    system: str, input_form: str, setup_dir: Path
) -> Path:
    """Place the source PDB at ``setup_dir/input.pdb``. Returns its path."""
    target = setup_dir / "input.pdb"
    if input_form == "pdb_file":
        shutil.copy2(system, target)
    elif input_form == "pdb_id":
        _fetch_pdb_from_rcsb(system, target)
    elif input_form == "sequence":
        # Sequence -> PDB generation needs a builder (PyRosetta, Modeller,
        # AlphaFold). Out of scope for v0.2; record and bail.
        raise NotImplementedError(
            "Sequence-to-structure not yet supported in the setup phase. "
            "Pass a PDB file or a 4-character PDB ID for now."
        )
    else:
        raise ValueError(f"Unknown input_form {input_form!r}")
    return target



def _keep_heterogens(params: dict, input_pdb) -> bool:
    """Resolve the heterogen policy, reporting what it implies.

    ``drop`` and ``keep`` are unconditional and match the historical options.
    ``auto`` inspects the structure and decides per component, refusing to
    proceed when the structure does not determine the answer.
    """
    from fastmdxplora.setup.heterogens import Action, resolve, summarize

    policy = str(params.get("heterogens", "auto")).strip().lower()
    if params.get("keep_heterogens"):
        policy = "keep"

    if policy == "keep":
        return True

    if policy == "auto":
        # Components worth simulating have already been written as ligand
        # files and routed through the small-molecule path, so the originals
        # are stripped here and re-added with real parameters.
        return False

    if policy != "drop":
        raise ValueError(
            f"heterogens: unknown policy {policy!r}; expected drop, keep, or auto"
        )
    return False




def _validate_ligand_forcefield(params: dict) -> None:
    """Ligand parameterization needs a force field that can do it.

    Checked before any OpenMM or OpenFF import, so the message is the same
    whichever backends happen to be installed.
    """
    from fastmdxplora.setup.forcefields import (
        available_forcefields,
        resolve_forcefield,
    )

    if params["force_field"]:
        raise ValueError(
            "A ligand was supplied with a raw `force_field` XML list. "
            "Protein-ligand parameterization needs the named `forcefield` "
            "selector (use forcefield='amber-openff'), which wires the OpenFF "
            "small-molecule generator."
        )
    choice = resolve_forcefield(params["forcefield"])
    if not choice.supports_ligand:
        valid = ", ".join(
            n for n in available_forcefields() if resolve_forcefield(n).supports_ligand
        )
        raise ValueError(
            f"Force field {choice.name!r} does not support ligands. For "
            f"protein-ligand systems use a ligand-capable force field: {valid}."
        )



def _repaired_complex(input_pdb, setup_dir, ph: float):
    """A structure with missing atoms rebuilt, for the pKa calculation.

    Crystal structures routinely leave surface side chains unmodelled where
    the density is poor. Computing a pKa in a structure with those gaps means
    computing it in the wrong electrostatic environment, so the ligand and the
    protein are repaired first. Falls back to the deposition if repair fails,
    since a less accurate pKa is better than none.
    """
    from fastmdxplora.setup.pdbfix import fix_pdb_with_pdbfixer

    # setup_dir is already the setup directory, as it is for the ligand
    # files: appending "setup" again buried this in setup/setup/.
    repaired = Path(setup_dir) / "complex_for_pka.pdb"
    repaired.parent.mkdir(parents=True, exist_ok=True)
    try:
        fix_pdb_with_pdbfixer(
            str(input_pdb), str(repaired), ph=ph,
            keep_heterogens=True,   # the ligand must be present to be assessed
            keep_water=False,
        )
        return repaired
    except Exception as exc:  # noqa: BLE001 - fall back rather than fail
        logger.warning(
            "Could not repair the structure for the pKa calculation (%s); "
            "using the deposition as given, which may have unmodelled side "
            "chains near the site.", exc,
        )
        return input_pdb



def _retain_in_structure(input_pdb, setup_dir, keep_decisions):
    """Write a structure holding the polymer plus the components to retain.

    PDBFixer keeps heterogens all or nothing, so selecting a few means
    filtering the input first. Everything the classifier discarded is dropped
    here; what remains is prepared by the protein force field, in place, under
    its own residue name, at its own coordinates.
    """
    keep = {d.resname.upper() for d in keep_decisions}
    target = Path(setup_dir) / "retained.pdb"
    target.parent.mkdir(parents=True, exist_ok=True)

    from fastmdxplora.setup.heterogens import ION_NAMES

    lines: list[str] = []
    dropped: set[str] = set()
    ions: set[str] = set()
    for raw in Path(input_pdb).read_text(
        encoding="utf-8", errors="ignore"
    ).splitlines():
        record = raw[:6].strip()
        if record == "HETATM" and raw[17:20].strip().upper() not in keep:
            dropped.add(raw[6:11].strip())
            continue
        if record in ("ATOM", "HETATM") and raw[17:20].strip().upper() in ION_NAMES:
            ions.add(raw[6:11].strip())
        lines.append(raw)

    # CONECT records name atoms by serial number, and OpenMM builds bonds from
    # them. Left pointing at atoms that are no longer here, they bond whatever
    # now holds those serials: a filtered structure produced a bond between
    # two zincs 9.6 A apart, which the force field rightly rejected. Any
    # record naming a removed atom is dropped with it.
    #
    # A record naming a retained ion is dropped too, and for a separate reason.
    # The legacy PDB format writes coordination and covalency with the same
    # record type -- only mmCIF separates them -- so a CONECT tying a metal to
    # the residue holding it asserts a bond that is not there. The classifier
    # above already applies that rule to LINK records; it holds here as well.
    # The force field could not honour the bond in any case: these ions are
    # modelled as Lennard-Jones plus a charge, with no bonded term to apply,
    # and OpenMM's ion templates permit no external bonds, so the bond makes
    # the ion unmatchable. 1ZNI reaches this, with two zincs and three
    # chlorides tied to each other.
    kept: list[str] = []
    coordination = 0
    for raw in lines:
        if raw[:6].strip() != "CONECT":
            kept.append(raw)
            continue
        serials = [raw[i:i + 5].strip() for i in range(6, min(len(raw), 31), 5)]
        present = [serial for serial in serials if serial]
        if any(serial in dropped for serial in present):
            continue
        if any(serial in ions for serial in present):
            coordination += 1
            continue
        kept.append(raw)

    if coordination:
        logger.info(
            "Dropped %d connectivity record(s) naming a retained ion: they "
            "describe coordination rather than a covalent bond, which the "
            "non-bonded ion parameters already represent.",
            coordination,
        )

    target.write_text("\n".join(kept) + "\n", encoding="utf-8")
    return target


def _auto_ligands(params: dict, input_pdb, setup_dir, entry_id: str | None) -> list[str]:
    """Turn the classifier's SIMULATE decisions into files the ligand path takes.

    The chemistry a structure omits is fetched from RCSB, completed with
    hydrogens, and checked for a determined protonation in the complex. What
    comes back is an ordinary ligand file, so everything downstream is the
    existing protein-ligand route rather than a parallel one.

    Returns the ligand files written, empty when the structure holds nothing
    worth simulating.
    """
    from fastmdxplora.setup.ccd import fetch_chemistry
    from fastmdxplora.setup.forcefields import (
        available_forcefields,
        resolve_forcefield,
    )
    from fastmdxplora.setup.heterogens import Action, resolve, summarize
    from fastmdxplora.setup.protonation import (
        POISED_MARGIN, apply_settled_state, settle,
    )

    decisions = resolve(input_pdb, keep_water=bool(params.get("keep_water")))
    logger.info("Heterogen decisions:\n%s", summarize(decisions))
    wanted = [d for d in decisions if d.action is Action.SIMULATE]
    if not wanted:
        return []

    names = ", ".join(
        f"{d.resname} x{d.count}" if d.count > 1 else d.resname for d in wanted
    )

    if entry_id is None:
        raise ValueError(
            f"heterogens: auto identified components to simulate ({names}), but "
            "their chemistry can only be retrieved for a structure given by PDB "
            "identifier: a local file carries no entry to look them up in. "
            "Supply the ligands as SDF or MOL2 files instead."
        )

    # A ligand needs a force field that can parameterize small molecules. The
    # protein force field is a scientific choice, so it is not changed here on
    # the user's behalf.
    choice = resolve_forcefield(params.get("forcefield"))
    if not choice.supports_ligand:
        capable = ", ".join(
            n for n in available_forcefields() if resolve_forcefield(n).supports_ligand
        )
        raise ValueError(
            f"heterogens: auto identified components to simulate ({names}), but "
            f"force field {choice.name!r} cannot parameterize small molecules. "
            f"Choose a ligand-capable force field ({capable}), or set the "
            "heterogens policy to 'drop' to exclude them."
        )

    # Monatomic ions are kept in the structure rather than extracted. The
    # protein force field already carries validated ion parameters, and the
    # small-molecule route is both unnecessary and destructive: a zinc pulled
    # out to an SDF and re-added came back renamed, in a different chain, and
    # at coordinates that were not its coordination site.
    from fastmdxplora.setup.heterogens import ION_NAMES

    ions = [d for d in wanted if d.resname in ION_NAMES]
    molecules = [d for d in wanted if d.resname not in ION_NAMES]
    if ions:
        params["_retained_pdb"] = str(
            _retain_in_structure(input_pdb, setup_dir, ions)
        )
        logger.info(
            "Keeping %s in the structure; the protein force field provides "
            "ion parameters.",
            ", ".join(f"{d.resname} x{d.count}" if d.count > 1 else d.resname
                      for d in ions),
        )
    if not molecules:
        return []

    copies = [(d, het) for d in molecules for het in d.instances]
    # setup_dir is already the setup directory: input.pdb and prepared.pdb
    # sit directly in it. Appending "setup" again buried the ligands in
    # setup/setup/ligands.
    ligand_dir = Path(setup_dir) / "ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)

    files: list[str] = []
    names: list[str] = []
    charges: list[int] = []
    for decision, instance in copies:
        chemistry = fetch_chemistry(
            entry_id,
            decision.resname,
            chain=instance.chain,
            resseq=instance.resseq,
            expected_heavy_atoms=sum(
                1 for atom in instance.atoms if atom.element.upper() != "H"
            ),
        )
        # A pKa is a property of the environment, so the environment should be
        # the one that will be simulated: repaired, with the missing side
        # chains rebuilt. The raw deposition has gaps that would bias the
        # answer. Only built when the ligand is actually titratable, since
        # repairing costs a second and otherwise changes nothing.
        pka_structure = input_pdb
        if chemistry.titratable_groups:
            pka_structure = _repaired_complex(input_pdb, setup_dir, float(params["ph"]))

        state = settle(
            pka_structure,
            decision.resname,
            chain=instance.chain,
            resseq=instance.resseq,
            ph=float(params["ph"]),
            expected_ionizable=bool(chemistry.titratable_groups),
            margin=float(params.get("protonation_margin") or POISED_MARGIN),
            known_groups=tuple(chemistry.titratable_groups),
        )
        logger.info("%s %s%s: %s", decision.resname, instance.chain,
                    instance.resseq, state.reason)

        # The settled state has to be the one that reaches the force field.
        sdf_text, net_charge = apply_settled_state(
            chemistry.path.read_text(encoding="utf-8"), chemistry, state
        )

        # Copies of one component need distinct residue names, or nothing
        # downstream can tell them apart.
        suffix = "" if len(copies) == 1 else f"_{instance.chain}{instance.resseq}"
        target = ligand_dir / f"{decision.resname}{suffix}.sdf"
        target.write_text(sdf_text, encoding="utf-8")
        files.append(str(target))
        names.append(decision.resname if len(copies) == 1
                     else f"{decision.resname[:1]}{len(names):02d}")
        charges.append(net_charge)

    # Remembered so preparation can report these as re-added with parameters
    # rather than warning that they were removed.
    params["_reinstated_heterogens"] = tuple(
        sorted({decision.resname for decision, _ in copies})
    )
    # Every component the classifier judged, so preparation does not warn
    # about ones whose fate has already been reported with a reason.
    params["_explained_heterogens"] = tuple(sorted({d.resname for d in decisions}))
    params["ligand_name"] = names[0] if len(names) == 1 else names
    if params.get("ligand_net_charge") is None:
        params["ligand_net_charge"] = charges[0] if len(charges) == 1 else charges
    return files


def run(
    *,
    orchestrator: "FastMDXplora",
    output_dir: Path,
    **options: Any,
) -> list[str]:
    """Run the setup phase.

    Parameters
    ----------
    orchestrator : FastMDXplora
        The parent orchestrator instance.
    setup_dir : pathlib.Path
        Where to write setup artifacts.
    **options
        Per-call overrides of the module-level :data:`DEFAULTS`.

    Returns
    -------
    list of str
        Paths (relative to ``setup_dir``) of artifacts produced.
    """
    # The orchestrator hands each phase its own directory, under the name every
    # phase uses. Here it is the setup directory, and calling it output_dir
    # invited "setup" being appended to it -- which happened to the ligand
    # files, to the structure built for the pKa calculation, and would have
    # happened again.
    setup_dir = output_dir

    params: dict[str, Any] = {**DEFAULTS, **options}

    # Force-field selection is either the named selector OR a raw XML list,
    # not both. Check the user-supplied options (not merged defaults), since
    # `forcefield` always has a default.
    if options.get("forcefield") is not None and options.get("force_field"):
        raise ValueError(
            "Specify either `forcefield` (a named force field) or "
            "`force_field` (a raw list of OpenMM XML files), not both. "
            "The named selector is recommended; the raw list is an escape "
            "hatch for combinations the named registry does not cover."
        )
    # Surface an unknown named force field early with a clear message
    # (only when a raw list is not overriding it).
    if not params["force_field"]:
        from fastmdxplora.setup.forcefields import resolve_forcefield

        resolve_forcefield(params["forcefield"])  # raises ValueError if unknown

    # Ligand / force-field coherence (pure logic — validate here, before any
    # OpenMM/OpenFF dependency, so the user gets a clear error regardless of
    # which backends are installed).
    if params.get("ligand"):
        _validate_ligand_forcefield(params)

    input_form = _classify_input(orchestrator.system)
    logger.debug("setup: input form detected as %s", input_form)

    presenter = getattr(orchestrator, "_presenter", None)
    artifacts: list[str] = []
    notes: list[str] = []

    # ---- Stage 1: resolve input ----------------------------------------
    try:
        input_pdb = _resolve_input(orchestrator.system, input_form, setup_dir)

        artifacts.append("input.pdb")
        if presenter:
            label = (
                orchestrator.system if input_form == "pdb_id"
                else Path(orchestrator.system).name
            )
            presenter.step(f"Loaded input: {label}")
    except NotImplementedError as exc:
        # Sequence input — manifest-only fallback
        (setup_dir / "input.sequence").write_text(f"{orchestrator.system}\n", encoding="utf-8")
        artifacts.append("input.sequence")
        notes.append(str(exc))
        if presenter:
            presenter.step(str(exc), status="warning")
        _write_manifest(setup_dir, orchestrator, input_form, params, artifacts, notes)
        artifacts.append("setup_parameters.json")
        return artifacts
    except Exception as exc:  # noqa: BLE001 -- network, IO, refusals
        # The manifest is written first, so whatever was learned before the
        # failure survives for diagnosis. Then the failure is raised.
        #
        # It used to be swallowed, on the reasoning that the manifest was the
        # source of truth. In practice that meant a mistyped PDB identifier
        # produced an empty output directory and exit code 0, and any refusal
        # raised here was reported as a successful run. A phase that could not
        # do its job must say so through the channel callers actually read.
        #
        # Note that this is a different situation from an absent optional
        # dependency, which is handled below and does degrade gracefully by
        # design: not having OpenMM installed is a choice, whereas failing to
        # fetch a structure is a failure.
        notes.append(f"Failed to resolve input ({input_form}): {exc}")
        _write_manifest(setup_dir, orchestrator, input_form, params, artifacts, notes)
        raise

    # With the auto policy, the structure decides what to simulate and the
    # chemistry it omits is retrieved before preparation runs, so the rest of
    # the phase sees an ordinary protein-ligand job.
    #
    # Deliberately outside the input-resolution try above: that handler
    # downgrades failures to a warning and carries on, which would turn a
    # refusal into a run that reports success while having simulated nothing
    # of the kind.
    if str(params.get("heterogens", "auto")).strip().lower() == "auto" \
            and not params.get("keep_heterogens") and not params.get("ligand"):
        discovered = _auto_ligands(
            params,
            input_pdb,
            setup_dir,
            orchestrator.system if input_form == "pdb_id" else None,
        )
        if discovered:
            params["ligand"] = discovered
            _validate_ligand_forcefield(params)
            logger.info(
                "Prepared %d ligand file(s) from the structure: %s",
                len(discovered), ", ".join(Path(d).name for d in discovered),
            )

    # ---- Stage 2: PDBFixer (or skip via fixed_pdb) ---------------------
    prepared_pdb = setup_dir / "prepared.pdb"
    fixed_pdb = params.get("fixed_pdb")
    if fixed_pdb:
        # User supplied an already-fixed PDB — skip PDBFixer entirely.
        fixed_src = Path(fixed_pdb)
        if not fixed_src.exists():
            notes.append(f"fixed_pdb not found: {fixed_src}")
            if presenter:
                presenter.step(f"fixed_pdb not found: {fixed_src}", status="warning")
            _write_manifest(setup_dir, orchestrator, input_form, params, artifacts, notes)
            artifacts.append("setup_parameters.json")
            return artifacts
        shutil.copy2(fixed_src, prepared_pdb)
        artifacts.append("prepared.pdb")
        if presenter:
            presenter.step(f"Using supplied fixed PDB: {fixed_src.name} (PDBFixer skipped)")
    else:
        try:
            from fastmdxplora.setup.pdbfix import fix_pdb_with_pdbfixer

            # When ions are being kept, PDBFixer runs on a structure already
            # filtered to the polymer plus those ions, so "keep heterogens"
            # retains exactly them and nothing else.
            retained = params.get("_retained_pdb")
            fix_pdb_with_pdbfixer(
                retained or str(input_pdb),
                str(prepared_pdb),
                ph=float(params["ph"]),
                keep_heterogens=(True if params.get("_retained_pdb")
                                 else _keep_heterogens(params, input_pdb)),
                keep_water=bool(params["keep_water"]),
                reinstated=tuple(params.get("_reinstated_heterogens", ())),
                explained=tuple(params.get("_explained_heterogens", ())),
                replace_nonstandard=bool(params["replace_nonstandard_residues"]),
            )
            artifacts.append("prepared.pdb")
            if presenter:
                presenter.step(f"Fixed PDB with PDBFixer (pH={params['ph']})")
        except ImportError as exc:
            missing = missing_dependencies()
            notes.append(
                dependency_error_message(missing)
                if missing
                else f"PDBFixer unavailable: {exc}"
            )
            if presenter:
                presenter.step(
                    notes[-1],
                    status="warning",
                )
            _write_manifest(setup_dir, orchestrator, input_form, params, artifacts, notes)
            artifacts.append("setup_parameters.json")
            return artifacts

    # ---- Stage 3: Solvate, ionize, parameterize, serialize -------------
    try:
        from fastmdxplora.setup.prepare import prepare_system

        produced = prepare_system(
            prepared_pdb,
            setup_dir,
            forcefield=params["forcefield"],
            force_field=params["force_field"],
            water_model=params["water_model"],
            ligand=params["ligand"],
            ligand_forcefield=params["ligand_forcefield"],
            # A list of names, one per ligand, must not be stringified: doing
            # so produced a single "name" of "['Z00', 'Z01', 'Z02']" and
            # defeated multi-ligand support one line after it was added.
            ligand_name=(
                params["ligand_name"]
                if isinstance(params["ligand_name"], (list, tuple))
                else str(params["ligand_name"])
            ),
            ligand_net_charge=params["ligand_net_charge"],
            check_ligand_clashes=bool(params["check_ligand_clashes"]),
            ligand_clash_threshold_nm=float(params["ligand_clash_threshold_nm"]),
            solvent_padding_nm=float(params["solvent_padding_nm"]),
            box_shape=str(params["box_shape"]),
            ion_positive=str(params["ion_positive"]),
            ion_negative=str(params["ion_negative"]),
            ion_concentration_M=float(params["ion_concentration_M"]),
            neutralize=bool(params["neutralize"]),
            nonbonded_method=str(params["nonbonded_method"]),
            nonbonded_cutoff_nm=float(params["nonbonded_cutoff_nm"]),
            ewald_error_tolerance=float(params["ewald_error_tolerance"]),
            use_switching_function=bool(params["use_switching_function"]),
            switch_distance_nm=params["switch_distance_nm"],
            dispersion_correction=bool(params["dispersion_correction"]),
            remove_cm_motion=bool(params["remove_cm_motion"]),
            constraints=str(params["constraints"]),
            rigid_water=bool(params["rigid_water"]),
            hydrogen_mass_amu=params["hydrogen_mass_amu"],
            temperature_K=float(params["temperature_K"]),
        )

        for _key, path in produced.items():
            artifacts.append(path.relative_to(setup_dir).as_posix())
        if presenter:
            if params["force_field"]:
                ff_label = ", ".join(params["force_field"])
            else:
                ff_label = str(params["forcefield"])
            presenter.step(f"Solvated and parameterized ({ff_label})")
            presenter.step("Wrote system.xml, state.xml, topology.pdb")
    except ImportError as exc:
        missing = missing_dependencies()
        notes.append(
            dependency_error_message(missing)
            if missing
            else f"OpenMM unavailable for parameterization: {exc}"
        )
        if presenter:
            presenter.step(
                notes[-1],
                status="warning",
            )
        _write_manifest(setup_dir, orchestrator, input_form, params, artifacts, notes)
        artifacts.append("setup_parameters.json")
        return artifacts

    # ---- Stage 4: Manifest --------------------------------------------
    _write_manifest(setup_dir, orchestrator, input_form, params, artifacts, notes)
    artifacts.append("setup_parameters.json")

    if presenter:
        presenter.step("Wrote setup_parameters.json")

    logger.debug("setup: wrote %d artifact(s) to %s", len(artifacts), setup_dir)
    return artifacts


def _write_manifest(
    setup_dir: Path,
    orchestrator: "FastMDXplora",
    input_form: str,
    params: dict[str, Any],
    artifacts: list[str],
    notes: list[str],
) -> None:
    """Write ``setup_parameters.json`` with the full provenance record."""
    # Record what the force-field selection actually resolved to, so the
    # manifest is reproducible regardless of whether the user picked a named
    # force field or passed a raw XML list.
    resolved_ff: dict[str, Any]
    if params.get("force_field"):
        resolved_ff = {
            "source": "explicit_xml_list",
            "xmls": list(params["force_field"]),
            "water_model": params.get("water_model"),
        }
    else:
        from fastmdxplora.setup.forcefields import resolve_forcefield

        try:
            choice = resolve_forcefield(params.get("forcefield"))
            resolved_ff = {
                "source": "named",
                "name": choice.name,
                "xmls": list(choice.xmls),
                "water_model": params.get("water_model") or choice.water_model,
                "supports_ligand": choice.supports_ligand,
                "small_molecule_forcefield": choice.small_molecule_forcefield,
            }
        except ValueError:
            resolved_ff = {"source": "named", "name": params.get("forcefield")}

    # Record ligand parameterization when a ligand was supplied.
    if params.get("ligand"):
        ligand_files = (
            params["ligand"] if isinstance(params["ligand"], list)
            else [params["ligand"]]
        )
        resolved_ff["ligand"] = {
            "files": [str(f) for f in ligand_files],
            "name": params.get("ligand_name", "LIG"),
            "small_molecule_forcefield": (
                params.get("ligand_forcefield")
                or resolved_ff.get("small_molecule_forcefield")
            ),
            "net_charge": params.get("ligand_net_charge"),
        }
    canonical = {
        "input_pdb": "input.pdb",
        "prepared_pdb": "prepared.pdb",
        "solvated_pdb": "solvated.pdb",
        "topology_pdb": "topology.pdb",
        "system_xml": "system.xml",
        "state_xml": "state.xml",
    }
    manifest = {
        "phase": "setup",
        "input": {
            "system": orchestrator.system,
            "form": input_form,
        },
        "parameters": params,
        "resolved_forcefield": resolved_ff,
        "artifacts_planned": canonical,
        "artifacts_written": list(artifacts),
        "notes": notes,
    }
    with (setup_dir / "setup_parameters.json").open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)
