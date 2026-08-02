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
    "ph": 7.0,
    "heterogens": "drop",
    "keep_heterogens": False,
    "keep_water": False,
    "fixed_pdb": None,             # skip PDBFixer; use this already-fixed PDB
    # System preparation options
    "forcefield": "charmm36",      # named selector (resolved to XMLs + water)
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
    "remove_cm_motion": False,
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
    system: str, input_form: str, output_dir: Path
) -> Path:
    """Place the source PDB at ``output_dir/input.pdb``. Returns its path."""
    target = output_dir / "input.pdb"
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

    policy = str(params.get("heterogens", "drop")).strip().lower()
    if params.get("keep_heterogens"):
        policy = "keep"

    if policy == "keep":
        return True

    if policy == "auto":
        # Raises AmbiguousStructureError when the structure is undetermined,
        # which is the intended outcome rather than a failure to handle.
        decisions = resolve(input_pdb, keep_water=bool(params.get("keep_water")))
        logger.info("Heterogen decisions:\n%s", summarize(decisions))
        wanted = [d for d in decisions if d.action is Action.SIMULATE]
        if wanted:
            names = ", ".join(f"{d.resname} x{d.count}" if d.count > 1 else d.resname
                              for d in wanted)
            raise ValueError(
                f"heterogens: auto identified components to simulate ({names}), "
                "but automatic ligand parameterization is not available yet. "
                "Supply them as SDF or MOL2 files with the ligand option "
                "(--setup-ligand under explore, --ligand under setup), or set "
                "the heterogens policy to 'drop' to exclude them."
            )
        # Everything was solvent or additive: dropping is the right call.
        return False

    if policy != "drop":
        raise ValueError(
            f"heterogens: unknown policy {policy!r}; expected drop, keep, or auto"
        )
    return False


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
    output_dir : pathlib.Path
        Where to write setup artifacts.
    **options
        Per-call overrides of the module-level :data:`DEFAULTS`.

    Returns
    -------
    list of str
        Paths (relative to ``output_dir``) of artifacts produced.
    """
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
        from fastmdxplora.setup.forcefields import (
            available_forcefields,
            resolve_forcefield,
        )

        if params["force_field"]:
            raise ValueError(
                "A ligand was supplied with a raw `force_field` XML list. "
                "Protein-ligand parameterization needs the named "
                "`forcefield` selector (use forcefield='amber-openff'), "
                "which wires the OpenFF small-molecule generator."
            )
        choice = resolve_forcefield(params["forcefield"])
        if not choice.supports_ligand:
            valid = ", ".join(
                n for n in available_forcefields()
                if resolve_forcefield(n).supports_ligand
            )
            raise ValueError(
                f"Force field {choice.name!r} does not support ligands. "
                f"For protein-ligand systems use a ligand-capable force "
                f"field: {valid}."
            )

    input_form = _classify_input(orchestrator.system)
    logger.debug("setup: input form detected as %s", input_form)

    presenter = getattr(orchestrator, "_presenter", None)
    artifacts: list[str] = []
    notes: list[str] = []

    # ---- Stage 1: resolve input ----------------------------------------
    try:
        input_pdb = _resolve_input(orchestrator.system, input_form, output_dir)
        artifacts.append("input.pdb")
        if presenter:
            label = (
                orchestrator.system if input_form == "pdb_id"
                else Path(orchestrator.system).name
            )
            presenter.step(f"Loaded input: {label}")
    except NotImplementedError as exc:
        # Sequence input — manifest-only fallback
        (output_dir / "input.sequence").write_text(f"{orchestrator.system}\n", encoding="utf-8")
        artifacts.append("input.sequence")
        notes.append(str(exc))
        if presenter:
            presenter.step(str(exc), status="warning")
        _write_manifest(output_dir, orchestrator, input_form, params, artifacts, notes)
        artifacts.append("setup_parameters.json")
        return artifacts
    except Exception as exc:  # noqa: BLE001 -- network, IO, etc.
        # Anything else (RCSB fetch failure, IO error) becomes a graceful
        # degradation: write the manifest with a note explaining what
        # went wrong and stop. The orchestrator's per-phase error handler
        # still sees a clean "ok" return; the manifest is the source of
        # truth for what actually happened.
        notes.append(f"Failed to resolve input ({input_form}): {exc}")
        if presenter:
            presenter.step(
                f"Could not resolve input ({input_form}): {exc}",
                status="warning",
            )
        _write_manifest(output_dir, orchestrator, input_form, params, artifacts, notes)
        artifacts.append("setup_parameters.json")
        return artifacts

    # ---- Stage 2: PDBFixer (or skip via fixed_pdb) ---------------------
    prepared_pdb = output_dir / "prepared.pdb"
    fixed_pdb = params.get("fixed_pdb")
    if fixed_pdb:
        # User supplied an already-fixed PDB — skip PDBFixer entirely.
        fixed_src = Path(fixed_pdb)
        if not fixed_src.exists():
            notes.append(f"fixed_pdb not found: {fixed_src}")
            if presenter:
                presenter.step(f"fixed_pdb not found: {fixed_src}", status="warning")
            _write_manifest(output_dir, orchestrator, input_form, params, artifacts, notes)
            artifacts.append("setup_parameters.json")
            return artifacts
        shutil.copy2(fixed_src, prepared_pdb)
        artifacts.append("prepared.pdb")
        if presenter:
            presenter.step(f"Using supplied fixed PDB: {fixed_src.name} (PDBFixer skipped)")
    else:
        try:
            from fastmdxplora.setup.pdbfix import fix_pdb_with_pdbfixer

            fix_pdb_with_pdbfixer(
                str(input_pdb),
                str(prepared_pdb),
                ph=float(params["ph"]),
                keep_heterogens=_keep_heterogens(params, input_pdb),
                keep_water=bool(params["keep_water"]),
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
            _write_manifest(output_dir, orchestrator, input_form, params, artifacts, notes)
            artifacts.append("setup_parameters.json")
            return artifacts

    # ---- Stage 3: Solvate, ionize, parameterize, serialize -------------
    try:
        from fastmdxplora.setup.prepare import prepare_system

        produced = prepare_system(
            prepared_pdb,
            output_dir,
            forcefield=params["forcefield"],
            force_field=params["force_field"],
            water_model=params["water_model"],
            ligand=params["ligand"],
            ligand_forcefield=params["ligand_forcefield"],
            ligand_name=str(params["ligand_name"]),
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
            artifacts.append(path.relative_to(output_dir).as_posix())
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
        _write_manifest(output_dir, orchestrator, input_form, params, artifacts, notes)
        artifacts.append("setup_parameters.json")
        return artifacts

    # ---- Stage 4: Manifest --------------------------------------------
    _write_manifest(output_dir, orchestrator, input_form, params, artifacts, notes)
    artifacts.append("setup_parameters.json")

    if presenter:
        presenter.step("Wrote setup_parameters.json")

    logger.debug("setup: wrote %d artifact(s) to %s", len(artifacts), output_dir)
    return artifacts


def _write_manifest(
    output_dir: Path,
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
    with (output_dir / "setup_parameters.json").open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)
