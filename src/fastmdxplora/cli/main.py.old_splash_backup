"""``fastmdx`` command-line entry point.

Subcommands
-----------

  explore / xplore   Run the full pipeline (setup → simulation → analysis → report)
  setup              Run only the setup phase
  simulate           Run only the simulation phase
  analyze            Run only the analysis phase
  report             Run only the report phase
  info               Print environment and component info

Each per-phase subcommand exposes the phase's own options (e.g.
``--ph`` for setup, ``--duration-ns`` for simulate). The ``explore``
subcommand additionally exposes those same options under per-phase
prefixes (``--setup-ph``, ``--simulate-duration-ns``) so a user can drive
the full pipeline from a single invocation.

Global flags:

  --version, -V
  --cite                 Print the citation and exit
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Sequence

from fastmdxplora import (
    __author__,
    __citation__,
    __doi__,
    __expansion__,
    __version__,
)
from fastmdxplora.orchestrator import FastMDXplora


# ---------------------------------------------------------------------------
# Phase-option definitions
# ---------------------------------------------------------------------------
# Each entry: (cli_flag_suffix, kwarg_name, argparse-kwargs)
#
# The cli_flag_suffix is appended to "--" for the per-phase subcommands
# (e.g. setup's ``--ph``) and to "--<prefix>-" for the `explore`
# subcommand (e.g. ``--setup-ph``).
#
# The kwarg_name is what gets passed to the phase's run() function.

_SETUP_OPTIONS: list[tuple[str, str, dict[str, Any]]] = [
    ("ph", "ph", {"type": float, "help": "pH for hydrogen placement (default 7.0)."}),
    ("keep-heterogens", "keep_heterogens", {"action": "store_true", "default": None,
        "help": "Retain non-standard residues (ligands, cofactors, ions)."}),
    ("keep-water", "keep_water", {"action": "store_true", "default": None,
        "help": "Retain crystallographic waters."}),
    ("fixed-pdb", "fixed_pdb", {"type": str, "metavar": "PATH",
        "help": "Use an already-fixed PDB and skip PDBFixer."}),
    ("forcefield", "forcefield", {"choices": ["charmm36", "amber14", "amber-fb15", "amber-openff"],
        "help": "Named force field (default 'charmm36'). Resolves to the "
                "right XMLs and water model. Use 'amber-openff' for ligands."}),
    ("force-field", "force_field", {"nargs": "+", "metavar": "XML",
        "help": "Raw OpenMM XML(s), overriding --forcefield (power users)."}),
    ("ligand", "ligand", {"nargs": "+", "metavar": "FILE",
        "help": "Ligand SDF/MOL2 file(s) (needs --setup-forcefield amber-openff)."}),
    ("ligand-forcefield", "ligand_forcefield", {"type": str, "metavar": "NAME",
        "help": "OpenFF small-molecule force field (e.g. openff-2.2.1)."}),
    ("ligand-name", "ligand_name", {"type": str, "metavar": "NAME",
        "help": "Residue/molecule name for the ligand (default 'LIG')."}),
    ("ligand-net-charge", "ligand_net_charge", {"type": int, "metavar": "INT",
        "help": "Ligand formal net charge (default: inferred from SDF)."}),
    ("no-ligand-clash-check", "check_ligand_clashes", {"action": "store_false",
        "default": None,
        "help": "Skip the ligand-protein clash check at setup."}),
    ("ligand-clash-threshold-nm", "ligand_clash_threshold_nm", {"type": float,
        "metavar": "NM",
        "help": "Min ligand-protein contact distance in nm (default 0.15)."}),
    ("water-model", "water_model", {"type": str, "metavar": "NAME",
        "help": "Water model for Modeller (e.g. 'tip3p', 'tip4pew')."}),
    ("solvent-padding-nm", "solvent_padding_nm", {"type": float,
        "help": "Min distance between solute and box wall in nm (default 1.0)."}),
    ("box-shape", "box_shape", {"choices": ["cube", "dodecahedron", "octahedron"],
        "help": "Periodic box geometry (default 'cube')."}),
    ("nonbonded-method", "nonbonded_method", {"choices": [
        "NoCutoff", "CutoffNonPeriodic", "CutoffPeriodic", "PME", "Ewald"],
        "help": "Nonbonded method (default 'PME')."}),
    ("ion-positive", "ion_positive", {"type": str, "metavar": "ION",
        "help": "Counter-ion cation (default 'Na+')."}),
    ("ion-negative", "ion_negative", {"type": str, "metavar": "ION",
        "help": "Counter-ion anion (default 'Cl-')."}),
    ("ion-concentration-M", "ion_concentration_M", {"type": float,
        "help": "Target ionic strength in M (default 0.15)."}),
    ("temperature-K", "temperature_K", {"type": float,
        "help": "Initial velocity temperature in K (default 300)."}),
]

_SIMULATION_OPTIONS: list[tuple[str, str, dict[str, Any]]] = [
    ("preset", "preset", {"choices": ["gentle"],
        "help": "Simulation preset. 'gentle' uses conservative smoke-test settings."}),
    ("duration-ns", "duration_ns", {"type": float,
        "help": "Production length in ns (standard MD convention; equilibration is independent)."}),
    ("nvt-duration-ns", "nvt_duration_ns", {"type": float,
        "help": "NVT equilibration in ns (default: fixed 500 ps regardless of production length)."}),
    ("npt-duration-ns", "npt_duration_ns", {"type": float,
        "help": "NPT equilibration in ns (default: fixed 1 ns regardless of production length)."}),
    ("nvt-steps", "nvt_steps", {"type": int,
        "help": "NVT step count (overrides --nvt-duration-ns)."}),
    ("npt-steps", "npt_steps", {"type": int,
        "help": "NPT step count (overrides --npt-duration-ns)."}),
    ("production-steps", "production_steps", {"type": int,
        "help": "Production step count (overrides --duration-ns)."}),
    ("timestep-fs", "timestep_fs", {"type": float,
        "help": "Integrator timestep in fs (default 2.0)."}),
    ("integrator", "integrator", {"choices": [
        "langevin_middle", "langevin", "brownian", "verlet",
        "variable_langevin", "variable_verlet"],
        "help": "Integrator (default 'langevin_middle')."}),
    ("temperature-K", "temperature_K", {"type": float,
        "help": "Production temperature in K (default 300)."}),
    ("pressure-bar", "pressure_bar", {"type": float,
        "help": "Barostat pressure in bar (OpenMM-native; default 1.0)."}),
    ("pressure-atm", "pressure_atm", {"type": float,
        "help": "Barostat pressure in atm (converted to bar)."}),
    ("friction-per-ps", "friction_per_ps", {"type": float,
        "help": "Langevin friction in 1/ps (default 1.0)."}),
    ("platform", "platform", {"choices": ["auto", "CUDA", "OpenCL", "CPU", "HIP"],
        "help": "OpenMM compute platform (default 'auto': CUDA → OpenCL → CPU)."}),
    ("precision", "precision", {"choices": ["single", "mixed", "double"],
        "help": "GPU precision (default 'mixed')."}),
    ("device-index", "device_index", {"type": str, "metavar": "IDX",
        "help": "GPU device index for multi-GPU machines (e.g. '0' or '0,1')."}),
    ("checkpoint-interval-steps", "checkpoint_interval_steps", {"type": int,
        "help": "Checkpoint (.chk) interval in steps; 0 disables (default 10000)."}),
    ("trajectory-interval-steps", "trajectory_interval_steps", {"type": int,
        "help": "Trajectory (.dcd) frame interval in steps (default: adaptive, ~2000 frames)."}),
    ("random-seed", "random_seed", {"type": int,
        "help": "Integrator random seed (default: not set)."}),
    ("plumed-script", "plumed_script", {"type": str, "metavar": "PATH",
        "help": "Enable PLUMED enhanced sampling with this script (path to a "
                ".dat file or inline text). Requires openmm-plumed."}),
    ("no-minimize", "minimize", {"action": "store_false", "default": None,
        "help": "Skip the energy minimization stage."}),
]

_ANALYSIS_OPTIONS: list[tuple[str, str, dict[str, Any]]] = [
    ("trajectory", "trajectory", {"type": str, "metavar": "PATH",
        "help": "Trajectory file (default: simulation/production.dcd)."}),
    ("topology", "topology", {"type": str, "metavar": "PATH",
        "help": "Topology file (default: simulation/topology.pdb)."}),
    ("analyses", "include", {"nargs": "+", "metavar": "NAME",
        "help": "Subset of analyses to run (e.g. rmsd rmsf rg). Default: all."}),
    ("exclude-analyses", "exclude", {"nargs": "+", "metavar": "NAME",
        "help": "Analyses to skip. Mutually exclusive with --analyses."}),
    ("selection", "selection", {"type": str, "metavar": "EXPR",
        "help": "Default MDTraj atom selection (e.g. 'name CA'). Overrides --scope."}),
    ("scope", "scope", {"choices": ["solute", "protein", "ligand", "all"],
        "help": "Atom scope for analyses (default 'solute' = protein+ligand)."}),
    ("stride", "stride", {"type": int,
        "help": "Frame stride for trajectory loading (default 1)."}),
    ("first", "first", {"type": int,
        "help": "First frame index to include (default 0)."}),
    ("last", "last", {"type": int,
        "help": "Last frame index (exclusive). Default: full trajectory."}),
]

_REPORT_OPTIONS: list[tuple[str, str, dict[str, Any]]] = [
    ("title", "title", {"type": str, "metavar": "STR",
        "help": "Report title (default auto-generated from system name)."}),
    ("author", "author", {"type": str, "metavar": "NAME",
        "help": "Author name for the report metadata."}),
    ("no-document", "document", {"action": "store_false", "default": None,
        "help": "Skip the Markdown document."}),
    ("no-slides", "slides", {"action": "store_false", "default": None,
        "help": "Skip the PPTX slide deck."}),
    ("no-bundle", "bundle", {"action": "store_false", "default": None,
        "help": "Skip the project_bundle.zip artifact."}),
    ("no-methods", "include_methods", {"action": "store_false", "default": None,
        "help": "Skip the Methods section in the document."}),
    ("no-reproducibility", "include_reproducibility", {"action": "store_false", "default": None,
        "help": "Skip the Reproducibility section."}),
]


# Map: phase -> (options-list, explore-prefix)
_PHASE_SPEC = {
    "setup":      (_SETUP_OPTIONS,      "setup"),
    "simulate":   (_SIMULATION_OPTIONS, "simulate"),
    "analyze":    (_ANALYSIS_OPTIONS,   "analyze"),
    "report":     (_REPORT_OPTIONS,     "report"),
}

# Phase-name aliases used when forwarding options to the orchestrator's
# `options=` dict. CLI says "simulate" / "analyze" but the orchestrator
# uses "simulation" / "analysis".
_PHASE_TO_ORCH = {
    "setup":    "setup",
    "simulate": "simulation",
    "analyze":  "analysis",
    "report":   "report",
}


def _attach_phase_options(
    parser: argparse.ArgumentParser,
    options: list[tuple[str, str, dict[str, Any]]],
    *,
    prefix: str = "",
    dest_prefix: str = "",
    group_title: str = "options",
) -> None:
    """Attach a phase's options to a parser under an argparse group.

    Parameters
    ----------
    parser : argparse.ArgumentParser
    options : list
        The per-phase option tuples (cli_suffix, kwarg, argparse_kwargs).
    prefix : str
        Prepended to the CLI flag with a dash (e.g. ``prefix="simulate"``
        → ``--simulate-duration-ns``). Empty for per-phase subcommands.
    dest_prefix : str
        Prepended to ``argparse.dest`` to avoid name collisions in
        ``explore``. Same convention as ``prefix`` but with underscores.
    group_title : str
        Title for the argument group (shown in --help).
    """
    group = parser.add_argument_group(group_title)
    for cli_suffix, kwarg, argparse_kwargs in options:
        if prefix:
            flag = f"--{prefix}-{cli_suffix}"
            dest = f"{dest_prefix}__{kwarg}"
        else:
            flag = f"--{cli_suffix}"
            dest = kwarg
        group.add_argument(flag, dest=dest, **argparse_kwargs)


def _harvest_phase_options(
    args: argparse.Namespace,
    options: list[tuple[str, str, dict[str, Any]]],
    *,
    dest_prefix: str = "",
) -> dict[str, Any]:
    """Pull phase options out of parsed args, dropping None values.

    Returns a kwargs-shaped dict ready to splat into the phase's run().
    None values are dropped so the phase falls back to its own DEFAULTS
    table — important because argparse uses None for unset args even
    when the phase's own default is something else.
    """
    out: dict[str, Any] = {}
    for _cli_suffix, kwarg, _ in options:
        dest = f"{dest_prefix}__{kwarg}" if dest_prefix else kwarg
        val = getattr(args, dest, None)
        if val is not None:
            out[kwarg] = val
    return out


# ---------------------------------------------------------------------------
# Parser construction
# ---------------------------------------------------------------------------
def _common_input_args(p: argparse.ArgumentParser) -> None:
    """Arguments shared by all subcommands that accept a system input.

    The ``system`` flag accepts three forms: ``-s`` (GNU short option),
    ``-system`` (single-dash long, the GROMACS / AMBER / NAMD convention
    MD researchers expect), and ``--system`` (GNU double-dash long). All
    three are equivalent.

    The system value is auto-classified downstream: a path ending in
    ``.pdb`` / ``.cif`` is loaded from disk, a 4-character alphanumeric
    string is fetched from RCSB as a PDB ID, and a longer alphabetic
    string is treated as a one-letter sequence. There is therefore no
    separate ``--pdb-id`` flag — ``--system 1L2Y`` does the right thing.
    """
    src = p.add_argument_group("input")
    src.add_argument(
        "-s", "-system", "--system",
        dest="system",
        metavar="SYSTEM",
        help=(
            "System input: a PDB/CIF file path, a 4-character PDB ID "
            "(e.g. 1L2Y, fetched from RCSB), or a one-letter sequence. "
            "May instead be supplied via --config."
        ),
    )
    src.add_argument(
        "-c", "-config", "--config",
        dest="config",
        metavar="FILE",
        help=(
            "YAML config file capturing the whole run (system, output, "
            "phase selection, per-phase options). Command-line flags "
            "override values in the file. See `fastmdx init-config`."
        ),
    )
    src.add_argument(
        "--output",
        dest="output_dir",
        metavar="DIR",
        help=(
            "Output directory for project artifacts "
            "(default: ./fastmdxplora_output_<UTC-timestamp>)."
        ),
    )
    src.add_argument(
        "--verbose",
        action="store_true",
        help="Also stream debug logging to the terminal.",
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="fastmdx",
        description=(
            f"FastMDXplora: {__expansion__}\n\n"
            "Project-level orchestrator for end-to-end molecular dynamics "
            "studies: setup → simulate → analyze → report."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  fastmdx explore -system protein.pdb\n"
            "  fastmdx xplore --system 1L2Y --simulate-duration-ns 50.0\n"
            "  fastmdx setup -system protein.pdb --ph 6.5\n"
            "  fastmdx simulate --duration-ns 100.0 --platform CUDA\n"
            "  fastmdx analyze --output run_001 --analyses rmsd rmsf rg\n"
            "  fastmdx report --output run_001 --no-slides\n"
            "\n"
            f"Citation: {__citation__}\n"
        ),
    )
    parser.add_argument(
        "--version",
        "-V",
        action="version",
        version=f"fastmdx {__version__} (FastMDXplora)",
    )
    parser.add_argument(
        "--cite",
        action="store_true",
        help="Print the citation and exit.",
    )

    sub = parser.add_subparsers(
        dest="command",
        metavar="<command>",
        title="commands",
    )

    # ---------- explore / xplore: full pipeline with prefixed flags ----------
    for verb in ("explore", "xplore"):
        ep = sub.add_parser(
            verb,
            help=(
                "Run the full pipeline (setup → simulate → analyze → report)."
                if verb == "explore"
                else "Alias of `explore` (matches the X branding)."
            ),
            description=(
                "Run the full FastMDXplora pipeline end-to-end on the given "
                "system. Use --include or --exclude to run a subset of phases. "
                "Each phase's flags are available under a per-phase prefix "
                "(--setup-*, --simulate-*, --analyze-*, --report-*)."
            ),
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        _common_input_args(ep)
        ep.add_argument(
            "--include",
            nargs="+",
            metavar="PHASE",
            help="Subset of phases to run: setup, simulation, analysis, report.",
        )
        ep.add_argument(
            "--exclude",
            nargs="+",
            metavar="PHASE",
            help="Phases to skip (mutually exclusive with --include).",
        )
        ep.add_argument(
            "--no-report",
            action="store_true",
            help="Skip the report phase even if it would otherwise run.",
        )
        ep.add_argument(
            "--dry-run",
            action="store_true",
            help="Print the plan (runs, systems, swept values, output dirs, "
                 "phases) without running anything.",
        )

        # Per-phase options under per-phase prefix
        for phase, (opts, prefix) in _PHASE_SPEC.items():
            _attach_phase_options(
                ep, opts,
                prefix=prefix,
                dest_prefix=phase,
                group_title=f"{phase} options",
            )

    # ---------- per-phase subcommands: phase-specific flags only ----------
    for phase, (opts, _) in _PHASE_SPEC.items():
        pp = sub.add_parser(
            phase,
            help=f"Run only the {phase} phase.",
            description=f"Run only the {phase} phase of the FastMDXplora pipeline.",
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        _common_input_args(pp)
        _attach_phase_options(pp, opts, group_title=f"{phase} options")

    sub.add_parser(
        "info",
        help="Print FastMDXplora environment information.",
        description=(
            "Print the installed FastMDXplora version, the detected backends "
            "for each phase, and the citation."
        ),
    )
    splash = sub.add_parser(
        "splash",
        help="Show or open the bundled FastMDXplora graphical abstract.",
        description=(
            "Print the terminal startup splash and optionally open the full "
            "graphical abstract image bundled with FastMDXplora."
        ),
    )
    splash.add_argument(
        "--open",
        action="store_true",
        help="Open the full graphical abstract image with the default viewer.",
    )
    splash.add_argument(
        "--path",
        action="store_true",
        help="Print the local path to the bundled graphical abstract image.",
    )
    splash.add_argument(
        "--no-preview",
        action="store_true",
        help="Do not print the terminal splash preview.",
    )
    # init-config: write a commented YAML template
    ic = sub.add_parser(
        "init-config",
        help="Write a commented YAML config template to edit.",
        description=(
            "Generate a FastMDXplora config template. By default writes a "
            "comprehensive, fully-commented template with every option, its "
            "default, and a description. Edit it and run with "
            "`fastmdx explore --config <file>`."
        ),
    )
    ic.add_argument(
        "-o", "--output",
        dest="config_output",
        metavar="FILE",
        default="fastmdxplora.yml",
        help="Where to write the template (default: fastmdxplora.yml).",
    )
    ic.add_argument(
        "--minimal",
        action="store_true",
        help="Write a short starter template with only the essentials.",
    )
    ic.add_argument(
        "--force",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )

    return parser


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------
def _infer_system_from_output(output_dir: str | None) -> str | None:
    """Best-effort system inference for report/analyze reruns on existing output."""
    if not output_dir:
        return None

    import json

    root = Path(output_dir)
    manifest = root / "manifest.json"
    try:
        data = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        data = {}
    system = data.get("system")
    if system:
        return str(system)

    topology = root / "simulation" / "topology.pdb"
    if topology.exists():
        return str(topology)
    return None


def _make_orchestrator(args: argparse.Namespace, *, phase: str | None = None) -> FastMDXplora:
    """Build a single-system orchestrator for the per-phase subcommands.

    The per-phase commands (setup/simulate/analyze/report) operate on one
    system directly, so they bypass the batch layer. `explore` always goes
    through BatchExplorer instead.
    """
    config = getattr(args, "config", None)
    inferred_system = (
        _infer_system_from_output(args.output_dir)
        if phase in {"analyze", "report"} else None
    )
    if not args.system and not config and not inferred_system:
        raise SystemExit(
            "fastmdx: this command requires a system input "
            "(-s / -system / --system) or a --config file."
        )
    # For per-phase commands with a config file, pull the first system out.
    if config and not args.system:
        from fastmdxplora.config import load_config_file
        from fastmdxplora.batch.sweep import normalize_systems

        raw = load_config_file(config)
        systems = normalize_systems(raw.get("systems") or [])
        system = systems[0]["system"]
    else:
        system = args.system or inferred_system
    return FastMDXplora(
        system=system,
        output_dir=args.output_dir,
        verbose=args.verbose,
    )


def _build_explore_config(args: argparse.Namespace) -> dict[str, Any]:
    """Assemble the config dict that drives an `explore` run.

    Two sources, in priority order (flags win):
      1. A ``--config`` YAML file (if given).
      2. Command-line flags: ``-s/--system`` builds a one-element
         ``systems`` list; per-phase prefixed flags become phase blocks.

    The result always has a ``systems`` list, so it flows through
    BatchExplorer like any other config (a single system is a batch of
    one, written with the flat output layout).
    """
    from fastmdxplora.config import load_config_file

    # Start from the file, if any.
    if getattr(args, "config", None):
        config = load_config_file(args.config)
    else:
        config = {}

    # Harvest per-phase option flags and merge them on top (flags win).
    for phase, (opts, _prefix) in _PHASE_SPEC.items():
        harvested = _harvest_phase_options(args, opts, dest_prefix=phase)
        if harvested:
            orch_phase = _PHASE_TO_ORCH[phase]
            block = dict(config.get(orch_phase, {}))
            block.update(harvested)
            config[orch_phase] = block

    # The flat --simulate-plumed-script flag maps to the nested `plumed` dict.
    sim_block = config.get("simulation")
    if isinstance(sim_block, dict) and "plumed_script" in sim_block:
        script = sim_block.pop("plumed_script")
        if script:
            sim_block["plumed"] = {"enabled": True, "script": script}

    # -s/--system builds (or replaces) a one-element systems list.
    if args.system:
        config["systems"] = [{"id": "s1", "system": args.system}]

    # Top-level scalars from flags.
    if args.output_dir:
        config["output"] = args.output_dir
    if getattr(args, "verbose", False):
        config["verbose"] = True
    if args.include:
        config["include"] = args.include
    if args.exclude:
        config["exclude"] = args.exclude

    return config


def _cmd_explore(args: argparse.Namespace) -> int:
    from fastmdxplora import FastMDXplora

    if args.include and args.exclude:
        print("fastmdx: --include and --exclude are mutually exclusive.", file=sys.stderr)
        return 2

    config = _build_explore_config(args)
    if not config.get("systems"):
        print(
            "fastmdx: explore requires a system — pass -s/--system PATH or a "
            "--config file with a `systems:` list.",
            file=sys.stderr,
        )
        return 2

    # --no-report removes report from the plan via exclude (unless the user
    # already constrained phases with include).
    if getattr(args, "no_report", False) and not config.get("include"):
        existing = config.get("exclude") or []
        if "report" not in existing:
            config["exclude"] = [*existing, "report"]

    fmdx = FastMDXplora(
        config_data=config,
        output_dir=args.output_dir,
        verbose=args.verbose,
    )
    results = fmdx.explore(dry_run=getattr(args, "dry_run", False))

    # Dry run: the plan was printed; nothing executed.
    if getattr(args, "dry_run", False):
        return 0

    # Single run -> flat layout; point at the project manifest.
    if len(results) == 1:
        print()
        print(f"Project output: {fmdx.output_dir}")
        print(f"Manifest:       {fmdx.output_dir / 'manifest.json'}")
    return 0 if all(r.status == "ok" for r in results) else 1


def _cmd_phase(phase: str, args: argparse.Namespace) -> int:
    fmdx = _make_orchestrator(args, phase=phase)
    opts_list, _ = _PHASE_SPEC[phase]
    kwargs = _harvest_phase_options(args, opts_list)

    method = {
        "setup":    fmdx.setup,
        "simulate": fmdx.simulate,
        "analyze":  fmdx.analyze,
        "report":   fmdx.report,
    }[phase]

    # Bracket the single-phase invocation with presenter output so the
    # user sees the same visual structure as during `fastmdx explore`.
    fmdx._presenter.phase_start(phase)  # noqa: SLF001 -- internal hook
    result = method(**kwargs)
    fmdx._presenter.phase_end(phase, status=result.status)
    fmdx._write_manifest()  # noqa: SLF001 -- single-phase still records
    print()
    print(f"Project output: {fmdx.output_dir}")
    return 0 if result.status == "ok" else 1


def _cmd_info() -> int:
    print(f"FastMDXplora {__version__}")
    print(f"  {__expansion__}")
    print(f"  Authors: {__author__}")
    print(f"  DOI:     {__doi__}")
    print()
    print("Phases:")
    for name in ("setup", "simulation", "analysis", "report"):
        try:
            module = __import__(f"fastmdxplora.{name}", fromlist=["run"])
            run_fn = getattr(module, "run", None)
            status = "available" if callable(run_fn) else "missing run()"
        except Exception as exc:  # noqa: BLE001
            status = f"import error: {exc}"
        print(f"  {name:<11} {status}")
    print()
    print("Optional backends (real chemistry):")
    for display_name, import_name, install_hint in (
        ("PDBFixer", "pdbfixer",
            "conda install -c conda-forge pdbfixer"),
        ("OpenMM",   "openmm",
            "conda install -c conda-forge openmm"),
    ):
        try:
            __import__(import_name)
            print(f"  {display_name:<10} installed")
        except ImportError:
            print(f"  {display_name:<10} not installed  ({install_hint})")
    print()
    print(f"Citation: {__citation__}")
    return 0

def _cmd_splash(args: argparse.Namespace) -> int:
    from fastmdxplora.utils.splash import (
        graphical_abstract_path,
        open_graphical_abstract,
        print_startup_splash,
    )

    path = graphical_abstract_path()

    if not getattr(args, "no_preview", False):
        print_startup_splash(show_asset_hint=False)

    if getattr(args, "path", False):
        print(path)

    if getattr(args, "open", False):
        if not path.exists():
            print(f"fastmdx: graphical abstract not found: {path}", file=sys.stderr)
            return 1

        opened = open_graphical_abstract()
        if not opened:
            print(f"fastmdx: could not open graphical abstract: {path}", file=sys.stderr)
            return 1

        print(f"Opened graphical abstract: {path}")

    elif not getattr(args, "path", False):
        print(f"Graphical abstract: {path}")
        print("Open it with: fastmdx splash --open")

    return 0

def _cmd_init_config(args: argparse.Namespace) -> int:
    from fastmdxplora.config import generate_template

    out_path = Path(args.config_output)
    if out_path.exists() and not args.force:
        print(
            f"fastmdx: {out_path} already exists. Use --force to overwrite, "
            f"or -o to choose a different path.",
            file=sys.stderr,
        )
        return 2

    text = generate_template(minimal=args.minimal)
    out_path.write_text(text, encoding="utf-8")
    kind = "minimal" if args.minimal else "comprehensive"
    print(f"Wrote {kind} config template to {out_path}")
    print(f"Edit it, then run:  fastmdx explore --config {out_path}")
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def main(argv: Sequence[str] | None = None) -> int:
    # Ensure the CLI can emit its Unicode output (box-drawing banner, "→",
    # "—") regardless of the platform's locale. On machines whose default
    # stdio encoding is ASCII, printing these would otherwise raise
    # UnicodeEncodeError. reconfigure() is available on Python 3.7+ text
    # streams; guard for unusual stream types.
    for _stream in (sys.stdout, sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8")  # type: ignore[union-attr]
        except (AttributeError, ValueError):
            pass

    # Initialize console logging on every CLI invocation. setup_console() is
    # idempotent (no duplicate handlers) and honors FASTMDX_LOG_STYLE /
    # FASTMDX_LOGLEVEL / NO_COLOR.
    from fastmdxplora.utils.logging import setup_console

    setup_console()

    parser = _build_parser()
    args = parser.parse_args(argv)
    command = getattr(args, "command", None)
    include_phases = set(getattr(args, "include", None) or [])
    exclude_phases = set(getattr(args, "exclude", None) or [])
    dry_run = bool(getattr(args, "dry_run", False))

    simulation_requested = (
        command == "simulate"
        or (
            command in {"explore", "xplore"}
            and (
                "simulation" in include_phases
                or (not include_phases and "simulation" not in exclude_phases)
            )
        )
    )

    if simulation_requested and not dry_run:
        try:
            from fastmdxplora.utils.splash import print_startup_splash

            print_startup_splash()
        except Exception:
            pass
    if args.cite:
        print(__citation__)
        return 0

    if args.command is None:
        parser.print_help()
        return 0

    if args.command == "init-config":
        return _cmd_init_config(args)

    if args.command == "splash":
        return _cmd_splash(args)

    # Commands that build an orchestrator can hit config-file errors;
    # surface those cleanly rather than as a traceback.
    from fastmdxplora.config import ConfigError

    try:
        if args.command in ("explore", "xplore"):
            return _cmd_explore(args)
        if args.command in ("setup", "simulate", "analyze", "report"):
            return _cmd_phase(args.command, args)
        if args.command == "info":
            return _cmd_info()
    except ConfigError as exc:
        print(f"fastmdx: config error: {exc}", file=sys.stderr)
        return 2

    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
