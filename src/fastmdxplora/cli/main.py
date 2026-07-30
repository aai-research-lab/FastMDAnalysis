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
import os
import sys
from datetime import datetime, timezone
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
    ("live-telemetry", "live_telemetry", {"action": "store_true", "default": None,
        "help": "Write lightweight live dashboard telemetry during simulation."}),
    ("telemetry-interval", "telemetry_interval", {"type": int,
        "help": "Minimum step interval for live telemetry updates (default 1000)."}),
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
    ("compat", "compat", {"choices": ["fastmdanalysis"],
        "help": "Use a compatibility profile for a published FastMDAnalysis workflow."}),
    ("dimred-methods", "dimred_methods", {"nargs": "+", "choices": ["pca", "tsne", "umap"],
        "metavar": "METHOD",
        "help": "Dimensionality-reduction methods (e.g. pca)."}),
    ("dimred-components", "dimred_components", {"type": int, "metavar": "N",
        "help": "Number of dimensionality-reduction components (default 2)."}),
    ("cluster-methods", "cluster_methods", {"nargs": "+",
        "choices": ["kmeans", "hierarchical", "dbscan"], "metavar": "METHOD",
        "help": "Clustering methods (e.g. hierarchical)."}),
    ("cluster-n-clusters", "cluster_n_clusters", {"type": int, "metavar": "N",
        "help": "Number of clusters for k-means/hierarchical clustering."}),
    ("cluster-linkage", "cluster_linkage",
        {"choices": ["ward", "complete", "average", "single"],
         "help": "Hierarchical clustering linkage method."}),
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


def _normalize_analysis_options(kwargs: dict[str, Any]) -> dict[str, Any]:
    """Move CLI analysis-method flags into the nested ``options`` mapping.

    The analysis orchestrator accepts method-specific settings under
    ``options`` (for example ``options.cluster.n_clusters``), while argparse
    naturally harvests flat flags. Keeping this conversion at the CLI
    boundary makes the same options work for both ``analyze`` and prefixed
    ``explore`` flags without changing the Python API.
    """
    out = dict(kwargs)
    compat = out.pop("compat", None)
    if compat == "fastmdanalysis":
        # Settings used by the paper's BPTI case study. Explicit CLI values
        # still win, so this is a reproducible opt-in profile rather than a
        # change to FastMDXplora's general-purpose defaults.
        out.setdefault("scope", "protein")
        # FastMDAnalysis uses the complete protein atom set when the BPTI
        # options file does not provide an atom selection.  Set this at the
        # compatibility boundary so RMSD/RMSF do not fall back to Xplora's
        # CA-only defaults.
        out.setdefault("selection", "protein")
        out.setdefault("stride", 2)
        out.setdefault(
            "include",
            ["rmsd", "rmsf", "rg", "hbonds", "sasa", "ss", "dimred", "cluster"],
        )
    nested = {
        "dimred": {
            key: out.pop(key)
            for key in ("dimred_methods", "dimred_components")
            if key in out
        },
        "cluster": {
            key: out.pop(key)
            for key in ("cluster_methods", "cluster_n_clusters", "cluster_linkage")
            if key in out
        },
    }
    options = dict(out.get("options") or {})
    if compat == "fastmdanalysis":
        paper_defaults = {
            "rmsd": {"ref": 0, "align": True},
            "rmsf": {"ref": 0, "per_residue": False},
            "sasa": {
                "mode": "total",
                "probe_radius": 0.14,
                "n_sphere_points": 960,
            },
            # FastMDAnalysis reports both donor/acceptor directions and
            # counts transient Baker-Hubbard bonds, rather than filtering
            # the per-frame series by the occupancy threshold first.
            "hbonds": {
                "method": "baker_hubbard",
                "freq": 0.1,
                "candidate_freq": 0.0,
                "count_multiplier": 2,
            },
            "dimred": {"methods": ["pca"], "n_components": 2},
            "cluster": {"methods": ["hierarchical"], "n_clusters": 6},
        }
        for name, values in paper_defaults.items():
            current = dict(options.get(name) or {})
            for key, value in values.items():
                current.setdefault(key, value)
            options[name] = current
    for name, values in nested.items():
        if values:
            current = dict(options.get(name) or {})
            current.update(
                {key.removeprefix(f"{name}_"): value for key, value in values.items()}
            )
            options[name] = current
    if options:
        out["options"] = options
    elif "options" in out:
        out.pop("options")
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
    dash = p.add_argument_group("dashboard")
    dash.add_argument(
        "--dashboard",
        "--live-dashboard",
        dest="dashboard",
        action="store_true",
        default=False,
        help=(
            "Launch the local live dashboard for this output folder before "
            "the workflow starts. Implies live telemetry when simulation runs."
        ),
    )
    dash.add_argument(
        "--dashboard-host",
        default="127.0.0.1",
        help="Dashboard bind address (default: 127.0.0.1).",
    )
    dash.add_argument(
        "--dashboard-port",
        type=int,
        default=8765,
        help="Dashboard port (default: 8765; next free port is used if busy).",
    )
    dash.add_argument(
        "--dashboard-stop-on-complete",
        action="store_true",
        default=False,
        help="Stop the dashboard automatically when the command completes.",
    )
    dash.add_argument(
        "--dashboard-refresh-seconds",
        type=float,
        default=None,
        metavar="SECONDS",
        help="Browser-side telemetry polling interval in seconds (default 3).",
    )
    dash.add_argument(
        "--dashboard-frame-interval",
        type=int,
        default=None,
        metavar="STEPS",
        help="Override simulation telemetry interval used by the live dashboard. "
             "Honored when the workflow is creating a telemetry writer; existing "
             "runs keep their stored value.",
    )
    dash.add_argument(
        "--dashboard-ligand-resname",
        type=str,
        default=None,
        metavar="RESNAME",
        help="Force a ligand residue name for the dashboard ligand tools pane. "
             "Auto-detection is used when omitted.",
    )
    dash.add_argument(
        "--dashboard-binding-pocket-cutoff-A",
        type=float,
        default=None,
        metavar="ANGSTROM",
        help="Default binding-pocket cutoff for the molecular viewer (default 5.0).",
    )
    dash.add_argument(
        "--dashboard-max-playback-frames",
        type=int,
        default=None,
        metavar="FRAMES",
        help="Maximum number of frames the molecular viewer will load for "
             "trajectory playback (default 200).",
    )
    dash.add_argument(
        "--dashboard-open-browser",
        action="store_true",
        default=False,
        help="Attempt to open the dashboard URL in the local browser. "
             "Disabled by default for headless / no-display environments.",
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

    health = sub.add_parser(
        "health",
        help="Run repository health checks and environment diagnostics.",
        description=(
            "Run the repository doctor checks to validate the local checkout, "
            "environment, and package readiness."
        ),
    )
    health.add_argument(
        "--no-fix",
        action="store_true",
        help="Only diagnose problems; do not install or modify anything.",
    )
    health.add_argument(
        "--yes",
        action="store_true",
        help="Accept all fixes automatically.",
    )

    sub.add_parser(
        "info",
        help="Print FastMDXplora environment information.",
        description=(
            "Print the installed FastMDXplora version, the detected backends "
            "for each phase, and the citation."
        ),
    )

    dashboard = sub.add_parser(
        "dashboard",
        help="Serve local dashboard views for an existing output directory.",
        description=(
            "Serve a local-only dashboard for completed outputs and live "
            "simulation telemetry. Binds to 127.0.0.1 by default."
        ),
    )
    dashboard_sub = dashboard.add_subparsers(dest="dashboard_command", metavar="<dashboard-command>")
    serve = dashboard_sub.add_parser(
        "serve",
        help="Serve the live dashboard for an output directory.",
        description="Start the local dashboard server for an existing FastMDXplora output.",
    )
    serve.add_argument(
        "--output",
        required=True,
        metavar="DIR",
        help="FastMDXplora output directory to watch.",
    )
    serve.add_argument(
        "--host",
        default="127.0.0.1",
        help="Bind address (default: 127.0.0.1).",
    )
    serve.add_argument(
        "--port",
        type=int,
        default=8765,
        help="Port to serve on (default: 8765).",
    )
    serve.add_argument(
        "--refresh-seconds",
        type=float,
        default=None,
        metavar="SECONDS",
        help="Browser polling interval hint surfaced to the dashboard.",
    )
    serve.add_argument(
        "--ligand-resname",
        type=str,
        default=None,
        metavar="RESNAME",
        help="Force a ligand residue name for the dashboard ligand tools.",
    )
    serve.add_argument(
        "--binding-pocket-cutoff-A",
        type=float,
        default=None,
        metavar="ANGSTROM",
        help="Binding-pocket cutoff in angstrom (default 5.0).",
    )
    serve.add_argument(
        "--frame-interval",
        type=int,
        default=None,
        metavar="STEPS",
        help="Simulation telemetry interval used by the live dashboard.",
    )
    serve.add_argument(
        "--max-playback-frames",
        type=int,
        default=None,
        metavar="FRAMES",
        help="Maximum frames the molecular viewer will load for playback.",
    )
    serve.add_argument(
        "--open-browser",
        action="store_true",
        default=False,
        help="Open the dashboard URL in the user's default browser.",
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

    def _add_bootstrap_parser(name: str, *, description: str, help_text: str) -> argparse.ArgumentParser:
        parser_obj = sub.add_parser(
            name,
            help=help_text,
            description=description,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )
        parser_obj.add_argument(
            "--env-name",
            default="fastmdxplora",
            help="Conda environment name to create (default: fastmdxplora).",
        )
        parser_obj.add_argument(
            "--python-version",
            default="3.10",
            help="Python version to install in the environment (3.9-3.12).",
        )
        parser_obj.add_argument(
            "--force",
            action="store_true",
            help="Recreate the environment if it already exists.",
        )
        parser_obj.add_argument(
            "--yes",
            "-y",
            action="store_true",
            help="Skip confirmation prompts if any.",
        )
        return parser_obj

    _add_bootstrap_parser(
        "install",
        description=(
            "Create a conda environment and install FastMDXplora. If you run it "
            "inside a repository checkout, the local checkout is installed; "
            "otherwise the published package is installed."
        ),
        help_text="Install FastMDXplora into a conda environment.",
    )
    _add_bootstrap_parser(
        "bootstrap",
        description=(
            "Alias for `install`; kept for compatibility with older docs and scripts."
        ),
        help_text="Alias for install.",
    )
    _add_bootstrap_parser(
        "install-e",
        description=(
            "Create a conda environment and install the current repository checkout "
            "in editable mode for contributors and editors."
        ),
        help_text="Install the current repository checkout in editable mode.",
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
        if phase == "analyze":
            harvested = _normalize_analysis_options(harvested)
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


def _dashboard_requested(args: argparse.Namespace) -> bool:
    return bool(getattr(args, "dashboard", False))


def _enable_dashboard_telemetry(
    config: dict[str, Any], args: argparse.Namespace | None = None
) -> None:
    simulation = dict(config.get("simulation", {}))
    simulation["live_telemetry"] = True
    if args is not None and getattr(args, "dashboard_frame_interval", None) is not None:
        simulation["telemetry_interval"] = int(args.dashboard_frame_interval)
    config["simulation"] = simulation


def _resolve_dashboard_output_dir(args: argparse.Namespace, config: dict[str, Any] | None = None) -> Path:
    raw_output = getattr(args, "output_dir", None)
    if not raw_output and config:
        raw_output = config.get("output")
    if not raw_output:
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        raw_output = f"fastmdxplora_output_{timestamp}"
    return Path(raw_output).expanduser().resolve()


def _start_dashboard_for_command(args: argparse.Namespace, output_dir: Path):
    # The orchestrator uses this process-local marker to publish setup,
    # analysis, and report phase transitions to the same live timeline as
    # the OpenMM simulation sub-stages.
    os.environ["FASTMDX_DASHBOARD_ACTIVE"] = "1"
    os.environ["FASTMDX_DASHBOARD_OUTPUT"] = str(output_dir)

    from fastmdxplora.live.server import (
        DashboardConfig,
        start_dashboard_session,
    )

    config = DashboardConfig(
        ligand_resname=getattr(args, "dashboard_ligand_resname", None),
        binding_pocket_cutoff_A=float(
            getattr(args, "dashboard_binding_pocket_cutoff_A", 5.0) or 5.0
        ),
        max_browser_frames=int(
            getattr(args, "dashboard_max_playback_frames", 200) or 200
        ),
        refresh_seconds=float(
            getattr(args, "dashboard_refresh_seconds", 3.0) or 3.0
        ),
    )
    session = start_dashboard_session(
        output=output_dir,
        host=args.dashboard_host,
        port=args.dashboard_port,
        config=config,
    )
    print(f"Live dashboard running at: {session.url}")
    if session.port_was_changed:
        print(
            f"Requested port {session.requested_port} was busy, "
            f"so FastMDXplora used {session.port}."
        )
    if args.dashboard_host == "0.0.0.0":
        print(
            "Warning: dashboard is bound to 0.0.0.0 and may be visible on your network."
        )
        print("Use --dashboard-host 127.0.0.1 for local-only access.")
    print(f"Watching output folder: {output_dir}")
    print("Open this URL in your browser to monitor the run.")
    if args.dashboard_stop_on_complete:
        print("Dashboard will stop automatically after the workflow completes.")
    else:
        print("Press Ctrl+C to stop the dashboard after the workflow completes.")
    print()
    return session


def _finish_dashboard_for_command(session, args: argparse.Namespace) -> None:
    if session is None:
        return
    if args.dashboard_stop_on_complete:
        session.stop()
        return
    print()
    print(f"Workflow complete. Live dashboard is still running at: {session.url}")
    print("Press Ctrl+C to stop the dashboard.")
    try:
        session.wait_forever()
    except KeyboardInterrupt:
        pass
    finally:
        session.stop()


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

    if _dashboard_requested(args):
        _enable_dashboard_telemetry(config, args)
    dashboard_output_dir: Path | None = None
    if _dashboard_requested(args):
        dashboard_output_dir = _resolve_dashboard_output_dir(args, config)
        config["output"] = str(dashboard_output_dir)

    fmdx = FastMDXplora(
        config_data=config,
        output_dir=args.output_dir,
        verbose=args.verbose,
    )
    session = None
    if _dashboard_requested(args) and not getattr(args, "dry_run", False):
        session = _start_dashboard_for_command(args, dashboard_output_dir)
    try:
        results = fmdx.explore(dry_run=getattr(args, "dry_run", False))
    except KeyboardInterrupt:
        if session is not None:
            session.stop()
        return 130
    except Exception:
        if session is not None:
            session.stop()
        raise

    # Dry run: the plan was printed; nothing executed.
    if getattr(args, "dry_run", False):
        return 0

    # Single run -> flat layout; point at the project manifest.
    if len(results) == 1:
        print()
        print(f"Project output: {fmdx.output_dir}")
        print(f"Manifest:       {fmdx.output_dir / 'manifest.json'}")
    rc = 0 if all(r.status == "ok" for r in results) else 1
    _finish_dashboard_for_command(session, args)
    return rc


def _cmd_phase(phase: str, args: argparse.Namespace) -> int:
    fmdx = _make_orchestrator(args, phase=phase)
    opts_list, _ = _PHASE_SPEC[phase]
    kwargs = _harvest_phase_options(args, opts_list)
    if phase == "analyze":
        kwargs = _normalize_analysis_options(kwargs)
    if _dashboard_requested(args) and phase == "simulate":
        kwargs["live_telemetry"] = True
        # Forward dashboard knobs when running live; ignored if the user
        # did not opt in to live telemetry.
        if getattr(args, "dashboard_frame_interval", None) is not None:
            kwargs["telemetry_interval"] = int(args.dashboard_frame_interval)
        if getattr(args, "dashboard_refresh_seconds", None) is not None:
            # The same value is embedded into the served dashboard HTML.
            print(
                f"  dashboard polling: every {args.dashboard_refresh_seconds}s"
            )

    method = {
        "setup":    fmdx.setup,
        "simulate": fmdx.simulate,
        "analyze":  fmdx.analyze,
        "report":   fmdx.report,
    }[phase]

    # Bracket the single-phase invocation with presenter output so the
    # user sees the same visual structure as during `fastmdx explore`.
    session = None
    if _dashboard_requested(args):
        output_dir = _resolve_dashboard_output_dir(args)
        if not getattr(args, "output_dir", None):
            output_dir = Path(fmdx.output_dir).expanduser().resolve()
        session = _start_dashboard_for_command(args, output_dir)
    try:
        fmdx._presenter.phase_start(phase)  # noqa: SLF001 -- internal hook
        result = method(**kwargs)
        fmdx._presenter.phase_end(phase, status=result.status)
        fmdx._write_manifest()  # noqa: SLF001 -- single-phase still records
    except KeyboardInterrupt:
        if session is not None:
            session.stop()
        return 130
    except Exception:
        if session is not None:
            session.stop()
        raise
    print()
    print(f"Project output: {fmdx.output_dir}")
    rc = 0 if result.status == "ok" else 1
    _finish_dashboard_for_command(session, args)
    return rc


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


def _cmd_bootstrap(args: argparse.Namespace, *, editable: bool = False, package_name: str = "fastmdxplora") -> int:
    from fastmdxplora.install import bootstrap_environment, BootstrapError

    repo_root = Path.cwd()
    repo_marker = (repo_root / "pyproject.toml").exists() and (repo_root / "src" / "fastmdxplora").exists()
    resolved_package_name = "." if package_name == "." and repo_marker else package_name
    resolved_editable = editable and repo_marker

    try:
        bootstrap_environment(
            env_name=args.env_name,
            python_version=args.python_version,
            yes=args.yes,
            force=args.force,
            package_name=resolved_package_name,
            editable=resolved_editable,
        )
        return 0
    except BootstrapError as exc:
        print(f"fastmdx: bootstrap failed: {exc}", file=sys.stderr)
        return 1


def _cmd_health(args: argparse.Namespace) -> int:
    from health import main as health_main

    argv: list[str] = []
    if getattr(args, "no_fix", False):
        argv.append("--no-fix")
    if getattr(args, "yes", False):
        argv.append("--yes")
    return health_main(argv)


# ---------------------------------------------------------------------------
# Self-healing prologue
# ---------------------------------------------------------------------------
# `explore` / `xplore` / `setup` / `simulate` need OpenMM (and PDBFixer
# for setup). Without the chemistry stack installed, the user would
# otherwise get a long stack trace. Detect this here and offer a one-stop
# install hint instead, while still letting `info` / `health` / etc.
# proceed normally.
_CHEMISTRY_PHASES = frozenset({"setup", "simulation"})


def _needs_chemistry(args: argparse.Namespace) -> bool:
    """Return True if this command will actually invoke a chemistry phase."""
    cmd = getattr(args, "command", None)
    if cmd in ("setup", "simulate"):
        return True
    if cmd in ("explore", "xplore"):
        # Dry runs only print the plan; they must work without chemistry so
        # users can use them as a teaching tool when explaining the install gap.
        if getattr(args, "dry_run", False):
            return False
        include = getattr(args, "include", None)
        exclude = set(getattr(args, "exclude", None) or ())
        if include is not None:
            return bool(set(include) & _CHEMISTRY_PHASES)
        if exclude >= _CHEMISTRY_PHASES:
            return False
        return True  # default plan runs every phase
    return False


def _missing_chemistry_backends() -> list[str]:
    """Return the chemistry backend modules that aren't importable here.

    Probes the *actual* import shape used by setup and simulation
    (e.g. ``from openmm.app import PDBFile``) so a broken partial install
    fails fast here instead of mid-phase.
    """
    probes: tuple[tuple[str, str], ...] = (
        ("openmm", None),                                # top-level package
        ("openmm.app", "from openmm.app import PDBFile"),
        ("openmm", "from openmm import unit"),
        ("pdbfixer", "from pdbfixer import PDBFixer"),
    )
    failing: list[str] = []
    for name, stmt in probes:
        try:
            if stmt is None:
                __import__(name)
            else:
                exec(stmt, {})  # noqa: S102 — string intentional, gated by probes tuple
        except ImportError:
            failing.append(name)
    # Reduce to the *top-level* packages the user has to install, so the
    # hint stays short and actionable.
    return sorted({("openmm" if name.startswith("openmm") else name) for name in failing})


def _cmd_dashboard(args: argparse.Namespace) -> int:
    if args.dashboard_command != "serve":
        print("fastmdx: dashboard requires a subcommand, e.g. `dashboard serve`.", file=sys.stderr)
        return 2
    from fastmdxplora.live.server import DashboardConfig, serve_dashboard

    config = DashboardConfig(
        ligand_resname=getattr(args, "ligand_resname", None),
        binding_pocket_cutoff_A=float(
            getattr(args, "binding_pocket_cutoff_A", 5.0) or 5.0
        ),
        max_browser_frames=int(
            getattr(args, "max_playback_frames", 200) or 200
        ),
    )
    serve_dashboard(
        output=args.output,
        host=args.host,
        port=args.port,
        config=config,
    )
    if getattr(args, "open_browser", False):
        import webbrowser
        try:
            webbrowser.open(f"http://{args.host}:{args.port}", new=2)
        except Exception:
            pass
    return 0


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def _startup_dashboard_details(argv: Sequence[str]) -> tuple[str, bool]:
    """Resolve the dashboard address shown by the startup wordmark."""
    host = "127.0.0.1"
    port = "8765"
    enabled = (
        "--dashboard" in argv
        or "--live-dashboard" in argv
        or ("dashboard" in argv and "serve" in argv)
        or os.getenv("FASTMDX_DASHBOARD_ACTIVE") == "1"
    )

    for index, token in enumerate(argv):
        if token in {"--dashboard-host", "--host"} and index + 1 < len(argv):
            host = str(argv[index + 1])
        elif token.startswith("--dashboard-host=") or token.startswith("--host="):
            host = token.split("=", 1)[1]

        if token in {"--dashboard-port", "--port"} and index + 1 < len(argv):
            port = str(argv[index + 1])
        elif token.startswith("--dashboard-port=") or token.startswith("--port="):
            port = token.split("=", 1)[1]

    if host in {"0.0.0.0", "::", "[::]"}:
        host = "127.0.0.1"

    url = os.getenv("FASTMDX_DASHBOARD_URL") or f"http://{host}:{port}"
    return url, enabled


def _cmd_dashboard_home() -> int:
    """Start the dashboard home screen for an empty CLI invocation."""
    from fastmdxplora.live.server import serve_dashboard

    serve_dashboard(
        output=Path.cwd(),
        host="127.0.0.1",
        port=8765,
    )
    return 0


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

    raw_argv = list(sys.argv[1:] if argv is None else argv)

    # Show the FastMDXplora identity as soon as the CLI starts. Keep version
    # and citation output machine-friendly; help and an empty invocation are
    # intentionally branded.
    if not any(flag in raw_argv for flag in ("--version", "-V", "--cite")):
        from fastmdxplora.utils.presenter import get_presenter

        dashboard_url, dashboard_enabled = _startup_dashboard_details(raw_argv)
        get_presenter().welcome(
            dashboard_url=dashboard_url,
            dashboard_enabled=dashboard_enabled,
        )

    parser = _build_parser()
    args = parser.parse_args(raw_argv)

    # Short-circuit cheap flags first so a missing chemistry backend never
    # *blocks* `--cite`, `--version`, or `--help`. These flags are how
    # users diagnose the install gap, so guarding them would defeat their
    # purpose.
    if args.cite:
        print(__citation__)
        return 0
    if args.command is None:
        return _cmd_dashboard_home()

    # Setup and simulation phases already handle missing optional chemistry
    # dependencies gracefully by recording the skipped work in their manifests.
    # Do not abort the CLI here: doing so prevents setup-only/config workflows
    # and the test matrix from exercising that documented fallback behavior.

    if args.command == "init-config":
        return _cmd_init_config(args)
    if args.command == "dashboard":
        return _cmd_dashboard(args)

    if args.command == "health":
        return _cmd_health(args)

    if args.command == "install":
        return _cmd_bootstrap(args, editable=False, package_name=".")
    if args.command == "bootstrap":
        return _cmd_bootstrap(args, editable=False, package_name=".")
    if args.command == "install-e":
        return _cmd_bootstrap(args, editable=True, package_name=".")

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
