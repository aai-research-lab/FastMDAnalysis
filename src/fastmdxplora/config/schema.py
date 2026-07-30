"""Configuration schema registry.

This module is the single source of truth for FastMDXplora's
configuration surface. Every option a user can set — top-level
(``system``, ``output``, ``include``/``exclude``, ``verbose``) and
per-phase (``setup``, ``simulation``, ``analysis``, ``report``) — is
declared here once, with its type, default, and a human-readable
description.

Four features read from this single registry, so they never drift apart:

  1. **Validation** (:mod:`fastmdxplora.config.loader`) — unknown keys
     are rejected with did-you-mean suggestions; values are type-checked.
  2. **Template generation** (``fastmdx init-config``) — a fully-commented
     YAML template is generated directly from the field descriptions and
     defaults.
  3. **Resolved-config dump** — after a run, the merged configuration is
     written to ``resolved_config.yml`` for reproducibility.
  4. **Documentation** — the field help strings are the canonical
     descriptions used in the template and (eventually) the docs.

The schema deliberately mirrors the keyword arguments accepted by each
phase's ``run()`` function and by :class:`fastmdxplora.FastMDXplora`,
so a config file and the equivalent flags/kwargs produce identical runs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


# ---------------------------------------------------------------------------
# Field descriptor
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Field:
    """One configurable option.

    Parameters
    ----------
    name : str
        The key as it appears in the YAML file and as the kwarg name.
    type : type | tuple[type, ...]
        Accepted Python type(s) after YAML parsing. Used for validation.
        ``list`` means "a YAML list"; element types are not deeply checked
        (MD option lists are heterogeneous enough that element-level
        checking causes more false positives than it's worth).
    default : Any
        The value used when the option is absent. ``None`` means "the
        phase supplies its own default" — we don't duplicate phase
        defaults here; we record ``None`` so the phase's DEFAULTS table
        remains the single source for the actual value.
    help : str
        One-line human-readable description (used in the template).
    example : Any, optional
        A representative value shown in the generated template when the
        default is ``None`` (so the template is illustrative, not blank).
    """

    name: str
    type: type | tuple[type, ...]
    default: Any
    help: str
    example: Any = None


@dataclass(frozen=True)
class PhaseSchema:
    """The schema for one phase (or the top-level block)."""

    name: str
    description: str
    fields: tuple[Field, ...]

    def field_names(self) -> set[str]:
        return {f.name for f in self.fields}

    def get(self, name: str) -> Field | None:
        for f in self.fields:
            if f.name == name:
                return f
        return None


# ---------------------------------------------------------------------------
# Top-level keys
# ---------------------------------------------------------------------------
TOP_LEVEL = PhaseSchema(
    name="(top-level)",
    description="Project-level settings.",
    fields=(
        Field("output", str, None,
              "Output directory for all artifacts. "
              "Default: ./fastmdxplora_output_<UTC-timestamp>.",
              example="./my_study"),
        Field("verbose", bool, False,
              "Stream debug logging to the terminal in addition to the log file."),
        Field("include", list, None,
              "Subset of phases to run, in order. "
              "Mutually exclusive with `exclude`.",
              example=["setup", "simulation", "analysis", "report"]),
        Field("exclude", list, None,
              "Phases to skip. Mutually exclusive with `include`.",
              example=["report"]),
    ),
)


# ---------------------------------------------------------------------------
# Setup phase
# ---------------------------------------------------------------------------
SETUP = PhaseSchema(
    name="setup",
    description="System preparation: fix structure, solvate, ionize, parameterize.",
    fields=(
        Field("ph", float, 7.0,
              "pH for hydrogen placement (sets protonation states)."),
        Field("keep_heterogens", bool, False,
              "Retain non-standard residues (ligands, cofactors, ions)."),
        Field("keep_water", bool, False,
              "Retain crystallographic waters."),
        Field("fixed_pdb", str, None,
              "Path to an already-fixed PDB to use directly, skipping "
              "PDBFixer. Default: run PDBFixer on the input.",
              example="prepared.pdb"),
        Field("forcefield", str, "charmm36",
              "Named force field: charmm36 (default), amber14, amber-fb15. "
              "Resolves to the right XML files and water model. For an "
              "unlisted combination, use `force_field` instead.",
              example="charmm36"),
        Field("force_field", list, None,
              "Raw OpenMM force-field XML file(s) — power-user escape hatch "
              "that overrides `forcefield`. Default: use `forcefield`.",
              example=["charmm36.xml", "charmm36/water.xml"]),
        Field("water_model", str, None,
              "Water model for Modeller (e.g. tip3p, tip4pew). "
              "Default: inferred from the force field.",
              example="tip3p"),
        Field("ligand", (str, list), None,
              "Ligand/cofactor SDF or MOL2 file(s) for protein-ligand "
              "systems. Requires a ligand-capable force field "
              "(forcefield: amber-openff).",
              example="ligand.sdf"),
        Field("ligand_forcefield", str, None,
              "OpenFF small-molecule force field for the ligand "
              "(e.g. openff-2.2.1, gaff-2.2.20). Default: per the chosen "
              "force field.",
              example="openff-2.2.1"),
        Field("ligand_name", str, "LIG",
              "Residue/molecule name assigned to the ligand."),
        Field("ligand_net_charge", int, None,
              "Ligand formal net charge. Default: inferred from the SDF."),
        Field("check_ligand_clashes", bool, True,
              "Fail setup if the ligand pose severely overlaps the protein "
              "(the provided coordinates must be a feasible bound pose)."),
        Field("ligand_clash_threshold_nm", float, 0.15,
              "Minimum allowed ligand-protein contact distance in nm; pairs "
              "closer than this count as a clash."),
        Field("solvent_padding_nm", float, 1.0,
              "Minimum distance (nm) between solute and the box wall."),
        Field("box_shape", str, "cube",
              "Periodic box geometry: cube, dodecahedron, or octahedron."),
        Field("ion_positive", str, "Na+",
              "Counter-ion cation."),
        Field("ion_negative", str, "Cl-",
              "Counter-ion anion."),
        Field("ion_concentration_M", float, 0.15,
              "Target ionic strength in molar (physiological is 0.15)."),
        Field("neutralize", bool, True,
              "Add ions to neutralize the net solute charge."),
        Field("nonbonded_method", str, "PME",
              "Nonbonded method: NoCutoff, CutoffNonPeriodic, "
              "CutoffPeriodic, PME, or Ewald."),
        Field("nonbonded_cutoff_nm", float, 1.0,
              "Real-space nonbonded cutoff in nm (cutoff/PME/Ewald methods)."),
        Field("ewald_error_tolerance", float, 0.0005,
              "Ewald/PME error tolerance."),
        Field("use_switching_function", bool, True,
              "Apply a switching function near the cutoff (cutoff methods)."),
        Field("switch_distance_nm", (int, float), None,
              "Switching-function turn-on distance in nm. "
              "Default: 0.9 × cutoff.",
              example=0.9),
        Field("dispersion_correction", bool, True,
              "Apply the long-range dispersion (vdW tail) correction."),
        Field("remove_cm_motion", bool, False,
              "Add a center-of-mass motion remover."),
        Field("constraints", str, "HBonds",
              "Bond constraints: None, HBonds, AllBonds, or HAngles."),
        Field("rigid_water", bool, True,
              "Constrain water bond lengths and angles."),
        Field("hydrogen_mass_amu", (int, float), None,
              "Hydrogen-mass-repartitioning mass in amu (enables longer "
              "timesteps). Default: off.",
              example=4.0),
        Field("temperature_K", (int, float), 300.0,
              "Temperature in K for initial velocity assignment."),
    ),
)


# ---------------------------------------------------------------------------
# Simulation phase
# ---------------------------------------------------------------------------
SIMULATION = PhaseSchema(
    name="simulation",
    description="Molecular dynamics: minimize, equilibrate (NVT, NPT), produce.",
    fields=(
        Field("preset", str, None,
              "Optional simulation preset. 'gentle' uses conservative "
              "smoke-test settings: 0.5 fs, 100 K, 5/ps friction, no NPT, "
              "and short NVT/production.",
              example="gentle"),
        Field("duration_ns", (int, float), None,
              "Production length in ns (standard MD convention — "
              "equilibration is independent). Default: 2 ns.",
              example=100.0),
        Field("nvt_duration_ns", (int, float), None,
              "NVT equilibration in ns. Default: fixed 500 ps regardless "
              "of production length.",
              example=1.0),
        Field("npt_duration_ns", (int, float), None,
              "NPT equilibration in ns. Default: fixed 1 ns regardless of "
              "production length.",
              example=2.0),
        Field("nvt_steps", int, None,
              "NVT step count (overrides nvt_duration_ns). Default: 250000.",
              example=250000),
        Field("npt_steps", int, None,
              "NPT step count (overrides npt_duration_ns). Default: 500000.",
              example=500000),
        Field("production_steps", int, None,
              "Production step count (overrides duration_ns). Default: 1000000.",
              example=1000000),
        Field("minimize", bool, True,
              "Run energy minimization before equilibration."),
        Field("integrator", str, "langevin_middle",
              "Integrator: langevin_middle, langevin, brownian, verlet, "
              "variable_langevin, or variable_verlet."),
        Field("integrator_error_tolerance", float, 0.001,
              "Error tolerance for the variable-timestep integrators "
              "(variable_langevin / variable_verlet)."),
        Field("minimize_tolerance_kjmol_per_nm", (int, float), 10.0,
              "Minimization force tolerance in kJ/mol/nm."),
        Field("minimize_max_iterations", int, 0,
              "Max minimization iterations (0 = until convergence)."),
        Field("timestep_fs", (int, float), 2.0,
              "Integrator timestep in fs."),
        Field("temperature_K", (int, float), 300.0,
              "Production temperature in K."),
        Field("pressure_bar", (int, float), 1.0,
              "Pressure for the Monte Carlo barostat in bar (OpenMM-native)."),
        Field("pressure_atm", (int, float), None,
              "Pressure in atm (converted to bar internally). Accepted as "
              "an alternative to pressure_bar; bar wins if both are given.",
              example=1.0),
        Field("friction_per_ps", (int, float), 1.0,
              "Langevin thermostat friction coefficient in 1/ps."),
        Field("barostat_frequency", int, 25,
              "MC barostat volume-move attempt interval in steps."),
        Field("random_seed", int, None,
              "Integrator random seed for reproducibility. Default: unset.",
              example=42),
        Field("platform", str, "auto",
              "OpenMM compute platform: auto, CUDA, OpenCL, CPU, or HIP."),
        Field("precision", str, "mixed",
              "GPU precision: single, mixed, or double."),
        Field("device_index", str, None,
              "GPU device index for multi-GPU machines (e.g. '0' or '0,1').",
              example="0"),
        Field("trajectory_interval_steps", int, None,
              "DCD reporter interval in steps. Default: adaptive "
              "(~2000 frames per production run).",
              example=1000),
        Field("state_interval_steps", int, 1000,
              "Energy/state reporter interval in steps."),
        Field("checkpoint_interval_steps", int, 10000,
              "Binary checkpoint (.chk) interval in steps, for restart / "
              "crash recovery. 0 disables checkpointing."),
        Field("live_telemetry", bool, False,
              "Write live_status.json, live_metrics.csv, and live_events.log "
              "for the local live dashboard. Also writes a live-frame PDB so "
              "the molecular viewer can refresh on the simulation."),
        Field("telemetry_interval", int, 1000,
              "Minimum step interval for live dashboard telemetry updates."),
        Field("dashboard_ligand_resname", str, None,
              "Override the ligand residue name surfaced in the dashboard's "
              "ligand tools pane. Auto-detected when omitted.",
              example="EPE"),
        Field("dashboard_binding_pocket_cutoff_A", float, 5.0,
              "Default binding-pocket cutoff in angstrom for the dashboard "
              "molecular viewer (used by the 'show pocket' tools).",
              example=5.0),
        Field("dashboard_max_playback_frames", int, 200,
              "Maximum frames the trajectory-playback panel will load. "
              "Frames are downsampled evenly from the full DCD.",
              example=200),
        Field("plumed", dict, None,
              "Optional PLUMED enhanced-sampling config: a dict with "
              "`enabled` (bool) and `script` (inline PLUMED text or a path "
              "to a .dat file). Requires the openmm-plumed package."),
    ),
)


# ---------------------------------------------------------------------------
# Analysis phase
# ---------------------------------------------------------------------------
ANALYSIS = PhaseSchema(
    name="analysis",
    description="Trajectory analysis: RMSD, RMSF, Rg, H-bonds, SS, SASA, etc.",
    fields=(
        Field("trajectory", str, None,
              "Trajectory file. Default: simulation/production.dcd.",
              example="simulation/production.dcd"),
        Field("topology", str, None,
              "Topology file. Default: simulation/topology.pdb.",
              example="simulation/topology.pdb"),
        Field("include", list, None,
              "Subset of analyses to run. Default: all ten. "
              "Mutually exclusive with `exclude`.",
              example=["rmsd", "rmsf", "rg", "cluster"]),
        Field("exclude", list, None,
              "Analyses to skip. Mutually exclusive with `include`.",
              example=["dimred"]),
        Field("selection", str, None,
              "Default MDTraj atom selection applied across analyses. "
              "Overrides `scope` when set.",
              example="name CA"),
        Field("scope", str, "solute",
              "Atom scope for analyses that don't set their own selection: "
              "solute (protein+ligand, default), protein, ligand, or all. "
              "Keeps analyses off solvent/ions.",
              example="solute"),
        Field("stride", int, None,
              "Load every Nth frame from the trajectory.",
              example=1),
        Field("first", int, None,
              "First frame index to include.",
              example=0),
        Field("last", int, None,
              "Last frame index (exclusive). Default: full trajectory.",
              example=10000),
        Field("options", dict, None,
              "Per-analysis option overrides, keyed by analysis name. "
              "E.g. {cluster: {methods: [kmeans], n_clusters: 5}}.",
              example={"cluster": {"methods": ["kmeans"], "n_clusters": 5}}),
    ),
)


# ---------------------------------------------------------------------------
# Report phase
# ---------------------------------------------------------------------------
REPORT = PhaseSchema(
    name="report",
    description="Generate the Markdown report, PPTX slides, and project bundle.",
    fields=(
        Field("title", str, None,
              "Report title. Default: auto-generated from the system name.",
              example="My MD Study"),
        Field("author", str, None,
              "Author name recorded in the report metadata.",
              example="A. Aina"),
        Field("document", bool, True,
              "Generate the Markdown study report."),
        Field("slides", bool, True,
              "Generate the PPTX slide deck."),
        Field("bundle", bool, True,
              "Generate the self-contained project_bundle.zip."),
        Field("include_methods", bool, True,
              "Include the Methods section in the document."),
        Field("include_reproducibility", bool, True,
              "Include the Reproducibility section in the document."),
        Field("region_highlights", list, None,
              "Optional user-defined residue ranges to highlight on RMSF "
              "report figures. Each item should include label, start, end, "
              "and optionally color.",
              example=[
                  {"label": "example region 1", "start": 3, "end": 7, "color": "#4E79A7"},
                  {"label": "example helix", "start": 10, "end": 14, "color": "#F28E2B"},
              ]),
        Field("comparison", bool, True,
              "For a multi-run study, build the cross-run comparison report "
              "(overlays + parameter trends) under comparison/."),
    ),
)


# ---------------------------------------------------------------------------
# Execution (batch run scheduling)
# ---------------------------------------------------------------------------
EXECUTION = PhaseSchema(
    name="execution",
    description="How the runs are scheduled: sequentially or in parallel.",
    fields=(
        Field("mode", str, "sequential",
              "Run scheduling: 'sequential' (one at a time) or 'parallel'."),
        Field("workers", int, None,
              "Parallel worker count. Default: one per device if `devices` "
              "is set, else the CPU count (capped at the number of runs).",
              example=2),
        Field("devices", list, None,
              "GPU device indices to distribute parallel runs across, "
              "round-robin (one run per device). GPU runs only.",
              example=[0, 1]),
        Field("continue_on_error", bool, True,
              "If a run fails, record it and continue (True) or stop the "
              "whole batch (False)."),
    ),
)


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------
PHASE_SCHEMAS: dict[str, PhaseSchema] = {
    "setup": SETUP,
    "simulation": SIMULATION,
    "analysis": ANALYSIS,
    "report": REPORT,
}

# Phase keys recognized at the top level of a config file (the per-phase
# option blocks).
PHASE_KEYS = tuple(PHASE_SCHEMAS.keys())

# Batch top-level keys. `systems` is the canonical (and only) way to
# specify input — always a list, even for a single system. `sweep`
# defines parameter axes. They have bespoke structure (a list of
# mappings; a mapping of dotted-axis -> value-list) validated by the
# batch layer rather than the per-field type checker.
BATCH_KEYS = ("systems", "sweep", "execution")

# All keys recognized at the top level: the scalar top-level fields plus
# the per-phase block names plus the batch keys.
TOP_LEVEL_KEYS = TOP_LEVEL.field_names() | set(PHASE_KEYS) | set(BATCH_KEYS)


def all_schemas() -> dict[str, PhaseSchema]:
    """Return every schema, including the top-level pseudo-phase."""
    return {"(top-level)": TOP_LEVEL, **PHASE_SCHEMAS}
