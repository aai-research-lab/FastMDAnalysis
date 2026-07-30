"""OpenMM simulation runner.

Takes the serialized ``System`` + ``State`` produced by the setup phase
and runs a standard four-stage MD pipeline:

  1. Minimization  -- local energy minimizer to a force-tolerance
  2. NVT equilibration  -- fixed volume, Langevin thermostat
  3. NPT equilibration  -- Monte Carlo barostat added, box equilibrates
  4. Production  -- the trajectory frames the analysis phase consumes

The runner is intentionally separate from the orchestrator-facing
:mod:`fastmdxplora.simulation.pipeline` so it can be exercised
directly from Python for tests and ad-hoc scripts.

OpenMM is a conda-forge package and is imported lazily — without it the
runner raises a helpful ImportError on first use, but importing this
module does not.
"""

from __future__ import annotations

import csv
import math
import shutil
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from fastmdxplora.utils.logging import get_logger

logger = get_logger("simulation.runner")


# ---------------------------------------------------------------------------
# Defaults — production-cadence reporters and standard
# stage step counts.
# ---------------------------------------------------------------------------
DEFAULT_TIMESTEP_FS = 2.0
DEFAULT_TEMPERATURE_K = 300.0
DEFAULT_FRICTION_PER_PS = 1.0          # Langevin collision frequency
DEFAULT_PRESSURE_BAR = 1.0
DEFAULT_BAROSTAT_FREQUENCY = 25        # MC barostat step interval
DEFAULT_INTEGRATOR = "langevin_middle"
DEFAULT_INTEGRATOR_ERROR_TOLERANCE = 0.001  # for the variable-step integrators

# Atmospheres to bar (OpenMM's barostat takes bar). 1 atm = 1.01325 bar.
ATM_TO_BAR = 1.01325

# Integrators we can construct. langevin_middle is the modern default
# (better configurational sampling than the legacy LangevinIntegrator).
SUPPORTED_INTEGRATORS = (
    "langevin",
    "langevin_middle",
    "brownian",
    "verlet",
    "variable_langevin",
    "variable_verlet",
)

# Standard stage step counts for general-purpose MD
DEFAULT_NVT_STEPS = 250_000            # 500 ps @ 2 fs
DEFAULT_NPT_STEPS = 500_000            # 1 ns  @ 2 fs
DEFAULT_PRODUCTION_STEPS = 1_000_000   # 2 ns  @ 2 fs
DEFAULT_MINIMIZE_TOLERANCE_KJMOL_PER_NM = 10.0
DEFAULT_MINIMIZE_MAX_ITERATIONS = 0     # 0 == until convergence

# Reporter cadence
DEFAULT_TRAJECTORY_INTERVAL_STEPS = 1000
DEFAULT_STATE_INTERVAL_STEPS = 1000
DEFAULT_CHECKPOINT_INTERVAL_STEPS = 10_000


# ---------------------------------------------------------------------------
# Result struct
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SimulationResult:
    """What the runner returns. All paths are absolute."""

    trajectory: Path
    topology: Path
    final_state: Path
    energy_csv: Path
    log_file: Path
    platform_used: str
    n_production_frames: int
    duration_ns_actual: float
    minimized_state: Path | None = None


# ---------------------------------------------------------------------------
# Lazy OpenMM import
# ---------------------------------------------------------------------------
def _import_openmm() -> dict:
    """Return a dict of the OpenMM symbols the runner needs.

    Raises a clean ImportError when OpenMM isn't installed. Bundling the
    imports into one function means the runner module itself imports
    cheaply, and tests can substitute the dict via monkeypatching.
    """
    try:
        import openmm
        from openmm import unit
        from openmm.app import (
            CheckpointReporter,
            DCDReporter,
            PDBFile,
            Simulation,
            StateDataReporter,
        )
    except ImportError as exc:
        raise ImportError(
            "Simulation phase requires OpenMM. Install via conda "
            "(recommended): conda install -c conda-forge openmm — or via "
            "pip with the optional [md] extras: "
            "pip install fastmdxplora[md]."
        ) from exc

    return {
        "openmm": openmm,
        "unit": unit,
        "CheckpointReporter": CheckpointReporter,
        "DCDReporter": DCDReporter,
        "PDBFile": PDBFile,
        "Simulation": Simulation,
        "StateDataReporter": StateDataReporter,
    }


# ---------------------------------------------------------------------------
# Platform selection
# ---------------------------------------------------------------------------
def select_platform(
    omm: dict,
    requested: str = "auto",
    precision: str = "mixed",
    device_index: str | int | None = None,
) -> tuple[Any, dict[str, str], str]:
    """Pick the best available OpenMM Platform.

    Parameters
    ----------
    omm : dict
        Output of :func:`_import_openmm`.
    requested : str, default "auto"
        One of ``"auto"``, ``"CUDA"``, ``"OpenCL"``, ``"CPU"``, ``"HIP"``.
        ``"auto"`` tries CUDA → OpenCL → CPU and uses the first that loads.
    precision : str, default "mixed"
        Numerical precision for GPU platforms. ``"single"``, ``"mixed"``,
        or ``"double"``. Ignored for CPU.
    device_index : str | int | None
        GPU device index for multi-GPU machines (e.g. ``"0"`` or ``"0,1"``).
        Maps to ``CudaDeviceIndex`` / ``OpenCLDeviceIndex``. Ignored for CPU.

    Returns
    -------
    platform : openmm.Platform
    properties : dict[str, str]
        Per-platform properties to pass to ``Simulation(...)``.
    name : str
        The platform's name as a string for logging.
    """
    Platform = omm["openmm"].Platform

    auto = requested == "auto"
    if auto:
        candidates = ["CUDA", "OpenCL", "CPU"]
    else:
        candidates = [requested]

    def _make_props(name: str) -> dict[str, str]:
        props: dict[str, str] = {}
        if name in ("CUDA", "OpenCL", "HIP"):
            props["Precision"] = precision
            if device_index is not None:
                prop_key = {
                    "CUDA": "CudaDeviceIndex",
                    "OpenCL": "OpenCLDeviceIndex",
                    "HIP": "HipDeviceIndex",
                }[name]
                props[prop_key] = str(device_index)
        return props

    def _platform_usable(platform: Any, name: str, props: dict[str, str]) -> bool:
        """A registered platform may still be unusable (e.g. OpenCL with no
        device). Probe it by building a trivial Context; only meaningful for
        the GPU platforms. CPU/Reference are always usable."""
        if name not in ("CUDA", "OpenCL", "HIP"):
            return True
        try:
            mm = omm["openmm"]
            unit = omm["unit"]
            sys_ = mm.System()
            sys_.addParticle(1.0)  # one dummy particle
            integ = mm.VerletIntegrator(0.001 * unit.picoseconds)
            ctx = mm.Context(sys_, integ, platform, props)
            del ctx, integ, sys_
            return True
        except Exception as exc:  # noqa: BLE001
            logger.info(
                "Platform %s is registered but not usable (%s); "
                "trying next candidate.", name, exc,
            )
            return False

    for name in candidates:
        try:
            platform = Platform.getPlatformByName(name)
        except Exception:  # noqa: BLE001 -- OpenMM raises different types
            continue
        props = _make_props(name)
        # For auto-selection, verify the platform actually works before
        # committing to it — a registered OpenCL/CUDA platform with no
        # usable device otherwise fails later at Context construction with
        # a confusing error. For an explicit request we honor it as-is so
        # the user sees the real error if their chosen platform is broken.
        if auto and not _platform_usable(platform, name, props):
            continue
        logger.info(
            "Selected OpenMM platform: %s (precision=%s%s)",
            name, precision,
            f", device={device_index}" if device_index is not None else "",
        )
        return platform, props, name

    raise RuntimeError(
        f"No usable OpenMM platform among {candidates}. "
        f"Install GPU drivers + a matching OpenMM build, or pass platform='CPU'."
    )


# ---------------------------------------------------------------------------
# Integrator + barostat helpers
# ---------------------------------------------------------------------------
def _make_integrator(
    omm: dict,
    *,
    name: str,
    temperature_K: float,
    friction_per_ps: float,
    timestep_fs: float,
    error_tolerance: float,
    random_seed: int | None,
):
    """Construct an OpenMM integrator by name.

    Supported names (see :data:`SUPPORTED_INTEGRATORS`):

      - ``langevin_middle`` -- LangevinMiddleIntegrator (default; best
        configurational sampling, the modern recommendation)
      - ``langevin``        -- LangevinIntegrator (legacy)
      - ``brownian``        -- BrownianIntegrator (overdamped)
      - ``verlet``          -- VerletIntegrator (NVE; no thermostat)
      - ``variable_langevin`` -- VariableLangevinIntegrator (adaptive dt)
      - ``variable_verlet``   -- VariableVerletIntegrator (adaptive dt)

    The fixed-timestep integrators take ``timestep_fs``; the variable
    ones take ``error_tolerance`` instead. The thermostatted integrators
    take ``temperature_K`` and ``friction_per_ps``; Verlet variants don't
    (a thermostat must come from a separate force / barostat coupling).
    """
    unit = omm["unit"]
    openmm = omm["openmm"]
    key = str(name).lower()

    T = temperature_K * unit.kelvin
    gamma = friction_per_ps / unit.picoseconds
    dt = timestep_fs * unit.femtoseconds

    if key == "langevin_middle":
        integ = openmm.LangevinMiddleIntegrator(T, gamma, dt)
    elif key == "langevin":
        integ = openmm.LangevinIntegrator(T, gamma, dt)
    elif key == "brownian":
        integ = openmm.BrownianIntegrator(T, gamma, dt)
    elif key == "verlet":
        integ = openmm.VerletIntegrator(dt)
    elif key == "variable_langevin":
        integ = openmm.VariableLangevinIntegrator(T, gamma, float(error_tolerance))
    elif key == "variable_verlet":
        integ = openmm.VariableVerletIntegrator(float(error_tolerance))
    else:
        raise ValueError(
            f"Unknown integrator {name!r}. Supported: "
            f"{', '.join(SUPPORTED_INTEGRATORS)}."
        )

    if random_seed is not None and hasattr(integ, "setRandomNumberSeed"):
        integ.setRandomNumberSeed(int(random_seed))
    return integ


def _add_barostat(omm: dict, system: Any, *, temperature_K: float,
                  pressure_bar: float, frequency: int) -> int:
    """Add a Monte Carlo barostat. Returns the force index for later removal."""
    unit = omm["unit"]
    barostat = omm["openmm"].MonteCarloBarostat(
        pressure_bar * unit.bar,
        temperature_K * unit.kelvin,
        int(frequency),
    )
    return system.addForce(barostat)


def _remove_force(system: Any, force_index: int) -> None:
    """Remove a force previously added by index."""
    system.removeForce(force_index)


# ---------------------------------------------------------------------------
# Reporter attach/detach helpers
# ---------------------------------------------------------------------------
def _attach_state_reporter(
    omm: dict,
    simulation: Any,
    csv_path: Path,
    *,
    interval: int,
    total_steps: int,
) -> Any:
    """Attach a CSV StateDataReporter with the standard observables."""
    # OpenMM's StateDataReporter writes a one-line header automatically.
    # We open with newline="" so the line endings are consistent
    # cross-platform and so the CSV opens cleanly in Excel.
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    reporter = omm["StateDataReporter"](
        str(csv_path),
        interval,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=total_steps,
        separator=",",
    )
    simulation.reporters.append(reporter)
    return reporter


def _attach_dcd_reporter(
    omm: dict, simulation: Any, dcd_path: Path, *, interval: int
) -> Any:
    dcd_path.parent.mkdir(parents=True, exist_ok=True)
    reporter = omm["DCDReporter"](str(dcd_path), interval)
    simulation.reporters.append(reporter)
    return reporter


def _attach_checkpoint_reporter(
    omm: dict, simulation: Any, chk_path: Path, *, interval: int
) -> Any:
    """Attach a binary CheckpointReporter for crash recovery / restart.

    OpenMM writes a portable ``.chk`` file every ``interval`` steps; the
    latest checkpoint can be loaded to resume a run from where it left
    off. Skipped cleanly when ``interval <= 0``.
    """
    if interval <= 0:
        return None
    chk_path.parent.mkdir(parents=True, exist_ok=True)
    reporter = omm["CheckpointReporter"](str(chk_path), interval)
    simulation.reporters.append(reporter)
    return reporter


def _detach_all_reporters(simulation: Any) -> None:
    # Closing the trajectory reporter is what flushes the DCD header.
    for r in simulation.reporters:
        close = getattr(r, "close", None)
        if callable(close):
            try:
                close()
            except Exception:  # noqa: BLE001
                pass
    simulation.reporters.clear()


# ---------------------------------------------------------------------------
# Stage runners
# ---------------------------------------------------------------------------
def _run_minimize(
    omm: dict,
    simulation: Any,
    *,
    tolerance_kjmol_per_nm: float,
    max_iterations: int,
    on_progress: Callable[[str], None] | None = None,
) -> None:
    unit = omm["unit"]
    if on_progress:
        on_progress("Minimizing energy...")
    simulation.minimizeEnergy(
        tolerance=tolerance_kjmol_per_nm * unit.kilojoules_per_mole / unit.nanometer,
        maxIterations=int(max_iterations),
    )


def _flatten_numbers(value: Any):
    """Yield numeric leaves from nested OpenMM/numpy/list containers."""
    if isinstance(value, (int, float)):
        yield float(value)
        return
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (str, bytes)):
        return
    try:
        iterator = iter(value)
    except TypeError:
        return
    for item in iterator:
        yield from _flatten_numbers(item)


def _value_in_unit(quantity: Any, unit_value: Any) -> Any:
    if hasattr(quantity, "value_in_unit"):
        return quantity.value_in_unit(unit_value)
    return quantity


def _validation_error(stage: str, detail: str) -> RuntimeError:
    return RuntimeError(
        f"Invalid simulation state after {stage}: {detail}. "
        "Try safer settings: lower --simulate-timestep-fs, lower "
        "--simulate-temperature-K, increase --simulate-friction-per-ps, use "
        "--simulate-precision double, or disable NPT for the first smoke test "
        "with --simulate-npt-steps 0."
    )


def _validate_state_finite(omm: dict, simulation: Any, *, stage: str) -> None:
    """Validate finite positions and potential energy at a stage boundary."""
    unit = omm["unit"]
    try:
        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
            enforcePeriodicBox=True,
        )
    except TypeError:
        try:
            state = simulation.context.getState()
        except Exception as exc:  # noqa: BLE001
            raise _validation_error(stage, f"could not evaluate state ({exc})") from exc
    except Exception as exc:  # noqa: BLE001
        raise _validation_error(stage, f"could not evaluate state ({exc})") from exc

    try:
        positions = state.getPositions(asNumpy=True)
    except TypeError:
        positions = state.getPositions()
    except AttributeError as exc:
        raise _validation_error(stage, "positions unavailable") from exc

    try:
        position_values = _value_in_unit(positions, unit.nanometer)
    except Exception as exc:  # noqa: BLE001
        raise _validation_error(stage, f"could not read positions ({exc})") from exc
    position_numbers = list(_flatten_numbers(position_values))
    if not position_numbers:
        raise _validation_error(stage, "no positions found")
    if not all(math.isfinite(x) for x in position_numbers):
        raise _validation_error(stage, "positions contain NaN or Inf")

    try:
        energy = state.getPotentialEnergy()
    except AttributeError as exc:
        raise _validation_error(stage, "potential energy unavailable") from exc
    try:
        energy_value = _value_in_unit(energy, unit.kilojoules_per_mole)
    except Exception as exc:  # noqa: BLE001
        raise _validation_error(stage, f"could not read potential energy ({exc})") from exc
    energy_numbers = list(_flatten_numbers(energy_value))
    if not energy_numbers:
        raise _validation_error(stage, "potential energy missing")
    if not all(math.isfinite(x) for x in energy_numbers):
        raise _validation_error(stage, "potential energy is NaN or Inf")


def _run_md_stage(
    simulation: Any,
    *,
    n_steps: int,
    label: str,
    on_progress: Callable[[str], None] | None = None,
) -> None:
    """Run ``n_steps`` of MD. Skips cleanly if ``n_steps <= 0``."""
    if n_steps <= 0:
        return
    if on_progress:
        on_progress(f"{label}: {n_steps:,} steps")
    try:
        simulation.step(int(n_steps))
    except Exception as exc:  # noqa: BLE001
        raise _validation_error(label, f"OpenMM integration failed ({exc})") from exc


def _run_md_stage_with_live_metrics(
    omm: dict,
    simulation: Any,
    telemetry: Any,
    *,
    n_steps: int,
    label: str,
    on_progress: Callable[[str], None] | None,
    current_step: int,
    total_steps: int,
    timestep_fs: float,
    telemetry_interval: int,
) -> int:
    """Run an MD stage in chunks so live telemetry can sample real state."""
    if n_steps <= 0:
        return current_step
    if on_progress:
        on_progress(f"{label}: {n_steps:,} steps")
    remaining = int(n_steps)
    interval = max(1, int(telemetry_interval or DEFAULT_STATE_INTERVAL_STEPS))
    while remaining > 0:
        chunk = min(interval, remaining)
        try:
            simulation.step(chunk)
        except Exception as exc:  # noqa: BLE001
            raise _validation_error(label, f"OpenMM integration failed ({exc})") from exc
        current_step += chunk
        remaining -= chunk
        _append_live_metric(
            omm,
            simulation,
            telemetry,
            stage=label,
            step=current_step,
            total_steps=total_steps,
            timestep_fs=timestep_fs,
        )
    return current_step


# ---------------------------------------------------------------------------
# Stage planning
# ---------------------------------------------------------------------------
def plan_stages(
    *,
    duration_ns: float | None,
    timestep_fs: float,
    nvt_steps: int | None,
    npt_steps: int | None,
    production_steps: int | None,
    nvt_duration_ns: float | None = None,
    npt_duration_ns: float | None = None,
) -> dict[str, int]:
    """Resolve per-stage step counts from a user's duration spec.

    Equilibration and production are independent. ``duration_ns`` sets
    *production* time only — standard MD convention, what people mean
    when they say "I ran a 10 ns simulation." Equilibration uses the
    the standard default lengths (500 ps NVT + 1 ns NPT) regardless
    of production length, because reaching a stable ensemble takes the
    same wall-time whether the production run is 10 ns or 1000 ns.

    Three ways to override the defaults:

      - ``nvt_steps`` / ``npt_steps`` / ``production_steps``: explicit
        step counts.
      - ``nvt_duration_ns`` / ``npt_duration_ns``: time-flavored
        equivalents for the equilibration stages.
      - ``duration_ns``: production time.

    Step-count overrides win over duration-ns overrides if both are
    supplied (lower-level wins; explicit beats inferred).
    """
    if production_steps is not None:
        auto_prod = int(production_steps)
    elif duration_ns is not None and duration_ns > 0:
        steps_per_ns = int(round(1_000_000.0 / float(timestep_fs)))
        auto_prod = int(round(duration_ns * steps_per_ns))
    else:
        auto_prod = DEFAULT_PRODUCTION_STEPS

    # NVT: fixed default, optionally overridden by ns-flavored kwarg
    if nvt_steps is not None:
        auto_nvt = int(nvt_steps)
    elif nvt_duration_ns is not None and nvt_duration_ns > 0:
        steps_per_ns = int(round(1_000_000.0 / float(timestep_fs)))
        auto_nvt = int(round(nvt_duration_ns * steps_per_ns))
    else:
        auto_nvt = DEFAULT_NVT_STEPS

    # NPT: same pattern
    if npt_steps is not None:
        auto_npt = int(npt_steps)
    elif npt_duration_ns is not None and npt_duration_ns > 0:
        steps_per_ns = int(round(1_000_000.0 / float(timestep_fs)))
        auto_npt = int(round(npt_duration_ns * steps_per_ns))
    else:
        auto_npt = DEFAULT_NPT_STEPS

    return {
        "nvt_steps": auto_nvt,
        "npt_steps": auto_npt,
        "production_steps": auto_prod,
    }


def trajectory_interval_for(
    production_steps: int,
    target_frames: int = 2000,
    min_interval: int = 100,
) -> int:
    """Compute a sensible DCD reporter interval.

    Aims for ~``target_frames`` frames in the production run, with a
    floor at ``min_interval`` to avoid absurd write rates on short runs.
    """
    if production_steps <= 0:
        return DEFAULT_TRAJECTORY_INTERVAL_STEPS
    interval = max(production_steps // target_frames, min_interval)
    return int(interval)


# ---------------------------------------------------------------------------
# Top-level runner
# ---------------------------------------------------------------------------
def run_simulation(
    *,
    system_xml: str | Path,
    state_xml: str | Path,
    topology_pdb: str | Path,
    output_dir: str | Path,
    # Stage controls
    minimize: bool = True,
    minimize_tolerance_kjmol_per_nm: float = DEFAULT_MINIMIZE_TOLERANCE_KJMOL_PER_NM,
    minimize_max_iterations: int = DEFAULT_MINIMIZE_MAX_ITERATIONS,
    nvt_steps: int | None = None,
    npt_steps: int | None = None,
    production_steps: int | None = None,
    duration_ns: float | None = None,
    nvt_duration_ns: float | None = None,
    npt_duration_ns: float | None = None,
    # Integrator
    integrator: str = DEFAULT_INTEGRATOR,
    integrator_error_tolerance: float = DEFAULT_INTEGRATOR_ERROR_TOLERANCE,
    timestep_fs: float = DEFAULT_TIMESTEP_FS,
    temperature_K: float = DEFAULT_TEMPERATURE_K,
    friction_per_ps: float = DEFAULT_FRICTION_PER_PS,
    pressure_bar: float | None = None,
    pressure_atm: float | None = None,
    barostat_frequency: int = DEFAULT_BAROSTAT_FREQUENCY,
    random_seed: int | None = None,
    # Hardware
    platform: str = "auto",
    precision: str = "mixed",
    device_index: str | int | None = None,
    # Reporters
    trajectory_interval_steps: int | None = None,
    state_interval_steps: int = DEFAULT_STATE_INTERVAL_STEPS,
    checkpoint_interval_steps: int = DEFAULT_CHECKPOINT_INTERVAL_STEPS,
    live_telemetry: bool = False,
    telemetry_interval: int = DEFAULT_STATE_INTERVAL_STEPS,
    # Hooks
    on_progress: Callable[[str], None] | None = None,
    # Enhanced sampling
    plumed: dict[str, Any] | None = None,
) -> SimulationResult:
    """Run minimize → NVT → NPT → production and return paths to outputs.

    This is the function the orchestrator-facing
    :mod:`fastmdxplora.simulation.pipeline` calls. It can also be used
    directly from Python:

    >>> from fastmdxplora.simulation.runner import run_simulation
    >>> result = run_simulation(                          # doctest: +SKIP
    ...     system_xml="setup/system.xml",
    ...     state_xml="setup/state.xml",
    ...     topology_pdb="setup/topology.pdb",
    ...     output_dir="simulation/",
    ...     duration_ns=10.0,
    ... )
    """
    omm = _import_openmm()
    unit = omm["unit"]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- Resolve pressure (OpenMM's barostat is in bar) ----------------
    # Accept either pressure_bar or pressure_atm; atm is the unit lab
    # scientists think in (and what AMBER uses), bar is
    # OpenMM-native. If both are given, bar wins (it's the native unit);
    # if neither, default to 1 bar.
    if pressure_bar is not None:
        resolved_pressure_bar = float(pressure_bar)
    elif pressure_atm is not None:
        resolved_pressure_bar = float(pressure_atm) * ATM_TO_BAR
    else:
        resolved_pressure_bar = DEFAULT_PRESSURE_BAR

    # ---- Resolve stage step counts -------------------------------------
    plan = plan_stages(
        duration_ns=duration_ns,
        timestep_fs=timestep_fs,
        nvt_steps=nvt_steps,
        npt_steps=npt_steps,
        production_steps=production_steps,
        nvt_duration_ns=nvt_duration_ns,
        npt_duration_ns=npt_duration_ns,
    )
    if trajectory_interval_steps is None:
        trajectory_interval_steps = trajectory_interval_for(plan["production_steps"])
    total_planned_steps = plan["nvt_steps"] + plan["npt_steps"] + plan["production_steps"]

    # ---- Deserialize System + State + Topology -------------------------
    if on_progress:
        on_progress("Loading System, State, and topology")

    system_xml_path = Path(system_xml)
    state_xml_path = Path(state_xml)
    topology_path = Path(topology_pdb)

    for label, path in [("system_xml", system_xml_path), ("state_xml", state_xml_path),
                        ("topology_pdb", topology_path)]:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    with system_xml_path.open(encoding="utf-8") as fh:
        system = omm["openmm"].XmlSerializer.deserialize(fh.read())
    with state_xml_path.open(encoding="utf-8") as fh:
        state = omm["openmm"].XmlSerializer.deserialize(fh.read())

    pdb = omm["PDBFile"](str(topology_path))
    topology = pdb.topology

    # PLUMED biasing (if enabled) is added just before the production stage,
    # not here — equilibration runs unbiased, matching standard enhanced-
    # sampling protocol. See Stage 4 below.

    # ---- Platform ------------------------------------------------------
    platform_obj, platform_props, platform_name = select_platform(
        omm, requested=platform, precision=precision, device_index=device_index
    )

    # ---- Integrator + Simulation ---------------------------------------
    integrator = _make_integrator(
        omm,
        name=integrator,
        temperature_K=temperature_K,
        friction_per_ps=friction_per_ps,
        timestep_fs=timestep_fs,
        error_tolerance=integrator_error_tolerance,
        random_seed=random_seed,
    )
    simulation = omm["Simulation"](
        topology, system, integrator, platform_obj, platform_props
    )
    simulation.context.setState(state)
    _validate_state_finite(omm, simulation, stage="loading state.xml")

    # Output paths
    traj_path = output_dir / "production.dcd"
    minimized_state_path = output_dir / "state_minimized.xml"
    final_state_path = output_dir / "state_final.xml"
    energy_csv = output_dir / "energy.csv"
    log_path = output_dir / "simulation.log"
    # Copy topology so the analysis phase finds it at the expected path.
    topo_out = output_dir / "topology.pdb"
    if not topo_out.exists() or topo_out.resolve() != topology_path.resolve():
        shutil.copy2(topology_path, topo_out)

    # Open the log file once and route the per-stage progress notes
    # through it in addition to whatever on_progress does.
    log_fh = log_path.open("w", encoding="utf-8")
    telemetry = None
    if live_telemetry:
        from fastmdxplora.live.telemetry import TelemetryWriter

        planned_frames = (
            plan["production_steps"] // trajectory_interval_steps
            if trajectory_interval_steps and plan["production_steps"] > 0 else 0
        )
        telemetry = TelemetryWriter(
            output_dir,
            enabled=True,
            total_steps=total_planned_steps,
            planned_frames=planned_frames,
            timestep_fs=timestep_fs,
            platform=platform_name,
            target_temperature_K=temperature_K,
        )
        telemetry.event("Live simulation telemetry started")
        telemetry.write_status(
            stage="loading",
            status="running",
            current_step=0,
            total_planned_steps=total_planned_steps,
            current_checkpoint_path=(output_dir / "checkpoint.chk").as_posix(),
        )

    def _log_step(msg: str) -> None:
        log_fh.write(msg + "\n")
        log_fh.flush()
        if telemetry is not None:
            telemetry.event(msg)
        if on_progress:
            on_progress(msg)

    try:
        # ---- Stage 1: Minimize ----------------------------------------
        if minimize:
            if telemetry is not None:
                telemetry.mark_stage("minimization", "current", status="running", current_step=0)
                telemetry.event("Minimization started")
            _run_minimize(
                omm,
                simulation,
                tolerance_kjmol_per_nm=minimize_tolerance_kjmol_per_nm,
                max_iterations=minimize_max_iterations,
                on_progress=_log_step,
            )
            _validate_state_finite(omm, simulation, stage="minimization")
            if random_seed is None:
                simulation.context.setVelocitiesToTemperature(
                    temperature_K * unit.kelvin
                )
            else:
                simulation.context.setVelocitiesToTemperature(
                    temperature_K * unit.kelvin,
                    int(random_seed),
                )
            _log_step(f"Reset velocities to {temperature_K:.1f} K after minimization")
            minimized_state = simulation.context.getState(
                getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=True,
            )
            with minimized_state_path.open("w", encoding="utf-8") as fh:
                fh.write(omm["openmm"].XmlSerializer.serialize(minimized_state))
            if telemetry is not None:
                _append_live_metric(
                    omm,
                    simulation,
                    telemetry,
                    stage="minimization",
                    step=0,
                    total_steps=total_planned_steps,
                    timestep_fs=timestep_fs,
                )
                telemetry.mark_stage("minimization", "completed", status="running", current_step=0)
                telemetry.event("Minimization completed")
        elif telemetry is not None:
            telemetry.mark_stage("minimization", "skipped", status="running", current_step=0)
            telemetry.event("Minimization skipped")

        # ---- Stage 2: NVT equilibration -------------------------------
        # Use the integrator's existing thermostat (Langevin). No barostat.
        _attach_state_reporter(
            omm, simulation, energy_csv,
            interval=state_interval_steps,
            total_steps=plan["nvt_steps"] + plan["npt_steps"] + plan["production_steps"],
        )
        current_step = 0
        if telemetry is not None:
            if plan["nvt_steps"] > 0:
                telemetry.mark_stage("nvt", "current", status="running", current_step=current_step)
                telemetry.event("NVT started")
            else:
                telemetry.mark_stage("nvt", "skipped", status="running", current_step=current_step)
                telemetry.event("NVT skipped (0 steps)")
        if telemetry is not None:
            current_step = _run_md_stage_with_live_metrics(
                omm,
                simulation,
                telemetry,
                n_steps=plan["nvt_steps"],
                label="NVT equilibration",
                on_progress=_log_step,
                current_step=current_step,
                total_steps=total_planned_steps,
                timestep_fs=timestep_fs,
                telemetry_interval=telemetry_interval,
            )
        else:
            _run_md_stage(
                simulation,
                n_steps=plan["nvt_steps"],
                label="NVT equilibration",
                on_progress=_log_step,
            )
            current_step += plan["nvt_steps"]
        _validate_state_finite(omm, simulation, stage="NVT equilibration")
        if telemetry is not None:
            _append_live_metric(
                omm,
                simulation,
                telemetry,
                stage="NVT",
                step=current_step,
                total_steps=total_planned_steps,
                timestep_fs=timestep_fs,
            )
            if plan["nvt_steps"] > 0:
                telemetry.mark_stage("nvt", "completed", status="running", current_step=current_step)
                telemetry.event("NVT completed")

        # ---- Stage 3: NPT equilibration -------------------------------
        # Add the barostat and reinitialize the context so the system picks up
        # the new force. Production then continues in NPT.
        if plan["npt_steps"] > 0:
            _add_barostat(
                omm, system,
                temperature_K=temperature_K,
                pressure_bar=resolved_pressure_bar,
                frequency=barostat_frequency,
            )
            simulation.context.reinitialize(preserveState=True)
            if telemetry is not None:
                telemetry.mark_stage("npt", "current", status="running", current_step=current_step)
                telemetry.event("NPT started")
            if telemetry is not None:
                current_step = _run_md_stage_with_live_metrics(
                    omm,
                    simulation,
                    telemetry,
                    n_steps=plan["npt_steps"],
                    label="NPT equilibration",
                    on_progress=_log_step,
                    current_step=current_step,
                    total_steps=total_planned_steps,
                    timestep_fs=timestep_fs,
                    telemetry_interval=telemetry_interval,
                )
            else:
                _run_md_stage(
                    simulation,
                    n_steps=plan["npt_steps"],
                    label="NPT equilibration",
                    on_progress=_log_step,
                )
                current_step += plan["npt_steps"]
            _validate_state_finite(omm, simulation, stage="NPT equilibration")
            if telemetry is not None:
                _append_live_metric(
                    omm,
                    simulation,
                    telemetry,
                    stage="NPT",
                    step=current_step,
                    total_steps=total_planned_steps,
                    timestep_fs=timestep_fs,
                )
                telemetry.mark_stage("npt", "completed", status="running", current_step=current_step)
                telemetry.event("NPT completed")
        elif telemetry is not None:
            telemetry.mark_stage("npt", "skipped", status="running", current_step=current_step)
            telemetry.event("NPT skipped (0 steps)")

        # ---- Stage 4: Production --------------------------------------
        # Production runs in NPT (the standard default ensemble).
        #
        # Enhanced sampling: add the PLUMED biasing force now (not during
        # equilibration), then reinitialize the context so it takes effect —
        # the standard protocol equilibrates unbiased and biases production.
        if plumed:
            from fastmdxplora.simulation.plumed import add_plumed_force
            plumed_force = add_plumed_force(omm, system, plumed, Path(output_dir))
            if plumed_force is not None:
                simulation.context.reinitialize(preserveState=True)
        _attach_dcd_reporter(
            omm, simulation, traj_path, interval=trajectory_interval_steps
        )
        # Checkpoint reporter for crash recovery / restart.
        _attach_checkpoint_reporter(
            omm, simulation, output_dir / "checkpoint.chk",
            interval=checkpoint_interval_steps,
        )
        if telemetry is not None:
            if plan["production_steps"] > 0:
                telemetry.mark_stage("production", "current", status="running", current_step=current_step)
                telemetry.event("Production started")
            else:
                telemetry.mark_stage("production", "skipped", status="running", current_step=current_step)
                telemetry.event("Production skipped (0 steps)")
        if telemetry is not None:
            current_step = _run_md_stage_with_live_metrics(
                omm,
                simulation,
                telemetry,
                n_steps=plan["production_steps"],
                label="Production",
                on_progress=_log_step,
                current_step=current_step,
                total_steps=total_planned_steps,
                timestep_fs=timestep_fs,
                telemetry_interval=telemetry_interval,
            )
        else:
            _run_md_stage(
                simulation,
                n_steps=plan["production_steps"],
                label="Production",
                on_progress=_log_step,
            )
            current_step += plan["production_steps"]

        # ---- Finalize -------------------------------------------------
        # Capture the final State for restarts.
        final_state = simulation.context.getState(
            getPositions=True,
            getVelocities=True,
            enforcePeriodicBox=True,
        )
        with final_state_path.open("w", encoding="utf-8") as fh:
            fh.write(omm["openmm"].XmlSerializer.serialize(final_state))
        if telemetry is not None:
            n_frames = (
                plan["production_steps"] // trajectory_interval_steps
                if plan["production_steps"] > 0 else 0
            )
            _append_live_metric(
                omm,
                simulation,
                telemetry,
                stage="production",
                step=current_step,
                total_steps=total_planned_steps,
                timestep_fs=timestep_fs,
                frame_count=n_frames,
            )
            if plan["production_steps"] > 0:
                telemetry.mark_stage(
                    "production",
                    "completed",
                    status="running",
                    current_step=current_step,
                    current_frame_count=n_frames,
                    simulation_time_completed_ns=(
                        float(current_step) * float(timestep_fs) / 1_000_000.0
                    ),
                )
                telemetry.event("Production completed")
            telemetry.write_status(
                status="running",
                current_step=current_step,
                current_frame_count=n_frames,
                simulation_time_completed_ns=(
                    float(current_step) * float(timestep_fs) / 1_000_000.0
                ),
            )

        _detach_all_reporters(simulation)
    except Exception as exc:
        if telemetry is not None:
            telemetry.event(f"error: {type(exc).__name__}: {exc}", level="error")
            try:
                from fastmdxplora.live.telemetry import read_status

                current = str(read_status(output_dir.parent).get("stage") or "production").lower()
                if current in {"minimization", "nvt", "npt", "production"}:
                    telemetry.mark_stage(
                        current,
                        "failed",
                        status="failed",
                        latest_error=f"{type(exc).__name__}: {exc}",
                    )
                else:
                    telemetry.write_status(
                        status="failed", latest_error=f"{type(exc).__name__}: {exc}"
                    )
            except Exception:
                telemetry.write_status(status="failed", latest_error=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        log_fh.close()

    # Frame-count estimate (DCD frame counts require reading the file
    # header; that's a hard dependency on mdtraj which we don't want here.
    # The step / interval division is exact assuming no early termination.)
    n_frames = (
        plan["production_steps"] // trajectory_interval_steps
        if plan["production_steps"] > 0 else 0
    )
    duration_ns_actual = (
        plan["production_steps"] * timestep_fs / 1_000_000.0
    )

    return SimulationResult(
        trajectory=traj_path,
        topology=topo_out,
        final_state=final_state_path,
        energy_csv=energy_csv,
        log_file=log_path,
        platform_used=platform_name,
        n_production_frames=int(n_frames),
        duration_ns_actual=float(duration_ns_actual),
        minimized_state=minimized_state_path if minimize else None,
    )


# ---------------------------------------------------------------------------
# Convenience reader -- a tiny CSV reader for the energy log (no pandas)
# ---------------------------------------------------------------------------
def read_energy_csv(path: str | Path) -> list[dict[str, str]]:
    """Read an energy.csv file and return a list of dict rows."""
    p = Path(path)
    with p.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return list(reader)


def _append_live_metric(
    omm: dict,
    simulation: Any,
    telemetry: Any,
    *,
    stage: str,
    step: int,
    total_steps: int,
    timestep_fs: float,
    frame_count: int | None = None,
) -> None:
    """Append one real state sample to live telemetry if OpenMM can provide it."""
    unit = omm["unit"]
    try:
        state = simulation.context.getState(getEnergy=True, enforcePeriodicBox=True)
    except Exception as exc:  # noqa: BLE001
        telemetry.event(f"warning: could not sample live metrics ({exc})", level="warning")
        return
    potential = _quantity_to_float(state.getPotentialEnergy(), unit.kilojoules_per_mole)
    kinetic = _quantity_to_float(state.getKineticEnergy(), unit.kilojoules_per_mole)
    total = None
    if potential is not None and kinetic is not None:
        total = potential + kinetic
    sim_time_ns = float(step) * float(timestep_fs) / 1_000_000.0
    progress = (float(step) / float(total_steps) * 100.0) if total_steps else None
    temperature = _temperature_from_kinetic_energy(simulation, kinetic)
    volume_nm3 = _periodic_volume_nm3(state, unit)
    density_g_ml = _system_density_g_ml(simulation, unit, volume_nm3)
    speed_ns_day = _telemetry_speed_ns_day(telemetry, sim_time_ns)
    telemetry.append_metric(
        stage=stage,
        step=step,
        simulation_time_ns=sim_time_ns,
        potential_energy=potential,
        kinetic_energy=kinetic,
        total_energy=total,
        temperature=temperature,
        volume=volume_nm3,
        density=density_g_ml,
        speed=speed_ns_day,
        current_frame_count=frame_count,
        progress_percent=progress,
    )
    stage_text = str(stage).lower()
    if "minim" in stage_text:
        stage_key = "minimization"
    elif "nvt" in stage_text:
        stage_key = "nvt"
    elif "npt" in stage_text:
        stage_key = "npt"
    elif "production" in stage_text:
        stage_key = "production"
    else:
        stage_key = stage_text
    telemetry.mark_stage(
        stage_key,
        "current",
        status="running",
        current_step=step,
        current_frame_count=frame_count,
        simulation_time_completed_ns=sim_time_ns,
    )
    _maybe_write_live_frame(
        omm,
        simulation,
        telemetry,
        step,
        stage=stage_key,
        simulation_time_ns=sim_time_ns,
    )


def _temperature_from_kinetic_energy(
    simulation: Any,
    kinetic_kj_per_mol: float | None,
) -> float | None:
    """Estimate instantaneous temperature using OpenMM's standard DOF rule."""

    if kinetic_kj_per_mol is None:
        return None
    try:
        system = simulation.system
        dof = 3 * int(system.getNumParticles()) - int(system.getNumConstraints())
        for index in range(int(system.getNumForces())):
            force = system.getForce(index)
            if force.__class__.__name__ == "CMMotionRemover":
                dof -= 3
                break
        if dof <= 0:
            return None
        # R in kJ mol^-1 K^-1.  OpenMM's StateDataReporter uses the same
        # equipartition relation for its Temperature column.
        gas_constant = 0.00831446261815324
        return (2.0 * float(kinetic_kj_per_mol)) / (float(dof) * gas_constant)
    except Exception:  # noqa: BLE001 - telemetry is best-effort
        return None


def _periodic_volume_nm3(state: Any, unit: Any) -> float | None:
    try:
        volume = state.getPeriodicBoxVolume()
    except Exception:  # noqa: BLE001
        return None
    try:
        return _quantity_to_float(volume, unit.nanometer**3)
    except Exception:  # noqa: BLE001
        return None


def _system_density_g_ml(
    simulation: Any,
    unit: Any,
    volume_nm3: float | None,
) -> float | None:
    if volume_nm3 is None or volume_nm3 <= 0:
        return None
    try:
        system = simulation.system
        mass_dalton = 0.0
        for index in range(int(system.getNumParticles())):
            mass_dalton += float(
                _value_in_unit(system.getParticleMass(index), unit.dalton)
            )
        # 1 Da / nm^3 = 1.66053906660e-3 g/mL.
        return mass_dalton * 1.66053906660e-3 / float(volume_nm3)
    except Exception:  # noqa: BLE001
        return None


def _telemetry_speed_ns_day(telemetry: Any, sim_time_ns: float) -> float | None:
    try:
        elapsed = (
            datetime.now(timezone.utc) - telemetry.start_time
        ).total_seconds()
    except Exception:  # noqa: BLE001
        return None
    if elapsed <= 0:
        return None
    return float(sim_time_ns) * 86400.0 / float(elapsed)


def _maybe_write_live_frame(
    omm: dict,
    simulation: Any,
    telemetry: Any,
    step: int,
    *,
    stage: str,
    simulation_time_ns: float,
) -> None:
    """Snapshot the current positions to ``simulation/live_frame.pdb``.

    Best-effort. Called from the telemetry loop so the dashboard's
    molecular viewer can refresh without disturbing the simulation.
    All exceptions are swallowed: live frames are a UI affordance and
    must never terminate the run.
    """
    try:
        from fastmdxplora.live.live_frames import write_openmm_live_frame
        positions_state = simulation.context.getState(
            getPositions=True, enforcePeriodicBox=True
        )
        write_openmm_live_frame(
            telemetry.root,
            pdbfile_writer=omm["PDBFile"].writeFile,
            topology=simulation.topology,
            positions=positions_state.getPositions(),
            frame_index=step,
            stage=stage,
            simulation_time_ns=simulation_time_ns,
            archive=True,
        )
    except Exception as exc:  # noqa: BLE001 - dashboard failures must not crash sim
        logger.debug("live frame snapshot skipped: %s", exc)


def _quantity_to_float(quantity: Any, unit_value: Any) -> float | None:
    try:
        value = _value_in_unit(quantity, unit_value)
    except Exception:  # noqa: BLE001
        return None
    numbers = list(_flatten_numbers(value))
    return numbers[0] if numbers else None
