"""Simulation pipeline: minimize, equilibrate, run production MD.

This module's public surface is :func:`run`, called by the FastMDXplora
orchestrator. Starting in v0.2.0 it performs real MD via OpenMM:

  1. Locate ``system.xml``, ``state.xml``, ``topology.pdb`` produced by
     the setup phase. (If they're missing, skip with a clear note.)
  2. Delegate to :func:`fastmdxplora.simulation.runner.run_simulation`
     for the actual minimize → NVT → NPT → production stages.
  3. Write ``simulation_parameters.json`` recording the full execution
     plan and what was produced.

Defaults
--------
Stage step counts and timestep follow standard
``build_auto_config`` so users moving between the two tools see the same
behavior:

  - Minimize: until convergence (10 kJ/mol/nm tolerance)
  - NVT: 250,000 steps (500 ps at 2 fs)
  - NPT: 500,000 steps (1 ns at 2 fs)
  - Production: 1,000,000 steps (2 ns at 2 fs)
  - Timestep: 2 fs; Temperature: 300 K; HBonds constraints (set in setup)

Pass ``duration_ns=`` to set the production length (standard MD
convention — "I ran a 10 ns simulation" means 10 ns of production).
Equilibration is independent: it uses fixed standard defaults
(500 ps NVT + 1 ns NPT) regardless of production length, because
reaching a stable ensemble takes the same wall-time whether the
production is 10 ns or 1000 ns.

To customize equilibration, pass ``nvt_duration_ns=`` /
``npt_duration_ns=`` (or the lower-level ``nvt_steps=`` /
``npt_steps=`` for exact control).

Graceful degradation
--------------------
When OpenMM isn't installed (the ``[setup]`` extras), or when the setup
phase didn't produce a ``system.xml`` (it was scaffolded only), the
simulation phase writes a manifest noting what was missing and returns
cleanly without raising. This keeps the project-level pipeline runnable
end-to-end for users who only want analysis or report.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING, Any

from fastmdxplora.dependencies import MissingBackendError, missing_dependencies
from fastmdxplora.config.schema import SIMULATION
from fastmdxplora.utils.logging import get_logger

if TYPE_CHECKING:
    from fastmdxplora.orchestrator import FastMDXplora

logger = get_logger("simulation")


# Default parameters — see runner.py for the precise step counts.
#: What each option starts as, read from the schema that declares it.
#: This was a second copy of those values, and the two drifted: the pH
#: default was 7.4 in the schema and 7.0 here. There is one declaration
#: now, and a phase reads it rather than restating it.
DEFAULTS: dict[str, Any] = SIMULATION.defaults()

PRESETS: dict[str, dict[str, Any]] = {
    "gentle": {
        "timestep_fs": 0.5,
        "temperature_K": 100.0,
        "friction_per_ps": 5.0,
        "nvt_steps": 1000,
        "npt_steps": 0,
        "production_steps": 1000,
        "duration_ns": 0.001,
        "precision": "double",
    },
}


def _resolve_params(options: dict[str, Any]) -> dict[str, Any]:
    """Merge defaults, an optional preset, and explicit user options."""
    preset = options.get("preset")
    if preset is None:
        return {**DEFAULTS, **options}
    preset_key = str(preset).lower()
    if preset_key not in PRESETS:
        valid = ", ".join(sorted(PRESETS))
        raise ValueError(f"Unknown simulation preset {preset!r}. Valid presets: {valid}.")
    explicit = dict(options)
    explicit["preset"] = preset_key
    return {**DEFAULTS, **PRESETS[preset_key], **explicit}


def _setup_outputs_present(setup_dir: Path) -> tuple[Path | None, Path | None, Path | None]:
    """Return paths to setup outputs that exist on disk, else (None, ...)."""
    system_xml = setup_dir / "system.xml"
    state_xml = setup_dir / "state.xml"
    topology = setup_dir / "topology.pdb"
    if not (system_xml.exists() and state_xml.exists() and topology.exists()):
        return None, None, None
    return system_xml, state_xml, topology


def run(
    *,
    orchestrator: "FastMDXplora",
    output_dir: Path,
    **options: Any,
) -> list[str]:
    """Run the simulation phase.

    Parameters
    ----------
    orchestrator : FastMDXplora
    output_dir : pathlib.Path
        Destination for simulation artifacts.
    **options
        Overrides of the module-level :data:`DEFAULTS`.

    Returns
    -------
    list of str
        Paths (relative to ``output_dir``) of artifacts produced.
    """
    params: dict[str, Any] = _resolve_params(options)
    presenter = getattr(orchestrator, "_presenter", None)
    artifacts: list[str] = []
    notes: list[str] = []

    # ---- Locate setup outputs ------------------------------------------
    setup_dir = orchestrator.output_dir / "setup"
    system_xml, state_xml, topology = _setup_outputs_present(setup_dir)
    if system_xml is None:
        notes.append(
            f"Setup outputs not found in {setup_dir} (system.xml / state.xml / "
            f"topology.pdb). Run the setup phase first, or skip simulation."
        )
        if presenter:
            presenter.step(
                "No setup outputs found — run setup first to produce "
                "system.xml/state.xml/topology.pdb",
                status="warning",
            )
        _write_manifest(output_dir, params, artifacts, notes, platform_used=None)
        artifacts.append("simulation_parameters.json")
        missing = missing_dependencies()
        if missing:
            raise MissingBackendError(missing)
        raise RuntimeError(
            f"Simulation cannot start because setup outputs are missing in {setup_dir}. "
            "Run the setup phase first, or choose an analysis-only workflow."
        )

    # ---- Run the simulation --------------------------------------------
    try:
        from fastmdxplora.simulation.runner import run_simulation

        if presenter:
            duration_label = (
                f"{params['duration_ns']} ns"
                if params["duration_ns"] is not None
                else "default stages"
            )
            presenter.step(
                f"Starting MD ({duration_label}, platform={params['platform']})"
            )

        def _progress(msg: str) -> None:
            if presenter:
                presenter.info(msg)

        result = run_simulation(
            system_xml=system_xml,
            state_xml=state_xml,
            topology_pdb=topology,
            output_dir=output_dir,
            minimize=bool(params["minimize"]),
            minimize_tolerance_kjmol_per_nm=float(
                params["minimize_tolerance_kjmol_per_nm"]
            ),
            minimize_max_iterations=int(params["minimize_max_iterations"]),
            nvt_steps=params["nvt_steps"],
            npt_steps=params["npt_steps"],
            production_steps=params["production_steps"],
            duration_ns=params["duration_ns"],
            nvt_duration_ns=params["nvt_duration_ns"],
            npt_duration_ns=params["npt_duration_ns"],
            integrator=str(params["integrator"]),
            integrator_error_tolerance=float(params["integrator_error_tolerance"]),
            timestep_fs=float(params["timestep_fs"]),
            temperature_K=float(params["temperature_K"]),
            friction_per_ps=float(params["friction_per_ps"]),
            pressure_bar=params["pressure_bar"],
            pressure_atm=params["pressure_atm"],
            barostat_frequency=int(params["barostat_frequency"]),
            random_seed=params["random_seed"],
            platform=str(params["platform"]),
            precision=str(params["precision"]),
            device_index=params["device_index"],
            trajectory_interval_steps=params["trajectory_interval_steps"],
            state_interval_steps=int(params["state_interval_steps"]),
            checkpoint_interval_steps=int(params["checkpoint_interval_steps"]),
            live_telemetry=bool(params["live_telemetry"]),
            telemetry_interval=int(params["telemetry_interval"]),
            on_progress=_progress,
            plumed=params.get("plumed"),
        )

        # Record artifacts relative to output_dir
        for path in (
            result.trajectory,
            result.topology,
            result.final_state,
            result.energy_csv,
            result.log_file,
        ):
            try:
                artifacts.append(path.relative_to(output_dir).as_posix())
            except ValueError:
                artifacts.append(str(path))
        if result.minimized_state is not None:
            try:
                artifacts.append(result.minimized_state.relative_to(output_dir).as_posix())
            except ValueError:
                artifacts.append(str(result.minimized_state))

        if presenter:
            presenter.step(
                f"Production complete: {result.n_production_frames:,} frames, "
                f"{result.duration_ns_actual:.3f} ns on {result.platform_used}"
            )

        _write_manifest(
            output_dir, params, artifacts, notes,
            platform_used=result.platform_used,
            pressure_bar_used=result.pressure_bar_used,
            n_frames=result.n_production_frames,
            duration_ns_actual=result.duration_ns_actual,
        )
    except ImportError as exc:
        notes.append(f"OpenMM unavailable: {exc}")
        if presenter:
            presenter.step(
                "OpenMM not installed — simulation skipped. "
                "Install via: conda install -c conda-forge openmm",
                status="warning",
            )
        _write_manifest(output_dir, params, artifacts, notes, platform_used=None)
        raise MissingBackendError(missing_dependencies()) from exc
    except Exception as exc:  # noqa: BLE001 -- runtime errors from OpenMM
        # Real runtime failure (numerical instability, bad topology, etc.).
        # Record it in the manifest so the project-level manifest still
        # picks up an actionable trace.
        notes.append(f"Simulation failed: {type(exc).__name__}: {exc}")
        if presenter:
            presenter.step(f"Simulation error: {exc}", status="error")
        _write_manifest(output_dir, params, artifacts, notes, platform_used=None)
        # Re-raise so the orchestrator marks the phase as errored.
        raise

    artifacts.append("simulation_parameters.json")

    if presenter:
        presenter.step("Wrote simulation_parameters.json")

    logger.debug("simulation: wrote %d artifact(s) to %s", len(artifacts), output_dir)
    return artifacts


def _write_manifest(
    output_dir: Path,
    params: dict[str, Any],
    artifacts: list[str],
    notes: list[str],
    *,
    platform_used: str | None,
    n_frames: int | None = None,
    duration_ns_actual: float | None = None,
    pressure_bar_used: float | None = None,
) -> None:
    """Write ``simulation_parameters.json`` with full provenance."""
    canonical = {
        "trajectory": "production.dcd",
        "topology": "topology.pdb",
        "state": "state_final.xml",
        "minimized_state": "state_minimized.xml",
        "energy_log": "energy.csv",
        "stdout_log": "simulation.log",
        "checkpoint": "checkpoint.chk",
    }
    manifest = {
        "phase": "simulation",
        "parameters": params,
        "platform_used": platform_used,
        # What the barostat ran at, rather than what was asked for: pressure
        # can be given in bar or atmospheres and unset means one bar, so the
        # number used was known only inside the runner.
        "pressure_bar_used": pressure_bar_used,
        "n_production_frames": n_frames,
        "duration_ns_actual": duration_ns_actual,
        "artifacts_planned": canonical,
        "artifacts_written": list(artifacts),
        "notes": notes,
    }
    with (output_dir / "simulation_parameters.json").open("w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)
