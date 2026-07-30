"""The FastMDXplora project-level orchestrator.

This module implements the central orchestrator class. Following the
phase-based orchestration pattern (Aina & Kwan, JCC 2026),
the orchestrator:

  1. Holds shared project state (system input, output directory, options)
  2. Knows its registered phases (setup, simulate, analyze, report)
  3. Applies intelligent defaults and validates per-phase options
  4. Executes phases in coordinated sequence
  5. Consolidates outputs into a single project directory

Unlike a generic workflow engine (Snakemake, Nextflow, Galaxy), the workflow
is built-in and the user expresses intent through include/exclude and option
overrides, not by describing a DAG (directed acyclic graph: the
task-and-dependency model those engines use).
"""

from __future__ import annotations

import json
import logging
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from fastmdxplora.utils.logging import get_logger

logger = get_logger("project")


# ---------------------------------------------------------------------------
# The canonical phase list. This is the only place phase order is defined.
# ---------------------------------------------------------------------------
PHASES: tuple[str, ...] = ("setup", "simulation", "analysis", "report")


@dataclass
class PhaseResult:
    """Lightweight record of a single phase invocation."""

    name: str
    status: str  # "ok" | "skipped" | "error"
    output_dir: Path | None = None
    started_at: str = ""
    finished_at: str = ""
    message: str = ""
    artifacts: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "status": self.status,
            "output_dir": str(self.output_dir) if self.output_dir else None,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "message": self.message,
            "artifacts": self.artifacts,
        }


@dataclass
class RunResult:
    """Result of one run within an exploration.

    ``explore()`` always returns a list of these — a single study is a
    list of one, a sweep is a list of many. Each carries the run's
    identity and the per-phase results inside ``phases``.
    """

    run_id: str
    system: str
    status: str  # "ok" | "error" | "skipped"
    output_dir: Path | None = None
    sweep_values: dict[str, Any] = field(default_factory=dict)
    phases: list[PhaseResult] = field(default_factory=list)
    message: str = ""
    error_type: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "system": self.system,
            "status": self.status,
            "output_dir": str(self.output_dir) if self.output_dir else None,
            "sweep_values": self.sweep_values,
            "phases": [p.to_dict() for p in self.phases],
            "message": self.message,
            "error_type": self.error_type,
        }

    # Convenience: treat a RunResult a bit like its phase list, so common
    # patterns (iterating phases, checking a phase status) stay ergonomic.
    def phase(self, name: str) -> PhaseResult | None:
        """Return the PhaseResult for ``name``, or None if it didn't run."""
        for p in self.phases:
            if p.name == name:
                return p
        return None


class FastMDXplora:
    """Project-level orchestrator for end-to-end MD studies.

    Parameters
    ----------
    system : str | os.PathLike
        Input for a single study. Accepted forms (auto-detected):

        - Path to a PDB / CIF file (e.g. ``"protein.pdb"``)
        - 4-character PDB ID (e.g. ``"1L2Y"``), fetched from RCSB
        - One-letter amino-acid sequence, if structure prediction is
          available (future)

        Mutually exclusive with ``config``.
    config : str | os.PathLike | None
        Path to a YAML config file. Drives one system or many (with an
        optional parameter sweep and parallel execution); the interface
        is the same either way. Mutually exclusive with ``system``.
    output_dir : str | os.PathLike | None
        Where to write project outputs. Defaults to
        ``./fastmdxplora_output_<timestamp>``.
    options : dict[str, dict] | None
        Per-phase keyword arguments, e.g.
        ``{"simulation": {"duration_ns": 100}}``.
    verbose : bool
        If True, log progress to stdout in addition to the project log file.
    include, exclude : list[str] | None
        Default phase selection (``explore()`` arguments still override).

    Examples
    --------
    >>> fmdx = FastMDXplora(system="protein.pdb")
    >>> fmdx.explore()                          # doctest: +SKIP

    >>> fmdx = FastMDXplora(system="1L2Y")     # PDB ID, fetched from RCSB
    >>> fmdx.explore(                            # doctest: +SKIP
    ...     include=["setup", "simulation"],
    ...     options={"simulation": {"duration_ns": 50}},
    ... )

    >>> # A config file: one system or many, same interface:
    >>> fmdx = FastMDXplora(config="study.yml")
    >>> fmdx.explore()                          # doctest: +SKIP
    """

    def __init__(
        self,
        system: str | os.PathLike | None = None,
        *,
        config: str | os.PathLike | None = None,
        config_data: dict[str, Any] | None = None,
        output_dir: str | os.PathLike | None = None,
        options: dict[str, dict[str, Any]] | None = None,
        verbose: bool = False,
        include: list[str] | None = None,
        exclude: list[str] | None = None,
    ) -> None:
        # FastMDXplora is the single user-facing entry point. Ways to
        # construct it:
        #   - config=...      : a YAML config path (one system or many, with
        #                       optional sweep / parallel execution).
        #   - config_data=... : the same, as an already-parsed dict (used by
        #                       the CLI, which assembles a config from flags).
        #   - system=...      : a single concrete study, run directly. This is
        #                       also the path the internal batch worker uses
        #                       for each run, so it must not recurse.
        # config/config_data execution is deferred to explore().
        n_config = sum(x is not None for x in (config, config_data))
        if n_config and system is not None:
            raise ValueError(
                "Pass either `system=` (a single study) or a config "
                "(`config=` / `config_data=`), not both."
            )

        self._config_path: str | None = (
            str(config) if config is not None else None
        )
        self._config_data: dict[str, Any] | None = config_data
        self._deferred_output_dir = output_dir
        self._deferred_verbose = verbose

        if n_config:
            # Config-driven: defer everything to explore(). We don't create
            # an output directory or banner here because the batch machinery
            # owns the layout (flat for one run, runs/<id>/ for many).
            self.system = None  # resolved per-run by the batch layer
            self.options = options or {}
            self.verbose = bool(verbose)
            self._config_include = include
            self._config_exclude = exclude
            self.results = []
            return

        # ---- Direct single-study path -----------------------------------
        if system is None:
            raise ValueError(
                "FastMDXplora requires either a `system` input (a PDB/CIF "
                "file path, a 4-character PDB ID, or a one-letter sequence) "
                "or a `config` file."
            )

        self.system: str = str(system)

        # Phase selection (the batch layer passes the config's include/exclude)
        self._config_include: list[str] | None = include
        self._config_exclude: list[str] | None = exclude

        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        self.output_dir: Path = (
            Path(output_dir) if output_dir
            else Path(f"fastmdxplora_output_{timestamp}")
        )
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.options: dict[str, dict[str, Any]] = options or {}
        self.verbose: bool = bool(verbose)

        # Per-phase output subdirectories (created lazily by each phase)
        self._phase_dirs: dict[str, Path] = {
            phase: self.output_dir / phase for phase in PHASES
        }

        # Record of phase executions in this session
        self.results: list[PhaseResult] = []

        self._configure_logging()
        self._presenter = self._configure_presenter()

        # Display the opening banner once the session is ready.
        from fastmdxplora import __version__

        self._presenter.banner(
            System=self.system,
            Output=str(self.output_dir),
            Version=__version__,
        )

        # The banner already shows system/output to the user; this log
        # is for the file/audit trail and verbose console only.
        logger.debug(
            "FastMDXplora initialized: system=%s output=%s", self.system, self.output_dir
        )

    # ------------------------------------------------------------------
    # Orchestration entry point
    # ------------------------------------------------------------------
    def explore(
        self,
        *,
        include: list[str] | None = None,
        exclude: list[str] | None = None,
        options: dict[str, dict[str, Any]] | None = None,
        report: bool = True,
        dry_run: bool = False,
    ) -> list[RunResult]:
        """Run the full pipeline, end to end.

        Parameters
        ----------
        include : list of str, optional
            Phases to run (subset of {"setup", "simulation", "analysis", "report"}).
            If omitted, all phases run.
        exclude : list of str, optional
            Phases to skip. Mutually exclusive with ``include``.
        options : dict, optional
            Per-phase option overrides applied on top of the orchestrator's
            ``options`` attribute.
        report : bool, default True
            Convenience flag. If False, skip the report phase even when
            ``include``/``exclude`` would otherwise enable it.
        dry_run : bool, default False
            If True, print the plan — every run, its system, swept values,
            output directory, and the phases that would execute — and
            return without running anything.

        Returns
        -------
        list[RunResult]
            One :class:`RunResult` per run, always. A single study is a
            list of one; a sweep is a list of many. Each ``RunResult``
            carries its per-phase :class:`PhaseResult` list in ``.phases``.

        Notes
        -----
        When ``include`` / ``exclude`` are omitted here but were set in a
        config file passed to the constructor, the config-file values are
        used. Explicit arguments to this method always win.
        """
        # Config-driven runs (one system or many) go through the batch
        # machinery internally. The user always sees the same FastMDXplora
        # interface; the batch layer is an implementation detail.
        if self._config_path is not None or self._config_data is not None:
            return self._explore_config(
                include=include, exclude=exclude, report=report, dry_run=dry_run,
            )

        # Config-file phase selection is the fallback when this call omits it.
        if include is None and exclude is None:
            include = self._config_include
            exclude = self._config_exclude

        plan = self._build_plan(include=include, exclude=exclude, want_report=report)

        # Dry run: report the plan and return without executing.
        if dry_run:
            self._print_dry_run_single(plan)
            return [RunResult(
                run_id="s1", system=self.system, status="planned",
                output_dir=self.output_dir, phases=[],
            )]

        merged_options = self._merge_options(options)
        dashboard_writer = self._dashboard_writer(merged_options)
        if dashboard_writer is not None:
            self._initialize_dashboard_timeline(dashboard_writer, plan)

        # Remember the resolved phase selection + merged options so the
        # resolved_config.yml dump reflects what *actually* ran (including
        # any per-call overrides), not just the construction-time config.
        self._resolved_include = include
        self._resolved_exclude = exclude
        self._resolved_options = {
            p: opts for p, opts in merged_options.items() if opts
        }

        # Plan goes to file/audit; the presenter shows headers visually.
        logger.debug("Plan: %s", " -> ".join(plan))
        for phase in plan:
            self._mark_dashboard_phase_start(
                dashboard_writer, phase, merged_options.get(phase, {})
            )
            self._presenter.phase_start(phase)
            result = self._run_phase(phase, merged_options.get(phase, {}))
            self.results.append(result)
            self._presenter.phase_end(phase, status=result.status)
            self._mark_dashboard_phase_end(dashboard_writer, phase, result)
            if result.status == "error":
                logger.error("Phase '%s' failed: %s", phase, result.message)
                break

        if dashboard_writer is not None:
            failed = next((item for item in self.results if item.status == "error"), None)
            dashboard_writer.write_status(
                status="failed" if failed else "completed",
                latest_error=failed.message if failed else None,
            )
            dashboard_writer.event(
                "FastMDXplora workflow failed" if failed else "FastMDXplora workflow completed",
                level="error" if failed else "info",
            )

        self._write_manifest()
        self._write_resolved_config()
        self._presenter.done()

        # Wrap the phase results into a single RunResult (a study of one).
        status = "error" if any(r.status == "error" for r in self.results) else "ok"
        return [RunResult(
            run_id="s1",
            system=self.system,
            status=status,
            output_dir=self.output_dir,
            phases=list(self.results),
        )]

    def _print_dry_run_single(self, plan: list[str]) -> None:
        """Print the plan for a single study without running it."""
        print("\nFastMDXplora dry run (no execution)")
        print("=" * 40)
        print(f"  system:  {self.system}")
        print(f"  output:  {self.output_dir}")
        print(f"  phases:  {' → '.join(plan) if plan else '(none)'}")

    def _explore_config(
        self,
        *,
        include: list[str] | None,
        exclude: list[str] | None,
        report: bool,
        dry_run: bool = False,
    ) -> list[RunResult]:
        """Run a config-driven study through the internal batch machinery.

        Handles one system or many identically. Exposed to the user only as
        ``FastMDXplora(config=...).explore()``; the batch layer underneath
        is private.
        """
        from fastmdxplora.batch import BatchExplorer

        batch = BatchExplorer(
            config=self._config_path,
            config_data=self._config_data,
            output_dir=self._deferred_output_dir,
            verbose=self._deferred_verbose,
        )
        # explore()-level phase overrides win over the config file.
        if include is not None:
            batch._raw["include"] = include
            batch._raw["exclude"] = None
        elif exclude is not None:
            batch._raw["exclude"] = exclude
            batch._raw["include"] = None
        if not report:
            existing = batch._raw.get("exclude") or []
            if "report" not in existing and not batch._raw.get("include"):
                batch._raw["exclude"] = [*existing, "report"]

        if dry_run:
            run_results = batch.dry_run()
            self.output_dir = batch.output_dir
            self.results = run_results
            return run_results

        run_results = batch.run()
        # Surface the resolved output location for callers that read it.
        self.output_dir = batch.output_dir
        self.results = run_results
        return run_results

    # Convenience: per-phase entry points (also called by the CLI)
    def setup(self, **kwargs: Any) -> PhaseResult:
        """Run only the setup phase."""
        return self._run_phase("setup", kwargs)

    def simulate(self, **kwargs: Any) -> PhaseResult:
        """Run only the simulation phase."""
        return self._run_phase("simulation", kwargs)

    def analyze(self, **kwargs: Any) -> PhaseResult:
        """Run only the analysis phase."""
        return self._run_phase("analysis", kwargs)

    def report(self, **kwargs: Any) -> PhaseResult:
        """Run only the report phase."""
        return self._run_phase("report", kwargs)

    def compare(self, *, output_dir: str | os.PathLike | None = None) -> Path | None:
        """(Re)build the cross-run comparison report for a multi-run study.

        A multi-run ``explore()`` builds this automatically; call this to
        regenerate it — for example after re-running some of the runs, or
        to produce it for a batch that finished earlier.

        Parameters
        ----------
        output_dir : str | os.PathLike, optional
            The batch output directory to read (the one containing
            ``batch_manifest.json``). Defaults to this object's
            ``output_dir`` — i.e. the study it just ran.

        Returns
        -------
        Path or None
            The ``comparison/`` directory, or None if there was nothing to
            compare (fewer than two successful runs, or no analysis
            outputs were found).
        """
        from fastmdxplora.batch.compare import build_comparison_report

        target = Path(output_dir) if output_dir is not None else getattr(
            self, "output_dir", None
        )
        if target is None:
            raise ValueError(
                "compare() needs an output directory — pass output_dir=, or "
                "call it after explore() so the run's output is known."
            )
        return build_comparison_report(target)

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------
    def _dashboard_writer(
        self,
        merged_options: dict[str, dict[str, Any]] | None = None,
    ):
        """Return the project-level telemetry writer when a dashboard is active."""

        simulation_options = (
            merged_options or self.options
        ).get("simulation", {})

        dashboard_active = (
            os.getenv("FASTMDX_DASHBOARD_ACTIVE") == "1"
        )
        dashboard_output = os.getenv("FASTMDX_DASHBOARD_OUTPUT")

        dashboard_matches_run = False

        if dashboard_active and dashboard_output:
            try:
                dashboard_matches_run = (
                    Path(dashboard_output).expanduser().resolve()
                    == self.output_dir.expanduser().resolve()
                )
            except (OSError, RuntimeError):
                dashboard_matches_run = False

        active = bool(
            simulation_options.get("live_telemetry")
            or dashboard_matches_run
        )

        if not active:
            return None

        try:
            from fastmdxplora.live.telemetry import TelemetryWriter

            return TelemetryWriter(
                self.output_dir / "simulation",
                enabled=True,
            )
        except Exception as exc:  # noqa: BLE001
            logger.debug(
                "Could not initialize dashboard phase telemetry: %s",
                exc,
            )
            return None

    @staticmethod
    def _initialize_dashboard_timeline(
        writer: Any,
        plan: list[str],
    ) -> None:
        states = {
            "setup": "waiting" if "setup" in plan else "skipped",
            "minimization": (
                "waiting" if "simulation" in plan else "skipped"
            ),
            "nvt": "waiting" if "simulation" in plan else "skipped",
            "npt": "waiting" if "simulation" in plan else "skipped",
            "production": (
                "waiting" if "simulation" in plan else "skipped"
            ),
            "analysis": "waiting" if "analysis" in plan else "skipped",
            "report": "waiting" if "report" in plan else "skipped",
        }

        writer.write_status(
            stage=plan[0] if plan else "completed",
            status="running" if plan else "completed",
            stage_states=states,
        )
        writer.event("Workflow timeline initialized")

    @staticmethod
    def _mark_dashboard_phase_start(
        writer: Any,
        phase: str,
        options: dict[str, Any],
    ) -> None:
        if writer is None:
            return

        if phase == "simulation":
            first_stage = (
                "minimization"
                if options.get("minimize", True)
                else "nvt"
            )
            writer.mark_stage(
                first_stage,
                "current",
                status="running",
            )
        else:
            writer.mark_stage(
                phase,
                "current",
                status="running",
            )

        writer.event(f"{phase.title()} phase started")

    @staticmethod
    def _mark_dashboard_phase_end(
        writer: Any,
        phase: str,
        result: PhaseResult,
    ) -> None:
        if writer is None:
            return

        if result.status == "ok":
            state = "completed"
        elif result.status == "skipped":
            state = "skipped"
        else:
            state = "failed"

        if phase == "simulation":
            from fastmdxplora.live.telemetry import read_status

            status = read_status(writer.root.parent)
            stages = (
                status.get("stage_states")
                if isinstance(status, dict)
                else {}
            )
            stages = stages if isinstance(stages, dict) else {}

            for name in (
                "minimization",
                "nvt",
                "npt",
                "production",
            ):
                current_state = str(
                    stages.get(name, "waiting")
                ).lower()

                if current_state in {"waiting", "current"}:
                    writer.mark_stage(
                        name,
                        state,
                        status=(
                            "failed"
                            if state == "failed"
                            else "running"
                        ),
                    )
        else:
            writer.mark_stage(
                phase,
                state,
                status=(
                    "failed"
                    if state == "failed"
                    else "running"
                ),
                latest_error=(
                    result.message
                    if state == "failed"
                    else None
                ),
            )

        writer.event(
            f"{phase.title()} phase {state}",
            level="error" if state == "failed" else "info",
        )
        
    def _build_plan(
        self,
        *,
        include: list[str] | None,
        exclude: list[str] | None,
        want_report: bool,
    ) -> list[str]:
        if include is not None and exclude is not None:
            raise ValueError("Specify either `include` or `exclude`, not both.")

        if include is not None:
            unknown = set(include) - set(PHASES)
            if unknown:
                raise ValueError(f"Unknown phase(s): {sorted(unknown)}. Valid: {PHASES}")
            plan = [p for p in PHASES if p in include]
        elif exclude is not None:
            unknown = set(exclude) - set(PHASES)
            if unknown:
                raise ValueError(f"Unknown phase(s): {sorted(unknown)}. Valid: {PHASES}")
            plan = [p for p in PHASES if p not in exclude]
        else:
            plan = list(PHASES)

        if not want_report and "report" in plan:
            plan.remove("report")

        return plan

    def _merge_options(
        self, override: dict[str, dict[str, Any]] | None
    ) -> dict[str, dict[str, Any]]:
        merged: dict[str, dict[str, Any]] = {p: dict(self.options.get(p, {})) for p in PHASES}
        if override:
            for phase, opts in override.items():
                if phase not in PHASES:
                    raise ValueError(
                        f"Unknown phase '{phase}' in options. Valid: {PHASES}"
                    )
                merged[phase].update(opts)
        return merged

    def _run_phase(self, phase: str, kwargs: dict[str, Any]) -> PhaseResult:
        phase_dir = self._phase_dirs[phase]
        phase_dir.mkdir(parents=True, exist_ok=True)

        started = datetime.now(timezone.utc).isoformat()
        logger.debug("--> Phase '%s' starting (output=%s)", phase, phase_dir)

        try:
            run_fn = self._resolve_phase_runner(phase)
            artifacts = run_fn(
                orchestrator=self,
                output_dir=phase_dir,
                **kwargs,
            )
            finished = datetime.now(timezone.utc).isoformat()
            return PhaseResult(
                name=phase,
                status="ok",
                output_dir=phase_dir,
                started_at=started,
                finished_at=finished,
                message=f"Phase '{phase}' completed.",
                artifacts=list(artifacts or []),
            )
        except Exception as exc:  # noqa: BLE001 -- we log and record
            finished = datetime.now(timezone.utc).isoformat()
            logger.exception("Phase '%s' raised an exception", phase)
            return PhaseResult(
                name=phase,
                status="error",
                output_dir=phase_dir,
                started_at=started,
                finished_at=finished,
                message=str(exc),
            )

    @staticmethod
    def _resolve_phase_runner(phase: str):
        """Look up the run() entry point for a given phase.

        Each phase package exposes a ``run(orchestrator, output_dir, **kwargs)``
        callable; the orchestrator imports it lazily so that an optional
        backend (e.g. OpenMM) is only required when its phase is invoked.
        """
        if phase == "setup":
            from fastmdxplora.setup.pipeline import run

            return run
        if phase == "simulation":
            from fastmdxplora.simulation.pipeline import run

            return run
        if phase == "analysis":
            from fastmdxplora.analysis.analyze import run

            return run
        if phase == "report":
            from fastmdxplora.report import run

            return run
        raise ValueError(f"Unknown phase: {phase}")

    def _write_manifest(self) -> None:
        """Write a single JSON manifest summarizing this session."""
        from fastmdxplora import __citation__, __doi__, __version__

        manifest = {
            "tool": "FastMDXplora",
            "version": __version__,
            "doi": __doi__,
            "citation": __citation__,
            "system": self.system,
            "output_dir": str(self.output_dir),
            "phases": [r.to_dict() for r in self.results],
            "options": self.options,
        }
        manifest_path = self.output_dir / "manifest.json"
        with manifest_path.open("w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2)
        logger.debug("Wrote manifest: %s", manifest_path)

    def _write_resolved_config(self) -> None:
        """Write the fully-merged configuration for reproducibility.

        Produces ``resolved_config.yml`` capturing the system, output,
        phase selection, and per-phase options actually used. The file is
        a valid FastMDXplora config — feeding it back to ``--config``
        reproduces the run.
        """
        from fastmdxplora.config import write_resolved_config

        resolved = {
            "system": self.system,
            "output": str(self.output_dir),
            "verbose": self.verbose,
            "include": getattr(self, "_resolved_include", None) or self._config_include,
            "exclude": getattr(self, "_resolved_exclude", None) or self._config_exclude,
            "options": getattr(self, "_resolved_options", None) or self.options,
        }
        try:
            path = write_resolved_config(resolved, self.output_dir)
            logger.debug("Wrote resolved config: %s", path)
        except Exception as exc:  # noqa: BLE001 -- never fail a run over this
            logger.debug("Could not write resolved config: %s", exc)

    def _configure_logging(self) -> None:
        """Wire up console and file logging for this project session.

        The root ``fastmdx`` logger is set to DEBUG so all records flow to
        the handlers; each handler then applies its own level filter. The
        file handler always captures at DEBUG (full audit trail). The
        console handler defaults to INFO, raised to DEBUG when
        ``verbose=True`` or ``FASTMDX_LOGLEVEL=DEBUG`` is set.
        """
        from fastmdxplora.utils.logging import attach_file_logger, set_level, setup_console

        console_level = logging.DEBUG if self.verbose else logging.INFO
        setup_console(level=console_level)
        attach_file_logger(self.output_dir / "fastmdxplora.log", level=logging.DEBUG)
        # Root logger must be at the lowest handler level so records flow.
        set_level(logging.DEBUG)
        # ...but the console handler still applies its own filter.
        # set_level above promoted ALL handlers to DEBUG; re-apply the
        # console-level filter so quiet mode stays quiet.
        from fastmdxplora.utils.logging import _console_handler

        if _console_handler is not None:
            _console_handler.setLevel(console_level)

    def _configure_presenter(self):
        """Create the session presenter for structural output.

        The presenter is silent when ``FASTMDX_LOG_STYLE=plain`` (handled
        internally by :class:`SessionPresenter`) or when stdout is not a
        TTY (handled by color auto-detection). Users wanting different
        behaviour can replace ``self._presenter`` after construction.
        """
        from fastmdxplora.utils.presenter import SessionPresenter

        return SessionPresenter()
