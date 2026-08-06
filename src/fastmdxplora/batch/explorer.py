"""Batch orchestration — the single execution path for all runs.

Every FastMDXplora run goes through :class:`BatchExplorer`. There is no
separate "single run" path: a one-system config is simply a batch of
one. This keeps one code path, one config shape (``systems:`` is always
a list), and one mental model.

Output layout adapts to the run count:

  - **One run** → flat, familiar layout written directly to the output
    directory: ``output/setup/``, ``output/simulation/``, etc., with the
    usual ``manifest.json`` and ``resolved_config.yml``. No ``runs/``
    wrapper, no ``batch_manifest.json``.
  - **Many runs** → each run in ``output/runs/<id>/`` (a complete study),
    plus a top-level ``batch_manifest.json`` indexing them all.

Execution modes (``execution:`` block):

  - **sequential** (default) — one run at a time, in process.
  - **parallel** — a process pool of ``workers`` runs at once. On GPU,
    set ``devices: [0, 1, ...]`` and each worker is pinned to a distinct
    device round-robin (one run per GPU), which is the only safe way to
    parallelize GPU MD — oversubscribing a single GPU is slower than
    sequential.

Process-based (not thread-based) parallelism is mandatory: OpenMM
contexts and the GIL don't share across threads. Each run is therefore
dispatched to a subprocess via a module-level worker function.
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any

from fastmdxplora.batch.sweep import (
    RunSpec,
    expand_runs,
    normalize_sweep,
    normalize_systems,
)
from fastmdxplora.config import load_config_file, validate_config
from fastmdxplora.utils.logging import get_logger

if TYPE_CHECKING:
    from fastmdxplora.orchestrator import RunResult

logger = get_logger("batch")


def _a_prepared_system_is_there(setup_dir: Path) -> bool:
    """Whether the three files a simulation starts from are on disk.

    Asked through the simulation phase's own check, so what counts as
    prepared is decided in one place. Its answer is a tuple of paths, and an
    empty answer is a tuple of Nones -- which is truthy, so it is the first
    element that has to be read.
    """
    from fastmdxplora.simulation.pipeline import _setup_outputs_present

    return _setup_outputs_present(setup_dir)[0] is not None


def _first_error_phase_message(phases: list[Any]) -> str:
    for phase in phases:
        if getattr(phase, "status", None) == "error":
            return getattr(phase, "message", "") or f"Phase '{phase.name}' failed."
    return ""


def _skipped_run_result(spec: RunSpec, run_out: Path, message: str) -> "RunResult":
    from fastmdxplora.orchestrator import RunResult

    return RunResult(
        run_id=spec.run_id,
        system=spec.system,
        status="skipped",
        output_dir=run_out,
        sweep_values=spec.sweep_values,
        phases=[],
        message=message,
    )


# ---------------------------------------------------------------------------
# Module-level worker (must be top-level so ProcessPoolExecutor can pickle it)
# ---------------------------------------------------------------------------
def _execute_run(
    spec_dict: dict[str, Any],
    run_out: str,
    include: list[str] | None,
    exclude: list[str] | None,
    verbose: bool,
    device_override: str | None,
) -> "RunResult":
    """Run one study and return a RunResult. Safe to call in a subprocess.

    Takes plain dicts/strings (picklable) rather than RunSpec objects so it
    works cleanly across the process boundary. Returns a RunResult, which
    is a dataclass and pickles back cleanly.
    """
    # Spawned workers (Windows uses 'spawn') do not pass through the CLI
    # entry point, so their stdout/stderr are not reconfigured to UTF-8 and
    # default to the platform codec (cp1252 on Windows). The presenter prints
    # non-ASCII status glyphs (✓, ▸, box-drawing), which then raise
    # UnicodeEncodeError and fail the run. Reconfigure here, as main() does.
    import sys as _sys

    for _stream in (_sys.stdout, _sys.stderr):
        try:
            _stream.reconfigure(encoding="utf-8")  # type: ignore[union-attr]
        except (AttributeError, ValueError):
            pass

    # Several workers share one terminal, so none of them should print a
    # banner to it. Three copies of the wordmark arriving a character at a
    # time is not three runs starting -- it is one unreadable screen. Each
    # run writes its own log file, and the batch prints the progress that
    # belongs to the whole study.
    # Restored when the run ends, because this function is called in-process
    # by the tests as well as in a subprocess by a real run, and a worker
    # that leaves a global changed poisons everything after it.
    import os as _os
    import sys as _sys

    class _QuietBanner:
        """Keep a worker's output out of the shared terminal.

        Two sources. Ours obeys an environment variable. PLUMED's comes from
        C++ and writes to the file descriptor directly, so redirecting
        Python's ``sys.stdout`` does not reach it -- the same lesson as a
        Markdown renderer that printed its installation guide from a dynamic
        loader. The descriptor itself has to move.

        It moves to the run's own log, so nothing is lost: a worker's output
        belongs beside its results rather than mixed with two others'.
        """

        def __init__(self, log_path):
            self.log_path = log_path

        def __enter__(self):
            self.previous = _os.environ.get("FASTMDX_LOG_STYLE")
            _os.environ["FASTMDX_LOG_STYLE"] = "plain"

            self.saved = None
            try:
                _os.makedirs(_os.path.dirname(self.log_path), exist_ok=True)
                self.sink = open(self.log_path, "a", encoding="utf-8")
                self.saved = (_os.dup(1), _os.dup(2))
                _os.dup2(self.sink.fileno(), 1)
                _os.dup2(self.sink.fileno(), 2)
            except OSError:
                # A worker that cannot redirect should still run.
                self.saved = None
            return self

        def __exit__(self, *_exc):
            if self.saved is not None:
                # Flush while the descriptors are still moved. PLUMED's
                # writes sit in the C runtime's buffer until something
                # empties it, and whatever is still there when the
                # descriptors go back lands in the shared terminal -- the
                # leak this guard exists to prevent, arriving late. The
                # trajectory loader learned the same thing from MDTraj's
                # plugin messages on macOS.
                from fastmdxplora.utils.native_output import _LIBC

                for stream in (_sys.stdout, _sys.stderr):
                    try:
                        stream.flush()
                    except (ValueError, OSError):
                        pass
                if _LIBC is not None:
                    try:
                        _LIBC.fflush(None)  # fflush(NULL): every open stream
                    except OSError:
                        pass

                out, err = self.saved
                _os.dup2(out, 1)
                _os.dup2(err, 2)
                _os.close(out)
                _os.close(err)
                self.sink.close()
            if self.previous is None:
                _os.environ.pop("FASTMDX_LOG_STYLE", None)
            else:
                _os.environ["FASTMDX_LOG_STYLE"] = self.previous
            return False

    # Imported here so the subprocess sets up its own logging/imports.
    from fastmdxplora import FastMDXplora
    from fastmdxplora.orchestrator import RunResult

    options = dict(spec_dict["options"])
    # Round-robin GPU pinning: stamp this worker's device onto the run.
    if device_override is not None:
        sim = dict(options.get("simulation", {}))
        sim["device_index"] = device_override
        options["simulation"] = sim

    try:
        with _QuietBanner(_os.path.join(run_out, "run.log")):
            fmdx = FastMDXplora(
                system=spec_dict["system"],
                output_dir=run_out,
                options=options,
                verbose=verbose,
            )
            # A single-system explore() returns a one-element list of
            # RunResult; take its phases and re-stamp the run's identity.
            inner = fmdx.explore(include=include, exclude=exclude, report=True)
        phases = inner[0].phases if inner else []
        status = "error" if any(p.status == "error" for p in phases) else "ok"
        message = _first_error_phase_message(phases) if status == "error" else ""
        return RunResult(
            run_id=spec_dict["run_id"],
            system=spec_dict["system"],
            status=status,
            output_dir=Path(run_out),
            sweep_values=spec_dict["sweep_values"],
            phases=phases,
            message=message,
            error_type="PhaseError" if status == "error" else None,
        )
    except Exception as exc:  # noqa: BLE001 -- isolate per-run failures
        return RunResult(
            run_id=spec_dict["run_id"],
            system=spec_dict["system"],
            status="error",
            output_dir=Path(run_out),
            sweep_values=spec_dict["sweep_values"],
            phases=[],
            message=f"{type(exc).__name__}: {exc}",
            error_type=type(exc).__name__,
        )


def _not_submitted(after: str, *, umbrella: bool) -> str:
    """Why a run was not started.

    Umbrella windows are one measurement, so a failure ends the study rather
    than costing it a data point -- and saying so is the difference between a
    reader thinking something crashed and knowing the rest was not attempted
    on purpose.
    """
    if umbrella:
        return (
            f"Not submitted: window '{after}' failed, and a free energy needs "
            "every window -- the rest would be hours spent on runs that get "
            "thrown away."
        )
    return (
        f"Not submitted because continue_on_error=False stopped after failed "
        f"run '{after}'."
    )


class BatchExplorer:
    """Run one or more FastMDXplora studies (systems × sweep).

    Parameters
    ----------
    config : str | os.PathLike
        Path to a YAML config with a ``systems:`` list (and optionally
        ``sweep:`` / ``execution:``).
    output_dir : str | os.PathLike | None
        Root output directory. One run → written here directly; many runs
        → each in ``runs/<id>/``. Defaults to a timestamped directory.
    verbose : bool
        Forwarded to each run.
    continue_on_error : bool | None
        Override the config's ``execution.continue_on_error``. If None,
        the config value (default True) is used.

    Examples
    --------
    >>> BatchExplorer(config="study.yml").run()       # doctest: +SKIP
    """

    def __init__(
        self,
        config: str | os.PathLike | None = None,
        *,
        config_data: dict[str, Any] | None = None,
        output_dir: str | os.PathLike | None = None,
        verbose: bool = False,
        continue_on_error: bool | None = None,
    ) -> None:
        if config is None and config_data is None:
            raise ValueError("BatchExplorer requires `config` (path) or `config_data` (dict).")

        if config_data is not None:
            self.config_path = str(config) if config is not None else "<in-memory>"
            raw = dict(config_data)
        else:
            self.config_path = str(config)
            raw = load_config_file(self.config_path)
        validate_config(raw, require_systems=True)
        self._raw = raw
        self.verbose = verbose

        # Execution settings
        execution = raw.get("execution") or {}
        self.mode = execution.get("mode", "sequential")
        self.workers = execution.get("workers")
        self.devices = execution.get("devices")
        # Umbrella windows are one measurement, not several systems. A
        # campaign should carry on when one system fails -- the others are
        # still results. A free energy cannot be computed at all if a window
        # is missing, so continuing spends hours producing runs that will be
        # thrown away.
        from fastmdxplora.simulation.umbrella import plan_from_expanded

        self._is_umbrella = plan_from_expanded(raw) is not None

        self.continue_on_error = (
            continue_on_error if continue_on_error is not None
            else execution.get(
                "continue_on_error", not self._is_umbrella)
        )
        # The cross-run comparison report is a reporting artifact, so it's
        # controlled by the `report` block (default on).
        report_block = raw.get("report") or {}
        self.comparison = report_block.get("comparison", True)

        # Output root
        resolved_output = output_dir or raw.get("output")
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
        self.output_dir = (
            Path(resolved_output) if resolved_output
            else Path(f"fastmdxplora_output_{timestamp}")
        )

        # Expand the run matrix
        self.run_specs = self._build_run_specs(raw)
        self.is_single = len(self.run_specs) == 1
        self.results: list[dict[str, Any]] = []

    # ------------------------------------------------------------------
    def _build_run_specs(self, raw: dict[str, Any]) -> list[RunSpec]:
        from fastmdxplora.config import phase_options

        base_options = phase_options(raw)

        systems = normalize_systems(raw["systems"])
        sweep = (
            normalize_sweep(raw["sweep"]) if raw.get("sweep") is not None else None
        )
        return expand_runs(systems=systems, sweep=sweep, base_options=base_options)

    # ------------------------------------------------------------------
    def _run_output_dir(self, spec: RunSpec) -> Path:
        """Flat output for a single run; runs/<id>/ for many."""
        if self.is_single:
            return self.output_dir
        return self.output_dir / "runs" / spec.run_id

    def _resolve_workers(self) -> int:
        """Decide the parallel worker count."""
        if self.workers:
            return max(1, int(self.workers))
        if self.devices:
            return max(1, len(self.devices))
        # CPU default: all cores, capped at the number of runs.
        return max(1, min(os.cpu_count() or 1, len(self.run_specs)))

    def _device_for_worker(self, worker_slot: int) -> str | None:
        """Round-robin GPU device for a worker slot, or None if no devices."""
        if not self.devices:
            return None
        dev = self.devices[worker_slot % len(self.devices)]
        return str(dev)

    # ------------------------------------------------------------------
    def run(self) -> list[dict[str, Any]]:
        """Execute every run. Returns per-run result records."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        n = len(self.run_specs)
        include = self._raw.get("include")
        exclude = self._raw.get("exclude")

        shared = self._maybe_prepare_once(include, exclude)
        if shared is not None:
            # The windows read the system that was just prepared, so none of
            # them prepares one of its own.
            exclude = list(exclude or []) + ["setup"]

        if self.is_single:
            logger.info("Exploring 1 molecular system in %s", self.output_dir)
        else:
            (self.output_dir / "runs").mkdir(exist_ok=True)
            logger.info(
                "Exploring %d molecular systems in %s (mode=%s)", n, self.output_dir, self.mode
            )
            print(f"\nFastMDXplora: exploring {n} systems ({self.mode})\n{'=' * 40}")

        if self.mode == "parallel" and not self.is_single:
            self.results = self._run_parallel(include, exclude)
        else:
            self.results = self._run_sequential(include, exclude)

        # Only write a batch manifest when there's actually a batch.
        if not self.is_single:
            self._write_batch_manifest()
            self._maybe_build_comparison()
            self._maybe_build_pmf()
            self._print_summary()
        return list(self.results)

    # ------------------------------------------------------------------
    def _maybe_prepare_once(self, include, exclude) -> Path | None:
        """Prepare the system once for a set of umbrella windows.

        Seven windows of one real study came out with 37,212, 37,254, 37,436
        and 37,445 atoms: four different systems for one measurement. The
        windows are the same molecule held at different points along a
        coordinate, so they should be the same molecule. Solvation does not
        place water the same way twice, and water arranged differently
        between windows is noise in the free energy rather than physics.

        A window is not a system of its own, so it does not get a system of
        its own. Preparing seven times was also seven times the work.

        Only for a study whose runs are one measurement. An ordinary campaign
        is untouched: different systems *are* different systems, and one
        preparation shared between them would be wrong rather than efficient.
        """
        if not self._is_umbrella:
            return None

        wanted = set(include) if include else {"setup", "simulation",
                                               "analysis", "report"}
        if "setup" not in wanted or "setup" in set(exclude or []):
            # Nothing is being prepared at all -- the windows are simulating
            # from something that already exists.
            return None

        # An umbrella study is one system by the time it gets here -- expansion
        # refuses several -- but a sweep can still ask for the windows to be
        # prepared differently from one another. Sharing then would quietly
        # ignore what the sweep asked for, which is worse than preparing
        # seven times.
        preparations = {
            json.dumps(spec.options.get("setup") or {}, sort_keys=True, default=str)
            for spec in self.run_specs
        }
        if len(preparations) != 1:
            logger.warning(
                "The windows are swept over %d different setup settings, so "
                "each prepares its own system.", len(preparations))
            return None

        shared = self.output_dir / "shared_setup"
        prepared = shared / "setup"
        if not _a_prepared_system_is_there(prepared):
            print("\nPreparing one system for every window\n" + "=" * 40)
            first = self.run_specs[0]
            result = _execute_run(
                first.to_dict(), str(shared), ["setup"], None,
                self.verbose, None,
            )
            if result.status == "error" or not _a_prepared_system_is_there(prepared):
                # Every window would fail the same way, one after another,
                # for hours. Say it once here instead.
                raise RuntimeError(
                    "The system could not be prepared, so there is nothing "
                    f"for the windows to simulate: {result.message or prepared}"
                )

        for spec in self.run_specs:
            simulation = dict(spec.options.get("simulation") or {})
            simulation["prepared_from"] = str(prepared)
            spec.options["simulation"] = simulation
        return prepared

    # ------------------------------------------------------------------
    def _maybe_build_pmf(self) -> None:
        """Recombine umbrella windows once they have all run.

        The windows are ordinary runs, so nothing before this point knows
        they belong together. Here is where they have all finished and the
        set is visible -- and where the overlap between them can be checked,
        which is the thing that decides whether a free energy exists at all.
        """
        import json

        from fastmdxplora.simulation.umbrella import (
            collect_samples,
            compute_pmf,
            plan_from_expanded,
        )

        # Rebuilt from the runs rather than carried in the config: a config is
        # what a user wrote, not a place to hide state, and validation
        # rightly refused an extra key when this was smuggled through one.
        plan = plan_from_expanded(self._raw or {})
        if plan is None:
            return

        # Where the runs actually went, from the same method that put them
        # there. Reconstructing the path from the output directory got both
        # the "runs" level and the slugged identifier wrong.
        directories = {}
        for spec in self.run_specs:
            block = (spec.options.get("simulation") or {}).get("umbrella") or {}
            if block.get("index") is not None:
                directories[int(block["index"])] = self._run_output_dir(spec)

        try:
            samples = collect_samples(directories)
        except FileNotFoundError as exc:
            # Some window did not produce sampling. Recorded rather than
            # raised: the runs that did work are still on disk and worth
            # keeping.
            payload = {"pmf": None, "refused": str(exc),
                       "plan": plan.as_record()}
        else:
            payload = compute_pmf(
                samples, plan,
                temperature_K=float(
                    ((self._raw or {}).get("simulation") or {})
                    .get("temperature_K", 300.0)))
            payload["plan"] = plan.as_record()

        destination = Path(self.output_dir) / "pmf.json"
        destination.write_text(json.dumps(payload, indent=2), encoding="utf-8")

        presenter = getattr(self, "_presenter", None) or getattr(self, "presenter", None)
        if presenter is not None:
            if payload.get("refused"):
                presenter.step(
                    "No free energy: " + payload["refused"].split(".")[0],
                    status="warn")
            else:
                presenter.step(f"Wrote {destination.name}")

    def _maybe_build_comparison(self) -> None:
        """Build the cross-run comparison report (best-effort).

        Controlled by ``execution.comparison`` (default True). A failure
        here must never fail the batch — the runs themselves succeeded.
        """
        if not self.comparison:
            return
        n_ok = sum(1 for r in self.results if r.status == "ok")
        if n_ok < 2:
            return
        try:
            from fastmdxplora.batch.compare import build_comparison_report

            path = build_comparison_report(self.output_dir)
            if path is not None:
                print(f"Comparison:     {path / 'comparison_report.md'}")
        except Exception as exc:  # noqa: BLE001 -- never break the batch
            logger.warning("Comparison report failed (runs are unaffected): %s", exc)

    # ------------------------------------------------------------------
    def dry_run(self) -> list["RunResult"]:
        """Report the plan without executing anything.

        Prints each run, its system, swept values, target output directory,
        and the phases that would execute, then returns a list of
        ``RunResult`` with status ``"planned"`` (no phases populated).
        """
        from fastmdxplora.orchestrator import RunResult, PHASES

        include = self._raw.get("include")
        exclude = self._raw.get("exclude")
        # Compute the phase plan the same way the orchestrator would.
        if include:
            plan = [p for p in PHASES if p in include]
        elif exclude:
            plan = [p for p in PHASES if p not in exclude]
        else:
            plan = list(PHASES)

        n = len(self.run_specs)
        layout = "flat" if self.is_single else "runs/<id>/"
        print("\nFastMDXplora dry run (no execution)")
        print("=" * 50)
        print(f"  runs:    {n}")
        print(f"  output:  {self.output_dir}  ({layout} layout)")
        print(f"  phases:  {' → '.join(plan) if plan else '(none)'}")
        if self.mode == "parallel":
            print(f"  mode:    parallel ({self._resolve_workers()} workers"
                  + (f", devices={self.devices}" if self.devices else "") + ")")
        else:
            print(f"  mode:    sequential")
        print("-" * 50)

        planned: list[RunResult] = []
        for i, spec in enumerate(self.run_specs, start=1):
            run_out = self._run_output_dir(spec)
            label = spec.run_id if not self.is_single else spec.system
            sweep = ""
            if spec.sweep_values:
                sweep = "  [" + ", ".join(
                    f"{k}={v}" for k, v in spec.sweep_values.items()) + "]"
            print(f"  [{i}/{n}] {label}{sweep}")
            print(f"          → {run_out}")
            planned.append(RunResult(
                run_id=spec.run_id, system=spec.system, status="planned",
                output_dir=run_out, sweep_values=spec.sweep_values, phases=[],
            ))
        print("=" * 50)
        return planned

    # ------------------------------------------------------------------
    def _run_sequential(self, include, exclude) -> list["RunResult"]:
        results: list[RunResult] = []
        n = len(self.run_specs)
        for i, spec in enumerate(self.run_specs, start=1):
            run_out = self._run_output_dir(spec)
            if not self.is_single:
                print(f"\n[{i}/{n}] {spec.run_id}")
                if spec.sweep_values:
                    pretty = ", ".join(f"{k}={v}" for k, v in spec.sweep_values.items())
                    print(f"        {pretty}")
            # Sequential: a single device (first listed) may still be pinned.
            device = self._device_for_worker(0) if self.devices else None
            result = _execute_run(
                spec.to_dict(), str(run_out), include, exclude,
                self.verbose, device,
            )
            results.append(result)
            if result.status == "error" and not self.continue_on_error:
                logger.error("Stopping after failed run '%s'.", spec.run_id)
                for skipped in self.run_specs[i:]:
                    results.append(_skipped_run_result(
                        skipped,
                        self._run_output_dir(skipped),
                        _not_submitted(
                            spec.run_id,
                            umbrella=getattr(self, "_is_umbrella", False),
                        ),
                    ))
                break
        return results

    # ------------------------------------------------------------------
    def _run_parallel(self, include, exclude) -> list["RunResult"]:
        from fastmdxplora.orchestrator import RunResult

        n_workers = self._resolve_workers()
        n = len(self.run_specs)
        print(f"Parallel execution: {n_workers} worker(s)"
              + (f", devices={self.devices}" if self.devices else ""))

        results: list[RunResult] = []
        futures = {}
        next_index = 0
        done = 0
        stopped_after: str | None = None

        def submit_next(pool: ProcessPoolExecutor) -> bool:
            nonlocal next_index
            if next_index >= n:
                return False
            spec = self.run_specs[next_index]
            run_out = self._run_output_dir(spec)
            device = self._device_for_worker(next_index)
            fut = pool.submit(
                _execute_run,
                spec.to_dict(), str(run_out), include, exclude,
                self.verbose, device,
            )
            futures[fut] = (next_index, spec)
            next_index += 1
            return True

        pool = ProcessPoolExecutor(max_workers=n_workers)
        try:
            for _ in range(min(n_workers, n)):
                submit_next(pool)

            while futures:
                fut = next(as_completed(futures))
                _idx, spec = futures.pop(fut)
                done += 1
                try:
                    result = fut.result()
                except Exception as exc:  # noqa: BLE001
                    result = RunResult(
                        run_id=spec.run_id, system=spec.system, status="error",
                        output_dir=self._run_output_dir(spec),
                        sweep_values=spec.sweep_values, phases=[],
                        message=f"{type(exc).__name__}: {exc}",
                        error_type=type(exc).__name__,
                    )
                mark = "✓" if result.status == "ok" else "✗"
                print(f"[{done}/{n}] {mark} {spec.run_id}")
                results.append(result)

                if result.status == "error" and not self.continue_on_error:
                    stopped_after = spec.run_id
                    logger.error("Stopping parallel batch after failed run '%s'.", spec.run_id)
                    break

                if stopped_after is None:
                    submit_next(pool)

            if stopped_after is not None:
                # Do not submit new work. Running tasks cannot always be
                # cancelled, so collect anything already in flight and mark
                # never-submitted runs explicitly as skipped.
                for fut, (_idx, spec) in list(futures.items()):
                    if fut.cancel():
                        results.append(_skipped_run_result(
                            spec,
                            self._run_output_dir(spec),
                            (
                                "Cancelled because continue_on_error=False "
                                f"stopped after failed run '{stopped_after}'."
                            ),
                        ))
                        futures.pop(fut, None)

                for fut in as_completed(futures):
                    _idx, spec = futures[fut]
                    done += 1
                    try:
                        result = fut.result()
                    except Exception as exc:  # noqa: BLE001
                        result = RunResult(
                            run_id=spec.run_id, system=spec.system, status="error",
                            output_dir=self._run_output_dir(spec),
                            sweep_values=spec.sweep_values, phases=[],
                            message=f"{type(exc).__name__}: {exc}",
                            error_type=type(exc).__name__,
                        )
                    mark = "✓" if result.status == "ok" else "✗"
                    print(f"[{done}/{n}] {mark} {spec.run_id}")
                    results.append(result)

                for spec in self.run_specs[next_index:]:
                    results.append(_skipped_run_result(
                        spec,
                        self._run_output_dir(spec),
                        _not_submitted(
                            stopped_after,
                            umbrella=getattr(self, "_is_umbrella", False),
                        ),
                    ))
        finally:
            pool.shutdown(wait=True, cancel_futures=True)

        # Preserve deterministic (submission) order in the manifest.
        order = {spec.run_id: i for i, spec in enumerate(self.run_specs)}
        results.sort(key=lambda r: order.get(r.run_id, 0))
        return results

    # ------------------------------------------------------------------
    def _write_batch_manifest(self) -> None:
        from fastmdxplora import __citation__, __doi__, __version__

        manifest = {
            "tool": "FastMDXplora",
            "kind": "batch",
            "version": __version__,
            "doi": __doi__,
            "citation": __citation__,
            "config": self.config_path,
            "output_dir": str(self.output_dir),
            "n_runs": len(self.run_specs),
            "execution": {
                "mode": self.mode,
                "workers": self._resolve_workers() if self.mode == "parallel" else 1,
                "devices": self.devices,
            },
            "systems": [
                {"id": s["id"], "system": s["system"]}
                for s in normalize_systems(self._raw["systems"])
            ],
            "sweep": self._raw.get("sweep") or {},
            "runs": [r.to_dict() for r in self.results],
        }
        path = self.output_dir / "batch_manifest.json"
        with path.open("w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2, default=str)
        logger.debug("Wrote batch manifest: %s", path)

    # ------------------------------------------------------------------
    def _print_summary(self) -> None:
        ok = sum(1 for r in self.results if r.status == "ok")
        err = sum(1 for r in self.results if r.status == "error")
        skipped = sum(1 for r in self.results if r.status == "skipped")
        print(f"\n{'=' * 40}")
        print(
            f"Batch complete: {ok} ok, {err} error(s), "
            f"{skipped} skipped, {len(self.results)} total"
        )
        print(f"Batch output:   {self.output_dir}")
        print(f"Manifest:       {self.output_dir / 'batch_manifest.json'}")
