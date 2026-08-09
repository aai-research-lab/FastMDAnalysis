"""Abstract base class for all FastMDXplora analysis modules.

Every analysis (RMSD, RMSF, Rg, ...) subclasses :class:`Analysis` and
implements ``compute()`` and ``plot()``. The base class handles everything
else: output directory creation, options-manifest serialization, atom
selection resolution, status tracking, and the ``run()`` convenience method
that does ``compute() -> save_data() -> plot() -> save_figure()`` in order.

The contract is deliberately small:

  - ``compute(traj)`` returns a Python object (usually a numpy array or
    pandas DataFrame). It must be deterministic, side-effect-free, and
    inexpensive to call again with the same input.
  - ``plot(result, ax)`` draws onto a matplotlib Axes. It must not call
    ``plt.show()`` or close the figure — the caller controls that.
  - ``save_data(result, path)`` writes the computed result to disk. The
    default implementation handles numpy arrays and DataFrames; analyses
    with non-tabular output can override it.

Outputs land in ``<output_dir>/<analysis_name>/``:

  ::

    <output_dir>/
    └── rmsd/
        ├── rmsd.dat       # the numerical data
        ├── rmsd.png       # the figure
        └── options.json   # parameter manifest for this analysis
"""

from __future__ import annotations

import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd

from fastmdxplora.analysis.plotting import new_figure, save_figure
from fastmdxplora.utils.logging import get_logger

logger = get_logger("analysis.base")


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
@dataclass
class AnalysisResult:
    """Lightweight record of one analysis invocation.

    Returned by :meth:`Analysis.run` and aggregated by the orchestrator.
    Includes both the computed data (for in-memory consumers) and the
    on-disk artifact paths (for report generation and provenance).
    """

    name: str
    status: str  # "ok" | "error" | "skipped"
    data: Any = None
    output_dir: Path | None = None
    figure_path: Path | None = None
    data_path: Path | None = None
    options_path: Path | None = None
    artifacts: list[Path] = field(default_factory=list)
    message: str = ""
    started_at: str = ""
    finished_at: str = ""

    def to_dict(self) -> dict[str, Any]:
        """Serialize for inclusion in the analysis manifest."""
        return {
            "name": self.name,
            "status": self.status,
            "output_dir": str(self.output_dir) if self.output_dir else None,
            "figure_path": str(self.figure_path) if self.figure_path else None,
            "data_path": str(self.data_path) if self.data_path else None,
            "options_path": str(self.options_path) if self.options_path else None,
            "artifacts": [str(path) for path in self.artifacts],
            "message": self.message,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
        }


# ---------------------------------------------------------------------------
# Abstract base class
# ---------------------------------------------------------------------------
class Analysis(ABC):
    """Base class for a single trajectory analysis module.

    Subclasses must:
      - Set the class attribute ``name`` to a short identifier (e.g. ``"rmsd"``).
      - Implement :meth:`compute` which takes an MDTraj trajectory and
        returns the analysis result.
      - Implement :meth:`plot` which renders ``self.result`` onto an Axes.

    Subclasses may override :meth:`save_data` and :attr:`default_selection`.
    """

    #: Short identifier used in output paths and the manifest.
    name: str = "analysis"

    #: Whether ``selection`` means anything for this analysis. An analysis
    #: that decides its own atoms -- a protein-ligand measure works out both
    #: sides from the ligand's residue name, dihedrals from the backbone --
    #: has nothing to apply it to, and offering a control that does nothing is
    #: worse than not offering one: it looks like it worked.
    honours_selection: bool = True

    #: Default atom selection (MDTraj selection language). ``None`` means
    #: "use the whole trajectory". Subclasses override when an analysis
    #: only makes sense on a subset of atoms (e.g. RMSF on CA atoms).
    default_selection: str | None = None

    #: The per-frame scalar whose ensemble average means something, as
    #: ``(column, label)``. ``column`` names a column when ``compute()``
    #: returns a DataFrame and is ``None`` when it returns a bare per-frame
    #: array. ``None`` for the attribute itself -- the default -- means no
    #: weighted average of this analysis is offered.
    #:
    #: On a metadynamics run the trajectory is not a Boltzmann ensemble, so
    #: every mean reported from it is an average over a distribution the bias
    #: flattened on purpose. Declaring this lets that mean be recomputed
    #: against the deposited bias and reported beside the raw one. Analyses
    #: whose result is not one number per frame leave it unset: reweighting a
    #: clustering or a projection is a harder question than a weighted mean,
    #: and the report says so rather than guessing.
    reweightable: tuple[str | None, str] | None = None

    #: Whether ``compute()`` returns per-frame categorical labels, either as
    #: an array or as a mapping of method name to array. A population is a
    #: weighted count of an indicator, so it reweights exactly as a mean
    #: does, and how often a state is visited is the thing a biased run
    #: distorts most.
    #:
    #: What this does not correct is which states exist. The clustering was
    #: performed on the biased frames, so the groupings themselves are shaped
    #: by where the bias sent the system; reweighting says how often each was
    #: really visited, not that the right ones were found.
    reweightable_populations: bool = False

    #: Human-readable description used in figure titles.
    description: str = ""

    #: True for analyses that only apply to protein-ligand complexes (e.g.
    #: ligand pose RMSD). The orchestrator runs these automatically when a
    #: ligand is present and skips them otherwise. They can still be
    #: explicitly requested via ``include``.
    requires_ligand: bool = False

    #: Whether ``compute`` returns one value per frame.
    #:
    #: A quantity measured every frame has a mean, and a mean is not a
    #: measurement until two things are known: whether the system had settled
    #: by the time the averaging started, and how many *independent*
    #: observations the average rests on. Both are recorded automatically for
    #: an analysis that says yes here.
    #:
    #: Declared rather than inferred from the array's length. A per-atom
    #: result on a trajectory that happens to have as many frames as atoms
    #: would otherwise be summarised as though it were a time series, and the
    #: numbers would look right.
    time_series: bool = False

    def __init__(
        self,
        *,
        selection: str | None = None,
        output_dir: str | Path | None = None,
        title: str | None = None,
        xlabel: str | None = None,
        ylabel: str | None = None,
        figsize: tuple[float, float] | None = None,
        xunit: str | None = None,
        **options: Any,
    ) -> None:
        """Initialize the analysis with user-supplied options.

        Parameters
        ----------
        selection : str, optional
            MDTraj atom selection string. Overrides ``default_selection``.
        output_dir : path, optional
            Where to write outputs. The analysis appends its own ``name``
            subdirectory. If omitted, defaults to the current directory.
        title, xlabel, ylabel : str, optional
            Figure customization hooks. If provided, these override the
            defaults used by :meth:`figure_title`, :meth:`default_xlabel`,
            and :meth:`default_ylabel` respectively.
        figsize : (float, float), optional
            Figure size in inches.
        xunit : {"ns", "ps", "frames", None}, optional
            X-axis unit for time-series analyses (RMSD, Rg, etc.). The
            default is ``"ns"`` when the trajectory carries a timestep,
            else ``"frames"``. Analyses without a time axis ignore this.
        **options
            Analysis-specific keyword arguments. Subclasses access these
            via ``self.options``.
        """
        if selection is not None and not self.honours_selection:
            raise ValueError(
                f"{type(self).__name__} works out its own atoms, so "
                f"`selection` would have no effect. Accepting it would let a "
                f"measurement look as though it had been restricted when it "
                f"had not."
            )
        #: What the analysis worked out while running, as opposed to what it
        #: was told. Recorded beside the options and kept out of them, because
        #: a report lists the options.
        self.findings: dict[str, Any] = {}
        self.selection: str | None = (
            selection if selection is not None else self.default_selection
        )
        self.output_dir: Path = (
            Path(output_dir) if output_dir is not None else Path.cwd()
        ) / self.name
        self.options: dict[str, Any] = dict(options)

        # User-overridable plot customizations
        self._user_title: str | None = title
        self._user_xlabel: str | None = xlabel
        self._user_ylabel: str | None = ylabel
        self._user_figsize: tuple[float, float] | None = figsize
        self._user_xunit: str | None = xunit

        self.result: Any = None

    # ------------------------------------------------------------------
    # Subclass hooks (mandatory)
    # ------------------------------------------------------------------
    @abstractmethod
    def compute(self, traj: md.Trajectory) -> Any:
        """Compute the analysis. Must be deterministic and side-effect-free.

        Parameters
        ----------
        traj : mdtraj.Trajectory
            The trajectory to analyze, already sliced/strided as the user
            requested. The analysis should respect ``self.selection`` if
            relevant.

        Returns
        -------
        Any
            The analysis result. Most commonly a 1-D or 2-D NumPy array or
            a pandas DataFrame. The same object is later passed to
            :meth:`plot` and :meth:`save_data`.
        """

    @abstractmethod
    def plot(self, result: Any, ax: plt.Axes) -> None:
        """Render ``result`` onto ``ax``. Do not call plt.show()."""

    # ------------------------------------------------------------------
    # Subclass hooks (optional)
    # ------------------------------------------------------------------
    def save_data(self, result: Any, path: Path) -> Path:
        """Write the computed result to ``path``.

        Default behaviour:
          - 1-D / 2-D numpy arrays: ``np.savetxt`` (whitespace-delimited).
          - pandas DataFrames: ``to_csv`` (comma-delimited).
          - Anything else: subclass must override.

        Returns the path actually written.
        """
        path.parent.mkdir(parents=True, exist_ok=True)
        if isinstance(result, np.ndarray):
            # %.8e preserves ~8 significant figures — well beyond what
            # MD trajectories ever resolve, and round-trips cleanly through
            # np.loadtxt.
            np.savetxt(path, result, fmt="%.8e")
        elif isinstance(result, pd.DataFrame):
            result.to_csv(path, index=False)
        else:
            raise NotImplementedError(
                f"{type(self).__name__}.save_data does not know how to "
                f"serialize {type(result).__name__}. Override save_data() "
                f"in the subclass."
            )
        return path

    def figure_title(self) -> str:
        """Title shown above the figure.

        If the user passed ``title=`` at construction, that wins. Otherwise
        the subclass's :attr:`description` is used; failing that, the
        analysis name in uppercase.
        """
        if self._user_title is not None:
            return self._user_title
        return self.description or self.name.upper()

    def default_xlabel(self) -> str | None:
        """X-axis label when the user has not overridden it.

        Override in subclasses to set a domain-specific default. Returning
        ``None`` means "leave whatever the plot() method set". User-supplied
        ``xlabel=`` at construction always wins regardless.
        """
        return None

    def default_ylabel(self) -> str | None:
        """Y-axis label when the user has not overridden it. See :meth:`default_xlabel`."""
        return None

    # ------------------------------------------------------------------
    # Public entry points
    # ------------------------------------------------------------------
    def run(self, traj: md.Trajectory) -> AnalysisResult:
        """Compute, plot, and save in one call.

        This is the orchestrator's standard entry point. Returns an
        :class:`AnalysisResult` regardless of success — check
        ``result.status`` for ``"ok"`` vs ``"error"``.
        """
        started = datetime.now(timezone.utc).isoformat()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        options_path = self._write_options_manifest()

        try:
            self.result = self.compute(traj)
            self._record_what_the_mean_is_worth(traj)
            # Written again, because an analysis can only record some things
            # once it has looked at the trajectory: how confidently it knew
            # the ligand's chemistry, which measurements it could not make and
            # why, what binding modes it found. Writing the manifest before
            # `compute` and never again meant that everything an analysis
            # learned about its own run was thrown away.
            options_path = self._write_options_manifest()
            data_path = self.save_data(self.result, self.output_dir / f"{self.name}.dat")
            figure_path = self._do_plot()
            svg_path = figure_path.with_suffix(".svg")
            figure_artifacts = [figure_path]
            if svg_path.is_file():
                figure_artifacts.append(svg_path)
            finished = datetime.now(timezone.utc).isoformat()
            return AnalysisResult(
                name=self.name,
                status="ok",
                data=self.result,
                output_dir=self.output_dir,
                figure_path=figure_path,
                data_path=data_path,
                options_path=options_path,
                artifacts=[data_path, *figure_artifacts, options_path],
                message=f"{self.name}: ok",
                started_at=started,
                finished_at=finished,
            )
        except Exception as exc:  # noqa: BLE001 -- captured, reported, propagated via status
            finished = datetime.now(timezone.utc).isoformat()
            logger.exception("Analysis '%s' failed", self.name)
            return AnalysisResult(
                name=self.name,
                status="error",
                output_dir=self.output_dir,
                options_path=options_path,
                message=f"{self.name}: {exc}",
                started_at=started,
                finished_at=finished,
            )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _record_what_the_mean_is_worth(self, traj: md.Trajectory) -> None:
        """Put the mean, its error and what stands behind it in the findings.

        Recorded rather than enforced: an analysis whose series does not
        support a mean still produces its figure and its data, and the reason
        travels with them. Withholding an RMSD plot because the run was short
        would hide the evidence that it was short.
        """
        if not self.time_series:
            return

        import numpy as np

        from fastmdxplora.statistics import summarise

        values = self.result
        if hasattr(values, "columns"):
            # A frame of (frame, value): the frame number is the axis, not a
            # quantity, and averaging it would give the middle of the run.
            columns = [c for c in values.columns if str(c).lower() != "frame"]
            if len(columns) != 1:
                return
            values = values[columns[0]]
        if hasattr(values, "to_numpy"):
            values = values.to_numpy()
        try:
            series = np.asarray(values, dtype=float).squeeze()
        except (TypeError, ValueError):
            return
        if series.ndim != 1 or series.size != traj.n_frames:
            return

        settled, reason = summarise(series)
        record: dict[str, Any] = {}
        if settled is not None:
            record.update(settled.as_record())
        if reason is not None:
            record["not_a_measurement"] = reason
        self.findings["mean"] = record

    def select_atoms(self, traj: md.Trajectory) -> np.ndarray:
        """Resolve :attr:`selection` to atom indices on a given trajectory.

        Returns the full atom index array when ``selection`` is ``None``.
        Raises ``ValueError`` if the selection matches zero atoms.
        """
        if self.selection is None:
            return np.arange(traj.n_atoms)
        idx = traj.topology.select(self.selection)
        if len(idx) == 0:
            raise ValueError(
                f"Atom selection {self.selection!r} matched zero atoms in "
                f"this trajectory."
            )
        return idx

    def frame_axis(self, traj: md.Trajectory) -> tuple[np.ndarray, str]:
        """Return ``(x_values, x_label)`` for a time-series plot.

        Resolution order for the unit:
          1. User-supplied ``xunit`` at construction (``"ns"``, ``"ps"``, or ``"frames"``).
          2. ``"ns"`` if the trajectory carries usable timing information
             (``traj.time`` or ``traj.timestep``).
          3. ``"frames"`` as the always-safe fallback.

        Parameters
        ----------
        traj : mdtraj.Trajectory

        Returns
        -------
        x : np.ndarray, shape (n_frames,)
            The numerical values for the x axis.
        label : str
            A pre-formatted axis label, e.g. ``"Time (ns)"`` or ``"Frame"``.

        Notes
        -----
        MDTraj stores ``traj.time`` in picoseconds. When the trajectory
        lacks timing (e.g., loaded from a PDB without a timestep), the
        method falls back to frame indices.
        """
        unit = self._user_xunit
        if unit is not None:
            unit = unit.lower()
            if unit not in {"ns", "ps", "frames"}:
                raise ValueError(
                    f"xunit must be one of 'ns', 'ps', or 'frames'; got {unit!r}"
                )

        # If user didn't specify, prefer ns when timing data is available.
        # MDTraj sets time=0 by default for trajectories without timestamps,
        # so check whether time actually varies.
        time_ps = np.asarray(traj.time, dtype=float)
        has_real_time = time_ps is not None and len(time_ps) > 1 and not np.allclose(
            time_ps, time_ps[0]
        )

        if unit is None:
            unit = "ns" if has_real_time else "frames"

        if unit == "frames" or not has_real_time:
            return np.arange(traj.n_frames), "Frame"
        if unit == "ps":
            return time_ps, "Time (ps)"
        # ns
        return time_ps / 1000.0, "Time (ns)"

    def _do_plot(self) -> Path:
        """Internal: create the figure, call self.plot, apply user overrides, save."""
        fig, ax = new_figure(title=self.figure_title(), figsize=self._user_figsize)
        self.plot(self.result, ax)

        # Apply user/default label overrides AFTER the subclass plot() runs
        # so users can replace anything the analysis set internally.
        xlabel = self._user_xlabel if self._user_xlabel is not None else self.default_xlabel()
        ylabel = self._user_ylabel if self._user_ylabel is not None else self.default_ylabel()
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)

        return save_figure(fig, self.output_dir / f"{self.name}.png")

    def _write_options_manifest(self) -> Path:
        """Record what this analysis was asked to do, and what it found out.

        The two are kept apart. ``options`` is what somebody set, and a report
        listing them is a record of the run. ``findings`` is what the analysis
        worked out along the way -- how confidently it knew the ligand's
        chemistry, which binding modes it saw, whether a rate was supportable.
        Both belong in the manifest and only the first belongs in a document:
        putting the findings in ``options`` filled two pages of a report with
        raw tuples of atom indices.
        """
        manifest = {
            "analysis": self.name,
            "class": type(self).__name__,
            "selection": self.selection,
            "options": self.options,
            "findings": self.findings,
        }
        path = self.output_dir / "options.json"
        with path.open("w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2)
        return path
