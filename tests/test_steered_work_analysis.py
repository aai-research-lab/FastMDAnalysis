"""The pull's work, drawn and recorded like every other analysis.

Written to `steered_work.json` and stopped there: no figure, no entry in the
analysis manifest, no mention in the report -- the same gap the umbrella
study and the metadynamics surface had, and the last of the three.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from fastmdxplora.analysis.steered_work import SteeredWork


def _record(run: Path, *, requested=1.0, reached=0.769, total=67.5,
            snap=False, samples=101) -> Path:
    """A record of the shape the simulation phase writes.

    The numbers are from a real pull: a tripeptide asked to 1.0 nm reached
    0.769 while 67.5 kJ/mol of work was done.
    """
    (run / "simulation").mkdir(parents=True, exist_ok=True)
    coordinate = list(np.linspace(0.73, reached, samples))
    if snap:
        # All the work in one step: a pull that broke past something.
        work = [0.0] * (samples // 2) + [total] * (samples - samples // 2)
    else:
        work = list(np.linspace(0.0, total, samples))
    (run / "simulation" / "steered_work.json").write_text(json.dumps({
        "work_kjmol": total,
        "pull_rate_per_ns": 1.35,
        "requested_to": requested,
        "from": coordinate[0],
        "to": coordinate[-1],
        "samples": samples,
        "trajectory": {"coordinate": coordinate, "work_kjmol": work},
    }), encoding="utf-8")
    return run


class TestItDrawsTheCurveRatherThanTheNumber:
    """A pull that accumulates work smoothly has met resistance all the way;
    one that accumulates it in a step has snapped past something. The total
    is the same in both cases."""

    def test_the_whole_trajectory_is_drawn(self, tmp_path) -> None:
        _record(tmp_path)
        analysis = SteeredWork(output_dir=tmp_path / "analysis")
        result = analysis.run(None)

        assert result.status == "ok"
        assert len(analysis.result["work_kjmol"]) == 101

    def test_a_smooth_pull_and_a_snap_differ(self, tmp_path) -> None:
        smooth = SteeredWork(output_dir=tmp_path / "smooth")
        _record(tmp_path)
        smooth.run(None)

        other = tmp_path / "snapped"
        other.mkdir()
        _record(other, snap=True)
        snapped = SteeredWork(output_dir=other / "analysis")
        snapped.run(None)

        # Same total, different shape -- which is the point of the curve.
        assert smooth.result["work_kjmol"][-1] == pytest.approx(
            snapped.result["work_kjmol"][-1])
        assert np.std(np.diff(smooth.result["work_kjmol"])) < np.std(
            np.diff(snapped.result["work_kjmol"]))

    def test_it_produces_the_same_three_files(self, tmp_path) -> None:
        _record(tmp_path)
        SteeredWork(output_dir=tmp_path / "analysis").run(None)
        produced = {p.suffix for p in (tmp_path / "analysis").rglob("*")
                    if p.is_file()}
        assert {".dat", ".png", ".svg"} <= produced


class TestItSaysWhatTheWorkIsNot:
    def test_the_table_carries_the_caveat(self, tmp_path) -> None:
        """Not only in a terminal somebody has since closed."""
        _record(tmp_path)
        SteeredWork(output_dir=tmp_path / "analysis").run(None)

        table = (tmp_path / "analysis" / "steered_work"
                 / "steered_work.dat").read_text(encoding="utf-8")
        assert "not a free energy" in table

    def test_the_figure_says_it_too(self) -> None:
        import inspect

        source = inspect.getsource(SteeredWork.plot)
        assert "a pathway, not a free energy" in source

    def test_the_module_explains_why(self) -> None:
        """Dissipated work does not cancel, and Jarzynski needs many pulls."""
        import inspect

        from fastmdxplora.analysis import steered_work

        source = inspect.getsource(steered_work)
        assert "Jarzynski" in source
        assert "how fast the anchor moved" in source


class TestTheGapBetweenAskedAndReached:
    """The anchor travels the whole way and the system lags behind it; that
    lag is the dissipation. A tripeptide asked to 1.0 nm reached 0.769."""

    def test_it_is_marked_where_they_differ(self, tmp_path) -> None:
        import inspect

        _record(tmp_path, requested=1.0, reached=0.769)
        analysis = SteeredWork(output_dir=tmp_path / "analysis")
        analysis.run(None)

        assert analysis._requested == 1.0
        source = inspect.getsource(SteeredWork.plot)
        assert "asked for" in source and "reached" in source

    def test_nothing_is_marked_where_the_pull_arrived(self, tmp_path) -> None:
        import inspect

        source = inspect.getsource(SteeredWork.plot)
        # Only drawn where they differ, so a pull that arrived is not
        # annotated with a line on top of its own endpoint.
        assert "abs(requested - reached) > 1e-6" in source


class TestItIsAnAnalysisLikeAnyOther:
    def test_it_is_registered(self) -> None:
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.orchestrator import available_analyses

        assert "steered_work" in available_analyses()

    def test_it_runs_only_where_a_pull_produced_one(self) -> None:
        assert SteeredWork.requires_steered is True

        import inspect

        from fastmdxplora.analysis import orchestrator

        assert "_steered_ok" in inspect.getsource(orchestrator)

    def test_a_run_without_a_pull_says_so(self, tmp_path) -> None:
        with pytest.raises(FileNotFoundError) as caught:
            SteeredWork(output_dir=tmp_path / "analysis").compute(None)
        assert "steered" in str(caught.value)
