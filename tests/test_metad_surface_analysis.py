"""The metadynamics surface, drawn and recorded like every other analysis.

Written to `metadynamics_surface.json` and stopped there: no figure, no entry
in the analysis manifest, no mention in the report. Ten analyses of the
trajectory each produced a curve and a plot, and the one result the run
existed to produce did not -- the same gap the umbrella study had.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from fastmdxplora.analysis.metad_surface import MetadynamicsSurface


def _record(run: Path, *, provisional=False, refused=None, empty=False) -> Path:
    (run / "simulation").mkdir(parents=True, exist_ok=True)
    grid = list(np.linspace(-np.pi, np.pi, 200))
    energy = [] if empty else [float(30.0 * (1 - np.cos(x))) for x in grid]
    (run / "simulation" / "metadynamics_surface.json").write_text(json.dumps({
        "refused": refused,
        "provisional": provisional,
        "evidence": {"drift_ceiling_kjmol": 20.0, "drift_kjmol": 2.3},
        "grid": [] if empty else grid,
        "free_energy_kjmol": energy,
    }), encoding="utf-8")
    return run


def _record_2d(run: Path) -> Path:
    (run / "simulation").mkdir(parents=True, exist_ok=True)
    first = np.linspace(-np.pi, np.pi, 8)
    second = np.linspace(-np.pi, np.pi, 6)
    energy = ((1.0 - np.cos(first[:, None]))
              + 2.0 * (1.0 - np.cos(second[None, :])))
    record = {
        "refused": None,
        "provisional": False,
        "evidence": {
            "per_dimension": [
                {"name": "cv1", "periodic": True},
                {"name": "cv2", "periodic": True},
            ],
        },
        "dimensions": 2,
        "grid": None,
        "axes": [first.tolist(), second.tolist()],
        "free_energy_kjmol": energy.tolist(),
    }
    (run / "simulation" / "metadynamics_surface.json").write_text(
        json.dumps(record), encoding="utf-8")
    return run


class TestItDrawsWhatTheRunComputed:
    def test_the_one_dimensional_result_keeps_its_existing_shape(
            self, tmp_path) -> None:
        _record(tmp_path)
        analysis = MetadynamicsSurface(output_dir=tmp_path / "analysis")

        result = analysis.compute(None)

        assert set(result) == {"coordinate", "free_energy_kjmol"}
        assert result["coordinate"].shape == result["free_energy_kjmol"].shape

    def test_a_settled_surface_is_drawn(self, tmp_path) -> None:
        _record(tmp_path)
        result = MetadynamicsSurface(output_dir=tmp_path / "analysis").run(None)

        assert result.status == "ok"
        produced = {p.suffix for p in (tmp_path / "analysis").rglob("*")
                    if p.is_file()}
        assert {".dat", ".png", ".svg"} <= produced

    def test_it_does_not_recompute_the_hills(self) -> None:
        """Summing them twice would invite two answers, and the trajectory
        of a metadynamics run is not a Boltzmann ensemble anyway."""
        import inspect

        from fastmdxplora.analysis import metad_surface

        # The module, not the class: `getsource` on a class does not reach
        # the module docstring where the reasoning lives.
        source = inspect.getsource(metad_surface)
        assert "recompute it" in source
        assert "not a Boltzmann ensemble" in source

    def test_a_run_without_metadynamics_says_so(self, tmp_path) -> None:
        analysis = MetadynamicsSurface(output_dir=tmp_path / "analysis")
        with pytest.raises(FileNotFoundError) as caught:
            analysis.compute(None)
        assert "metadynamics" in str(caught.value)


class TestATwoDimensionalSurfaceKeepsBothCoordinates:
    def test_grid_may_be_null_when_axes_describe_the_surface(
            self, tmp_path) -> None:
        _record_2d(tmp_path)
        analysis = MetadynamicsSurface(output_dir=tmp_path / "analysis")

        result = analysis.compute(None)

        assert result["dimensions"] == 2
        assert tuple(len(axis) for axis in result["axes"]) == (8, 6)
        assert result["free_energy_kjmol"].shape == (8, 6)

    def test_plotting_and_export_use_both_cv_axes(self, tmp_path) -> None:
        _record_2d(tmp_path)
        analysis = MetadynamicsSurface(output_dir=tmp_path / "analysis")

        completed = analysis.run(None)

        assert completed.status == "ok"
        assert completed.figure_path.is_file()
        lines = completed.data_path.read_text(encoding="utf-8").splitlines()
        assert lines[0] == "# cv1 cv2 free_energy_kjmol"
        values = np.loadtxt(completed.data_path)
        assert values.shape == (8 * 6, 3)
        assert len(np.unique(values[:, 0])) == 8
        assert len(np.unique(values[:, 1])) == 6


class TestAProvisionalSurfaceIsStillDrawn:
    """A run whose bias has not settled still describes the landscape it has
    filled so far, and the refusal beside it explains what is missing.
    Withholding the picture as well leaves a reader with a sentence where
    they could have had both."""

    def test_it_is_drawn(self, tmp_path) -> None:
        _record(tmp_path, provisional=True,
                refused="the surface moved 8.8 kJ/mol")
        result = MetadynamicsSurface(output_dir=tmp_path / "analysis").run(None)
        assert result.status == "ok"

    def test_the_caveat_travels_with_the_numbers(self, tmp_path) -> None:
        """Not only in a terminal somebody has since closed."""
        _record(tmp_path, provisional=True, refused="not settled")
        MetadynamicsSurface(output_dir=tmp_path / "analysis").run(None)

        table = (tmp_path / "analysis" / "metad_surface"
                 / "metad_surface.dat").read_text(encoding="utf-8")
        assert "provisional" in table.splitlines()[1]

    def test_a_settled_one_carries_no_such_note(self, tmp_path) -> None:
        _record(tmp_path, provisional=False)
        MetadynamicsSurface(output_dir=tmp_path / "analysis").run(None)

        table = (tmp_path / "analysis" / "metad_surface"
                 / "metad_surface.dat").read_text(encoding="utf-8")
        assert "provisional" not in table

    def test_a_coordinate_that_never_moved_has_nothing_to_draw(
        self, tmp_path
    ) -> None:
        """The one refusal where there is genuinely no surface."""
        _record(tmp_path, empty=True,
                refused="Every hill was deposited at the same value.")
        analysis = MetadynamicsSurface(output_dir=tmp_path / "analysis")
        with pytest.raises(ValueError) as caught:
            analysis.compute(None)
        assert "no surface" in str(caught.value)


class TestItIsAnAnalysisLikeAnyOther:
    def test_it_is_registered(self) -> None:
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.orchestrator import available_analyses

        assert "metad_surface" in available_analyses()

    def test_it_runs_only_where_a_run_produced_one(self) -> None:
        assert MetadynamicsSurface.requires_metadynamics is True

        import inspect

        from fastmdxplora.analysis import orchestrator

        assert "_metadynamics_ok" in inspect.getsource(orchestrator)

    def test_it_marks_where_the_surface_is_trustworthy(self, tmp_path) -> None:
        """Above the drift ceiling the estimate moves by several kJ/mol
        however long the run, so the shape is read rather than the points."""
        import inspect

        source = inspect.getsource(MetadynamicsSurface.plot)
        assert "drift_ceiling_kjmol" in source
        assert "axhline" in source
