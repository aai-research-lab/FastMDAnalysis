"""The umbrella result, drawn and recorded like every other analysis.

FastMDXplora's claim is an end-to-end study, and for an umbrella run the free
energy along the coordinate is the study. It was written to `pmf.json` and
stopped there: no figure, no entry in the analysis manifest, no mention in the
report. Sixteen analyses of the trajectory each produced a curve and a plot,
and the one result the run existed for did not.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from fastmdxplora.analysis.pmf import PMF


def _written(directory: Path, *, unsampled: int = 0) -> Path:
    """A `pmf.json` of the shape the umbrella phase writes."""
    coordinate = list(np.linspace(0.55, 0.95, 40))
    energy = [float(30.0 * (x - 0.73) ** 2 * 100) for x in coordinate]
    for index in range(unsampled):
        energy[10 + index] = None
    (directory / "pmf.json").write_text(json.dumps({
        "pmf": {"coordinate": coordinate, "free_energy_kjmol": energy},
        "overlaps": [
            {"between": [0, 1], "centres": [0.6, 0.7], "overlap": 0.43},
            {"between": [1, 2], "centres": [0.7, 0.8], "overlap": 0.54},
        ],
        "refused": None,
        "unsampled_bins": unsampled,
        "n_windows": 5,
    }), encoding="utf-8")
    return directory / "pmf.json"


class TestItReadsWhatTheStudyComputed:
    def test_it_finds_the_result_beside_the_run(self, tmp_path) -> None:
        """`pmf.json` sits beside the windows rather than inside any one of
        them, because it is the answer from all of them together."""
        _written(tmp_path)
        analysis = PMF(output_dir=tmp_path / "analysis")
        result = analysis.compute(None)

        assert len(result["coordinate"]) == 40
        assert np.isfinite(result["free_energy_kjmol"]).all()

    def test_it_does_not_recompute_the_stitching(self) -> None:
        """Stitching windows is delicate -- they have to overlap, and where
        they do not the gap is reported rather than bridged. Doing it twice
        would invite two answers to one question."""
        import inspect

        source = inspect.getsource(PMF)
        assert "does not recompute" in source
        assert "WHAM" not in source

    def test_a_run_without_an_umbrella_study_says_so(self, tmp_path) -> None:
        analysis = PMF(output_dir=tmp_path / "analysis")
        with pytest.raises(FileNotFoundError) as caught:
            analysis.compute(None)
        assert "umbrella" in str(caught.value)


class TestAnUnsampledBinStaysAGap:
    """`null` turned into zero would read as a minimum where there is no
    data -- the same defect the umbrella phase had when it wrote 1724 kJ/mol
    into bins nobody visited."""

    def test_it_is_not_a_zero(self, tmp_path) -> None:
        _written(tmp_path, unsampled=3)
        result = PMF(output_dir=tmp_path / "analysis").compute(None)
        energy = result["free_energy_kjmol"]

        assert np.isnan(energy).sum() == 3
        assert not (energy == 0).any()

    def test_the_table_says_nan_rather_than_nothing(self, tmp_path) -> None:
        """A blank field makes a row two columns in some places and one in
        others, and every reader of a whitespace table then disagrees."""
        _written(tmp_path, unsampled=2)
        analysis = PMF(output_dir=tmp_path / "analysis")
        result = analysis.compute(None)
        target = tmp_path / "pmf.dat"
        analysis.save_data(result, target)

        rows = [line.split() for line in
                target.read_text(encoding="utf-8").splitlines()
                if not line.startswith("#")]
        assert all(len(row) == 2 for row in rows)
        assert sum(1 for row in rows if row[1] == "nan") == 2


class TestItIsAnAnalysisLikeAnyOther:
    def test_it_is_registered(self) -> None:
        import fastmdxplora.analysis  # noqa: F401
        from fastmdxplora.analysis.orchestrator import available_analyses

        assert "pmf" in available_analyses()

    def test_it_runs_only_where_a_study_produced_one(self) -> None:
        """An analysis that failed on every unbiased trajectory would turn a
        missing study into a failed phase."""
        assert PMF.requires_umbrella is True

        import inspect

        from fastmdxplora.analysis import orchestrator

        source = inspect.getsource(orchestrator)
        assert "_umbrella_ok" in source
        assert 'requires_umbrella' in source

    def test_it_does_not_take_a_selection(self) -> None:
        """The coordinate was chosen when the windows were planned; a
        selection here would describe a different question."""
        assert PMF.honours_selection is False

    def test_it_is_not_a_time_series(self) -> None:
        """The x axis is the collective variable, not time, so a mean over
        it would be a mean over positions rather than over a run."""
        assert PMF.time_series is False

    def test_it_draws_the_minimum_and_the_windows(self) -> None:
        import inspect

        source = inspect.getsource(PMF.plot)
        assert "nanargmin" in source
        # The window centres, so a feature sitting on a seam can be judged.
        assert "centres" in source


class TestTheStudyDrawsItsOwnFreeEnergy:
    """Putting the drawing in the per-run analysis phase left it undrawn.
    Each window analysed itself before the last one finished, so no window
    found a result to draw, and by the time `pmf.json` existed no analysis
    phase remained. The windows had each produced ten figures of their own
    restrained trajectory, and the one result the study existed for had
    none."""

    def test_the_study_draws_it_after_writing_it(self) -> None:
        import inspect

        from fastmdxplora.batch import explorer

        source = inspect.getsource(explorer.BatchExplorer._maybe_build_pmf)
        assert "_draw_pmf(destination)" in source
        # After the numbers are on disk, so a drawing failure cannot lose them.
        assert source.index("destination.write_text") < source.index("_draw_pmf")

    def test_drawing_cannot_fail_the_study(self) -> None:
        """A study whose windows all succeeded must not fail at the last
        step, for the same reason the comparison report beside it cannot."""
        import inspect

        from fastmdxplora.batch import explorer

        source = inspect.getsource(explorer.BatchExplorer._draw_pmf)
        assert "except Exception" in source
        assert "The numbers are unaffected" in source

    def test_it_produces_the_same_three_files_as_any_analysis(
        self, tmp_path
    ) -> None:
        _written(tmp_path)
        analysis = PMF(output_dir=tmp_path / "free_energy")
        result = analysis.run(None)

        assert result.status == "ok"
        produced = {p.suffix for p in (tmp_path / "free_energy").rglob("*")
                    if p.is_file()}
        assert {".dat", ".png", ".svg"} <= produced

    def test_it_runs_without_a_trajectory(self, tmp_path) -> None:
        """The trajectory of any one window is not the study: it is a system
        held at one position by a spring."""
        _written(tmp_path)
        result = PMF(output_dir=tmp_path / "free_energy").run(None)
        assert result.status == "ok"
