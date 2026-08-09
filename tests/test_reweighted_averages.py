"""Reporting a metadynamics run's analyses as corrected averages.

The machinery that computes weights is tested in `test_reweight.py` against a
two-state system with a known answer. What is tested here is the wiring: that
the weights reach the analyses that have a meaningful average, that they stay
away from the ones that do not, and that a run nobody biased is left alone.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fastmdxplora.analysis.reweight import KB_KJ_PER_MOL_K
from fastmdxplora.analysis.reweighted_averages import (
    frame_series,
    read_colvar,
    reweight_results,
    weights_for_run,
)


# ---------------------------------------------------------------------------
# A metadynamics run on disk, small enough to reason about
# ---------------------------------------------------------------------------
def _write_run(directory: Path, *, n_hills: int = 200, settled: bool = True,
               temperature_K: float = 300.0) -> Path:
    """A run directory holding what PLUMED and the simulation phase leave.

    The collective variable oscillates, so hills land across its range and
    the bias each frame feels genuinely differs between frames -- a constant
    coordinate would make every weight equal and test nothing.
    """
    simulation = directory / "simulation"
    simulation.mkdir(parents=True, exist_ok=True)

    times = np.arange(1.0, n_hills + 1.0)
    centres = np.sin(times / 7.0)
    sigma, height = 0.2, 1.2

    simulation.joinpath("HILLS").write_text(
        "#! FIELDS time cv sigma_cv height biasf\n"
        + "\n".join(
            f"{t:.4f} {c:.6f} {sigma:.4f} {height:.4f} 10.0"
            for t, c in zip(times, centres))
        + "\n", encoding="utf-8")

    # COLVAR is written more often than frames are, as PACE is finer than the
    # reporting interval in any real run.
    colvar_times = np.arange(0.5, n_hills + 1.0, 0.5)
    colvar_cv = np.sin(colvar_times / 7.0)
    bias = _bias_from_hills(times, centres, sigma, height,
                            colvar_times, colvar_cv)
    simulation.joinpath("COLVAR").write_text(
        "#! FIELDS time cv metad.bias\n"
        + "\n".join(
            f"{t:.4f} {c:.6f} {b:.6f}"
            for t, c, b in zip(colvar_times, colvar_cv, bias))
        + "\n", encoding="utf-8")

    simulation.joinpath("metadynamics_surface.json").write_text(
        json.dumps({"provisional": not settled}), encoding="utf-8")
    simulation.joinpath("simulation_parameters.json").write_text(
        json.dumps({"temperature_K": temperature_K}), encoding="utf-8")

    analysis = directory / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    return analysis


def _bias_from_hills(hill_times, centres, sigma, height, at_times, at_values):
    """The bias felt, computed independently of the code under test."""
    out = np.zeros(at_times.size)
    for index, (t, value) in enumerate(zip(at_times, at_values)):
        laid = hill_times <= t
        out[index] = np.sum(
            height * np.exp(-((value - centres[laid]) ** 2) / (2 * sigma ** 2)))
    return out


class _Result:
    """Stands in for an AnalysisResult."""

    def __init__(self, data, status="ok"):
        self.data = data
        self.status = status


class _Cls:
    def __init__(self, reweightable):
        self.reweightable = reweightable


# ---------------------------------------------------------------------------
class TestItReadsWhatPlumedWrote:
    def test_colvar_columns_come_back_by_name(self, tmp_path: Path) -> None:
        """A steered COLVAR's third column is work and a metadynamics one's
        is bias. Reading by position would silently swap them."""
        path = tmp_path / "COLVAR"
        path.write_text(
            "#! FIELDS time cv pull.work pull.bias\n"
            "0.0 1.0 2.0 3.0\n1.0 1.5 2.5 3.5\n", encoding="utf-8")
        colvar = read_colvar(path)
        assert set(colvar) == {"time", "cv", "pull.work", "pull.bias"}
        assert colvar["pull.work"].tolist() == [2.0, 2.5]

    def test_the_reconstructed_bias_matches_plumed(self, tmp_path: Path) -> None:
        """The one independent check available: PLUMED computed the same
        quantity as the run went, and the two must agree."""
        analysis = _write_run(tmp_path)
        _, provenance = weights_for_run(analysis, np.arange(1.0, 200.0, 4.0))
        assert provenance["applies"] is True
        assert provenance["largest_disagreement_with_plumed"] < 1e-6


class TestAnUnbiasedRunIsLeftAlone:
    def test_no_hills_means_no_reweighting(self, tmp_path: Path) -> None:
        """Plain MD is already a Boltzmann ensemble. Producing a 'corrected'
        table for it would imply its averages needed correcting."""
        analysis = tmp_path / "analysis"
        analysis.mkdir(parents=True)
        weights, provenance = weights_for_run(analysis, np.arange(10.0))
        assert weights is None
        assert provenance["applies"] is False
        assert "no bias was deposited" in provenance["reason"]

    def test_the_pass_produces_nothing(self, tmp_path: Path) -> None:
        analysis = tmp_path / "analysis"
        analysis.mkdir(parents=True)
        record = reweight_results(
            {"rmsd": _Result(np.linspace(0.1, 0.3, 20))},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=20, frame_times_ps=np.arange(20.0), output_dir=analysis)
        assert record is None
        assert not (analysis / "reweighted").exists()


class TestOnlyPerFrameScalarsAreReweighted:
    """A per-frame weight lines up with a per-frame value and nothing else."""

    def test_a_bare_per_frame_array_is_taken(self) -> None:
        series = frame_series(_Result(np.arange(10.0)), (None, "x"), 10)
        assert series is not None and series.size == 10

    def test_a_breakdown_table_uses_its_total(self) -> None:
        """Rg returns the total in column zero and components beside it."""
        data = np.column_stack([np.arange(10.0), np.ones(10), np.ones(10)])
        series = frame_series(_Result(data), (None, "Rg"), 10)
        assert series is not None
        assert series.tolist() == list(np.arange(10.0))

    def test_a_named_column_is_taken(self) -> None:
        frame = pd.DataFrame({"frame": np.arange(6), "n_hbonds": np.arange(6)})
        series = frame_series(_Result(frame), ("n_hbonds", "H-bonds"), 6)
        assert series is not None and series.size == 6

    def test_a_per_residue_table_is_refused(self) -> None:
        """SASA per residue is a row per residue per frame. Weighting it by
        position would produce a number that is not a mean of anything."""
        frame = pd.DataFrame({
            "frame": np.repeat(np.arange(4), 3),
            "residue": np.tile(np.arange(3), 4),
            "sasa_nm2": np.linspace(0.0, 1.0, 12),
        })
        assert frame_series(_Result(frame), ("sasa_nm2", "SASA"), 4) is None

    def test_an_analysis_that_did_not_declare_one_is_skipped(
            self, tmp_path: Path) -> None:
        analysis = _write_run(tmp_path)
        record = reweight_results(
            {"cluster": _Result(np.arange(50.0)),
             "rmsd": _Result(np.linspace(0.1, 0.4, 50))},
            {"cluster": _Cls(None), "rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=50, frame_times_ps=np.linspace(1.0, 200.0, 50),
            output_dir=analysis)
        assert record is not None
        assert [q["analysis"] for q in record["quantities"]] == ["rmsd"]

    def test_a_failed_analysis_contributes_nothing(self, tmp_path: Path) -> None:
        analysis = _write_run(tmp_path)
        record = reweight_results(
            {"rmsd": _Result(np.linspace(0.1, 0.4, 50), status="error")},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=50, frame_times_ps=np.linspace(1.0, 200.0, 50),
            output_dir=analysis)
        assert record is None


class TestTheCorrectionIsActuallyApplied:
    def test_the_reweighted_mean_differs_from_the_raw_one(
            self, tmp_path: Path) -> None:
        """If it did not, nothing here would be doing anything."""
        analysis = _write_run(tmp_path)
        times = np.linspace(1.0, 200.0, 100)
        # A quantity that tracks the coordinate, so the bias distorts it.
        values = np.sin(times / 7.0) + 2.0
        record = reweight_results(
            {"rmsd": _Result(values)}, {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=100, frame_times_ps=times, output_dir=analysis)

        item = record["quantities"][0]
        assert item["raw_mean"] == pytest.approx(float(np.mean(values)))
        assert item["reweighted_mean"] != pytest.approx(item["raw_mean"])
        assert np.isfinite(item["shift_percent"])

    def test_the_files_land_where_the_manifest_points(
            self, tmp_path: Path) -> None:
        analysis = _write_run(tmp_path)
        times = np.linspace(1.0, 200.0, 60)
        reweight_results(
            {"rg": _Result(np.sin(times / 7.0) + 1.5)},
            {"rg": _Cls((None, "Radius of gyration (nm)"))},
            n_frames=60, frame_times_ps=times, output_dir=analysis)

        directory = analysis / "reweighted"
        assert (directory / "reweighted_averages.json").is_file()
        assert (directory / "reweighted_averages.dat").is_file()
        assert (directory / "reweighted_averages.png").is_file()


class TestItSaysWhatTheAnswerRestsOn:
    def test_the_effective_sample_size_is_recorded(self, tmp_path: Path) -> None:
        analysis = _write_run(tmp_path)
        times = np.linspace(1.0, 200.0, 80)
        record = reweight_results(
            {"rmsd": _Result(np.sin(times / 7.0) + 2.0)},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=80, frame_times_ps=times, output_dir=analysis)
        assert 0.0 < record["effective_sample_size"] <= 80.0
        assert record["usable_fraction"] <= 1.0

    def test_an_unsettled_bias_carries_its_warning(self, tmp_path: Path) -> None:
        """The simulation phase already judged the surface provisional.
        Judging it a second time here could disagree with the report."""
        analysis = _write_run(tmp_path, settled=False)
        times = np.linspace(1.0, 200.0, 60)
        record = reweight_results(
            {"rmsd": _Result(np.sin(times / 7.0) + 2.0)},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=60, frame_times_ps=times, output_dir=analysis)
        assert record["settled"] is False
        assert any("approximate" in w for w in record["warnings"])

    def test_a_missing_temperature_is_said_rather_than_assumed_silently(
            self, tmp_path: Path) -> None:
        analysis = _write_run(tmp_path)
        (analysis.parent / "simulation" / "simulation_parameters.json").unlink()
        times = np.linspace(1.0, 200.0, 60)
        record = reweight_results(
            {"rmsd": _Result(np.sin(times / 7.0) + 2.0)},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=60, frame_times_ps=times, output_dir=analysis)
        assert any("300 K" in w for w in record["warnings"])

    def test_the_temperature_that_was_recorded_is_the_one_used(
            self, tmp_path: Path) -> None:
        """Weights go as exp(V/RT); the wrong T is wrong by an exponent."""
        analysis = _write_run(tmp_path, temperature_K=350.0)
        _, provenance = weights_for_run(analysis, np.linspace(1.0, 200.0, 40))
        assert provenance["temperature_K"] == pytest.approx(350.0)
        assert provenance["temperature_from_record"] is True


class TestTheWeightsMeasureTheCoordinateAndNotTheClock:
    """The bias grows as hills accumulate, so exp(V/RT) grows with time
    wherever the system is. Without the c(t) offset the weights rank frames
    by when they were written: on a converged well-tempered run below, the
    last fifth of the frames carried the whole weight and a mean over five
    hundred rested on seven of them.
    """

    @staticmethod
    def _well_tempered(directory: Path, n: int = 3000):
        """A run whose surface fills and whose coordinate diffuses."""
        rng = np.random.default_rng(0)
        times = np.arange(1.0, n + 1.0)
        cv = np.clip(np.cumsum(rng.normal(0, 0.05, n)), -1.5, 1.5)
        heights = 1.2 * np.exp(-times / 1200.0)  # well-tempered decay
        sigma = 0.2

        simulation = directory / "simulation"
        simulation.mkdir(parents=True, exist_ok=True)
        simulation.joinpath("HILLS").write_text(
            "#! FIELDS time cv sigma_cv height biasf\n"
            + "\n".join(f"{t:.4f} {c:.6f} {sigma} {h:.6f} 10.0"
                        for t, c, h in zip(times, cv, heights)) + "\n",
            encoding="utf-8")
        colvar_times = np.arange(0.5, n + 1.0, 0.5)
        simulation.joinpath("COLVAR").write_text(
            "#! FIELDS time cv metad.bias\n"
            + "\n".join(f"{t:.4f} {c:.6f} 0.0" for t, c in
                        zip(colvar_times, np.interp(colvar_times, times, cv)))
            + "\n", encoding="utf-8")
        simulation.joinpath("metadynamics_surface.json").write_text(
            json.dumps({"provisional": False}), encoding="utf-8")
        simulation.joinpath("simulation_parameters.json").write_text(
            json.dumps({"temperature_K": 300.0}), encoding="utf-8")
        analysis = directory / "analysis"
        analysis.mkdir(parents=True, exist_ok=True)
        return analysis, times, n

    def test_the_weight_does_not_pile_up_at_the_end_of_the_run(
            self, tmp_path: Path) -> None:
        analysis, _, n = self._well_tempered(tmp_path)
        frames = np.linspace(1.0, n, 500)
        weights, _ = weights_for_run(analysis, frames)

        values = weights.values / np.sum(weights.values)
        assert np.sum(values[400:]) < 0.75, (
            "the last fifth of the frames carries almost all the weight, "
            "which is the c(t) failure")
        assert weights.effective_sample_size > 50.0

    def test_the_offset_grows_with_the_surface(self, tmp_path: Path) -> None:
        """c(t) tracks the filling. It should rise over the run and it is
        what cancels the growth of the bias."""
        from fastmdxplora.analysis.reweighted_averages import c_of_t
        from fastmdxplora.simulation.metad_surface import read_hills

        analysis, times, _ = self._well_tempered(tmp_path)
        hills = read_hills(analysis.parent / "simulation" / "HILLS")
        offset = c_of_t(hills, times, 300.0)
        assert offset[0] < offset[-1]
        assert np.all(np.isfinite(offset))

    def test_an_untempered_run_still_gets_an_offset(self, tmp_path: Path) -> None:
        """A bias factor of one would divide by zero in the well-tempered
        form, so that case takes the limit instead of failing."""
        from fastmdxplora.analysis.reweighted_averages import c_of_t
        from fastmdxplora.simulation.metad_surface import read_hills

        analysis = _write_run(tmp_path)
        path = analysis.parent / "simulation" / "HILLS"
        path.write_text(path.read_text().replace(" 10.0", " 1.0"),
                        encoding="utf-8")
        offset = c_of_t(read_hills(path), np.arange(1.0, 200.0), 300.0)
        assert np.all(np.isfinite(offset))
        assert offset[-1] > offset[0]


class TestTheDeclarationsAreWiredUp:
    """The pass finds analyses by a class attribute, so a typo in one of them
    is a quantity that silently stops being corrected."""

    @pytest.mark.parametrize("module,name,column", [
        ("rmsd", "RMSD", None),
        ("rg", "RadiusOfGyration", None),
        ("hbonds", "HBonds", "n_hbonds"),
        ("sasa", "SASA", "sasa_nm2"),
    ])
    def test_the_per_frame_analyses_declare_a_scalar(
            self, module: str, name: str, column: str | None) -> None:
        import importlib
        imported = importlib.import_module(f"fastmdxplora.analysis.{module}")
        cls = next(
            value for key, value in vars(imported).items()
            if isinstance(value, type)
            and getattr(value, "time_series", False)
            and getattr(value, "name", "") == module)
        assert cls.reweightable is not None, f"{module} declares no scalar"
        assert cls.reweightable[0] == column

    def test_an_analysis_without_one_defaults_to_none(self) -> None:
        from fastmdxplora.analysis.cluster import Cluster
        assert getattr(Cluster, "reweightable", None) is None


class _Populations:
    reweightable = None
    reweightable_populations = True


class TestPopulationsAreCounts:
    """A population is the mean of an indicator, so it reweights exactly as
    any other average does -- and it is the quantity a bias distorts most,
    since escaping a well is what the bias is for."""

    def test_a_known_population_is_recovered(self) -> None:
        """The two-state case with the answer known independently, asked at
        the level a clustering reports it."""
        from fastmdxplora.analysis.reweight import (
            KB_KJ_PER_MOL_K, weights_from_bias)
        from fastmdxplora.analysis.reweighted_averages import populations

        rng = np.random.default_rng(1)
        n = 4000
        # What metadynamics samples once the bias has levelled the states.
        in_b = rng.random(n) < 0.5
        labels = in_b.astype(int)
        kT = KB_KJ_PER_MOL_K * 300.0
        bias = np.where(labels == 0, kT * np.log(0.8 / 0.2), 0.0)
        weights = weights_from_bias(bias, temperature_K=300.0)

        rows = {row["label"]: row for row in populations(labels, weights)}
        assert rows[0]["raw_fraction"] == pytest.approx(0.5, abs=0.03)
        assert rows[0]["reweighted_fraction"] == pytest.approx(0.8, abs=0.03)
        assert rows[1]["reweighted_fraction"] == pytest.approx(0.2, abs=0.03)

    def test_the_populations_still_sum_to_one(self) -> None:
        from fastmdxplora.analysis.reweight import weights_from_bias
        from fastmdxplora.analysis.reweighted_averages import populations

        rng = np.random.default_rng(2)
        labels = rng.integers(0, 4, 500)
        weights = weights_from_bias(rng.normal(0, 5, 500), temperature_K=300.0)
        rows = populations(labels, weights)
        assert sum(r["reweighted_fraction"] for r in rows) == pytest.approx(1.0)

    def test_an_unbiased_run_leaves_them_alone(self) -> None:
        from fastmdxplora.analysis.reweight import weights_from_bias
        from fastmdxplora.analysis.reweighted_averages import populations

        labels = np.array([0, 0, 0, 1, 1, 1, 1, 1, 1, 1])
        weights = weights_from_bias(np.zeros(10), temperature_K=300.0)
        rows = {r["label"]: r for r in populations(labels, weights)}
        assert rows[0]["reweighted_fraction"] == pytest.approx(0.3)
        assert rows[1]["reweighted_fraction"] == pytest.approx(0.7)


class TestOnlyLabelsBecomePopulations:
    def test_a_mapping_of_methods_is_read(self) -> None:
        from fastmdxplora.analysis.reweighted_averages import frame_labels

        labels = frame_labels(
            _Result({"kmeans": np.zeros(8, int),
                     "hierarchical": np.ones(8, int)}), 8)
        assert set(labels) == {"kmeans", "hierarchical"}

    def test_a_float_series_is_not_a_set_of_states(self) -> None:
        """RMSD per frame is a measurement. Treated as labels it would make
        every frame its own state and produce a row for each."""
        from fastmdxplora.analysis.reweighted_averages import frame_labels

        assert frame_labels(_Result(np.linspace(0.0, 1.0, 20)), 20) == {}

    def test_labels_of_the_wrong_length_are_refused(self) -> None:
        from fastmdxplora.analysis.reweighted_averages import frame_labels

        assert frame_labels(_Result({"kmeans": np.zeros(5, int)}), 20) == {}

    def test_a_single_state_is_not_a_population(self, tmp_path: Path) -> None:
        """Nothing was distinguished, and a table saying one state was
        visited all of the time before and after tells a reader nothing."""
        analysis = _write_run(tmp_path)
        times = np.linspace(1.0, 200.0, 40)
        record = reweight_results(
            {"cluster": _Result({"kmeans": np.zeros(40, int)})},
            {"cluster": _Populations()},
            n_frames=40, frame_times_ps=times, output_dir=analysis)
        assert record is None

    def test_populations_alone_still_produce_a_record(
            self, tmp_path: Path) -> None:
        """A run whose only reweightable result is a clustering has no means
        to plot, and an empty figure would read as one that failed."""
        analysis = _write_run(tmp_path)
        times = np.linspace(1.0, 200.0, 60)
        labels = (np.sin(times / 7.0) > 0).astype(int)
        record = reweight_results(
            {"cluster": _Result({"kmeans": labels})},
            {"cluster": _Populations()},
            n_frames=60, frame_times_ps=times, output_dir=analysis)
        assert record is not None
        assert record["quantities"] == []
        assert len(record["populations"][0]["states"]) == 2
        assert not (analysis / "reweighted" / "reweighted_averages.png").exists()

    def test_the_clustering_declares_itself(self) -> None:
        from fastmdxplora.analysis.cluster import Cluster
        assert Cluster.reweightable_populations is True

    def test_a_measurement_analysis_does_not(self) -> None:
        from fastmdxplora.analysis.rmsd import RMSD
        assert RMSD.reweightable_populations is False


class TestABiasedRunThatCannotBeCorrectedStillSaysSo:
    """Umbrella and steered runs have no weights that recover an equilibrium
    average. Writing nothing left their averages looking exactly like an
    ordinary run's -- which is worse than the metadynamics case, not better,
    because there is no corrected number to put beside them."""

    @staticmethod
    def _method_run(tmp_path: Path, script: str) -> Path:
        simulation = tmp_path / "simulation"
        simulation.mkdir(parents=True, exist_ok=True)
        (simulation / script).write_text("# PLUMED\n", encoding="utf-8")
        (simulation / "COLVAR").write_text(
            "#! FIELDS time cv\n0.0 1.0\n1.0 1.1\n", encoding="utf-8")
        analysis = tmp_path / "analysis"
        analysis.mkdir(parents=True, exist_ok=True)
        return analysis

    @pytest.mark.parametrize("script,method", [
        ("umbrella.plumed", "umbrella"),
        ("steered.plumed", "steered"),
    ])
    def test_the_method_is_recognised_and_recorded(
            self, tmp_path: Path, script: str, method: str) -> None:
        analysis = self._method_run(tmp_path, script)
        record = reweight_results(
            {"rmsd": _Result(np.linspace(0.1, 0.4, 20))},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=20, frame_times_ps=np.arange(20.0), output_dir=analysis)

        assert record is not None
        assert record["applies"] is False
        assert record["biasing_method"] == method
        assert (analysis / "reweighted" / "reweighted_averages.json").is_file()

    def test_no_table_and_no_figure_are_invented(self, tmp_path: Path) -> None:
        """There is nothing to correct, so a corrected column would be a
        column of the same numbers under a heading claiming otherwise."""
        analysis = self._method_run(tmp_path, "steered.plumed")
        record = reweight_results(
            {"rmsd": _Result(np.linspace(0.1, 0.4, 20))},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=20, frame_times_ps=np.arange(20.0), output_dir=analysis)

        assert record["quantities"] == []
        assert not (analysis / "reweighted" / "reweighted_averages.png").exists()

    def test_the_reason_is_the_one_the_methods_text_gives(
            self, tmp_path: Path) -> None:
        """Two statements of the same claim that could drift apart is worse
        than one, so these say the same thing the report already says."""
        analysis = self._method_run(tmp_path, "umbrella.plumed")
        record = reweight_results(
            {"rmsd": _Result(np.linspace(0.1, 0.4, 20))},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=20, frame_times_ps=np.arange(20.0), output_dir=analysis)
        assert "potential of mean force" in record["reason"]
        assert "not measurements of the unrestrained" in record["reason"]

    def test_an_ordinary_run_is_still_left_entirely_alone(
            self, tmp_path: Path) -> None:
        """No marker, no bias, nothing to label: a record here would imply
        the averages had needed one."""
        analysis = tmp_path / "analysis"
        analysis.mkdir(parents=True)
        record = reweight_results(
            {"rmsd": _Result(np.linspace(0.1, 0.4, 20))},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=20, frame_times_ps=np.arange(20.0), output_dir=analysis)
        assert record is None
        assert not (analysis / "reweighted").exists()

    def test_metadynamics_is_not_caught_by_this_path(
            self, tmp_path: Path) -> None:
        """It has weights, and must still take the reweighting route."""
        analysis = _write_run(tmp_path)
        (analysis.parent / "simulation" / "metadynamics.plumed").write_text(
            "# PLUMED\n", encoding="utf-8")
        times = np.linspace(1.0, 200.0, 50)
        record = reweight_results(
            {"rmsd": _Result(np.sin(times / 7.0) + 2.0)},
            {"rmsd": _Cls((None, "RMSD (nm)"))},
            n_frames=50, frame_times_ps=times, output_dir=analysis)
        assert record.get("applies", True) is not False
        assert record["quantities"]
