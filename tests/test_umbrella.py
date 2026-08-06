"""A free energy along a coordinate, from equilibrium sampling at each point.

Steered MD drags a system and reports work that depends on how fast it was
dragged. Umbrella sampling holds it at a series of positions, lets each
equilibrate, and recombines the sampling -- so nothing is hurried and nothing
is dissipated.
"""

from __future__ import annotations

import numpy as np
import pytest

from fastmdxplora.simulation.umbrella import KB_KJ

#: A double well with a barrier of exactly 10 kJ/mol, so the recovered value
#: can be checked against a number rather than against plausibility.
def _true_free_energy(x):
    return 10.0 * (np.asarray(x) ** 2 - 1.0) ** 2


def _sample(centre, force, n=6000, seed=0, temperature_K=300.0, width=0.8):
    """Draw from the biased distribution the window would sample."""
    kT = KB_KJ * temperature_K
    rng = np.random.RandomState(seed)
    drawn = []
    while len(drawn) < n:
        x = rng.uniform(centre - width, centre + width, 400)
        energy = _true_free_energy(x) + 0.5 * force * (x - centre) ** 2
        accept = np.exp(-(energy - energy.min()) / kT)
        drawn.extend(x[rng.uniform(size=x.size) < accept])
    return np.array(drawn[:n])


class TestItRecoversAFreeEnergyItWasGiven:
    """The check that matters: sample a landscape whose barrier is known and
    see whether the number comes back."""

    def test_the_barrier_is_recovered(self) -> None:
        from fastmdxplora.simulation.umbrella import compute_pmf, plan_windows

        plan = plan_windows({
            "collective_variable": "distance", "from": -1.6, "to": 1.6,
            "n_windows": 17, "force_constant": 200.0})
        samples = {w.index: _sample(w.centre, w.force_constant, seed=w.index)
                   for w in plan.windows}

        result = compute_pmf(samples, plan, temperature_K=300.0)
        assert result["refused"] is None

        x = np.array(result["pmf"]["coordinate"])
        f = np.array(result["pmf"]["free_energy_kjmol"])
        inside = (x > -1.4) & (x < 1.4) & np.isfinite(f)

        truth = _true_free_energy(x[inside])
        truth = truth - truth.min()
        error = np.sqrt(((f[inside] - truth) ** 2).mean())
        assert error < 1.0, f"recovered PMF is off by {error:.2f} kJ/mol RMS"

    def test_the_barrier_height_is_about_right(self) -> None:
        from fastmdxplora.simulation.umbrella import compute_pmf, plan_windows

        plan = plan_windows({
            "collective_variable": "distance", "from": -1.6, "to": 1.6,
            "n_windows": 17, "force_constant": 200.0})
        result = compute_pmf(
            {w.index: _sample(w.centre, w.force_constant, seed=w.index)
             for w in plan.windows}, plan)

        x = np.array(result["pmf"]["coordinate"])
        f = np.array(result["pmf"]["free_energy_kjmol"])
        at_barrier = f[np.argmin(np.abs(x))]
        assert at_barrier == pytest.approx(10.0, abs=1.5)


class TestWindowsThatDoNotOverlap:
    """Recombination stitches histograms together. Where two neighbours never
    visit the same value there is nothing to stitch, and a curve through the
    gap is interpolation presented as a measurement.
    """

    @staticmethod
    def _separated():
        from fastmdxplora.simulation.umbrella import compute_pmf, plan_windows

        plan = plan_windows({
            "collective_variable": "distance", "from": -1.5, "to": 1.5,
            "n_windows": 4, "force_constant": 3000.0})
        samples = {w.index: _sample(w.centre, w.force_constant, n=4000,
                                    seed=w.index, width=0.5)
                   for w in plan.windows}
        return compute_pmf(samples, plan)

    def test_no_free_energy_is_returned(self) -> None:
        assert self._separated()["pmf"] is None

    def test_the_gap_is_named(self) -> None:
        refusal = self._separated()["refused"]
        assert "do not overlap" in refusal
        assert "windows 0 and 1" in refusal

    def test_and_what_would_close_it(self) -> None:
        """More windows or a softer constant. Not more sampling."""
        refusal = self._separated()["refused"]
        assert "More windows" in refusal
        assert "Sampling for longer will not" in refusal

    def test_the_overlaps_are_reported_either_way(self) -> None:
        """So somebody can see how close the windows came."""
        result = self._separated()
        assert result["overlaps"]
        assert all("overlap" in o for o in result["overlaps"])


class TestPlanningWindows:
    def test_a_range_and_a_count_becomes_windows(self) -> None:
        from fastmdxplora.simulation.umbrella import plan_windows

        plan = plan_windows({"collective_variable": "distance", "from": 0.4,
                             "to": 2.0, "n_windows": 9,
                             "force_constant": 500.0})
        assert len(plan.windows) == 9
        assert plan.windows[0].centre == pytest.approx(0.4)
        assert plan.windows[-1].centre == pytest.approx(2.0)

    def test_explicit_centres_are_taken_as_given(self) -> None:
        from fastmdxplora.simulation.umbrella import plan_windows

        plan = plan_windows({"collective_variable": "distance",
                             "centres": [0.4, 0.5, 0.7, 1.2],
                             "force_constant": 500.0})
        assert [w.centre for w in plan.windows] == [0.4, 0.5, 0.7, 1.2]

    def test_the_force_constant_has_no_default(self) -> None:
        """It sets how far a window wanders and therefore whether neighbours
        overlap. There is no value right for an arbitrary coordinate."""
        from fastmdxplora.simulation.umbrella import plan_windows

        with pytest.raises(ValueError, match="force_constant"):
            plan_windows({"collective_variable": "distance", "from": 0.0,
                          "to": 1.0, "n_windows": 5})

    def test_one_window_is_not_umbrella_sampling(self) -> None:
        from fastmdxplora.simulation.umbrella import plan_windows

        with pytest.raises(ValueError, match="at least two windows"):
            plan_windows({"collective_variable": "distance",
                          "centres": [1.0], "force_constant": 500.0})

    def test_windows_become_a_sweep(self) -> None:
        """One system at many restraint positions is the shape the batch
        machinery already runs, so it schedules them rather than a second
        one."""
        from fastmdxplora.simulation.umbrella import plan_windows, windows_as_sweep

        plan = plan_windows({"collective_variable": "distance", "from": 0.0,
                             "to": 1.0, "n_windows": 5,
                             "force_constant": 500.0})
        runs = windows_as_sweep(plan)
        assert len(runs) == 5
        assert runs[0]["id"] == "window_00"
        assert runs[-1]["umbrella_centre"] == pytest.approx(1.0)


class TestOverlapItself:
    def test_identical_sampling_overlaps_completely(self) -> None:
        from fastmdxplora.simulation.umbrella import overlap_between

        x = np.random.RandomState(0).normal(size=2000)
        assert overlap_between(x, x) == pytest.approx(1.0, abs=0.05)

    def test_separated_sampling_does_not_overlap(self) -> None:
        from fastmdxplora.simulation.umbrella import overlap_between

        rng = np.random.RandomState(0)
        assert overlap_between(rng.normal(0, 0.1, 2000),
                               rng.normal(5, 0.1, 2000)) < 0.01

    def test_empty_sampling_overlaps_nothing(self) -> None:
        from fastmdxplora.simulation.umbrella import overlap_between

        assert overlap_between(np.array([]), np.array([1.0, 2.0])) == 0.0


class TestExpandingIntoRuns:
    """An umbrella job is one system held at many positions, which is the
    shape the batch machinery already runs -- so the scheduling, parallelism
    and per-GPU pinning come from the code that already does those things."""

    @staticmethod
    def _config():
        return {
            "output": "runs/pmf",
            "systems": [{"system": "181L"}],
            "simulation": {"duration_ns": 5, "umbrella": {
                "collective_variable": "ligand_distance",
                "site_selection": "resid 84 to 121 and name CA",
                "from": 0.4, "to": 2.0, "n_windows": 9,
                "force_constant": 1000.0}},
        }

    def test_one_block_becomes_a_run_per_window(self) -> None:
        from fastmdxplora.simulation.umbrella import expand_umbrella

        expanded = expand_umbrella(self._config())
        assert len(expanded["systems"]) == 9
        assert expanded["systems"][0]["id"] == "window_00"

    def test_each_run_holds_its_own_position(self) -> None:
        from fastmdxplora.simulation.umbrella import expand_umbrella

        expanded = expand_umbrella(self._config())
        centres = [e["simulation"]["umbrella"]["centre"]
                   for e in expanded["systems"]]
        assert centres == sorted(centres)
        assert centres[0] == pytest.approx(0.4)
        assert centres[-1] == pytest.approx(2.0)

    def test_shared_settings_survive(self) -> None:
        from fastmdxplora.simulation.umbrella import expand_umbrella

        expanded = expand_umbrella(self._config())
        assert expanded["simulation"]["duration_ns"] == 5

    def test_the_block_is_consumed(self) -> None:
        """Leaving it would have every run try to expand it again."""
        from fastmdxplora.simulation.umbrella import expand_umbrella

        assert "umbrella" not in expand_umbrella(self._config())["simulation"]

    def test_a_config_without_one_is_untouched(self) -> None:
        from fastmdxplora.simulation.umbrella import expand_umbrella

        plain = {"output": "runs/x", "systems": [{"system": "1UBQ"}]}
        assert expand_umbrella(plain) == plain

    def test_several_systems_at_once_is_refused(self) -> None:
        """A window set for each would be several separate free energies,
        each needing its own overlap check."""
        from fastmdxplora.simulation.umbrella import expand_umbrella

        config = self._config()
        config["systems"].append({"system": "1UBQ"})
        with pytest.raises(ValueError, match="one system at many positions"):
            expand_umbrella(config)


class TestReadingTheSamplingBack:
    @staticmethod
    def _written(tmp_path, plan, n=3000):
        for window in plan.windows:
            directory = tmp_path / f"window_{window.index:02d}" / "simulation"
            directory.mkdir(parents=True)
            values = _sample(window.centre, window.force_constant,
                             n=n, seed=window.index)
            directory.joinpath("COLVAR").write_text(
                "#! FIELDS time cv restraint.bias\n"
                + "\n".join(f"{i * 0.1:.1f} {v:.6f} 0.0"
                            for i, v in enumerate(values)),
                encoding="utf-8")

    def test_a_free_energy_comes_back_from_files_on_disk(self, tmp_path) -> None:
        from fastmdxplora.simulation.umbrella import (
            collect_samples,
            compute_pmf,
            plan_windows,
        )

        plan = plan_windows({
            "collective_variable": "distance", "from": -1.6, "to": 1.6,
            "n_windows": 17, "force_constant": 200.0})
        self._written(tmp_path, plan)

        result = compute_pmf(collect_samples(tmp_path, plan), plan)
        x = np.array(result["pmf"]["coordinate"])
        f = np.array(result["pmf"]["free_energy_kjmol"])
        assert f[np.argmin(np.abs(x))] == pytest.approx(10.0, abs=1.5)

    def test_the_approach_to_each_window_is_discarded(self) -> None:
        """A window begins away from where it will settle, and counting the
        approach biases the histogram towards where the run started."""
        import inspect

        from fastmdxplora.simulation import umbrella

        source = inspect.getsource(umbrella.collect_samples)
        assert "equilibration_fraction" in source
        # Flattened, because the sentence wraps and a substring check
        # across a line break tests the layout rather than the meaning.
        prose = " ".join(inspect.getdoc(umbrella.collect_samples).split())
        assert "biases the histogram" in prose

    def test_a_missing_window_stops_it(self, tmp_path) -> None:
        """The windows either side of a missing one have nothing between them
        to stitch through."""
        from fastmdxplora.simulation.umbrella import collect_samples, plan_windows

        plan = plan_windows({
            "collective_variable": "distance", "from": -1.0, "to": 1.0,
            "n_windows": 5, "force_constant": 200.0})
        self._written(tmp_path, plan, n=500)
        import shutil
        shutil.rmtree(tmp_path / "window_02")

        with pytest.raises(FileNotFoundError, match="window_02"):
            collect_samples(tmp_path, plan)


class TestOneWindowInTheRunner:
    def test_the_setting_is_declared(self) -> None:
        from fastmdxplora.config.schema import PHASE_SCHEMAS

        field = PHASE_SCHEMAS["simulation"].get("umbrella")
        assert field is not None
        assert "overlap" in field.help

    def test_it_reaches_the_command_line(self) -> None:
        from fastmdxplora.cli.main import _PHASE_SPEC

        table, _prefix = _PHASE_SPEC["simulate"]
        assert "umbrella" in {dest for _flag, dest, _kw in table}

    def test_the_pipeline_passes_it(self) -> None:
        import inspect

        from fastmdxplora.simulation import pipeline

        assert "umbrella=params.get" in inspect.getsource(pipeline)

    def test_only_one_way_of_moving_a_coordinate_at_a_time(self) -> None:
        """Steering, metadynamics and an umbrella restraint would add their
        forces."""
        import inspect

        from fastmdxplora.simulation import runner

        source = inspect.getsource(runner.run_simulation)
        assert "not more than one" in source or "not more than one" in source.replace("\n", " ")

    def test_a_window_needs_a_position(self) -> None:
        import inspect

        from fastmdxplora.simulation import runner

        assert "needs a `centre`" in inspect.getsource(runner.run_simulation)


class TestItIsActuallyReached:
    """The expansion was written, tested and called by nothing -- the second
    time in a row a piece of this feature worked and was unreachable.

    A config file is the one place every route into the software passes
    through, so expanding there means the command line, Python and the
    browser all get it without each remembering to ask.
    """

    def test_a_config_file_expands_into_windows(self, tmp_path) -> None:
        from fastmdxplora.config.loader import load_config_file

        path = tmp_path / "pmf.yml"
        path.write_text(
            "output: runs/pmf\n"
            "systems:\n"
            "  - system: 181L\n"
            "simulation:\n"
            "  duration_ns: 5\n"
            "  umbrella:\n"
            "    collective_variable: ligand_distance\n"
            '    site_selection: "resid 84 to 121 and name CA"\n'
            "    from: 0.4\n"
            "    to: 2.0\n"
            "    n_windows: 9\n"
            "    force_constant: 1000\n",
            encoding="utf-8")

        config = load_config_file(path)
        assert len(config["systems"]) == 9
        assert config["systems"][0]["id"] == "window_00"
        assert config["simulation"]["duration_ns"] == 5

    def test_a_bad_umbrella_block_reads_as_a_config_error(self, tmp_path) -> None:
        """Rather than a ValueError from somewhere further in."""
        from fastmdxplora.config.loader import ConfigError, load_config_file

        path = tmp_path / "bad.yml"
        path.write_text(
            "output: runs/x\n"
            "systems:\n"
            "  - system: 181L\n"
            "simulation:\n"
            "  umbrella:\n"
            "    collective_variable: distance\n"
            "    from: 0.0\n"
            "    to: 1.0\n"
            "    n_windows: 5\n",
            encoding="utf-8")

        with pytest.raises(ConfigError, match="force_constant"):
            load_config_file(path)

    def test_every_public_piece_of_this_module_is_called(self) -> None:
        """Twice now something here worked and nothing reached it. Checked by
        searching the package rather than by remembering."""
        from pathlib import Path

        from fastmdxplora.simulation import umbrella

        package = Path(__file__).resolve().parents[1] / "src" / "fastmdxplora"
        sources = "\n".join(
            path.read_text(encoding="utf-8") for path in package.rglob("*.py")
            if path.name != "umbrella.py")

        # Entry points a user's run should pass through.
        for name in ("expand_umbrella", "collect_samples", "compute_pmf"):
            assert name in umbrella.__all__, f"{name} is not public"
            assert name in sources, (
                f"{name} is exported and nothing in the package calls it"
            )
