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

        directories = {w.index: tmp_path / f"window_{w.index:02d}"
                       for w in plan.windows}
        result = compute_pmf(collect_samples(directories), plan)
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

        directories = {w.index: tmp_path / f"window_{w.index:02d}"
                       for w in plan.windows}
        with pytest.raises(FileNotFoundError, match="window 2"):
            collect_samples(directories)


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


class TestNothingIsSmuggledThroughTheConfig:
    """The expansion stashed the plan under `_umbrella_plan`, and validation
    refused the unknown key -- correctly, since a config is what somebody
    wrote and not a place to hide state.

    Each window carries its own block, so the set is rebuilt from the runs.
    """

    @staticmethod
    def _written(tmp_path):
        path = tmp_path / "pmf.yml"
        path.write_text(
            "output: pmftest\n"
            "systems:\n"
            "  - system: 181L\n"
            "simulation:\n"
            "  umbrella:\n"
            "    collective_variable: ligand_distance\n"
            '    site_selection: "resid 84 to 121 and name CA"\n'
            "    from: 0.3\n"
            "    to: 1.5\n"
            "    n_windows: 7\n"
            "    force_constant: 1000\n",
            encoding="utf-8")
        return path

    def test_an_expanded_config_validates(self, tmp_path) -> None:
        from fastmdxplora.config.loader import load_config_file, validate_config

        config = load_config_file(self._written(tmp_path))
        validate_config(config, require_systems=True)   # raises if not

    def test_no_private_key_is_added(self, tmp_path) -> None:
        from fastmdxplora.config.loader import load_config_file

        config = load_config_file(self._written(tmp_path))
        assert not [k for k in config if k.startswith("_")], (
            "the config should hold nothing but what a user could have written"
        )

    def test_the_plan_is_rebuilt_from_the_windows(self, tmp_path) -> None:
        from fastmdxplora.config.loader import load_config_file
        from fastmdxplora.simulation.umbrella import plan_from_expanded

        plan = plan_from_expanded(load_config_file(self._written(tmp_path)))
        assert plan is not None
        assert len(plan.windows) == 7
        assert plan.windows[0].centre == pytest.approx(0.3)
        assert plan.windows[-1].centre == pytest.approx(1.5)
        assert plan.collective_variable == "ligand_distance"

    def test_a_plain_study_rebuilds_nothing(self) -> None:
        from fastmdxplora.simulation.umbrella import plan_from_expanded

        assert plan_from_expanded(
            {"systems": [{"system": "1UBQ"}, {"system": "181L"}]}) is None

    def test_the_explorer_rebuilds_rather_than_reads_a_stash(self) -> None:
        import inspect

        from fastmdxplora.batch import explorer

        source = inspect.getsource(explorer.BatchExplorer._maybe_build_pmf)
        assert "plan_from_expanded" in source
        assert "_umbrella_plan" not in source


class TestAWindowInheritsTheStudysSettings:
    """A per-system block replaces the top-level one rather than merging, so a
    block holding only the umbrella settings discarded the step counts. A run
    asking for five thousand production steps ran a million.
    """

    @staticmethod
    def _expanded():
        from fastmdxplora.simulation.umbrella import expand_umbrella

        return expand_umbrella({
            "output": "x",
            "systems": [{"system": "181L"}],
            "simulation": {
                "nvt_steps": 500, "npt_steps": 500, "production_steps": 5000,
                "timestep_fs": 1.0,
                "umbrella": {
                    "collective_variable": "ligand_distance",
                    "site_selection": "resid 84 to 121 and name CA",
                    "from": 0.3, "to": 1.5, "n_windows": 7,
                    "force_constant": 1000}},
        })

    def test_the_step_counts_reach_every_window(self) -> None:
        for entry in self._expanded()["systems"]:
            simulation = entry["simulation"]
            assert simulation["production_steps"] == 5000
            assert simulation["nvt_steps"] == 500
            assert simulation["timestep_fs"] == 1.0

    def test_and_each_still_holds_its_own_position(self) -> None:
        centres = [e["simulation"]["umbrella"]["centre"]
                   for e in self._expanded()["systems"]]
        assert len(set(centres)) == 7


class TestFindingTheLigand:
    """A ligand variable in an umbrella window refused for want of something
    the run knew: the runner never received the residue name, and neither did
    the metadynamics path -- it had simply never been run with one.
    """

    @staticmethod
    def _system(extra):
        import mdtraj as md

        top = md.Topology()
        chain = top.add_chain()
        for index in range(10):
            residue = top.add_residue("ALA", chain, resSeq=index + 1)
            for name in ("N", "CA", "C"):
                top.add_atom(name, md.element.carbon, residue)
        for name in extra:
            other = top.add_chain()
            residue = top.add_residue(name, other, resSeq=900)
            for index in range(6):
                top.add_atom(f"C{index}", md.element.carbon, residue)
        waters = top.add_chain()
        for index in range(3):
            residue = top.add_residue("HOH", waters, resSeq=1000 + index)
            top.add_atom("O", md.element.oxygen, residue)
        return top

    def test_one_candidate_is_the_ligand(self) -> None:
        from fastmdxplora.simulation.metadynamics import detect_ligand

        assert detect_ligand(self._system(["BNZ"])) == "BNZ"

    def test_no_candidate_is_not_a_ligand(self) -> None:
        from fastmdxplora.simulation.metadynamics import detect_ligand

        assert detect_ligand(self._system([])) is None

    def test_two_candidates_is_a_question_the_topology_cannot_answer(self) -> None:
        from fastmdxplora.simulation.metadynamics import detect_ligand

        assert detect_ligand(self._system(["BNZ", "LIG"])) is None

    def test_the_refusal_says_to_name_it(self) -> None:
        from fastmdxplora.simulation.metadynamics import plan_from_config

        with pytest.raises(ValueError, match="ligand_resname"):
            plan_from_config({"collective_variable": "ligand_distance",
                              "site_selection": "name CA", "sigma": 0.05,
                              "unbounded": True},
                             self._system(["BNZ", "LIG"]))

    def test_lipids_are_not_mistaken_for_a_ligand(self) -> None:
        from fastmdxplora.simulation.metadynamics import detect_ligand

        assert detect_ligand(self._system(["POP"])) is None


class TestAFailedWindowStopsTheStudy:
    """Umbrella windows are one measurement, not several systems.

    A campaign should carry on when one system fails -- the others are still
    results. A free energy cannot be computed at all if a window is missing,
    so continuing spends hours producing runs that get thrown away.
    """

    @staticmethod
    def _umbrella():
        from fastmdxplora.simulation.umbrella import expand_umbrella

        return expand_umbrella({
            "output": "/tmp/u", "systems": [{"system": "181L"}],
            "simulation": {"umbrella": {
                "collective_variable": "ligand_distance",
                "site_selection": "name CA", "from": 0.3, "to": 1.5,
                "n_windows": 5, "force_constant": 1000}}})

    def test_an_umbrella_study_stops_on_the_first_failure(self, tmp_path) -> None:
        from fastmdxplora.batch.explorer import BatchExplorer

        explorer = BatchExplorer(config_data=self._umbrella(),
                                 output_dir=str(tmp_path))
        assert explorer.continue_on_error is False

    def test_an_ordinary_campaign_carries_on(self, tmp_path) -> None:
        from fastmdxplora.batch.explorer import BatchExplorer

        explorer = BatchExplorer(
            config_data={"output": "x",
                         "systems": [{"system": "1UBQ"}, {"system": "181L"}]},
            output_dir=str(tmp_path))
        assert explorer.continue_on_error is True

    def test_it_can_still_be_asked_for_explicitly(self, tmp_path) -> None:
        """Somebody who wants the surviving windows anyway should get them."""
        from fastmdxplora.batch.explorer import BatchExplorer

        explorer = BatchExplorer(config_data=self._umbrella(),
                                 output_dir=str(tmp_path),
                                 continue_on_error=True)
        assert explorer.continue_on_error is True

    def test_the_message_says_why_the_rest_were_not_run(self) -> None:
        """The difference between a reader thinking something crashed and
        knowing the rest was not attempted on purpose."""
        from fastmdxplora.batch.explorer import _not_submitted

        umbrella = _not_submitted("window_02", umbrella=True)
        assert "a free energy needs every window" in umbrella
        assert "thrown away" in umbrella

        campaign = _not_submitted("run_02", umbrella=False)
        assert "continue_on_error=False" in campaign

    def test_the_expanded_config_is_still_json(self) -> None:
        """The window centres came out of numpy, and a numpy scalar
        serialises as an opaque object in YAML and not at all in JSON -- so a
        config written back out would not reload."""
        import json

        json.dumps(self._umbrella())
        for entry in self._umbrella()["systems"]:
            centre = entry["simulation"]["umbrella"]["centre"]
            assert type(centre) is float


class TestItLooksWhereTheRunsActuallyWent:
    """The first version guessed the path -- `<output>/window_00/simulation/
    COLVAR` -- and the runs go to `<output>/runs/window-00/simulation/COLVAR`:
    under a `runs` directory, with the identifier slugged.

    Two mistakes in one path, neither visible until a real study finished and
    the recombination found nothing. The caller knows where it put things.
    """

    @staticmethod
    def _explorer(tmp_path):
        from fastmdxplora.batch.explorer import BatchExplorer

        config = tmp_path / "c.yml"
        config.write_text(
            "output: out\n"
            "include: [setup, simulation]\n"
            "systems:\n"
            "  - system: 181L\n"
            "simulation:\n"
            "  umbrella:\n"
            "    collective_variable: distance\n"
            '    selection_a: "name CA"\n'
            '    selection_b: "name CB"\n'
            "    from: -1.6\n"
            "    to: 1.6\n"
            "    n_windows: 17\n"
            "    force_constant: 200\n",
            encoding="utf-8")
        return BatchExplorer(config=config, output_dir=str(tmp_path / "out"))

    def test_a_free_energy_comes_out_of_a_real_layout(self, tmp_path) -> None:
        import json

        explorer = self._explorer(tmp_path)
        for spec in explorer.run_specs:
            block = spec.options["simulation"]["umbrella"]
            simulation = explorer._run_output_dir(spec) / "simulation"
            simulation.mkdir(parents=True)
            values = _sample(block["centre"], block["force_constant"],
                             n=3000, seed=block["index"])
            simulation.joinpath("COLVAR").write_text(
                "#! FIELDS time cv bias\n"
                + "\n".join(f"{i * 0.1:.1f} {v:.6f} 0.0"
                            for i, v in enumerate(values)),
                encoding="utf-8")

        explorer._maybe_build_pmf()

        written = json.loads(
            (tmp_path / "out" / "pmf.json").read_text(encoding="utf-8"))
        assert written["refused"] is None
        x = np.array(written["pmf"]["coordinate"])
        f = np.array(written["pmf"]["free_energy_kjmol"])
        assert f[np.argmin(np.abs(x))] == pytest.approx(10.0, abs=1.5)

    def test_the_directories_come_from_the_explorer(self) -> None:
        """Not from reconstructing a path, which is what got it wrong."""
        import inspect

        from fastmdxplora.batch import explorer
        from fastmdxplora.simulation import umbrella

        assert "_run_output_dir(spec)" in inspect.getsource(
            explorer.BatchExplorer._maybe_build_pmf)
        # And the collector no longer builds one.
        # The body, not the docstring -- which explains the mistake and
        # therefore contains the very string this looks for.
        source = inspect.getsource(umbrella.collect_samples)
        body = source[source.index('"""', source.index('"""') + 3):]
        assert "window_" not in body, (
            "the collector should not know how a run directory is named"
        )

    def test_a_missing_window_names_the_path_it_looked_at(self, tmp_path) -> None:
        """So the next path mistake is one line of output rather than an
        empty result."""
        from fastmdxplora.simulation.umbrella import collect_samples

        with pytest.raises(FileNotFoundError, match="COLVAR"):
            collect_samples({0: tmp_path / "nowhere"})


class TestParallelRunsShareOneTerminal:
    """Three workers printing the wordmark at once produced a screen of
    interleaved fragments -- characters from three banners arriving mixed
    together, which is not three runs starting but one unreadable page.

    Each run writes its own log file. The terminal belongs to the study.
    """

    def test_a_worker_prints_no_banner(self) -> None:
        import io
        import os

        from fastmdxplora.utils.presenter import SessionPresenter

        previous = os.environ.get("FASTMDX_LOG_STYLE")
        os.environ["FASTMDX_LOG_STYLE"] = "plain"
        try:
            out = io.StringIO()
            SessionPresenter(stream=out).banner(
                System="181L", Output="x", Version="2.3.0")
            assert out.getvalue() == ""
        finally:
            if previous is None:
                os.environ.pop("FASTMDX_LOG_STYLE", None)
            else:
                os.environ["FASTMDX_LOG_STYLE"] = previous

    def test_an_ordinary_run_still_prints_one(self) -> None:
        import io
        import os

        from fastmdxplora.utils.presenter import SessionPresenter

        previous = os.environ.pop("FASTMDX_LOG_STYLE", None)
        try:
            out = io.StringIO()
            SessionPresenter(stream=out).banner(
                System="181L", Output="x", Version="2.3.0")
            assert "181L" in out.getvalue()
        finally:
            if previous is not None:
                os.environ["FASTMDX_LOG_STYLE"] = previous

    def test_the_worker_asks_for_it(self) -> None:
        import inspect

        from fastmdxplora.batch import explorer

        source = inspect.getsource(explorer._execute_run)
        assert "FASTMDX_LOG_STYLE" in source
        assert "share one terminal" in source

    def test_plumeds_own_output_goes_to_the_run_log(self) -> None:
        """PLUMED prints its banner from C++, straight to the file
        descriptor, so silencing ours left three copies of its own
        interleaved in the terminal. Python-level redirection does not reach
        a C write -- the descriptor has to move, which is the same lesson a
        Markdown renderer taught by printing its installation guide from a
        dynamic loader.
        """
        import ctypes
        import inspect
        import os
        import pathlib
        import tempfile
        import textwrap

        from fastmdxplora.batch import explorer

        source = inspect.getsource(explorer._execute_run)
        block = source[source.index("    class _QuietBanner:"):
                       source.index("    # Imported here")]
        namespace: dict = {}
        exec("import os as _os\n" + textwrap.dedent(block), namespace)

        libc = ctypes.CDLL(None)
        with tempfile.TemporaryDirectory() as directory:
            log = os.path.join(directory, "run", "run.log")
            with namespace["_QuietBanner"](log):
                libc.puts(b"PLUMED: pretending to start")
                libc.fflush(None)
            assert "PLUMED" in pathlib.Path(log).read_text(encoding="utf-8")

    def test_and_puts_it_back(self) -> None:
        """The worker is called in-process by the tests as well as in a
        subprocess by a real run. Setting a global and not restoring it made
        seventeen unrelated tests fail."""
        import inspect

        from fastmdxplora.batch import explorer

        source = inspect.getsource(explorer._execute_run)
        assert "__exit__" in source
        # The call, whatever it now takes as an argument -- pinning the exact
        # text made this fail the moment the guard learned where to write.
        assert "with _QuietBanner(" in source
