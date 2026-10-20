"""A number a setting cannot be is refused with the Config, before a run.

Found by the weekly test of 10-20: a timestep of zero passed every check and
stopped the run on a division by zero; a negative production ran "0 frames"
and exited 0; an interval of zero was reported as an unstable simulation; a
negative hydrogen mass built a system of negative masses; `ph: .nan` was
prepared at a pH of NaN.
"""

from __future__ import annotations

import math
import unittest

from fastmdxplora.config.loader import ConfigError, validate_config


def refusal(config: dict) -> ConfigError:
    try:
        validate_config(config)
    except ConfigError as exc:
        return exc
    raise AssertionError(f"accepted: {config}")


class ImpossibleNumbersAreRefused(unittest.TestCase):
    CASES = (
        ("simulation", "timestep_fs", 0),
        ("simulation", "timestep_fs", -1.0),
        ("simulation", "production_steps", -5000),
        ("simulation", "nvt_steps", -5),
        ("simulation", "npt_steps", -1),
        ("simulation", "trajectory_interval_steps", 0),
        ("simulation", "trajectory_interval_steps", -10),
        ("simulation", "state_interval_steps", 0),
        ("simulation", "pressure_atm", -1),
        ("simulation", "friction_per_ps", -1),
        ("simulation", "barostat_frequency", 0),
        ("simulation", "integrator_error_tolerance", 0),
        ("setup", "hydrogen_mass_amu", -1),
        ("setup", "hydrogen_mass_amu", 0),
        ("setup", "ewald_error_tolerance", 0),
        ("analysis", "stride", 0),
        ("analysis", "first", -1),
        ("execution", "workers", 0),
        ("execution", "workers", -2),
    )

    def test_each_is_refused_as_out_of_range_naming_the_setting(self):
        for block, key, value in self.CASES:
            exc = refusal({block: {key: value}})
            self.assertEqual(exc.code, "config.option.out_of_range")
            self.assertIn(key, str(exc))

    def test_a_value_that_is_not_finite_is_refused_whatever_the_bounds(self):
        for block, key, value in (("setup", "ph", math.nan),
                                  ("simulation", "temperature_K", math.inf),
                                  ("setup", "switch_distance_nm", math.nan),
                                  ("simulation", "pressure_atm", math.inf)):
            exc = refusal({block: {key: value}})
            self.assertEqual(exc.code, "config.option.out_of_range")
            self.assertIn("finite", str(exc))

    def test_last_at_or_before_first_is_refused(self):
        exc = refusal({"analysis": {"first": 50, "last": 10}})
        self.assertEqual(exc.code, "config.option.conflicting")
        self.assertIn("first", str(exc))
        validate_config({"analysis": {"first": 10, "last": 50}})
        validate_config({
            "setup": {"hydrogen_mass_amu": 4.0, "ph": 7.0},
            "simulation": {"timestep_fs": 4, "nvt_steps": 0, "npt_steps": 0,
                           "production_steps": 1000,
                           "trajectory_interval_steps": 100,
                           "checkpoint_interval_steps": 0,
                           "friction_per_ps": 0, "pressure_atm": 1.0},
            "analysis": {"stride": 1, "first": 0},
            "execution": {"workers": 1},
        })

    def test_a_trajectory_interval_longer_than_production_is_refused(self):
        exc = refusal({"simulation": {"production_steps": 500,
                                      "trajectory_interval_steps": 5000}})
        self.assertEqual(exc.code, "config.option.conflicting")
        self.assertIn("trajectory_interval_steps", str(exc))
        # From a duration at the timestep: 0.001 ns at 2 fs is 500 steps.
        exc = refusal({"simulation": {"duration_ns": 0.001,
                                      "trajectory_interval_steps": 5000}})
        self.assertEqual(exc.code, "config.option.conflicting")
        # As long as production, one frame: accepted. No production at all
        # writes no trajectory and asks for none.
        validate_config({"simulation": {"production_steps": 500,
                                        "trajectory_interval_steps": 500}})
        validate_config({"simulation": {"production_steps": 0,
                                        "trajectory_interval_steps": 5000}})

    def test_the_interval_is_checked_for_each_system_and_sweep_value(self):
        exc = refusal({
            "simulation": {"trajectory_interval_steps": 1000},
            "systems": [{"system": "1L2Y"},
                        {"system": "1UBQ", "id": "short",
                         "simulation": {"production_steps": 200}}],
        })
        self.assertEqual(exc.code, "config.option.conflicting")
        self.assertIn("short", str(exc))
        exc = refusal({
            "simulation": {"trajectory_interval_steps": 1000},
            "systems": [{"system": "1L2Y"}],
            "sweep": {"simulation.production_steps": [100000, 300]},
        })
        self.assertEqual(exc.code, "config.option.conflicting")


class AnErrorThatIsNotInstabilityIsNotCalledOne(unittest.TestCase):
    def test_a_modulo_by_zero_is_said_as_it_is(self):
        from fastmdxplora.refusals import UnstableRun
        from fastmdxplora.simulation import runner

        made = runner._integration_failed(
            "production", ZeroDivisionError("integer modulo by zero"), object())
        self.assertNotIsInstance(made, UnstableRun)
        self.assertEqual(made.refusal.code, "simulation.run.engine_failed")
        self.assertIn("modulo by zero", str(made))
        self.assertNotIn("Minimize", str(made))

    def test_a_nan_is_still_an_instability(self):
        from fastmdxplora.refusals import UnstableRun
        from fastmdxplora.simulation import runner

        class NoState:
            def __getattr__(self, name):
                raise RuntimeError("no context")

        made = runner._integration_failed(
            "production", Exception("Particle coordinate is NaN."), NoState())
        self.assertIsInstance(made, UnstableRun)


if __name__ == "__main__":
    unittest.main()


def test_one_worker_is_said_as_one(tmp_path, capsys) -> None:
    """The dry run said "parallel (1 workers)"."""
    from fastmdxplora.batch.explorer import BatchExplorer

    structure = tmp_path / "protein.pdb"
    structure.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "END\n", encoding="utf-8")
    BatchExplorer(config_data={
        "output": str(tmp_path / "out"),
        "systems": [{"id": "a", "system": str(structure)}],
        "sweep": {"simulation.random_seed": [1, 2]},
        "execution": {"mode": "parallel", "workers": 1},
    }).dry_run()
    said = capsys.readouterr().out
    assert "parallel (1 worker)" in said, said


class RoundTwo(unittest.TestCase):
    """What the reviews of the first fixes found still let through."""

    def test_a_hydrogen_mass_at_hydrogen_s_or_over_4_amu_is_refused(self):
        for value in (1.0, 1.008, 4.5, 6):
            exc = refusal({"setup": {"hydrogen_mass_amu": value}})
            self.assertEqual(exc.code, "config.option.out_of_range", value)
        for value in (1.5, 3, 4.0):
            validate_config({"setup": {"hydrogen_mass_amu": value}})

    def test_a_setting_not_finite_inside_an_analysis_option_is_refused(self):
        exc = refusal({"analysis": {"options": {"hbonds": {"distance_cutoff": math.nan}}}})
        self.assertEqual(exc.code, "config.option.out_of_range")
        self.assertIn("finite", str(exc))

    def test_last_at_the_first_frame_with_no_first_is_refused(self):
        exc = refusal({"analysis": {"last": 0}})
        self.assertEqual(exc.code, "config.option.conflicting")
        validate_config({"analysis": {"last": 5}})

    def test_a_stage_given_as_steps_and_a_duration_that_disagree_is_refused(self):
        for steps, ns in (("production_steps", "duration_ns"),
                          ("nvt_steps", "nvt_duration_ns"),
                          ("npt_steps", "npt_duration_ns")):
            exc = refusal({"simulation": {steps: 300, ns: 10}})
            self.assertEqual(exc.code, "config.option.conflicting")
            self.assertIn(steps, str(exc))
            # The same length said twice is one length.
            validate_config({"simulation": {steps: 500, ns: 0.001}})
        exc = refusal({"systems": [{"system": "1L2Y",
                                    "simulation": {"duration_ns": 0.002}}],
                       "simulation": {"production_steps": 300}})
        self.assertEqual(exc.code, "config.option.conflicting")

    def test_an_exclusive_minimum_is_said_as_more_than(self):
        from fastmdxplora.config.describe import describe_field, schema_as_json
        from fastmdxplora.config.phase_settings import _bounds
        from fastmdxplora.config.schema import SIMULATION

        field = SIMULATION.get("timestep_fs")
        said = describe_field(field)
        self.assertIn("more than 0", said)
        self.assertNotIn("at least 0", said)
        self.assertIn("more than", _bounds(field))
        entry = schema_as_json(phases=["simulation"])["simulation"]["timestep_fs"]
        self.assertEqual(entry.get("exclusiveMinimum"), 0.0)
        self.assertNotIn("minimum", entry)



def _system_of(masses):
    """A system and topology whose atoms carry these masses, as setup's
    `createSystem` hands them over: a methyl carbon at -2.966 amu is what a
    hydrogen mass of 6 built."""
    import pytest

    pytest.importorskip("openmm")
    import openmm
    from openmm import app

    topology = app.Topology()
    chain = topology.addChain()
    residue = topology.addResidue("ALA", chain)
    system = openmm.System()
    for i, mass in enumerate(masses):
        topology.addAtom(f"A{i}", app.element.carbon, residue)
        system.addParticle(mass)
    return system, topology


def test_setup_refuses_a_system_built_with_a_mass_at_or_below_zero() -> None:
    """The validator's bound is the first guard; setup's holds whatever
    reaches it."""
    import pytest

    from fastmdxplora.refusals import StudyError
    from fastmdxplora.setup import prepare

    system, topology = _system_of([12.011 - 3 * 5.0, 6.0, 14.007])
    with pytest.raises(StudyError) as caught:
        prepare._refuse_masses_not_positive(system, topology, 6.0)
    assert caught.value.code == "setup.system.mass_not_positive"
    assert "ALA1:A0" in str(caught.value) and "4.0 or less" in str(caught.value)
    system, topology = _system_of([3.011, 4.0, 14.007, 0.0])
    with pytest.raises(StudyError):
        prepare._refuse_masses_not_positive(system, topology, None)
    # A virtual site has no mass by design.
    import openmm

    system.setVirtualSite(3, openmm.TwoParticleAverageSite(0, 1, 0.5, 0.5))
    prepare._refuse_masses_not_positive(system, topology, None)


def test_an_interval_that_does_not_divide_production_is_worth_knowing() -> None:
    from fastmdxplora.advisories import advise

    said = [a for a in advise({}, {"production_steps": 500,
                                   "trajectory_interval_steps": 300})
            if a.setting == "trajectory_interval_steps"]
    assert len(said) == 1, said
    text = said[0].as_text()
    assert "1 frame will be written" in text and "last 200 steps" in text, text
    assert "such as 250" in text, text
    assert not [a for a in advise({}, {"production_steps": 1000,
                                       "trajectory_interval_steps": 100})
                if a.setting == "trajectory_interval_steps"]


def test_a_study_whose_phase_spellings_disagree_is_read_as_naming_none(tmp_path) -> None:
    """The timeline read the config's phases through the loader's settling,
    which now refuses two spellings that disagree; it raised in the GUI."""
    from fastmdxplora.gui.telemetry import phases_asked

    (tmp_path / "exploration.yml").write_text(
        "systems: [{system: 1L2Y}]\ninclude: [setup]\n"
        "include_phase: [setup, simulation]\n", encoding="utf-8")
    assert phases_asked(tmp_path) == []
