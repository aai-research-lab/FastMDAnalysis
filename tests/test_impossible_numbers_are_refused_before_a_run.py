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
