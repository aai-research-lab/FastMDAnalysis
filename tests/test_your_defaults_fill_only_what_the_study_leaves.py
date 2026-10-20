"""Your defaults fill what a study leaves unset, in either spelling of it,
and are said only of the runs that take them.

Found by the weekly test of 10-20: a lab default `pressure_bar: 1.0` beside a
study's `pressure_atm: 1.2` ran at 1.0 bar (bar wins), and the default
`production_steps: 300` beside a study's `duration_ns` ran 300 steps; a lab
default of 310 K was recorded as the reason for the temperature of a system
that set 320 K of its own.
"""

from __future__ import annotations

import unittest
from pathlib import Path

from fastmdxplora.config.defaults_file import YourDefaults, refused_with, with_defaults

WHY = "Our assays run at body temperature"


def lab(**simulation) -> YourDefaults:
    return YourDefaults(Path("/lab/fastmdx-defaults.yml"),
                        {"simulation": dict(simulation)},
                        {"simulation.temperature_K": WHY})


class TheOtherSpellingIsTheStudys(unittest.TestCase):
    def test_pressure_in_atm_keeps_the_file_s_bar_out(self):
        study = {"systems": [{"system": "1L2Y"}],
                 "simulation": {"pressure_atm": 1.2}}
        out, filled = with_defaults(study, lab(pressure_bar=1.0))
        self.assertNotIn("pressure_bar", out["simulation"])
        self.assertEqual(filled, [])
        study = {"systems": [{"system": "1L2Y"}], "simulation": {"pressure_atm": 1.2}}
        out, filled = with_defaults(study, lab(temperature_K=310))
        self.assertEqual(out["simulation"]["temperature_K"], 310)
        self.assertEqual(filled, ["simulation.temperature_K"])

    def test_a_duration_keeps_the_file_s_steps_out(self):
        study = {"systems": [{"system": "1L2Y"}],
                 "simulation": {"duration_ns": 0.0002, "nvt_duration_ns": 0.0001,
                                "npt_duration_ns": 0.0001}}
        out, filled = with_defaults(
            study, lab(production_steps=300, nvt_steps=100, npt_steps=100))
        for key in ("production_steps", "nvt_steps", "npt_steps"):
            self.assertNotIn(key, out["simulation"])
        self.assertEqual(filled, [])

    def test_steps_keep_the_file_s_duration_out(self):
        study = {"systems": [{"system": "1L2Y"}],
                 "simulation": {"production_steps": 500}}
        out, _ = with_defaults(study, lab(duration_ns=10))
        self.assertNotIn("duration_ns", out["simulation"])

    def test_a_named_force_field_keeps_the_file_s_list_out(self):
        defaults = YourDefaults(Path("/lab/fastmdx-defaults.yml"),
                                {"setup": {"forcefield": "charmm36"}})
        study = {"systems": [{"system": "1L2Y"}],
                 "setup": {"force_field": ["amber14-all.xml", "amber14/tip3p.xml"]}}
        out, filled = with_defaults(study, defaults)
        self.assertNotIn("forcefield", out["setup"])
        self.assertIsNone(refused_with(study, out, defaults, require_systems=True))

    def test_a_swept_spelling_keeps_the_file_s_other_out(self):
        study = {"systems": [{"system": "1L2Y"}],
                 "sweep": {"simulation.pressure_atm": [1.0, 2.0]}}
        out, _ = with_defaults(study, lab(pressure_bar=1.0))
        self.assertNotIn("pressure_bar", out.get("simulation", {}))


class ASystemsOwnValueIsItsOwn(unittest.TestCase):
    def test_every_system_setting_it_leaves_the_file_s_value_unused_and_unsaid(self):
        study = {"systems": [{"system": "1L2Y", "id": "a",
                              "simulation": {"temperature_K": 320}}]}
        out, filled = with_defaults(study, lab(temperature_K=310))
        self.assertNotIn("simulation.temperature_K", filled)
        self.assertNotIn("simulation.temperature_K", out.get("decisions") or {})
        self.assertNotIn("temperature_K", out.get("simulation") or {})

    def test_some_systems_setting_it_give_it_to_the_others_and_say_so(self):
        study = {"systems": [{"system": "1L2Y", "id": "a",
                              "simulation": {"temperature_K": 320}},
                             {"system": "1UBQ", "id": "b"}]}
        out, filled = with_defaults(study, lab(temperature_K=310))
        self.assertEqual(out["systems"][0]["simulation"], {"temperature_K": 320})
        self.assertEqual(out["systems"][1]["simulation"], {"temperature_K": 310})
        self.assertNotIn("temperature_K", out.get("simulation") or {})
        why = out["decisions"]["simulation.temperature_K"]["why"]
        self.assertIn("For b", why)
        self.assertIn(WHY, why)
        # The study given is not changed.
        self.assertNotIn("simulation", study["systems"][1])

    def test_a_system_s_other_spelling_keeps_the_file_s_value_from_it(self):
        study = {"systems": [{"system": "1L2Y", "id": "a",
                              "simulation": {"pressure_atm": 2.0}},
                             {"system": "1UBQ", "id": "b"}]}
        out, _ = with_defaults(study, lab(pressure_bar=1.0))
        self.assertNotIn("pressure_bar", out["systems"][0]["simulation"])
        self.assertEqual(out["systems"][1]["simulation"], {"pressure_bar": 1.0})


class EachRunRecordsOnlyItsOwnReason(unittest.TestCase):
    """The study-wide reason "For plain; the other systems give their own."
    was copied into every run's record, so warm's report, at 320 K of its
    own, named the lab's 310 K as why."""

    def _specs(self, systems):
        import tempfile

        from fastmdxplora.batch.explorer import BatchExplorer

        folder = Path(tempfile.mkdtemp())
        pdb = folder / "p.pdb"
        pdb.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00"
                       "           N\nEND\n", encoding="utf-8")
        study = {"output": str(folder / "out"),
                 "systems": [{**entry, "system": str(pdb)} for entry in systems]}
        out, _ = with_defaults(study, lab(temperature_K=310))
        return {spec.system_id: spec for spec in BatchExplorer(config_data=out).run_specs}

    def test_the_system_that_took_it_records_it_and_the_other_does_not(self):
        specs = self._specs([{"id": "plain"},
                             {"id": "warm", "simulation": {"temperature_K": 320}}])
        plain = specs["plain"].study["decisions"]["simulation.temperature_K"]
        self.assertEqual(plain, {"why": f"{WHY}.", "source": "fastmdx-defaults.yml"})
        self.assertNotIn("simulation.temperature_K",
                         specs["warm"].study.get("decisions") or {})
        self.assertEqual(specs["plain"].options["simulation"]["temperature_K"], 310)
        self.assertEqual(specs["warm"].options["simulation"]["temperature_K"], 320)

    def test_systems_with_no_id_are_told_apart_as_their_folders_are(self):
        specs = self._specs([{}, {"simulation": {"temperature_K": 320}}])
        took, own = sorted(specs)
        self.assertIn("simulation.temperature_K", specs[took].study["decisions"])
        self.assertNotIn("simulation.temperature_K", specs[own].study.get("decisions") or {})


if __name__ == "__main__":
    unittest.main()
