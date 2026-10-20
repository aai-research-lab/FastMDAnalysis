"""A study with a budget runs the phases it asks for, and records them.

A study carrying `budget_hours` runs in two stages: setup, a price, then
the rest. Each stage set one phase key on top of the study's own, so a
study that named its phases (`include_phase` from the GUI or the template,
`exclude_phase: [analysis, report]` by hand) was refused after setup for
"both 'include' and 'exclude'", a conflict the person never wrote. A study
that named none ran, and its `resolved_config.yml` said `exclude_phase:
[setup]` with the stage's own `setup_from`, so replaying the record skipped
the setup the study asked for.
"""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import yaml

from fastmdxplora.agent import run_in_stages


class TestTheStudySPhases(unittest.TestCase):

    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.sent: list[dict] = []
        from fastmdxplora.cost import calibrate

        calibrate(particles=30_000, steps=5_000, seconds=42.0,
                  platform_name="CUDA", precision="mixed")

    def explore(self, *, config, output_dir):
        """What the real run does with the config it is given: validates it
        (where the conflict was refused), writes the record, and for setup
        the count the price is read from."""
        from fastmdxplora.config import validate_config, write_resolved_config

        checked = json.loads(json.dumps(config))
        validate_config(checked)
        self.sent.append(checked)
        out = Path(output_dir)
        if checked.get("include_phase") == ["setup"]:
            (out / "setup").mkdir(parents=True, exist_ok=True)
            (out / "setup" / "setup_parameters.json").write_text(
                json.dumps({"n_atoms_solvated": 60_000}))
        write_resolved_config({
            "systems": checked["systems"], "output": str(out),
            "include_phase": checked.get("include_phase"),
            "exclude_phase": checked.get("exclude_phase"),
            "options": {"simulation": dict(checked.get("simulation") or {})},
        }, out)

    def study(self, **phases):
        return {"systems": [{"id": "a", "system": "x.pdb"}],
                "simulation": {"production_steps": 1000,
                               "nvt_steps": 0, "npt_steps": 0},
                **phases}

    def run_it(self, config):
        return run_in_stages(config, self.root, budget_hours=1000,
                             platform_name="CUDA", explore=self.explore)

    def recorded(self):
        return yaml.safe_load((self.root / "resolved_config.yml").read_text())

    def test_a_study_naming_its_phases_runs_through(self):
        for phases in ({"include_phase": ["setup", "simulation", "analysis", "report"]},
                       {"include_phase": ["setup", "simulation"]},
                       {"exclude_phase": ["analysis", "report"]},
                       {"exclude_phase": ["report"]},
                       {"include": ["setup", "simulation", "analysis"]},
                       {"exclude": ["report"]},
                       {}):
            self.root = Path(tempfile.mkdtemp())
            self.sent = []
            staged = self.run_it(self.study(**phases))
            self.assertIsNone(staged.refusal,
                              f"{phases}: {staged.refusal and staged.refusal.message}")
            self.assertTrue(staged.simulated, phases)

    def test_each_stage_runs_only_its_own_share_of_the_study_s_phases(self):
        self.run_it(self.study(exclude_phase=["analysis", "report"]))
        self.assertEqual(self.sent[0]["include_phase"], ["setup"])
        self.assertEqual(self.sent[1]["include_phase"], ["simulation"])
        for sent in self.sent:
            self.assertFalse(sent.get("exclude_phase"))

    def test_the_record_says_the_study_s_phases(self):
        self.run_it(self.study(exclude_phase=["report"]))
        recorded = self.recorded()
        self.assertEqual(recorded.get("exclude_phase"), ["report"])
        self.assertNotIn("include_phase", recorded)

    def test_the_record_of_a_study_naming_none_names_none(self):
        self.run_it(self.study())
        recorded = self.recorded()
        self.assertNotIn("exclude_phase", recorded)
        self.assertNotIn("include_phase", recorded)

    def test_the_record_keeps_no_setup_from_the_stage_gave_it(self):
        self.run_it(self.study(include_phase=["setup", "simulation"]))
        recorded = self.recorded()
        self.assertEqual(recorded["include_phase"], ["setup", "simulation"])
        simulation = recorded.get("simulation") or {}
        self.assertFalse(simulation.get("setup_from"))
        self.assertFalse(simulation.get("prepared_from"))
        self.assertTrue(recorded and (self.root / "resolved_config.yml")
                        .read_text().startswith("# FastMDXplora resolved configuration"))
        # A prepared system the study itself named stays named.
        elsewhere = Path(tempfile.mkdtemp()) / "prepared"
        elsewhere.mkdir()
        (elsewhere / "setup_parameters.json").write_text(
            json.dumps({"n_atoms_solvated": 60_000}))
        self.root = Path(tempfile.mkdtemp())
        config = self.study(include_phase=["setup", "simulation"])
        config["simulation"]["setup_from"] = str(elsewhere)
        self.run_it(config)
        self.assertEqual(self.recorded()["simulation"]["setup_from"], str(elsewhere))

    def test_a_study_of_setup_alone_is_not_priced_or_run_again(self):
        staged = self.run_it(self.study(include_phase=["setup"]))
        self.assertIsNone(staged.refusal)
        self.assertEqual(len(self.sent), 1)


class TestAnAIAppIsToldWhatTheGUIWarnsOf(unittest.TestCase):
    """check_study and save_study say what the GUI's Check it says under
    Worth knowing, from the same function: a 4 fs timestep with no
    constraints was warned of there and passed here as it stands."""

    RISKY = ("systems:\n  - system: ghg.pdb\nsimulation:\n  timestep_fs: 4\n"
             "setup:\n  constraints: None\n")

    def setUp(self):
        import os

        from fastmdxplora.mcp import App, Workspace
        from tests._mcp_wire import Wire
        from tests.test_an_ai_app_reads_and_checks_studies import _structure

        self.tmp = Path(tempfile.mkdtemp())
        self._env = {k: os.environ.get(k) for k in ("FASTMDXPLORA_CACHE_DIR",
                                                    "FASTMDXPLORA_CONFIG_DIR")}
        os.environ["FASTMDXPLORA_CACHE_DIR"] = str(self.tmp / "cache")
        os.environ["FASTMDXPLORA_CONFIG_DIR"] = str(self.tmp / "settings")
        self.root = self.tmp / "work"
        self.root.mkdir()
        _structure(self.root / "ghg.pdb")
        self.wire = Wire(App(Workspace.at(self.root)).server())

    def tearDown(self):
        import os

        self.wire.close()
        for key, value in self._env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

    def said(self, tool, **arguments):
        result = self.wire.request("tools/call", {"name": tool, "arguments": arguments})["result"]
        self.assertFalse(result["isError"], result)
        return result["content"][0]["text"]

    def gui_says(self):
        from fastmdxplora.gui.config_builder import check_config

        knowing = check_config(yaml.safe_load(self.RISKY))["worth_knowing"]
        self.assertTrue(knowing, "the GUI warns of nothing here: the case is not the one")
        return knowing

    def test_check_study_says_it(self):
        said = self.said("check_study", config=self.RISKY)
        self.assertIn("\nWorth knowing:\n", said)
        for line in self.gui_says():
            self.assertIn(f"  {line}", said)
        # And a config with nothing worth knowing has no heading for it.
        plain = self.said("check_study", config="systems:\n  - system: ghg.pdb\n")
        self.assertNotIn("Worth knowing", plain)

    def test_save_study_says_it(self):
        said = self.said("save_study", config=self.RISKY, name="risky")
        self.assertIn("\nWorth knowing:\n", said)
        for line in self.gui_says():
            self.assertIn(f"  {line}", said)


if __name__ == "__main__":
    unittest.main()
