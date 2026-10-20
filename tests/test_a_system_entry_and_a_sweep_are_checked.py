"""A system entry's blocks and a sweep's values are held to the schema, as a
top-level block is.

Found by the weekly test of 10-20: `simulation: {temprature_K: 350}` in a
system entry ran at 300 K and exited 0; a swept pH of 25 was prepared at 25;
a swept temperature of -50 K ran setup before OpenMM refused it; a repeated
sweep value made two identical runs called replicas; a system's own value of
a swept setting was dropped without a word; a misspelt analysis option was
refused only once the simulation had finished.
"""

from __future__ import annotations

import unittest

from fastmdxplora.config.loader import ConfigError, validate_config


def refusal(config: dict) -> ConfigError:
    try:
        validate_config(config, require_systems=True)
    except ConfigError as exc:
        return exc
    raise AssertionError(f"accepted: {config}")


ONE = [{"system": "1L2Y"}]


class ASystemEntryIsChecked(unittest.TestCase):
    def test_a_misspelt_setting_in_an_entry_is_refused_with_the_nearest(self):
        exc = refusal({"systems": [{"system": "1L2Y",
                                    "simulation": {"temprature_K": 350}}]})
        self.assertEqual(exc.code, "config.option.unknown")
        self.assertIn("temperature_K", str(exc))
        self.assertIn("systems[0]", str(exc))
        validate_config({"systems": [
            {"system": "1UBQ", "id": "wild_type"},
            {"system": "mutant.pdb", "id": "L50A",
             "simulation": {"temperature_K": 320}, "setup": None}]},
            require_systems=True)

    def test_a_misspelt_block_or_unknown_key_in_an_entry_is_refused(self):
        for key in ("simulaton", "sim", "options"):
            exc = refusal({"systems": [{"system": "1L2Y", key: {}}]})
            self.assertEqual(exc.code, "config.option.unknown")

    def test_a_value_out_of_range_in_an_entry_is_refused(self):
        exc = refusal({"systems": [{"system": "1L2Y", "setup": {"ph": 25}}]})
        self.assertEqual(exc.code, "config.option.out_of_range")

    def test_a_block_that_is_not_a_mapping_is_refused_not_a_traceback(self):
        exc = refusal({"systems": [{"system": "1L2Y", "simulation": 5}]})
        self.assertEqual(exc.code, "config.option.wrong_type")

    def test_an_entry_s_analysis_names_and_options_are_checked(self):
        exc = refusal({"systems": [{"system": "1L2Y",
                                    "analysis": {"include": ["rmfs"]}}]})
        self.assertIn(exc.code, ("analysis.unknown", "config.option.not_permitted"))
        exc = refusal({"systems": [{"system": "1L2Y", "analysis": {
            "options": {"hbonds": {"distance_cutof": 0.3}}}}]})
        self.assertEqual(exc.code, "analysis.option.inapplicable")

    def test_a_missing_system_and_a_shared_id_have_their_own_codes(self):
        self.assertEqual(refusal({"systems": [{"id": "a"}]}).code,
                         "config.option.missing_companion")
        exc = refusal({"systems": [{"system": "1L2Y", "id": "a"},
                                   {"system": "1UBQ", "id": "a"}]})
        self.assertEqual(exc.code, "config.option.conflicting")
        self.assertEqual(refusal({"systems": [5]}).code, "config.option.wrong_type")


class ASweepIsChecked(unittest.TestCase):
    def test_each_value_is_held_to_its_setting(self):
        for axis, values, code in (
                ("setup.ph", [7, 25], "config.option.out_of_range"),
                ("simulation.temperature_K", [300, "hot"], "config.option.wrong_type"),
                ("simulation.temperature_K", [-50], "config.option.out_of_range"),
                ("setup.box_shape", ["cube", "sphere"], "config.option.not_permitted"),
                ("simulation.timestep_fs", [2, 0], "config.option.out_of_range")):
            exc = refusal({"systems": ONE, "sweep": {axis: values}})
            self.assertEqual(exc.code, code)

    def test_a_swept_list_of_analyses_is_checked_by_name(self):
        exc = refusal({"systems": ONE,
                       "sweep": {"analysis.include": [["rmsd"], ["rmfs"]]}})
        self.assertIn("rmfs", str(exc))

    def test_a_value_given_twice_is_refused(self):
        for axis, values in (("simulation.random_seed", [1, 1]),
                             ("setup.solvent_padding_nm", [1.0, 1.00])):
            exc = refusal({"systems": ONE, "sweep": {axis: values}})
            self.assertEqual(exc.code, "batch.sweep.invalid")
            self.assertIn("more than once", str(exc))

    def test_an_axis_a_system_also_sets_is_refused_naming_both(self):
        exc = refusal({"systems": [{"system": "1L2Y", "id": "hot",
                                    "simulation": {"temperature_K": 350}}],
                       "sweep": {"simulation.temperature_K": [300, 310]}})
        self.assertEqual(exc.code, "config.option.conflicting")
        self.assertIn("hot", str(exc))
        self.assertIn("350", str(exc))
        validate_config({"systems": [{"system": "1L2Y", "id": "a"},
                                     {"system": "1UBQ", "id": "b",
                                      "setup": {"ph": 7.0}}],
                         "sweep": {"simulation.temperature_K": [300, 310],
                                   "simulation.random_seed": [1, 2, 3]}},
                        require_systems=True)

    def test_an_axis_misspelt_is_refused_with_the_nearest(self):
        exc = refusal({"systems": ONE, "sweep": {"simulation.temprature_K": [300]}})
        self.assertEqual(exc.code, "config.option.unknown")
        self.assertEqual(exc.refusal.details["suggestion"], "simulation.temperature_K")
        self.assertIn("did you mean", str(exc))


class TheRestOfTheConfig(unittest.TestCase):
    def test_a_misspelt_analysis_option_is_refused_before_the_run(self):
        exc = refusal({"systems": ONE, "analysis": {
            "options": {"hbonds": {"distance_cutof": 0.3}}}})
        self.assertEqual(exc.code, "analysis.option.inapplicable")
        self.assertIn("distance_cutoff", str(exc))
        exc = refusal({"systems": ONE, "analysis": {
            "options": {"hbondz": {"distance_cutoff": 0.3}}}})
        self.assertEqual(exc.code, "analysis.unknown")
        self.assertIn("hbonds", str(exc))
        validate_config({"systems": ONE, "analysis": {"options": {
            "cluster": {"methods": ["kmeans"], "n_clusters": 5}}}},
            require_systems=True)

    def test_a_named_force_field_beside_a_list_is_refused(self):
        exc = refusal({"systems": ONE, "setup": {
            "forcefield": "charmm36", "force_field": ["amber14-all.xml"]}})
        self.assertEqual(exc.code, "config.option.conflicting")
        self.assertIn("force_field", str(exc))
        # A resolved_config.yml records both; fed back, it runs.
        from fastmdxplora.setup.forcefields import resolve_forcefield

        xmls = list(resolve_forcefield("amber14").xmls)
        validate_config({"systems": ONE, "setup": {
            "forcefield": "amber14", "force_field": xmls}}, require_systems=True)

    def test_include_beside_include_phase_saying_otherwise_is_refused(self):
        exc = refusal({"systems": ONE, "include": ["setup"],
                       "include_phase": ["setup", "simulation"]})
        self.assertEqual(exc.code, "config.option.conflicting")
        validate_config({"systems": ONE, "include": ["setup"],
                         "include_phase": ["setup"]}, require_systems=True)

    def test_the_phase_conflict_names_the_keys_as_written_now(self):
        exc = refusal({"systems": ONE, "include_phase": ["setup"],
                       "exclude_phase": ["report"]})
        self.assertIn("include_phase", str(exc))
        self.assertIn("exclude_phase", str(exc))

    def test_constraints_are_checked_in_any_letter_case(self):
        exc = refusal({"systems": ONE, "setup": {"constraints": "HBond"}})
        self.assertEqual(exc.code, "config.option.not_permitted")
        self.assertIn("HBonds", str(exc))
        for value in ("HBonds", "hbonds", "None", "AllBonds", "hangles"):
            validate_config({"systems": ONE, "setup": {"constraints": value}},
                            require_systems=True)


class RoundTwo(unittest.TestCase):
    """What the reviews of the first fixes found still let through, or broke."""

    def test_an_axis_whose_other_spelling_the_study_gives_is_refused(self):
        for given, axis in ((("simulation", "pressure_bar", 1.0), "simulation.pressure_atm"),
                            (("simulation", "pressure_atm", 1.0), "simulation.pressure_bar"),
                            (("simulation", "production_steps", 300), "simulation.duration_ns"),
                            (("simulation", "nvt_duration_ns", 0.1), "simulation.nvt_steps")):
            block, key, value = given
            values = [1, 2] if "steps" in axis else [0.001, 0.002]
            exc = refusal({"systems": ONE, block: {key: value}, "sweep": {axis: values}})
            self.assertEqual(exc.code, "config.option.conflicting", axis)
            self.assertIn(key, str(exc))
        exc = refusal({"systems": ONE, "sweep": {"simulation.pressure_bar": [1, 2],
                                                 "simulation.pressure_atm": [1, 2]}})
        self.assertEqual(exc.code, "config.option.conflicting")
        validate_config({"systems": ONE, "sweep": {"simulation.pressure_atm": [1, 2]}},
                        require_systems=True)

    def test_a_swept_value_refused_is_named_by_its_place_in_the_sweep(self):
        from fastmdxplora.remedies import remedy_for

        exc = refusal({"systems": ONE, "sweep": {"setup.ph": [7, 25]}})
        self.assertEqual(exc.refusal.details["option"], "sweep.setup.ph")
        fix = remedy_for(exc.refusal, price=False)
        self.assertIn("`sweep.setup.ph`", str(getattr(fix, "fix", fix)))

    def test_a_record_of_a_study_given_a_force_field_list_replays(self):
        import tempfile
        from pathlib import Path

        from fastmdxplora.cli.main import main
        from fastmdxplora.config import load_config_file, write_resolved_config

        where = Path(tempfile.mkdtemp())
        structure = where / "protein.pdb"
        structure.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000"
                             "  1.00  0.00           N\nEND\n", encoding="utf-8")
        recorded = write_resolved_config({
            "system": str(structure), "output": str(where / "study"),
            "options": {"setup": {"force_field": ["charmm36.xml", "charmm36/water.xml"]}},
        }, where / "study")
        loaded = load_config_file(recorded)
        self.assertEqual(loaded["setup"]["forcefield"], "auto")  # as every record writes it
        validate_config(loaded, require_systems=True)
        self.assertEqual(main(["explore", "-c", str(recorded), "--include-phase",
                               "analysis", "report", "--dry-run"]), 0)


def test_setup_reads_the_default_beside_a_list_of_other_files_as_unset(tmp_path) -> None:
    """Setup refused such a record whenever it ran setup again, as a resume
    that has to start a run over does."""
    from unittest.mock import MagicMock, patch

    from fastmdxplora.setup import prepare
    from fastmdxplora.setup.pipeline import run as setup_run

    pdb = tmp_path / "mini.pdb"
    pdb.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00"
                   "           N\nEND\n", encoding="utf-8")
    orchestrator = MagicMock()
    orchestrator.system = str(pdb)
    orchestrator.output_dir = tmp_path / "proj"
    orchestrator.output_dir.mkdir()
    orchestrator._presenter = None
    out = orchestrator.output_dir / "setup"
    out.mkdir()
    listed = ["charmm36.xml", "charmm36/water.xml"]
    with patch.object(prepare, "prepare_system", return_value={}) as prepared:
        try:
            setup_run(orchestrator=orchestrator, output_dir=out,
                      forcefield="auto", force_field=listed)
        except Exception as exc:  # noqa: BLE001 - only the pair's refusal is the defect
            assert "not both" not in str(exc), exc
    if prepared.called:  # not reached without PDBFixer
        assert prepared.call_args.kwargs["force_field"] == listed


if __name__ == "__main__":
    unittest.main()
