"""The Config Builder runs what it shows.

Found by the weekly test of 10-20 (area 2, New study): with a
fastmdx-defaults.yml beside the studies, the Builder's plan said 300 K
where the run used 310 K, and a value typed to set FastMDXplora's own
default back was left out of the file; `fastmdx gui` started in a folder put
its studies, read its relative paths and looked for its defaults one folder
up; Check it said "Runs." of a Config with no systems; Run as it stands went
to the Results box's old value, not the Config's own `output:`; the plan
named the pressure in atm where the run used bar; a Run without PDBFixer was
refused even for a structure already fixed, and left a folder behind; Full
config ticked still ran the short file; and the Changed card said "not set"
of a setting whose default is 2 ns. Each test here says what a person sees
or what a run is given.
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pytest
import yaml

ATOM = "ATOM  {:5d}  {:<3s} {:>3s} {}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00          {:>2s}\n"


def _peptide(n: int = 3) -> str:
    lines, serial = [], 1
    for i in range(n):
        for k, name in enumerate(("N", "CA", "C", "O")):
            lines.append(ATOM.format(serial, name, "GLY", "A", i + 1, 3.6 * i + 0.4 * k,
                                     0.0, 0.0, name[0]))
            serial += 1
    return "".join(lines) + "END\n"


@pytest.fixture(autouse=True)
def _own_folders(monkeypatch, tmp_path):
    monkeypatch.setenv("FASTMDXPLORA_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    # Your defaults are looked for up to the home folder: kept to this test.
    monkeypatch.setenv("HOME", str(tmp_path))


def _post(url: str, body: dict) -> dict:
    request = urllib.request.Request(url, data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"},
                                     method="POST")
    with urllib.request.urlopen(request) as response:
        return json.load(response)


def _get(url: str) -> dict:
    with urllib.request.urlopen(url) as response:
        return json.load(response)


@pytest.fixture
def home_gui(tmp_path, monkeypatch):
    """`fastmdx gui` typed in a folder of studies, with no study open,
    from a process whose own folder is elsewhere."""
    from fastmdxplora.gui import preview
    from fastmdxplora.gui.server import start_dashboard_session

    monkeypatch.setattr(preview, "_platform", lambda requested, precision: ("CPU", precision))
    studies = tmp_path / "lab" / "studies"
    studies.mkdir(parents=True)
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    session = start_dashboard_session(output=str(studies), host="127.0.0.1", port=0,
                                      home_mode=True)
    try:
        yield session, studies, f"http://127.0.0.1:{session.port}"
    finally:
        session.stop()


class TestTheFolderItWasStartedIn:

    def test_new_studies_go_in_it(self, home_gui) -> None:
        session, studies, url = home_gui
        assert Path(session.runtime.exploration_root) == studies.resolve()
        assert Path(_get(url + "/api/schema")["workspace"]) == studies.resolve()

    def test_its_defaults_are_the_form_s(self, home_gui) -> None:
        _, studies, url = home_gui
        (studies / "fastmdx-defaults.yml").write_text(
            "simulation:\n  temperature_K: 310\n", encoding="utf-8")
        schema = _get(url + "/api/schema")
        fields = {f["name"]: f for f in schema["phases"]["simulation"]["fields"]}
        assert fields["temperature_K"]["default"] == 310
        assert fields["temperature_K"].get("default_from")

    def test_a_structure_beside_it_is_previewed_as_run_will_read_it(self, home_gui) -> None:
        _, studies, url = home_gui
        (studies / "only_here.pdb").write_text(_peptide(), encoding="utf-8")
        answer = _post(url + "/api/preview-system",
                       {"system": "only_here.pdb", "include_phase": ["setup", "simulation"]})
        assert answer["ok"], answer
        assert answer["estimate"]["residues"] == 3


def _yours(folder: Path, text: str):
    from fastmdxplora.config.defaults_file import DEFAULTS_FILE, read_defaults

    path = folder / DEFAULTS_FILE
    path.write_text(text, encoding="utf-8")
    return read_defaults(path)


_LAB = ("simulation:\n  temperature_K: 310\n"
        "decisions:\n  simulation.temperature_K:\n    why: Our assays run at body temperature.\n")


def _state(**simulation) -> dict:
    state = {"system": "1L2Y", "include_phase": ["setup", "simulation"]}
    if simulation:
        state["simulation"] = simulation
    return state


class TestWithYourDefaults:

    def test_fastmdxplora_s_own_value_chosen_over_yours_is_kept(self, tmp_path) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        built = config_yaml(_state(temperature_K="300"), defaults=_yours(tmp_path, _LAB))
        assert built["ok"], built
        written = yaml.safe_load(built["yaml"])
        assert written["simulation"]["temperature_K"] == 300

    def test_a_value_left_at_yours_is_left_to_your_file(self, tmp_path) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        built = config_yaml(_state(temperature_K="310"), defaults=_yours(tmp_path, _LAB))
        assert "temperature_K" not in (yaml.safe_load(built["yaml"]).get("simulation") or {})

    def test_the_plan_says_the_value_the_run_takes_marked_as_yours(self, tmp_path) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        built = config_yaml(_state(), defaults=_yours(tmp_path, _LAB))
        line = next(x for x in built["plan"] if x["label"] == "Conditions")
        assert line["value"].startswith("310 K"), line
        assert line.get("source") == "yours"

    def test_the_full_config_names_yours_and_where_it_came_from(self, tmp_path) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        built = config_yaml(_state(), full=True, defaults=_yours(tmp_path, _LAB))
        written = yaml.safe_load(built["yaml"])
        assert written["simulation"]["temperature_K"] == 310
        decision = written["decisions"]["simulation.temperature_K"]
        assert decision["source"] == "fastmdx-defaults.yml"
        assert decision["why"] == "Our assays run at body temperature."

    def test_the_builder_s_route_reads_the_workspace_s_defaults(self, home_gui) -> None:
        _, studies, url = home_gui
        (studies / "fastmdx-defaults.yml").write_text(_LAB, encoding="utf-8")
        built = _post(url + "/api/config", _state())
        line = next(x for x in built["plan"] if x["label"] == "Conditions")
        assert line["value"].startswith("310 K"), line
        built = _post(url + "/api/config", _state(temperature_K="300"))
        assert yaml.safe_load(built["yaml"])["simulation"]["temperature_K"] == 300


class TestFullConfigTicked:

    def test_run_writes_the_full_file_the_preview_showed(self, tmp_path) -> None:
        from fastmdxplora.gui.config_builder import config_yaml
        from fastmdxplora.gui.run_from_config import prepare_run

        state = {**_state(duration_ns="0.002"), "full": True}
        prepared = prepare_run(state, tmp_path / "run")
        assert prepared["ok"], prepared
        written = yaml.safe_load(Path(prepared["config_path"]).read_text(encoding="utf-8"))
        preview = yaml.safe_load(config_yaml(state, full=True)["yaml"])
        assert written == preview
        assert written["simulation"]["timestep_fs"] == 2.0

    def test_a_run_given_your_defaults_writes_them_into_the_full_file(self, tmp_path) -> None:
        from fastmdxplora.gui.run_from_config import prepare_run

        yours = _yours(tmp_path, _LAB)
        prepared = prepare_run({**_state(), "full": True}, tmp_path / "run", defaults=yours)
        written = yaml.safe_load(Path(prepared["config_path"]).read_text(encoding="utf-8"))
        assert written["simulation"]["temperature_K"] == 310


class TestCheckIt:

    def test_a_config_with_no_systems_does_not_run(self, home_gui) -> None:
        _, studies, url = home_gui
        path = studies / "no_systems.yml"
        path.write_text("simulation:\n  duration_ns: 0.001\n", encoding="utf-8")
        answer = _post(url + "/api/check-config", {"path": str(path)})
        assert answer["ok"] is False
        assert answer["refusal"]["code"] == "config.option.missing_companion"
        # Opened in the form, it is a start to fill in, not refused.
        from fastmdxplora.gui.config_builder import load_config_into_state

        assert load_config_into_state(str(path))["ok"]
        # Given a system, it runs.
        path.write_text("systems:\n  - system: 1L2Y\n", encoding="utf-8")
        answer = _post(url + "/api/check-config", {"path": str(path)})
        assert answer["ok"] and answer["systems"] == 1, answer

    def test_a_stopping_rule_on_an_analysis_that_is_not_one_does_not_run(self, home_gui) -> None:
        _, studies, url = home_gui
        path = studies / "stop.yml"
        path.write_text(
            "systems:\n  - system: 1L2Y\n"
            "simulation:\n  stop_when:\n    measures: [{analysis: rmsdd, standard_error: 0.01}]\n",
            encoding="utf-8")
        answer = _post(url + "/api/check-config", {"path": str(path)})
        assert answer["ok"] is False and "rmsdd" in answer["error"], answer


class TestRunAsItStands:

    def test_it_goes_to_the_config_s_own_output(self, home_gui, monkeypatch) -> None:
        session, studies, url = home_gui
        runtime = session.runtime
        spawned: list = []

        def spawn(command, output_dir, dashboard_url, *, config_path=None):
            spawned.append((command, Path(output_dir)))
            return {"ok": True}

        monkeypatch.setattr(runtime, "_spawn", spawn)
        path = studies / "asis.yml"
        path.write_text("systems:\n  - system: 1L2Y\noutput: out\n", encoding="utf-8")
        answer = _post(url + "/api/run-config", {"path": str(path)})
        assert answer["ok"], answer
        command, output_dir = spawned[0]
        assert output_dir == (studies / "out").resolve()
        assert command[command.index("--output") + 1] == str((studies / "out").resolve())

    def test_without_one_it_goes_beside_the_config(self, home_gui) -> None:
        session, studies, _ = home_gui
        path = studies / "plain.yml"
        path.write_text("systems:\n  - system: 1L2Y\n", encoding="utf-8")
        assert session.runtime.output_of_existing_config(str(path)) == (
            studies / "plain_output").resolve()

    def test_a_folder_asked_for_still_wins(self, home_gui) -> None:
        session, studies, _ = home_gui
        path = studies / "asis.yml"
        path.write_text("systems:\n  - system: 1L2Y\noutput: out\n", encoding="utf-8")
        assert session.runtime.output_of_existing_config(str(path), "elsewhere") == (
            studies / "elsewhere").resolve()
