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

    def test_a_value_typed_at_yours_is_kept_so_the_file_runs_alone(self, tmp_path) -> None:
        """Download config and Copy the command left out 310 K typed where
        your defaults say 310 K: taken to a machine with no defaults file,
        the file ran at 300 K. The file names what the run takes from your
        defaults; only FastMDXplora's own value, where yours say nothing
        else, goes without saying."""
        from fastmdxplora.gui.config_builder import config_yaml

        lab = _LAB.replace("simulation:\n", "setup:\n  ph: 7.0\nsimulation:\n")
        built = config_yaml({**_state(temperature_K="310", timestep_fs="2"),
                             "setup": {"ph": "7.0", "box_shape": "cube"}},
                            defaults=_yours(tmp_path, lab))
        assert built["ok"], built
        written = yaml.safe_load(built["yaml"])
        assert written["simulation"] == {"temperature_K": 310}, written
        assert written["setup"]["ph"] == 7.0
        assert "--simulate-temperature-K 310" in built["command"], built["command"]
        # Left at yours, not typed: still the value the file gives, with
        # where it came from.
        built = config_yaml(_state(), defaults=_yours(tmp_path, _LAB))
        written = yaml.safe_load(built["yaml"])
        assert written["simulation"] == {"temperature_K": 310}, written
        assert written["decisions"]["simulation.temperature_K"]["source"] == "fastmdx-defaults.yml"

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


_TWO = {"systems": [{"id": "plain", "system": "1L2Y"},
                    {"id": "warm", "system": "1L2Y", "simulation": {"temperature_K": 320}}],
        "include_phase": ["setup", "simulation"]}


class TestTwoSystemsWithYourDefaults:
    """Your defaults 310 K, one system at 320 K of its own and one with none:
    the runs are at 310 K (plain) and 320 K (warm)."""

    def test_the_plan_names_each_system_s_temperature(self, tmp_path) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        built = config_yaml(dict(_TWO), defaults=_yours(tmp_path, _LAB))
        assert built["ok"], built
        line = next(x for x in built["plan"] if x["label"] == "Conditions")
        assert line["value"] == "310 K, 1 bar (plain); 320 K, 1 bar (warm)", line
        assert "300 K" not in json.dumps(built["plan"])

    def test_the_full_config_gives_yours_only_to_the_system_that_takes_it(self, tmp_path) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        built = config_yaml(dict(_TWO), full=True, defaults=_yours(tmp_path, _LAB))
        assert built["ok"], built
        written = yaml.safe_load(built["yaml"])
        plain, warm = written["systems"]
        assert plain["simulation"]["temperature_K"] == 310
        assert warm["simulation"]["temperature_K"] == 320
        # The study's own block is not your 310 K: no run would take it
        # there, and warm's record would name your reason for its 320 K.
        assert written["simulation"]["temperature_K"] == 300
        why = written["decisions"]["simulation.temperature_K"]["why"]
        assert "For plain" in why and "warm" not in why, why

    def test_the_full_config_names_no_pressure_beside_the_one_a_system_gives(
            self, tmp_path) -> None:
        """A full config wrote the default `pressure_bar: 1.0` beside a
        system's `pressure_atm: 2`, and bar is the one the run takes."""
        from fastmdxplora.gui.config_builder import config_yaml

        state = {"systems": [{"id": "a", "system": "1L2Y",
                              "simulation": {"pressure_atm": 2}}],
                 "include_phase": ["setup", "simulation"]}
        built = config_yaml(state, full=True)
        assert built["ok"], built
        written = yaml.safe_load(built["yaml"])
        assert "pressure_bar" not in written["simulation"], written["simulation"]
        line = next(x for x in built["plan"] if x["label"] == "Conditions")
        assert "2 atm" in line["value"], line

    def test_phase_lists_given_twice_and_differing_are_a_refusal_not_a_failure(
            self, home_gui) -> None:
        _, _, url = home_gui
        built = _post(url + "/api/config", {**_state(), "include": ["setup"],
                                              "include_phase": ["setup", "simulation"]})
        assert built["ok"] is False
        assert built["refusal"]["code"] == "config.option.conflicting", built


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


class TestThePressure:

    def test_bar_given_beside_atm_is_the_pressure_said(self) -> None:
        from fastmdxplora.gui.plan import plan_of

        plan = plan_of({"systems": [{"system": "1L2Y"}],
                        "simulation": {"pressure_bar": 1.5, "pressure_atm": 2}})
        line = next(x for x in plan if x["label"] == "Conditions")
        assert line["value"] == "300 K, 1.5 bar", line
        plan = plan_of({"systems": [{"system": "1L2Y"}], "simulation": {"pressure_atm": 2}})
        assert next(x for x in plan if x["label"] == "Conditions")["value"] == "300 K, 2 atm"


def _no_pdbfixer(monkeypatch) -> None:
    from fastmdxplora.dependencies import CHEMISTRY_DEPENDENCIES

    pdbfixer = next(d for d in CHEMISTRY_DEPENDENCIES if d.module == "pdbfixer")
    for where in ("fastmdxplora.gui.exploration.missing_dependencies",
                  "fastmdxplora.dependencies.missing_dependencies"):
        monkeypatch.setattr(where, lambda include_analysis: [pdbfixer])


class TestWithoutPDBFixer:

    def test_a_structure_already_fixed_runs(self, monkeypatch) -> None:
        from fastmdxplora.gui.exploration import exploration_environment_error

        _no_pdbfixer(monkeypatch)
        assert exploration_environment_error(
            {"systems": [{"system": "a.pdb"}], "setup": {"fixed_pdb": "a_fixed.pdb"}}) is None
        assert exploration_environment_error(
            {"systems": [{"system": "a.pdb", "setup": {"fixed_pdb": "a_fixed.pdb"}},
                         {"system": "b.pdb", "setup": {"fixed_pdb": "b_fixed.pdb"}}]}) is None
        # One system still to repair needs it.
        assert "PDBFixer" in exploration_environment_error(
            {"systems": [{"system": "a.pdb", "setup": {"fixed_pdb": "a_fixed.pdb"}},
                         {"system": "b.pdb"}]})

    def test_a_refused_run_leaves_its_folder_free(self, tmp_path, monkeypatch) -> None:
        from fastmdxplora.gui.exploration import DashboardRuntime

        _no_pdbfixer(monkeypatch)
        runtime = DashboardRuntime(workspace_root=tmp_path / "workspace",
                                   exploration_root=tmp_path / "runs")
        state = {**_state(), "output": "trial"}
        first = runtime.launch_from_config(state)
        assert first["ok"] is False and "PDBFixer" in first["error"]
        assert not (tmp_path / "runs" / "trial").exists()
        assert first["refusal"]["code"] == "environment.backend.missing"
        second = runtime.launch_from_config(state)
        assert "not empty" not in second["error"] and "PDBFixer" in second["error"]


# ---------------------------------------------------------------- in a browser

@pytest.fixture(scope="module")
def browser():
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


def _builder(browser, url: str, width: int = 1440, scheme: str = "light",
             system: str | None = "1L2Y"):
    context = browser.new_context(viewport={"width": width, "height": 900 if width > 600 else 844})
    context.add_init_script(f"localStorage.setItem('fmx.theme', '{scheme}');")
    page = context.new_page()
    page.set_default_timeout(60000)
    errors: list[str] = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.goto(url + "/#run", wait_until="domcontentloaded")
    page.wait_for_selector("#run-start [data-start='structure']", state="attached")
    page.evaluate("() => window.FastMDXRun.setStart('structure')")
    if system:
        _system(page, system)
    page.errors = errors
    return page


def _placeholders_cut(page) -> list[dict]:
    """The boxes whose placeholder is drawn short of its end, read from the
    pixels as a person sees them. Measuring the text against the box's width
    passed while Chromium kept the right of a number box for its step
    arrows, and the last word was cut there; and fields in a closed section
    were never measured."""
    import io

    from PIL import Image

    page.evaluate("() => document.querySelectorAll('#run-settings details')"
                  ".forEach(d => { d.open = true; })")
    page.evaluate("() => document.querySelectorAll('#run-settings [aria-expanded=false]')"
                  ".forEach(b => b.click())")
    page.wait_for_timeout(500)
    cut = []
    for box in page.query_selector_all("#run-settings input[placeholder]"):
        said = page.evaluate("""box => {
            if (!box.placeholder || !box.offsetParent || box.value) return null;
            const style = getComputedStyle(box);
            const canvas = document.createElement('canvas').getContext('2d');
            canvas.font = style.font;
            return {id: box.id, text: box.placeholder,
                    ends: parseFloat(style.paddingLeft) + parseFloat(style.borderLeftWidth)
                          + canvas.measureText(box.placeholder).width};
        }""", box)
        if not said:
            continue
        box.scroll_into_view_if_needed()
        shot = Image.open(io.BytesIO(box.screenshot())).convert("L")
        width, height = shot.size
        ground = shot.getpixel((width - 3, height // 2))
        inked = [x for x in range(width) for y in range(4, height - 4)
                 if abs(shot.getpixel((x, y)) - ground) > 40]
        drawn_to = max(inked) if inked else 0
        if drawn_to < said["ends"] - 4:
            cut.append({**said, "drawn_to": drawn_to})
    return cut


def _system(page, system: str) -> None:
    page.fill("#run-system", system)
    page.evaluate("() => document.getElementById('run-system')"
                  ".dispatchEvent(new Event('change', {bubbles: true}))")


def _type(page, selector: str, value: str) -> None:
    page.evaluate("""([selector, value]) => {
        const box = document.querySelector(selector);
        box.value = value;
        box.dispatchEvent(new Event('input', {bubbles: true}));
        box.dispatchEvent(new Event('change', {bubbles: true}));
    }""", [selector, value])


class TestInTheBuilder:

    def test_a_length_changed_says_what_it_was(self, browser, home_gui) -> None:
        _, _, url = home_gui
        page = _builder(browser, url)
        _type(page, "#set-simulation-duration_ns", "0.002")
        page.wait_for_function(
            "() => document.getElementById('run-changed').textContent.includes('0.002')")
        said = page.text_content("#run-changed")
        assert "2 ns → 0.002 ns" in said, said
        assert "not set" not in said
        assert not page.errors
        page.context.close()

    @pytest.mark.parametrize("width", [1440, 390])
    @pytest.mark.parametrize("scheme", ["light", "dark"])
    def test_the_verdict_and_the_empty_boxes_are_read_whole(
            self, browser, home_gui, width, scheme) -> None:
        _, _, url = home_gui
        page = _builder(browser, url, width, scheme, system=None)

        def cut() -> list:
            page.wait_for_timeout(300)
            shown = page.evaluate("""() => ['run-status', 'run-phone-status'].map(id => {
                const chip = document.getElementById(id);
                if (!chip || !chip.offsetParent) return null;
                const box = chip.getBoundingClientRect();
                return {id, text: chip.textContent, cut: chip.scrollWidth > chip.clientWidth + 1,
                        off: box.right > document.documentElement.clientWidth + 1};
            }).filter(Boolean)""")
            assert shown, "a verdict is shown"
            return [chip for chip in shown if chip["cut"] or chip["off"]]

        # The two longest verdicts the form says: what to give first, and
        # checks passed with something worth knowing (a 4 fs step with no
        # constraints).
        page.wait_for_function("() => document.getElementById('run-status')"
                               ".textContent.startsWith('Choose a structure')")
        assert not cut()
        page.evaluate("""() => {
            const s = window.FastMDXRun.state;
            s.values.simulation = Object.assign({}, s.values.simulation,
                {timestep_fs: 4, constraints: 'None'});
        }""")
        _system(page, "1L2Y")
        _type(page, "#set-simulation-duration_ns", "0.002")
        page.wait_for_function("() => document.getElementById('run-status')"
                               ".textContent.includes('worth knowing')")
        assert not cut()
        # Each box, empty, says its whole placeholder, every section open.
        short = _placeholders_cut(page)
        assert not short, short
        assert page.evaluate(
            "() => document.documentElement.scrollWidth <= document.documentElement.clientWidth")
        assert not page.errors
        page.context.close()

    def test_run_as_it_stands_sends_no_results_folder(self, browser, home_gui, tmp_path) -> None:
        _, studies, url = home_gui
        path = studies / "asis.yml"
        path.write_text("systems:\n  - system: 1L2Y\noutput: out\n", encoding="utf-8")
        page = _builder(browser, url)
        _type(page, "#run-output", "stale_box")
        sent: list[dict] = []

        def run_config(route):
            sent.append(json.loads(route.request.post_data or "{}"))
            route.fulfill(status=200, content_type="application/json",
                          body=json.dumps({"ok": False, "error": "held by the test"}))

        page.route("**/api/run-config", run_config)
        page.evaluate("path => { window.FastMDXRun.setStart('config'); "
                      "const box = document.getElementById('run-config-path'); "
                      "box.value = path; box.dispatchEvent(new Event('change')); }", str(path))
        page.wait_for_function("() => !document.getElementById('run-as-is').disabled")
        page.click("#run-as-is")
        page.wait_for_function("() => document.getElementById('run-note').textContent"
                               ".includes('held by the test')")
        assert sent == [{"path": str(path)}], sent
        page.context.close()

    def test_a_timestep_of_zero_is_refused_on_its_box(self, browser, home_gui) -> None:
        """The schema refuses a timestep of 0; the box let it through to the server."""
        _, _, url = home_gui
        page = _builder(browser, url)
        _type(page, "#set-simulation-timestep_fs", "0")
        page.wait_for_function(
            "() => !document.querySelector('#set-simulation-timestep_fs')"
            ".closest('.builder-field, label, div').querySelector('.builder-limit').hidden")
        said = page.evaluate(
            "() => document.querySelector('#set-simulation-timestep_fs')"
            ".closest('.builder-field, label, div').querySelector('.builder-limit').textContent")
        assert "must be more than 0" in said, said
        _type(page, "#set-simulation-timestep_fs", "2")
        page.wait_for_function(
            "() => document.querySelector('#set-simulation-timestep_fs')"
            ".closest('.builder-field, label, div').querySelector('.builder-limit').hidden")
        assert not page.errors
        page.context.close()

    @pytest.mark.parametrize("scheme", ["light", "dark"])
    def test_a_long_verdict_at_the_foot_of_a_phone_is_three_lines_and_leads_to_the_rest(
            self, browser, home_gui, scheme) -> None:
        """A long refusal grew the bar at a phone's foot to half the screen."""
        _, _, url = home_gui
        page = _builder(browser, url, 390, scheme)
        _type(page, "#set-simulation-temperature_K", "-5")
        page.wait_for_function(
            "() => document.getElementById('run-phone-status').textContent.includes('at least 0')")
        # The longest a refusal says: one naming every setting it could be.
        page.evaluate("""() => { const chip = document.getElementById('run-phone-status');
            chip.textContent += ' Valid options: '
              + Array.from({length: 45}, (_, i) => 'option_' + i).join(', '); }""")
        held = page.evaluate("""() => {
            const chip = document.getElementById('run-phone-status');
            const bar = document.getElementById('run-phone-bar');
            return {chip: chip.getBoundingClientRect().height,
                    line: parseFloat(getComputedStyle(chip).lineHeight),
                    bar: bar.getBoundingClientRect().height,
                    run: document.getElementById('run-phone-start').getBoundingClientRect().height};
        }""")
        assert held["chip"] <= 3 * held["line"] + 12, held
        assert held["bar"] < 844 / 6, held
        assert held["run"] < 50, held
        page.click("#run-phone-status")
        top = page.evaluate("() => document.getElementById('run-status').getBoundingClientRect().top")
        assert 0 <= top < 844, top
        assert not page.errors
        page.context.close()

    def test_a_study_opened_with_workers_says_the_mode_it_runs_in(
            self, browser, home_gui) -> None:
        """Opened with `execution: {workers: 2}`, the Mode box said "sequential
        (default)" of a study that runs in parallel."""
        _, studies, url = home_gui
        path = studies / "two.yml"
        path.write_text("systems:\n  - system: 1L2Y\n  - system: 1UBQ\n"
                        "execution:\n  workers: 2\n", encoding="utf-8")
        loaded = _post(url + "/api/load-config", {"path": str(path)})
        page = _builder(browser, url)
        page.evaluate("from => window.FastMDXRun.applyLoadedState(from)", loaded.get("state", loaded))
        page.wait_for_selector("#set-execution-mode", state="attached")
        box = page.evaluate("""() => { const box = document.getElementById('set-execution-mode');
            return {value: box.value, said: box.options[box.selectedIndex].textContent}; }""")
        assert box["value"] == "parallel", box
        assert "default" not in box["said"], box
        assert not page.errors
        page.context.close()


class TestCheckItReadsWhereTheRunReads:
    """Check it read a relative structure from the process's own folder,
    and Run as it stands from the folder the GUI was started in."""

    def test_a_relative_structure_is_read_from_the_studies(self, home_gui) -> None:
        _, studies, url = home_gui
        Path("only_there.pdb").write_text(_peptide(), encoding="utf-8")
        path = studies / "rel2.yml"
        path.write_text("systems:\n  - system: only_there.pdb\n", encoding="utf-8")
        answer = _post(url + "/api/check-config", {"path": str(path)})
        assert answer["ok"] is False, answer
        assert answer["refusal"]["code"] == "environment.path.not_found"
        # One beside the studies, which the run reads, is found, the
        # Config being kept in a folder of its own.
        (studies / "only_here.pdb").write_text(_peptide(), encoding="utf-8")
        (studies / "configs").mkdir()
        path = studies / "configs" / "rel.yml"
        path.write_text("systems:\n  - system: only_here.pdb\n", encoding="utf-8")
        answer = _post(url + "/api/check-config", {"path": str(path)})
        assert answer["ok"], answer

    def test_a_config_with_no_systems_is_told_what_to_add(self, home_gui) -> None:
        _, studies, url = home_gui
        path = studies / "none.yml"
        path.write_text("simulation:\n  duration_ns: 0.001\n", encoding="utf-8")
        answer = _post(url + "/api/check-config", {"path": str(path)})
        assert answer["ok"] is False
        assert "Add one under `systems:`" in answer["error"], answer["error"]
        assert "not invented" not in answer["error"]


class TestTheThirdRound:
    """What the second review of these fixes found."""

    def test_a_full_config_of_a_swept_atm_runs(self, tmp_path) -> None:
        """Full config wrote FastMDXplora's `pressure_bar` beside a swept
        `pressure_atm`, and the study was refused as giving both."""
        from fastmdxplora.gui.config_builder import config_yaml

        state = {**_state(), "full": True,
                 "sweep": [{"axis": "simulation.pressure_atm", "values": "1, 2"}]}
        made = config_yaml(state, full=True)
        assert made["ok"], made
        written = yaml.safe_load(made["yaml"])
        assert "pressure_bar" not in (written.get("simulation") or {})

    def test_your_defaults_written_in_are_said_and_kept_out_of_the_command(
            self, tmp_path) -> None:
        from fastmdxplora.gui.config_builder import config_yaml

        yours = _yours(tmp_path, "simulation:\n  temperature_K: 310\n"
                                 "decisions:\n  simulation.temperature_K:\n"
                                 "    why: Our assays run at body temperature.\n")
        made = config_yaml(_state(), defaults=yours)
        assert made["ok"], made
        assert "fastmdx-defaults.yml gives written in" in made["yaml"]
        assert "--decisions" not in (made.get("command") or ""), made.get("command")
        assert "--temperature" in (made.get("command") or "") or "310" in (made.get("command") or "")

    def test_the_mode_box_says_parallel_for_a_study_opened_with_workers(
            self, browser, home_gui) -> None:
        _, studies, url = home_gui
        path = studies / "par.yml"
        path.write_text("systems:\n  - system: 1L2Y\nexecution:\n  workers: 2\n",
                        encoding="utf-8")
        # Opened as Open for editing opens it, on a page that has drawn
        # nothing else first: the review's boxes are drawn again after the
        # study arrives, and the new Mode box said "sequential (default)".
        context = browser.new_context(viewport={"width": 1440, "height": 900})
        page = context.new_page()
        page.errors = []
        page.on("pageerror", lambda error: page.errors.append(str(error)))
        page.goto(url + "/#run", wait_until="domcontentloaded")
        page.wait_for_selector("#run-start [data-start='structure']", state="attached")
        loaded = _post(url + "/api/load-config", {"path": str(path)})
        page.evaluate("from => window.FastMDXRun.applyLoadedState(from)",
                      loaded.get("state", loaded))
        page.wait_for_timeout(1500)
        said = page.evaluate("() => { const box = document.getElementById('set-execution-mode');"
                             " return [box.value, box.options[box.selectedIndex].text]; }")
        assert said[0] == "parallel", said
        assert "default" not in said[1], said
        assert not page.errors
        page.context.close()
