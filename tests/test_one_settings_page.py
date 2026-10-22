"""Every setting on one Settings page (IMDE 2.9, user 10-21).

Settings were four places: the gear's popup (theme, accent, the Agent),
the Preferences dialog (the Viewer, time, charts, notices), the Agent's
own dialog, and the workspace's `fastmdx-defaults.yml`, which the GUI read
for new studies and could not change. They are one page, found by their
words, with two scopes: You, saved in this browser, and This folder, the
defaults file every new study made here starts from. Each applies as it
is changed. The Display size (the ideas session's large-screen entry,
Auto by default, user 10-21) is one of You's: the whole page at once,
text included, Auto from the screen.
"""

from __future__ import annotations

import pytest
import yaml

from tests.test_the_workspace_says_its_studies import _study


# --------------------------------------------------------------------------
# The defaults file, written
# --------------------------------------------------------------------------


def test_a_defaults_file_is_written_checked_and_taken_out_again(tmp_path):
    from fastmdxplora.config.defaults_file import read_defaults, write_defaults

    where = tmp_path / "fastmdx-defaults.yml"
    write_defaults(where, given={"simulation.temperature_K": 310, "analysis.include": ["rmsd"]},
                   whys={"simulation.temperature_K": "body temperature"})
    assert yaml.safe_load(where.read_text(encoding="utf-8")) == {
        "simulation": {"temperature_K": 310},
        "analysis": {"include": ["rmsd"]},
        "decisions": {"simulation.temperature_K": {"why": "body temperature"}},
    }
    # Another setting leaves the others as they are.
    write_defaults(where, given={"setup.ph": 7.0})
    found = read_defaults(where)
    assert dict(found.settings()) == {"setup.ph": 7.0, "simulation.temperature_K": 310,
                                      "analysis.include": ["rmsd"]}
    # Taken out, its reason goes with it; an empty reason is none.
    write_defaults(where, unset=["simulation.temperature_K"], whys={"setup.ph": ""})
    found = read_defaults(where)
    assert dict(found.settings()) == {"setup.ph": 7.0, "analysis.include": ["rmsd"]}
    assert dict(found.whys) == {}
    # Nothing is left beside it, and it is readable as it was.
    assert sorted(p.name for p in tmp_path.iterdir()) == ["fastmdx-defaults.yml"]
    where.chmod(0o644)
    write_defaults(where, given={"setup.ph": 7.4})
    assert where.stat().st_mode & 0o777 == 0o644


@pytest.mark.parametrize("given", [
    {"simulation.temperature_K": "hot"},
    {"simulation.timestep_fs": -2},
    {"simulation.no_such_setting": 1},
    {"systems.x": 1},
    {"simulation.umbrella.windows": 3},
])
def test_a_value_any_defaults_file_refuses_leaves_the_file_as_it_was(tmp_path, given):
    from fastmdxplora.config.defaults_file import write_defaults
    from fastmdxplora.config.loader import ConfigError

    where = tmp_path / "fastmdx-defaults.yml"
    write_defaults(where, given={"setup.ph": 7.0})
    before = where.read_text(encoding="utf-8")
    with pytest.raises(ConfigError) as refused:
        write_defaults(where, given=given)
    assert where.read_text(encoding="utf-8") == before
    # Said of the file, never of the file written beside it to check.
    assert ".fastmdx-defaults." not in str(refused.value)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["fastmdx-defaults.yml"]


def test_a_file_wrong_already_is_not_written_over(tmp_path):
    from fastmdxplora.config.defaults_file import write_defaults
    from fastmdxplora.config.loader import ConfigError

    where = tmp_path / "fastmdx-defaults.yml"
    wrong = "setup:\n  ph: acid\n"
    where.write_text(wrong, encoding="utf-8")
    with pytest.raises(ConfigError):
        write_defaults(where, given={"simulation.temperature_K": 310})
    assert where.read_text(encoding="utf-8") == wrong


def test_a_file_with_comments_is_changed_by_hand(tmp_path):
    """Written again from its values, a file would lose its comments: Settings
    says so and leaves it, and shows it as it is."""
    from fastmdxplora.gui.settings_page import defaults_changed, defaults_said, has_comments

    commented = "# the lab's\nsetup:\n  ph: 7.0  # as the buffer\n"
    where = tmp_path / "fastmdx-defaults.yml"
    where.write_text(commented, encoding="utf-8")
    said = defaults_said(tmp_path)
    assert said["comments"] and said["values"] == {"setup.ph": 7.0}
    refused = defaults_changed(tmp_path, {"set": {"setup.ph": 6.5}})
    assert refused["ok"] is False and "comments" in refused["reason"]
    assert where.read_text(encoding="utf-8") == commented
    # A colour is no comment.
    assert not has_comments("report:\n  region_highlights:\n  - {label: a, start: 1, end: 3, "
                            "color: '#4E79A7'}\n")


# --------------------------------------------------------------------------
# This folder, as Settings says it
# --------------------------------------------------------------------------


def test_this_folder_says_every_setting_with_fastmdxplora_s_own(tmp_path):
    from fastmdxplora.gui.settings_page import defaults_said

    said = defaults_said(tmp_path)
    assert said["ok"] and not said["exists"]
    assert said["path"] == str(tmp_path / "fastmdx-defaults.yml")
    fields = {f["name"]: f for f in said["fields"]}
    assert {f["block"] for f in said["fields"]} == {"setup", "simulation", "analysis", "report",
                                                    "execution"}
    assert fields["simulation.temperature_K"]["default"] == 300
    assert fields["simulation.temperature_K"]["control"] == "number"
    # A list of numbers is written from numbers between commas; one of
    # mappings is not drawn, and is changed in the file by hand.
    assert fields["simulation.restraint_release"]["items"] == "number"
    assert fields["setup.force_field"]["items"] == "text"
    assert fields["report.region_highlights"]["control"] == ""
    # Who chose a phase's settings is the study's record, and a setting of
    # one study (what it continues or is prepared from, the files it
    # analyses, what it simulates) is never a folder's: not offered.
    for one_study_s in ("setup.agent", "simulation.resume_from", "simulation.resume_unsealed",
                        "simulation.setup_from", "simulation.prepared_from",
                        "simulation.extra_ns", "analysis.trajectory", "analysis.topology",
                        "setup.fixed_pdb", "setup.ligand", "setup.ligand_pose",
                        "setup.mutations", "setup.chains", "setup.model"):
        assert one_study_s not in fields, one_study_s


def test_a_file_setting_what_is_simulated_is_refused_and_said(tmp_path):
    """`L99A` in the folder's file would mutate every new study's protein."""
    from fastmdxplora.gui.settings_page import defaults_changed, defaults_said

    where = tmp_path / "fastmdx-defaults.yml"
    where.write_text("setup:\n  mutations: [L99A]\n", encoding="utf-8")
    said = defaults_said(tmp_path)
    # Said once, under the file's own name: not its path again.
    assert said["refused"].startswith("It sets `setup.mutations`, which says what one study "
                                      "simulates"), said["refused"]
    assert str(where) not in said["refused"] and said["values"] == {}
    refused = defaults_changed(tmp_path, {"set": {"setup.ph": 7.0}})
    assert refused["ok"] is False and str(where) not in refused["reason"]
    assert where.read_text(encoding="utf-8") == "setup:\n  mutations: [L99A]\n"
    # Asked of a file that is right, it is refused the same, the file as it was.
    where.write_text("setup:\n  ph: 7.0\n", encoding="utf-8")
    refused = defaults_changed(tmp_path, {"set": {"simulation.resume_from": "runs/old"}})
    assert refused["ok"] is False
    assert refused["reason"].startswith("It sets `simulation.resume_from`, which names the study")
    assert str(where) not in refused["reason"]
    assert where.read_text(encoding="utf-8") == "setup:\n  ph: 7.0\n"


def test_a_setting_of_several_values_the_file_holds_is_shown(tmp_path):
    """A stopping rule in the file makes every new study run until it is
    determined, up to its ceiling: shown as the file holds it, and counted."""
    from fastmdxplora.gui.settings_page import defaults_said

    (tmp_path / "fastmdx-defaults.yml").write_text(
        "simulation:\n  temperature_K: 310\n  stop_when:\n    measures:\n"
        "    - analysis: rmsd\n      standard_error: 0.01\n    max_duration_ns: 50\n"
        "setup:\n  dispersion_correction: null\n", encoding="utf-8")
    said = defaults_said(tmp_path)
    fields = {f["name"]: f for f in said["fields"]}
    assert set(said["values"]) <= set(fields)
    for several in ("simulation.stop_when", "simulation.umbrella", "simulation.metadynamics",
                    "setup.residue_states", "analysis.options"):
        assert fields[several]["control"] == "", several
    assert said["values"]["simulation.stop_when"]["max_duration_ns"] == 50


def test_two_changes_at_once_are_both_written(tmp_path, monkeypatch):
    """Each request is answered in its own thread: two changes that read
    the file before either was written lost the first."""
    import os
    import threading
    import time

    from fastmdxplora.config.defaults_file import read_defaults
    from fastmdxplora.gui.settings_page import defaults_changed

    replace = os.replace

    def slowly(source, target):
        time.sleep(0.4)
        replace(source, target)

    monkeypatch.setattr(os, "replace", slowly)
    asked = [{"set": {"setup.ph": 7.0}}, {"set": {"simulation.temperature_K": 310}}]
    threads = [threading.Thread(target=defaults_changed, args=(tmp_path, body)) for body in asked]
    for thread in threads:
        thread.start()
        time.sleep(0.1)
    for thread in threads:
        thread.join()
    assert dict(read_defaults(tmp_path / "fastmdx-defaults.yml").settings()) == {
        "setup.ph": 7.0, "simulation.temperature_K": 310}


def test_the_nearest_file_above_is_the_one_shown_and_changed(tmp_path):
    from fastmdxplora.gui.settings_page import defaults_changed, defaults_said

    (tmp_path / "fastmdx-defaults.yml").write_text("setup:\n  ph: 7.0\n", encoding="utf-8")
    workspace = tmp_path / "studies"
    workspace.mkdir()
    said = defaults_said(workspace)
    assert said["exists"] and said["path"] == str(tmp_path / "fastmdx-defaults.yml")
    assert said["values"] == {"setup.ph": 7.0}
    changed = defaults_changed(workspace, {"set": {"simulation.temperature_K": 310},
                                           "why": {"simulation.temperature_K": "warm"}})
    assert changed["ok"] and changed["values"]["simulation.temperature_K"] == 310
    assert changed["whys"] == {"simulation.temperature_K": "warm"}
    assert not (workspace / "fastmdx-defaults.yml").exists()
    refused = defaults_changed(workspace, {"set": {"simulation.temperature_K": -5}})
    assert refused["ok"] is False
    assert refused["reason"].startswith("Simulation option 'temperature_K'")
    assert str(tmp_path) not in refused["reason"]
    assert refused["values"]["simulation.temperature_K"] == 310


# --------------------------------------------------------------------------
# In a browser
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def browser():
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


@pytest.fixture
def session(tmp_path):
    from fastmdxplora.gui.server import start_dashboard_session

    study = _study(tmp_path / "ubiquitin", system="1UBQ")
    started = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    yield started
    started.server.shutdown()


def _open(browser, url, where="#settings", screen=None, folder_404=False):
    options = {"viewport": {"width": 1440, "height": 900}}
    if screen:
        options["screen"] = {"width": screen[0], "height": screen[1]}
    context = browser.new_context(**options)
    page = context.new_page()
    page.set_default_timeout(60000)
    page.errors = []
    page.on("pageerror", lambda error: page.errors.append(str(error)))
    if folder_404:
        # A server without This folder's route (the GUI before it).
        page.route("**/api/defaults-file", lambda route: route.fulfill(status=404, body="{}"))
    page.goto(url + where, wait_until="domcontentloaded")
    page.wait_for_function("() => window.FastMDXDashboard && window.FastMDXSettings")
    page.wait_for_selector("body:not(.state-loading)")
    return page


def test_one_page_two_scopes_one_search(browser, session):
    page = _open(browser, session.url, folder_404=True)
    page.wait_for_selector('.page[data-page="settings"]:not([hidden])')
    shown = page.evaluate("""() => ['settings-you', 'settings-folder']
        .map((id) => !document.getElementById(id).hidden)""")
    # Every preference that was a dialog's, and the popup's theme and
    # accent, are rows of You.
    rows = page.evaluate("""() => [...document.querySelectorAll('#settings-you .pref-row > span:first-child')]
        .map((s) => s.textContent.trim())""")
    page.fill("#settings-find", "water")
    found = page.evaluate("""() => [...document.querySelectorAll('#settings-you .pref-row')]
        .filter((r) => !r.hidden).map((r) => r.textContent.trim())""")
    counted = page.evaluate("""() => [...document.querySelectorAll('.settings-count')]
        .map((c) => c.hidden ? null : c.textContent)""")
    page.fill("#settings-find", "no such words")
    none = page.text_content("#settings-none-you")
    page.keyboard.press("Escape")
    cleared = page.input_value("#settings-find")
    # The keyboard goes from scope to scope with the arrows.
    page.focus("#settings-scope-you")
    page.keyboard.press("ArrowRight")
    moved = page.evaluate("""() => [document.activeElement.id,
        document.getElementById('settings-folder').hidden, document.getElementById('settings-you').hidden]""")
    page.wait_for_selector("#settings-folder-body .settings-cannot")
    cannot = page.text_content("#settings-folder-body")
    errors = page.errors
    page.context.close()
    assert errors == []
    assert shown == [True, False]
    for name in ("Theme", "Accent", "Display size", "Representation", "Ground", "Water", "Ions",
                 "Spin", "Pocket cutoff (Å)", "Time", "Chart history", "Notify when a run ends",
                 "Engine", "Mode"):
        assert any(r.startswith(name) for r in rows), name
    assert len(found) == 1 and found[0].startswith("Water")
    # This folder's count is not said where its file could not be read.
    assert counted == ["1", None]
    assert none.startswith("No setting here has those words.")
    assert cleared == ""
    assert moved == ["settings-scope-folder", False, True]
    assert "cannot read the folder's defaults yet" in cannot


def test_the_popup_the_agent_s_dialog_and_the_palette_lead_here(browser, session):
    page = _open(browser, session.url, "#analysis", folder_404=True)
    page.click("#settings-open")
    page.click("#settings-page-link")
    page.wait_for_selector('.page[data-page="settings"]:not([hidden])')
    from_popup = (page.is_hidden("#settings-popup"), page.evaluate("() => location.hash"))
    page.evaluate("() => window.FastMDXDashboard.navigate('overview')")
    # The Agent's dialog, to its part of the page.
    page.evaluate("() => window.FastMDXAgent.openSettings()")
    page.wait_for_selector("#agent-settings:not([hidden])")
    page.click('#agent-settings [data-settings-part-link="agent"]')
    page.wait_for_selector('.page[data-page="settings"]:not([hidden])')
    page.wait_for_function("() => document.activeElement.id === 'settings-page-agent'")
    from_agent = (page.is_hidden("#agent-settings"),
                  page.text_content("#settings-page-mode") == page.text_content("#settings-mode"))
    # And back to the Agent's dialog from there.
    page.click("#settings-page-agent")
    page.wait_for_selector("#agent-settings:not([hidden])")
    page.keyboard.press("Escape")
    # What still opens Preferences by the dialog's name (the palette's
    # Preferences) is taken here.
    page.evaluate("() => window.FastMDXDashboard.navigate('overview')")
    page.evaluate("() => window.FastMDXDialog.open('prefs-dialog')")
    page.wait_for_selector('.page[data-page="settings"]:not([hidden])')
    errors = page.errors
    page.context.close()
    assert errors == []
    assert from_popup == (True, "#settings")
    assert from_agent == (True, True)


def test_the_theme_here_is_the_popup_s(browser, session):
    page = _open(browser, session.url, folder_404=True)
    page.click('#settings-you [data-theme="light"]')
    page.click('#settings-you [data-accent="violet"]')
    said = page.evaluate("""() => [document.documentElement.dataset.theme,
        document.documentElement.dataset.accent,
        document.querySelector('#settings-popup [data-theme="light"]').getAttribute('aria-pressed'),
        document.querySelector('#settings-popup [data-accent="violet"]').getAttribute('aria-pressed')]""")
    page.reload(wait_until="domcontentloaded")
    page.wait_for_selector("body:not(.state-loading)")
    again = page.evaluate("""() => [document.documentElement.dataset.theme,
        document.querySelector('#settings-you [data-accent="violet"]').getAttribute('aria-pressed')]""")
    page.context.close()
    assert said == ["light", "violet", "true", "true"]
    assert again == ["light", "true"]


@pytest.mark.parametrize("screen, auto", [((1920, 1080), 100), ((2560, 1440), 125),
                                          ((3440, 1440), 125), ((3840, 2160), 150),
                                          ((1366, 768), 100)])
def test_auto_is_chosen_from_the_screen(browser, session, screen, auto):
    page = _open(browser, session.url, screen=screen, folder_404=True)
    said = page.evaluate("""() => [window.FastMDXDisplay.chosen(), window.FastMDXDisplay.zoom(),
        document.getElementById('setting-display-size').selectedOptions[0].textContent]""")
    page.context.close()
    assert said == ["auto", auto / 100, f"Auto: {auto}%"]


def test_a_display_size_chosen_sizes_everything_and_is_saved(browser, session):
    page = _open(browser, session.url, screen=(1920, 1080), folder_404=True)
    page.select_option("#setting-display-size", "125")
    page.wait_for_function("() => window.FastMDXDisplay.zoom() === 1.25")
    sized = page.evaluate("""() => {
        const shell = document.querySelector('.app-shell');
        const title = document.querySelector('.page[data-page="settings"] .page-title');
        const box = title.getBoundingClientRect();
        return {
            zoom: getComputedStyle(document.documentElement).zoom,
            saved: localStorage.getItem('fmx.displaySize'),
            // The page is told of the screen in its own pixels: the
            // shell fills the window it is told of down to the status
            // bar, which ends at its foot, and a box is where the page
            // put it.
            fits: Math.abs(shell.getBoundingClientRect().bottom
                           - document.querySelector('.status-bar').getBoundingClientRect().top) < 2
                && Math.abs(document.querySelector('.status-bar').getBoundingClientRect().bottom
                            - innerHeight) < 2,
            told: innerWidth,
            box: Math.abs(box.width - title.offsetWidth) < 1.5,
        };
    }""")
    page.reload(wait_until="domcontentloaded")
    page.wait_for_selector("body:not(.state-loading)")
    again = page.evaluate("() => [window.FastMDXDisplay.zoom(), "
                          "document.getElementById('setting-display-size').value]")
    # Restore defaults puts Auto back.
    page.click("#prefs-restore")
    page.wait_for_function("() => window.FastMDXDisplay.zoom() === 1")
    restored = page.evaluate("() => [localStorage.getItem('fmx.displaySize'), "
                             "document.getElementById('setting-display-size').value]")
    errors = page.errors
    page.context.close()
    assert errors == []
    assert sized["zoom"] == "1.25" and sized["saved"] == "125"
    assert sized["fits"] and sized["told"] == 1440 / 1.25 and sized["box"]
    assert again == [1.25, "125"]
    assert restored == ["auto", "auto"]


def _at(browser, url, size, viewport, where="#overview"):
    context = browser.new_context(viewport={"width": viewport[0], "height": viewport[1]})
    context.add_init_script(f"try {{ localStorage.setItem('fmx.displaySize', '{size}'); }} "
                            "catch (e) {}")
    page = context.new_page()
    page.set_default_timeout(60000)
    page.errors = []
    page.on("pageerror", lambda error: page.errors.append(str(error)))
    page.route("**/api/defaults-file", lambda route: route.fulfill(status=404, body="{}"))
    page.goto(url + where, wait_until="domcontentloaded")
    page.wait_for_function("() => window.FastMDXDashboard && window.FastMDXSettings")
    page.wait_for_selector("body:not(.state-loading)")
    return page


def _box(page, selector):
    return page.evaluate(f"""() => {{ const r = document.querySelector('{selector}')
        .getBoundingClientRect(); return [r.left, r.top, r.right, r.bottom]; }}""")


@pytest.mark.parametrize("size", ["90", "110", "125", "150"])
def test_at_every_display_size_tips_and_menus_land_where_they_are_seen(browser, session, size):
    """A tooltip was placed from the window's height read from the root
    element, which said the screen's pixels: at 125% the status bar's tip
    was put wholly below the window, and one near the right edge ran past
    it. The gear's menu, opened by a click where its button is seen, opens
    over it; and the columns end at the status bar at every size."""
    page = _at(browser, session.url, size, (1440, 900))
    zoom = int(size) / 100
    said = {"window": page.evaluate("() => [innerWidth, innerHeight]")}
    shell, bar = _box(page, ".app-shell"), _box(page, ".status-bar")
    said["columns"] = abs(shell[3] - bar[1]) < 2 and abs(bar[3] - said["window"][1]) < 2
    tips = {}
    for target in ("#sb-machine", "#refresh-now"):
        at = _box(page, target)
        page.mouse.move(0, 0)
        page.wait_for_timeout(350)
        page.mouse.move((at[0] + at[2]) / 2 * zoom, (at[1] + at[3]) / 2 * zoom)
        page.wait_for_selector("#fmx-tooltip:not([hidden])")
        tip = _box(page, "#fmx-tooltip")
        tips[target] = {"inside": tip[0] >= 0 and tip[1] >= 0 and tip[2] <= said["window"][0] + 0.5
                        and tip[3] <= said["window"][1] + 0.5,
                        "beside": min(abs(at[1] - tip[3]), abs(tip[1] - at[3])) < 12}
    page.mouse.move(0, 0)
    gear = _box(page, "#settings-open")
    page.mouse.click((gear[0] + gear[2]) / 2 * zoom, (gear[1] + gear[3]) / 2 * zoom)
    page.wait_for_selector("#settings-popup:not([hidden])")
    menu = _box(page, "#settings-popup")
    errors = page.errors
    page.context.close()
    assert errors == []
    assert abs(said["window"][0] - 1440 / zoom) < 0.01, said
    assert said["columns"], said
    for target, placed in tips.items():
        assert placed == {"inside": True, "beside": True}, (target, placed)
    assert menu[0] >= 0 and menu[1] >= 0 and menu[2] <= said["window"][0] + 0.5
    assert menu[3] <= said["window"][1] + 0.5
    # Over the button it was opened from, its foot on the button.
    assert menu[0] < gear[2] and gear[0] < menu[2], (menu, gear)
    assert gear[1] - 12 < menu[3] <= gear[3] + 12, (menu, gear)


def test_the_docs_say_the_settings_page(tmp_path):
    from pathlib import Path

    text = (Path(__file__).resolve().parent.parent / "docs" / "gui.md").read_text(encoding="utf-8")
    said = text[text.index("**Settings**"):]
    assert "This folder" in said and "fastmdx-defaults.yml" in said
    assert "Display size" in said and "Auto" in said
