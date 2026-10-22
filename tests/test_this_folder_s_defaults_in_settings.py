"""Settings' This folder changes the workspace's defaults file, and the
Agent's mode and ceiling are saved (IMDE 2.9, the parts in `server.py` and
`agent-panel.js`).

The page is `test_one_settings_page.py`'s; here are the server's route to
the file (`/api/defaults-file`, read and written by `settings_page.py`) and
the Agent's choices saved over a reload, which until now went back to Seeker
with no ceiling.
"""

from __future__ import annotations

import json

import pytest
import yaml

from tests.test_the_workspace_says_its_studies import _study


def _ask(study, hosting, method="GET", body=None, secret="s" * 40):
    import threading
    import urllib.error
    import urllib.request
    from http.server import ThreadingHTTPServer

    from fastmdxplora.gui.exploration import DashboardRuntime
    from fastmdxplora.gui.hosting import SECRET_HEADER
    from fastmdxplora.gui.server import make_handler

    runtime = DashboardRuntime(workspace_root=study.parent, exploration_root=study.parent)
    runtime.hosting = hosting
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(
        study.parent, runtime=runtime, allow_control=True, hosting=hosting))
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    headers = {SECRET_HEADER: secret, "Host": "app.example.org"} if hosting else {}
    data = None
    if body is not None:
        data = json.dumps(body).encode("utf-8")
        headers["content-type"] = "application/json"
        if hosting:
            headers["Origin"] = "https://app.example.org"
    try:
        request = urllib.request.Request(
            f"http://127.0.0.1:{httpd.server_address[1]}/api/defaults-file",
            data=data, method=method, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                return response.status, json.loads(response.read())
        except urllib.error.HTTPError as refused:
            return refused.code, json.loads(refused.read() or b"{}")
    finally:
        httpd.shutdown()
        httpd.server_close()


def test_the_route_reads_and_writes_the_workspace_s_file(tmp_path):
    study = _study(tmp_path / "person" / "study")
    where = study.parent / "fastmdx-defaults.yml"
    status, said = _ask(study, None)
    assert status == 200 and said["path"] == str(where) and not said["exists"]
    status, said = _ask(study, None, "POST", {"set": {"simulation.temperature_K": 310},
                                              "why": {"simulation.temperature_K": "warm"}})
    assert status == 200 and said["ok"] and said["values"] == {"simulation.temperature_K": 310}
    assert yaml.safe_load(where.read_text(encoding="utf-8"))["decisions"] == {
        "simulation.temperature_K": {"why": "warm"}}
    status, said = _ask(study, None, "POST", {"set": {"simulation.temperature_K": "warm"}})
    assert status == 200 and said["ok"] is False and "temperature_K" in said["reason"]
    assert said["values"] == {"simulation.temperature_K": 310}


def test_a_service_s_guests_neither_read_nor_write_it(tmp_path):
    from fastmdxplora.gui.hosting import Hosting

    study = _study(tmp_path / "person" / "study")
    where = study.parent / "fastmdx-defaults.yml"
    where.write_text("setup:\n  ph: 7.0\n", encoding="utf-8")
    hosting = Hosting(workspace=study.parent.resolve(),
                      allowed_hosts=frozenset({"app.example.org"}), secret="s" * 40)
    read = _ask(study, hosting)
    written = _ask(study, hosting, "POST", {"set": {"setup.ph": 6.0}})
    assert read[0] == 403 and read[1]["ok"] is False and "values" not in read[1]
    assert written[0] == 403 and written[1]["ok"] is False
    assert where.read_text(encoding="utf-8") == "setup:\n  ph: 7.0\n"


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
    yield started, study.parent / "fastmdx-defaults.yml"
    started.server.shutdown()


def _open(browser, url, where="#settings"):
    page = browser.new_context(viewport={"width": 1440, "height": 900}).new_page()
    page.set_default_timeout(60000)
    page.errors = []
    page.on("pageerror", lambda error: page.errors.append(str(error)))
    page.goto(url + where, wait_until="domcontentloaded")
    page.wait_for_function("() => window.FastMDXDashboard && window.FastMDXSettings")
    page.wait_for_selector("body:not(.state-loading)")
    return page


ROW = '.settings-field[data-name="simulation.temperature_K"]'


def test_a_setting_of_this_folder_is_set_given_its_reason_undone_and_taken_out(browser, session):
    started, where = session
    page = _open(browser, started.url)
    page.click("#settings-scope-folder")
    page.wait_for_selector(f"{ROW} .settings-control")
    lead = page.text_content("#settings-folder-lead")
    theirs = page.text_content(f"{ROW} .settings-theirs")
    # Set, as it is changed: no Save.
    page.fill(f"{ROW} .settings-control", "310")
    page.press(f"{ROW} .settings-control", "Enter")
    page.wait_for_selector(f"{ROW}.is-set")
    set_said = page.text_content(f"{ROW} .settings-status")
    page.fill(f"{ROW} .settings-why", "the body's temperature")
    page.press(f"{ROW} .settings-why", "Tab")
    page.wait_for_function(f"() => /Reason saved/.test(document.querySelector('{ROW} .settings-status')"
                           ".textContent)")
    written = yaml.safe_load(where.read_text(encoding="utf-8"))
    # The New study form offers it as its default.
    offered = page.evaluate("""async () => {
        const s = await (await fetch('/api/schema')).json();
        return s.phases.simulation.fields.find((f) => f.name === 'temperature_K').default; }""")
    # Refused: said beside it, the file as it was.
    page.fill(f"{ROW} .settings-control", "-4")
    page.press(f"{ROW} .settings-control", "Enter")
    page.wait_for_selector(f"{ROW} .settings-status.is-refused")
    refused = (page.text_content(f"{ROW} .settings-status"),
               page.get_attribute(f"{ROW} .settings-control", "aria-invalid"))
    after_refusal = yaml.safe_load(where.read_text(encoding="utf-8"))
    # Found by its words, This folder's count beside its name.
    page.fill("#settings-find", "production temperature")
    count = page.text_content("#settings-scope-folder .settings-count")
    page.fill("#settings-find", "")
    # FastMDXplora's again, and Undo puts it back with its reason.
    page.click(f"{ROW} .settings-use-theirs")
    page.wait_for_selector(f"{ROW}:not(.is-set)")
    taken_out = yaml.safe_load(where.read_text(encoding="utf-8")) if where.read_text(
        encoding="utf-8").strip() else {}
    page.click(f"{ROW} .settings-undo")
    page.wait_for_selector(f"{ROW}.is-set")
    put_back = yaml.safe_load(where.read_text(encoding="utf-8"))
    why_back = page.input_value(f"{ROW} .settings-why")
    errors = page.errors
    page.context.close()
    assert errors == []
    assert str(where) in lead and "Not there yet" in lead
    assert theirs.startswith("FastMDXplora's: 300")
    assert set_said.startswith("Saved in the file.")
    assert written == {"simulation": {"temperature_K": 310},
                       "decisions": {"simulation.temperature_K": {"why": "the body's temperature"}}}
    assert offered == 310
    assert refused[1] == "true" and "temperature_K" in refused[0] and str(where) not in refused[0]
    assert after_refusal == written
    assert count == "1"
    assert taken_out == {}
    assert put_back == written and why_back == "the body's temperature"


def test_a_file_with_comments_is_shown_and_not_written(browser, session):
    started, where = session
    where.write_text("simulation:\n  temperature_K: 305  # the lab's\n", encoding="utf-8")
    page = _open(browser, started.url)
    page.click("#settings-scope-folder")
    page.wait_for_selector(f"{ROW} .settings-fixed")
    said = (page.text_content("#settings-folder-body .settings-cannot"),
            page.text_content(f"{ROW} .settings-fixed"),
            page.locator(f"{ROW} input, {ROW} select").count())
    page.context.close()
    assert "comments" in said[0] and said[1] == "305" and said[2] == 0


def test_the_agent_s_mode_and_ceiling_are_saved(browser, session):
    started, _ = session
    page = _open(browser, started.url, "#agent")
    page.evaluate("() => window.FastMDXAgent.openSettings()")
    page.select_option("#agent-mode", "autonomous")
    page.fill("#agent-budget", "12")
    page.dispatch_event("#agent-budget", "change")
    page.keyboard.press("Escape")
    page.reload(wait_until="domcontentloaded")
    page.wait_for_selector("body:not(.state-loading)")
    page.wait_for_function("() => document.getElementById('sb-agent')"
                           " && !document.getElementById('sb-agent').hidden")
    said = page.evaluate("""() => [document.getElementById('agent-mode').value,
        document.getElementById('agent-budget').value,
        document.getElementById('settings-mode').textContent,
        document.getElementById('settings-page-mode').textContent,
        document.getElementById('sb-agent').textContent]""")
    errors = page.errors
    page.context.close()
    assert errors == []
    assert said[:2] == ["autonomous", "12"]
    assert said[2] == said[3] == "Ranger: built-in, autonomous"
    assert "Ranger" in said[4]
