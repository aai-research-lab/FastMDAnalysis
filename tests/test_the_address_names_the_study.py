"""The address names the open study, so Back returns to the study before.

Opening a study switched the server's study and reloaded the page at
`#overview`: the address never said which study it showed, so Back could
not return to the study open before and a bookmark could not hold one. The
address now names the study by the server's key for it (a hash of its
folder, never its path): a study opened from the page is a step of its own,
Back and Forward open the study each step names, and an address naming
another study opens it, or says it is not among the studies here.
"""

from __future__ import annotations

import pytest

from fastmdxplora.gui.workspace import study_key, study_keyed
from tests.test_the_workspace_says_its_studies import _replicas, _study


def test_a_study_is_found_by_its_key_among_those_near(tmp_path) -> None:
    workspace = tmp_path / "workspace"
    first = _study(workspace / "lab" / "ubiquitin", system="1UBQ")
    second = _study(workspace / "chignolin", system="5AWL")
    sweep = _replicas(workspace / "seeds")
    run = next(p for p in (sweep / "runs").iterdir() if p.is_dir())
    (run / "manifest.json").write_text("{}", encoding="utf-8")
    near, under = [first, None], [workspace]
    assert study_keyed(study_key(second), near, under) == second.resolve()
    # A run inside a study of several runs, as the Viewer opens one.
    assert study_keyed(study_key(run), near, under) == run.resolve()
    # The folders around the one open are not looked through.
    assert study_keyed(study_key(second), near) is None
    # The key names the folder and nothing of its path.
    assert len(study_key(first)) == 16 and "ubiquitin" not in study_key(first)
    assert study_keyed("0" * 16, near, under) is None
    assert study_keyed("../" + study_key(second), near, under) is None


@pytest.fixture(scope="module")
def session(tmp_path_factory):
    pytest.importorskip("playwright.sync_api")
    from fastmdxplora.gui.server import start_dashboard_session

    workspace = tmp_path_factory.mktemp("address")
    first = _study(workspace / "ubiquitin", system="1UBQ", started="2026-09-01T10:00:00+00:00")
    second = _study(workspace / "chignolin", system="5AWL", started="2026-09-02T10:00:00+00:00")
    started = start_dashboard_session(output=str(first), host="127.0.0.1", port=0)
    started.first, started.second = first.resolve(), second.resolve()
    yield started
    started.server.shutdown()


@pytest.fixture(scope="module")
def browser():
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


def _open(browser, session, address="#overview"):
    page = browser.new_context(viewport={"width": 1440, "height": 900}).new_page()
    page.set_default_timeout(60000)
    page.errors = []
    page.on("pageerror", lambda error: page.errors.append(str(error)))
    page.goto("about:blank")
    page.goto(session.url.rstrip("/") + "/" + address, wait_until="domcontentloaded")
    return page


def _named(page, folder, where="overview"):
    """The address names ``folder``'s study, on ``where``, the server shows
    that study, and the page has read it (a page shown again by the browser
    may load itself again first, so each look may meet a page going)."""
    from playwright.sync_api import Error

    key = study_key(folder)
    said = None
    for _ in range(300):
        try:
            said = page.evaluate("""([key, where, folder]) => fetch('/api/app-state')
                .then(r => r.json()).then(s => [new URL(location.href).searchParams.get('study'),
                    location.hash, s.active_run, document.readyState,
                    document.body.classList.contains('state-ready')])""",
                                 [key, where, str(folder)])
            if said == [key, "#" + where, str(folder), "complete", True]:
                return
        except Error:
            said = "a page going"
        page.wait_for_timeout(100)
    raise AssertionError(said)


def _switch(page, folder):
    page.evaluate("""(folder) => fetch('/api/explore/switch', {method: 'POST',
        headers: {'content-type': 'application/json'}, body: JSON.stringify({folder})})
        .then(r => r.json()).then(answer => window.FastMDXDashboard.studyOpened(answer))""",
                  str(folder))


def test_back_and_forward_open_the_study_each_step_names(browser, session) -> None:
    page = _open(browser, session)
    _named(page, session.first)
    address = page.url
    _switch(page, session.second)
    _named(page, session.second)
    page.go_back()
    _named(page, session.first)
    page.go_forward()
    _named(page, session.second)
    page.go_back()
    _named(page, session.first)
    errors = page.errors
    page.context.close()
    assert str(session.first) not in address and "ubiquitin" not in address
    assert errors == []


def test_an_address_naming_another_study_opens_it_on_its_page(browser, session) -> None:
    page = _open(browser, session, f"?study={study_key(session.second)}#viewer")
    _named(page, session.second, "viewer")
    page.context.close()
    # And back to the first, for the tests after.
    page = _open(browser, session, f"?study={study_key(session.first)}#overview")
    _named(page, session.first)
    page.context.close()


def test_an_address_naming_no_study_here_says_so(browser, session) -> None:
    page = _open(browser, session, "?study=0123456789abcdef#overview")
    page.wait_for_function("() => /not among the studies here/.test("
                           "(document.getElementById('dashboard-toast') || {}).textContent || '')")
    _named(page, session.first)
    page.context.close()


def test_back_to_a_step_not_yet_naming_the_study_opens_its_page(browser, session) -> None:
    """The address is put right to name the study once the page has read
    it; a page opened before that left a step without the study's name,
    and Back to it changed the address but left the page shown."""
    page = _open(browser, session)
    _named(page, session.first)
    # The step as it was before the study was named in it.
    page.evaluate("() => history.replaceState(history.state, '', location.pathname + '#overview')")
    page.click("[data-view-link='viewer']")
    page.wait_for_function("() => location.hash === '#viewer' && "
                           "new URL(location.href).searchParams.get('study')")
    page.go_back()
    page.wait_for_function("() => location.hash === '#overview' && "
                           "!document.querySelector(`.page[data-page='overview']`).hidden")
    errors = page.errors
    page.context.close()
    assert errors == []


def test_a_reload_keeps_the_study_opened_elsewhere(browser, session) -> None:
    """Another study's conversation opened from the Agent's list switches
    the server and reloads the page: the address, naming the study before,
    opened that one again."""
    page = _open(browser, session)
    _named(page, session.first)
    page.evaluate("""(folder) => fetch('/api/explore/switch', {method: 'POST',
        headers: {'content-type': 'application/json'}, body: JSON.stringify({folder})})""",
                  str(session.second))
    page.reload(wait_until="domcontentloaded")
    _named(page, session.second)
    errors = page.errors
    page.context.close()
    assert errors == []


def test_an_address_naming_the_open_study_asks_nothing(browser, session) -> None:
    """A dashboard watched from another machine may not switch: a page
    opened on the study already open asked to, was refused, and said so."""
    import json
    import urllib.request

    with urllib.request.urlopen(session.url.rstrip("/") + "/api/status", timeout=30) as answer:
        key = json.loads(answer.read())["study"]
    page = browser.new_context(viewport={"width": 1440, "height": 900}).new_page()
    page.set_default_timeout(60000)
    asked = []
    page.on("request", lambda request: asked.append(request.url)
            if "/api/explore/switch" in request.url else None)
    page.goto(session.url.rstrip("/") + f"/?study={key}#overview", wait_until="domcontentloaded")
    page.wait_for_function("() => document.body.classList.contains('state-ready')")
    page.wait_for_timeout(1500)
    page.context.close()
    assert asked == []
