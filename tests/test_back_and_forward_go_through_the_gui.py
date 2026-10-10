"""Back and Forward go through the GUI's pages and dialogs.

Pages were changed with the address replaced, not stepped: Overview, Viewer,
Analysis, then Back left the GUI for whatever came before it. Each page is a
step of its own now, and each dialog opened is one over its page: Back
closes it before it leaves the page, Forward opens it again, and a dialog
closed by its own Close takes its step back with it, so Back then goes to
the page before. Where the GUI opens is put right in its address, not a step.
An analysis opened from elsewhere is a step to its place on the Analysis
page (`#analysis/rmsd`), which a reload opens again.
"""

from __future__ import annotations

import pytest

pytest.importorskip("playwright.sync_api")
pytest.importorskip("mdtraj")

from tests.test_the_drawing_scripts_run_in_a_browser import _write_study  # noqa: E402


@pytest.fixture(scope="module")
def session(tmp_path_factory):
    from fastmdxplora.gui.server import start_dashboard_session

    from tests.test_an_analysis_is_drawn_from_its_numbers import _figure
    from tests.test_the_analysis_page_reads_as_a_whole import _analysis, _ar1, _manifest

    study = _write_study(tmp_path_factory.mktemp("back") / "study")
    names = ["rmsd", "rg", "sasa", "rmsf", "hbonds", "ligand_rmsd"]
    _manifest(study, {name: {"status": "ok"} for name in names})
    for seed, name in enumerate(names):
        _analysis(study, name, 1.0 + 0.01 * _ar1(4000, 0.5, seed + 1))
        _figure(study / "analysis" / name / f"{name}.png")
    started = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    yield started
    started.server.shutdown()


@pytest.fixture(scope="module")
def browser():
    from playwright.sync_api import sync_playwright

    with sync_playwright() as pw:
        launched = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
        yield launched
        launched.close()


def _open(browser, session):
    page = browser.new_context(viewport={"width": 1440, "height": 900}).new_page()
    page.set_default_timeout(60000)
    page.errors = []
    page.on("pageerror", lambda error: page.errors.append(str(error)))
    page.goto("about:blank")
    page.goto(session.url + "#overview", wait_until="domcontentloaded")
    page.wait_for_function("() => document.body.classList.contains('state-ready')")
    return page


def _on(page, name):
    page.wait_for_function("(name) => location.hash === '#' + name && "
                           "!document.querySelector(`.page[data-page='${name}']`).hidden", arg=name)


def _shown(page, dialog):
    return page.evaluate("(id) => !document.getElementById(id).hidden", dialog)


def _settled(page):
    # The dialog's step taken back: the step now is the page's own.
    page.wait_for_function("() => !(history.state && history.state.fmxDialog)")


def test_back_and_forward_go_between_pages(browser, session) -> None:
    page = _open(browser, session)
    page.click("[data-view-link='viewer']")
    _on(page, "viewer")
    page.click("[data-view-link='analysis']")
    _on(page, "analysis")
    page.go_back()
    _on(page, "viewer")
    page.go_back()
    _on(page, "overview")
    page.go_forward()
    _on(page, "viewer")
    page.go_back()
    _on(page, "overview")
    # Where it opened is not a step of its own: one more Back leaves it.
    page.go_back()
    left = page.url
    errors = page.errors
    page.context.close()
    assert left == "about:blank"
    assert errors == []


def test_back_closes_a_dialog_before_it_leaves_the_page(browser, session) -> None:
    page = _open(browser, session)
    page.click("[data-view-link='viewer']")
    _on(page, "viewer")
    page.click("#settings-open")
    page.click("[data-dialog-open='prefs-dialog']")
    assert _shown(page, "prefs-dialog")
    page.go_back()
    page.wait_for_function("() => document.getElementById('prefs-dialog').hidden")
    _on(page, "viewer")
    page.go_forward()
    page.wait_for_function("() => !document.getElementById('prefs-dialog').hidden")
    _on(page, "viewer")
    # Closed by its own Close: its step goes with it, so Back is the page
    # before, not the dialog's step again.
    page.click("#prefs-dialog [data-dialog-close]")
    _settled(page)
    page.go_back()
    _on(page, "overview")
    shown = _shown(page, "prefs-dialog")
    errors = page.errors
    page.context.close()
    assert not shown
    assert errors == []


def test_a_link_to_a_dialog_opens_it_over_the_page(browser, session) -> None:
    page = _open(browser, session)
    page.click("[data-view-link='files']")
    _on(page, "files")
    page.evaluate("() => { location.hash = '#cite'; }")
    page.wait_for_function("() => !document.getElementById('cite-dialog').hidden")
    _on(page, "files")
    page.go_back()
    page.wait_for_function("() => document.getElementById('cite-dialog').hidden")
    _on(page, "files")
    page.go_back()
    _on(page, "overview")
    page.context.close()


def test_a_page_opened_as_a_dialog_closes_is_not_taken_back_with_it(browser, session) -> None:
    """Go to anything closes and opens a page at once: the page's step was
    taken back with the dialog's, and Back went nowhere."""
    page = _open(browser, session)
    page.evaluate("""() => { window.FastMDXDialog.open('cite-dialog');
                             window.FastMDXDialog.close('cite-dialog');
                             window.FastMDXDashboard.navigate('files'); }""")
    _settled(page)
    _on(page, "files")
    page.go_back()
    _on(page, "overview")
    page.go_forward()
    _on(page, "files")
    errors = page.errors
    page.context.close()
    assert errors == []


def test_an_analysis_opened_is_a_place_a_reload_opens_again(browser, session) -> None:
    page = _open(browser, session)
    page.evaluate("() => window.FastMDXDashboard.showAnalysis('ligand_rmsd')")
    page.wait_for_function("() => location.hash === '#analysis/ligand_rmsd'")
    _on_place = ("() => { const c = document.querySelector("
                 "'.page[data-page=\"analysis\"] .analysis-card[data-analysis=\"ligand_rmsd\"]');"
                 " if (!c || c.offsetParent === null) return false;"
                 " const r = c.getBoundingClientRect(); return r.bottom > 0 && r.top < innerHeight; }")
    page.wait_for_function(_on_place)
    page.go_back()
    _on(page, "overview")
    page.go_forward()
    page.wait_for_function("() => location.hash === '#analysis/ligand_rmsd' && "
                           "!document.querySelector('.page[data-page=\"analysis\"]').hidden")
    page.reload(wait_until="domcontentloaded")
    page.wait_for_function(_on_place)
    kept = page.evaluate("() => location.hash")
    errors = page.errors
    page.context.close()
    assert kept == "#analysis/ligand_rmsd"
    assert errors == []
