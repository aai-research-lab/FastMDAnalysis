"""A file asked for on the Files page stays open in the side panel.

The panel closes by itself on a study's first status when the study is
finished, or silent (no live record, nothing in its log): the log is then a
record the Files and Report pages hold better. View on the Files page opens
the panel and the file in it; asked before that first status came, the
panel was then closed under the file, as if nobody had opened it (CI on
99b21cc6, a browser test waiting 20 s for its View button).
"""

from __future__ import annotations

import pytest

from tests.test_hostile_text_renders_inert import _hostile_study

#: The study's own status is never answered, so the first status is the one
#: a test sends, at the moment it chooses.
HOLD_THE_STATUS = """
(() => {
  const real = window.fetch;
  window.fetch = function (url, ...rest) {
    if (String(url).includes("/api/status")) return new Promise(() => {});
    return real.call(this, url, ...rest);
  };
})();
"""
FIRST_STATUS = {
    "finished": {"status": {"status": "completed"}, "health": {}, "study": ""},
    "silent": {"status": {}, "health": {}, "study": ""},
}
VIEW = "#side-preview-modes:not([hidden]) [data-mode=view]"


@pytest.fixture()
def page():
    sync_playwright = pytest.importorskip("playwright.sync_api").sync_playwright
    from fastmdxplora.gui.server import start_dashboard_session

    session = start_dashboard_session(output=str(_hostile_study()), host="127.0.0.1",
                                      port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page()
            page.set_default_timeout(20000)
            page.add_init_script(HOLD_THE_STATUS)
            page.goto(session.url, wait_until="domcontentloaded")
            page.wait_for_function("!!(window.FastMDXFrame && window.FastMDXFrame.previewPath)")
            # The log read, empty: what a silent study needs besides its status.
            page.wait_for_selector(".side-log-empty", state="attached")
            yield page
            browser.close()
    finally:
        session.server.shutdown()
        session.server.server_close()


def _first_status(page, kind: str) -> None:
    page.evaluate("d => window.dispatchEvent(new CustomEvent('dashboard:status-updated', "
                  "{detail: d}))", FIRST_STATUS[kind])


def _panel_open(page) -> bool:
    return not page.evaluate("document.body.classList.contains('panel-collapsed')")


@pytest.mark.parametrize("kind", sorted(FIRST_STATUS))
def test_a_file_opened_before_the_first_status_stays_open(page, kind) -> None:
    page.evaluate("window.FastMDXFrame.previewPath('analysis/notes.md')")
    page.wait_for_selector(VIEW)
    _first_status(page, kind)
    assert page.is_visible(VIEW)
    assert _panel_open(page)


@pytest.mark.parametrize("kind", sorted(FIRST_STATUS))
def test_with_nothing_asked_the_first_status_still_closes_it(page, kind) -> None:
    page.evaluate("window.FastMDXFrame.setCollapsed(false, false)")
    assert _panel_open(page)
    _first_status(page, kind)
    assert not _panel_open(page)
