"""The study environment as a person uses it, after its first review.

Option and a digit typed in a field went to another page and lost the
character; the run menu closed under the pointer, and the runs table lost
the keyboard, each time a run's progress moved; a run of a study of runs was
named by its folder in the Runs compared table, its tab and the sidebar; the
Overview of a study of runs said "What the run did"; the runs table hid that
it scrolls; the Agent beside the page went on under the status bar; the
Viewer said a 40 ps run's frames in ns under a ruler in ps; a track's note
was drawn over its name.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

pytest.importorskip("playwright.sync_api")

from tests.test_the_ide_bugs_in_a_browser import DONE, _serve  # noqa: E402


@pytest.fixture(scope="module")
def sweep(tmp_path_factory) -> Path:
    from tests.test_a_study_of_runs_is_shown_as_one import _a_sweep_study

    root = _a_sweep_study(tmp_path_factory.mktemp("used"))
    for folder in sorted((root / "runs").iterdir()):
        (folder / "simulation" / "live_status.json").write_text(json.dumps(DONE), encoding="utf-8")
    return root


def test_option_and_a_digit_typed_in_a_field_stay_in_it(sweep):
    def look(page, url):
        page.add_init_script("Object.defineProperty(navigator, 'platform', {get: () => 'MacIntel'})")
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector(".page[data-page=overview] .place-run")
        page.evaluate("""() => {
            const box = document.createElement('input');
            box.id = 'typed-here';
            document.querySelector('.page[data-page=overview]').appendChild(box);
            box.focus();
        }""")
        page.keyboard.press("Alt+Digit3")
        typed = page.evaluate("() => [document.documentElement.dataset.page,"
                              " document.activeElement.id]")
        page.locator("body").click(position={"x": 5, "y": 5})
        page.keyboard.press("Alt+Digit3")
        page.wait_for_function("() => document.documentElement.dataset.page !== 'overview'")
        return typed

    assert _serve(sweep, look) == ["overview", "typed-here"]


def test_the_run_menu_stays_open_while_runs_move_and_the_table_keeps_the_keyboard(sweep):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector(".page[data-page=overview] .place-run")
        page.wait_for_selector("#overview-runs:not([hidden]) .overview-runs-open")
        page.click(".page[data-page=overview] .place-run")
        page.wait_for_selector(".place-menu-item:focus")
        page.evaluate("""() => {
            const d = window.FastMDXDashboard;
            const next = d.state.appState;
            next.runs[0].state = 'running';
            next.runs[0].fraction = 0.42;
            d.applyAppState(next);
        }""")
        menu = page.evaluate("""() => [window.FastMDXPlace.menuOpen(),
            document.activeElement.className,
            document.querySelectorAll('.place-menu-item .place-menu-said')[1].textContent]""")
        page.keyboard.press("Escape")
        # The runs table, drawn again, keeps the run the keyboard is on.
        page.locator(".overview-runs-open").nth(1).focus()
        page.evaluate("""() => {
            const s = window.FastMDXOverviewRuns;
            return fetch('/api/runs-compared').then((r) => r.json()).then((d) => {
                d.runs[1].production_ns = 0.123;
                s.show(d);
            });
        }""")
        table = page.evaluate("""() => [document.activeElement.className,
            document.activeElement.closest('.overview-runs-row').dataset.run,
            document.querySelector('.overview-runs-row:nth-child(2) td:nth-child(3)').textContent]""")
        return menu, table

    menu, table = _serve(sweep, look)
    assert menu == [True, "place-menu-item", "42%"]
    assert table[0] == "overview-runs-open" and table[2].startswith("123 ps")


def test_a_run_of_a_study_is_named_by_what_sets_it_apart_and_the_study_by_its_runs(sweep):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector("#overview-runs:not([hidden]) .overview-runs-open")
        head = page.text_content("#overview-subtitle")
        page.goto(url + "#analysis", wait_until="domcontentloaded")
        page.wait_for_selector("#runs-compared-table tbody th")
        table = page.evaluate("""() => ({
            runs: [...document.querySelectorAll('#runs-compared-table tbody th')]
                .map((th) => [th.textContent, th.title]),
            axes: [...document.querySelectorAll('#runs-compared-table th.is-varied')]
                .map((th) => th.textContent)})""")
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector("#overview-runs:not([hidden]) .overview-runs-open")
        page.locator(".overview-runs-open").first.click()
        # The page loads again as the run opens: asked of the new one.
        page.wait_for_function("() => window.FastMDXDashboard"
                               " && (window.FastMDXDashboard.state.appState || {}).run_of")
        page.wait_for_function("() => document.getElementById('sidebar-run-name')"
                               ".textContent === 'temperature 300 K'")
        # Its tab, named beside another of its system by its study and run.
        page.wait_for_function("""() => Object.keys(sessionStorage)
            .filter((k) => k.startsWith('fmx.studyTabs:'))
            .some((k) => (JSON.parse(sessionStorage.getItem(k)).tabs || [])
                .some((t) => t.folder.endsWith('\u203a temperature 300 K')))""")
        return head, table

    head, table = _serve(sweep, look)
    assert head == "What the study's runs did, and how they went."
    assert [label for label, _ in table["runs"]] == ["temperature 300 K", "temperature 310 K"]
    assert all(folder and folder != label for label, folder in table["runs"])
    assert table["axes"] == ["temperature (K)"]


def test_the_runs_table_shows_it_scrolls_and_the_agent_ends_above_the_status_bar(sweep):
    def look(page, url):
        page.goto(url + "#overview", wait_until="domcontentloaded")
        page.wait_for_selector("#overview-runs:not([hidden]) .overview-runs-open")
        more = page.evaluate("""() => {
            const scroll = document.querySelector('#overview-runs .runs-compared-scroll');
            scroll.style.width = '200px';
            return new Promise((done) => requestAnimationFrame(() => requestAnimationFrame(
                () => done(scroll.dataset.more || ''))));
        }""")
        page.keyboard.press("Control+KeyJ")
        page.wait_for_selector(".agent-drawer:not([hidden])")
        edges = page.evaluate("""() => [
            document.querySelector('.agent-drawer').getBoundingClientRect().bottom,
            document.querySelector('.status-bar').getBoundingClientRect().top]""")
        return more, edges

    for width in (1250, 1440):
        more, edges = _serve(sweep, look, width=width)
        assert more == "right"
        assert edges[0] <= edges[1] + 0.5
