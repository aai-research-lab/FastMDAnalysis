"""One confirmation, one empty, loading and error state, one list of
conversations (IMDE 2.10).

Stopping a run, stopping a job on another machine and clearing the scratch
each asked in a way of their own (a line, another line, a dialog); they ask
with one component now, beside what was pressed, Not now focused, Escape for
no. What can be undone (a selection forgotten) asks nothing and offers Undo.
A page with nothing to show, still reading or not able to read says it with
one block: what, why, and the one thing to do; the Files page said "Reading
the study's files..." for ever when the server did not answer. And the
sidebar's Conversations and Chats are one list, the study's first.
"""

from __future__ import annotations

import pytest

from tests.test_the_workspace_says_its_studies import _study

ASK = """async () => {
    const host = document.createElement('div');
    host.id = 'an-ask';
    document.body.appendChild(host);
    const answer = window.FastMDXConfirm.ask(host, {title: 'Stop the run?',
        detail: 'It ends at its next frame.', yes: 'Yes, stop it', list: ['seed 1', 'seed 2']});
    const shown = {
        text: host.querySelector('.fmx-ask').innerText.replace(/\\s+/g, ' ').trim(),
        focus: document.activeElement.textContent,
        role: host.querySelector('.fmx-ask').getAttribute('role'),
        list: [...host.querySelectorAll('.fmx-ask-list li')].map((li) => li.textContent),
    };
    document.activeElement.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
    shown.escaped = await answer;
    const again = window.FastMDXConfirm.ask(host, {title: 'Clear it?', yes: 'Clear',
        detail: Promise.resolve({detail: '3 files, 1.2 MB.', ready: true})});
    shown.waiting = host.querySelector('.fmx-ask-yes').disabled;
    await new Promise((r) => setTimeout(r, 50));
    shown.ready = !host.querySelector('.fmx-ask-yes').disabled;
    shown.read = host.querySelector('.fmx-ask-detail').textContent;
    host.querySelector('.fmx-ask-yes').click();
    shown.yes = await again;
    window.FastMDXConfirm.done(host);
    shown.gone = host.hidden && !host.childNodes.length;
    return shown;
}"""


def test_one_ask_one_state_and_one_list(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    study = _study(tmp_path / "alpha")
    _study(tmp_path / "beta")
    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    post = """async ([url, body]) => (await fetch(url, {method: 'POST',
        headers: {'content-type': 'application/json'}, body: JSON.stringify(body)})).json()"""
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            # The Files page's read refused: said, with Try again.
            refused = {"now": True}

            def files(route):
                if refused["now"]:
                    route.fulfill(status=500, body="no")
                else:
                    route.continue_()
            page.route("**/api/files-page*", files)
            page.goto(session.url + "#files", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.wait_for_selector('#files-page .fmx-state[data-kind="error"]')
            files_said = page.inner_text("#files-page .fmx-state")
            refused["now"] = False
            page.click("#files-page .fmx-state-act")
            page.wait_for_function("() => !document.querySelector('#files-page .fmx-state')")
            asked = page.evaluate(ASK)
            # Nothing to compare yet.
            page.goto(session.url + "#compare", wait_until="domcontentloaded")
            page.wait_for_selector('#compare-empty .fmx-state[data-kind="empty"]')
            compare = page.inner_text("#compare-empty")
            # All studies: none matching, then not listed with Try again.
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            page.fill("#studies-search", "nothing like it")
            none = page.inner_text("#studies-empty")
            page.fill("#studies-search", "")
            page.route("**/api/studies*", lambda route: route.fulfill(
                json={"ok": False, "reason": "No such folder: /gone"}))
            page.evaluate("() => window.FastMDXStudies.load('/gone')")
            page.wait_for_selector('#studies-empty .fmx-state[data-kind="error"]')
            listed = page.inner_text("#studies-empty")
            page.unroute("**/api/studies*")
            page.click("#studies-empty .fmx-state-act")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            # One list of conversations: the study's, then the chats.
            run = page.evaluate("() => window.FastMDXDashboard.state.appState.active_run")
            for where, text in ((None, "A chat about nothing open"), (run, "About this study")):
                begun = page.evaluate(post, ["/api/agent/conversation/new", {"study": where}])
                page.evaluate(post, ["/api/agent/conversation", {
                    "id": begun["id"], "study": where, "entries": [{"role": "user", "text": text}]}])
            page.evaluate("() => window.FastMDXChats.show()")
            page.wait_for_selector("#sidebar-chats-list .sidebar-conv")
            talks = page.evaluate("""() => ({
                head: document.querySelector('#sidebar-chats > summary').innerText.trim(),
                add: document.getElementById('sidebar-chat-new').getAttribute('aria-label'),
                rows: [...document.querySelectorAll('#sidebar-chats-list > *')]
                    .filter((n) => !n.classList.contains('sidebar-conv-new'))
                    .map((n) => n.classList.contains('sidebar-conv-part') ? '# ' + n.textContent
                        : n.dataset.group + ': ' + n.querySelector('.sidebar-conv-title').textContent),
                two: !!document.getElementById('sidebar-study-convs'),
            })""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert "The study's files could not be read" in files_said and "Try again" in files_said
    assert asked == {
        "text": "Stop the run? It ends at its next frame. seed 1 seed 2 Yes, stop it Not now",
        "focus": "Not now", "role": "group", "list": ["seed 1", "seed 2"], "escaped": False,
        "waiting": True, "ready": True, "read": "3 files, 1.2 MB.", "yes": True, "gone": True}
    assert compare.startswith("Nothing to compare yet")
    assert none.startswith("No study matches that")
    assert listed.startswith("The studies could not be listed") and "No such folder: /gone" in listed
    assert talks == {
        "head": "Conversations", "add": "A new conversation about this study",
        "rows": ["# This study", "study: About this study",
                 "# Not about a study", "chats: A chat about nothing open"],
        "two": False}
