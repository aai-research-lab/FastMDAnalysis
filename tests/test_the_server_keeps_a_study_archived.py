"""The server keeps a study archived where All studies' bar asks.

Archive on All studies' bar posts `/api/study-tags` with ``archived``
(studies.js `setArchived`): the server passes it on to `set_tags`, which
keeps it in the study's own `study_tags.json`, nothing of the study moved.
A server that dropped it kept the tags alone, and the page said this GUI's
server does not keep Archive yet.
"""

from __future__ import annotations

import json
import urllib.request

import pytest

from tests.test_the_workspace_says_its_studies import _study


def _post(url: str, body: dict) -> dict:
    request = urllib.request.Request(url + "/api/study-tags", data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.loads(response.read())


def test_the_server_passes_archived_on_and_back(tmp_path):
    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.gui.workspace import studies_in
    from fastmdxplora.study_tags import tags_of

    _study(tmp_path / "wild_type", started="2026-09-02T10:00:00+00:00")
    study = _study(tmp_path / "mutant", started="2026-09-03T10:00:00+00:00")
    records = {p.name: p.read_bytes() for p in study.iterdir() if p.is_file()}
    session = start_dashboard_session(output=str(tmp_path / "wild_type"), host="127.0.0.1",
                                      port=0)
    try:
        archived = _post(session.url, {"path": str(study), "tags": ["L67S"], "archived": True})
        on_disk = json.loads((study / "study_tags.json").read_text(encoding="utf-8"))
        listed = {card["name"]: card for card in studies_in(tmp_path)["studies"]}
        # Tags set by a card, with no word on it: still archived.
        kept = _post(session.url, {"path": str(study), "tags": ["L67S", "redo"]})
        back = _post(session.url, {"path": str(study), "tags": ["L67S"], "archived": False})
    finally:
        session.server.shutdown()
    assert archived == {"ok": True, "tags": ["L67S"], "note": "", "archived": True}
    assert on_disk["archived"].endswith("+00:00")
    assert listed["mutant"]["archived"] is True and "archived" not in listed["wild_type"]
    assert kept["archived"] is True
    assert back == {"ok": True, "tags": ["L67S"], "note": ""}
    assert "archived" not in tags_of(study)
    # Nothing of the study moved or changed but its tags.
    assert {p.name: p.read_bytes() for p in study.iterdir()
            if p.is_file() and p.name != "study_tags.json"} == records


def test_archive_on_the_bar_is_kept_and_undone(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.study_tags import tags_of

    for k, name in enumerate(("alpha", "beta", "gamma")):
        _study(tmp_path / name, started=f"2026-09-0{k + 1}T10:00:00+00:00")
    session = start_dashboard_session(output=str(tmp_path / "alpha"), host="127.0.0.1", port=0)
    rows = ("() => [...document.querySelectorAll('#studies-table tbody tr.studies-row')]"
            ".map((r) => r.dataset.path.split('/').pop()).sort()")
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1400, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            page.click('#studies-table tr.studies-row[data-path$="/beta"] .studies-pick input')
            page.click("#studies-archive")
            page.wait_for_selector(".studies-bar-said:not([hidden])")
            said = page.inner_text(".studies-bar-said")
            listed = page.evaluate(rows)
            on_disk = tags_of(tmp_path / "beta").get("archived")
            page.click('.studies-view-tab[data-view="archived"]')
            page.wait_for_selector('#studies-table tr.studies-row[data-path$="/beta"]')
            in_archived = page.evaluate(rows)
            page.click('#studies-table tr.studies-row[data-path$="/beta"] .studies-pick input')
            page.click("#studies-archive")
            page.wait_for_function("() => document.querySelector('.studies-bar-said')"
                                   ".textContent.includes('back in All studies')")
            after = page.evaluate(rows)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert "archived" in said and "does not keep Archive" not in said
    assert listed == ["alpha", "gamma"] and on_disk is True
    assert in_archived == ["beta"]
    assert after == ["alpha", "beta", "gamma"]
    assert "archived" not in tags_of(tmp_path / "beta")
