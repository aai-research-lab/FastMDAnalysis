"""The studies chosen on All studies have one bar, at the foot of the list.

Compare and Clear came beside the pointer once a second study was chosen
(user, 10-09), and no more than two could be. Choosing is for more than
comparing now: tagging and archiving many at once (user, 10-18: the mock's
bar "seems better"; "good ..... built it"). So while any study is chosen,
one bar at the foot of the list says how many and offers Compare (with
two), Tag (a tag on all of them or off all of them, a study's other tags
kept), Files (with one), Archive and Clear; Undo follows Tag and Archive.
The box in the list's head chooses every study shown, Shift and a click a
range, Esc clears and t opens Tag. Archived is kept beside the tags, in the
study's own `study_tags.json`: nothing of the study moves.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.test_the_workspace_says_its_studies import _study


def _tags(folder: Path) -> list[str]:
    return json.loads((folder / "study_tags.json").read_text(encoding="utf-8"))["tags"]


def test_archived_is_kept_beside_the_tags_until_it_is_undone(tmp_path):
    from fastmdxplora.gui.workspace import card_of
    from fastmdxplora.study_tags import add_tags, set_tags, tags_of

    study = _study(tmp_path / "wild_type")
    assert tags_of(study) == {"tags": [], "note": ""}
    said = set_tags(study, ["wild type"], "first try", archived=True)
    assert said["ok"] and said["archived"] is True
    kept = json.loads((study / "study_tags.json").read_text(encoding="utf-8"))
    when = kept["archived"]
    assert isinstance(when, str) and when.endswith("+00:00")
    assert card_of(study)["archived"] is True
    # Tags changed by a card, or added by an AI app, leave it archived, when
    # it was.
    set_tags(study, ["wild type", "JCIM Fig. 4"])
    add_tags(study, ["redo"])
    kept = json.loads((study / "study_tags.json").read_text(encoding="utf-8"))
    assert kept["archived"] == when and kept["note"] == "first try"
    assert set_tags(study, kept["tags"], archived=True)["archived"] is True
    assert json.loads((study / "study_tags.json").read_text())["archived"] == when
    # Taken back: the key is gone, the tags and note kept.
    said = set_tags(study, kept["tags"], archived=False)
    assert said["ok"] and "archived" not in said
    kept = json.loads((study / "study_tags.json").read_text(encoding="utf-8"))
    assert "archived" not in kept and kept["tags"] == ["wild type", "JCIM Fig. 4", "redo"]
    # Anything but True or False leaves it as it was.
    set_tags(study, [], archived="yes")
    assert "archived" not in tags_of(study) and "archived" not in card_of(study)


LOOK = """() => {
    const bar = document.getElementById('studies-bar');
    const r = bar.querySelector('.studies-bar-in:not([hidden])');
    const box = r ? r.getBoundingClientRect() : null;
    const rows = [...document.querySelectorAll('#studies-table tbody tr.studies-row')];
    const head = document.querySelector('#studies-table thead .studies-pick-all');
    return {
        shown: !bar.hidden && !!box,
        said: r ? r.innerText.replace(/\\s+/g, ' ').trim() : '',
        compare: document.getElementById('studies-compare').disabled,
        files: document.getElementById('studies-files-chosen').disabled,
        chosen: rows.filter((t) => t.querySelector('.studies-pick input').checked)
            .map((t) => t.dataset.path.split('/').pop()),
        lit: rows.filter((t) => t.classList.contains('is-chosen')).length,
        head: head ? [head.checked, head.indeterminate] : null,
        heard: document.getElementById('studies-chosen-heard').textContent,
        tagger: (document.querySelector('.studies-tagger') || {}).innerText || '',
        bottom: box ? Math.round(box.bottom) : null,
        height: innerHeight,
        tags: Object.fromEntries(rows.map((t) => [t.dataset.path.split('/').pop(),
            [...t.querySelectorAll('.studies-row-tag')].map((g) => g.textContent)])),
    };
}"""


def test_the_bar_chooses_tags_compares_and_clears(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.study_tags import set_tags

    names = ["alpha", "beta", "gamma", "delta"]
    for k, name in enumerate(names):
        _study(tmp_path / name, means={"rmsd": (0.1 + 0.01 * k, 0.004)},
               started=f"2026-09-0{k + 1}T10:00:00+00:00")
    set_tags(tmp_path / "alpha", ["wild type"], "kept as it is")
    set_tags(tmp_path / "beta", ["wild type", "JCIM Fig. 4"])
    session = start_dashboard_session(output=str(tmp_path / "alpha"), host="127.0.0.1", port=0)
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
            order = page.eval_on_selector_all(
                "#studies-table tbody tr.studies-row",
                "(all) => all.map((r) => r.dataset.path.split('/').pop())")

            def box(name):
                return page.locator(f'#studies-table tr.studies-row[data-path$="/{name}"] '
                                    ".studies-pick input")

            none = page.evaluate(LOOK)
            box(order[0]).click()
            one = page.evaluate(LOOK)
            box(order[1]).click()
            two = page.evaluate(LOOK)
            # Shift and a click: every study from the one chosen last.
            box(order[3]).click(modifiers=["Shift"])
            ranged = page.evaluate(LOOK)
            # The head's box: with every study shown chosen, none; again,
            # every one.
            page.click("#studies-table thead .studies-pick-all")
            cleared_by_head = page.evaluate(LOOK)
            page.click("#studies-table thead .studies-pick-all")
            all_by_head = page.evaluate(LOOK)
            page.click("#studies-table thead .studies-pick-all")
            box("alpha").click()
            box("beta").click()
            box("gamma").click()
            page.keyboard.press("t")
            page.wait_for_selector(".studies-tagger")
            opened = page.evaluate(LOOK)
            # "wild type" on 2 of 3: ticked, on all three. A new tag too.
            page.click('.studies-tagger-item:has-text("wild type") input')
            page.fill(".studies-tagger-new", "batch one")
            page.press(".studies-tagger-new", "Enter")
            asked = page.evaluate(LOOK)
            page.click('.studies-tagger-foot button:has-text("Apply")')
            page.wait_for_selector(".studies-bar-said:not([hidden])")
            tagged = page.evaluate(LOOK)
            on_disk = {name: _tags(tmp_path / name) for name in ("alpha", "beta", "gamma")}
            note = json.loads((tmp_path / "alpha" / "study_tags.json").read_text())["note"]
            page.click(".studies-bar-undo")
            page.wait_for_function("() => document.querySelector('.studies-bar-said')"
                                   ".textContent.includes('back as they were')")
            undone = {name: _tags(tmp_path / name) for name in ("alpha", "beta", "gamma")}
            page.click(".studies-bar-said .is-icon")
            # Esc clears them, the keyboard put back where it can be.
            page.keyboard.press("Escape")
            escaped = page.evaluate(LOOK)
            # Two chosen and c: the Compare page.
            box("alpha").click()
            box("delta").click()
            page.keyboard.press("c")
            page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
            compared = page.url
            page.go_back()
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            # On a phone, the bar across the window's foot.
            page.set_viewport_size({"width": 390, "height": 844})
            page.wait_for_timeout(300)
            box("beta").click()
            phone = page.evaluate(LOOK)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert not none["shown"]
    assert one["shown"] and one["said"].startswith("1 chosen")
    assert one["compare"] and not one["files"]
    assert one["heard"] == f"{order[0]} chosen." and one["head"] == [False, True]
    assert two["said"].startswith("2 chosen") and not two["compare"] and two["files"]
    assert two["heard"] == (f"{order[0]} and {order[1]} chosen. "
                            "Compare is in the bar under the list.")
    assert sorted(ranged["chosen"]) == sorted(names) and ranged["lit"] == 4
    assert ranged["head"] == [True, False] and ranged["compare"]
    assert not cleared_by_head["shown"] and cleared_by_head["chosen"] == []
    assert sorted(all_by_head["chosen"]) == sorted(names)
    assert all_by_head["said"].startswith("4 chosen")
    assert "Tags of the 3 chosen" in opened["tagger"]
    assert "wild type\n2 of 3" in opened["tagger"] and "JCIM Fig. 4\n1 of 3" in opened["tagger"]
    assert "wild type\nall 3" in asked["tagger"] and "batch one\nall 3" in asked["tagger"]
    assert tagged["said"].startswith("Tags changed on 3 studies.")
    assert on_disk == {"alpha": ["wild type", "batch one"],
                       "beta": ["wild type", "JCIM Fig. 4", "batch one"],
                       "gamma": ["wild type", "batch one"]}
    assert note == "kept as it is"
    assert tagged["tags"]["beta"] == ["wild type", "JCIM Fig. 4", "batch one"]
    assert undone == {"alpha": ["wild type"], "beta": ["wild type", "JCIM Fig. 4"], "gamma": []}
    assert not escaped["shown"] and escaped["chosen"] == []
    assert "#compare" in compared
    assert phone["shown"] and phone["bottom"] == phone["height"]


def test_beyond_loopback_nothing_that_posts_is_offered(tmp_path):
    """Bound beyond loopback the server refuses every POST but the Config's
    and `/api/compared`: Tag, Archive, Open and Files, which POST, and
    Compare were offered on a row and in the bar all the same. Its list is
    given here as the page would be given it (the route is loopback's)."""
    import re

    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.exploration import DashboardRuntime
    from fastmdxplora.gui.workspace import studies_in
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _serving

    for k, name in enumerate(("alpha", "beta")):
        _study(tmp_path / name, means={"rmsd": (0.1 + 0.01 * k, 0.004)},
               started=f"2026-09-0{k + 1}T10:00:00+00:00")
    listed = json.dumps(studies_in(tmp_path))
    runtime = DashboardRuntime(workspace_root=tmp_path / "alpha", exploration_root=tmp_path,
                               active_root=tmp_path / "alpha")
    look = """() => ({
        row: [...document.querySelectorAll('#studies-table tr.studies-row button')]
            .filter((b) => b.offsetParent !== null).map((b) => b.className),
        bar: [...document.querySelectorAll('#studies-bar button')]
            .filter((b) => b.offsetParent !== null).map((b) => b.id || b.className),
        tagger: !!document.querySelector('.studies-tagger'),
        page: location.hash})"""
    with _serving(runtime, allow_control=False) as url, sync_playwright() as pw:
        browser = pw.chromium.launch()
        page = browser.new_page(viewport={"width": 1440, "height": 900})
        page.set_default_timeout(60000)
        errors: list[str] = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.route(re.compile(r".*/api/studies(\?.*)?$"), lambda route: route.fulfill(
            status=200, content_type="application/json", body=listed))
        page.goto(url + "#studies", wait_until="domcontentloaded")
        page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
        page.wait_for_selector("#studies-table tbody tr.studies-row")
        bound = page.evaluate("() => document.body.dataset.beyondLoopback")
        # Nothing to choose: a bar of "2 chosen" and Clear alone was all
        # choosing gave there. Nor by keys.
        boxes = page.locator("#studies-table input[type=checkbox], .studies-pick-all").count()
        page.focus('#studies-table tr.studies-row[data-path$="/beta"] .studies-name-text')
        page.keyboard.press("x")
        page.keyboard.press("t")
        page.keyboard.press("c")
        seen = page.evaluate(look)
        bar = page.is_hidden("#studies-bar")
        page.click('[data-studies-view="cards"]')
        cards = page.locator(".study-card input[type=checkbox]").count()
        browser.close()
    assert errors == []
    assert bound == "true"
    assert boxes == 0 and cards == 0 and bar
    assert seen["row"] == [] and not seen["tagger"] and seen["page"] == "#studies"


def test_back_from_compare_is_where_the_two_were_chosen(tmp_path):
    """Compare opened from the bar, then Back: All studies where it was
    left, the two still chosen and in view, and the bar with them (kept
    from the tests of the offer beside the pointer this bar replaced)."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    for k in range(30):
        _study(tmp_path / f"study{k:02d}", means={"rmsd": (0.1 + 0.001 * k, 0.004)},
               started=f"2026-09-{k + 1:02d}T10:00:00+00:00")
    session = start_dashboard_session(output=str(tmp_path / "study00"), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 800})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            rows = page.locator("#studies-table tbody tr.studies-row .studies-pick input")
            rows.nth(24).scroll_into_view_if_needed()
            rows.nth(24).click()
            rows.nth(26).click()
            scrolled = page.evaluate("() => document.querySelector('.main').scrollTop")
            page.click("#studies-compare")
            page.wait_for_selector('.page[data-page="compare"][data-compared="yes"]')
            page.go_back()
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            # Given back once the list is read again.
            page.wait_for_function("() => document.querySelector('.main').scrollTop > 0",
                                   timeout=15000)
            page.wait_for_timeout(300)
            back = page.evaluate("() => document.querySelector('.main').scrollTop")
            chosen = page.evaluate("""() => [...document.querySelectorAll(
                '#studies-table tbody .studies-pick input')].filter((b) => b.checked).length""")
            box = rows.nth(26).bounding_box()
            bar = page.evaluate("() => !document.getElementById('studies-bar').hidden")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert scrolled > 0 and abs(back - scrolled) <= 2, (back, scrolled)
    assert chosen == 2 and bar
    assert 0 <= box["y"] <= 800 - box["height"], box


def test_no_rule_is_left_for_the_offer_beside_the_pointer():
    """The bar replaced Compare and Clear beside the pointer: the page's
    rules for that offer, left in `dashboard.css`, are gone with it."""
    from fastmdxplora.gui import server

    css = (Path(server.__file__).parent / "static" / "dashboard.css").read_text(encoding="utf-8")
    assert "studies-compare-bar" not in css


@pytest.mark.parametrize("width", [1440, 390])
def test_the_bar_leaves_the_last_row_in_view(tmp_path, width):
    """The bar was drawn up into the list from under it, over its last
    row, wherever the list ended in view. It has room of its own under the
    last row (on a phone, the page room under the bar across its foot)."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    for k, name in enumerate(("alpha", "beta", "gamma")):
        _study(tmp_path / name, means={"rmsd": (0.1 + 0.01 * k, 0.004)},
               started=f"2026-09-0{k + 1}T10:00:00+00:00")
    session = start_dashboard_session(output=str(tmp_path / "alpha"), host="127.0.0.1", port=0)
    look = """() => {
        const rows = [...document.querySelectorAll('#studies-table tbody tr.studies-row')];
        const last = rows[rows.length - 1].getBoundingClientRect();
        const bar = document.getElementById('studies-bar-in').getBoundingClientRect();
        return {lastBottom: last.bottom, barTop: bar.top, barBottom: bar.bottom, height: innerHeight};
    }"""
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": width, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            page.locator("#studies-table tbody tr.studies-row .studies-pick input").first.click()
            page.wait_for_selector("#studies-bar-in:not([hidden])")
            # As far down as the page goes.
            page.evaluate("""() => { document.querySelectorAll('*').forEach((n) => {
                if (n.scrollHeight > n.clientHeight + 1) n.scrollTop = n.scrollHeight; }); }""")
            page.wait_for_timeout(300)
            seen = page.evaluate(look)
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert seen["lastBottom"] <= seen["barTop"] + 1, seen
    assert seen["barBottom"] <= seen["height"] + 1, seen


def test_a_change_keeps_a_tag_kept_since_the_list_was_read(tmp_path):
    """Archive and the bar's Tag sent the tags as the page had read them,
    and the server put them in place of those kept: a tag an AI app added
    through `fastmdx mcp` since (or another window) was lost without a
    word. They send the change alone, applied to the tags as kept."""
    from fastmdxplora.study_tags import (A_RUN_IS_ARCHIVED_WITH_ITS_STUDY, add_tags,
                                         change_tags, set_tags, tags_of)
    from tests.test_the_workspace_says_its_studies import _replicas

    study = _study(tmp_path / "alpha")
    set_tags(study, ["reference"], "kept as it is")
    add_tags(study, ["from the agent"])
    said = change_tags(study, add=["batch one"], remove=["REFERENCE"])
    assert said["tags"] == ["from the agent", "batch one"] and said["note"] == "kept as it is"
    archived = change_tags(study, archived=True)
    assert archived["archived"] is True and archived["tags"] == ["from the agent", "batch one"]
    assert change_tags(study, add="x")["reason"] == "Tags are given as a list of words."
    # A run of a study is archived with its study, everywhere.
    run = _replicas(tmp_path / "reps") / "runs" / "seed1"
    (run / "manifest.json").write_text('{"phases": []}', encoding="utf-8")
    assert change_tags(run, archived=True)["reason"] == A_RUN_IS_ARCHIVED_WITH_ITS_STUDY
    assert set_tags(run, [], archived=True)["reason"] == A_RUN_IS_ARCHIVED_WITH_ITS_STUDY
    assert "archived" not in tags_of(run)


def test_the_bar_keeps_tags_kept_elsewhere_and_archives_no_run_alone(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.study_tags import add_tags, set_tags
    from tests.test_the_workspace_says_its_studies import _replicas

    for k, name in enumerate(("alpha", "beta")):
        _study(tmp_path / name, started=f"2026-09-0{k + 1}T10:00:00+00:00")
    reps = _replicas(tmp_path / "reps")
    for n in (1, 2, 3):
        (reps / "runs" / f"seed{n}" / "manifest.json").write_text(
            '{"phases": [{"name": "simulation", "status": "ok"}]}', encoding="utf-8")
    set_tags(tmp_path / "alpha", ["reference"])
    session = start_dashboard_session(output=str(tmp_path / "alpha"), host="127.0.0.1", port=0)

    def box(page, name):
        return page.locator(f'#studies-table tr.studies-row[data-path$="/{name}"] .studies-pick input')

    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            # Kept elsewhere once the list is read, as `fastmdx mcp` adds one.
            add_tags(tmp_path / "alpha", ["from the agent"])
            box(page, "alpha").click()
            page.keyboard.press("t")
            page.fill(".studies-tagger-new", "batch one")
            page.press(".studies-tagger-new", "Enter")
            page.click('.studies-tagger-foot button:has-text("Apply")')
            page.wait_for_selector(".studies-bar-said:not([hidden])")
            tagged = json.loads((tmp_path / "alpha" / "study_tags.json").read_text())["tags"]
            add_tags(tmp_path / "alpha", ["again"])
            page.click(".studies-bar-said .is-icon")
            page.click("#studies-archive")
            page.wait_for_function("() => document.querySelector('.studies-bar-said')"
                                   ".textContent.includes('archived')")
            archived = json.loads((tmp_path / "alpha" / "study_tags.json").read_text())
            # A run of a study chosen: Archive is not offered, and says why.
            page.click('#studies-table tr.studies-row[data-path$="/reps"] .studies-fold')
            page.wait_for_selector('#studies-table tr.studies-child[data-path$="/seed1"]')
            box(page, "seed1").click()
            run = page.evaluate("""() => { const b = document.getElementById('studies-archive');
                return [b.disabled, b.title]; }""")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert tagged == ["reference", "from the agent", "batch one"]
    assert archived["tags"] == ["reference", "from the agent", "batch one", "again"]
    assert archived["archived"]
    assert run == [True, "A run of a study is archived with its study: archive the study."]


def test_a_chip_ticked_keeps_the_keyboard_on_its_box(tmp_path):
    """Ticking a State chip's third box put the keyboard on its first: a
    keyboard user ticking two had to go back down each time."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    _study(tmp_path / "alpha")
    _study(tmp_path / "beta", failed="simulation.run.stopped")
    _study(tmp_path / "gamma", failed="simulation.run.failed")
    session = start_dashboard_session(output=str(tmp_path / "alpha"), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.wait_for_selector("body.state-ready")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            page.click('.studies-chip[data-key="states"] > summary')
            boxes = page.locator('.studies-chip[data-key="states"] input')
            last = boxes.count() - 1
            boxes.nth(last).focus()
            page.keyboard.press("Space")
            on = page.evaluate("""() => [document.activeElement.dataset.value,
                [...document.querySelectorAll('.studies-chip[data-key="states"] input')]
                    .map((b) => b.dataset.value)]""")
            browser.close()
    finally:
        session.server.shutdown()
    assert on[0] == on[1][-1] and on[0] != on[1][0]
