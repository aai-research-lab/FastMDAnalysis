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
