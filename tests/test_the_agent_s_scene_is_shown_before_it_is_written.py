"""The Agent's scene is shown in the Viewer before anything is written.

To see a scene the Agent proposed, it had to be written first as
`scenes/<name>.mvsx` in the study (Write this scene). The Viewer now takes a
scene's parts (`showScene`): its frame, representation, colour,
superposition and highlighted atoms, read by MDTraj on the server as every
selection is, with their labels, the camera kept and nothing written. A part
the study cannot show is said in the notice over the molecule and the rest
shown; Back to my view puts back the view from before the first scene shown.
"""

from __future__ import annotations

import pytest

pytest.importorskip("mdtraj")

from tests.test_the_frames_are_fitted_to_a_structure import made_study  # noqa: E402

import re  # noqa: E402
import time  # noqa: E402

VIEWER = "window.FastMDXMoleculeViewer"
STATE = VIEWER + ".STATE"
HIGHLIGHT = 0xE69F00


@pytest.fixture(scope="module")
def browser():
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    pw = sync_playwright().start()
    browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
    try:
        yield browser
    finally:
        browser.close()
        pw.stop()


def _opened(browser, tmp_path_factory, label):
    from fastmdxplora.gui.server import start_dashboard_session

    study = made_study(tmp_path_factory.mktemp(label) / "study")
    (study / "setup" / "topology.pdb").write_text(
        (study / "simulation" / "trajectory_topology.pdb").read_text(encoding="utf-8"),
        encoding="utf-8")
    session = start_dashboard_session(output=str(study), host="127.0.0.1", port=0)
    page = browser.new_page(viewport={"width": 1400, "height": 900})
    try:
        page.set_default_timeout(60000)
        errors: list[str] = []
        page.on("pageerror", lambda error: errors.append(str(error)))
        page.goto(session.url + "#viewer", wait_until="domcontentloaded")
        if not page.evaluate("() => !!document.createElement('canvas').getContext('webgl')"):
            pytest.skip("this browser has no WebGL, so the viewer cannot render")
        page.wait_for_function(f"() => {VIEWER} && {STATE}.model && {STATE}.engine")
        yield {"page": page, "study": study, "errors": errors}
    finally:
        page.close()
        session.server.shutdown()


@pytest.fixture(scope="module")
def viewer(browser, tmp_path_factory):
    yield from _opened(browser, tmp_path_factory, "scene")


@pytest.fixture(scope="module")
def pieces(browser, tmp_path_factory):
    """A study whose frames come two at a time, as a big study's do."""
    from fastmdxplora.gui import trajectory_frames

    was = trajectory_frames.PIECE_ATOM_FRAMES
    trajectory_frames.PIECE_ATOM_FRAMES = 130
    try:
        yield from _opened(browser, tmp_path_factory, "pieces")
    finally:
        trajectory_frames.PIECE_ATOM_FRAMES = was


def _forNow(page):
    """The highlight drawn for now, as the engine is given it."""
    return page.evaluate("""() => (window.FastMDXSelections.forScene() || [])
        .filter((one) => one.colour === %d)
        .map((one) => ({atoms: one.atoms.length, labelled: one.labelled,
                        representation: one.representation}))""" % HIGHLIGHT)


def _written(study):
    return sorted(str(path.relative_to(study)) for path in study.rglob("*")
                  if path.name in ("viewer_views.json", "viewer_selections.json")
                  or "scenes" in path.parts)


def test_a_scene_is_shown_and_its_view_put_back(viewer):
    page, study = viewer["page"], viewer["study"]
    before = page.evaluate(f"() => {VIEWER}.viewNow()")
    shown = page.evaluate(f"""() => {VIEWER}.showScene({{frame: 2, representation: 'sticks',
        colour: 'chain', highlight: 'resSeq 3 to 5', labels: true}})""")
    now = page.evaluate(f"""() => ({{representation: {STATE}.representation,
        colour: {STATE}.colorMode, of: {STATE}.model.of, frame: {STATE}.engine.frame(),
        rep: document.getElementById('viewer-rep').value,
        notice: document.getElementById('viewer-scene-notice').hidden,
        said: document.getElementById('viewer-scene-said').textContent}})""")
    highlighted = _forNow(page)
    page.click("#viewer-scene-back")
    page.wait_for_function("() => document.getElementById('viewer-scene-notice').hidden"
                           f" && {STATE}.model && {STATE}.model.of !== 'frames'")
    page.wait_for_function(f"() => {STATE}.representation === '{before['representation']}'")
    back = page.evaluate(f"""() => ({{representation: {STATE}.representation,
        colour: {STATE}.colorMode, of: {STATE}.model.of}})""")
    assert shown == {"shown": True, "notes": []}
    assert now == {"representation": "sticks", "colour": "chain", "of": "frames", "frame": 2,
                   "rep": "sticks", "notice": False, "said": "The Agent's scene, not written."}
    # Three alanines of five atoms each, orange sticks with their labels.
    assert highlighted == [{"atoms": 15, "labelled": True, "representation": "sticks"}]
    assert back == {"representation": before["representation"], "colour": before["colour"],
                    "of": "structure"}
    assert _forNow(page) == []
    assert _written(study) == []
    assert viewer["errors"] == []


def test_what_the_study_cannot_show_is_said_and_the_rest_shown(viewer):
    page, study = viewer["page"], viewer["study"]
    shown = page.evaluate(f"""() => {VIEWER}.showScene({{frame: 99, colour: 'result:rmsf',
        representation: 'spacefill', highlight: 'resname and', labels: true}})""")
    said = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    representation = page.evaluate(f"() => {STATE}.representation")
    of = page.evaluate(f"() => {STATE}.model.of")
    # The structure shown: no frame of the frames is said over it, and the
    # notice is read out through the page's live region.
    overlay = page.evaluate("() => document.getElementById('overlay-frame').textContent")
    read_out = page.evaluate("() => document.getElementById('sr-live').textContent")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert overlay == ""
    assert read_out == said
    assert shown["shown"] is True
    # The study may have an RMSF of its own; none by residue to colour by.
    assert shown["notes"][0] == ("This study has no rmsf result by residue to colour by, so the "
                                 "colours are as they were.")
    assert shown["notes"][1] == "The study has no frame 99 to show; the structure is shown."
    assert shown["notes"][2].startswith("The atoms of resname and are not highlighted: ")
    assert len(shown["notes"]) == 3
    assert said == " ".join(["The Agent's scene, not written.", *shown["notes"]])
    assert representation == "spacefill"
    assert of == "structure"
    assert _forNow(page) == []
    assert _written(study) == []
    assert viewer["errors"] == []


def test_a_second_scene_goes_back_to_the_view_before_the_first(viewer):
    page = viewer["page"]
    before = page.evaluate(f"() => {VIEWER}.viewNow()")
    page.evaluate(f"() => {VIEWER}.showScene({{representation: 'sticks'}})")
    page.evaluate(f"() => {VIEWER}.showScene({{representation: 'lines', highlight: 'resSeq 1'}})")
    highlighted = _forNow(page)
    page.evaluate(f"() => {VIEWER}.showScene({{representation: 'surface'}})")
    after_plain = _forNow(page)
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    representation = page.evaluate(f"() => {STATE}.representation")
    assert highlighted == [{"atoms": 5, "labelled": False, "representation": "sticks"}]
    # A scene with no highlight leaves none from the one before.
    assert after_plain == []
    assert representation == before["representation"]
    assert viewer["errors"] == []


def test_written_says_it_is_kept(viewer):
    page = viewer["page"]
    page.evaluate(f"() => {VIEWER}.showScene({{representation: 'sticks'}})")
    page.evaluate(f"() => {VIEWER}.sceneWritten('Core')")
    said = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    hidden = page.evaluate("() => document.getElementById('viewer-scene-notice').hidden")
    assert said == "The Agent's scene, written as scenes/Core.mvsx."
    assert hidden is True


def test_a_superposition_is_shown_only_as_a_written_scene_shows_it(viewer):
    page = viewer["page"]
    # No frame named: the structure is shown, so nothing is superposed.
    plain = page.evaluate(f"() => {VIEWER}.showScene({{superposed: 'backbone'}})")
    plain_of = page.evaluate(f"() => {STATE}.model.of")
    # The pocket is not offered in a study with no ligand.
    pocket = page.evaluate(f"() => {VIEWER}.showScene({{frame: 1, superposed: 'pocket'}})")
    pocket_now = page.evaluate(f"() => [{STATE}.model.of, {STATE}.superposed || 'none']")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert plain == {"shown": True, "notes": [
        "The frames are not shown, so nothing is superposed on the backbone."]}
    assert plain_of == "structure"
    assert pocket == {"shown": True, "notes": [
        "This structure has no ligand, so nothing is superposed on the pocket."]}
    assert pocket_now == ["frames", "none"]
    assert viewer["errors"] == []


def test_a_scene_overtaken_leaves_nothing_of_itself(viewer):
    page = viewer["page"]

    # The second scene, with no highlight, is asked for as the first asks
    # for its frames superposed, which come late.
    answers = page.evaluate(f"""async () => {{
        const real = window.fetch;
        let second = null;
        window.fetch = (url, ...rest) => {{
            if (!second && String(url).includes('/api/frames-superposed')) {{
                second = {VIEWER}.showScene({{representation: 'lines'}});
                return new Promise((done) => setTimeout(done, 1500)).then(() => real(url, ...rest));
            }}
            return real(url, ...rest);
        }};
        try {{
            const first = await {VIEWER}.showScene({{frame: 3, superposed: 'backbone',
                representation: 'spacefill', highlight: 'resSeq 1 to 4', labels: true}});
            return [first, await second];
        }} finally {{
            window.fetch = real;
        }} }}""")
    now = page.evaluate(f"() => {STATE}.representation")
    shown = page.evaluate(f"() => [{STATE}.mode, {STATE}.model.of]")
    highlighted = _forNow(page)
    said = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert answers[0]["shown"] is False
    assert answers[1] == {"shown": True, "notes": []}
    assert now == "lines"
    # The second scene names no frame: the first's frames are not left.
    assert shown[0] != "playback" and shown[1] == "structure"
    # Nothing of the first scene's highlight is left under the second.
    assert highlighted == []
    assert said == "The Agent's scene, not written."
    assert viewer["errors"] == []


def test_back_while_a_scene_is_shown_leaves_the_view_before(viewer):
    page = viewer["page"]
    before = page.evaluate(f"() => {VIEWER}.viewNow()")
    page.evaluate(f"""() => Promise.all([
        {VIEWER}.showScene({{frame: 3, representation: 'sticks', highlight: 'resSeq 2 to 4'}}),
        {VIEWER}.backToMyView()])""")
    after = page.evaluate(f"""() => ({{representation: {STATE}.representation,
        of: {STATE}.model.of,
        notice: document.getElementById('viewer-scene-notice').hidden}})""")
    assert after == {"representation": before["representation"], "of": "structure",
                     "notice": True}
    assert _forNow(page) == []
    assert viewer["errors"] == []


def test_back_asked_as_a_scene_is_shown_waits_for_it_and_undoes_it(viewer):
    """Back asked while a scene's frames were being superposed: run at
    once, it put back the view before, and the scene, finishing after it,
    left the frames superposed under the structure Back showed."""
    page = viewer["page"]
    page.evaluate("() => document.querySelector('[data-action=\"live-toggle\"]').click()")
    page.wait_for_function(f"() => {STATE}.mode !== 'playback' && {STATE}.model"
                           f" && {STATE}.model.of === 'structure'")
    answers = page.evaluate(f"""async () => {{
        const real = window.fetch;
        let back = null;
        window.fetch = (url, ...rest) => {{
            if (!back && String(url).includes('/api/frames-superposed')) {{
                back = {VIEWER}.backToMyView();
                return new Promise((done) => setTimeout(done, 1500)).then(() => real(url, ...rest));
            }}
            return real(url, ...rest);
        }};
        try {{
            const shown = await {VIEWER}.showScene({{frame: 3, superposed: 'backbone'}});
            return [shown.shown, await back];
        }} finally {{
            window.fetch = real;
        }} }}""")
    page.wait_for_timeout(1000)
    left = page.evaluate(f"""() => [{STATE}.superposed || 'none',
        document.getElementById('traj-superpose').value, {STATE}.mode !== 'playback',
        {STATE}.model.of]""")
    assert answers == [False, True]
    assert left == ["none", "none", True, "structure"]
    assert viewer["errors"] == []


def test_focus_goes_to_the_molecule_after_back(viewer):
    page = viewer["page"]
    page.evaluate(f"() => {VIEWER}.showScene({{representation: 'sticks'}})")
    page.focus("#viewer-scene-back")
    page.keyboard.press("Enter")
    page.wait_for_function("() => document.getElementById('viewer-scene-notice').hidden")
    page.wait_for_function("() => document.activeElement"
                           " && document.activeElement.id === 'viewer-canvas'")
    assert viewer["errors"] == []


def test_a_scene_with_nothing_to_show_is_said_and_not_shown(viewer):
    page = viewer["page"]
    empty = page.evaluate(f"() => {VIEWER}.showScene(null)")
    hidden = page.evaluate("() => document.getElementById('viewer-scene-notice').hidden")
    odd = page.evaluate(f"() => {VIEWER}.showScene({{frame: '1e0'}})")
    long = page.evaluate(f"() => {VIEWER}.showScene({{highlight: 'resSeq 1 or '.repeat(30) + 'resSeq 2'}})")
    named = page.evaluate(f"() => {VIEWER}.showScene({{representation: 'Ball and stick', colour: 'result:'}})")
    representation = page.evaluate(f"() => {STATE}.representation")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert empty == {"shown": False, "notes": ["The scene has no parts to show."]}
    assert hidden is True
    assert odd["notes"] == ["The study has no frame 1e0 to show; the structure is shown."]
    assert long["notes"] == ["The atoms to highlight are named in more than 200 characters,"
                             " so none are highlighted."]
    assert named["notes"] == ["The scene names no result to colour by, so the colours are as they were."]
    assert representation == "ballAndStick"
    assert viewer["errors"] == []


def test_the_scene_s_frame_is_held_while_playback_runs(viewer):
    page = viewer["page"]
    page.evaluate(f"""async () => {{ await {VIEWER}.movie.frames();
        {STATE}.playbackLoop = true;
        document.querySelector('[data-action="play-trajectory"]').click(); }}""")
    page.wait_for_function(f"() => {STATE}.playbackPlaying")
    page.evaluate(f"() => {VIEWER}.showScene({{frame: 2}})")
    time.sleep(1.5)
    held = page.evaluate(f"""() => [{STATE}.playbackPlaying,
        document.getElementById('traj-slider').value]""")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    # Back to the frames playing, as they were before the scene.
    playing = page.evaluate(f"() => {STATE}.playbackPlaying")
    page.evaluate("() => document.querySelector('[data-action=\"play-trajectory\"]').click()")
    page.wait_for_function(f"() => !{STATE}.playbackPlaying")
    assert held == [False, "2"]
    assert playing is True
    assert viewer["errors"] == []


def test_back_reads_on_the_notice_in_every_scheme(viewer):
    page = viewer["page"]
    page.evaluate(f"() => {VIEWER}.showScene({{representation: 'sticks'}})")
    colours = []
    for scheme in ("light", "dark"):
        page.evaluate("(name) => { document.documentElement.dataset.theme = name;"
                      " document.body.dataset.theme = name; }", scheme)
        page.hover("#viewer-scene-back")
        page.wait_for_timeout(500)
        colours.append(page.evaluate("""() => getComputedStyle(
            document.getElementById('viewer-scene-back')).color"""))
    page.mouse.move(0, 0)
    page.evaluate(f"() => {VIEWER}.backToMyView()")

    def light(colour):
        red, green, blue = (int(part) for part in re.findall(r"\d+", colour)[:3])
        return 0.2126 * red + 0.7152 * green + 0.0722 * blue > 140

    # Light text on the notice's dark ground, hovered, in both schemes.
    assert [light(colour) for colour in colours] == [True, True], colours
    assert viewer["errors"] == []


def test_the_ground_and_the_person_s_changes_over_a_scene(viewer):
    """The ground chosen as the page's scheme was set to white or dark by
    a scene and by Back, and no longer followed the scheme; Now and the
    ground chosen over a scene were not said as changes."""
    page = viewer["page"]
    page.evaluate("() => { const ground = document.getElementById('viewer-ground');"
                  " ground.value = 'scheme'; ground.dispatchEvent(new Event('change')); }")
    page.evaluate(f"() => {VIEWER}.showScene({{frame: 2}})")
    during = page.evaluate(f"() => {STATE}.ground")
    page.evaluate("() => document.querySelector('[data-action=\"live-toggle\"]').click()")
    page.wait_for_function(f"() => {STATE}.mode !== 'playback'")
    said = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    after = page.evaluate(f"() => {STATE}.ground")
    assert (during, after) == ("scheme", "scheme")
    assert said == "The Agent's scene, not written, changed here since."
    assert viewer["errors"] == []


def test_the_highlight_is_read_again_when_the_frames_are_shown(viewer):
    page = viewer["page"]
    page.evaluate(f"() => {VIEWER}.showScene({{highlight: 'resSeq 3 to 5'}})")
    on_structure = _forNow(page)
    page.evaluate(f"() => {VIEWER}.movie.frames()")
    page.wait_for_function(f"() => {STATE}.model.of === 'frames'")
    page.wait_for_function("() => (window.FastMDXSelections.forScene() || []).length")
    on_frames = _forNow(page)
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert on_structure == on_frames == [{"atoms": 15, "labelled": False, "representation": "sticks"}]
    assert viewer["errors"] == []


def test_a_change_made_over_the_scene_is_said(viewer):
    page = viewer["page"]
    before = page.evaluate(f"() => {VIEWER}.viewNow()")
    page.evaluate(f"() => {VIEWER}.showScene({{representation: 'sticks'}})")
    page.evaluate("""() => { const rep = document.getElementById('viewer-rep');
        rep.value = 'lines'; rep.dispatchEvent(new Event('change')); }""")
    said = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    page.evaluate(f"() => {VIEWER}.sceneWritten('Kept')")
    written = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert said == "The Agent's scene, not written, changed here since."
    assert written == "The Agent's scene, written as scenes/Kept.mvsx, changed here since."
    assert page.evaluate(f"() => {STATE}.representation") == before["representation"]
    assert viewer["errors"] == []


def test_the_scene_s_notes_stay_once_it_is_written(viewer):
    page = viewer["page"]
    shown = page.evaluate(f"() => {VIEWER}.showScene({{representation: 'ribbon'}})")
    page.evaluate(f"() => {VIEWER}.sceneWritten('Core')")
    said = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert shown["notes"] == ["The Viewer has no representation called ribbon."]
    assert said == ("The Agent's scene, written as scenes/Core.mvsx."
                    " The Viewer has no representation called ribbon.")


def test_the_frames_are_fitted_as_the_written_scene_fits_them(viewer):
    """Fitted to the deposited structure and smoothed by the person, the
    frames were shown so under a scene that, written, fits them to the
    first frame unsmoothed."""
    page = viewer["page"]
    page.evaluate(f"""() => {{ {STATE}.superposedTo = 'deposited'; {STATE}.smoothedOver = 3;
        {STATE}.pocketCutoff = 8; }}""")
    page.evaluate(f"() => {VIEWER}.showScene({{frame: 2, superposed: 'backbone'}})")
    fitted = page.evaluate(f"() => [{STATE}.superposedTo, {STATE}.smoothedOver, {STATE}.pocketCutoff]")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert fitted == ["first", 1, 5]
    assert viewer["errors"] == []


def test_a_frame_stepped_to_over_the_scene_is_said(viewer):
    page = viewer["page"]
    page.evaluate(f"() => {VIEWER}.showScene({{frame: 2}})")
    unchanged = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    page.evaluate(f"() => {VIEWER}.movie.showFrame(3)")
    said = page.evaluate("() => document.getElementById('viewer-scene-said').textContent")
    page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert unchanged == "The Agent's scene, not written."
    assert said == "The Agent's scene, not written, changed here since."
    assert viewer["errors"] == []


def test_the_notice_on_a_phone_leaves_the_colour_bar_clear(viewer):
    page = viewer["page"]
    page.set_viewport_size({"width": 390, "height": 844})
    try:
        page.evaluate(f"() => {VIEWER}.showScene({{representation: 'sticks'}})")
        # However much it says.
        page.evaluate("""() => { document.getElementById('viewer-scene-said').textContent =
            'The atoms are not highlighted. '.repeat(60); }""")
        boxes = page.evaluate("""() => {
            const box = (id) => document.getElementById(id).getBoundingClientRect();
            const notice = box('viewer-scene-notice'), frame = box('viewer-canvas-frame');
            const back = box('viewer-scene-back');
            return {bottom: frame.bottom - notice.bottom,
                    back: back.top >= notice.top && back.bottom <= notice.bottom + 0.5};
        }""")
        page.evaluate(f"() => {VIEWER}.backToMyView()")
    finally:
        page.set_viewport_size({"width": 1400, "height": 900})
    # The colour bar takes the 100 px at the foot of the canvas, and Back to
    # my view stays in sight however much is said above it.
    assert boxes == {"bottom": boxes["bottom"], "back": True}
    assert boxes["bottom"] >= 99
    assert viewer["errors"] == []


def test_the_person_s_own_view_ends_the_scene(viewer):
    page, study = viewer["page"], viewer["study"]
    # The structure shown, not the frames, as the view saved.
    page.evaluate("() => document.querySelector('[data-action=\"live-toggle\"]').click()")
    page.wait_for_function(f"() => {STATE}.mode !== 'playback' && {STATE}.model"
                           f" && {STATE}.model.of === 'structure'")
    page.evaluate("""() => { document.getElementById('viewer-view-save').click();
        document.getElementById('viewer-view-name').value = 'Mine';
        document.getElementById('viewer-view-keep').click(); }""")
    page.wait_for_function("() => [...document.getElementById('viewer-views').options]"
                           ".some((o) => o.value === 'Mine')")
    page.evaluate(f"""() => {VIEWER}.showScene({{frame: 3, representation: 'sticks',
        highlight: 'resSeq 1'}})""")
    page.evaluate("""() => { const views = document.getElementById('viewer-views');
        views.value = 'Mine'; views.dispatchEvent(new Event('change')); }""")
    hidden = page.evaluate("() => document.getElementById('viewer-scene-notice').hidden")
    # Saved on the structure: the frames the scene loaded are left.
    page.wait_for_function(f"() => {STATE}.mode !== 'playback' && {STATE}.model"
                           f" && {STATE}.model.of === 'structure'")
    back = page.evaluate(f"() => {VIEWER}.backToMyView()")
    assert hidden is True
    assert _forNow(page) == []
    assert back is False
    assert _written(study) == ["viewer_views.json"]
    assert viewer["errors"] == []


def test_a_saved_view_chosen_as_a_scene_is_shown_is_the_one_left(viewer):
    """A saved view chosen while a scene was still being shown: the scene,
    finishing after it, left its own camera over the view chosen."""
    page = viewer["page"]
    camera = "{STATE}.engine && {VIEWER}.viewNow().camera.position.map((x) => Math.round(x * 10) / 10)"
    page.evaluate(f"""async () => {{ const view = {VIEWER}.viewNow();
        view.representation = 'lines';
        view.camera.position = view.camera.position.map((x) => x * 1.5);
        await {VIEWER}.showView(view); }}""")
    page.evaluate("""() => { document.getElementById('viewer-view-save').click();
        document.getElementById('viewer-view-name').value = 'Far';
        document.getElementById('viewer-view-keep').click(); }""")
    page.wait_for_function("() => [...document.getElementById('viewer-views').options]"
                           ".some((o) => o.value === 'Far')")
    far = page.evaluate(f"() => {camera.format(STATE=STATE, VIEWER=VIEWER)}")
    page.evaluate(f"""async () => {{ const view = {VIEWER}.viewNow();
        view.representation = 'sticks';
        view.camera.position = view.camera.position.map((x) => x / 1.5);
        await {VIEWER}.showView(view);
        document.getElementById('viewer-views').value = ''; }}""")
    page.evaluate(f"""async () => {{
        const real = window.fetch;
        let chosen = false;
        window.fetch = (url, ...rest) => {{
            if (!chosen && String(url).includes('/api/frames-superposed')) {{
                chosen = true;
                setTimeout(() => {{ const views = document.getElementById('viewer-views');
                    views.value = 'Far'; views.dispatchEvent(new Event('change')); }}, 50);
                return new Promise((done) => setTimeout(done, 1500)).then(() => real(url, ...rest));
            }}
            return real(url, ...rest);
        }};
        try {{
            return await {VIEWER}.showScene({{frame: 3, superposed: 'backbone',
                representation: 'spacefill', highlight: 'resSeq 1'}});
        }} finally {{
            window.fetch = real;
        }} }}""")
    page.wait_for_function(f"() => {STATE}.representation === 'lines'")
    page.wait_for_timeout(1500)
    left = page.evaluate(f"() => [{STATE}.representation, {camera.format(STATE=STATE, VIEWER=VIEWER)}]")
    assert left == ["lines", far]
    assert _forNow(page) == []
    assert viewer["errors"] == []


def test_a_study_changed_under_a_scene_leaves_nothing_of_it(viewer):
    page = viewer["page"]
    first = page.evaluate(f"() => {VIEWER}.showScene({{highlight: 'resSeq 2'}})")
    answer = page.evaluate(f"""async () => {{
        const shown = {VIEWER}.showScene({{frame: 3, highlight: 'resSeq 1'}});
        await new Promise((done) => setTimeout(done, 0));
        window.dispatchEvent(new CustomEvent('dashboard:run-changed', {{detail: {{first: false}}}}));
        return shown; }}""")
    left = page.evaluate("""() => [window.FastMDXSelections.state.forNow,
        document.getElementById('viewer-scene-notice').hidden]""")
    assert first == {"shown": True, "notes": []}
    assert answer["shown"] is False
    # Neither scene's highlight is left, nor its notice.
    assert left == [None, True]


def test_back_asked_as_the_study_changes_does_nothing_to_the_next(viewer):
    """Back went on over the next study: live updates paused before the
    scene were paused again there, and its view was said to be back."""
    page = viewer["page"]
    page.wait_for_function(f"() => {STATE}.model && {STATE}.engine && {STATE}.mode !== 'playback'")
    page.evaluate(f"() => {{ {STATE}.liveUpdates = false; }}")
    page.evaluate(f"() => {VIEWER}.showScene({{frame: 2, representation: 'sticks'}})")
    back = page.evaluate(f"""async () => {{
        const back = {VIEWER}.backToMyView();
        await new Promise((done) => setTimeout(done, 0));
        window.dispatchEvent(new CustomEvent('dashboard:run-changed', {{detail: {{first: false}}}}));
        return back; }}""")
    left = page.evaluate(f"""() => [{STATE}.liveUpdates,
        document.getElementById('sr-live').textContent]""")
    assert back is False
    assert left[0] is True
    assert left[1] != "Back to the view before the Agent's scene."


def _fresh(page):
    """The Viewer opened again, no frame loaded yet."""
    page.reload(wait_until="domcontentloaded")
    page.wait_for_function(f"() => {VIEWER} && {STATE}.model && {STATE}.engine")


def _held(page, *, first_too=False):
    """The pieces of the frames held, unanswered, until let go."""
    held = []

    def hold(route):
        if first_too or "k=0" not in route.request.url:
            held.append(route)
        else:
            route.continue_()

    page.route("**/structure/frames-piece.xtc*", hold)
    return held


def _let_go(page, held, *, only=None):
    from playwright.sync_api import Error

    for route in [one for one in held if only is None or only in one.request.url]:
        held.remove(route)
        try:
            route.continue_()
        except Error:
            pass  # answered as the page let its route go


def _until(page, test):
    for _ in range(600):
        if test():
            return
        page.wait_for_timeout(100)
    raise AssertionError("not reached")


def test_back_does_not_wait_for_the_frames_still_to_come(pieces):
    page = pieces["page"]
    _fresh(page)
    held = _held(page)
    try:
        page.evaluate(f"() => {{ window.__shown = {VIEWER}.showScene({{frame: 4}}); }}")
        # The first piece is here; the scene waits for the one with frame 4.
        _until(page, lambda: held)
        back = page.evaluate(f"() => {VIEWER}.backToMyView()")
        still = [len(held) > 0, page.evaluate(f"() => {STATE}.playbackFrames < 5")]
        shown = page.evaluate("() => window.__shown")
        hidden = page.evaluate("() => document.getElementById('viewer-scene-notice').hidden")
    finally:
        _let_go(page, held)
        page.unroute("**/structure/frames-piece.xtc*")
    assert still == [True, True]
    assert back is True
    assert shown["shown"] is False
    assert hidden is True
    assert pieces["errors"] == []


def test_a_scene_ended_as_its_frames_load_says_nothing_after(pieces):
    """A saved view chosen while the scene's first piece of frames was on
    its way: the scene, coming out of the load, said it was waiting for its
    frame over the view chosen, with Back to my view, and left it there."""
    page = pieces["page"]
    _fresh(page)
    held = _held(page, first_too=True)
    try:
        page.evaluate(f"() => {{ window.__shown = {VIEWER}.showScene({{frame: 4}}); }}")
        _until(page, lambda: held)
        page.evaluate(f"() => {{ window.__ended = {VIEWER}.endTheScene(null); }}")
        _let_go(page, held, only="k=0")
        shown = page.evaluate("() => window.__shown")
        page.evaluate("() => window.__ended")
        notice = page.evaluate("""() => [document.getElementById('viewer-scene-notice').hidden,
            document.getElementById('viewer-scene-said').textContent]""")
    finally:
        _let_go(page, held)
        page.unroute("**/structure/frames-piece.xtc*")
    assert shown["shown"] is False
    assert notice[0] is True, notice
    assert pieces["errors"] == []


def test_a_frame_still_arriving_is_waited_for(pieces):
    """The page asks how many frames there are every few seconds and was
    told all of them: the wait ended, and frame 4 was said shown while the
    engine, holding fewer, showed its last."""
    page = pieces["page"]
    _fresh(page)
    held = _held(page)
    try:
        page.evaluate(f"""() => {{ window.__done = false;
            window.__shown = {VIEWER}.showScene({{frame: 4}});
            window.__shown.then(() => {{ window.__done = true; }}); }}""")
        _until(page, lambda: held)
        # Long enough for the page to ask for the frames' count again.
        page.wait_for_timeout(4500)
        waiting = page.evaluate("() => !window.__done")
    finally:
        _let_go(page, held)
        page.unroute("**/structure/frames-piece.xtc*")
        _let_go(page, held)
    shown = page.evaluate("() => window.__shown")
    frame = page.evaluate(f"() => {STATE}.engine.frame()")
    assert waiting is True
    assert shown == {"shown": True, "notes": []}
    assert frame == 4
    assert pieces["errors"] == []
