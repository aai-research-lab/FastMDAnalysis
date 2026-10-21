"""All studies opens on its list, each study a row that shows its molecule.

The page opened on two cards a row, each a plot of the study's first
measure, and its Table said a study's name, system, state, dates and its
first three means. A workspace of forty studies is read as rows, and a
study is known by its molecule before its RMSD (user, 10-18: "why not just
show the molecule?"). So the list comes first: views as tabs (All,
Running, Needs more sampling, Not started, and the person's own, kept in
this browser), chips for state, system and tags, rows grouped (replicas
folded under their study, folders where nothing ran gathered at the foot),
and in each row the molecule, where the study stands (a run going, how far
and how long it has left), the ns simulated and every mean, a filled dot
where it is determined and a ring where it is not, with why. A study of
replicas gives the mean of its replicas' means, its error their spread
over the square root of their number, as the Compare page does; and where
some replicas could not determine their own mean, the spread, not an
error. A trajectory analysed with no run of its own says "analysed", not
"not started". The cards are a gallery of the same molecules. Keys: / to
search, j and k to move, x to choose, m for every mean.
"""

from __future__ import annotations

import json
import math

import pytest

from tests.test_a_study_without_a_figure_shows_its_backbone import _helix, _pdb
from tests.test_the_workspace_says_its_studies import _replicas, _study


def _live(folder, *, status="running", current=6000, elapsed=381.555):
    from datetime import datetime, timedelta, timezone

    (folder / "simulation").mkdir(parents=True, exist_ok=True)
    (folder / "simulation" / "live_status.json").write_text(json.dumps({
        "stage": "production", "status": status, "current_step": current,
        "total_planned_steps": 12000, "timestep_fs": 2.0, "elapsed_wall_time_s": elapsed,
        "nvt_steps_planned": 1000, "npt_steps_planned": 1000, "production_steps_planned": 10000,
        "stage_states": {"production": "running" if status == "running" else "completed"},
    }), encoding="utf-8")
    # Production's own samples: 4,000 steps in 200 s, 20 steps a second.
    # The run's whole time (381.555 s for 6,000 steps) has setup in it.
    start = datetime(2026, 10, 8, 12, 0, tzinfo=timezone.utc)
    lines = ["timestamp,stage,step,simulation_time_ns,temperature,speed"]
    for k in range(5):
        lines.append(f"{(start + timedelta(seconds=50 * k)).isoformat()},Production,"
                     f"{2000 + 1000 * k},{(2000 + 1000 * k) * 2e-6},300,1.0")
    (folder / "simulation" / "live_metrics.csv").write_text("\n".join(lines) + "\n")
    return folder


def _members(folder, means, *, short=()):
    rows = [{"member": f"seed{n}", "mean": m,
             "standard_error": None if f"seed{n}" in short else 0.01,
             "not_a_measurement": "Too short." if f"seed{n}" in short else None}
            for n, m in enumerate(means, start=1)]
    average = sum(means) / len(means)
    spread = (sum((m - average) ** 2 for m in means) / (len(means) - 1)) ** 0.5
    entry = {"members": rows, "mean_of_means": sum(means) / len(means), "spread_of_means": spread,
             "spread_is": "the error on one estimate", "observed_standard_error": spread}
    if short:
        entry["qualified_by_their_own_runs"] = list(short)
    (folder / "members.json").write_text(json.dumps(
        {"members": [r["member"] for r in rows], "replicas": True, "analyses": {"rmsd": entry}}),
        encoding="utf-8")


def test_a_run_going_says_how_far_and_how_long_it_has_left(tmp_path):
    from fastmdxplora.gui.workspace import card_of

    card = card_of(_live(_study(tmp_path / "going")))
    # 6,000 steps taken, 2,000 of them equilibration: 4,000 of production's
    # 10,000 at 2 fs, 8 ps of 20.
    assert card["state"] == "running"
    assert card["simulated_ns"] == pytest.approx(0.008)
    assert card["planned_ns"] == pytest.approx(0.02)
    assert card["progress"] == pytest.approx(0.4)
    # 6,000 steps more at production's own 20 steps a second, as the
    # Overview's time left reads it; not at the whole time's 381.555 s a
    # 6,000, setup in it, which said 382.
    assert card["seconds_left"] == 300
    ended = card_of(_live(_study(tmp_path / "ended"), status="completed", current=12000))
    assert ended["simulated_ns"] == pytest.approx(0.02) and "progress" not in ended
    # Stopped once its quantities were determined, short of the ceiling:
    # what it ran, as its Overview says, not the 20 ps it might have.
    early = _live(_study(tmp_path / "early"), status="completed", current=8000)
    (early / "simulation" / "simulation_parameters.json").write_text(
        json.dumps({"duration_ns_actual": 0.012}), encoding="utf-8")
    early = card_of(early)
    assert early["simulated_ns"] == pytest.approx(0.012)
    assert early["planned_ns"] == pytest.approx(0.02)


def test_a_trajectory_analysed_is_analysed_and_every_mean_is_sent(tmp_path):
    from fastmdxplora.gui.workspace import card_of

    analysed = tmp_path / "helical"
    (analysis := analysed / "analysis").mkdir(parents=True)
    (analysis / "analysis_manifest.json").write_text("{}", encoding="utf-8")
    card = card_of(analysed)
    assert (card["state"], card["kind"]) == ("analysed", "a trajectory analysed")
    (tmp_path / "empty").mkdir()
    assert card_of(tmp_path / "empty")["state"] == "not started"
    many = _study(tmp_path / "many",
                  means={n: (1.0, 0.01) for n in ("rmsd", "rg", "sasa", "hbonds")})
    assert len(card_of(many)["means"]) == 4


def test_a_study_of_replicas_gives_its_mean_of_means(tmp_path):
    from fastmdxplora.gui.workspace import card_of

    study = _replicas(tmp_path / "chignolin")
    _members(study, [0.10, 0.12, 0.14])
    card = card_of(study)
    rmsd = card["means"][0]
    assert card["runs"] == 3
    assert rmsd["mean"] == pytest.approx(0.12) and rmsd["spread"] == pytest.approx(0.02)
    # The error of the mean of three: their spread over the root of three.
    assert rmsd["error"] == pytest.approx(0.02 / math.sqrt(3)) and rmsd["withheld"] is None
    _members(study, [0.10, 0.12, 0.14], short=("seed2",))
    rmsd = card_of(study)["means"][0]
    assert rmsd["error"] is None and rmsd["spread"] == pytest.approx(0.02)
    assert rmsd["withheld"] == ("The spread of 3 replicas' means; 1 of them could not determine "
                                "their own mean, so the spread is not yet an error.")


def test_the_replicas_mean_and_error_are_the_runs_tables(tmp_path, monkeypatch):
    """One derivation: the replicas' mean and its error as the runs' table
    gives them (`runs_compared.mean_of_replicas`), not a second formula."""
    from fastmdxplora.gui import runs_compared
    from fastmdxplora.gui.workspace import card_of

    study = _replicas(tmp_path / "chignolin")
    _members(study, [0.10, 0.12, 0.14])
    entry = json.loads((study / "members.json").read_text(encoding="utf-8"))["analyses"]["rmsd"]
    together = runs_compared.mean_of_replicas(entry)
    rmsd = card_of(study)["means"][0]
    assert (rmsd["mean"], rmsd["error"], rmsd["spread"]) == (
        together["mean"], together["error"], together["spread"])
    monkeypatch.setattr(runs_compared, "mean_of_replicas", lambda found: {**together, "error": 0.5})
    assert card_of(study)["means"][0]["error"] == 0.5


def test_the_molecule_is_its_backbone_in_its_own_axes(tmp_path):
    from fastmdxplora.gui import backbone_picture
    from fastmdxplora.gui.backbone_picture import LIST_POINTS, backbone_points

    study = _study(tmp_path / "helix")
    assert backbone_points(study) is None
    _pdb(study / "setup" / "prepared.pdb", _helix(20, "A", 1))
    drawn = backbone_points(study)
    (piece,) = drawn["pieces"]
    assert drawn["from"] == "structure" and drawn["residues"] == 20 and len(piece) == 20
    # The N terminus on the left, coloured from 0 to 1 along the chain,
    # the height filled, the width its own; depth from 0 to 1.
    assert piece[0][0] < piece[-1][0]
    assert [p[3] for p in piece] == sorted(p[3] for p in piece) and piece[-1][3] == 1
    assert min(p[1] for p in piece) == 0 and max(p[1] for p in piece) == 1
    assert max(p[0] for p in piece) == pytest.approx(drawn["aspect"], abs=1e-3)
    assert drawn["aspect"] >= 1
    assert {min(p[2] for p in piece), max(p[2] for p in piece)} == {0, 1}
    # Kept until the file changes.
    assert backbone_points(study) == drawn
    assert any(key[0].endswith("prepared.pdb") for key in backbone_picture._POINTS_KEPT)
    # The latest frame a run wrote is drawn before the structure it began from.
    _pdb(study / "simulation" / "live_frame.pdb", _helix(500, "A", 1))
    latest = backbone_points(study)
    assert latest["from"] == "latest frame" and latest["residues"] == 500
    assert sum(len(p) for p in latest["pieces"]) <= LIST_POINTS


LOOK = """() => ({
    rows: [...document.querySelectorAll('#studies-table tbody tr.studies-row')]
        .map((r) => r.querySelector('.studies-name-text').textContent),
    paths: [...document.querySelectorAll('#studies-table tbody tr.studies-row')]
        .map((r) => r.dataset.path),
    children: document.querySelectorAll('#studies-table tr.studies-child').length,
    groups: [...document.querySelectorAll('#studies-table tr.studies-group b')]
        .map((b) => b.textContent),
    tabs: [...document.querySelectorAll('#studies-views .studies-view-tab')]
        .map((t) => [...t.childNodes].map((n) => n.textContent.trim()).filter(Boolean).join(' ')),
    selected: (document.querySelector('#studies-views [aria-selected="true"]') || {}).textContent,
    grid: !document.getElementById('studies-grid').hidden,
    focused: (document.activeElement.closest('tr') || {}).dataset?.path || '',
    pop: (document.querySelector('.studies-means-pop') || {}).textContent || '',
    wide: document.documentElement.scrollWidth > innerWidth,
})"""


def test_the_list_first_with_views_molecules_means_and_keys(tmp_path):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    helix = _study(tmp_path / "helix", system="1L2Y", started="2026-09-03T10:00:00+00:00",
                   means={"rmsd": (0.05, 0.01), "rg": (0.7, 0.01), "sasa": (12.0, 0.1)})
    _pdb(helix / "setup" / "prepared.pdb", _helix(20, "A", 1))
    _live(_study(tmp_path / "going", system="1UAO", started="2026-09-02T10:00:00+00:00"))
    replicas = _replicas(tmp_path / "chignolin")
    _members(replicas, [0.10, 0.12, 0.14], short=("seed2",))
    for n in (1, 2, 3):
        run = replicas / "runs" / f"seed{n}"
        (run / "resolved_config.yml").write_text("systems: [{system: 1UAO}]\n", encoding="utf-8")
        (run / "manifest.json").write_text(json.dumps(
            {"phases": [{"name": "simulation", "status": "ok"}]}), encoding="utf-8")
    (tmp_path / "begun").mkdir()
    (tmp_path / "begun" / "resolved_config.yml").write_text("systems: [{system: 1ABC}]\n")
    session = start_dashboard_session(output=str(helix), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            first = page.evaluate(LOOK)
            row = '#studies-table tr.studies-row[data-path$="helix"]'
            molecule = page.eval_on_selector(row, "(r) => [r.querySelector('.studies-mol-frame')"
                                                  ".dataset.picture, r.querySelectorAll("
                                                  "'svg.study-molecule line').length]")
            going = page.text_content(
                '#studies-table tr.studies-row[data-path$="going"] .studies-state')
            bar = page.get_attribute('#studies-table tr[data-path$="going"] .studies-progress',
                                     "aria-valuenow")
            means = page.eval_on_selector(
                '#studies-table tr.studies-row[data-path$="chignolin"]',
                "(r) => [...r.querySelectorAll('.studies-mean-item')].map((m) => "
                "[m.dataset.determined, m.querySelector('.studies-mean-v').textContent])")
            # The replicas unfolded: read when asked.
            page.click('#studies-table tr[data-path$="chignolin"] .studies-fold')
            page.wait_for_selector("#studies-table tr.studies-child")
            unfolded = page.evaluate(LOOK)
            # A view: Running alone.
            page.click('#studies-views .studies-view-tab:has-text("Running")')
            running = page.evaluate(LOOK)
            page.click('#studies-views .studies-view-tab:has-text("All")')
            # The person's own view: a search kept under a name, and kept.
            page.fill("#studies-search", "1L2Y")
            page.click("#studies-views .studies-view-save")
            page.fill("#studies-views .studies-view-name", "Trp-cage")
            page.keyboard.press("Enter")
            page.wait_for_selector('#studies-views .studies-view-tab.is-saved:has-text("Trp-cage")')
            page.reload(wait_until="domcontentloaded")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            page.click('#studies-views .studies-view-tab:has-text("Trp-cage")')
            saved = page.evaluate(LOOK)
            page.click('#studies-views .studies-view-tab:has-text("All")')
            # The keys: j moves down the rows, x chooses, m every mean.
            page.click("#studies-views .studies-view-tab[aria-selected='true']")
            page.keyboard.press("j")
            page.keyboard.press("j")
            moved = page.evaluate(LOOK)
            page.keyboard.press("x")
            chosen = page.eval_on_selector_all("#studies-table .studies-pick input:checked",
                                               "(all) => all.map((b) => b.dataset.path)")
            page.keyboard.press("m")
            page.wait_for_selector(".studies-means-pop")
            popped = page.evaluate(LOOK)
            page.keyboard.press("Escape")
            page.keyboard.press("/")
            searching = page.evaluate("() => document.activeElement.id")
            page.keyboard.press("Escape")
            # The cards: a gallery of the same molecules.
            page.click('[data-studies-view="cards"]')
            page.wait_for_selector(".study-card")
            gallery = page.eval_on_selector(
                '.study-card[data-path$="helix"] .study-thumb',
                "(t) => [t.dataset.picture, Math.round(t.getBoundingClientRect().height)]")
            page.click('[data-studies-view="table"]')
            page.set_viewport_size({"width": 390, "height": 844})
            page.wait_for_timeout(300)
            phone = page.evaluate(LOOK)
            phone_row = page.eval_on_selector(row, "(r) => getComputedStyle(r).display")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert not first["grid"]
    # Most recently active first; the folder where nothing ran at the foot.
    assert first["rows"][-1] == "begun" and first["groups"] == ["Not started"]
    assert first["tabs"][:4] == ["All 4", "Running 1", "Needs more sampling 1", "Not started 1"]
    assert molecule == ["molecule", 19]
    assert going.strip().startswith("running, 40%") and "8 ps of 20 ps · 5 min left" in going
    assert bar == "40"
    assert means == [["false", "0.120 nm, spread 0.020"]]
    assert unfolded["children"] == 3
    assert running["rows"] == ["going"] and running["selected"].startswith("Running")
    assert saved["rows"] == ["helix"] and saved["selected"].startswith("Trp-cage")
    assert moved["focused"] == moved["paths"][1]
    assert len(chosen) == 1 and chosen[0] == moved["focused"]
    assert "every mean recorded" in popped["pop"]
    assert searching == "studies-search"
    assert gallery == ["molecule", 130]
    assert phone_row == "grid" and not phone["wide"]


def test_why_a_mean_is_not_determined_is_reached_from_the_keyboard(tmp_path):
    """A row of one or two means, one not determined, said why only on
    hover over its ring: no button to reach it from the keyboard, as +N is
    where there are more means than a row shows."""
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    replicas = _replicas(tmp_path / "chignolin")
    _members(replicas, [0.10, 0.12, 0.14], short=("seed2",))
    _study(tmp_path / "helix", means={"rmsd": (0.05, 0.01)})
    session = start_dashboard_session(output=str(tmp_path / "helix"), host="127.0.0.1", port=0)
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            page.wait_for_selector("#studies-table tbody tr.studies-row")
            # The page settled: drawn again, a list closes what is open on it.
            page.wait_for_selector("body.state-ready")
            page.evaluate(f"() => window.FastMDXStudies.load({str(tmp_path)!r})")
            row = '#studies-table tr.studies-row[data-path$="/chignolin"]'
            offered = page.locator(f"{row} .studies-means-more").count()
            determined = page.locator(
                '#studies-table tr.studies-row[data-path$="/helix"] .studies-means-more').count()
            page.focus(f"{row} .studies-means-more")
            page.keyboard.press("Enter")
            page.wait_for_selector(".studies-means-pop")
            said = page.inner_text(".studies-means-pop")
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert offered == 1 and determined == 0
    assert "could not determine their own mean" in said
