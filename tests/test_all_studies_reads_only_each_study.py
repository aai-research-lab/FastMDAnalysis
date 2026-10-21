"""All studies reads each study only, as the study's own pages do.

A study's card (`workspace.card_of`, behind `/api/studies`) read its
records through any link: a downloaded study holding a link to a file
outside it had that file read and shown as its means, tags, state and
molecule, and `/api/study-thumbnail` sent a figure it linked to. Each card
is now read held to its study (`server._confined`, `record_access`): a
record linked out of it, or a file of it that is not one of its results
(a hidden one), is not opened, and the card names it and says what it
leaves unknown, never a value taken from its absence. Checked by an audit
of the opens, on loopback and on a public bind.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from tests.test_the_study_pages_read_only_the_study import _OPENED, _WATCHED, _link, _watch

#: A Cα trace of ten residues, 3.8 Å apart.
PDB = "".join(
    f"ATOM  {i:5d}  CA  ALA A{i:4d}    {3.8 * i:8.3f}{1.5 * (i % 2):8.3f}{0.5 * (i % 3):8.3f}"
    "  1.00  0.00           C\n" for i in range(1, 11)) + "END\n"
PIXEL = bytes.fromhex(
    "89504e470d0a1a0a0000000d4948445200000001000000010806000000"
    "1f15c4890000000d49444154789c6360000002000154a24f5f0000000049454e44ae426082")
#: What a study's records say once they are linked out: none of it may
#: reach a card.
LEAKED = ("0.4321", "2.2222", "leaked tag", "leaked note", "Leaked windows")
OUT = dict(rmsd=0.4321, rg=2.2222, tag="leaked tag", note="leaked note", pmf="Leaked windows")


def _write(root: Path, relative: str, value: object) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(value, bytes):
        path.write_bytes(value)
    else:
        path.write_text(value if isinstance(value, str) else json.dumps(value), encoding="utf-8")
    return path


def _mean(name: str, value: float) -> dict:
    return {"analysis": name, "findings": {"mean": {
        "mean": value, "standard_error": 0.002, "unit": "nm", "effective_samples": 40}}}


def _a_study(root: Path, *, rmsd=0.1234, rg=1.0567, tag="wild type", note="kept with care",
             pmf="Too few windows") -> Path:
    """A study run and analysed, with every record a card reads."""
    _write(root, "manifest.json", {"phases": [
        {"name": "setup", "status": "ok", "started_at": "2026-09-01T10:00:00+00:00"},
        {"name": "simulation", "status": "ok"}, {"name": "analysis", "status": "ok"}]})
    _write(root, "resolved_config.yml", yaml.safe_dump({
        "systems": [{"system": "1UBQ"}], "simulation": {"duration_ns": 1.0}}))
    _write(root, "analysis/analysis_manifest.json", {"plan": ["rmsd", "rg"]})
    _write(root, "analysis/rmsd/options.json", _mean("rmsd", rmsd))
    _write(root, "analysis/rg/options.json", _mean("rg", rg))
    _write(root, "analysis/rmsd/rmsd.png", PIXEL)
    _write(root, "analysis/rmsd/rmsd.dat", "".join(f"{0.1 + 0.001 * i:.4f}\n" for i in range(20)))
    _write(root, "simulation/live_status.json", {
        "status": "completed", "stage": "production", "timestep_fs": 2.0,
        "nvt_steps_planned": 0, "npt_steps_planned": 0,
        "production_steps_planned": 500000, "current_step": 500000,
        "total_planned_steps": 500000})
    _write(root, "setup/prepared.pdb", PDB)
    _write(root, "study_tags.json", {"tags": [tag], "note": note,
                                     "updated": "2026-09-02T10:00:00+00:00"})
    _write(root, "pmf.json", {"refused": f"{pmf}. More are needed."})
    return root


def _linked_out(root: Path, outside: Path) -> None:
    """Each regular file of the study moved outside it, a link left in its
    place: what the link leads to is the study's own record."""
    for path in sorted(root.rglob("*")):
        if path.is_file() and not path.is_symlink():
            target = outside / path.relative_to(root)
            target.parent.mkdir(parents=True, exist_ok=True)
            path.rename(target)
            _link(path, target)


@pytest.fixture
def workspace(tmp_path) -> tuple[Path, Path]:
    space, outside = tmp_path / "ws", tmp_path / "outside"
    _a_study(space / "own")
    _linked_out(_a_study(space / "linked", **OUT), outside / "linked")
    # A record linked to a hidden file of the study: not one of its results.
    hidden = _a_study(space / "hidden", rg=3.3333)
    (hidden / "analysis" / "rg" / "options.json").rename(hidden / ".rg.json")
    _link(hidden / "analysis" / "rg" / "options.json", Path("../../.rg.json"))
    # Replicas whose comparison is linked out, and one of whose runs is.
    reps = space / "reps"
    planned = [{"run_id": f"seed{n}", "system": "1UAO", "options": {},
                "sweep_values": {"simulation.random_seed": n}} for n in (1, 2)]
    _write(reps, "batch_manifest.json", {
        "n_runs": 2, "planned": planned, "sweep": {"simulation.random_seed": [1, 2]},
        "systems": [{"id": "s", "system": "1UAO"}],
        "runs": [{"run_id": p["run_id"], "status": "ok"} for p in planned]})
    _write(reps, "resolved_config.yml", yaml.safe_dump({"systems": [{"system": "1UAO"}]}))
    _a_study(reps / "runs" / "seed1")
    _a_study(outside / "seed2", **OUT)
    _link(reps / "runs" / "seed2", outside / "seed2")
    _write(outside, "members.json", {"replicas": True, "analyses": {"rmsd": {
        "mean_of_means": 0.4321, "spread_of_means": 0.01,
        "members": [{"member": "seed1", "mean": 0.12}, {"member": "seed2", "mean": 0.13}]}}})
    _link(reps / "members.json", outside / "members.json")
    return space, outside


def _opened_outside(outside: Path, read) -> tuple[object, list[str]]:
    _watch(outside)
    try:
        answer = read()
        return answer, sorted(set(_OPENED))
    finally:
        _WATCHED.clear()


def test_a_card_reads_nothing_through_a_link_and_says_what_it_did_not_read(workspace):
    from fastmdxplora.gui.workspace import studies_in

    space, outside = workspace
    found, opened = _opened_outside(outside, lambda: studies_in(space))
    cards = {card["name"]: card for card in found["studies"]}
    runs, opened_runs = _opened_outside(outside, lambda: studies_in(space / "reps" / "runs"))
    assert opened == [] and opened_runs == []
    # The study read whole: all of it, nothing said not read.
    own = cards["own"]
    assert own["state"] == "completed" and own["kind"] == "one run"
    assert [m["mean"] for m in own["means"]] == [0.1234, 1.0567]
    assert own["tags"] == ["wild type"] and own["simulated_ns"] == 1.0
    assert own["molecule"]["residues"] == 10 and own["free_energy"]
    assert "unavailable" not in own and "unread" not in own and "not_read" not in own
    # Every record linked out: none of it on the card, each named.
    linked = cards["linked"]
    text = json.dumps(linked)
    for value in LEAKED:
        assert value not in text, value
    assert linked["state"] == "not known" and linked["kind"] == ""
    assert linked["means"] == [] and linked["tags"] == [] and linked["molecule"] is None
    assert "simulated_ns" not in linked and linked["last_active"] == ""
    assert set(linked["unread"]) == {"state", "length", "means", "molecule"}
    assert {"manifest.json", "analysis/rmsd/options.json", "analysis/rg/options.json",
            "study_tags.json", "simulation/live_status.json",
            "setup/prepared.pdb"} <= set(linked["unavailable"])
    assert "analysis/rg/options.json" in linked["not_read"]
    assert "a link to a file outside the study" in linked["not_read"]
    # A record linked to a hidden file of the study: not read, the rest is.
    hidden = cards["hidden"]
    assert [m["analysis"] for m in hidden["means"]] == ["rmsd"]
    assert hidden["unread"] == ["means"] and hidden["unavailable"] == ["analysis/rg/options.json"]
    assert "not one of its results" in hidden["not_read"]
    assert hidden["state"] == "completed" and "3.3333" not in json.dumps(hidden)
    # Replicas: their comparison linked out is not their mean, and a run
    # linked out is not counted in how long they ran.
    reps = cards["reps"]
    assert reps["means"] == [] and "means" in reps["unread"]
    assert "simulated_ns" not in reps and "length" in reps["unread"]
    assert reps["state"] == "not known"
    assert "0.4321" not in json.dumps(reps)
    named = {card["name"]: card for card in runs["studies"]}
    assert named["seed1"]["means"] and "unavailable" not in named["seed1"]
    assert named["seed2"]["state"] == "not known" and named["seed2"]["means"] == []
    assert named["seed2"]["not_read"] == ("runs/seed2 was not read: "
                                          "it is a link to a folder outside the study.")


def test_a_start_outside_the_study_is_not_drawn(tmp_path):
    from fastmdxplora.gui.workspace import card_of

    start = _write(tmp_path / "structures", "start.pdb", PDB)
    study = tmp_path / "ws" / "waiting"
    _write(study, "resolved_config.yml", yaml.safe_dump({"systems": [{"system": str(start)}]}))
    _write(study, "manifest.json", {"phases": []})
    card, opened = _opened_outside(tmp_path / "structures", lambda: card_of(study))
    assert opened == []
    assert card["molecule"] is None and card["molecule_outside"] == "start.pdb"
    assert card["state"] == "not started"
    # Where it is in the study, it is drawn.
    _write(study, "start.pdb", PDB)
    (study / "resolved_config.yml").write_text(
        yaml.safe_dump({"systems": [{"system": "start.pdb"}]}), encoding="utf-8")
    assert card_of(study)["molecule"]["from"] == "start"


@pytest.mark.parametrize("public", [False, True], ids=["loopback", "public"])
def test_the_routes_read_nothing_through_a_link(workspace, public):
    import urllib.error
    import urllib.request
    from urllib.parse import urlencode

    from fastmdxplora.gui.exploration import DashboardRuntime
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _serving

    space, outside = workspace
    # A study whose only structure is linked out: no figure, no backbone.
    drawn = _a_study(space / "drawn")
    (drawn / "analysis" / "rmsd" / "rmsd.png").unlink()
    (outside / "drawn").mkdir()
    (drawn / "setup" / "prepared.pdb").rename(outside / "drawn" / "prepared.pdb")
    _link(drawn / "setup" / "prepared.pdb", outside / "drawn" / "prepared.pdb")
    runtime = DashboardRuntime(workspace_root=space / "own", exploration_root=space,
                               active_root=space / "own")
    asked = {"studies": f"/api/studies?{urlencode({'path': str(space)})}",
             "runs": f"/api/studies?{urlencode({'path': str(space / 'reps' / 'runs')})}",
             "figure": f"/api/study-thumbnail?{urlencode({'path': str(space / 'linked')})}",
             "series": f"/api/study-thumbnail?{urlencode({'path': str(space / 'linked'), 'series': 1})}",
             "backbone": f"/api/study-thumbnail?{urlencode({'path': str(drawn)})}",
             "own": f"/api/study-thumbnail?{urlencode({'path': str(space / 'own')})}"}
    said: dict[str, tuple[int, bytes]] = {}

    def ask_all() -> None:
        with _serving(runtime, allow_control=not public) as url:
            for key, route in asked.items():
                try:
                    with urllib.request.urlopen(url + route, timeout=60) as response:
                        said[key] = (response.status, response.read())
                except urllib.error.HTTPError as error:
                    said[key] = (error.code, error.read())

    _, opened = _opened_outside(outside, ask_all)
    assert opened == []
    for key, (_, body) in said.items():
        if key == "own":
            continue
        for value in LEAKED:
            assert value.encode() not in body, (key, value)
        assert PIXEL not in body and b"<svg" not in body, key
    if public:
        # Beyond loopback none of these is answered.
        assert {status for status, _ in said.values()} == {403}
        return
    assert said["own"] == (200, PIXEL)
    cards = {card["name"]: card for card in json.loads(said["studies"][1])["studies"]}
    assert cards["linked"]["state"] == "not known" and cards["own"]["means"]
    assert said["figure"][0] == 404
    assert "analysis/rmsd/rmsd.png" in json.loads(said["figure"][1])["unavailable"]
    assert said["backbone"][0] == 404
    assert json.loads(said["backbone"][1])["unavailable"] == ["setup/prepared.pdb"]
    series = json.loads(said["series"][1])
    assert series["ok"] is False and "unavailable" in series


def test_the_page_says_what_a_card_did_not_read(workspace):
    pytest.importorskip("playwright.sync_api")
    from playwright.sync_api import sync_playwright

    from fastmdxplora.gui.server import start_dashboard_session

    space, _ = workspace
    session = start_dashboard_session(output=str(space / "own"), host="127.0.0.1", port=0)
    look = """(name) => {
        const row = document.querySelector(`#studies-table tr.studies-row[data-path$="/${name}"]`);
        const frame = row.querySelector('.studies-mol-frame');
        return {state: row.querySelector('.studies-state-line').textContent.trim(),
                length: row.querySelector('.studies-length').textContent.trim(),
                means: row.querySelector('.studies-means-cell').textContent.trim(),
                said: (row.querySelector('.studies-not-read') || {}).textContent || '',
                why: (row.querySelector('.studies-not-read') || {}).title || '',
                molecule: frame ? frame.title : ''};
    }"""
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch(args=["--enable-unsafe-swiftshader"])
            page = browser.new_page(viewport={"width": 1440, "height": 900})
            page.set_default_timeout(60000)
            errors: list[str] = []
            page.on("pageerror", lambda error: errors.append(str(error)))
            page.goto(session.url + "#studies", wait_until="domcontentloaded")
            page.evaluate(f"() => window.FastMDXStudies.load({str(space)!r})")
            page.wait_for_selector('#studies-table tr.studies-row[data-path$="/linked"]')
            linked = page.evaluate(look, "linked")
            hidden = page.evaluate(look, "hidden")
            own = page.evaluate(look, "own")
            page.click('[data-studies-view="cards"]')
            card = page.inner_text('.study-card[data-path$="/linked"]')
            browser.close()
    finally:
        session.server.shutdown()
    assert errors == []
    assert linked["state"] == "not known"
    assert linked["length"] == "not read" and linked["means"] == "not read"
    assert linked["said"].startswith("Not read: ") and "manifest.json" in linked["said"]
    assert "a link to a file outside the study" in linked["why"]
    assert linked["molecule"].startswith("Not drawn: its structure was not read")
    assert "Rg" not in hidden["means"] and "not read" in hidden["means"]
    assert own["said"] == "" and "not read" not in own["means"]
    assert "Not read: " in card
