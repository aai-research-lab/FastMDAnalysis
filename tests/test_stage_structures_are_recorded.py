"""A stage is not a promise that coordinates were retained for it."""
from __future__ import annotations

import builtins
import contextlib
import json
import math
import shutil
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest


DATA = Path(__file__).parent / "data" / "stage_states"


def _api():
    from fastmdxplora.gui import stage_structures
    return stage_structures


@pytest.fixture
def study(tmp_path, monkeypatch):
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    root = tmp_path / "study"
    (root / "setup").mkdir(parents=True)
    (root / "simulation").mkdir()
    (root / "manifest.json").write_text(json.dumps({"phases": [
        {"name": "setup", "status": "completed"}, {"name": "simulation", "status": "completed"}]}))
    for name in ["input.pdb", "prepared.pdb", "solvated.pdb"]:
        shutil.copy(DATA / "topology.pdb", root / "setup" / name)
    shutil.copy(DATA / "topology.pdb", root / "simulation" / "topology.pdb")
    for name in ["state_minimized.xml", "state_final.xml"]:
        shutil.copy(DATA / "state.xml", root / "simulation" / name)
    (root / "simulation" / "live_status.json").write_text(json.dumps({"status": "completed",
        "stage": "production", "stage_states": {"minimization": "completed", "nvt": "completed",
        "npt": "completed", "production": "completed"}}))
    return root


def test_openmm_state_coordinates_are_converted_from_nm_in_atom_order():
    pdb, reason = _api().state_coordinates((DATA / "topology.pdb").read_text(), (DATA / "state.xml").read_text())
    assert reason == "" and pdb
    atoms = [line for line in pdb.splitlines() if line.startswith("ATOM")]
    assert [line[12:16].strip() for line in atoms] == ["N", "CA", "C", "O"]
    assert [[float(line[a:b]) for a, b in [(30,38),(38,46),(46,54)]] for line in atoms] == [
        [3,4,5],[4.4,4,5],[5.7,4,5],[6.7,4.5,5]]


def test_the_box_is_the_states_triclinic_box_in_angstroms():
    pdb, _ = _api().state_coordinates((DATA / "topology.pdb").read_text(), (DATA / "state.xml").read_text())
    line = pdb.splitlines()[0]
    assert line.startswith("CRYST1")
    edges = [float(line[a:b]) for a,b in [(6,15),(15,24),(24,33)]]
    assert edges == pytest.approx([30, math.sqrt(4**2+28**2), math.sqrt(2**2+3**2+32**2)], abs=.001)
    angles = [float(line[a:b]) for a,b in [(33,40),(40,47),(47,54)]]
    assert all(0 < x < 180 for x in angles)
    assert angles[2] == pytest.approx(math.degrees(math.acos(4 / edges[1])), abs=.01)


@pytest.mark.parametrize("old,new", [
    ('x=".3"','x="nan"'), ('x=".3"','x="inf"'), ('x=".3"','x="10000000"'),
    ('version="1"','version="2"'), ('<Positions>','<Missing>'),
    ('<A x="3" y="0" z="0"/>','<A x="3" y="1" z="0"/>'),
    ('<A x="3" y="0" z="0"/>','<A x="0" y="0" z="0"/>')])
def test_invalid_state_coordinates_are_withheld(old, new):
    result, reason = _api().state_coordinates((DATA / "topology.pdb").read_text(),
        (DATA / "state.xml").read_text().replace(old,new))
    assert result is None and reason


def test_atom_count_mismatch_is_withheld():
    topology = (DATA / "topology.pdb").read_text()
    result, reason = _api().state_coordinates("\n".join(line for line in topology.splitlines()
        if not line.startswith("ATOM      4")), (DATA / "state.xml").read_text())
    assert result is None and "atom count" in reason


def test_xml_entities_are_not_expanded():
    result, reason = _api().state_coordinates((DATA / "topology.pdb").read_text(),
        '<!DOCTYPE State [<!ENTITY x "large">]><State>&x;</State>')
    assert result is None and "entities" in reason


def test_completion_and_retained_coordinates_are_separate(study):
    entries = {e["id"]:e for e in _api().timeline_payload(study)["stages"]}
    assert list(entries) == ["input","prepared","solvated","minimized","equilibrated","production"]
    assert entries["minimized"]["state"] == "completed" and entries["minimized"]["available"]
    assert entries["equilibrated"]["state"] == "completed" and not entries["equilibrated"]["available"]
    assert "not retained" in entries["equilibrated"]["reason"]


def test_a_missing_state_never_uses_the_current_topology_as_its_snapshot(study):
    (study / "simulation" / "state_minimized.xml").unlink()
    assert not next(e for e in _api().timeline_payload(study)["stages"] if e["id"] == "minimized")["available"]
    data, entry = _api().stage_snapshot(study, "minimized")
    assert data is None and entry["reason"]


def test_an_analysis_only_study_has_no_unperformed_simulation_stages(study):
    shutil.rmtree(study / "setup")
    for name in ["state_minimized.xml","state_final.xml","live_status.json"]:
        (study / "simulation" / name).unlink()
    (study / "manifest.json").write_text('{"phases":[{"name":"analysis","status":"completed"}]}')
    assert _api().timeline_payload(study)["stages"] == []


def test_interruption_and_skipped_equilibration_stay_named(study):
    (study / "simulation" / "live_status.json").write_text(json.dumps({"stage":"npt",
        "status":"interrupted", "stage_states":{"nvt":"skipped","npt":"current"}}))
    entry = next(e for e in _api().timeline_payload(study)["stages"] if e["id"] == "equilibrated")
    assert entry["state"] == "interrupted" and not entry["available"]


def test_the_current_stage_is_named_without_a_stage_states_record(study):
    (study / "simulation" / "live_status.json").write_text(json.dumps({
        "stage":"production","status":"running"}))
    entry=next(e for e in _api().timeline_payload(study)["stages"] if e["id"]=="production")
    assert entry["state"]=="current"


def test_a_changed_source_version_is_not_silently_displayed(study):
    entry = next(e for e in _api().timeline_payload(study)["stages"] if e["id"] == "minimized")
    target = study / "simulation" / "state_minimized.xml"
    target.write_text(target.read_text()+"\n")
    data, refused = _api().stage_snapshot(study,"minimized",entry["version"])
    assert data is None and "changed" in refused["reason"]


def test_inspection_does_not_import_the_simulation_engine(study, monkeypatch):
    api = _api()
    original = builtins.__import__
    def checked(name, *args, **kwargs):
        assert not name.startswith("openmm")
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,"__import__",checked)
    data, _ = api.stage_snapshot(study,"minimized",with_solvent=True)
    assert b"   3.000   4.000   5.000" in data


@pytest.mark.parametrize("target", ["setup/input.pdb","simulation/state_minimized.xml","simulation/topology.pdb"])
def test_links_outside_the_study_are_never_opened(study, tmp_path, target):
    outside=tmp_path/"private"; outside.write_text("private marker")
    path=study/target; path.unlink()
    try: path.symlink_to(outside)
    except OSError: pytest.skip("symlinks unavailable")
    payload=_api().timeline_payload(study)
    data, refused=_api().stage_snapshot(study,"minimized")
    assert not payload["available"] and not payload["stages"]
    assert data is None and "private marker" not in str(refused)


def test_oversized_sources_are_named_without_reading_them(study, monkeypatch):
    monkeypatch.setattr(_api(),"MAX_SOURCE_BYTES",10)
    entry=next(e for e in _api().timeline_payload(study)["stages"] if e["id"]=="minimized")
    assert not entry["available"] and "size limit" in entry["reason"]


@contextlib.contextmanager
def _serving(root, public=False):
    from fastmdxplora.gui.server import make_handler
    server=ThreadingHTTPServer(("127.0.0.1",0),make_handler(root,allow_control=not public))
    threading.Thread(target=server.serve_forever,daemon=True).start()
    try: yield f"http://127.0.0.1:{server.server_port}"
    finally: server.shutdown(); server.server_close()


def _get(url):
    try:
        with urllib.request.urlopen(url) as r: return r.status,r.read()
    except urllib.error.HTTPError as e: return e.code,e.read()


@pytest.mark.parametrize("public",[False,True])
def test_read_only_routes_return_recorded_coordinates(study,public):
    with _serving(study,public) as url:
        code,body=_get(url+"/api/stage-timeline")
        assert code==200
        payload=json.loads(body)
        assert payload["available"] and len(payload["stages"])==6
        code,pdb=_get(url+"/structure/stage.pdb?stage=minimized&solvent=1")
        assert code==200 and b"   3.000   4.000   5.000" in pdb


def test_study_mismatch_is_refused_before_serving_coordinates(study):
    with _serving(study) as url:
        code,body=_get(url+"/structure/stage.pdb?stage=minimized&study=0000000000000000")
        assert code==409 and "changed" in json.loads(body)["reason"]


def test_stage_names_cannot_name_arbitrary_files(study):
    with _serving(study) as url:
        code,body=_get(url+"/structure/stage.pdb?stage=../../private")
        assert code==422 and "no retained coordinates" in json.loads(body)["reason"]


def test_a_study_of_runs_uses_the_viewed_runs_state_and_names_it(study,monkeypatch):
    from fastmdxplora.gui.workspace import study_key
    child=study/"runs"/"replica";child.mkdir(parents=True)
    shutil.copytree(study/"simulation",child/"simulation")
    shutil.copy(study/"manifest.json",child/"manifest.json")
    state=child/"simulation"/"state_minimized.xml"
    state.write_text(state.read_text().replace('x=".3"','x=".9"'))
    monkeypatch.setattr("fastmdxplora.gui.runs_together.run_shown",lambda root:child)
    with _serving(study) as url:
        code,body=_get(url+"/api/stage-timeline")
        payload=json.loads(body)
        assert code==200 and payload["study"]==study_key(study)
        assert payload["run"]==study_key(child) and payload["run_name"]=="runs/replica"
        code,pdb=_get(url+"/structure/stage.pdb?stage=minimized&solvent=1")
        assert code==200 and b"   9.000   4.000   5.000" in pdb
        code,body=_get(url+"/structure/stage.pdb?stage=minimized&run="+study_key(study))
        assert code==409 and "changed" in json.loads(body)["reason"]


@pytest.mark.parametrize("run_id",[".private","../agent/conversations/private"])
@pytest.mark.parametrize("public",[False,True])
def test_a_nested_run_cannot_remove_its_private_parent_boundary(study,run_id,public):
    child=study/"runs"/run_id;child.mkdir(parents=True)
    shutil.copytree(study/"simulation",child/"simulation")
    (child/"simulation"/"production.dcd").write_bytes(b"recorded trajectory")
    with _serving(study,public) as url:
        code,body=_get(url+"/api/stage-timeline")
        assert code==200 and json.loads(body)["available"]
        (study/"batch_manifest.json").write_text(json.dumps({"planned":[{"run_id":run_id}]}))
        code,body=_get(url+"/api/stage-timeline")
        assert code==403 and "available" not in json.loads(body)
        code,body=_get(url+"/structure/stage.pdb?stage=minimized&solvent=1")
        assert code==403 and b"ATOM" not in body


@pytest.mark.parametrize("step,expected",[(10,"interrupted"),(100,"completed")])
def test_a_silent_productions_state_is_reconciled_from_recorded_steps(study,step,expected):
    (study/"simulation"/"live_status.json").write_text(json.dumps({
        "stage":"production","status":"running","current_step":step,"total_planned_steps":100,
        "last_update_timestamp":"2020-01-01T00:00:00Z","stage_states":{"production":"current"}}))
    payload=_api().timeline_payload(study)
    assert payload["available"]
    entry=next(e for e in payload["stages"] if e["id"]=="production")
    assert entry["state"]==expected and entry["available"]


def test_public_telemetry_does_not_attempt_the_optional_private_process_record(study,monkeypatch):
    from datetime import datetime,timezone
    from fastmdxplora import record_access
    from fastmdxplora.gui.server import _private
    from fastmdxplora.gui.telemetry import status_as_it_stands
    from fastmdxplora.orchestrator import RUN_PROCESS_FILE

    (study/"simulation"/"live_status.json").write_text(json.dumps({"status":"running",
        "stage":"production","last_update_timestamp":datetime.now(timezone.utc).isoformat()}))
    original=record_access.read_text
    def checked(path,**kwargs):
        assert Path(path).name != RUN_PROCESS_FILE
        return original(path,**kwargs)
    monkeypatch.setattr(record_access,"read_text",checked)
    with record_access.confined_records(study,lambda path:not _private(path.relative_to(study).parts)) as scope:
        assert status_as_it_stands(study,include_campaign=False)["status"]=="running"
    assert not scope.denied
    assert _api().timeline_payload(study)["available"]


@pytest.mark.parametrize("route", ["/api/stage-timeline", "/structure/stage.pdb?stage=minimized&solvent=1"])
@pytest.mark.parametrize("run_id", ["replica", "../../outside", ".private"])
@pytest.mark.parametrize("adaptive", [False, True])
def test_run_selection_admits_sources_before_reading_them(study, monkeypatch, route, run_id, adaptive):
    from datetime import datetime, timezone
    from fastmdxplora import record_access
    from fastmdxplora.orchestrator import RUN_PROCESS_FILE

    child = study / "runs" / run_id
    child.mkdir(parents=True)
    shutil.copytree(study / "simulation", child / "simulation")
    (child / "simulation" / "production.dcd").write_bytes(b"recorded trajectory")
    (child / "simulation" / "live_status.json").write_text(json.dumps({
        "stage": "production", "status": "running",
        "last_update_timestamp": datetime.now(timezone.utc).isoformat()}))
    (child / RUN_PROCESS_FILE).write_text("{}")
    (study / "batch_manifest.json").write_text(json.dumps({"planned": [{"run_id": run_id}]}))
    if adaptive:
        (study / "stopping.json").write_text(json.dumps({"outcome": "running", "runs": [run_id],
            "rounds": [{"decision": "extend", "production_ns": 1, "more_ns": 1}]}))
        (study / RUN_PROCESS_FILE).write_text("{}")
    reads = []
    original = record_access.read_text

    def checked(path, **kwargs):
        resolved = Path(path).resolve()
        assert resolved.is_relative_to(study), "Selection read outside the study."
        assert not any(part.startswith(".") for part in resolved.relative_to(study).parts)
        reads.append(resolved)
        return original(path, **kwargs)

    # Server construction may adopt a local process. This public request
    # must discover its run entirely through admitted public records.
    with _serving(study, public=True) as url:
        monkeypatch.setattr(record_access, "read_text", checked)
        code, body = _get(url + route)
    assert code == (200 if run_id == "replica" else 403)
    if run_id == "replica":
        assert child / "simulation" / "live_status.json" in reads
    else:
        assert b"ATOM" not in body


def test_a_selected_runs_joined_record_cannot_read_an_outside_file(study, monkeypatch, tmp_path):
    child = study / "runs" / "replica"
    child.mkdir(parents=True)
    shutil.copytree(study / "simulation", child / "simulation")
    (child / "joined").mkdir()
    outside = tmp_path / "private.json"
    outside.write_text('{"secret":"PRIVATE_STAGE_SELECTION"}')
    (child / "joined" / "joined.json").symlink_to(outside)
    (study / "batch_manifest.json").write_text('{"planned":[{"run_id":"replica"}]}')
    original = Path.open
    def checked(path, *args, **kwargs):
        assert path.resolve() != outside, "Selection opened an outside joined record."
        return original(path, *args, **kwargs)
    with _serving(study, public=True) as url:
        monkeypatch.setattr(Path, "open", checked)
        code, body = _get(url + "/api/stage-timeline")
    assert code == 403 and b"PRIVATE_STAGE_SELECTION" not in body


def test_the_latest_continuation_status_replaces_the_first_pieces_status(study):
    (study/"simulation"/"live_status.json").write_text(json.dumps({
        "stage":"production","status":"stopped","last_update_timestamp":"2020-01-01T00:00:00Z",
        "stage_states":{"minimization":"completed","nvt":"completed","npt":"completed","production":"stopped"}}))
    piece=study/"segment-001"/"simulation";piece.mkdir(parents=True)
    (piece/"live_status.json").write_text(json.dumps({"stage":"production","status":"failed",
        "last_update_timestamp":"2020-01-02T00:00:00Z","stage_states":{"production":"failed"}}))
    payload=_api().timeline_payload(study)
    assert payload["available"]
    assert next(e for e in payload["stages"] if e["id"]=="production")["state"]=="failed"


def test_an_unavailable_process_query_does_not_call_a_live_process_dead(monkeypatch):
    import os
    import subprocess
    from datetime import datetime,timezone
    from fastmdxplora.gui import exploration,telemetry
    from fastmdxplora.gui.exploration import _process_alive
    from fastmdxplora.orchestrator import this_machine

    monkeypatch.setattr(subprocess,"run",lambda *args,**kwargs:subprocess.CompletedProcess(
        args[0],1,stdout="",stderr="fatal library error, lookup self"))
    assert _process_alive(os.getpid())
    # A failed query alone is insufficient, but a known newer start time
    # still proves PID reuse even when this runtime cannot expose /proc.
    now=datetime.now(timezone.utc)
    monkeypatch.setattr(exploration,"_process_alive",lambda pid:True)
    monkeypatch.setattr(exploration,"_identify_run",lambda *args:None)
    monkeypatch.setattr(telemetry,"_process_started_at",lambda pid:now.timestamp()+20)
    record={**this_machine(),"pid":os.getpid()+100000,"started_at":now.isoformat(),"argv":[]}
    assert telemetry._what_its_process_says(Path("/an-admitted-study"),record)=="gone"
