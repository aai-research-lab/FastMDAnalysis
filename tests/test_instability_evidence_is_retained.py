"""Actual Reference-platform forces, and their confined frozen correspondence."""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest
import yaml


@pytest.fixture
def observed(tmp_path):
    openmm = pytest.importorskip("openmm")
    from openmm import app, unit
    from fastmdxplora.simulation.instability import FiniteForceObserver

    root = tmp_path / "study"; (root / "simulation").mkdir(parents=True)
    config = {"systems": [{"system": "1UBQ"}], "output": "failed_run",
        "setup": {"ph": 7.4}, "simulation": {"timestep_fs": 2, "integrator": "verlet",
            "nvt_steps": 10, "npt_steps": 0, "production_steps": 100,
            "state_interval_steps": 10, "trajectory_interval_steps": 10},
        "decisions": {"setup.ph": {"why": "Assay pH.", "source": "person"},
            "simulation.timestep_fs": {"why": "Original choice.", "source": "person"}}}
    (root / "resolved_config.yml").write_text(yaml.safe_dump(config, sort_keys=False))
    (root / "manifest.json").write_text(json.dumps({"phases": [{"phase": "simulation", "status": "failed", "message": "Particle coordinate is NaN"}]}))
    topology = app.Topology()
    chain = topology.addChain("A"); residue = topology.addResidue("GLY", chain, "7", "B")
    topology.addAtom("CA", app.element.carbon, residue)
    topology.addAtom("N", app.element.nitrogen, residue)
    water = topology.addResidue("HOH", topology.addChain("W"), "1")
    topology.addAtom("O", app.element.oxygen, water); topology.addAtom("H1", app.element.hydrogen, water)
    topology.setPeriodicBoxVectors(np.eye(3) * 3 * unit.nanometer)
    system = openmm.System(); force = openmm.CustomExternalForce("k*x*x/2")
    force.addPerParticleParameter("k")
    for i in range(4):
        system.addParticle(12); force.addParticle(i, [100 * (i + 1)])
    system.addForce(force)
    integrator = openmm.VerletIntegrator(.002)
    simulation = app.Simulation(topology, system, integrator, openmm.Platform.getPlatformByName("Reference"))
    simulation.context.setPositions(np.array([[.1, .1, .1], [.2, .1, .1], [.3, .2, .1], [.4, .2, .1]]) * unit.nanometer)
    simulation.context.setVelocities(np.zeros((4, 3)) * unit.nanometer / unit.picosecond)
    observer = FiniteForceObserver(simulation, root / "simulation", 10)
    observer.begin(simulation, "production")
    observer.failed(simulation, RuntimeError("Particle coordinate is NaN"))
    yield {"root": root, "simulation": simulation, "observer": observer, "config": config}
    del simulation, integrator


def _read(root):
    from fastmdxplora.gui.instability_view import instability_payload
    return instability_payload(root)


def _record(root):
    path = root / "simulation/instability.json"
    return path, json.loads(path.read_text())


def test_known_forces_keep_solvent_and_insertion_codes_from_the_same_state(observed):
    data = _read(observed["root"])
    assert data["available"], data
    assert [row["force_kj_mol_nm"] for row in data["ranking"]] == pytest.approx([160, 90, 40, 10])
    assert data["ranking"][0]["identity"] == ["W", "1", "", "HOH", "H1"]
    assert data["ranking"][-1]["identity"] == ["A", "7", "B", "GLY", "CA"]
    assert data["ranking"][0]["position_nm"] == pytest.approx([.4, .2, .1])
    assert data["atom_count"] == 4 and "not necessarily" in data["notice"]
    assert "unknown" in data["failure_step_note"] and data["interval_steps"] == 10
    path, record = _record(observed["root"])
    assert len(record["system_sha256"]) == 64
    assert record["context_parameters"] == {}
    assert record["snapshot_sha256"] == hashlib.sha256((path.parent / record["snapshot"]).read_bytes()).hexdigest()


def test_a_nonfinite_observation_never_replaces_the_last_finite_state(observed):
    from openmm import unit
    sim, observer = observed["simulation"], observed["observer"]
    previous = observer.last
    sim.context.setPositions(np.full((4, 3), np.nan) * unit.nanometer)
    observer.report(sim, sim.context.getState(getPositions=True, getForces=True, getEnergy=True))
    assert observer.last is previous
    observer.failed(sim, RuntimeError("potential energy is NaN"))
    assert _read(observed["root"])["ranking"][0]["force_kj_mol_nm"] == 160


def test_the_reporter_observes_during_integration_without_changing_dynamics(observed):
    from openmm import unit
    sim, observer = observed["simulation"], observed["observer"]
    checkpoint = sim.context.createCheckpoint()
    sim.step(30)
    expected = sim.context.getState(getPositions=True, getVelocities=True)
    sim.context.loadCheckpoint(checkpoint)
    sim.reporters.insert(0, observer); sim.step(30)
    actual = sim.context.getState(getPositions=True, getVelocities=True)
    assert observer.last["step"] == 30 and observer.last["time_ps"] == pytest.approx(.06)
    np.testing.assert_array_equal(expected.getPositions(asNumpy=True).value_in_unit(unit.nanometer), actual.getPositions(asNumpy=True).value_in_unit(unit.nanometer))
    np.testing.assert_array_equal(expected.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond), actual.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond))
    observer.failed(sim, RuntimeError("NaN")); assert _read(observed["root"])["step"] == 30


def test_the_rank_is_bounded_with_deterministic_atom_order_for_ties(observed, monkeypatch):
    from openmm import unit
    from fastmdxplora.simulation import instability
    monkeypatch.setattr(instability, "RANKS", 2)
    sim, observer = observed["simulation"], observed["observer"]
    sim.context.setPositions(np.zeros((4, 3)) * unit.nanometer)
    observer.report(sim, sim.context.getState(getPositions=True, getForces=True, getEnergy=True))
    assert [row["index"] for row in observer.last["ranking"]] == [0, 1]


def test_a_new_run_archives_the_previous_failure_instead_of_reusing_it(observed):
    from fastmdxplora.simulation.instability import FiniteForceObserver
    root = observed["root"]; path, _ = _record(root); old = path.read_bytes()
    FiniteForceObserver(observed["simulation"], root / "simulation", 10)
    assert not path.exists()
    archived = list((root / "simulation").glob("instability-record-*.json"))
    assert len(archived) == 1 and archived[0].read_bytes() == old
    assert not _read(root)["available"]


def test_unverified_enhanced_sampling_is_not_subject_to_extra_force_evaluations(observed, monkeypatch):
    from fastmdxplora.simulation.instability import FiniteForceObserver
    sim, root = observed["simulation"], observed["root"]
    observer = FiniteForceObserver(sim, root / "simulation", 10, supported=False)
    monkeypatch.setattr(sim.context, "getState", lambda **kw: pytest.fail("Unsupported methods must not receive an extra force evaluation."))
    observer.begin(sim, "production")
    observer.failed(sim, RuntimeError("NaN"))
    assert not observer.enabled
    data = _read(root)
    assert not data["available"] and "enhanced sampling" in data["reason"]


def test_a_previous_evidence_link_is_archived_without_reading_its_target(observed, tmp_path):
    from fastmdxplora.simulation.instability import FiniteForceObserver
    root = observed["root"]; path, _ = _record(root)
    outside = tmp_path / "outside.json"; path.rename(outside); original = outside.read_bytes()
    path.symlink_to(outside)
    FiniteForceObserver(observed["simulation"], root / "simulation", 10)
    assert not path.exists() and outside.read_bytes() == original
    archives = list((root / "simulation").glob("instability-unreadable-*.json"))
    assert len(archives) == 1 and archives[0].is_symlink()


def test_stage_transition_drops_the_old_force_model_observation(observed, monkeypatch):
    observer, sim = observed["observer"], observed["simulation"]
    monkeypatch.setattr(sim.context, "getState", lambda **kw: (_ for _ in ()).throw(RuntimeError("unreadable")))
    observer.begin(sim, "npt")
    observer.failed(sim, RuntimeError("Particle coordinate is NaN"))
    data = _read(observed["root"])
    assert not data["available"] and "stage's start" in data["reason"]


def test_failure_evidence_never_masks_or_retries_the_original_failure(observed, monkeypatch):
    observer = observed["observer"]
    def broken(*args): raise OSError("disk full")
    monkeypatch.setattr(observer, "_write", broken)
    observer.failed(observed["simulation"], RuntimeError("NaN"))
    # No integration or retry occurs when saving the evidence fails.
    assert observed["simulation"].currentStep == 0


def test_unrelated_errors_do_not_publish_force_evidence(observed):
    root = observed["root"]; path, _ = _record(root); path.unlink()
    observed["observer"].failed(observed["simulation"], OSError("file not found"))
    assert not path.exists()


@pytest.mark.parametrize("change", ["snapshot", "identity", "coordinate", "magnitude", "index", "units", "step", "interval", "stage", "timestep", "counter"])
def test_tampered_or_unmappable_evidence_is_refused(observed, change):
    root = observed["root"]; path, record = _record(root)
    if change == "snapshot":
        with (path.parent / record["snapshot"]).open("a") as f: f.write("REMARK changed\n")
    elif change == "identity": record["ranking"][0]["identity"][0] = "X"
    elif change == "coordinate": record["ranking"][0]["position_nm"][0] += .1
    elif change == "magnitude": record["ranking"][0]["force_kj_mol_nm"] += 1
    elif change == "index": record["ranking"][0]["index"] = -1
    elif change == "units": record["force_units"] = "N"
    elif change == "step": record["step"] = True
    elif change == "interval": record["interval_steps"] = 0
    elif change == "stage": record["failure_stage"] = "npt"
    elif change == "timestep": record["timestep_fs"] = float("nan")
    elif change == "counter": record["step"] = record["failure_reported_step"] + 1
    path.write_text(json.dumps(record))
    assert not _read(root)["available"]


@pytest.mark.parametrize("target", ["record", "snapshot", "config"])
@pytest.mark.parametrize("where", ["outside", "private"])
def test_outside_and_private_sources_are_not_read(observed, tmp_path, target, monkeypatch, where):
    import os
    from fastmdxplora.gui.instability_view import timestep_proposal
    root = observed["root"]; path, record = _record(root)
    evidence = _read(root)
    source = {"record": path, "snapshot": path.parent / record["snapshot"], "config": root / "resolved_config.yml"}[target]
    folder = tmp_path if where == "outside" else root / ".private"
    folder.mkdir(exist_ok=True)
    outside = folder / ("outside-" + source.name); source.rename(outside); source.symlink_to(outside)
    reads = []
    original = os.fdopen
    identity = outside.stat()
    def watch(descriptor, *a, **kw):
        actual = os.fstat(descriptor)
        if (actual.st_dev, actual.st_ino) == (identity.st_dev, identity.st_ino): reads.append(descriptor)
        return original(descriptor, *a, **kw)
    monkeypatch.setattr(os, "fdopen", watch)
    result = timestep_proposal(root, evidence["version"]) if target == "config" else _read(root)
    assert not result.get("available", result.get("ok")) and reads == []


def test_legacy_run_does_not_evaluate_forces_or_invent_a_snapshot(tmp_path):
    data = _read(tmp_path)
    assert not data["available"] and "Older runs" in data["reason"]
    assert list(tmp_path.iterdir()) == []


def test_a_moved_study_keeps_self_contained_snapshot_identity(observed, tmp_path):
    old = observed["root"]; new = tmp_path / "moved"; old.rename(new)
    assert _read(new)["available"]


def test_proposal_preserves_physical_times_unrelated_settings_and_original_records(observed):
    from fastmdxplora.gui.instability_view import timestep_proposal
    from fastmdxplora.config.loader import validate_config
    root = observed["root"]; before = {p.relative_to(root): p.read_bytes() for p in root.rglob("*") if p.is_file()}
    evidence = _read(root); proposal = timestep_proposal(root, evidence["version"])
    assert proposal["ok"], proposal
    new = yaml.safe_load(proposal["yaml"]); validate_config(new, require_systems=True)
    assert new["simulation"]["timestep_fs"] == 1 and new["simulation"]["production_steps"] == 200
    assert new["simulation"]["state_interval_steps"] == 20
    assert new["simulation"]["checkpoint_interval_steps"] == 20000
    assert new["setup"] == observed["config"]["setup"] and new["systems"] == observed["config"]["systems"]
    assert new["decisions"]["setup.ph"] == observed["config"]["decisions"]["setup.ph"]
    assert new["decisions"]["simulation.timestep_fs"]["alternatives"][0]["earlier_decision"] == observed["config"]["decisions"]["simulation.timestep_fs"]
    assert new["output"] != observed["config"]["output"]
    assert {p.relative_to(root): p.read_bytes() for p in root.rglob("*") if p.is_file()} == before


@pytest.mark.parametrize("change", ["stale", "minimization", "changed-config", "enhanced", "variable", "invalid", "missing-dt"])
def test_unsupported_or_stale_proposal_is_visibly_refused(observed, change):
    from fastmdxplora.gui.instability_view import timestep_proposal
    root = observed["root"]; path, record = _record(root)
    config = observed["config"]
    if change == "minimization": record["stage"] = record["failure_stage"] = "minimization"
    if change in {"changed-config", "enhanced", "variable", "invalid", "missing-dt"}:
        if change == "changed-config": config["setup"]["ph"] = 6
        elif change == "enhanced": config["simulation"]["plumed"] = {"enabled": True}
        elif change == "variable": config["simulation"]["integrator"] = "variable_verlet"
        elif change == "invalid": config["simulation"]["timestep_fs"] = -1
        else: config["simulation"].pop("timestep_fs")
        raw = yaml.safe_dump(config); (root / "resolved_config.yml").write_text(raw)
        if change != "changed-config": record["config_sha256"] = hashlib.sha256(raw.encode()).hexdigest()
    path.write_text(json.dumps(record))
    data = _read(root)
    result = timestep_proposal(root, "stale" if change == "stale" else data["version"])
    assert not result["ok"] and result["reason"]


@pytest.mark.parametrize("live", [False, True])
def test_the_real_runner_retains_the_observation_and_preserves_the_failure(observed, monkeypatch, live):
    from openmm import XmlSerializer, app
    from fastmdxplora.simulation import runner
    import io

    sim, root = observed["simulation"], observed["root"]
    setup = root / "setup"; setup.mkdir()
    (setup / "system.xml").write_text(XmlSerializer.serialize(sim.system))
    (setup / "state.xml").write_text(XmlSerializer.serialize(sim.context.getState(getPositions=True, getVelocities=True)))
    out = io.StringIO(); app.PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(), out, keepIds=True)
    (setup / "topology.pdb").write_text(out.getvalue())
    error = RuntimeError("Particle coordinate is NaN")
    def fail_at_production(simulation, **kw):
        if kw["label"] == "Production":
            simulation.step(10)
            raise error
        return kw.get("current_step", 0)
    monkeypatch.setattr(runner, "_run_md_stage", fail_at_production)
    monkeypatch.setattr(runner, "_run_md_stage_with_live_metrics", lambda omm, simulation, telemetry, **kw: fail_at_production(simulation, **kw))
    with pytest.raises(RuntimeError) as raised:
        runner.run_simulation(system_xml=setup / "system.xml", state_xml=setup / "state.xml",
            topology_pdb=setup / "topology.pdb", output_dir=root / "simulation",
            minimize=False, nvt_steps=0, npt_steps=0, production_steps=20,
            platform="Reference", integrator="verlet", timestep_fs=2, state_interval_steps=10,
            trajectory_interval_steps=10, live_telemetry=live, telemetry_interval=10, save_selection=None)
    assert raised.value is error
    data = _read(root)
    assert data["available"] and data["stage"] == "production" and data["step"] == 10, data


@pytest.mark.parametrize("public", [False, True])
def test_http_inspection_uses_the_selected_nested_run_and_public_views_cannot_propose(observed, tmp_path, monkeypatch, public):
    import http.client
    from urllib.parse import urlencode
    from fastmdxplora.gui.server import start_dashboard_session
    from fastmdxplora.gui.workspace import study_key
    root = tmp_path / "batch"; child = root / "runs/failure"; child.parent.mkdir(parents=True)
    observed["root"].rename(child)
    (root / "batch_manifest.json").write_text(json.dumps({"planned": [{"run_id": "empty"}, {"run_id": "failure"}]}))
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    session = start_dashboard_session(output=root, host="0.0.0.0" if public else "127.0.0.1", port=0)
    def get(path, query):
        connection = http.client.HTTPConnection("127.0.0.1", session.port)
        connection.request("GET", path + "?" + urlencode(query)); response = connection.getresponse()
        result = response.status, json.loads(response.read()); connection.close(); return result
    monkeypatch.setattr(observed["simulation"].context, "getState", lambda **kw: pytest.fail("Inspection must never query forces."))
    try:
        query = {"study": study_key(root), "run": study_key(child)}
        code, data = get("/api/instability", query)
        assert code == 200 and data["available"] and data["run_name"] == "failure", data
        assert [one["label"] for one in data["runs"]] == ["empty", "failure"]
        code, proposal = get("/api/instability-proposal", {**query, "version": data["version"]})
        assert code == (403 if public else 200)
        assert not data["can_propose"] if public else proposal["ok"]
        code, refused = get("/api/instability", {**query, "study": "other"})
        assert code == 409 and "changed" in refused["reason"]
        code, refused = get("/api/instability", {**query, "run": "unadmitted"})
        assert not refused["available"]
    finally: session.stop()


def test_the_run_list_is_bounded_and_never_uses_an_unsafe_run_id(observed):
    from fastmdxplora.gui.instability_view import evidence_runs
    from fastmdxplora.record_access import confined_records
    root = observed["root"]
    manifest = root / "batch_manifest.json"
    for value in ({"planned": [{"run_id": "../outside"}]}, {"planned": [{"run_id": str(i)} for i in range(129)]}):
        manifest.write_text(json.dumps(value))
        with confined_records(root, lambda path: True): assert evidence_runs(root) == []


def test_retargeted_run_cannot_replace_the_study_boundary_before_inspection(observed, tmp_path):
    from fastmdxplora.gui.instability_view import instability_payload, timestep_proposal
    boundary = tmp_path / "batch"; boundary.mkdir()
    child = boundary / "runs/one"; child.parent.mkdir()
    # Admit a real directory, then replace it by a link before the read.
    child.mkdir(); assert child.resolve().is_relative_to(boundary)
    child.rmdir(); child.symlink_to(observed["root"], target_is_directory=True)
    assert not instability_payload(child, boundary=boundary)["available"]
    assert not timestep_proposal(child, "stale", boundary=boundary)["ok"]
