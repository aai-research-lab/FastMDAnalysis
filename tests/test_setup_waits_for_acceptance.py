"""Setup-only launch and explicit acceptance keep the same prepared study."""
from __future__ import annotations

import json
import hashlib
from pathlib import Path
import shutil

import pytest
import yaml

from fastmdxplora.gui.run_from_config import prepare_run


DATA = Path(__file__).parent / "data" / "preparation_changes"
PLAN = "simulation_plan.yml"
RECORD = ".preparation-review.json"


def _gate():
    from fastmdxplora.gui import preparation_gate
    return preparation_gate


@pytest.fixture
def prepared(tmp_path):
    root = tmp_path / "study"
    root.mkdir()
    plan = yaml.safe_dump({"systems": [{"system": "1UBQ"}], "include_phase": ["setup", "simulation"],
                           "simulation": {"temperature_K": 310, "duration_ns": 7}})
    (root / PLAN).write_bytes(plan.encode("utf-8"))
    (root / RECORD).write_text(json.dumps({"version": 1, "status": "waiting",
        "plan_sha256": hashlib.sha256(plan.encode()).hexdigest()}))
    shutil.copytree(DATA / "setup", root / "setup")
    (root / "setup" / "system.xml").write_text("<System/>")
    (root / "setup" / "state.xml").write_text("<State/>")
    (root / "setup" / "setup_parameters.json").write_text('{"phase":"setup"}')
    (root / "resolved_config.yml").write_text("include_phase: [setup]\n")
    (root / "manifest.json").write_text('{"phases":[{"name":"setup","status":"ok"}]}')
    return root


def test_first_launch_runs_only_audited_setup_and_keeps_the_full_plan(tmp_path):
    root = tmp_path / "new"
    result = prepare_run({"system": "1UBQ", "include_phase": ["setup", "simulation"],
        "review_after_setup": True, "simulation": {"duration_ns": 7, "temperature_K": 310}}, root)
    assert result["ok"] and result["review_after_setup"]
    config = yaml.safe_load(Path(result["config_path"]).read_text())
    assert config["include_phase"] == ["setup"]
    assert config["setup"]["preparation_audit"] is True
    assert config["systems"][0]["setup"]["preparation_audit"] is True
    plan = yaml.safe_load((root / PLAN).read_text())
    assert hashlib.sha256((root / PLAN).read_bytes()).hexdigest() == json.loads((root / RECORD).read_text())["plan_sha256"]
    assert plan["include_phase"] == ["setup", "simulation"]
    assert plan["simulation"]["temperature_K"] == 310
    assert plan["simulation"]["duration_ns"] == 7
    assert "--no-defaults" in result["command"]
    assert not _gate().review_state(root)["can_accept"]


@pytest.mark.parametrize("change", [
    {"systems": [{"system": "1UBQ"}, {"system": "1L2Y"}]},
    {"include_phase": ["simulation"]}, {"sweep": {"simulation.temperature_K": [300, 310]}},
    {"simulation": {"resume_from": "old"}}, {"simulation": {"setup_from": "old"}},
    {"budget_hours": 2},
])
def test_unsupported_workflows_do_not_silently_bypass_the_stop(change):
    config = {"systems": [{"system": "1UBQ"}], "include_phase": ["setup", "simulation"], **change}
    assert _gate().scope_reason(config)
    assert _gate().scope_reason(config, hosted=True)


def test_the_completed_audit_is_ready_but_never_starts_a_process_on_read(prepared):
    first = _gate().review_state(prepared)
    assert first["requested"] and first["can_accept"] and first["status"] == "waiting"
    assert _gate().review_state(prepared) == first
    assert json.loads((prepared / RECORD).read_text())["status"] == "waiting"


def test_acceptance_uses_the_fixed_plan_and_existing_simulate_command(prepared):
    from fastmdxplora.cli.main import _build_parser
    calls = []
    def spawn(command, plan):
        calls.append((command, plan))
        args = _build_parser().parse_args(command[3:])
        assert args.command == "simulate"
        assert args.config == str(plan) and plan.parent == prepared
        kept = yaml.safe_load(plan.read_text())
        assert kept["simulation"]["temperature_K"] == 310
        assert kept["simulation"]["duration_ns"] == 7
        assert Path(kept["simulation"]["setup_from"]).parent.parent == prepared
        assert args.output_dir == str(prepared)
        return {"launched": True, "pid": 123}
    revision = _gate().review_state(prepared)["revision"]
    result = _gate().accept(prepared, revision, spawn)
    assert result["ok"] and result["review"] == "accepted" and len(calls) == 1
    assert not _gate().review_state(prepared)["can_accept"]
    assert not _gate().accept(prepared, revision, spawn)["ok"]
    assert len(calls) == 1


@pytest.mark.parametrize("source", ["simulation_plan.yml", "setup/topology.pdb", "setup/system.xml",
    "setup/state.xml", "setup/preparation_audit.json", "setup/audit/source-0002.pdb", "resolved_config.yml"])
def test_changed_evidence_or_plan_refuses_a_stale_acceptance(prepared, source):
    before = _gate().review_state(prepared)
    target = prepared / source
    target.write_text(target.read_text() + "\n")
    calls = []
    result = _gate().accept(prepared, before["revision"], lambda *args: calls.append(args))
    assert not result["ok"] and not calls
    assert json.loads((prepared / RECORD).read_text())["status"] == "waiting"
    assert not _gate().review_state(prepared)["can_accept"]


@pytest.mark.parametrize("status", ["partial", "incomplete", "stopped", "recording"])
def test_incomplete_evidence_never_offers_accept(prepared, status):
    path = prepared / "setup" / "preparation_audit.json"
    record = json.loads(path.read_text()); record["status"] = status
    path.write_text(json.dumps(record))
    assert not _gate().review_state(prepared)["can_accept"]


def test_existing_simulation_output_refuses_a_second_start(prepared):
    folder = prepared / "simulation"; folder.mkdir()
    (folder / "exploration.log").write_text("Previous failed attempt")
    assert not _gate().review_state(prepared)["can_accept"]


def test_a_failed_launch_does_not_record_acceptance(prepared):
    revision = _gate().review_state(prepared)["revision"]
    result = _gate().accept(prepared, revision, lambda *args: {"ok": False, "error": "No GPU available."})
    assert not result["ok"] and result["error"] == "No GPU available."
    assert _gate().review_state(prepared)["can_accept"]
    assert json.loads((prepared / RECORD).read_text())["status"] == "waiting"


def test_an_outside_snapshot_is_not_read_or_accepted(prepared, tmp_path, monkeypatch):
    private = tmp_path / "private.pdb"; private.write_text("PRIVATE_MARKER")
    target = prepared / "setup" / "audit" / "source-0002.pdb"
    target.unlink(); target.symlink_to(private)
    original = Path.open
    def checked(path, *args, **kwargs):
        assert path.resolve() != private
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "open", checked)
    result = _gate().review_state(prepared)
    assert not result["can_accept"] and "PRIVATE_MARKER" not in str(result)


def _runtime(prepared, tmp_path, monkeypatch):
    from fastmdxplora.gui import exploration
    monkeypatch.setattr(exploration, "exploration_environment_error", lambda config: None)
    runtime = exploration.DashboardRuntime(workspace_root=tmp_path / "workspace",
        exploration_root=tmp_path, active_root=prepared)
    return runtime


def test_runtime_refuses_a_changed_study_before_using_the_acceptance(prepared, tmp_path, monkeypatch):
    runtime = _runtime(prepared, tmp_path, monkeypatch)
    calls = []
    monkeypatch.setattr(runtime, "_spawn", lambda *args, **kwargs: calls.append(args))
    answer = runtime.accept_preparation({"study": "wrong", "revision": _gate().review_state(prepared)["revision"]})
    assert not answer["ok"] and "changed" in answer["error"] and not calls


def test_duplicate_acceptance_requests_share_one_launch_lock(prepared, tmp_path, monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    from fastmdxplora.gui.workspace import study_key

    runtime = _runtime(prepared, tmp_path, monkeypatch)
    calls = []
    def spawn(*args, **kwargs):
        calls.append(args)
        return {"launched": True, "pid": 123}
    monkeypatch.setattr(runtime, "_spawn", spawn)
    request = {"study": study_key(prepared), "revision": _gate().review_state(prepared)["revision"]}
    with ThreadPoolExecutor(max_workers=2) as pool:
        answers = list(pool.map(lambda _: runtime.accept_preparation(request), range(2)))
    assert sum(answer["ok"] for answer in answers) == 1
    assert len(calls) == 1


def test_an_active_workflow_holds_acceptance_without_changing_it(prepared, tmp_path, monkeypatch):
    from types import SimpleNamespace
    from fastmdxplora.gui.workspace import study_key
    runtime = _runtime(prepared, tmp_path, monkeypatch)
    runtime.process = SimpleNamespace(poll=lambda: None)
    answer = runtime.accept_preparation({"study": study_key(prepared), "revision": _gate().review_state(prepared)["revision"]})
    assert not answer["ok"] and "already running" in answer["error"]
    assert json.loads((prepared / RECORD).read_text())["status"] == "waiting"


def test_hosted_runtime_cannot_accept_or_launch_a_local_review(prepared, tmp_path, monkeypatch):
    from types import SimpleNamespace
    runtime = _runtime(prepared, tmp_path, monkeypatch)
    runtime.hosting = SimpleNamespace()
    assert not runtime.accept_preparation({})["ok"]


def test_public_views_inspect_only_public_preparation_evidence(prepared, monkeypatch):
    from tests.test_stage_structures_are_recorded import _get, _serving
    from fastmdxplora import record_access
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _ask
    original = record_access.open_text
    def checked(path, **kwargs):
        assert not any(part.startswith(".") for part in Path(path).relative_to(prepared).parts)
        return original(path, **kwargs)
    monkeypatch.setattr("fastmdxplora.gui.preparation_review.open_text", checked)
    monkeypatch.setattr("fastmdxplora.gui.preparation_gate.review_state",
                        lambda root: pytest.fail("A public view read the private acceptance ledger."))
    with _serving(prepared, public=True) as url:
        # The named imported helper reads via the same admission scope.
        code, body = _get(url + "/api/preparation-review")
        answer = json.loads(body)
        assert code == 200 and answer["available"] and not answer["control"]
        assert not answer["gate"]["requested"]
        code, body = _get(url + "/api/preparation-change?event=event-0002&version=" + answer["version"])
        assert code == 200 and len(json.loads(body)["added"]) == 1
        code, body = _ask(url + "/api/accept-preparation", body={"study": answer["study"],
            "revision": "A public view cannot read the private acceptance revision."})
        assert code == 403
    assert json.loads((prepared / RECORD).read_text())["status"] == "waiting"


def test_running_setup_does_not_bind_partial_outputs(prepared):
    answer = _gate().review_state(prepared, running=True)
    assert answer["requested"] and not answer["can_accept"]
    assert "setup_revision" not in json.loads((prepared / RECORD).read_text())


def test_an_unresolved_spawn_exception_does_not_allow_another_launch(prepared):
    def failed(*args):
        raise RuntimeError("Spawn failed before starting a child.")
    revision = _gate().review_state(prepared)["revision"]
    with pytest.raises(RuntimeError, match="Spawn failed"):
        _gate().accept(prepared, revision, failed)
    assert not _gate().review_state(prepared)["can_accept"]
    assert json.loads((prepared / RECORD).read_text())["status"] == "starting"


def test_a_missing_simulation_backend_does_not_accept_or_launch(prepared, tmp_path, monkeypatch):
    from fastmdxplora.gui.workspace import study_key
    runtime = _runtime(prepared, tmp_path, monkeypatch)
    monkeypatch.setattr("fastmdxplora.gui.exploration.exploration_environment_error", lambda config: "OpenMM is unavailable.")
    calls = []
    monkeypatch.setattr(runtime, "_spawn", lambda *args, **kwargs: calls.append(args))
    result = runtime.accept_preparation({"study": study_key(prepared), "revision": _gate().review_state(prepared)["revision"]})
    assert not result["ok"] and "OpenMM is unavailable" in result["error"] and not calls
    assert json.loads((prepared / RECORD).read_text())["status"] == "waiting"


def test_a_remote_plan_refuses_the_requested_pause_before_reading_a_machine(prepared, tmp_path, monkeypatch):
    from fastmdxplora.gui.remote_routes import RemoteDesk
    runtime = _runtime(prepared, tmp_path, monkeypatch)
    desk = RemoteDesk(runtime)
    monkeypatch.setattr("fastmdxplora.remote.inputs.home_folders", lambda: pytest.fail("Remote planning read a machine."))
    result = desk.plan(None, "remote", state={"review_after_setup": True})
    assert not result["ok"] and "locally" in str(result)


def test_a_hosted_save_cannot_drop_the_requested_pause(prepared, tmp_path, monkeypatch):
    from types import SimpleNamespace
    runtime = _runtime(prepared, tmp_path, monkeypatch)
    runtime.hosting = SimpleNamespace()
    result = runtime.save_config_for_elsewhere({"review_after_setup": True})
    assert not result["ok"] and "locally" in result["error"]


def test_actual_post_popen_response_failure_cannot_launch_a_second_child(prepared, tmp_path, monkeypatch):
    from types import SimpleNamespace
    from fastmdxplora.gui.workspace import study_key

    runtime = _runtime(prepared, tmp_path, monkeypatch)
    calls = []
    monkeypatch.setattr(runtime, "_others_running", lambda: None)
    monkeypatch.setattr(runtime, "_spawn", lambda command, root, url, **kwargs:
                        runtime._spawn_now(command, root, url))
    def popen(*args, **kwargs):
        calls.append(args[0])
        return SimpleNamespace(pid=987654, poll=lambda: 0)
    monkeypatch.setattr("fastmdxplora.gui.exploration.subprocess.Popen", popen)
    request = {"study": study_key(prepared), "revision": _gate().review_state(prepared)["revision"]}
    with monkeypatch.context() as patch:
        patch.setattr(runtime, "snapshot", lambda: (_ for _ in ()).throw(OSError("Post-Popen snapshot failed")))
        with pytest.raises(OSError, match="Post-Popen"):
            runtime.accept_preparation(request)
    assert len(calls) == 1
    assert json.loads((prepared / RECORD).read_text())["status"] == "starting"
    assert not runtime.accept_preparation(request)["ok"]
    assert len(calls) == 1


@pytest.mark.parametrize("source", ["setup/system.xml", "setup/state.xml", PLAN])
def test_changes_during_backend_preflight_refuse_before_spawn(prepared, tmp_path, monkeypatch, source):
    from fastmdxplora.gui.workspace import study_key

    runtime = _runtime(prepared, tmp_path, monkeypatch)
    calls = []
    def probe(config):
        path = prepared / source
        path.write_text(path.read_text() + "\n")
        return None
    monkeypatch.setattr("fastmdxplora.gui.exploration.exploration_environment_error", probe)
    monkeypatch.setattr(runtime, "_spawn", lambda *args, **kwargs: calls.append(args))
    result = runtime.accept_preparation({"study": study_key(prepared), "revision": _gate().review_state(prepared)["revision"]})
    assert not result["ok"] and "changed" in result["error"] and not calls
    assert json.loads((prepared / RECORD).read_text())["status"] == "waiting"


def test_child_inputs_keep_the_checked_bytes_after_original_sources_change(prepared):
    from fastmdxplora.simulation.pipeline import _where_the_system_was_prepared

    original = {name: (prepared / "setup" / name).read_bytes()
                for name in ("system.xml", "state.xml", "topology.pdb")}
    original_plan = (prepared / PLAN).read_bytes()
    def spawn(command, plan):
        assert Path(command[command.index("--config") + 1]) == plan
        for name in original:
            (prepared / "setup" / name).write_text("changed while the child was starting")
        (prepared / PLAN).write_text("changed: original plan")
        kept = yaml.safe_load(plan.read_text())
        location, named = _where_the_system_was_prepared(None, kept["systems"][0]["simulation"]["setup_from"])
        assert named and location.parent.parent == prepared
        assert all((location / name).read_bytes() == data for name, data in original.items())
        retained = yaml.safe_load(original_plan)
        for key, value in retained["simulation"].items():
            assert kept["simulation"][key] == value
        assert all((location / name).stat().st_mode & 0o222 == 0 for name in original)
        return {"launched": True, "pid": 123}
    answer = _gate().accept(prepared, _gate().review_state(prepared)["revision"], spawn)
    assert answer["ok"]


def test_changes_during_input_copy_refuse_without_spawning(prepared, monkeypatch):
    read = _gate()._read
    calls = []
    def changing(path, limit):
        data = read(path, limit)
        if path.name == "system.xml":
            (prepared / "setup" / "state.xml").write_text("changed during copy")
        return data
    monkeypatch.setattr(_gate(), "_read", changing)
    result = _gate().accept(prepared, _gate().review_state(prepared)["revision"], lambda *args: calls.append(args))
    assert not result["ok"] and not calls
