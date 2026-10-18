"""Problems admits actual record reads without inspecting unrelated trees."""

from __future__ import annotations

import json
import logging
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier

import pytest


def _write(root, name, value):
    path = root / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def _study(root):
    from fastmdxplora.gui.telemetry import TelemetryWriter

    _write(root, "simulation/live_status.json", {"status": "completed", "stage": "production"})
    _write(root, "manifest.json", {"phases": [{"name": "simulation", "status": "ok"}]})
    TelemetryWriter(root / "simulation").event("Recorded warning.", level="warning")
    return root


def _link(path, target):
    try:
        path.symlink_to(target, target_is_directory=target.is_dir())
    except OSError:
        pytest.skip("symbolic links are unavailable")


def _public(root):
    from fastmdxplora.gui.server import _private

    return lambda path: path.is_relative_to(root) and not _private(path.relative_to(root).parts)


@pytest.mark.parametrize("moved_checkpoint", [False, True])
def test_unrelated_links_and_files_do_not_require_a_tree_scan(tmp_path, monkeypatch, moved_checkpoint):
    from fastmdxplora.gui.problems_view import problems_payload

    root = _study(tmp_path / "study")
    outside = _write(tmp_path, "large.dcd", "not a record")
    _link(root / "simulation/production.dcd", outside)
    if moved_checkpoint:
        _write(root, "simulation/live_status.json", {
            "status": "completed", "stage": "production",
            "current_checkpoint_path": str(tmp_path / "old-machine/study/simulation/checkpoint.chk")})
        (root / "simulation/checkpoint.chk").write_bytes(b"checkpoint")
    for index in range(100):
        _write(root, f"scratch/{index}/unused.json", {"unused": True})

    def no_recursive_walk(*args, **kwargs):
        pytest.fail("Problems recursively walked the study")

    monkeypatch.setattr(Path, "rglob", no_recursive_walk)
    payload = problems_payload(root)
    assert payload["available"]
    assert payload["items"][0]["reason"] == "Recorded warning."


def test_selected_segment_keeps_own_records_without_guessing_about_parent_checkpoints(tmp_path):
    from fastmdxplora import record_access
    from fastmdxplora.gui.problems_view import problems_payload
    from fastmdxplora.remedies import remedy_for

    root = _study(tmp_path / "segment-001")
    refusal = {"code": "simulation.run.stopped", "message": "Stopped."}
    _write(root, "simulation/live_status.json", {"status": "stopped", "stage": "production"})
    payload = problems_payload(root)
    assert payload["available"]
    assert "owning study" in payload["items"][0]["fix"]
    assert "no checkpoint" not in payload["items"][0]["fix"]
    with record_access.confined_records(root, _public(root)) as scope:
        remedy = remedy_for(refusal, study=root, price=False)
        assert "owning study" in remedy.fix
        assert not scope.denied


def test_endpoint_admits_only_candidate_evidence(tmp_path, monkeypatch):
    from fastmdxplora.gui import server
    from fastmdxplora.gui.exploration import DashboardRuntime
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _ask, _serving

    root = _study(tmp_path / "study")
    runtime = DashboardRuntime(workspace_root=root, exploration_root=tmp_path, active_root=root)

    def no_artifact_listing(*args, **kwargs):
        pytest.fail("Problems enumerated all artifacts")

    monkeypatch.setattr(server, "_artifact_records", no_artifact_listing)
    with _serving(runtime, allow_control=False) as url:
        status, body = _ask(url + "/api/problems")
    assert status == 200
    payload = json.loads(body)
    assert payload["items"][0]["file"] == "simulation/live_events.log"


def test_unexpected_reader_failure_is_logged_without_public_exception_details(
        tmp_path, monkeypatch, caplog):
    from fastmdxplora.gui import problems_view

    root = _study(tmp_path)

    def broken(*args, **kwargs):
        raise RuntimeError("PRIVATE_READER_CONTEXT")

    monkeypatch.setattr(problems_view, "overview_of", broken)
    with caplog.at_level(logging.DEBUG, logger=problems_view.__name__):
        payload = problems_view.problems_payload(root)
    assert payload["items"][0]["reason"] == "Recorded warning."
    assert any(record.exc_info and "analysis" in record.message for record in caplog.records)
    assert "PRIVATE_READER_CONTEXT" in caplog.text
    assert "PRIVATE_READER_CONTEXT" not in json.dumps(payload)


def test_selected_run_does_not_consult_parent_campaign(tmp_path, monkeypatch):
    from fastmdxplora.gui import telemetry
    from fastmdxplora.gui.problems_view import problems_payload

    root = _study(tmp_path / "campaign/runs/one")

    def parent_is_not_selected(*args, **kwargs):
        pytest.fail("Selected-run Problems consulted a parent campaign")

    monkeypatch.setattr(telemetry, "between_rounds", parent_is_not_selected)
    payload = problems_payload(root)
    assert payload["available"]
    assert payload["items"][0]["reason"] == "Recorded warning."


def test_unused_axis_and_data_discovery_preserve_recorded_verdicts(tmp_path, monkeypatch):
    from fastmdxplora.gui import analysis_overview
    from fastmdxplora.gui.problems_view import problems_payload
    from fastmdxplora.statistics import MEAN_RULES

    _study(tmp_path)
    _write(tmp_path, "analysis/analysis_manifest.json", {"results": {"rmsd": {"status": "ok"}}})
    _write(tmp_path, "analysis/rmsd/options.json", {"findings": {"mean": {
        "mean": 1, "standard_error": 0.1, "n_frames": 100, "discard": 10,
        "effective_samples": 4, "rules": MEAN_RULES,
        "not_a_measurement": "Constant series."}}})

    def unused(*args, **kwargs):
        pytest.fail("Problems computed an unused clock or listed data files")

    monkeypatch.setattr(analysis_overview, "analysed_axis", unused)
    monkeypatch.setattr(analysis_overview, "_data_files", unused)
    payload = problems_payload(tmp_path)
    assert any(item["reason"] == "Constant series." for item in payload["items"])


@pytest.mark.parametrize("checkpoint", [False, True])
def test_checkpoint_wording_does_not_need_price_reads(tmp_path, monkeypatch, checkpoint):
    from fastmdxplora import remedies
    from fastmdxplora.gui.problems_view import problems_payload

    _study(tmp_path)
    _write(tmp_path, "simulation/live_status.json", {"status": "stopped", "stage": "production"})
    if checkpoint:
        (tmp_path / "simulation/checkpoint.chk").write_bytes(b"checkpoint")

    def unused(*args, **kwargs):
        pytest.fail("Problems priced an undisplayed remedy")

    monkeypatch.setattr(remedies, "_what_remains", unused)
    monkeypatch.setattr(remedies, "_speed", unused)
    payload = problems_payload(tmp_path)
    fix = payload["items"][0]["fix"]
    assert ("no checkpoint" in fix) is not checkpoint
    assert ("Carry it on" in fix) is checkpoint


@pytest.mark.parametrize("source", ["simulation/live_metrics.csv",
                                    "analysis/reweighted/reweighted_averages.json",
                                    "segment-001/simulation/live_status.json"])
def test_nested_readers_refuse_a_link_installed_after_the_first_record(
        tmp_path, monkeypatch, source):
    from fastmdxplora import record_access
    from fastmdxplora.gui.problems_view import problems_payload

    root = _study(tmp_path / "study")
    secret = _write(root, ".private/record.json", {"latest_error": "PRIVATE_CONTEXT"})
    original = record_access.open_text
    installed = False

    def install_after_manifest(path, **kwargs):
        nonlocal installed
        stream = original(path, **kwargs)
        if Path(path) == root / "manifest.json" and not installed:
            link = root / source
            link.parent.mkdir(parents=True, exist_ok=True)
            _link(link, secret)
            installed = True
        return stream

    monkeypatch.setattr(record_access, "open_text", install_after_manifest)
    payload = problems_payload(root)
    assert installed
    assert payload["available"] is False and payload["items"] == []
    assert "PRIVATE_CONTEXT" not in json.dumps(payload)
    assert any("confined" in note for note in payload["notes"])


def test_confined_open_uses_the_admitted_destination_if_a_link_is_retargeted(
        tmp_path, monkeypatch):
    from fastmdxplora import record_access

    root = tmp_path / "study"
    safe = _write(root, "public.json", {"value": "public"})
    private = _write(root, ".private/secret.json", {"value": "PRIVATE"})
    link = root / "record.json"
    _link(link, safe)
    checked = record_access.ReadScope.checked

    def retarget(scope, path, **kwargs):
        result = checked(scope, path, **kwargs)
        if path == link:
            link.unlink()
            _link(link, private)
        return result

    monkeypatch.setattr(record_access.ReadScope, "checked", retarget)
    with record_access.confined_records(root, _public(root)):
        assert "PRIVATE" not in record_access.read_text(link)


def test_retargeting_an_admitted_parent_directory_cannot_redirect_the_open(tmp_path, monkeypatch):
    from fastmdxplora import record_access

    if os.open not in os.supports_dir_fd or not hasattr(os, "O_NOFOLLOW"):
        pytest.skip("this regression exercises POSIX directory-descriptor opening")
    root = tmp_path / "study"
    source = _write(root, "records/record.json", {"value": "public"})
    external = _write(tmp_path, "external/record.json", {"value": "PRIVATE"})
    checked = record_access.ReadScope.checked
    changed = False

    def retarget(scope, path, **kwargs):
        nonlocal changed
        result = checked(scope, path, **kwargs)
        if path == source and not changed:
            source.parent.rename(root / "previous")
            _link(root / "records", external.parent)
            changed = True
        return result

    monkeypatch.setattr(record_access.ReadScope, "checked", retarget)
    with record_access.confined_records(root, _public(root)) as scope:
        with pytest.raises(OSError):
            record_access.read_text(source)
        assert scope.denied


def test_read_scopes_are_nested_and_isolated_between_threads(tmp_path):
    from fastmdxplora import record_access

    one = _write(tmp_path / "one", "record.json", {"value": "one"})
    two = _write(tmp_path / "two", "record.json", {"value": "two"})
    barrier = Barrier(2)

    def read(own, other):
        with record_access.confined_records(own.parent, _public(own.parent)) as outer:
            barrier.wait(timeout=10)
            assert json.loads(record_access.read_text(own))["value"] == own.parent.name
            with pytest.raises(PermissionError):
                record_access.read_text(other)
            with record_access.confined_records(other.parent, _public(other.parent)):
                assert record_access.read_text(other) == other.read_text()
            assert record_access.read_text(own) == own.read_text()
            assert outer.denied
        assert record_access.read_text(other) == other.read_text()

    with ThreadPoolExecutor(max_workers=2) as executor:
        first = executor.submit(read, one, two)
        second = executor.submit(read, two, one)
        first.result(timeout=15)
        second.result(timeout=15)


def test_a_record_denied_by_a_nested_reader_is_not_silently_accepted(tmp_path):
    from fastmdxplora import record_access
    from fastmdxplora.gui.telemetry import read_status

    root = tmp_path / "study"
    secret = _write(root, ".private/status.json", {"status": "failed"})
    (root / "simulation").mkdir()
    _link(root / "simulation/live_status.json", secret)
    with record_access.confined_records(root, _public(root)) as scope:
        assert read_status(root) == {}
        assert "simulation/live_status.json" in scope.denied


def test_nonregular_sources_are_refused_without_waiting_for_bytes(tmp_path):
    from fastmdxplora import record_access

    if not hasattr(os, "mkfifo"):
        pytest.skip("named pipes are not supported")
    pipe = tmp_path / "pipe.json"
    os.mkfifo(pipe)
    with record_access.confined_records(tmp_path, _public(tmp_path)) as scope:
        with pytest.raises(PermissionError):
            record_access.read_text(pipe)
        assert scope.denied
