"""The Problems panel gathers recorded evidence without changing the study."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fastmdxplora.gui.workspace import study_key
from fastmdxplora.remedies import remedy_for
from fastmdxplora.statistics import MEAN_RULES


def problems_payload(root, **kwargs):
    # Import per test so the parent-commit check runs every case, even
    # though that commit has no Problems reader yet.
    from fastmdxplora.gui.problems_view import problems_payload as read

    return read(root, **kwargs)


def _write(root: Path, relative: str, value: object) -> None:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def _status(root: Path, state: str, **extra: object) -> None:
    _write(root, "simulation/live_status.json", {"status": state, "stage": "production",
                                                **extra})


def _mean(**extra: object) -> dict:
    return {"mean": 1.4, "standard_error": 0.1, "n_frames": 200,
            "effective_samples": 100, "rules": MEAN_RULES, **extra}


def test_every_quantity_is_checked_including_after_the_first_nine(tmp_path):
    results = {f"quantity_{index}": {"status": "ok"} for index in range(12)}
    _write(tmp_path, "analysis/analysis_manifest.json", {"results": results})
    for name in results:
        _write(tmp_path, f"analysis/{name}/options.json", {"findings": {
            "mean": _mean(), "secondary": _mean(not_a_measurement="Not equilibrated.")}})
    _status(tmp_path, "completed")
    payload = problems_payload(tmp_path)
    assert payload["available"] is True
    assert payload["study"] == study_key(tmp_path)
    assert len(payload["items"]) == 12
    assert all(item["reason"] == "Not equilibrated." for item in payload["items"])
    assert all("longer" not in item["fix"] for item in payload["items"])


def test_old_rules_are_distinct_from_a_withheld_mean(tmp_path):
    _write(tmp_path, "analysis/rmsd/options.json", {"findings": {
        "mean": _mean(rules=0), "other": _mean(not_a_measurement="Periodic image.")}})
    by_id = {item["id"]: item for item in problems_payload(tmp_path)["items"]}
    assert "earlier rules" in by_id["analysis:rmsd:mean"]["title"]
    assert "Analyze again" in by_id["analysis:rmsd:mean"]["fix"]
    assert by_id["analysis:rmsd:other"]["reason"] == "Periodic image."
    assert "code" not in by_id["analysis:rmsd:other"]


def test_nested_thermodynamic_quantities_are_included(tmp_path):
    _write(tmp_path, "analysis/thermodynamics/options.json", {"findings": {
        "thermodynamics": {"samples": 40, "temperature": {
            "units": "K", "mean": 310, "not_a_measurement": "Temperature is drifting."}}}})
    items = problems_payload(tmp_path)["items"]
    assert len(items) == 1
    assert items[0]["id"] == "analysis:thermodynamics:thermodynamics.temperature"
    assert items[0]["reason"] == "Temperature is drifting."


def test_all_refusals_use_the_registry_without_exposing_commands(tmp_path):
    first = {"code": "config.option.unknown", "message": "Unknown setting.",
             "details": {"context": "simulation", "option": "temprature_K"},
             "permitted": ["temperature_K"]}
    second = {"code": "unclassified", "message": "A separate failure.",
              "details": {"install_command": "SECRET_COMMAND"}}
    _write(tmp_path, "manifest.json", {"phases": [
        {"name": "setup", "status": "error", "refusal": first},
        {"name": "analysis", "status": "error", "refusal": second}],
        "argv": ["PRIVATE_ARGV"], "config": {"private": "PRIVATE_CONFIG"}})
    _status(tmp_path, "failed", latest_error="Unknown setting.")
    payload = problems_payload(tmp_path)
    assert len(payload["items"]) == 2
    assert payload["items"][0]["fix"] == remedy_for(first).fix
    assert payload["items"][1]["reason"] == "A separate failure."
    assert set(payload["items"][0]) == {"id", "title", "reason", "fix", "code"}
    assert "PRIVATE" not in json.dumps(payload)
    assert "SECRET_COMMAND" not in json.dumps(payload)


@pytest.mark.parametrize("state", ["stopped", "interrupted"])
def test_stops_and_interruptions_are_not_called_failures(tmp_path, state):
    _status(tmp_path, state)
    items = problems_payload(tmp_path)["items"]
    assert len(items) == 1
    assert items[0]["id"] == f"health:{state}"
    assert items[0]["title"] == state.capitalize()


def test_health_warning_uses_the_existing_telemetry_interpretation(tmp_path):
    _status(tmp_path, "completed", target_temperature_K=300)
    (tmp_path / "simulation/live_metrics.csv").write_text(
        "stage,temperature\nproduction,500\n", encoding="utf-8")
    items = problems_payload(tmp_path)["items"]
    assert len(items) == 1
    assert items[0]["id"] == "health:warning"
    assert "500.0 K" in items[0]["reason"]


def test_a_manifest_stop_does_not_duplicate_or_promise_a_checkpoint(tmp_path):
    _write(tmp_path, "manifest.json", {"phases": [{"name": "simulation", "status": "error",
        "refusal": {"code": "simulation.run.stopped", "message": "Stopped during NVT."}}]})
    _write(tmp_path, "simulation/live_status.json", {"status": "stopped", "stage": "nvt"})
    items = problems_payload(tmp_path)["items"]
    assert len(items) == 1
    assert "before its production began" in items[0]["fix"]
    assert "from the start" in items[0]["fix"]


def test_failed_analysis_keeps_its_reason_and_only_admitted_evidence_links(tmp_path):
    _write(tmp_path, "analysis/analysis_manifest.json", {"results": {
        "rmsd": {"status": "error", "message": "rmsd: Atoms do not match."}}})
    path = "analysis/analysis_manifest.json"
    no_links = problems_payload(tmp_path)["items"]
    assert len(no_links) == 1 and no_links[0]["reason"] == "Atoms do not match."
    assert "file" not in no_links[0]
    admitted = problems_payload(tmp_path, artifacts=[{"path": path}])["items"]
    assert admitted[0]["file"] == path
    denied = problems_payload(tmp_path, artifacts=[{"path": "unrelated.json"}])["items"]
    assert "file" not in denied[0]


@pytest.mark.parametrize("content", ["{", "[]", "{}"])
def test_malformed_and_empty_records_are_not_healthy(tmp_path, content):
    (tmp_path / "manifest.json").write_text(content, encoding="utf-8")
    payload = problems_payload(tmp_path)
    assert payload["available"] is False
    assert payload["items"] == []
    assert any("manifest.json" in note for note in payload["notes"])
    assert any("unknown" in note for note in payload["notes"])


def test_empty_folder_and_campaign_scope_are_stated(tmp_path):
    payload = problems_payload(tmp_path)
    assert payload["available"] is False and payload["notes"]
    _write(tmp_path, "batch_manifest.json", {"runs": []})
    campaign = problems_payload(tmp_path)
    assert any("campaign run" in note for note in campaign["notes"])
    assert campaign["available"] is False


def test_malformed_findings_cannot_be_reported_as_checked(tmp_path):
    _write(tmp_path, "analysis/rmsd/options.json", {"findings": "unreadable"})
    payload = problems_payload(tmp_path)
    assert payload["available"] is False
    assert any("no readable findings" in note for note in payload["notes"])


def test_successful_current_records_have_no_problems_and_are_unchanged(tmp_path):
    _status(tmp_path, "completed")
    _write(tmp_path, "analysis/analysis_manifest.json", {"results": {"rmsd": {"status": "ok"}}})
    _write(tmp_path, "analysis/rmsd/options.json", {"findings": {"mean": _mean()}})
    before = {str(p): (p.read_bytes(), p.stat().st_mtime_ns) for p in tmp_path.rglob("*")
              if p.is_file()}
    payload = problems_payload(tmp_path)
    assert payload["available"] is True and payload["items"] == []
    assert payload["notes"] == []
    assert before == {str(p): (p.read_bytes(), p.stat().st_mtime_ns) for p in tmp_path.rglob("*")
                      if p.is_file()}


def test_analysis_read_failure_does_not_erase_a_stopped_run(tmp_path, monkeypatch):
    from fastmdxplora.gui import problems_view

    _status(tmp_path, "stopped")

    def unreadable(root):
        raise ValueError("private exception context")

    monkeypatch.setattr(problems_view, "overview_of", unreadable)
    payload = problems_payload(tmp_path)
    assert payload["items"][0]["id"] == "health:stopped"
    assert any("could not be read" in note for note in payload["notes"])
    assert "private exception" not in json.dumps(payload)


def test_analysis_name_cannot_follow_a_path_outside_the_study(tmp_path):
    _write(tmp_path, "analysis/analysis_manifest.json", {"results": {
        "../../outside": {"status": "ok"}}})
    payload = problems_payload(tmp_path)
    assert payload["items"] == []
    assert any("unreadable result" in note for note in payload["notes"])
