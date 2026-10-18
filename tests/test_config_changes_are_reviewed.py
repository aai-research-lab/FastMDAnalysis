"""A reviewed version contains only the selected values and their reasons."""
from copy import deepcopy

import pytest


BASE = {
    "systems": [{"id": "ubq", "system": "1UBQ"}],
    "simulation": {"temperature_K": 300, "duration_ns": 2},
    "decisions": {"simulation.temperature_K": {"why": "Original temperature.", "source": "person"}},
}


def _review(before, after, accepted=None):
    from fastmdxplora.gui.config_review import review_config

    body = {"before": before, "after": after}
    if accepted is not None:
        body["accepted"] = accepted
    return review_config(body)


def _choose(before, after, settings):
    rows = _review(before, after)
    assert rows["ok"], rows
    ids = [r["id"] for r in rows["changes"] if r["setting"] in settings]
    return _review(before, after, ids)


def test_review_derives_changes_and_keeps_both_originals():
    after = deepcopy(BASE)
    after["simulation"].update(temperature_K=330, duration_ns=10)
    after["decisions"]["simulation.temperature_K"] = {"why": "Heat test.", "source": "agent"}
    old, new = deepcopy(BASE), deepcopy(after)
    result = _choose(BASE, after, {"simulation.duration_ns"})
    assert result["ok"]
    assert result["config"]["simulation"] == {"temperature_K": 300, "duration_ns": 10}
    assert result["config"]["decisions"] == BASE["decisions"]
    assert result["review"] == {"accepted": ["simulation.duration_ns"], "rejected": ["simulation.temperature_K"],
                                "said": "Accepted: Production 2 ns \u2192 10 ns. Rejected: Temperature 300 K \u2192 330 K."}
    assert BASE == old and after == new


def test_accepted_provenance_travels_with_its_value():
    after = deepcopy(BASE)
    after["simulation"]["temperature_K"] = 330
    after["decisions"]["simulation.temperature_K"] = {"why": "As requested.", "source": "person"}
    result = _choose(BASE, after, {"simulation.temperature_K"})
    assert result["config"] == after
    assert "330 K" in str(result["plan"])


def test_rejecting_an_added_setting_drops_its_decision():
    after = deepcopy(BASE)
    after["setup"] = {"ph": 6.5}
    after["decisions"]["setup.ph"] = {"why": "A test.", "source": "agent"}
    result = _choose(BASE, after, set())
    assert result["ok"] and result["config"] == BASE


def test_accepting_a_removed_setting_removes_its_decision():
    after = deepcopy(BASE)
    del after["simulation"]["temperature_K"]
    after.pop("decisions")
    result = _choose(BASE, after, {"simulation.temperature_K"})
    assert result["ok"]
    assert "temperature_K" not in result["config"]["simulation"]
    assert "decisions" not in result["config"]


def test_a_reason_only_change_can_be_rejected():
    after = deepcopy(BASE)
    after["decisions"]["simulation.temperature_K"]["why"] = "A corrected explanation."
    rows = _review(BASE, after)
    assert len(rows["changes"]) == 1
    assert rows["changes"][0]["before_decision"]["why"] == "Original temperature."
    assert _choose(BASE, after, set())["config"] == BASE
    assert _choose(BASE, after, {"simulation.temperature_K"})["config"] == after


def test_an_unset_default_can_have_its_reason_reviewed():
    after = deepcopy(BASE)
    after["decisions"]["setup.ph"] = {"why": "Use the default pH.", "source": "person"}
    rows = _review(BASE, after)
    assert [r["setting"] for r in rows["changes"]] == ["setup.ph"]
    assert _choose(BASE, after, {"setup.ph"})["config"]["decisions"] == after["decisions"]


def test_systems_and_sweeps_are_atomic_choices():
    after = deepcopy(BASE)
    after["systems"] = [{"system": "1L2Y", "setup": {"ph": 6}}]
    after["sweep"] = {"simulation.temperature_K": [300, 310]}
    result = _choose(BASE, after, {"sweep"})
    assert result["ok"]
    assert result["config"]["systems"] == BASE["systems"]
    assert result["config"]["sweep"] == after["sweep"]


def test_compound_settings_are_not_split_into_unsafe_fragments():
    before = deepcopy(BASE)
    before["simulation"]["stop_when"] = {"max_duration_ns": 10, "independent_starts": "not_required", "measures": [{"analysis": "rg", "standard_error": 0.1}]}
    after = deepcopy(before)
    after["simulation"]["stop_when"]["measures"][0]["standard_error"] = 0.2
    rows = _review(before, after)
    assert rows["ok"] and [r["setting"] for r in rows["changes"]] == ["simulation.stop_when"]
    assert _choose(before, after, set())["config"] == before


def test_a_combination_is_checked_not_just_its_two_valid_parents():
    before = deepcopy(BASE)
    before["analysis"] = {"include": ["rmsd"]}
    after = deepcopy(BASE)
    after["analysis"] = {"exclude": ["rg"]}
    result = _choose(before, after, {"analysis.exclude"})
    assert not result["ok"]
    assert result["refusal"]["code"] == "config.option.conflicting"
    assert "config" not in result and "yaml" not in result


@pytest.mark.parametrize("accepted", [["simulation.temperature_K"], ["0", "0"], [0], "0"])
def test_choices_cannot_supply_paths_or_unknown_ids(accepted):
    after = deepcopy(BASE)
    after["simulation"]["temperature_K"] = 330
    result = _review(BASE, after, accepted)
    assert not result["ok"] and "config" not in result


def test_agent_metadata_is_kept_without_becoming_a_scientific_choice():
    after = deepcopy(BASE)
    after.update(agent="assisted", agent_model="provider/name")
    after["simulation"]["temperature_K"] = 330
    result = _choose(BASE, after, set())
    assert result["ok"] and result["config"]["simulation"] == BASE["simulation"]
    assert result["config"]["agent"] == "assisted"
    assert result["config"]["agent_model"] == "provider/name"
    assert result["review"]["rejected"] == ["simulation.temperature_K"]


@pytest.mark.parametrize("payload", [None, [], {}, {"before": [], "after": BASE},
                                    {"before": BASE, "after": {"simulation": {"temperature_K": "hot"}}}])
def test_review_refuses_invalid_configs_before_offering_choices(payload):
    from fastmdxplora.gui.config_review import review_config

    result = review_config(payload)
    assert not result["ok"] and "config" not in result


def test_the_endpoint_checks_without_writing_or_starting_a_study(tmp_path, monkeypatch):
    import json
    import urllib.request

    from fastmdxplora.gui.exploration import DashboardRuntime
    from fastmdxplora.gui.server import start_dashboard_session

    def no_run(*args, **kwargs):
        raise AssertionError("Review started a study")

    monkeypatch.setattr(DashboardRuntime, "launch_from_config", no_run)
    monkeypatch.setenv("FASTMDXPLORA_CONFIG_DIR", str(tmp_path / "config"))
    session = start_dashboard_session(output=str(tmp_path / "study"), host="127.0.0.1", port=0)
    try:
        before_files = {p.relative_to(tmp_path) for p in tmp_path.rglob("*") if p.is_file()}
        after = deepcopy(BASE)
        after["simulation"]["temperature_K"] = 330
        request = urllib.request.Request(session.url + "/api/agent/review-config",
            data=json.dumps({"before": BASE, "after": after, "accepted": ["0"]}).encode(),
            headers={"Origin": session.url, "Content-Type": "application/json"})
        with urllib.request.urlopen(request, timeout=10) as response:
            result = json.load(response)
        assert result["ok"] and result["config"]["simulation"]["temperature_K"] == 330
        assert {p.relative_to(tmp_path) for p in tmp_path.rglob("*") if p.is_file()} == before_files
    finally:
        session.stop()


def test_accepting_a_setting_in_a_new_phase_creates_that_phase():
    after = deepcopy(BASE)
    after["setup"] = {"ph": 6.5}
    result = _choose(BASE, after, {"setup.ph"})
    assert result["ok"] and result["config"] == after


def test_accepting_the_last_removed_setting_removes_its_absent_phase():
    before = deepcopy(BASE)
    before["setup"] = {"ph": 6.5}
    result = _choose(before, BASE, {"setup.ph"})
    assert result["ok"] and result["config"] == BASE


def test_a_writer_record_is_not_a_rejectable_scientific_change():
    before = deepcopy(BASE)
    before.update(agent="assisted")
    before["decisions"]["agent"] = {"why": "Assisted draft.", "source": "person"}
    after = deepcopy(before)
    after.update(agent="autonomous")
    after["decisions"]["agent"]["why"] = "Autonomous draft."
    rows = _review(before, after)
    assert rows["ok"] and rows["changes"] == []
    result = _review(before, after, [])
    assert result["config"]["agent"] == "autonomous"
    assert result["config"]["decisions"]["agent"]["why"] == "Autonomous draft."


@pytest.mark.parametrize("phase, established, general, old, new", [
    ("simulation", "prepared_from", "setup_from", "runs/old", "runs/new"),
    ("setup", "ligand_name", "ligand_resname", "OLD", "NEW"),
    ("analysis", "selection", "select_atoms", "protein", "backbone"),
])
def test_mirrored_names_are_one_choice_that_survives_validation(phase, established, general, old, new):
    from fastmdxplora.config.loader import validate_config

    before = deepcopy(BASE)
    before.setdefault(phase, {}).update({established: old, general: old})
    after = deepcopy(before)
    after[phase].update({established: new, general: new})
    after["decisions"][f"{phase}.{general}"] = {"why": "The new choice.", "source": "person"}
    rows = _review(before, after)
    assert rows["ok"] and len(rows["changes"]) == 1
    row = rows["changes"][0]
    assert row["aliases"] == [f"{phase}.{established}", f"{phase}.{general}"]
    assert len(row["alias_labels"]) == 2
    accepted = _review(before, after, [row["id"]])
    rejected = _review(before, after, [])
    assert accepted["config"][phase][established] == accepted["config"][phase][general] == new
    assert rejected["config"][phase][established] == rejected["config"][phase][general] == old
    validate_config(accepted["config"], require_systems=True)
    assert accepted["config"][phase][established] == new
    assert accepted["config"]["decisions"][f"{phase}.{general}"]["why"] == "The new choice."
    assert f"{phase}.{general}" not in rejected["config"]["decisions"]


def test_earlier_phase_list_spelling_cannot_override_a_review_choice():
    before = dict(deepcopy(BASE), include=["simulation"])
    after = dict(deepcopy(BASE), include_phase=["analysis"])
    rows = _review(before, after)
    assert rows["ok"] and [r["setting"] for r in rows["changes"]] == ["include_phase"]
    result = _review(before, after, ["0"])
    assert result["config"]["include_phase"] == ["analysis"]
    assert "include" not in result["config"]


def test_a_proposal_without_writer_metadata_does_not_retain_the_previous_writers():
    before = dict(deepcopy(BASE), agent="autonomous", agent_model="provider/old")
    before["decisions"]["agent"] = {"why": "Old draft.", "source": "person"}
    result = _review(before, BASE, [])
    assert result["ok"] and result["config"] == BASE


def test_review_requires_a_complete_study_as_an_agent_proposal_does():
    result = _review(BASE, {}, ["0"])
    assert not result["ok"] and "config" not in result
    assert result["refusal"]["code"] == "config.option.missing_companion"


@pytest.mark.parametrize("full", [False, True])
def test_export_keeps_chosen_values_and_decisions_from_an_inactive_phase(full):
    from fastmdxplora.gui.config_builder import export_config
    import yaml

    config = deepcopy(BASE)
    config.update(include_phase=["simulation"], report={"title": "Keep this accepted title."})
    config["decisions"]["report.title"] = {"why": "Chosen while planning.", "source": "person"}
    original = deepcopy(config)
    result = export_config(config, full=full)
    assert result["ok"], result
    exported = yaml.safe_load(result["yaml"])
    assert exported["report"] == config["report"]
    assert exported["decisions"] == config["decisions"]
    assert exported["simulation"]["temperature_K"] == 300
    assert exported["include_phase"] == ["simulation"]
    assert config == original
    if full:
        assert "timestep_fs" in exported["simulation"]


@pytest.mark.parametrize("spelling", ["prepared_from", "setup_from"])
def test_a_reviewed_version_holds_only_the_spelling_each_version_used(spelling):
    before = deepcopy(BASE)
    before["simulation"][spelling] = "runs/old"
    after = deepcopy(before)
    after["simulation"].update({spelling: "runs/new", "temperature_K": 330})
    rows = _review(before, after)["changes"]
    shown = {r["label"] for r in rows}
    assert shown == {"Temperature", "Prepared from"}
    assert all("alias_labels" not in r for r in rows)
    rejected = _review(before, after, [])
    assert rejected["config"] == before
    assert {c["label"] for c in rejected["changes"]} == shown
    accepted = _review(before, after, [r["id"] for r in rows])
    assert accepted["config"] == after and accepted["changes"] == []
    mixed = _choose(before, after, {"simulation.prepared_from"})
    assert mixed["config"]["simulation"] == {"temperature_K": 300, "duration_ns": 2, spelling: "runs/new"}
    assert [c["label"] for c in mixed["changes"]] == ["Temperature"]


def test_values_and_reasons_are_said_as_the_builder_says_them():
    after = deepcopy(BASE)
    after["simulation"]["temperature_K"] = 330
    after["sweep"] = {"simulation.random_seed": [1, 2, 3]}
    after["decisions"]["simulation.temperature_K"] = {
        "why": "Warmer.", "source": "agent", "alternatives": [310, 350]}
    rows = {r["setting"]: r for r in _review(BASE, after)["changes"]}
    sweep = rows["sweep"]
    assert (sweep["label"], sweep["before"], sweep["after"]) == ("Varied random seed", "not set", "1, 2, 3")
    temperature = rows["simulation.temperature_K"]
    assert temperature["before_reason"] == "Original temperature. (you asked)"
    assert temperature["after_reason"] == "Warmer. (Agent’s choice; alternatives: 310 K, 350 K)"


def test_a_compound_change_is_said_part_by_part():
    after = deepcopy(BASE)
    after["sweep"] = {"simulation.random_seed": [1, 2], "simulation.temperature_K": [300, 310]}
    row = _review(BASE, after)["changes"][0]
    assert row["label"] == "Sweep"
    assert row["parts"] == [
        {"label": "Varied random seed", "before": "not set", "after": "1, 2"},
        {"label": "Varied temperature", "before": "not set", "after": "300, 310"}]
    assert "{" not in row["before"] + row["after"]


def test_a_reason_only_change_says_the_value_is_the_same():
    after = deepcopy(BASE)
    after["decisions"]["simulation.temperature_K"]["why"] = "A corrected explanation."
    row = _review(BASE, after)["changes"][0]
    assert row["same_value"] is True and row["before"] == row["after"] == "300 K"


def test_the_review_record_says_the_choices_as_they_were_shown():
    after = deepcopy(BASE)
    after["simulation"].update(temperature_K=330, duration_ns=10)
    result = _choose(BASE, after, {"simulation.duration_ns"})
    assert result["review"]["said"] == "Accepted: Production 2 ns → 10 ns. Rejected: Temperature 300 K → 330 K."
    assert "simulation." not in result["review"]["said"]
    assert _choose(BASE, after, set())["review"]["said"] == (
        "Rejected: Temperature 300 K → 330 K; Production 2 ns → 10 ns.")


def test_a_refused_combination_names_the_rows_to_change():
    before = deepcopy(BASE)
    before["analysis"] = {"include": ["rmsd"]}
    after = deepcopy(BASE)
    after["analysis"] = {"exclude": ["rg"]}
    after["simulation"]["temperature_K"] = 330
    result = _choose(before, after, {"analysis.exclude", "simulation.temperature_K"})
    assert not result["ok"]
    involved = [r["label"] for r in result["changes"] if r["id"] in result["involved"]]
    assert involved == ["Analyses left out", "Analyses"]
