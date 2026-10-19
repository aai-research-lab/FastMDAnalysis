"""Edited YAML locates actual source keys and uses the shared Config rules."""
from __future__ import annotations

from copy import deepcopy

import pytest
import yaml


def checked(text, **extra):
    from fastmdxplora.gui.config_yaml_edit import validate_yaml
    return validate_yaml({"text": text, **extra})


BASE = "systems:\n  - system: 1UBQ\ninclude_phase: [setup, simulation]\n"


@pytest.mark.parametrize("text,line,code", [
    (BASE + "simulation:\n  temperatur_K: 310\n", 5, "config.option.unknown"),
    (BASE + "simulation:\n  temperature_K: -1\n", 5, "config.option.out_of_range"),
    (BASE + "simulation:\n  temperature_K: cold\n", 5, "config.option.wrong_type"),
    (BASE + "unknown_key: 4\n", 4, "config.option.unknown"),
    (BASE + "execution:\n  workers: many\n", 5, "config.option.wrong_type"),
    (BASE + "simulation:\n  temperature_K: 300\n  temperature_K: 310\n", 6, "config.file.unparseable"),
    ("systems: [\n", 2, "config.file.unparseable"),
])
def test_parser_and_schema_refusals_name_the_real_source_line(text, line, code):
    answer = checked(text)
    assert not answer["ok"]
    diagnostic = answer["diagnostics"][0]
    assert diagnostic["code"] == code
    assert line in [place["line"] for place in diagnostic["locations"]]


def test_nested_system_override_is_checked_without_guessing_a_global_line():
    text = "systems:\n  - system: 1UBQ\n    simulation:\n      temperature_K: -1\nsimulation:\n  temperature_K: 310\n"
    answer = checked(text)
    assert not answer["ok"]
    assert answer["diagnostics"][0]["locations"] == [{"line":4,"column":7,"path":"systems.0.simulation.temperature_K"}]


def test_missing_required_input_is_document_level_instead_of_an_invented_line():
    answer = checked("simulation: {temperature_K: 310}\n")
    assert not answer["ok"] and answer["diagnostics"][0]["locations"] == []


@pytest.mark.parametrize("text", ["[]", "null", "12", "!!python/object/apply:os.system ['false']",
    "systems: &loop [*loop]", "systems: " + "[" * 65 + "1" + "]" * 65,
    "x: " + "a" * 1_000_001])
def test_unsafe_recursive_or_oversized_yaml_is_refused(text):
    assert not checked(text)["ok"]


def test_acyclic_yaml_merges_keep_explicit_overrides():
    text = BASE + "simulation: &settings\n  temperature_K: 310\n  nvt_steps: 100\n"
    assert checked(text)["ok"]
    text = "systems:\n  - system: 1UBQ\n    simulation:\n      <<: &defaults {temperature_K: 300}\n      temperature_K: 315\ninclude_phase: [setup, simulation]\n"
    answer = checked(text)
    assert answer["ok"] and answer["state"]["systems"][0]["simulation"]["temperature_K"] == 315


def _form_of(state):
    # The same request shape as the existing browser mapping; browser tests
    # additionally call its actual snapshotOfLoaded/bodyOf functions.
    result = {"systems": state["systems"], "system": state["system"],
              "include_phase": state["include_phase"], "output": state["output"],
              "__run__": state["study"], "__execution__": state["execution"]}
    result.update(state["phases"])
    return result


def test_roundtrip_preserves_scientific_values_and_per_system_overrides():
    data = {"systems":[{"system":"1UBQ","id":"ubq","simulation":{"temperature_K":315}}],
            "include_phase":["setup","simulation"],"simulation":{"temperature_K":310,"random_seed":42}}
    text = yaml.safe_dump(data)
    answer = checked(text)
    assert answer["ok"]
    assert checked(text, form=_form_of(answer["state"]))["apply_ok"]
    assert data["systems"][0]["simulation"]["temperature_K"] == 315


def test_a_valid_document_that_the_form_would_lose_is_not_applied():
    text = BASE + "analysis:\n  stride: 5\n"
    answer = checked(text)
    assert answer["ok"]
    mapped = checked(text, form=_form_of(answer["state"]))
    assert mapped["ok"] and not mapped["apply_ok"]
    assert "analysis" in mapped["mapping_reason"]


def test_prospective_form_changes_are_detected_before_the_builder_is_mutated():
    text = BASE + "simulation:\n  temperature_K: 310\n"
    answer = checked(text)
    form = _form_of(answer["state"])
    form["simulation"]["temperature_K"] = 330
    mapped = checked(text, form=form)
    assert mapped["ok"] and not mapped["apply_ok"]
    assert "simulation.temperature_K" in mapped["mapping_reason"]
    assert mapped["diagnostics"][0]["locations"][0]["line"] == 5


def test_validation_never_changes_the_request_or_reads_a_named_structure(monkeypatch):
    from pathlib import Path
    from fastmdxplora.gui.config_yaml_edit import validate_yaml
    request = {"text": BASE.replace("1UBQ", "/never-read-this.pdb")}
    before = deepcopy(request)
    original = Path.open
    def guarded(path, *args, **kwargs):
        assert str(path) != "/never-read-this.pdb"
        return original(path, *args, **kwargs)
    monkeypatch.setattr(Path, "open", guarded)
    assert validate_yaml(request)["ok"]
    assert request == before


def test_phase_alias_refusal_points_at_its_original_yaml_key():
    answer = checked("systems: [{system: 1UBQ}]\ninclude: [unknown_phase]\n")
    assert not answer["ok"]
    assert answer["diagnostics"][0]["locations"][0]["line"] == 2


def test_a_builder_conversion_cannot_introduce_unrequested_run_options():
    text = BASE + "simulation: {temperature_K: 310}\n"
    answer = checked(text)
    form = _form_of(answer["state"])
    form["__run__"]["budget_hours"] = 4
    result = checked(text, form=form)
    assert not result["apply_ok"] and "budget_hours" in result["mapping_reason"]


def test_unencodable_yaml_is_refused_instead_of_becoming_an_internal_error():
    assert not checked(BASE + "# \ud800")["ok"]


@pytest.mark.parametrize("public", [False, True])
def test_yaml_validation_route_is_read_only_and_local_control_is_required(tmp_path, public):
    from tests.test_stage_structures_are_recorded import _serving
    from tests.test_a_public_dashboard_answers_only_what_it_lists import _ask
    import json
    if public:
        # A positive control distinguishes this route's public refusal from
        # the blanket refusal that an entirely absent POST route also gets.
        with _serving(tmp_path, public=False) as local:
            status, body = _ask(local + "/api/config-yaml", body={"text":BASE})
            assert status == 200 and json.loads(body)["ok"]
    with _serving(tmp_path, public=public) as url:
        before = set(tmp_path.rglob("*"))
        status, body = _ask(url + "/api/config-yaml", body={"text":BASE})
        assert status == (403 if public else 200)
        if not public:
            assert json.loads(body)["ok"]
        assert set(tmp_path.rglob("*")) == before


@pytest.mark.parametrize("suffix", ["execution: null\n", "simulation: {temperature_K: 310}\n"])
def test_optional_null_execution_can_roundtrip(suffix):
    text = BASE + suffix
    answer = checked(text)
    assert answer["ok"]
    assert checked(text, form=_form_of(answer["state"]))["apply_ok"]


def test_nontext_config_keys_are_refused_without_a_server_error():
    answer = checked(BASE + "100: unknown\n")
    assert not answer["ok"] and answer["diagnostics"][0]["code"]


@pytest.mark.parametrize("value", [".nan", ".inf", "-.inf"])
def test_nonfinite_values_have_finite_json_refusals_and_source_lines(value):
    import json
    answer = checked(BASE + "simulation:\n  temperature_K: " + value + "\n")
    assert not answer["ok"]
    assert answer["diagnostics"][0]["locations"][0]["line"] == 5
    json.dumps(answer, allow_nan=False)


@pytest.mark.parametrize("value", ["2026-20-01", "!!float invalid", "!!int invalid"])
def test_failed_yaml_scalar_construction_is_marked_and_refused(value):
    answer = checked(BASE + "simulation:\n  temperature_K: " + value + "\n")
    assert not answer["ok"] and answer["diagnostics"][0]["locations"][0]["line"] == 5


@pytest.mark.parametrize("value", ["!!omap [{x: .nan}]", "!!pairs [{x: .inf}]", "!!set {x: null}",
                                 "!!binary YQ==", "2026-10-19"])
def test_yaml_only_types_cannot_hide_nonfinite_values_in_a_json_response(value):
    import json
    answer = checked(BASE + "simulation:\n  metadynamics:\n    x: " + value + "\n")
    assert not answer["ok"]
    assert answer["diagnostics"][0]["code"] == "config.option.wrong_type"
    assert answer["diagnostics"][0]["locations"][0]["line"] == 6
    json.dumps(answer, allow_nan=False)


def test_duplicate_merge_keys_are_refused_but_one_merge_with_override_works():
    answer = checked(BASE + "simulation:\n  <<: {temperature_K: 300}\n  <<: {temperature_K: 330}\n")
    assert not answer["ok"]
    assert answer["diagnostics"][0]["locations"][0]["line"] == 6
    assert checked(BASE + "simulation:\n  <<: {temperature_K: 300}\n  temperature_K: 330\n")["ok"]


@pytest.mark.parametrize("location", ["global", "system", "form"])
def test_window_expansion_is_refused_before_any_array_is_allocated(monkeypatch, location):
    import numpy
    monkeypatch.setattr(numpy, "linspace", lambda *args, **kwargs: pytest.fail("An unbounded window array was requested."))
    spec = {"collective_variable":"distance", "selection_a":"name CA", "selection_b":"name N",
            "force_constant":1000,"from":0.3,"to":1.5,"n_windows":1000000000}
    config = {"systems":[{"system":"1UBQ"}],"include_phase":["setup","simulation"]}
    extra = {}
    if location == "system":
        config["systems"][0]["simulation"] = {"umbrella":spec}
    elif location == "global":
        config["simulation"] = {"umbrella":spec}
    else:
        initial = checked(yaml.safe_dump(config))
        form = _form_of(initial["state"])
        form["simulation"] = {"umbrella":spec}
        extra["form"] = form
    answer = checked(yaml.safe_dump(config), **extra)
    assert not answer["ok"] and "1000 umbrella windows" in answer["error"]


def test_a_single_expanded_umbrella_window_is_not_reexpanded_during_override_checks():
    text = BASE + ""
    data = yaml.safe_load(text)
    data["systems"][0]["simulation"] = {"umbrella":{"collective_variable":"distance",
        "selection_a":"name CA","selection_b":"name N","force_constant":1000,"centre":0.4,"index":0}}
    assert checked(yaml.safe_dump(data))["ok"]


@pytest.mark.parametrize("key", [".nan", "2026-10-01", "!!binary YQ=="])
def test_mapping_keys_that_cannot_travel_as_safe_json_are_refused(key):
    import json
    answer = checked(BASE + "simulation:\n  metadynamics:\n    " + key + ": 1\n")
    assert not answer["ok"]
    json.dumps(answer, allow_nan=False)


@pytest.mark.parametrize("name", ["centres", "centers", "n_windows"])
def test_quoted_counts_and_centres_cannot_bypass_the_expansion_bound(monkeypatch, name):
    import numpy
    monkeypatch.setattr(numpy, "linspace", lambda *args, **kwargs: pytest.fail("An array was allocated before refusing invalid windows."))
    spec = {"collective_variable":"distance", "selection_a":"name CA", "selection_b":"name N",
            "force_constant":1000, name:"1" * 5000}
    data = yaml.safe_load(BASE)
    data["simulation"] = {"umbrella":spec}
    assert not checked(yaml.safe_dump(data))["ok"]


@pytest.mark.parametrize("field,value", [("centres",["1e9999","2"]), ("from","1e9999"),
    ("force_constant","1e9999")])
def test_umbrella_reader_cannot_turn_numeric_strings_into_infinities(field, value):
    spec = {"collective_variable":"distance", "selection_a":"name CA", "selection_b":"name N",
            "force_constant":1000,"from":0.3,"to":1.5,"n_windows":3}
    spec[field] = value
    data = yaml.safe_load(BASE)
    data["simulation"] = {"umbrella":spec}
    assert not checked(yaml.safe_dump(data))["ok"]


@pytest.mark.parametrize("analysis", [{"include":["not_an_analysis"]}, {"include":["rmsd"],"exclude":["rmsf"]}])
def test_per_system_analysis_names_and_conflicting_lists_use_the_full_validator(analysis):
    data = yaml.safe_load(BASE)
    data["systems"][0]["analysis"] = analysis
    assert not checked(yaml.safe_dump(data))["ok"]


@pytest.mark.parametrize("phase_line", ["include: [setup, simulation]", "include_phase: setup, simulation", "exclude: [analysis, report]"])
def test_phase_aliases_are_normalized_before_comparing_full_defaults(phase_line):
    text = "systems: [{system: 1UBQ}]\n" + phase_line + "\nsimulation: {temperature_K: 310}\n"
    initial = checked(text)
    assert initial["ok"]
    assert checked(text, form=_form_of(initial["state"]))["apply_ok"]
