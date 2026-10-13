"""A protocol's settings reach only the studies the paper says them of.

A study takes its protocol's settings where it states none of its own. A
paper may describe some runs apart, after the others ("run with a similar
setup as described above; however ..."): a setting whose words are where
the paper describes another protocol's studies, and not where it describes
this one, is not this study's on the paper's word. It is not used; the
person chooses, shown the words of both. A study described where its
settings are, or a paper with one protocol, takes them as before.
"""

from __future__ import annotations

import copy
import json

import pytest
import yaml

from fastmdxplora.config.loader import ConfigError, validate_config
from fastmdxplora.paper.extract import STATUSES, check_reading, kept_reading, tie_to_studies
from fastmdxplora.paper.mapping import plan_study
from fastmdxplora.paper.studies import config_text
from fastmdxplora.paper.text import read_paper
from fastmdxplora.validation.paper_reading import score

from tests._a_paper import PROTOCOL, RESTRAINED, STUDIES, field, jats

MAIN = "Methods: Molecular dynamics simulations"
APART = "Methods: Restrained runs"
AT_310 = "harmonic restraints of 1000 kJ/mol/nm2, at 310 K"

#: The restrained runs, read as their own study under their own protocol,
#: whose settings the reading took from the part on the other runs.
RESTRAINED_STUDY = {
    "id": "S3", "label": "Wild type, restrained", "protocol": "P2",
    "fields": {"system": field("wild type", "Restrained runs of the wild type"),
               "production": field(50, "Restrained runs of the wild type (50 ns)", "ns"),
               "temperature": field(310, AT_310, "K")}}
RESTRAINED_PROTOCOL = {
    "pressure": PROTOCOL["pressure"],
    "timestep": PROTOCOL["timestep"],
    "thermostat": PROTOCOL["thermostat"],
    "engine": PROTOCOL["engine"],
    "equilibration_restraints": field("protein heavy atoms",
                                      "the protein heavy atoms were held by harmonic restraints"),
}


@pytest.fixture(autouse=True)
def _own_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("FASTMDXPLORA_CACHE_DIR", str(tmp_path / "cache"))


@pytest.fixture
def paper(tmp_path):
    path = tmp_path / "paper.xml"
    path.write_bytes(jats(restrained=RESTRAINED))
    return read_paper(path)


def _raw(studies=None, protocols=None):
    return {"studies": copy.deepcopy(STUDIES + [RESTRAINED_STUDY] if studies is None else studies),
            "protocols": [{"id": "P1", "label": "unbiased"}, {"id": "P2", "label": "restrained"}],
            "protocol_fields": copy.deepcopy({"P1": PROTOCOL, "P2": RESTRAINED_PROTOCOL}
                                             if protocols is None else protocols)}


def _study(reading, sid):
    return next(study for study in reading["studies"] if study["id"] == sid)


def test_the_runs_described_apart_do_not_take_the_other_runs_settings(paper):
    reading = check_reading(paper, _raw(), model="x")
    apart = _study(reading, "S3")["fields"]
    for name in ("pressure", "timestep", "thermostat", "engine"):
        assert apart[name]["status"] == "of_another", name
        assert apart[name]["where"] == MAIN
        assert apart[name]["its_study"] == {"quote": "Restrained runs of the wild type",
                                            "where": APART}
    # Its own settings, and its protocol's said where it is described, are used.
    assert apart["production"]["status"] == apart["temperature"]["status"] == "stated"
    assert apart["equilibration_restraints"]["status"] == "stated"
    assert apart["equilibration_restraints"]["from"] == "protocol"
    # The other runs keep their protocol's, said where they are described.
    for sid in ("S1", "S2"):
        fields = _study(reading, sid)["fields"]
        assert fields["pressure"]["status"] == "stated"
        assert fields["pressure"]["from"] == "protocol"
    assert "of_another" in STATUSES


def test_the_person_chooses_shown_the_words_of_both(paper, tmp_path):
    reading = check_reading(paper, _raw(), model="x")
    plan = plan_study(_study(reading, "S3"), reading)
    assert plan["state"] == "needs_you"
    asked = {c["field"]: c for c in plan["choices"] if c["label"] == "needs_you"}
    for name in ("pressure", "timestep", "thermostat"):
        why = asked[name]["why"]
        assert MAIN in why and APART in why, why
        assert "another protocol's studies" in why
        assert "if it holds for this study, set it in `" in why
        assert "if not, set this study's value instead" in why
        assert asked[name]["quote"] and asked[name]["where"] == MAIN
    assert "1 bar" in asked["pressure"]["why"]
    assert "`simulation.pressure_bar`" in asked["pressure"]["why"]
    assert "`simulation.timestep_fs`" in asked["timestep"]["why"]
    # The study's own words are said in full once, not with every setting.
    whys = [c["why"] for c in plan["choices"] if APART in c["why"]]
    assert sum("Restrained runs of the wild type" in why for why in whys) == 1
    # Each said once: not also as this software's own value, not stated.
    assert [c["field"] for c in plan["choices"]].count("pressure") == 1
    # The engine is not set here: said, and not asked.
    engine = next(c for c in plan["choices"] if c["field"] == "engine")
    assert engine["label"] == "not_stated" and "not used" in engine["why"]
    # The config says why, and is refused until the person has chosen.
    text = config_text(plan)
    assert "where it describes another protocol's studies" in text
    with pytest.raises(ConfigError):
        validate_config(yaml.safe_load(text), require_systems=True)
    # The runs described with their settings need nothing for them.
    other = plan_study(_study(reading, "S1"), reading)
    assert not [c for c in other["choices"] if c["label"] == "needs_you"]


def _record(value, quote, where, given_by):
    return {"status": "stated", "value": value, "quote": quote, "where": where,
            "from": given_by}


def _apart(**fields):
    """The other runs, and a study described apart whose ``fields`` are
    given as checked records, marked as a reading is."""
    other = {"id": "S1", "label": "Wild type", "protocol": "P1",
             "fields": {"system": _record("wild type", "wild type", MAIN, "study"),
                        "pdb_id": _record("1UBQ", "PDB ID: 1UBQ", MAIN, "study")}}
    apart = {"id": "S2", "label": "Wild type, apart", "protocol": "P2",
             "fields": {"system": _record("wild type", "wild type", APART, "study"),
                        "pdb_id": _record("1UBQ", "PDB ID: 1UBQ", APART, "study"), **fields}}
    studies = tie_to_studies([other, apart])
    return studies[1]


def test_a_choice_the_study_s_own_words_led_to_is_kept():
    # Its own details name a method this software does not run; its
    # protocol's method, said of the other runs, does not make it runnable.
    details = _record("temperature replica exchange MD with 24 replicas",
                      "we ran temperature replica exchange MD (REMD) with 24 replicas",
                      APART, "study")
    method = _record("plain MD", "molecular dynamics", MAIN, "protocol")
    study = _apart(method_details=details, method=method)
    assert study["fields"]["method"]["status"] == "of_another"
    plan = plan_study(study, {})
    assert plan["state"] == "cannot_run"
    assert ("method", "not_possible") in {(c["field"], c["label"]) for c in plan["choices"]}


def test_a_question_its_own_step_raises_is_kept():
    # Its own 4 fs step needs heavier hydrogens; the mass, given where the
    # other runs are described, is asked, and why the step needs it stays.
    study = _apart(timestep=_record(4, "a 4 fs time step", APART, "study"),
                   hydrogen_mass=_record(3.024, "hydrogen masses of 3.024 amu", MAIN,
                                         "protocol"))
    assert study["fields"]["hydrogen_mass"]["status"] == "of_another"
    plan = plan_study(study, {})
    whys = [c["why"] for c in plan["choices"]
            if c["field"] == "hydrogen_mass" and c["label"] == "needs_you"]
    assert len(whys) == 2
    assert any("A 4 fs step needs hydrogen mass repartitioning" in why for why in whys)
    assert any("`setup.hydrogen_mass_amu`" in why and MAIN in why for why in whys)


def test_a_lipid_force_field_given_for_other_runs_is_asked():
    study = _apart(membrane=_record("POPC bilayer", "a POPC bilayer", APART, "study"),
                   protein_forcefield=_record("ff14SB", "AMBER ff14SB", APART, "study"),
                   lipid_forcefield=_record("Lipid21", "Lipid21 for the lipids", MAIN,
                                            "protocol"))
    assert study["fields"]["lipid_forcefield"]["status"] == "of_another"
    plan = plan_study(study, {})
    asked = [c for c in plan["choices"] if c["field"] == "lipid_forcefield"]
    assert [c["label"] for c in asked] == ["needs_you"]
    assert "`setup.force_field` (the force field's files)" in asked[0]["why"]
    assert "Lipid21" in asked[0]["why"]
    assert plan["state"] == "needs_you"


def test_the_study_s_own_reference_is_said_in_place_of_another_s_setting(paper):
    study = copy.deepcopy(RESTRAINED_STUDY)
    study["fields"]["timestep"] = {"status": "by_reference",
                                   "quote": "a similar setup as described above"}
    reading = check_reading(paper, _raw(STUDIES + [study]), model="x")
    timestep = _study(reading, "S3")["fields"]["timestep"]
    assert timestep["status"] == "by_reference" and timestep["from"] == "study"
    assert timestep["where"] == APART
    assert tie_to_studies(copy.deepcopy(reading["studies"])) == reading["studies"]


def test_a_setting_said_where_the_study_is_described_reaches_it(paper):
    protocols = {"P1": PROTOCOL,
                 "P2": {**RESTRAINED_PROTOCOL, "temperature": field(310, AT_310, "K")}}
    study = copy.deepcopy(RESTRAINED_STUDY)
    del study["fields"]["temperature"]
    reading = check_reading(paper, _raw(STUDIES + [study], protocols), model="x")
    temperature = _study(reading, "S3")["fields"]["temperature"]
    assert temperature["status"] == "stated" and temperature["from"] == "protocol"
    assert temperature["where"] == APART


def test_with_one_protocol_its_settings_reach_every_study(paper):
    # The study told apart in the results, its settings in the methods:
    # no other protocol's studies are described there, so it takes them.
    told = {"id": "S3", "label": "Wild type, again", "protocol": "P1",
            "fields": {"system": field("wild type",
                                       "The root-mean-square deviation of the wild type")}}
    raw = _raw(STUDIES + [told], {"P1": PROTOCOL})
    raw["protocols"] = raw["protocols"][:1]
    reading = check_reading(paper, raw, model="x")
    fields = _study(reading, "S3")["fields"]
    assert fields["system"]["where"] == "Results"
    assert all(record["status"] == "stated" for name, record in fields.items()
               if record.get("from") == "protocol")
    assert fields["pressure"]["status"] == "stated"


def test_a_study_with_none_of_its_own_words_takes_its_protocol_s(paper):
    bare = {"id": "S3", "label": "Wild type, restrained", "protocol": "P2", "fields": {}}
    reading = check_reading(paper, _raw(STUDIES + [bare]), model="x")
    assert _study(reading, "S3")["fields"]["pressure"]["status"] == "stated"


def test_a_reading_kept_from_before_is_marked_when_used(paper):
    reading = check_reading(paper, _raw(), model="x")
    # As a reading kept before the marking had it.
    before = copy.deepcopy(reading)
    for study in before["studies"]:
        for record in study["fields"].values():
            if record.get("status") == "of_another":
                record["status"] = "stated"
                del record["its_study"]
    from fastmdxplora.paper.extract import _key
    path = _key(paper, "x")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(before), encoding="utf-8")
    kept = kept_reading(paper, "x")
    assert _study(kept, "S3")["fields"]["pressure"]["status"] == "of_another"
    assert kept["studies"] == reading["studies"]
    # Marking twice changes nothing.
    again = copy.deepcopy(kept["studies"])
    assert tie_to_studies(again) == kept["studies"]


def test_a_setting_of_another_study_is_not_counted_as_used(paper):
    reading = check_reading(paper, _raw(), model="x")
    apart = _study(reading, "S3")
    truth = {"key": "test", "studies": [
        {"label": apart["label"],
         "fields": {"production": {"status": "stated", "value": 50},
                    "temperature": {"status": "stated", "value": 310},
                    "pressure": {"status": "not_stated"},
                    "timestep": {"status": "not_stated"},
                    "thermostat": {"status": "not_stated"}}}]}
    counted = score(truth, {"studies": [apart]})
    assert counted["fields"]["invented"] == 0
    assert counted["fields"]["correct"] == 2
