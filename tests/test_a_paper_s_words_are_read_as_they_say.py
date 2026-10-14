"""A paper's words are read as they say, in the ways the first counted
reading (check 5) found they were not.

Check 5's readings, read value by value, gave words that hold a value this
software did not read from them: "particle-mesh Ewald" taken as two
methods, "CHARMM22 with CMAP" as CHARMM22, ten PDB entries for one study,
"set the net system charge to zero", "in pure water", "400 ns per Voronoi
cell", "forty-eight trajectories". Each is read here as it says. A reading
kept from before is checked again when it is used, so a value now read is
used without asking anything again.
"""

from __future__ import annotations

import copy
import json

import pytest

from fastmdxplora.paper.extract import _key, check_field, check_reading, kept_reading
from fastmdxplora.paper.mapping import plan_study, recognized
from fastmdxplora.paper.quotes import QuoteIndex
from fastmdxplora.paper.text import read_paper
from fastmdxplora.validation.paper_reading import same_value

from tests._a_paper import PROTOCOL, STUDIES, field, jats

#: Methods made of words check 5's papers used, one sentence each.
METHODS = (
    "Long range electrostatic interactions were calculated using the particle-mesh "
    "Ewald method. Long-range electrostatics were treated with the Gaussian split Ewald "
    "method. The protein was described by CHARMM22 with the CMAP correction. The initial "
    "conformations of the complexes were taken from holo crystal structures (3U5J, ; 3U5L, "
    "; 4OGI). One sodium ion was added to set the net system charge to zero. We elected "
    "to perform all MD simulations in pure water. MMVT simulations were then run for 400 "
    "ns per Voronoi cell. Four trajectories of 1 us and forty-eight trajectories of 2 us "
    "were started."
)


@pytest.fixture(autouse=True)
def _own_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("FASTMDXPLORA_CACHE_DIR", str(tmp_path / "cache"))


@pytest.fixture
def paper(tmp_path):
    path = tmp_path / "paper.xml"
    path.write_bytes(jats(METHODS))
    return read_paper(path)


def _checked(paper, name, value, quote, unit=None):
    return check_field(QuoteIndex(paper), name, field(value, quote, unit))


def test_particle_mesh_ewald_is_one_method():
    assert recognized("electrostatics", "particle-mesh Ewald") == {"PME"}
    assert recognized("electrostatics", "smooth particle mesh Ewald") == {"PME"}
    assert recognized("electrostatics", "Gaussian split Ewald") == {"Gaussian split Ewald"}
    assert recognized("electrostatics", "Ewald summation") == {"Ewald"}
    # The harness counts the same method by either name.
    assert same_value("electrostatics", "PME", "particle-mesh Ewald")


def test_charmm22_with_cmap_is_charmm27(paper):
    for said in ("CHARMM22 with the CMAP correction", "CHARMM22/CMAP", "CHARMM27"):
        assert recognized("protein_forcefield", said) == {"CHARMM27"}, said
    assert recognized("protein_forcefield", "CHARMM22") == {"CHARMM22"}
    assert recognized("protein_forcefield", "CHARMM22*") == {"CHARMM22*"}
    record = _checked(paper, "protein_forcefield", "CHARMM27",
                      "described by CHARMM22 with the CMAP correction")
    assert record["status"] == "stated"
    record = _checked(paper, "protein_forcefield", "CHARMM22",
                      "described by CHARMM22 with the CMAP correction")
    assert record["status"] == "unread"


def test_a_study_s_several_pdb_entries_are_each_read(paper):
    quote = "taken from holo crystal structures (3U5J, ; 3U5L, ; 4OGI)"
    record = _checked(paper, "pdb_id", "3U5J, 3U5L, 4OGI", quote)
    assert record["status"] == "stated" and record["value"] == "3U5J, 3U5L, 4OGI"
    assert _checked(paper, "pdb_id", ["3U5J", "3U5L"], quote)["value"] == "3U5J, 3U5L"
    # One entry its words do not hold, and none is used.
    assert _checked(paper, "pdb_id", "3U5J, 1UBQ", quote)["status"] == "unread"
    # Run, it asks which: one config runs one structure.
    study = {"id": "S1", "label": "complexes", "fields": {"pdb_id": record}}
    plan = plan_study(study, {})
    choice = next(c for c in plan["choices"] if c["field"] == "pdb_id")
    assert choice["label"] == "needs_you"
    assert "3 PDB entries (3U5J, 3U5L, 4OGI)" in choice["why"]


def test_a_charge_set_to_zero_is_neutralised(paper):
    record = _checked(paper, "neutralized", "yes",
                      "One sodium ion was added to set the net system charge to zero")
    assert record["status"] == "stated" and record["value"] is True


def test_pure_water_is_no_salt(paper):
    record = _checked(paper, "salt_concentration", "pure water",
                      "We elected to perform all MD simulations in pure water")
    assert record["status"] == "stated" and record["value"] == 0.0


def test_a_length_per_cell_is_a_length(paper):
    record = _checked(paper, "production", 400, "MMVT simulations were then run for 400 ns "
                      "per Voronoi cell", "ns per Voronoi cell")
    assert record["status"] == "stated" and record["value"] == 400.0


def test_a_count_in_words_with_tens_and_units(paper):
    quote = "Four trajectories of 1 us and forty-eight trajectories of 2 us were started"
    record = _checked(paper, "replicas", "forty-eight", quote)
    assert record["status"] == "stated" and record["value"] == 48
    assert _checked(paper, "replicas", 4, quote)["value"] == 4
    # Two counts in the words, none given: neither is chosen.
    assert _checked(paper, "replicas", None, quote)["status"] == "unread"


def test_a_kept_reading_is_checked_again_when_used(paper):
    studies = copy.deepcopy(STUDIES)
    studies[0]["fields"]["neutralized"] = field(
        "yes", "One sodium ion was added to set the net system charge to zero")
    reading = check_reading(paper, {"studies": studies, "protocol_fields": {"P1": PROTOCOL}},
                            model="x")
    # As kept before the words were read so: unread, with what was said.
    before = copy.deepcopy(reading)
    before["studies"][0]["fields"]["neutralized"] = {
        "status": "unread", "quote": "One sodium ion was added to set the net system charge "
        "to zero", "where": "Methods: Molecular dynamics simulations", "in": "paper",
        "said": {"value": "yes"}, "from": "study"}
    before["studies"][0]["fields"]["water_model"] = {
        "status": "not_found", "said": {"value": "TIP3P", "quote": "TIP3P 44"},
        "from": "protocol"}
    path = _key(paper, "x")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(before), encoding="utf-8")
    kept = kept_reading(paper, "x")
    fields = kept["studies"][0]["fields"]
    assert fields["neutralized"]["status"] == "stated" and fields["neutralized"]["value"] is True
    assert fields["neutralized"]["from"] == "study"
    # Words still not found stay as they were.
    assert fields["water_model"]["status"] == "not_found"
    # Checking again changes nothing more.
    path.write_text(json.dumps(kept), encoding="utf-8")
    assert kept_reading(paper, "x") == kept


def test_what_only_looks_like_these_is_not_read_as_them(paper):
    from fastmdxplora.paper.extract import _count_read, _flag_read
    from fastmdxplora.paper.values import read_value

    # A charge that is not zero, or said not to be.
    assert _flag_read("the net charge was 0.5 e", "neutralized") is None
    assert _flag_read("a net charge of zero was not achieved", "neutralized") is False
    # A "not" about something else leaves the neutralising as it was.
    for words in ("Na+ counterions were added to neutralize the charge; ions were not "
                  "restrained", "counterions were added to neutralize the system charge, "
                  "which was not zero initially", "Counter-ions were added to compensate the "
                  "protein charge (not shown)"):
        assert _flag_read(words, "neutralized") is True, words
    # Twenty runs of one microsecond are not twenty-one runs.
    assert _count_read("twenty one-microsecond runs were performed", "replicas", None) != 21
    assert _count_read("twenty one-microsecond runs were performed", "replicas", 20) == 20
    # Pure water with salt added is not without salt.
    for words in ("solvated in pure water and 150 mM NaCl was added",
                  "in pure water containing 0.15 M NaCl", "in pure water with 0.1 M KCl",
                  "in pure water and 0.15 M sodium chloride",
                  "solvated in pure water, then 0.15 M KCl added"):
        assert read_value(words, "concentration", 0, "M") is None, words
    # Ewald and CMAP however the words are broken or put.
    for words in ("particle mesh\nEwald", "particle\u00a0mesh\u00a0Ewald"):
        assert recognized("electrostatics", words) == {"PME"}
    assert recognized("electrostatics", "Gaussian split\nEwald") == {"Gaussian split Ewald"}
    for words in ("CHARMM22 (with CMAP)", "CHARMM22 force field with the CMAP correction",
                  "CHARMM22 plus CMAP"):
        assert recognized("protein_forcefield", words) == {"CHARMM27"}, words
    assert recognized("protein_forcefield", "CHARMM22 without CMAP") == {"CHARMM22"}


def test_coupled_replicas_and_walkers_are_not_independent_runs():
    from fastmdxplora.paper.extract import _count_read

    # Walkers share one bias; bias-exchange and temperature replicas
    # exchange: each is a part of one run, not a run.
    for words, given in (("we perform an OPES simulation using 32 walkers, biasing z", 32),
                         ("Simulations were run with a BEMD scheme of three replicas", 3),
                         ("we used R = 40 replicas distributed exponentially in the "
                          "temperature range from 300 to 440 K", 40),
                         ("REMD with 32 replicas was run at 300-450 K", 32),
                         ("Multiple-walker metadynamics with 32 walkers was run", 32),
                         ("32 replicas were run with exchanges every 2 ps", 32),
                         ("three replicas each with exchange attempts every 2 ps", 3),
                         ("a 24-replica REMD simulation", 24),
                         ("REST2 with 8 replicas", 8),
                         ("HREX with 16 replicas was run in triplicate", 16)):
        assert _count_read(words, "replicas", given) is None, words
    assert _count_read("HREX with 16 replicas was run in triplicate", "replicas", 3) == 3
    # "Rest", or an exchange of something else, is no coupled method.
    for words in ("Production consisted of three replicas, the rest of the protocol as above",
                  "simulated in three replicas; the ligand was at rest in the pocket",
                  "three replicas of 200 ns, with water exchange monitored",
                  "Three replicas were heated from 0 to 300 K over 1 ns"):
        assert _count_read(words, "replicas", 3) == 3, words
    # Separate runs of such a method are runs.
    assert _count_read("For each simulation system we have produced four REMD trajectories",
                       "replicas", 4) == 4
    assert _count_read("we ran three independent bias-exchange simulations", "replicas", 3) == 3
    assert _count_read("Production consisted of three replicas of 100 ns", "replicas", 3) == 3


def test_a_kept_quote_not_found_is_not_found_again(paper):
    studies = copy.deepcopy(STUDIES)
    reading = check_reading(paper, {"studies": studies, "protocol_fields": {"P1": PROTOCOL}},
                            model="x")
    # Kept cut short: its start is the paper's, its end was not.
    reading["studies"][0]["fields"]["production"] = {
        "status": "not_found", "from": "study",
        "said": {"value": 100, "unit": "ns", "quote": "three independent 100-ns simulations"}}
    path = _key(paper, "x")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(reading), encoding="utf-8")
    kept = kept_reading(paper, "x")
    assert kept["studies"][0]["fields"]["production"]["status"] == "not_found"
