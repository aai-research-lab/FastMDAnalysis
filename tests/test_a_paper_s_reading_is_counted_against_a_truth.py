"""How well a paper's MD studies are read, counted against a truth made by hand.

Two readers read each registered paper apart; what both state with one value
is the truth, the rest contested and left out of the count. A reading is
counted field by field: the value used is the truth's, another one, none
where the truth states one, or one where the truth says the paper gives none.
The truth keeps values and where they are, never the paper's words.
"""

from __future__ import annotations

import json

import pytest

from fastmdxplora.paper.fields import BY_NAME
from fastmdxplora.validation import paper_reading
from fastmdxplora.validation.paper_reading import (PAPERS, TRUTH, main, match_studies, reconcile,
                                                   same_value, score)


def _stated(value, unit=""):
    return {"status": "stated", "value": value, "unit": unit, "where": "Methods",
            "quote": "the paper's words"}


def _study(sid, label, **fields):
    return {"id": sid, "label": label, "fields": fields, "claims": []}


@pytest.mark.parametrize("name, a, b, same", [
    ("temperature", 300, 300.0, True),
    ("temperature", 300, 310, False),
    ("timestep", 2, "2.0", True),
    ("protein_forcefield", "AMBER ff14SB", "ff14SB", True),
    ("protein_forcefield", "ff14SB", "CHARMM36m", False),
    ("water_model", "TIP3P", "tip3p", True),
    ("mutations", "L50A, K63R", ["K63R", "Leu50Ala"], True),
    ("mutations", "L50A", "L50G", False),
    ("pdb_id", "1UBQ", "PDB 1ubq", True),
    ("pdb_id", "1UBQ", "3PTB", False),
    ("neutralized", True, "yes", True),
])
def test_two_readings_of_one_field_are_one_value_or_not(name, a, b, same):
    assert same_value(name, a, b) is same


def test_studies_are_paired_by_what_tells_them_apart():
    first = [_study("S1", "wild type", pdb_id=_stated("1UBQ")),
             _study("S2", "L50A mutant", pdb_id=_stated("1UBQ"), mutations=_stated("L50A"))]
    second = [_study("S1", "the L50A mutant", pdb_id=_stated("1UBQ"),
                     mutations=_stated("L50A")),
              _study("S2", "wild-type ubiquitin", pdb_id=_stated("1UBQ"))]
    assert match_studies(first, second) == [(0, 1), (1, 0)]


def test_what_two_readers_agree_on_is_the_truth_and_the_rest_contested():
    a = {"key": "x", "doi": "10.9999/x", "studies": [_study(
        "S1", "wild type", temperature=_stated(300, "K"), timestep=_stated(2, "fs"),
        water_model=_stated("TIP3P"), replicas={"status": "in_si"})]}
    b = {"key": "x", "studies": [_study(
        "S1", "wild type", temperature=_stated(300, "K"), timestep=_stated(4, "fs"),
        replicas={"status": "in_si"})], "claims": []}
    a["studies"][0]["claims"] = [{"quantity": "rmsd", "value": 0.12, "unit": "nm",
                                  "quote": "0.12 nm"}]
    b["studies"][0]["claims"] = [{"quantity": "rmsd", "value": 0.120, "unit": "nm"}]
    truth = reconcile(a, b)
    study = truth["studies"][0]
    assert study["fields"]["temperature"]["value"] == 300
    assert study["fields"]["replicas"] == {"status": "in_si"}
    assert set(study["contested"]) == {"timestep", "water_model"}
    assert study["claims"][0]["value"] == 0.12
    assert "quote" not in json.dumps(truth)


def test_a_reading_is_counted_against_the_agreed_fields_only():
    truth = {"key": "x", "studies": [{
        "id": "S1/S1", "label": "wild type", "contested": ["timestep"],
        "fields": {"temperature": {"status": "stated", "value": 300},
                   "water_model": {"status": "stated", "value": "TIP3P"},
                   "pressure": {"status": "stated", "value": 1.0},
                   "system": {"status": "stated", "value": "ubiquitin"}},
        "claims": [{"quantity": "rmsd", "value": 0.12}]}]}
    reading = {"studies": [_study(
        "S1", "wild type", temperature=_stated(300), water_model=_stated("TIP4P-Ew"),
        timestep=_stated(9), thermostat=_stated("Langevin"),
        system=_stated("something else entirely"))]}
    reading["studies"][0]["claims"] = [{"quantity": "rmsd", "value": 0.12, "status": "stated"}]
    counted = score(truth, reading)
    fields = counted["fields"]
    assert (fields["correct"], fields["wrong"], fields["missed"], fields["invented"]) == (1, 1, 1, 1)
    assert "timestep" not in counted["by_field"]  # contested
    assert "system" not in counted["by_field"]  # a description
    assert counted["by_field"]["thermostat"]["invented"] == 1
    assert (counted["claims_agreed"], counted["claims_found"]) == (1, 1)


def test_every_registered_paper_has_its_truth_and_the_truth_keeps_no_words():
    assert len(PAPERS) == 20 and len({paper["key"] for paper in PAPERS}) == 20
    for paper in PAPERS:
        truth = json.loads((TRUTH / f"{paper['key']}.json").read_text(encoding="utf-8"))
        assert truth["studies"], paper["key"]
        for study in truth["studies"]:
            assert set(study["fields"]) <= set(BY_NAME), paper["key"]
            assert set(study["contested"]) <= set(BY_NAME), paper["key"]
            for record in study["fields"].values():
                assert "quote" not in record
    assert sorted(path.stem for path in TRUTH.glob("*.json")) == sorted(p["key"] for p in PAPERS)


def test_the_harness_reads_each_paper_and_writes_the_counts(tmp_path, monkeypatch):
    truth = json.loads((TRUTH / "fischer_2024.json").read_text(encoding="utf-8"))
    reading = {"model": "test/model", "studies": [
        {"id": study["id"], "label": study["label"], "claims": [],
         "fields": {name: dict(record) for name, record in study["fields"].items()}}
        for study in truth["studies"]]}
    asked = []

    def studies_in(source, use_kept=True):
        asked.append(source)
        return None, reading

    monkeypatch.setattr("fastmdxplora.paper.studies.studies_in", studies_in)
    out = tmp_path / "counts.json"
    assert main(["--papers", "fischer_2024", "--file", "fischer_2024=here.pdf",
                 "--out", str(out)]) == 0
    assert asked == ["here.pdf"]
    totals = json.loads(out.read_text(encoding="utf-8"))["totals"]
    assert totals["wrong"] == totals["missed"] == totals["invented"] == 0
    assert totals["correct"] > 0


def test_the_paper_not_open_is_tried_by_its_doi_first(tmp_path, monkeypatch):
    from fastmdxplora.paper import PaperRefused

    def refuse(source, *args, **kwargs):
        raise PaperRefused("not open", code="environment.paper.not_open")

    monkeypatch.setattr("fastmdxplora.paper.fetch.open_paper", refuse)
    monkeypatch.setattr("fastmdxplora.paper.studies.studies_in",
                        lambda source, use_kept=True: (None, {"studies": []}))
    out = tmp_path / "counts.json"
    assert main(["--papers", "buch_2011", "--file", "buch_2011=buch.pdf", "--out", str(out)]) == 0
    result = json.loads(out.read_text(encoding="utf-8"))["results"][0]
    assert result["refused_by_doi"] == "environment.paper.not_open"


def test_the_truth_is_made_again_from_two_readers_folders(tmp_path, capsys):
    first, second, out = tmp_path / "a", tmp_path / "b", tmp_path / "truth"
    first.mkdir()
    second.mkdir()
    one = {"key": "x", "studies": [_study("S1", "wild type", temperature=_stated(300))]}
    (first / "x.json").write_text(json.dumps(one), encoding="utf-8")
    (second / "x.json").write_text(json.dumps(one), encoding="utf-8")
    assert paper_reading.main(["--reconcile", str(first), str(second), str(out)]) == 0
    assert "1 studies matched" in capsys.readouterr().out
    assert json.loads((out / "x.json").read_text())["studies"][0]["fields"]["temperature"][
        "value"] == 300


# -- sittings: a paper not read, not found, and the claims judged at the end --
def _reading_of(key):
    truth = json.loads((TRUTH / f"{key}.json").read_text(encoding="utf-8"))
    return {"model": "test/model", "studies": [
        {"id": study["id"], "label": study["label"], "claims": [],
         "fields": {name: dict(record) for name, record in study["fields"].items()}}
        for study in truth["studies"]]}


_KEYS = [paper["key"] for paper in PAPERS]
_LIMIT = ("Anthropic refused the request (400): {\"type\":\"error\",\"error\":{\"type\":"
          "\"invalid_request_error\",\"message\":\"You have reached your specified API usage "
          "limits. You will regain access on 2026-11-01 at 00:00 UTC.\"}}")


def _reader(monkeypatch, refuse=None):
    """The harness's reading, faked: each paper read as its truth, but those
    ``refuse`` names, which raise what it gives them. Buch 2011's DOI is
    refused as not open."""
    from fastmdxplora.paper import PaperRefused
    from fastmdxplora.refusals import StudyError

    asked = []
    route = {paper["route"]: paper["key"] for paper in PAPERS} | {"buch.pdf": "buch_2011"}

    def studies_in(source, use_kept=True):
        key = route.get(source, source)
        asked.append(key)
        if refuse and key in refuse:
            code, said = refuse[key]
            if isinstance(code, BaseException):
                raise code
            raise (PaperRefused if code.startswith("environment.paper.") else StudyError)(
                said, code=code)
        return None, _reading_of(key)

    def not_open(source, *args, **kwargs):
        raise PaperRefused("not open", code="environment.paper.not_open")

    monkeypatch.setattr("fastmdxplora.paper.studies.studies_in", studies_in)
    monkeypatch.setattr("fastmdxplora.paper.fetch.open_paper", not_open)
    return asked


def test_a_refusal_is_printed_with_its_reason(tmp_path, monkeypatch, capsys):
    _reader(monkeypatch, {"ilter_2025": ("environment.service.unusable_response", _LIMIT)})
    main(["--papers", "ilter_2025", "--out", str(tmp_path / "s.json")])
    printed = capsys.readouterr().out
    assert "ilter_2025: not read (environment.service.unusable_response)" in printed
    assert "You have reached your specified API usage limits" in printed


def test_the_file_is_written_after_each_paper(tmp_path, monkeypatch):
    asked = _reader(monkeypatch)
    out = tmp_path / "s.json"
    from fastmdxplora.paper import studies as studies_module

    reading = studies_module.studies_in

    def stopped_by_hand(source, use_kept=True):
        if len(asked) == 2:
            raise KeyboardInterrupt
        return reading(source, use_kept)

    monkeypatch.setattr("fastmdxplora.paper.studies.studies_in", stopped_by_hand)
    with pytest.raises(KeyboardInterrupt):
        main(["--out", str(out), "--file", "buch_2011=buch.pdf"])
    kept = json.loads(out.read_text(encoding="utf-8"))
    assert [r["key"] for r in kept["results"]] == _KEYS[:2]
    assert kept["claims"]["values_used_right"]["met"] is None  # not judged: 18 not read


def test_the_account_s_limit_stops_the_sitting_and_the_rest_are_not_read(tmp_path, monkeypatch,
                                                                         capsys):
    asked = _reader(monkeypatch, {"ilter_2025": ("environment.service.unusable_response",
                                                 _LIMIT)})
    out = tmp_path / "s.json"
    assert main(["--out", str(out), "--file", "buch_2011=buch.pdf"]) == 0
    assert asked[-1] == "ilter_2025" and len(asked) == _KEYS.index("ilter_2025") + 1
    kept = json.loads(out.read_text(encoding="utf-8"))
    rest = _KEYS[_KEYS.index("ilter_2025"):]
    assert kept["claims"]["not_read"] == rest
    assert kept["claims"]["studies_found"] == {"value": None, "met": None}
    printed = capsys.readouterr().out
    assert "the sitting stops here" in printed and f"--earlier {out}" in printed


@pytest.mark.parametrize("code, said", [
    ("environment.service.unusable_response",
     "Anthropic refused the request (529): {\"type\":\"overloaded_error\"}"),
    ("environment.service.unreachable", "Could not reach https://api.anthropic.com"),
])
def test_a_refusal_for_this_paper_alone_does_not_stop_the_sitting(tmp_path, monkeypatch,
                                                                  code, said):
    asked = _reader(monkeypatch, {"ilter_2025": (code, said)})
    out = tmp_path / "s.json"
    main(["--out", str(out), "--file", "buch_2011=buch.pdf"])
    assert asked == _KEYS
    kept = json.loads(out.read_text(encoding="utf-8"))
    assert kept["claims"]["not_read"] == ["ilter_2025"]


def test_a_paper_fastmdxplora_could_not_read_is_not_found_and_the_claims_judged(
        tmp_path, monkeypatch, capsys):
    _reader(monkeypatch, {"love_2023": (
        "environment.paper.answer_cut_off", "FastMDXplora could not read the paper: ...")})
    out = tmp_path / "s.json"
    main(["--out", str(out), "--file", "buch_2011=buch.pdf"])
    claims = json.loads(out.read_text(encoding="utf-8"))["claims"]
    assert claims["not_found"] == ["love_2023"] and claims["not_read"] == []
    assert claims["studies_found"]["value"] == pytest.approx((245 - 48) / 245, abs=1e-4)
    assert claims["values_used_right"]["met"] is True
    assert claims["not_open_refused"]["met"] is True
    assert "love_2023: not found (environment.paper.answer_cut_off)" in capsys.readouterr().out


def test_a_later_sitting_reads_only_the_papers_not_read(tmp_path, monkeypatch):
    first, second = tmp_path / "first.json", tmp_path / "second.json"
    _reader(monkeypatch, {"ilter_2025": ("environment.service.unusable_response", _LIMIT)})
    main(["--out", str(first), "--file", "buch_2011=buch.pdf"])
    asked = _reader(monkeypatch)
    main(["--earlier", str(first), "--out", str(second)])
    rest = _KEYS[_KEYS.index("ilter_2025"):]
    assert asked == rest  # the papers read before are not read again
    kept = json.loads(second.read_text(encoding="utf-8"))
    assert kept["sittings"] == 2 and kept["read_again"] == rest
    assert kept["claims"]["not_read"] == []
    assert kept["claims"]["studies_found"] == {"value": 1.0, "met": True}
    assert {r["key"]: r["sitting"] for r in kept["results"]}["seekr2_2022"] == 1


def test_a_paper_not_found_stays_not_found_whatever_a_later_file_holds(tmp_path, monkeypatch):
    # check 5's files: FastMDXplora's own failure then carried a service's
    # code; pantsar_2018 was read whole in a third sitting, against the rule.
    second = {"sittings": 1, "results": [
        {"key": "pantsar_2018", "code": "environment.service.unusable_response",
         "refused": "FastMDXplora could not read the paper: its reading came back cut off",
         "studies_truth": 14},
        {"key": "siwy_2017", "code": "environment.service.unusable_response",
         "refused": _LIMIT, "studies_truth": 4}]}
    third = {"sittings": 2, "results": [
        {"key": "pantsar_2018", "studies_truth": 14, "studies_read": 14,
         "studies_matched": 14, "fields": {"correct": 320}}]}
    paths = []
    for number, body in enumerate((second, third)):
        path = tmp_path / f"sitting{number}.json"
        path.write_text(json.dumps(body), encoding="utf-8")
        paths += ["--earlier", str(path)]
    asked = _reader(monkeypatch)
    out = tmp_path / "s.json"
    main(paths + ["--papers", "pantsar_2018,siwy_2017", "--out", str(out)])
    assert asked == ["siwy_2017"]
    kept = json.loads(out.read_text(encoding="utf-8"))
    assert kept["sittings"] == 3
    assert "pantsar_2018" in kept["claims"]["not_found"]


def test_buch_2011_is_not_read_without_its_pdf(tmp_path, monkeypatch, capsys):
    asked = _reader(monkeypatch)
    assert main(["--out", str(tmp_path / "s.json")]) == 2
    assert asked == [] and "--file buch_2011=PATH" in capsys.readouterr().err
    first = tmp_path / "first.json"
    assert main(["--out", str(first), "--file", "buch_2011=buch.pdf",
                 "--papers", "buch_2011"]) == 0
    assert main(["--earlier", str(first), "--papers", "buch_2011",
                 "--out", str(tmp_path / "second.json")]) == 0  # read before: no file needed


@pytest.mark.parametrize("refusal", [
    ("environment.paper.not_open", "Europe PMC says it is not open access"),
    ("environment.path.not_found", "No file here.pdf"),
    ("environment.paper.unreadable", "here.pdf holds no text (a scanned PDF)"),
    (ConnectionResetError(54, "Connection reset by peer"), ""),
    ("environment.service.unusable_response",
     "Gemini refused the request (429): Quota exceeded for requests per minute"),
])
def test_a_paper_not_reached_is_not_read_and_the_sitting_goes_on(tmp_path, monkeypatch,
                                                                 refusal):
    asked = _reader(monkeypatch, {"ilter_2025": refusal})
    out = tmp_path / "s.json"
    assert main(["--out", str(out), "--file", "buch_2011=buch.pdf"]) == 0
    assert asked == _KEYS
    claims = json.loads(out.read_text(encoding="utf-8"))["claims"]
    assert claims["not_read"] == ["ilter_2025"] and claims["not_found"] == []


def test_found_no_study_is_a_reading_and_counts_as_not_found(tmp_path, monkeypatch):
    _reader(monkeypatch, {"atlas_2024": (
        "environment.paper.unreadable", "FastMDXplora found no MD study in ATLAS. If it ...")})
    out = tmp_path / "s.json"
    main(["--out", str(out), "--file", "buch_2011=buch.pdf"])
    assert json.loads(out.read_text(encoding="utf-8"))["claims"]["not_found"] == ["atlas_2024"]


def test_a_paper_not_read_and_not_tried_now_is_kept_in_the_next_file(tmp_path, monkeypatch):
    first, second, third = (tmp_path / f"{name}.json" for name in ("a", "b", "c"))
    _reader(monkeypatch, {"ilter_2025": ("environment.service.unreachable", "no answer"),
                          "siwy_2017": ("environment.service.unreachable", "no answer")})
    main(["--out", str(first), "--file", "buch_2011=buch.pdf"])
    _reader(monkeypatch)
    main(["--earlier", str(first), "--papers", "ilter_2025", "--out", str(second)])
    kept = json.loads(second.read_text(encoding="utf-8"))
    assert kept["claims"]["not_read"] == ["siwy_2017"]
    assert any(r["key"] == "siwy_2017" and r["status"] == "not_read" for r in kept["results"])
    asked = _reader(monkeypatch)
    main(["--earlier", str(second), "--out", str(third)])
    assert asked == ["siwy_2017"]
    last = json.loads(third.read_text(encoding="utf-8"))
    assert last["sittings"] == 3 and last["read_again"] == ["siwy_2017"]
    assert last["claims"]["not_read"] == []


def test_an_earlier_file_that_is_no_sitting_s_is_refused(tmp_path, capsys):
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"sittings": "two", "results": []}), encoding="utf-8")
    assert main(["--earlier", str(bad), "--out", str(tmp_path / "s.json")]) == 2
    assert "not a sitting's file" in capsys.readouterr().err
