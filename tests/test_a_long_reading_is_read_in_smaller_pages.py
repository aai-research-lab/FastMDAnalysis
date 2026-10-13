"""A long paper's reading is read whole, in smaller pages where it must be.

An answer that comes back cut off (at the length FastMDXplora reads at a
time) is read again smaller: half as many studies or results, down to one,
and a protocol's settings in two parts, down to one setting. Still cut off
at one, the reading is refused with a code of its own; an answer whole but
not in the form asked for is asked for once more, then refused with
another. Neither is the code a provider's refusal carries.
"""

from __future__ import annotations

import json
import re

import pytest

from fastmdxplora.paper import PaperRefused
from fastmdxplora.paper.extract import (
    CLAIMS_PER_PAGE,
    CUT_OFF,
    MALFORMED,
    STUDIES_PER_PAGE,
    _cut_off,
    _json_from,
    _studies_prompt,
    read_studies,
)
from fastmdxplora.paper.fields import PROTOCOL_FIELDS
from fastmdxplora.paper.text import Part, PaperText


def _paper() -> PaperText:
    return PaperText(parts=[Part("Methods", "Each system was simulated for 100 ns at 300 K.",
                                 "paper")])


def _cut(answer: dict) -> str:
    whole = json.dumps(answer)
    return whole[: len(whole) * 2 // 3]


class _Long:
    """Answers as a provider does whose answers stop at a fixed length: a
    page of more than ``fits`` studies or results, or more than ``fits``
    settings, comes back cut off."""

    def __init__(self, studies: int = 12, claims: int = 0, fits: int = 2,
                 settings_fit: int = len(PROTOCOL_FIELDS)) -> None:
        self.studies = [{"id": f"S{n}", "label": f"system {n}", "protocol": "P1",
                         "fields": {"temperature": {"value": 300, "unit": "K",
                                                    "quote": "simulated for 100 ns at 300 K"}}}
                        for n in range(1, studies + 1)]
        self.claims = [{"study": "S1", "quantity": "rmsd", "value": n / 10, "unit": "nm",
                        "quote": f"RMSD {n}"} for n in range(1, claims + 1)]
        self.fits, self.settings_fit = fits, settings_fit
        self.sizes: dict[str, list[int]] = {"studies": [], "claims": [], "settings": []}
        self.settings_asked: list[tuple[str, ...]] = []

    def __call__(self, prompt: str) -> str:
        if "list the molecular dynamics (MD) studies" in prompt:
            most = int(re.search(r"list at most\s+(\d+) per answer", prompt).group(1))
            self.sizes["studies"].append(most)
            listed = re.search(r"The studies already listed are: ([^.]*)\.", prompt)
            done = len(listed.group(1).split(", ")) if listed else 0
            page = self.studies[done:done + most]
            answer = {"title": "A long paper", "studies": page,
                      "protocols": [{"id": "P1", "label": "every study"}],
                      "more": done + most < len(self.studies)}
            return _cut(answer) if most > self.fits else json.dumps(answer)
        if "settings of one of its MD protocols" in prompt:
            names = tuple(re.findall(r'^- "([a-z_0-9]+)":', prompt, re.M))
            self.sizes["settings"].append(len(names))
            self.settings_asked.append(names)
            answer = {"fields": {name: {"value": 300, "unit": "K", "quote": "simulated for 100 ns at 300 K"}
                                 for name in names if name == "temperature"}}
            return _cut(answer) if len(names) > self.settings_fit else json.dumps(answer)
        if "numerical results" in prompt:
            most = int(re.search(r"List at most (\d+) per answer", prompt).group(1))
            self.sizes["claims"].append(most)
            listed = re.search(r"(\d+) results are already listed", prompt)
            done = int(listed.group(1)) if listed else 0
            answer = {"claims": self.claims[done:done + most],
                      "more": done + most < len(self.claims)}
            return _cut(answer) if most > self.fits else json.dumps(answer)
        return "{}"


@pytest.fixture(autouse=True)
def _cache(tmp_path, monkeypatch):
    monkeypatch.setenv("FASTMDXPLORA_CACHE_DIR", str(tmp_path))


# -- telling a cut answer from one in another form ---------------------------
@pytest.mark.parametrize("reply, cut", [
    ('{"studies": [{"id": "S1", "label": "wild', True),
    ('{"studies": [{"id": "S1"}, {"id": "S2"', True),
    ('Here it is:\n```json\n{"fields": {"temperature": {"value": 30', True),
    ('{"studies": [], "more": false}', False),
    ('Here it is:\n```json\n{"studies": [], "more": false}\n```\nThat is all.', False),
    ('{"label": "a brace } and a quote \\" inside"}', False),
    ("I could not find any MD studies in this paper.", False),
    ('[{"id": "S1"}, {"id": "S2"', True),
    ('Note: I wrote { as a placeholder, and found no study.', False),
    ('{"studies": []}  (an unmatched { in a note)', False),
    ('{ "studies": [{"id": "S1", "label": "wild', True),
])
def test_an_answer_ending_inside_its_json_is_cut_off(reply, cut):
    assert _cut_off(reply) is cut


def test_a_cut_answer_and_one_in_another_form_are_refused_apart():
    with pytest.raises(PaperRefused) as cut:
        _json_from('{"studies": [{"id": "S1", "lab', "studies")
    with pytest.raises(PaperRefused) as malformed:
        _json_from('{"studies": "S1 and S2"}', "studies")
    assert cut.value.code == CUT_OFF == "environment.paper.answer_cut_off"
    assert malformed.value.code == MALFORMED == "environment.paper.answer_malformed"
    assert "cut off" in str(cut.value) and "cut off" not in str(malformed.value)
    for refused in (cut.value, malformed.value):
        assert refused.code != "environment.service.unusable_response"


def test_both_codes_are_registered_and_documented():
    from pathlib import Path

    from fastmdxplora.refusals import known

    page = (Path(__file__).resolve().parents[1] / "docs" / "refusals.md").read_text(
        encoding="utf-8")
    for code in (CUT_OFF, MALFORMED):
        assert known(code)
        assert f"`{code}`" in page


# -- smaller pages -----------------------------------------------------------
def test_the_first_page_is_asked_for_as_before():
    # A reading kept from before was asked with these words, so it is kept.
    assert "list at most\n5 per answer" in _studies_prompt("the paper", [])
    assert STUDIES_PER_PAGE == 5 and CLAIMS_PER_PAGE == 12


def test_studies_cut_off_are_read_in_smaller_pages_and_all_are_kept():
    model = _Long(studies=12, fits=2)
    said: list[str] = []
    reading = read_studies(_paper(), model, model="test/long", said=said.append)
    assert [s["id"] for s in reading["studies"]] == [f"S{n}" for n in range(1, 13)]
    assert model.sizes["studies"][:2] == [5, 2]
    assert set(model.sizes["studies"][1:]) == {2}  # the rest at the size that fits
    assert any("reading the studies 2 at a time" in line for line in said)


def test_a_page_cut_off_even_at_two_is_read_one_study_at_a_time():
    model = _Long(studies=3, fits=1)
    said: list[str] = []
    reading = read_studies(_paper(), model, model="test/long", said=said.append)
    assert len(reading["studies"]) == 3
    assert model.sizes["studies"] == [5, 2, 1, 1, 1]
    assert any("one at a time" in line for line in said)


def test_results_cut_off_are_read_in_smaller_pages():
    model = _Long(studies=1, claims=20, fits=5)
    reading = read_studies(_paper(), model, model="test/long")
    assert sum(len(study["claims"]) for study in reading["studies"]) == 20
    assert model.sizes["claims"][0] == 12 and model.sizes["claims"][1] == 6
    assert max(model.sizes["claims"][2:]) <= 3


def test_a_protocol_s_settings_cut_off_are_read_in_parts_down_to_what_fits():
    model = _Long(studies=1, settings_fit=7)
    said: list[str] = []
    reading = read_studies(_paper(), model, model="test/long", said=said.append)
    leaves = [names for names in model.settings_asked if len(names) <= 7]
    asked = [name for names in leaves for name in names]
    assert sorted(asked) == sorted(PROTOCOL_FIELDS)  # each setting asked once, none left out
    assert model.sizes["settings"][0] == len(PROTOCOL_FIELDS)
    assert reading["studies"][0]["fields"]["temperature"]["status"] == "stated"
    assert any("two parts" in line for line in said)


def test_a_study_cut_off_even_alone_is_refused_with_its_own_code(tmp_path):
    model = _Long(studies=3, fits=0)
    with pytest.raises(PaperRefused) as refused:
        read_studies(_paper(), model, model="test/long")
    assert refused.value.code == CUT_OFF
    assert model.sizes["studies"] == [5, 2, 1]
    assert not list(tmp_path.rglob("reading-*.json"))  # nothing kept


def test_a_setting_cut_off_even_alone_is_refused_with_its_own_code():
    model = _Long(studies=1, settings_fit=0)
    with pytest.raises(PaperRefused) as refused:
        read_studies(_paper(), model, model="test/long")
    assert refused.value.code == CUT_OFF
    assert min(model.sizes["settings"]) == 1


def test_an_answer_not_in_the_form_is_asked_for_once_more_then_refused():
    tries: list[str] = []

    def prose(prompt: str) -> str:
        tries.append(prompt)
        return "The paper reports one study of ubiquitin."

    with pytest.raises(PaperRefused) as refused:
        read_studies(_paper(), prose, model="test/prose")
    assert refused.value.code == MALFORMED
    assert len(tries) == 2 and tries[1].endswith("Answer with the JSON object only, whole.")


def test_an_answer_in_prose_once_is_read_when_asked_again():
    model = _Long(studies=2, fits=5)
    first = {"done": False}

    def once_in_prose(prompt: str) -> str:
        if "list the molecular dynamics (MD) studies" in prompt and not first["done"]:
            first["done"] = True
            return "Sure: there are two studies."
        return model(prompt)

    reading = read_studies(_paper(), once_in_prose, model="test/once")
    assert len(reading["studies"]) == 2


def test_listing_stops_after_its_pages_counted_at_their_first_size():
    class Endless(_Long):
        def __call__(self, prompt):
            if "list the molecular dynamics (MD) studies" in prompt:
                most = int(re.search(r"list at most\s+(\d+) per answer", prompt).group(1))
                self.sizes["studies"].append(most)
                if most > 1:
                    return '{"studies": [{"id": "S'
                listed = re.search(r"The studies already listed are: ([^.]*)\.", prompt)
                done = len(listed.group(1).split(", ")) if listed else 0
                return json.dumps({"studies": [{"id": f"S{done + 1}", "label": "x"}],
                                   "protocols": [], "more": True})
            return super().__call__(prompt)

    model = Endless(studies=0)
    said: list[str] = []
    reading = read_studies(_paper(), model, model="test/endless", said=said.append)
    # Twelve pages of five: sixty studies at most, one at a time here.
    assert len(reading["studies"]) == 60
    assert any(line.startswith("Stopped after 60 rounds") for line in said)
