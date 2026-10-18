"""A remedy is not priced from a record a held read did not read.

What remains of a stopped run, and how much longer a thin window must run,
are read from the run's own records (`resume.production_done_ns`,
`resume.extension_of`, its `resolved_config.yml`). In a read held to the
study (`record_access`), a record linked out of it is refused; the remedy
caught the refusal with every other error and priced the run as though
nothing had been done (its whole plan), the planned length read through the
very link refused. A remedy now says its price, or its length, is not given
and which record was not read. Outside a held read nothing changes.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fastmdxplora.remedies import remedies_of
from tests.test_every_refusal_says_its_fix_and_price import (
    STOPPED,
    _batch,
    _config,
    _manifest,
    _pmf,
    _run,
    _speed,
    _stopped_at,
    _windows,
)

#: What a remedy says of a record a held read refused.
NOT_READ = ("was not read: it leads out of the study, or to a file of it that is not one "
            "of its results.")


def _linked_out(root: Path, tmp_path: Path, relative: str) -> None:
    path = root / relative
    target = tmp_path / "outside" / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    path.rename(target)
    try:
        path.symlink_to(target)
    except OSError:
        pytest.skip("this system cannot create symbolic links")


def _stopped(tmp_path: Path) -> Path:
    """A study stopped in production at 0.2 of its 0.5 ns."""
    root = tmp_path / "study"
    _config(root)
    _speed(root)
    _stopped_at(root, 100_000)
    _manifest(root, "simulation", STOPPED)
    return root


def _thin(tmp_path: Path) -> Path:
    """An umbrella study whose window 2 recorded 150 of the 200 values."""
    root = tmp_path / "umbrella"
    planned = _windows(root)
    _batch(root, [_run(root, p["run_id"], "ok") for p in planned], planned)
    _pmf(root, refused="1 of 4 windows recorded fewer than 200 values.",
         thin=[{"window": 2, "samples": 150}], plan={"minimum_samples": 200})
    return root


def _held(root: Path) -> dict:
    """The study's fixes, read as a confined route reads (`server._confined`)."""
    from fastmdxplora.gui.fixes_view import fixes_payload
    from fastmdxplora.gui.server import _confined

    return _confined(root, fixes_payload)


@pytest.mark.parametrize("record", ["resolved_config.yml", "simulation/checkpoint.chk.json"])
def test_a_stopped_run_s_record_not_read_gives_no_price(tmp_path, record):
    root = _stopped(tmp_path)
    (fix,) = _held(root)["fixes"]
    assert fix["price"]["production_ns"] == pytest.approx(0.3)

    _linked_out(root, tmp_path, record)
    answer = _held(root)
    (fix,) = answer["fixes"]
    assert record in answer["unavailable"]
    # Not its whole plan, as a run that did nothing.
    assert fix["price"] is None and fix["price_said"] == ""
    assert fix["fix"].endswith(f" The price is not given, since {record} {NOT_READ}")
    assert fix["command"] == f"fastmdx resume {root.resolve()}"


def test_a_thin_window_s_record_not_read_gives_no_length(tmp_path):
    root = _thin(tmp_path)
    (fix,) = _held(root)["fixes"]
    assert "--simulate-duration-ns 1.4" in fix["command"]

    _linked_out(root, tmp_path, "runs/w02/resolved_config.yml")
    answer = _held(root)
    (fix,) = answer["fixes"]
    assert "runs/w02/resolved_config.yml" in answer["unavailable"]
    # No length read through the link, and no price.
    assert "--simulate-duration-ns" not in fix["command"]
    assert fix["command"].endswith("--rerun-window 2")
    assert fix["price"] is None
    assert fix["fix"] == ("Run it again for longer; how much longer is not given, since "
                          f"runs/w02/resolved_config.yml {NOT_READ}")


@pytest.mark.parametrize("study,record", [
    (_stopped, "resolved_config.yml"), (_stopped, "simulation/checkpoint.chk.json"),
    (_thin, "runs/w02/resolved_config.yml")])
def test_outside_a_held_read_the_remedy_is_as_before(tmp_path, study, record):
    """The command line follows a link as it always did: the same remedy,
    to the word, and the same console lines."""
    from fastmdxplora.batch.explorer import _what_would_fix_it

    root = study(tmp_path)
    before = [remedy.as_text() for remedy in remedies_of(root)]
    lines = _what_would_fix_it(root)
    _linked_out(root, tmp_path, record)
    assert [remedy.as_text() for remedy in remedies_of(root)] == before
    assert _what_would_fix_it(root) == lines
    assert "Costs " in before[0]


def test_any_other_error_is_priced_as_before(tmp_path, monkeypatch):
    """A record that cannot be made sense of is not a refusal: what remains
    is the whole plan, as it was."""
    from fastmdxplora.simulation import resume

    root = _stopped(tmp_path)

    def broken(*_args, **_kwargs):
        raise RuntimeError("a record that makes no sense")

    monkeypatch.setattr(resume, "extension_of", broken)
    monkeypatch.setattr(resume, "production_done_ns", broken)
    (fix,) = _held(root)["fixes"]
    assert (fix["price"]["production_ns"], fix["price"]["equilibration_ns"]) == (
        pytest.approx(0.5), pytest.approx(0.2))
    assert json.dumps(fix).count("not given") == 0


def test_a_run_again_from_the_top_is_not_priced_from_a_plan_not_read(tmp_path):
    """Run again from the top, a study costs its whole plan: with its plan
    not read, no price is given, not one read through the link."""
    from fastmdxplora.gui.server import _confined
    from fastmdxplora.remedies import remedy_for
    from tests.test_every_refusal_says_its_fix_and_price import PROTONATION

    root = tmp_path / "study"
    _config(root, duration_ns=10.0)
    _manifest(root, "setup", PROTONATION)

    def fix(inside: Path) -> dict:
        return remedy_for(PROTONATION, study=inside).as_record()

    assert _confined(root, fix)["price"]["production_ns"] == 10.0
    _linked_out(root, tmp_path, "resolved_config.yml")
    said = _confined(root, fix)
    assert said["price"] is None and said["unavailable"] == ["resolved_config.yml"]
    assert said["fix"].endswith(f" The price is not given, since resolved_config.yml {NOT_READ}")
    # Outside a held read, as before.
    assert remedy_for(PROTONATION, study=root).price.production_ns == 10.0


@pytest.mark.parametrize("record", ["resolved_config.yml", "simulation/checkpoint.chk.json"])
def test_the_fixes_route_prices_no_remedy_from_a_record_not_read(tmp_path, record):
    """The Overview's card asks `/api/fixes`: a stopped run's config, or
    its checkpoint's record, linked out of the study prices nothing there,
    where the page says the record was not read."""
    from tests.test_the_study_pages_read_only_the_study import _ask

    root = _stopped(tmp_path)
    before = json.loads(_ask(root, tmp_path, "/api/fixes", public=False)[1])
    assert before["fixes"][0]["price"]["production_ns"] == pytest.approx(0.3)
    assert "unavailable" not in before

    _linked_out(root, tmp_path, record)
    status, body = _ask(root, tmp_path, "/api/fixes", public=False)
    assert status == 200
    answer = json.loads(body)
    (fix,) = answer["fixes"]
    assert fix["price"] is None and fix["price_said"] == ""
    assert fix["fix"].endswith(f" The price is not given, since {record} {NOT_READ}")
    assert answer["unavailable"] == [record]
    assert answer["study"] == str(root)
