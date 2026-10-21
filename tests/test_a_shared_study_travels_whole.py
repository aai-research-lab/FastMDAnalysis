"""A shared study reaches Zenodo, and comes back from it, whole.

A study of a paper is gigabytes, and a connection is not always there for
as long: Zenodo limits how often it is asked (429, saying how long to wait)
and is at times briefly unable to answer (502, 503, 504). A download cut
short goes on from where it stopped, an answer that says to ask again is
asked again after a wait, and the archive is checked against Zenodo's MD5
as it lands. A request that makes something is never asked twice where
the first may have made it.
"""

from __future__ import annotations

import hashlib
import http.server
import json
import os
import threading

import pytest

from fastmdxplora.sharing import ShareRefused, zenodo
from fastmdxplora.sharing.opening import _download, open_shared
from tests.test_a_shared_study_goes_to_zenodo import _share_to_zenodo, _Zenodo
from tests.test_a_study_is_shared import _study


class _Unsteady(_Zenodo):
    """Zenodo as above, made to fail as asked: ``faults`` maps a request
    (``"GET file"``, ``"PUT file"``, ``"POST"``, ``"GET record"``) to the
    answers it gives first, in order: an HTTP status, ``"cut"`` (the file's
    first bytes, then the connection closed), ``"whole"`` (the whole file
    with 200, whatever range was asked), ``"elsewhere"`` (206 with the file
    from byte 100, whatever range was asked), ``"stall"`` (the file's first
    bytes, then nothing for 3 s), ``"drop"`` (the connection closed with no
    answer) or ``"made"`` (a draft made, then a 504 in place of the
    answer), ``"made late"`` (as ``"made"``, the draft listed only from the
    third look) or ``"504 cut"`` (a 504 whose own words are cut short). It
    lists the person's drafts as Zenodo does, and its drafts' files carry
    their names (``files``, set by a test)."""

    faults: dict = {}
    ranges: list = []
    payload: bytes | None = None
    drafts: list = []
    files: list | None = None
    late = 2

    def _deposition(self, number):
        made = super()._deposition(number)
        if type(self).files is not None:
            made["files"] = type(self).files
        return made

    def _fault(self, which):
        queue = type(self).faults.get(which) or []
        return queue.pop(0) if queue else None

    def _list(self):
        fault = self._fault("GET list")
        type(self).calls.append(("GET", self.path, self.headers.get("Authorization")))
        if isinstance(fault, int):
            return self._status(fault)
        shown = []
        for draft in type(self).drafts:
            if draft["hidden"]:
                draft["hidden"] -= 1
                continue
            shown.append({"id": draft["id"], "metadata": {"title": draft["title"]}})
        return self._json(shown)

    def _drop(self):
        self.close_connection = True
        self.connection.shutdown(2)

    def _status(self, status):
        if status == "504 cut":
            self.send_response(504)
            self.send_header("Content-Length", "400")
            self.end_headers()
            self.wfile.write(b'{"message": ')
            self.wfile.flush()
            return self._drop()
        self.send_response(status)
        if status == 429:
            self.send_header("Retry-After", "7")
        self.send_header("Content-Length", "0")
        self.end_headers()

    def do_POST(self):  # noqa: N802
        fault = self._fault("POST")
        if self.path == "/api/deposit/depositions" and fault in (None, "made", "made late"):
            body = self._body()
            title = (json.loads(body or b"{}").get("metadata") or {}).get("title")
            self.rfile = __import__("io").BytesIO(body)
            self.headers.replace_header("Content-Length", str(len(body))) \
                if "Content-Length" in self.headers else None
            # Made at once, or seen among the drafts only from a later look.
            type(self).drafts.append({"id": 7 + len(type(self).drafts), "title": title,
                                      "hidden": type(self).late if fault == "made late"
                                      else 0})
        if fault in ("made", "made late"):
            type(self).calls.append(("POST", self.path, self.headers.get("Authorization")))
            return self._status(504)
        if fault == "504 cut":
            type(self).calls.append(("POST", self.path, self.headers.get("Authorization")))
            self._body()
            return self._status("504 cut")
        if fault == "drop":
            type(self).calls.append(("POST", self.path, self.headers.get("Authorization")))
            self._body()
            return self._drop()
        if isinstance(fault, int):
            type(self).calls.append(("POST", self.path, self.headers.get("Authorization")))
            self._body()
            return self._status(fault)
        return super().do_POST()

    def do_DELETE(self):  # noqa: N802
        fault = self._fault("DELETE")
        if isinstance(fault, int):
            type(self).calls.append(("DELETE", self.path, self.headers.get("Authorization")))
            return self._status(fault)
        return super().do_DELETE()

    def do_PUT(self):  # noqa: N802
        if self.path.startswith("/api/deposit/depositions/"):
            fault = self._fault("PUT metadata")
            if isinstance(fault, int):
                type(self).calls.append(("PUT", self.path, self.headers.get("Authorization")))
                self._body()
                return self._status(fault)
        if self.path.startswith("/files/"):
            fault = self._fault("PUT file")
            if isinstance(fault, int):
                type(self).calls.append(("PUT", self.path, self.headers.get("Authorization")))
                self._body()
                return self._status(fault)
        return super().do_PUT()

    def do_GET(self):  # noqa: N802
        archive = type(self).payload or type(self).stored.get("archive", b"")
        name = type(self).stored.get("name", "shared.zip")
        if self.path.startswith("/api/deposit/depositions?"):
            return self._list()
        if self.path.startswith("/api/deposit/depositions/"):
            fault = self._fault("GET deposition")
            if isinstance(fault, int):
                return self._status(fault)
        if self.path == "/api/records/7":
            fault = self._fault("GET record")
            if fault == "drop":
                type(self).calls.append(("GET", self.path, None))
                return self._drop()
            if isinstance(fault, int):
                return self._status(fault)
            if type(self).payload is not None:
                return self._json({"id": 7, "doi": "10.5072/zenodo.7",
                                   "files": [{"key": name, "size": len(archive),
                                              "checksum": "md5:"
                                              + hashlib.md5(archive).hexdigest()}]})
        if self.path.startswith(f"/records/7/files/{name}"):
            asked = self.headers.get("Range")
            type(self).ranges.append(asked)
            fault = self._fault("GET file")
            if fault == 416 and not asked:
                fault = None
            if isinstance(fault, int):
                return self._status(fault)
            start = int(asked.split("=")[1].split("-")[0]) if asked and fault != "whole" else 0
            if fault == "elsewhere":
                start = 100
            body = archive[start:]
            self.send_response(206 if start else 200)
            if start:
                self.send_header("Content-Range", f"bytes {start}-{len(archive) - 1}/"
                                                  f"{len(archive)}")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            if fault == "stall":
                import time

                self.wfile.write(body[:len(body) // 3])
                self.wfile.flush()
                time.sleep(3)
                self.close_connection = True
                return
            if fault == "cut":
                self.wfile.write(body[:len(body) // 3])
                self.wfile.flush()
                self.close_connection = True
                self.connection.shutdown(2)
                return
            self.wfile.write(body)
            return
        return super().do_GET()


@pytest.fixture
def unsteady(monkeypatch):
    _Unsteady.calls, _Unsteady.stored, _Unsteady.submitted = [], {}, False
    _Unsteady.wrong_md5, _Unsteady.faults, _Unsteady.ranges = False, {}, []
    _Unsteady.payload, _Unsteady.drafts, _Unsteady.files = None, [], None
    _Unsteady.late = 2
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _Unsteady)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    monkeypatch.setenv("FASTMDXPLORA_ZENODO_URL", f"http://127.0.0.1:{server.server_port}")
    monkeypatch.setenv("ZENODO_SANDBOX_TOKEN", "sandbox-token")
    waits: list[float] = []
    monkeypatch.setattr(zenodo, "_sleep", waits.append)
    _Unsteady.waits = waits
    yield _Unsteady
    server.shutdown()


def _shared(tmp_path, unsteady) -> None:
    assert _share_to_zenodo(_study(tmp_path / "study"), tmp_path / "shared.zip") == 0
    unsteady.calls.clear()


def test_a_download_cut_short_goes_on_from_where_it_stopped(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET file": ["cut", "cut"]}
    told: list[str] = []
    opened = open_shared("10.5072/zenodo.7", workspace=tmp_path / "w", said=told.append)
    assert (opened / "shared_from.json").is_file()
    assert unsteady.ranges[0] is None
    assert unsteady.ranges[1].startswith("bytes=") and unsteady.ranges[1] != "bytes=0-"
    assert len(unsteady.ranges) == 3
    assert any("stopped at" in line and "trying again" in line for line in told)


def test_a_server_that_sends_the_whole_file_again_is_read_from_its_start(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET file": ["cut", "whole"]}
    told: list[str] = []
    opened = open_shared("10.5072/zenodo.7", workspace=tmp_path / "w", said=told.append)
    assert (opened / "shared_from.json").is_file()
    assert any("sent the whole file again, not the rest of it" in line for line in told)


def test_too_many_requests_waits_as_long_as_zenodo_asks(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET record": [429], "GET file": [503, 429]}
    open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert unsteady.waits == [7.0, 1.0, 7.0]


def test_zenodo_unable_for_longer_than_the_tries_is_said(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET file": [503] * zenodo.TRIES}
    with pytest.raises(ShareRefused) as refused:
        open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert refused.value.refusal.code == "environment.service.unreachable"
    assert f"HTTP 503, still after {zenodo.TRIES} tries" in str(refused.value)
    assert not list((tmp_path / "w").glob("*/shared_from.json"))


def test_a_resumed_download_is_still_checked_against_zenodo_s_md5(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET file": ["cut"]}
    unsteady.wrong_md5 = True
    with pytest.raises(ShareRefused) as refused:
        open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert refused.value.refusal.code == "environment.share.unverified"
    # Checked once whole: the rest was fetched before the MD5 was compared.
    assert len(unsteady.ranges) == 2 and unsteady.ranges[1].startswith("bytes=")


def test_an_upload_zenodo_could_not_take_is_sent_again(tmp_path, unsteady):
    unsteady.faults = {"PUT file": [503, 502]}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.stored["archive"] == (tmp_path / "shared.zip").read_bytes()
    assert [c[0] for c in unsteady.calls if c[1].startswith("/files/")] == ["PUT"] * 3


def test_a_draft_is_never_asked_for_twice_where_the_first_may_have_been_made(tmp_path,
                                                                              unsteady,
                                                                              capsys):
    unsteady.faults = {"POST": [503]}
    assert _share_to_zenodo(_study(tmp_path / "study"), tmp_path / "shared.zip") == 2
    assert [c for c in unsteady.calls if c[0] == "POST"] == [
        ("POST", "/api/deposit/depositions", "Bearer sandbox-token")]
    assert "HTTP 503" in capsys.readouterr().err
    # Too many requests: Zenodo made nothing, and the draft is asked for again.
    unsteady.calls.clear()
    unsteady.faults = {"POST": [429]}
    assert _share_to_zenodo(_study(tmp_path / "again"), tmp_path / "again.zip") == 0
    assert [c[1] for c in unsteady.calls if c[0] == "POST"] == [
        "/api/deposit/depositions"] * 2


def test_a_large_transfer_says_how_far_it_has_gone(tmp_path, unsteady):
    unsteady.payload = os.urandom(12 << 20)
    told: list[str] = []
    _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
              tmp_path / "got.zip", most=1 << 30,
              md5=hashlib.md5(unsteady.payload).hexdigest(), said=told.append,
              size=len(unsteady.payload))
    assert (tmp_path / "got.zip").read_bytes() == unsteady.payload
    progress = [line for line in told if line.startswith("Downloaded ") and " of " in line]
    assert 5 <= len(progress) <= 9 and progress[-1].endswith("of 13 MB.")
    told.clear()
    big = tmp_path / "big.zip"
    big.write_bytes(unsteady.payload)
    zenodo._call("PUT", f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/files/b7/big.zip",
                 upload=big, said=told.append)
    assert unsteady.stored["archive"] == unsteady.payload
    assert 5 <= len([line for line in told if line.startswith("Uploaded ")]) <= 9


def test_the_wait_is_what_zenodo_asks_within_a_minute():
    assert zenodo.wait_before(1, "7") == 7.0
    assert zenodo.wait_before(3, None) == 4.0
    assert zenodo.wait_before(2, "Wed, 21 Oct 2026 07:28:00 GMT") == 2.0
    assert zenodo.wait_before(1, "3600") == 60.0
    assert json.dumps(sorted(zenodo.AGAIN)) == "[429, 502, 503, 504]"



def test_a_part_from_anywhere_but_where_it_stopped_is_asked_for_again(tmp_path, unsteady):
    """Review of round 1: a part of the file from another place was added
    where the download had stopped."""
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET file": ["cut", "elsewhere"]}
    opened = open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert (opened / "shared_from.json").is_file()
    assert len(unsteady.ranges) == 3 and unsteady.ranges[2] == unsteady.ranges[1]


def test_a_disk_that_fills_is_said_not_taken_for_a_lost_connection(tmp_path, unsteady,
                                                                   monkeypatch):
    """Review of round 1: a piece this computer could not write was taken
    for a connection cut, skipped, and the MD5 passed over the gap."""
    import errno
    import io

    unsteady.payload = os.urandom(3 << 20)
    real = io.open
    writes = []
    failing = [errno.ENOSPC, "No space left on device"]

    class Full:
        def __init__(self, handle):
            self.handle = handle

        def write(self, piece):
            writes.append(len(piece))
            if len(writes) == 2:
                raise OSError(*failing)
            return self.handle.write(piece)

        def __getattr__(self, name):
            return getattr(self.handle, name)

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.handle.close()

    target = tmp_path / "got.zip"

    def opened(file, mode="r", *args, **kwargs):
        handle = real(file, mode, *args, **kwargs)
        return Full(handle) if str(file) == str(target) and "w" in mode else handle

    monkeypatch.setattr(io, "open", opened)
    with pytest.raises(ShareRefused) as failed:
        _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
                  target, most=1 << 30, md5=hashlib.md5(unsteady.payload).hexdigest(),
                  said=None, size=len(unsteady.payload))
    assert failed.value.refusal.code == "environment.share.too_large"
    assert f"could not be written to {tmp_path}: the disk there has no room for it" \
        in str(failed.value)
    assert unsteady.ranges == [None]
    # Another error writing is the system's own, never asked for again.
    writes.clear()
    failing[:] = [errno.EIO, "Input/output error"]
    with pytest.raises(OSError):
        _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
                  target, most=1 << 30, md5="", said=None, size=len(unsteady.payload))
    assert unsteady.ranges == [None, None]


def test_a_request_that_lost_its_connection_is_asked_again_unless_it_makes_something(
        tmp_path, unsteady):
    """Review of round 1: asking again after a connection lost, and not for a
    draft, were untested."""
    _shared(tmp_path, unsteady)
    unsteady.calls.clear()
    unsteady.faults = {"GET record": ["drop"]}
    open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert [c for c in unsteady.calls if c[1] == "/api/records/7"] == [
        ("GET", "/api/records/7", None)] * 2
    unsteady.calls.clear()
    unsteady.faults = {"POST": ["drop"]}
    assert _share_to_zenodo(_study(tmp_path / "again"), tmp_path / "again.zip") == 2
    assert [c[1] for c in unsteady.calls if c[0] == "POST"] == ["/api/deposit/depositions"]
    base = os.environ["FASTMDXPLORA_ZENODO_URL"]
    assert zenodo._call("DELETE", f"{base}/api/deposit/depositions/99", missing_ok=True) \
        is None


def test_the_wait_is_until_zenodo_s_limit_resets():
    import time

    reset = time.time() + 5
    assert 5.0 <= zenodo.wait_before(1, {"X-RateLimit-Reset": str(reset)}) <= 6.5
    assert zenodo.wait_before(2, {"X-RateLimit-Reset": "0"}) == 2.0
    assert zenodo.wait_before(1, {"Retry-After": "3", "X-RateLimit-Reset": str(reset)}) == 3.0


def test_a_draft_made_whose_answer_was_lost_is_taken_not_made_again(tmp_path, unsteady):
    """Hands-on review of round 1: Zenodo made the draft and its answer was
    lost (a 504), and sharing again made a second, the first left empty."""
    unsteady.faults = {"POST": ["made"]}
    told: list[str] = []
    from fastmdxplora.sharing.pack import share_study

    study = _study(tmp_path / "study")
    zenodo.draft(study, tmp_path / "shared.zip", site="zenodo-sandbox",
                 pack=lambda doi: share_study(study, tmp_path / "shared.zip",
                                              authors=[("A B", "")], identifier=doi),
                 authors=[("A B", "")], license_id="CC-BY-4.0", system="1ABC",
                 said=told.append)
    assert [c[1] for c in unsteady.calls if c[0] == "POST"] == ["/api/deposit/depositions"]
    assert [d["id"] for d in unsteady.drafts] == [7]
    assert unsteady.stored["archive"] == (tmp_path / "shared.zip").read_bytes()
    assert any("Zenodo made draft 7 though its answer was lost; taking that one" in line
               for line in told)


def test_a_draft_keeps_its_archive_until_the_new_one_is_there(tmp_path, unsteady, capsys):
    """Hands-on review of round 1: shared again while Zenodo would take no
    upload, the draft lost the archive it had, and the message said neither
    that nor what to do."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    capsys.readouterr()
    unsteady.calls.clear()
    unsteady.stored["archive"] = b"the archive it had"
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert not any(c[0] == "DELETE" for c in unsteady.calls)
    assert unsteady.stored["archive"] == b"the archive it had"
    said = capsys.readouterr().err
    assert "/uploads/7 is kept, with the archive it had: share again once Zenodo answers" in said
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert not any(c[0] == "POST" for c in unsteady.calls)


def test_a_stalled_download_is_asked_for_again_from_where_it_stopped(tmp_path, unsteady,
                                                                     monkeypatch):
    """Hands-on review of round 1: a server that stalled took two minutes a
    try before the download went on, twelve to give up."""
    from fastmdxplora.sharing import opening

    _shared(tmp_path, unsteady)
    monkeypatch.setattr(opening, "STALLED_S", 1)
    unsteady.faults = {"GET file": ["stall"]}
    told: list[str] = []
    opened = open_shared("10.5072/zenodo.7", workspace=tmp_path / "w", said=told.append)
    assert (opened / "shared_from.json").is_file()
    assert unsteady.ranges[1].startswith("bytes=")
    assert opening.STALLED_S == 1 and zenodo.TRIES == 5


def test_a_refused_open_leaves_no_folder_it_made(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.wrong_md5 = True
    with pytest.raises(ShareRefused):
        open_shared("10.5072/zenodo.7", into=tmp_path / "new" / "place")
    assert not (tmp_path / "new" / "place").exists()


def test_asking_zenodo_again_for_a_record_is_said(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET record": [503]}
    told: list[str] = []
    open_shared("10.5072/zenodo.7", workspace=tmp_path / "w", said=told.append)
    assert any("/api/records/7 with HTTP 503; asking again in 1 s (2 of 5)" in line
               for line in told)


def test_giving_up_says_how_far_the_download_got(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET file": ["cut"] + [503] * zenodo.TRIES}
    with pytest.raises(ShareRefused) as refused:
        open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    size = len(unsteady.stored["archive"])
    assert str(refused.value).endswith(
        f"HTTP 503, still after {zenodo.TRIES} tries, with "
        f"{(size - size * 2 // 3) / 1e6:.1f} of {size / 1e6:.1f} MB downloaded. Nothing was "
        "kept; open it again later.")


def test_a_draft_the_person_made_is_never_taken(tmp_path, unsteady, capsys):
    """Review of round 2: a draft that appeared among the person's after a
    failed request was taken as this share's, theirs or not, and its files
    and description replaced."""
    unsteady.drafts = [{"id": 8, "title": "My own draft", "hidden": 0}]
    for fault in (400, 503, "drop"):
        unsteady.calls.clear()
        unsteady.faults = {"POST": [fault]}
        assert _share_to_zenodo(_study(tmp_path / str(fault)), tmp_path / f"{fault}.zip") == 2
        assert not any("/8" in c[1] for c in unsteady.calls), fault
    assert "may still make the draft (titled 'FastMDXplora: sharing" in capsys.readouterr().err


def test_a_draft_listed_late_is_still_the_one_taken(tmp_path, unsteady):
    """Review of round 2: a draft Zenodo listed only a moment after its
    gateway gave up was not found, and sharing again made a second."""
    unsteady.faults = {"POST": ["made late"]}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert [c[1] for c in unsteady.calls if c[0] == "POST"] == ["/api/deposit/depositions"]
    assert unsteady.waits[:2] == [2.0, 4.0]
    # Listed only after every look: found at the next share, not made again.
    unsteady.calls.clear()
    unsteady.drafts = []
    unsteady.faults = {"POST": ["made late"]}
    study = _study(tmp_path / "later")
    (study / "shared_to.json").unlink(missing_ok=True)
    zenodo._LOOKS_S, looks = (1.0,), zenodo._LOOKS_S
    try:
        assert _share_to_zenodo(study, tmp_path / "later.zip") == 2
        assert "asked_as" in json.loads((study / "shared_to.json").read_text())[
            "zenodo-sandbox"]
        # Zenodo has since finished making it.
        unsteady.drafts[0]["hidden"] = 0
        assert _share_to_zenodo(study, tmp_path / "later.zip") == 0
    finally:
        zenodo._LOOKS_S = looks
    assert [c[1] for c in unsteady.calls if c[0] == "POST"] == ["/api/deposit/depositions"]
    assert json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"]["id"] == "7"


def test_an_answer_whose_own_words_are_cut_is_a_refusal(tmp_path, unsteady, capsys):
    """Review of round 2: Zenodo's 504 cut short in its own words escaped as
    a traceback, the draft it may have made never looked for."""
    unsteady.faults = {"POST": ["504 cut"]}
    assert _share_to_zenodo(_study(tmp_path / "study"), tmp_path / "shared.zip") == 2
    assert "HTTP 504" in capsys.readouterr().err
    assert any(c[1].startswith("/api/deposit/depositions?") for c in unsteady.calls)


def test_a_cut_is_seen_where_zenodo_records_no_size(tmp_path, unsteady):
    unsteady.payload = os.urandom(3 << 20)
    unsteady.faults = {"GET file": ["cut"]}
    _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
              tmp_path / "got.zip", most=1 << 30, md5="", said=None, size=0)
    assert (tmp_path / "got.zip").read_bytes() == unsteady.payload
    assert unsteady.ranges[1].startswith("bytes=")


def test_a_kept_draft_is_said_as_it_is(tmp_path, unsteady, capsys):
    """Review of round 2: a draft whose new archive was up, its description
    not, was said to keep the archive it had."""
    unsteady.faults = {"PUT metadata": [400]}
    unsteady.files = []
    assert _share_to_zenodo(_study(tmp_path / "study"), tmp_path / "shared.zip") == 2
    assert "is kept, with the new archive, its description not yet written" in \
        capsys.readouterr().err


def test_a_draft_s_archive_of_the_same_name_is_replaced_and_others_go(tmp_path, unsteady):
    """Review of round 2: the stand-in's files had no names, so a test
    passed taking away the archive just uploaded."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.calls.clear()
    unsteady.files = [{"id": "same", "filename": "shared.zip"},
                      {"id": "other", "filename": "old-name.zip"}]
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    deleted = [c[1] for c in unsteady.calls if c[0] == "DELETE"]
    assert deleted == ["/api/deposit/depositions/7/files/other"]
    put = next(i for i, c in enumerate(unsteady.calls) if c[1].startswith("/files/"))
    assert put < next(i for i, c in enumerate(unsteady.calls) if c[0] == "DELETE")


def test_every_folder_made_to_open_into_is_taken_away(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.wrong_md5 = True
    with pytest.raises(ShareRefused):
        open_shared("10.5072/zenodo.7", into=tmp_path / "new" / "deeper" / "place")
    assert not (tmp_path / "new").exists()


def _studied(tmp_path, name):
    return _study(tmp_path / name), tmp_path / f"{name}.zip"


def test_a_draft_asked_for_again_keeps_its_version(tmp_path, unsteady, capsys):
    """Review of round 3: a first share that failed made the next one say
    "version 2", one more for each try."""
    study, out = _studied(tmp_path, "study")
    for fault in (400, "drop"):
        unsteady.faults = {"POST": [fault]}
        assert _share_to_zenodo(study, out) == 2
    assert _share_to_zenodo(study, out) == 0
    assert unsteady.stored["metadata"]["version"] == "1"
    assert "Draft made:" in capsys.readouterr().out


def test_a_draft_not_listed_at_the_next_share_is_not_made_again(tmp_path, unsteady,
                                                                capsys):
    """Review of round 3: at the next share, a listing Zenodo did not answer,
    or one that did not yet show the draft, made a second."""
    study, out = _studied(tmp_path, "study")
    unsteady.late = 99
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, out) == 2
    unsteady.calls.clear()
    unsteady.faults = {"GET list": [503]}
    assert _share_to_zenodo(study, out) == 2
    assert "Zenodo did not list your drafts" in capsys.readouterr().err
    assert not any(c[0] == "POST" for c in unsteady.calls)
    # Listed, still without it: asked for again under the same title, and
    # once both are seen, named for the person rather than a third made.
    unsteady.faults = {}
    assert _share_to_zenodo(study, out) == 0
    titles = {d["title"] for d in unsteady.drafts}
    assert len(unsteady.drafts) == 2 and len(titles) == 1
    unsteady.drafts[0]["hidden"] = 0
    (study / "shared_to.json").write_text(json.dumps(
        {"zenodo-sandbox": {"asked_as": titles.pop(), "version": 1}}))
    assert _share_to_zenodo(study, out) == 2
    assert "Zenodo holds 2 drafts titled" in capsys.readouterr().err


def test_the_looks_for_a_lost_draft_are_few_and_short(tmp_path, unsteady, monkeypatch):
    """Review of round 3: each look asked Zenodo five times, waiting what it
    asked, and the last wait was followed by no look."""
    unsteady.late = len(zenodo._LOOKS_S) - 1
    unsteady.faults = {"POST": ["made late"]}
    timeouts = []
    real = zenodo.urllib.request.urlopen

    def urlopen(request, timeout=None):
        if "?status=draft" in request.full_url:
            timeouts.append(timeout)
        return real(request, timeout=timeout)

    monkeypatch.setattr(zenodo.urllib.request, "urlopen", urlopen)
    study, out = _studied(tmp_path, "study")
    assert _share_to_zenodo(study, out) == 0
    assert len(timeouts) == len(zenodo._LOOKS_S)
    assert set(timeouts) == {zenodo._LOOK_TIMEOUT_S}
    assert sum(unsteady.waits) == sum(zenodo._LOOKS_S)
    unsteady.calls.clear()
    unsteady.faults = {"POST": ["made late"], "GET list": [429] * 9}
    unsteady.late = 99
    assert _share_to_zenodo(*_studied(tmp_path, "slow")) == 2
    looks = [c for c in unsteady.calls if "?status=draft" in c[1]]
    assert len(looks) == len(zenodo._LOOKS_S)


def test_a_connection_never_made_is_not_taken_for_a_lost_answer(tmp_path, unsteady,
                                                                monkeypatch):
    """Review of round 3: a connection that could not be made, nothing sent,
    was looked after for a draft for a minute and a half."""
    import urllib.error

    real = zenodo.urllib.request.urlopen

    def urlopen(request, timeout=None):
        if request.get_method() == "POST":
            raise urllib.error.URLError(OSError(101, "Network is unreachable"))
        return real(request, timeout=timeout)

    monkeypatch.setattr(zenodo.urllib.request, "urlopen", urlopen)
    assert _share_to_zenodo(*_studied(tmp_path, "study")) == 2
    assert not any("?status=draft" in c[1] for c in unsteady.calls)
    assert unsteady.waits == []


def test_a_disk_that_fills_unpacking_or_after_is_said(tmp_path, unsteady, monkeypatch):
    """Hands-on review of round 3: a disk that filled unpacking was said to
    be an unsafe archive, and one that filled writing the study's records
    ended in a traceback."""
    import errno

    from fastmdxplora.sharing import opening

    _shared(tmp_path, unsteady)
    archive = tmp_path / "shared.zip"

    def full(*_args, **_kwargs):
        raise OSError(errno.ENOSPC, "No space left on device")

    for step in ("_unpacked", "_put_back"):
        monkeypatch.setattr(opening, step, full)
        with pytest.raises(ShareRefused) as refused:
            open_shared(str(archive), into=tmp_path / "here")
        assert refused.value.refusal.code == "environment.share.too_large", step
        assert str(refused.value).startswith(
            f"The shared study could not be written to {tmp_path / 'here'}: the disk there "
            "has no room for it"), step
        monkeypatch.undo()
        monkeypatch.setattr(zenodo, "_sleep", unsteady.waits.append)


def test_two_archives_left_on_a_draft_are_said_before_it_is_published(tmp_path, unsteady,
                                                                       capsys):
    """Hands-on review of round 3: an earlier archive Zenodo would not take
    away was left beside the new one, unsaid; published so, the record
    opens as neither."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    capsys.readouterr()
    unsteady.files = [{"id": "other", "filename": "old-name.zip"}]
    unsteady.faults = {"DELETE": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "with the new archive and the earlier old-name.zip, which Zenodo did not take away" \
        in said and "before you publish" in said


def test_a_range_refused_is_read_from_the_start(tmp_path, unsteady):
    _shared(tmp_path, unsteady)
    unsteady.faults = {"GET file": ["cut", 416]}
    opened = open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert (opened / "shared_from.json").is_file()
    assert unsteady.ranges[-1] is None


def test_what_to_do_is_said_when_zenodo_stays_down(tmp_path, unsteady, capsys):
    """Hands-on review of round 3: giving up on the record, or on the draft
    a study was shared as, said what failed but not what to do."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.faults = {"GET record": [503] * zenodo.TRIES}
    with pytest.raises(ShareRefused) as refused:
        open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert str(refused.value).endswith("Nothing was kept; open it again later.")
    capsys.readouterr()
    unsteady.faults = {"GET deposition": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "The draft or record 7 the study was shared as is as it was; share again once " \
        "Zenodo answers." in capsys.readouterr().err
