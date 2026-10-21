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
    answer). It lists the person's drafts as Zenodo does."""

    faults: dict = {}
    ranges: list = []
    payload: bytes | None = None
    drafts: list = []

    def _fault(self, which):
        queue = type(self).faults.get(which) or []
        return queue.pop(0) if queue else None

    def _drop(self):
        self.close_connection = True
        self.connection.shutdown(2)

    def _status(self, status):
        self.send_response(status)
        if status == 429:
            self.send_header("Retry-After", "7")
        self.send_header("Content-Length", "0")
        self.end_headers()

    def do_POST(self):  # noqa: N802
        fault = self._fault("POST")
        if self.path == "/api/deposit/depositions" and fault in (None, "made"):
            type(self).drafts.append(7)
        if fault == "made":
            type(self).calls.append(("POST", self.path, self.headers.get("Authorization")))
            self._body()
            return self._status(504)
        if fault == "drop":
            type(self).calls.append(("POST", self.path, self.headers.get("Authorization")))
            self._body()
            return self._drop()
        if isinstance(fault, int):
            type(self).calls.append(("POST", self.path, self.headers.get("Authorization")))
            self._body()
            return self._status(fault)
        return super().do_POST()

    def do_PUT(self):  # noqa: N802
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
            type(self).calls.append(("GET", self.path, self.headers.get("Authorization")))
            return self._json([{"id": n} for n in type(self).drafts])
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
    _Unsteady.payload, _Unsteady.drafts = None, []
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
    opened = open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert (opened / "shared_from.json").is_file()


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

    class Full:
        def __init__(self, handle):
            self.handle = handle

        def write(self, piece):
            writes.append(len(piece))
            if len(writes) == 2:
                raise OSError(errno.ENOSPC, "No space left on device")
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
    with pytest.raises(OSError) as failed:
        _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
                  target, most=1 << 30, md5=hashlib.md5(unsteady.payload).hexdigest(),
                  said=None, size=len(unsteady.payload))
    assert failed.value.errno == errno.ENOSPC
    assert unsteady.ranges == [None]


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
    assert unsteady.drafts == [7]
    assert unsteady.stored["archive"] == (tmp_path / "shared.zip").read_bytes()
    assert any("it made draft 7; taking that one" in line for line in told)


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
        f"{(size - size * 2 // 3) / 1e6:.1f} of {size / 1e6:.1f} MB downloaded.")
