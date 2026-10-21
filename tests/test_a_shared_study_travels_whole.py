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
import re
import threading
from pathlib import Path

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

    gone: set = set()
    uploaded_to: dict = {}
    edited: dict = {}
    latest: dict | None = None
    files_of: dict = {}
    # Drafts of a new version of record 7 there is (Zenodo keeps one).
    versions_begun: list = []
    # Set by a test: each draft keeps the metadata it was made and written
    # with, as Zenodo does, and only the drafts named are published.
    meta: dict | None = None
    published_ids: set | None = None

    def _deposition(self, number):
        made = super()._deposition(number)
        # A draft holds what was uploaded to its own bucket (one made since
        # holds none), where a test does not say.
        uploaded = type(self).uploaded_to.get(number)
        made["files"] = [{"id": "f1", "filename": uploaded}] if uploaded else []
        if type(self).files is not None:
            made["files"] = type(self).files
        if type(self).meta is not None:
            made["metadata"].update(type(self).meta.get(number, {}))
        if type(self).published_ids is not None:
            made["submitted"] = number in type(self).published_ids
        made["metadata"].update(type(self).edited.get(number, {}))
        if number in type(self).files_of:
            made["files"] = type(self).files_of[number]
        if type(self).latest is not None:
            made["links"]["latest"] = f"{self._base()}/api/records/latest"
        made["conceptrecid"] = "c7"
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
        shown += [self._deposition(n) for n in type(self).versions_begun]
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
        latest = type(self).latest
        if self.path.endswith("/actions/newversion") and latest is not None \
                and f"/{latest['id']}/actions/" not in self.path:
            # As Zenodo documents: a new version only of the latest one.
            type(self).calls.append(("POST", self.path, None))
            self._body()
            return self._json({"message": "not the latest version"}, 400)
        if self.path.endswith("/actions/newversion"):
            # A new version's draft begun again, the one deleted replaced;
            # with the files of the version it is of.
            type(self).gone.clear()
            if 7 in type(self).uploaded_to:
                type(self).uploaded_to[8] = type(self).uploaded_to[7]
            if type(self).meta is not None:
                # As Zenodo begins one: the latest version's description.
                type(self).meta[8] = {k: v for k, v in type(self).meta.get(7, {}).items()
                                      if k not in ("doi", "publication_date", "version")}
            if self._fault("POST newversion") == "made":
                type(self).versions_begun.append(8)
                type(self).calls.append(("POST", self.path, None))
                self._body()
                return self._status(504)
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
        if self.path == "/api/deposit/depositions" and fault is None:
            type(self).calls.append(("POST", self.path, self.headers.get("Authorization")))
            made = type(self).drafts[-1]
            if type(self).meta is not None:
                type(self).meta[made["id"]] = {"title": made["title"]}
            return self._json(self._deposition(made["id"]), 201)
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
            if type(self).meta is not None:
                body = self._body()
                number = int(self.path.rsplit("/", 1)[1])
                sent = dict(json.loads(body)["metadata"])
                sent.pop("prereserve_doi", None)
                type(self).meta[number] = {**type(self).meta.get(number, {}), **sent}
                type(self).stored["metadata"] = json.loads(body)["metadata"]
                type(self).calls.append(("PUT", self.path, self.headers.get("Authorization")))
                if fault == "applied":
                    # Kept, and its answer lost.
                    return self._status(500)
                return self._json(self._deposition(number))
        if self.path.startswith("/files/"):
            fault = self._fault("PUT file")
            if isinstance(fault, int):
                type(self).calls.append(("PUT", self.path, self.headers.get("Authorization")))
                self._body()
                return self._status(fault)
            bucket, name = self.path.split("/")[2:4]
            type(self).uploaded_to[int(bucket.lstrip("b"))] = name
            # A file of that name a test lists as the upload begins holds what
            # was uploaded (one a test lists once it is answered is the test's).
            listed = [*(type(self).files or []),
                      *type(self).files_of.get(int(bucket.lstrip("b")), [])]
            super().do_PUT()
            sent = hashlib.md5(type(self).stored["archive"]).hexdigest()  # noqa: S324
            for entry in listed:
                if entry.get("filename") == name and "checksum" in entry:
                    entry["checksum"] = sent
            return None
        return super().do_PUT()

    def do_GET(self):  # noqa: N802
        archive = type(self).payload or type(self).stored.get("archive", b"")
        name = type(self).stored.get("name", "shared.zip")
        if self.path.startswith("/api/deposit/depositions?"):
            return self._list()
        if self.path == "/api/records/latest" and type(self).latest is not None:
            return self._json(type(self).latest)
        latest = type(self).latest
        if latest is not None and self.path == f"/api/deposit/depositions/{latest['id']}":
            made = self._deposition(7)
            made.update(id=latest["id"], submitted=True,
                        files=[{"id": f["key"], "filename": f["key"]} for f in latest["files"]])
            made["metadata"].update(latest["metadata"])
            return self._json(made)
        if self.path.startswith("/api/deposit/depositions/"):
            fault = self._fault("GET deposition")
            if isinstance(fault, int):
                return self._status(fault)
            if self.path.rsplit("/", 1)[1] in {str(n) for n in type(self).gone}:
                type(self).calls.append(("GET", self.path, None))
                return self._json({"message": "not found"}, 404)
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
    _Unsteady.late, _Unsteady.gone, _Unsteady.edited = 2, set(), {}
    _Unsteady.latest, _Unsteady.files_of, _Unsteady.versions_begun = None, {}, []
    _Unsteady.meta, _Unsteady.published_ids = None, None
    _Unsteady.uploaded_to = {}
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
    assert any("Zenodo made draft 7, asked for by this study; taking that one" in line
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
    assert "/uploads/7 is kept, with the archive it had" in said
    assert "share again once Zenodo answers" in said
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
    # Another study's draft, not the first one's (the stand-in has but one).
    unsteady.stored.clear()
    unsteady.uploaded_to.clear()
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
    _packed_before_as(study, "old-name.zip")
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
    unsteady.faults = {"GET list": [503] * zenodo.TRIES}
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
    # The looks, and one more once the share is done, for a draft made twice.
    assert len(timeouts) == len(zenodo._LOOKS_S) + 1
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
    _packed_before_as(study, "old-name.zip")
    unsteady.faults = {"DELETE": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "with the new archive and the earlier old-name.zip, which Zenodo may not have " \
        "taken away" in said and "before you publish" in said


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


def test_a_record_that_will_never_open_is_not_said_to_open_later(tmp_path, unsteady):
    """Review of the fourth commit: a DOI naming no record, or a record of no
    one study, was told to open it again later."""
    with pytest.raises(ShareRefused) as refused:
        open_shared("10.5072/zenodo.12345", workspace=tmp_path / "w")
    assert "again later" not in str(refused.value)
    assert str(refused.value).startswith("Zenodo has no published record 12345: look at "
                                         "the DOI or address given")
    assert refused.value.refusal.code == "environment.share.not_a_study"


def test_a_draft_zenodo_listed_only_late_is_named_to_delete(tmp_path, unsteady):
    """Review of the fourth commit: a draft Zenodo listed only after it was
    asked for again was left beside the second, never named."""
    from fastmdxplora.sharing.pack import share_study

    study = _study(tmp_path / "study")
    unsteady.late = 99
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.late = 2
    told: list[str] = []
    zenodo.draft(study, tmp_path / "shared.zip", site="zenodo-sandbox",
                 pack=lambda doi: share_study(study, tmp_path / "shared.zip",
                                              authors=[("A B", "")], identifier=doi),
                 authors=[("A B", "")], license_id="CC-BY-4.0", system="1ABC",
                 said=told.append)
    unsteady.drafts[0]["hidden"] = 0
    told.clear()
    zenodo.draft(study, tmp_path / "shared.zip", site="zenodo-sandbox",
                 pack=lambda doi: share_study(study, tmp_path / "shared.zip",
                                              authors=[("A B", "")], identifier=doi),
                 authors=[("A B", "")], license_id="CC-BY-4.0", system="1ABC",
                 said=told.append)
    named = [line for line in told if line.startswith("Zenodo also holds draft 7 ")]
    assert named and "delete it among your uploads" in named[0]


def test_a_file_too_large_for_the_disk_is_said(tmp_path, unsteady, monkeypatch):
    import errno
    import io

    unsteady.payload = os.urandom(1 << 20)
    real = io.open
    target = tmp_path / "got.zip"

    class Fat:
        def __init__(self, handle):
            self.handle = handle

        def write(self, piece):
            raise OSError(errno.EFBIG, "File too large")

        def __getattr__(self, name):
            return getattr(self.handle, name)

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.handle.close()

    monkeypatch.setattr(io, "open", lambda file, mode="r", *a, **k: Fat(real(file, mode, *a, **k))
                        if str(file) == str(target) and "w" in mode else real(file, mode, *a, **k))
    with pytest.raises(ShareRefused) as refused:
        _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
                  target, most=1 << 30, md5="", said=None, size=len(unsteady.payload))
    assert "takes no file this large" in str(refused.value)


def test_a_token_refused_is_not_told_to_wait(tmp_path, unsteady, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    capsys.readouterr()
    unsteady.faults = {"GET deposition": [401]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "the token was refused" in said and "once Zenodo answers" not in said


def test_a_file_over_the_limit_is_refused_before_it_comes(tmp_path, unsteady):
    unsteady.payload = os.urandom(3 << 20)
    got = tmp_path / "got.zip"
    with pytest.raises(ShareRefused) as refused:
        _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
                  got, most=1 << 20, md5="", said=None, size=0)
    assert refused.value.refusal.code == "environment.share.too_large"
    assert got.stat().st_size == 0


def test_a_new_version_s_draft_deleted_is_begun_again_from_its_record(tmp_path, unsteady):
    """Hands-on review of the fourth commit: the draft of a new version
    deleted, the next share made a record of its own, version 1."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    kept = json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"]
    assert {k: kept[k] for k in ("id", "version", "record", "doi", "archive")} == {
        "id": "8", "version": "2", "record": "7", "doi": "10.5072/zenodo.8",
        "archive": "shared.zip"}
    unsteady.gone = {8}
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert any(c[1].endswith("/7/actions/newversion") for c in unsteady.calls)
    assert not any(c[1] == "/api/deposit/depositions" for c in unsteady.calls if c[0] == "POST")
    assert unsteady.stored["metadata"]["version"] == "2"


def test_a_new_version_the_person_began_and_changed_is_not_taken(tmp_path, unsteady, capsys):
    """Hands-on review of the fourth commit: a draft of a new version the
    person began and changed by hand had its title and files replaced."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    # Begun by the person on Zenodo's site, a co-author added: not the
    # study's, whatever in it was changed.
    unsteady.versions_begun = [8]
    unsteady.edited = {8: {"creators": [{"name": "Their, Co-author"}]}}
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "a draft of a new version of record 7 (8), begun on its site or by a share" in said
    assert not any(c[0] in ("PUT", "DELETE") or c[1].endswith("newversion")
                   for c in unsteady.calls)


def test_a_file_the_person_put_beside_the_archive_stays(tmp_path, unsteady):
    """Hands-on review of the fourth commit: a paper the person added to the
    draft was taken away at the next share; only another archive is."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.calls.clear()
    unsteady.files = [{"id": "p", "filename": "paper.pdf"}, {"id": "z", "filename": "old.zip"}]
    _packed_before_as(study, "old.zip")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert [c[1] for c in unsteady.calls if c[0] == "DELETE"] == [
        "/api/deposit/depositions/7/files/z"]


def test_a_refusal_says_what_to_do(tmp_path, unsteady, capsys, monkeypatch):
    """Hands-on review of the fourth commit: a refused token, or a draft
    never asked for, said what failed and not what to do."""
    import urllib.error

    unsteady.faults = {"POST": [401]}
    assert _share_to_zenodo(_study(tmp_path / "a"), tmp_path / "a.zip") == 2
    assert "with the deposit:write scope, and give it in ZENODO_TOKEN" in \
        capsys.readouterr().err
    real = zenodo.urllib.request.urlopen

    def urlopen(request, timeout=None):
        if request.get_method() == "POST":
            raise urllib.error.URLError(ConnectionRefusedError(111, "Connection refused"))
        return real(request, timeout=timeout)

    monkeypatch.setattr(zenodo.urllib.request, "urlopen", urlopen)
    assert _share_to_zenodo(_study(tmp_path / "b"), tmp_path / "b.zip") == 2
    assert "Nothing was made; share again once Zenodo answers." in capsys.readouterr().err


def _packed_before_as(study, name: str) -> None:
    """The study's last share named its archive ``name``."""
    path = study / "shared_to.json"
    kept = json.loads(path.read_text())
    kept["zenodo-sandbox"]["archive"] = name
    path.write_text(json.dumps(kept))


def test_an_archive_the_person_put_on_the_draft_is_never_taken_away(tmp_path, unsteady,
                                                                     capsys):
    """Review of the fifth commit: a zip the person put on the draft was
    taken away without a word as if it were an earlier archive."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.calls.clear()
    unsteady.files = [{"id": "m", "filename": "supplementary-movies.zip"}]
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "holds supplementary-movies.zip, which this study did not put there" in \
        capsys.readouterr().err
    assert not any(c[0] in ("DELETE",) or c[1].startswith("/files/") for c in unsteady.calls)


def test_a_new_version_is_compared_with_and_counted_from_the_latest_published(
        tmp_path, unsteady):
    """Review of the fifth commit: a version the person published by hand
    since made every new version look changed, and its number was the
    study's own count, not the record's."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    # Published by hand since as version 2, with a paper: the new draft
    # Zenodo begins from it carries the paper too.
    unsteady.latest = {"id": 9, "metadata": {"version": "2"},
                       "files": [{"key": "shared.zip"}, {"key": "paper.pdf"}]}
    unsteady.files_of = {8: [{"id": "z", "filename": "shared.zip"},
                             {"id": "p", "filename": "paper.pdf"}]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.stored["metadata"]["version"] == "3"


def test_a_refusal_after_a_deleted_new_version_names_its_record(tmp_path, unsteady,
                                                                 capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.gone = {8}
    unsteady.faults = {"POST": [504]}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "Zenodo may have begun the new version of record 7" in said
    assert "record 8" not in said


def test_a_token_refused_at_the_next_share_is_said_as_such(tmp_path, unsteady, capsys):
    study = _study(tmp_path / "study")
    unsteady.late = 99
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {"GET list": [401]}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "the token was refused" in said and "did not list your drafts" not in said


def test_a_file_too_large_unpacking_is_said(tmp_path, unsteady, monkeypatch):
    import errno

    from fastmdxplora.sharing import opening

    _shared(tmp_path, unsteady)

    def big(*_args, **_kwargs):
        raise OSError(errno.EFBIG, "File too large")

    monkeypatch.setattr(opening, "_unpacked", big)
    with pytest.raises(ShareRefused) as refused:
        open_shared(str(tmp_path / "shared.zip"), into=tmp_path / "fat")
    assert refused.value.refusal.code == "environment.share.too_large"
    assert "takes no file this large" in str(refused.value)


def test_a_new_version_begun_whose_answer_was_lost_is_left_for_the_person(tmp_path,
                                                                         unsteady, capsys):
    """Hands-on review of the sixth commit: a draft of a new version a share
    began and lost the answer to was taken later, with the person's edits
    to it written over. It cannot be told from one the person began: it is
    left for them, and once they delete it the next share begins one."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.faults = {"POST newversion": ["made"]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "the next share finds its draft and asks you to delete it there" in \
        capsys.readouterr().err
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "begun on its site or by a share whose answer was lost" in capsys.readouterr().err
    assert not any(c[0] in ("PUT", "DELETE") or c[1].endswith("newversion")
                   for c in unsteady.calls)
    unsteady.versions_begun = []
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0


def test_the_person_s_edits_to_a_draft_not_yet_published_are_kept(tmp_path, unsteady,
                                                                   capsys):
    """Hands-on review of the fifth commit: shared again, the description the
    person had changed on the draft was written over."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.edited = {7: {"description": "<p>As I want it said.</p>"}}
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "it is kept as it is, and only the archive is new" in capsys.readouterr().out
    assert unsteady.stored["metadata"]["description"] == "<p>As I want it said.</p>"
    assert any(c[1].startswith("/files/") for c in unsteady.calls)


def test_a_full_disk_while_sharing_is_said(tmp_path, unsteady, monkeypatch, capsys):
    """Hands-on review of the fifth commit: a disk that filled writing the
    archive, the study's record of where it is shared, or the archive named
    with its DOI ended in a traceback, a draft left unsaid."""
    import errno

    from fastmdxplora.sharing import pack

    study = _study(tmp_path / "study")
    real = pack._packed
    calls = []

    def full(*args, **kwargs):
        calls.append(kwargs.get("identifier"))
        if len(calls) == 2:
            raise OSError(errno.ENOSPC, "No space left on device", str(tmp_path / "x.zip"))
        return real(*args, **kwargs)

    monkeypatch.setattr(pack, "_packed", full)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert f"could not be written to {tmp_path}: the disk there has no room" in said
    assert "/uploads/7 is kept, with no archive yet: share again" in said
    # Room again (the packing fails only the once): that draft finished.
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert not any(c[0] == "POST" for c in unsteady.calls)


def test_the_person_s_edits_stay_theirs_at_every_later_share(tmp_path, unsteady):
    """Review of the sixth commit: an edit kept at one share was written
    over at the next."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.edited = {7: {"description": "<p>As I want it said.</p>"}}
    for _ in range(2):
        unsteady.calls.clear()
        assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
        assert unsteady.stored["metadata"]["description"] == "<p>As I want it said.</p>"


def test_a_new_version_s_draft_not_listed_and_changed_is_left_alone(tmp_path, unsteady,
                                                                     capsys):
    """Review of the sixth commit: a draft of a new version the person began
    and changed, left out of Zenodo's list, was given back to the share's
    request and written over."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.edited = {8: {"title": "Theirs"}}
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "(8) that has been changed (its title)" in capsys.readouterr().err
    assert not any(c[0] in ("PUT", "DELETE") for c in unsteady.calls)


def test_a_new_version_refused_outright_is_not_taken_as_asked(tmp_path, unsteady):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.faults = {"POST": [403]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    kept = json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"]
    assert "asked_version_of" not in kept


def test_a_temporary_folder_with_no_room_is_named_as_such(tmp_path, monkeypatch):
    import errno
    import tempfile

    from fastmdxplora.sharing import pack

    def full(*_args, **_kwargs):
        raise OSError(errno.ENOSPC, "No space left on device",
                      str(Path(tempfile.gettempdir()) / "fastmdx-share-x1" / "0.txt"))

    monkeypatch.setattr(pack, "_packed", full)
    with pytest.raises(ShareRefused) as refused:
        pack.share_study(tmp_path, tmp_path / "o.zip", authors=[("A B", "")])
    assert str(refused.value).startswith(f"The temporary folder {tempfile.gettempdir()}, ")
    assert "set TMPDIR" in str(refused.value)


def test_a_token_refused_naming_a_stray_draft_fails_no_share(tmp_path, unsteady):
    study = _study(tmp_path / "study")
    unsteady.late = 99
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {"GET list": [None, 401]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"]["archive"] \
        == "shared.zip"


def test_the_person_s_archive_is_kept_after_a_share_that_failed(tmp_path, unsteady, capsys):
    """Hands-on review of the sixth commit: a share that failed part way
    forgot the archive it had put there, and the next took the person's
    supplement.zip away without a word."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"][
        "archive"] == "shared.zip"
    unsteady.files = [{"id": "s", "filename": "shared.zip"},
                      {"id": "m", "filename": "supplement.zip"}]
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "holds supplement.zip, which this study did not put there" in \
        capsys.readouterr().err
    assert not any(c[0] == "DELETE" for c in unsteady.calls)


def test_a_draft_the_share_had_not_yet_described_keeps_the_person_s_words(tmp_path,
                                                                           unsteady):
    """Hands-on review of the sixth commit: a draft whose description the
    share could not write had the words the person wrote there since
    written over."""
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT metadata": [503] * (zenodo.TRIES + 1)}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.edited = {7: {"title": "My title", "description": "<p>Mine.</p>"}}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.stored["metadata"]["description"] == "<p>Mine.</p>"


def test_a_published_record_s_edits_are_carried_into_its_new_version(tmp_path, unsteady,
                                                                      capsys):
    """Hands-on review of the sixth commit: what the person changed on a
    record before publishing it was dropped from every new version."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.edited = {7: {"description": "<p>Theirs.</p>"},
                       8: {"description": "<p>Theirs.</p>"}}
    unsteady.submitted = True
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "it is kept as it is, and only the archive is new" in capsys.readouterr().out
    assert unsteady.stored["metadata"]["description"] == "<p>Theirs.</p>"
    assert unsteady.stored["metadata"]["version"] == "2"


def test_a_disk_with_room_for_one_archive_makes_it_again(tmp_path, monkeypatch):
    """Hands-on review of the sixth commit: packed again with its DOI, an
    archive needed room for two, the one before still there as it wrote."""
    import zipfile

    from fastmdxplora.sharing import pack

    study = _study(tmp_path / "study")
    out = tmp_path / "shared.zip"
    out.write_bytes(b"an archive from before")
    there = []
    real = zipfile.ZipFile

    def opened(file, mode="r", *args, **kwargs):
        if mode == "w":
            there.append(out.exists())
        return real(file, mode, *args, **kwargs)

    monkeypatch.setattr(pack.zipfile, "ZipFile", opened)
    pack.share_study(study, out, authors=[("A B", "")])
    assert there == [False]


def test_a_draft_the_study_made_is_described_once_a_share_stopped_after_making_it(
        tmp_path, unsteady):
    """Review of the seventh commit: a draft the study made, its share
    stopped before describing it, was left titled as it was asked for,
    with no description, at every later share."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert unsteady.meta[7]["title"].startswith("FastMDXplora: sharing study (")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[7]["title"] == "The 1ABC study"
    assert unsteady.meta[7]["description"]


def test_a_new_draft_after_one_deleted_is_described(tmp_path, unsteady):
    """Review of the seventh commit: a new draft, made after the person
    deleted the one the study was shared as, was held against the deleted
    one's description and never described."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.gone = {7}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[8]["title"] == "The 1ABC study"


def test_a_new_version_shared_again_and_deleted_is_begun_again_from_its_record(tmp_path,
                                                                               unsteady):
    """Review of the seventh commit: a draft of a new version shared again
    forgot its record, so once deleted the next share made a record of its
    own."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.published_ids = {7}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"][
        "record"] == "7"
    unsteady.gone = {8}
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert any(c[1].endswith("/7/actions/newversion") for c in unsteady.calls)
    assert not any(c[0] == "POST" and c[1] == "/api/deposit/depositions"
                   for c in unsteady.calls)


def _analysed_again(study: Path, words: str) -> None:
    """The study's report sums it up in ``words`` now, and so does its
    description."""
    report = study / "report" / "report.md"
    text = report.read_text()
    assert "The system was simulated for 4 ps." in text
    report.write_text(text.replace("The system was simulated for 4 ps.", words))


def test_a_description_written_whose_answer_was_lost_is_still_the_study_s(tmp_path,
                                                                           unsteady):
    """Hands-on review of the seventh commit: a description Zenodo kept,
    its answer a 500, was taken for the person's at every later share."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT metadata": ["applied"] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    _analysed_again(study, "Analysed again, a longer run.")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Analysed again, a longer run." in unsteady.meta[7]["description"]
    assert json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"]["written"]


def test_a_new_version_whose_share_failed_is_described_at_the_next(tmp_path, unsteady):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.published_ids = {7}
    _analysed_again(study, "Analysed again, a longer run.")
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[8].get("version") == "2"
    assert "Analysed again, a longer run." in unsteady.meta[8]["description"]


def test_an_archive_being_replaced_stays_the_study_s_until_it_is_gone(tmp_path, unsteady):
    """Hands-on review of the seventh commit: the old archive Zenodo would
    not take away refused every later share as not the study's."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.files = [{"id": "a", "filename": "a.zip"}]
    unsteady.faults = {"DELETE": [500] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    unsteady.faults = {}
    unsteady.files = [{"id": "a", "filename": "a.zip"}, {"id": "b", "filename": "b.zip"}]
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 0
    assert [c[1] for c in unsteady.calls if c[0] == "DELETE"] == [
        "/api/deposit/depositions/7/files/a"]


def test_a_version_labelled_with_letters_is_counted_by_its_number(tmp_path, unsteady):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.latest = {"id": 9, "metadata": {"version": "v5"},
                       "files": [{"key": "shared.zip"}]}
    unsteady.files_of = {8: [{"id": "z", "filename": "shared.zip"}]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.stored["metadata"]["version"] == "v6"


def test_where_a_study_is_shared_is_written_whole_or_not_at_all(tmp_path, monkeypatch):
    import errno

    study = tmp_path / "study"
    study.mkdir()
    (study / "shared_to.json").write_text('{"zenodo": {"id": "7"}}')
    real = Path.write_text

    def full(self, *args, **kwargs):
        if self.name.endswith(".tmp"):
            raise OSError(errno.ENOSPC, "No space left on device")
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", full)
    with pytest.raises(ShareRefused):
        zenodo._keep(study, {"zenodo": {"id": "7", "record": "5"}})
    assert json.loads((study / "shared_to.json").read_text()) == {"zenodo": {"id": "7"}}
    assert not list(study.glob(".*.tmp"))


def test_a_new_version_begun_again_is_described_by_the_study(tmp_path, unsteady):
    """Hands-on review of the seventh commit: a new version begun again
    after the person deleted the share's kept the record's description, its
    own held against the deleted draft's."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.published_ids = {7}
    _analysed_again(study, "Analysed again, a longer run.")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.gone = {8}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[8].get("version") == "2"
    assert "Analysed again, a longer run." in unsteady.meta[8]["description"]


def test_a_lost_answer_then_published_leaves_no_later_version_stale(tmp_path, unsteady):
    """Review of the eighth commit: a description written whose answer was
    lost, the draft then published, left every later version with the
    record's old description."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT metadata": ["applied"] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.published_ids = {7}
    _analysed_again(study, "Words of version two.")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Words of version two." in unsteady.meta[8]["description"]


def test_the_person_s_edits_after_a_lost_answer_are_kept(tmp_path, unsteady):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT metadata": ["applied"] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.meta[7]["creators"] = [{"name": "Aina, Adekunle"}, {"name": "Their, Co-author"}]
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert len(unsteady.meta[7]["creators"]) == 2


def test_a_version_label_is_counted_on_in_its_own_form():
    assert [zenodo._bumped(v) for v in ("1", "1.0", "2.1.5", "v5", "draft")] == [
        "2", "1.1", "2.1.6", "v6", ""]


def test_every_archive_being_replaced_stays_the_study_s(tmp_path, unsteady):
    """Review of the eighth commit: two earlier archives Zenodo would not
    take away, and only the last was remembered as the study's."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.files = [{"id": "a", "filename": "a.zip"}]
    unsteady.faults = {"DELETE": [500] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    unsteady.files = [{"id": "a", "filename": "a.zip"}, {"id": "b", "filename": "b.zip"}]
    unsteady.faults = {"DELETE": [None] + [500] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 2
    unsteady.faults = {}
    unsteady.files = [{"id": "a", "filename": "a.zip"}, {"id": "c", "filename": "c.zip"}]
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 0
    assert [c[1] for c in unsteady.calls if c[0] == "DELETE"] == [
        "/api/deposit/depositions/7/files/a"]


def test_an_archive_whose_upload_answer_was_lost_is_the_study_s(tmp_path, unsteady):
    """Review of the eighth commit: an archive Zenodo kept whose answer was
    lost was, at a later share under another name, said not the study's."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.files = [{"id": "a", "filename": "a.zip"}]
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    unsteady.faults = {}
    unsteady.files = [{"id": "a", "filename": "a.zip"}, {"id": "b", "filename": "b.zip"}]
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 0
    assert sorted(c[1] for c in unsteady.calls if c[0] == "DELETE") == [
        "/api/deposit/depositions/7/files/a", "/api/deposit/depositions/7/files/b"]


def test_a_description_sent_twice_and_lost_is_still_the_study_s(tmp_path, unsteady,
                                                                 capsys):
    """Hands-on review of the eighth commit: a description that landed, its
    answer lost, then another that never landed, left the draft taken for
    the person's for good."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _analysed_again(study, "Words of the second share.")
    unsteady.faults = {"PUT metadata": ["applied"] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    report = study / "report" / "report.md"
    report.write_text(report.read_text().replace("Words of the second share.",
                                                 "Words of the third share."))
    unsteady.faults = {"PUT metadata": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "changed there since this study last described it" not in capsys.readouterr().out
    assert "Words of the third share." in unsteady.meta[7]["description"]


def test_one_429_on_the_listing_before_a_new_version_is_asked_again(tmp_path, unsteady):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.faults = {"GET list": [429]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0


def test_names_the_person_gave_a_draft_not_yet_described_are_kept(tmp_path, unsteady,
                                                                   capsys):
    """Review of the ninth commit: a co-author and keywords the person added
    to a draft whose first upload failed were written over at the next
    share, the draft taken as still the study's to describe."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.meta[7]["creators"] = [{"name": "Aina, Adekunle"}, {"name": "Their, Co-author"}]
    unsteady.meta[7]["keywords"] = ["their keyword"]
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert len(unsteady.meta[7]["creators"]) == 2
    assert unsteady.meta[7]["keywords"] == ["their keyword"]
    assert "changed there before this study described it" in capsys.readouterr().out


def test_a_new_version_begun_again_holds_the_record_s_archive_as_the_study_s(tmp_path,
                                                                            unsteady):
    """Review of the ninth commit: a new version shared under another name,
    then deleted, was begun again holding the record's archive, which the
    share refused as not the study's."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.published_ids = {7}
    # As Zenodo begins a new version: with the files of the one before.
    unsteady.files_of = {n: [{"id": "a", "filename": "a.zip"}] for n in (7, 8)}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 0
    unsteady.gone = {8}
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 0
    assert any(c[1].endswith("/7/actions/newversion") for c in unsteady.calls)
    assert [c[1] for c in unsteady.calls if c[0] == "DELETE"] == [
        "/api/deposit/depositions/8/files/a"]


def test_a_draft_holding_only_the_person_s_files_is_not_said_to_hold_an_archive(
        tmp_path, unsteady, capsys):
    """Review of the ninth commit: a draft holding only the person's paper,
    its upload failing, was said kept "with the archive it had"."""
    study = _study(tmp_path / "study")
    unsteady.files = [{"id": "p", "filename": "paper.pdf"}]
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "with the archive it had" not in said
    assert "with no archive yet, or the new one if Zenodo kept it" in said


def test_one_429_on_the_listing_after_a_lost_draft_is_asked_again(tmp_path, unsteady):
    """Review of the ninth commit: at the share after a draft whose answer
    was lost, one 429 on the listing refused the share."""
    study, out = _studied(tmp_path, "study")
    unsteady.late = len(zenodo._LOOKS_S)
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, out) == 2
    unsteady.calls.clear()
    unsteady.faults = {"GET list": [429]}
    assert _share_to_zenodo(study, out) == 0
    assert not any(c[0] == "POST" for c in unsteady.calls)
    assert len(unsteady.drafts) == 1


def test_a_description_that_landed_is_the_study_s_whatever_zenodo_added(tmp_path,
                                                                        unsteady):
    """Review of the ninth commit: a description written whose answer was
    lost was told by every field of the draft, so one Zenodo added there
    made it the person's for good."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT metadata": ["applied"] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.meta[7]["imprint_publisher"] = "Zenodo"
    _analysed_again(study, "Words of the next share.")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Words of the next share." in unsteady.meta[7]["description"]
    assert unsteady.meta[7]["imprint_publisher"] == "Zenodo"


def test_a_full_disk_refuses_before_a_new_version_is_begun(tmp_path, unsteady, monkeypatch,
                                                            capsys):
    """Hands-on review of the ninth commit: a published study whose disk
    was full began a new version it could not remember, then refused it at
    the next share as one the person began."""
    import errno

    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    real = Path.write_text

    def full(self, *args, **kwargs):
        if self.name.startswith(f".{zenodo.SHARED_TO}."):
            raise OSError(errno.ENOSPC, "No space left on device")
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", full)
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "the disk has no room" in capsys.readouterr().err
    assert not any(c[1].endswith("/actions/newversion") for c in unsteady.calls)
    monkeypatch.setattr(Path, "write_text", real)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0


def test_a_record_zenodo_did_not_answer_for_is_not_said_begun(tmp_path, unsteady,
                                                             monkeypatch, capsys):
    """Hands-on review of the ninth commit: a record's GET answering 503 was
    said to have perhaps begun a new version, though none was asked for."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True

    def unanswered(*_args, **_kwargs):
        refused = ShareRefused("Zenodo answered HTTP 503, still after 5 tries.",
                               code="environment.service.unreachable")
        refused.unclear = refused.later = True
        raise refused

    monkeypatch.setattr(zenodo, "_latest_published", unanswered)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "may have begun" not in said
    assert "Record 7 is as it was; share again once Zenodo answers." in said


def test_a_token_refused_part_way_is_not_told_to_wait(tmp_path, unsteady, capsys):
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT file": [401]}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "once Zenodo answers" not in said and "share again once that is mended" in said


def test_a_small_limit_is_said_in_megabytes(tmp_path, unsteady):
    unsteady.payload = os.urandom(12 << 20)
    with pytest.raises(ShareRefused) as refused:
        _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
                  tmp_path / "got.zip", most=10 ** 7, md5="", said=None, size=0)
    assert "over 10 MB" in str(refused.value) and "0.0 GB" not in str(refused.value)


def test_names_the_person_gave_a_draft_found_by_its_title_are_kept(tmp_path, unsteady):
    """Review of the tenth commit: a draft whose answer was lost, found by
    its title at the next share, was taken as made with the co-author and
    keywords the person had given it since, and both written over."""
    unsteady.meta = {}
    study, out = _studied(tmp_path, "study")
    unsteady.late = len(zenodo._LOOKS_S)
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, out) == 2
    unsteady.faults = {}
    unsteady.meta[7] = {"title": unsteady.drafts[0]["title"],
                        "creators": [{"name": "Aina, Adekunle"}, {"name": "Their, Co-author"}],
                        "keywords": ["their keyword"]}
    assert _share_to_zenodo(study, out) == 0
    assert len(unsteady.meta[7]["creators"]) == 2
    assert unsteady.meta[7]["keywords"] == ["their keyword"]


def test_a_new_version_begun_whose_draft_was_not_read_is_said_begun(tmp_path, unsteady,
                                                                   capsys):
    """Review of the eleventh commit: a new version Zenodo began, its draft
    then refused to be read, was said not begun ("Record 7 is as it was")."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.faults = {"GET deposition": [None, 403]}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "is as it was" not in said
    assert "Zenodo began the new version of record 7 (8)" in said
    assert "holding the record as it was and asks you to delete it there" in said


def test_a_new_version_begun_again_holds_whichever_archive_the_record_was_published_with(
        tmp_path, unsteady):
    """Review of the eleventh commit: a record published with the archive
    whose upload answer was lost, the old one deleted, refused that archive
    as not the study's once a new version was begun again."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    unsteady.faults = {}
    # The person deleted a.zip, kept b.zip and published.
    unsteady.published_ids = {7}
    unsteady.files_of = {n: [{"id": "b", "filename": "b.zip"}] for n in (7, 8)}
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 0
    unsteady.gone = {8}
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 0
    assert [c[1] for c in unsteady.calls if c[0] == "DELETE"] == [
        "/api/deposit/depositions/8/files/b"]


def test_an_answer_that_is_no_draft_is_not_taken_for_one_gone(tmp_path, unsteady,
                                                             monkeypatch, capsys):
    """Review of the eleventh commit: a draft's GET answered with an empty
    list was taken as the draft gone, and a second made beside it."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    real = zenodo._call

    def odd(method, url, **kwargs):
        if method == "GET" and url.endswith("/api/deposit/depositions/7"):
            return []
        return real(method, url, **kwargs)

    monkeypatch.setattr(zenodo, "_call", odd)
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "is no draft or record; nothing was made" in capsys.readouterr().err
    assert not any(c[0] == "POST" for c in unsteady.calls)


def test_a_new_version_begun_on_a_disk_that_then_fills_is_named(tmp_path, unsteady,
                                                                monkeypatch, capsys):
    """Hands-on review of the eleventh commit: a new version begun, the
    study's disk full just after, was not named, and the next share took it
    for one whose answer was lost."""
    import errno

    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    real, writes = Path.write_text, []

    def fills(self, *args, **kwargs):
        if self.name.startswith(f".{zenodo.SHARED_TO}."):
            writes.append(self)
            if len(writes) > 1:
                raise OSError(errno.ENOSPC, "No space left on device")
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", fills)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "Zenodo had begun a new version of record 7, draft 8" in said
    assert "delete that draft there, make room, then share again" in said


def test_a_new_version_left_as_begun_is_not_told_to_be_published(tmp_path, unsteady,
                                                                 capsys):
    """Hands-on review of the eleventh commit: a draft of a new version
    holding the record as it was was said to be published or deleted;
    published, it is the record again."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.published_ids = {7}
    unsteady.meta[8] = {k: v for k, v in unsteady.meta[7].items() if k != "version"}
    unsteady.uploaded_to[8] = unsteady.uploaded_to[7]
    unsteady.versions_begun = [8]
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "delete it there (it holds the record as it was)" in said
    assert "publish it" not in said


def test_names_changed_on_a_described_draft_are_not_said_a_description(tmp_path, unsteady,
                                                                       capsys):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.meta[7]["creators"] = [{"name": "Aina, Adekunle"}, {"name": "Their, Co-author"}]
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    said = capsys.readouterr().out
    assert "holds a description changed" not in said
    assert "changed there since this study last described it (its words, names, " \
        "keywords, access or licence)" in said


def test_a_new_version_shared_again_holds_the_record_s_archives_as_the_study_s(tmp_path,
                                                                              unsteady):
    """Review of the thirteenth commit: a new version's draft shared again
    after a failure took the archive the record was published with for one
    the study did not put there."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    unsteady.published_ids = {7}
    unsteady.files_of = {n: [{"id": "b", "filename": "b.zip"}] for n in (7, 8)}
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 2
    unsteady.faults = {}
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 0
    assert [c[1] for c in unsteady.calls if c[0] == "DELETE"] == [
        "/api/deposit/depositions/8/files/b"]


def test_a_new_version_with_a_file_replaced_by_the_person_is_theirs(tmp_path, unsteady,
                                                                    capsys):
    """Review of the thirteenth commit: a new version the person began and
    put a corrected a.zip on was said to hold the record as it was, to be
    deleted."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.published_ids = {7}
    unsteady.files_of = {7: [{"id": "a", "filename": "a.zip", "checksum": "md5:1"}],
                         8: [{"id": "a2", "filename": "a.zip", "checksum": "md5:2"}]}
    unsteady.meta[8] = {k: v for k, v in unsteady.meta[7].items() if k != "version"}
    unsteady.versions_begun = [8]
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 2
    said = capsys.readouterr().err
    assert "publish it, or delete it there" in said and "as it was" not in said


def test_the_person_s_choice_of_access_and_licence_is_kept(tmp_path, unsteady):
    """Review of the thirteenth commit: an embargo, licence and kind the
    person chose for a draft not yet described were written over."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.meta[7].update(access_right="embargoed", embargo_date="2027-01-01",
                            license="cc0-1.0", upload_type="software")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[7]["access_right"] == "embargoed"
    assert unsteady.meta[7]["license"] == "cc0-1.0"
    assert unsteady.meta[7]["upload_type"] == "software"


def test_a_doi_zenodo_did_not_reserve_says_the_draft_is_kept(tmp_path, unsteady, capsys):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.meta[7]["prereserve_doi"] = {}
    unsteady.faults = {"PUT metadata": [503] * zenodo.TRIES}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "is kept as it was, its DOI not yet reserved: share again once Zenodo answers" \
        in capsys.readouterr().err


def test_what_the_person_writes_while_the_archive_uploads_is_kept(tmp_path, unsteady,
                                                                  monkeypatch):
    """Hands-on review of the thirteenth commit: words the person wrote on
    the draft while its archive was uploading were written over."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    real = zenodo._call

    def meanwhile(method, url, **kwargs):
        answer = real(method, url, **kwargs)
        if kwargs.get("upload") is not None:
            unsteady.meta[7].update(description="<p>Typed by me during the upload.</p>",
                                    keywords=["mine"])
        return answer

    monkeypatch.setattr(zenodo, "_call", meanwhile)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[7]["description"] == "<p>Typed by me during the upload.</p>"
    assert unsteady.meta[7]["keywords"] == ["mine"]


def test_an_answer_that_is_no_json_is_said_to_be_asked_again(tmp_path, unsteady,
                                                             monkeypatch, capsys):
    import urllib.request

    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    real = urllib.request.urlopen

    class _Page:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def read(self, *_args):
            return b"<html>Bad gateway</html>"

    def gateway(request, timeout=None):
        if request.get_method() == "GET" and request.full_url.endswith("/depositions/7"):
            return _Page()
        return real(request, timeout=timeout)

    monkeypatch.setattr(urllib.request, "urlopen", gateway)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "was not JSON." in said and "share again once Zenodo answers" in said


def test_a_download_refused_is_not_told_to_be_asked_later(tmp_path, unsteady):
    unsteady.payload = os.urandom(1 << 16)
    unsteady.faults = {"GET file": [403]}
    with pytest.raises(ShareRefused) as refused:
        _download(f"{os.environ['FASTMDXPLORA_ZENODO_URL']}/records/7/files/shared.zip",
                  tmp_path / "got.zip", most=1 << 30, md5="", said=None, size=0)
    assert "ask its owner for access" in str(refused.value)
    assert "open it again later" not in str(refused.value)


def test_a_draft_found_by_its_title_with_zenodo_s_defaults_is_the_study_s(tmp_path,
                                                                         unsteady):
    """Review of the fifteenth commit: a draft whose answer was lost, found
    by its title holding the licence Zenodo gives a draft it makes, was
    taken for the person's and never described."""
    unsteady.meta = {}
    study, out = _studied(tmp_path, "study")
    unsteady.late = len(zenodo._LOOKS_S)
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, out) == 2
    unsteady.faults = {}
    unsteady.meta[7] = {"title": unsteady.drafts[0]["title"], "license": "cc-by",
                        "access_right": "open"}
    assert _share_to_zenodo(study, out) == 0
    assert unsteady.meta[7]["title"] == "The 1ABC study"
    assert unsteady.meta[7]["description"]


@pytest.mark.parametrize("chosen", [{"license": "cc0-1.0"}, {"upload_type": "software"},
                                    {"access_right": "closed"}])
def test_a_licence_kind_or_access_chosen_alone_is_kept(tmp_path, unsteady, chosen):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.meta[7].update(chosen)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert all(unsteady.meta[7][key] == value for key, value in chosen.items())


def test_a_draft_deleted_during_the_share_is_not_said_kept(tmp_path, unsteady, monkeypatch,
                                                           capsys):
    """Review of the fifteenth commit: a draft the person deleted while its
    archive uploaded was said kept, with the new archive."""
    study = _study(tmp_path / "study")
    real = zenodo._call

    def deleting(method, url, **kwargs):
        answer = real(method, url, **kwargs)
        if kwargs.get("upload") is not None:
            unsteady.gone = {7}
        return answer

    monkeypatch.setattr(zenodo, "_call", deleting)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "was deleted on zenodo-sandbox during this share" in said and "is kept" not in said
    monkeypatch.setattr(zenodo, "_call", real)
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert any(c[0] == "POST" and c[1] == "/api/deposit/depositions" for c in unsteady.calls)


def test_a_new_version_begun_is_read_whole_before_it_is_said_changed(tmp_path, unsteady,
                                                                    monkeypatch, capsys):
    """Review of the fifteenth commit: a draft of a new version, listed
    without its files, was said changed, to be published."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.published_ids = {7}
    unsteady.files_of = {n: [{"id": "a", "filename": "a.zip", "checksum": "md5:1"}]
                         for n in (7, 8)}
    unsteady.meta[8] = {k: v for k, v in unsteady.meta[7].items() if k != "version"}
    unsteady.versions_begun = [8]
    real = zenodo._concept_drafts

    def listed(*args, **kwargs):
        return [{k: v for k, v in entry.items() if k != "files"}
                for entry in real(*args, **kwargs) or []]

    monkeypatch.setattr(zenodo, "_concept_drafts", listed)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 2
    assert "delete it there (it holds the record as it was)" in capsys.readouterr().err


def test_a_file_the_person_names_as_the_record_s_archive_is_not_taken_away(tmp_path,
                                                                           unsteady,
                                                                           capsys):
    """Review of the fifteenth commit: a file the person put on the study's
    draft of a new version, named as the record's earlier archive, was
    deleted without a word."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.published_ids = {7}
    unsteady.files_of = {n: [{"id": "a", "filename": "a.zip", "checksum": "md5:1"}]
                         for n in (7, 8)}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 0
    sent = hashlib.md5((tmp_path / "b.zip").read_bytes()).hexdigest()  # noqa: S324
    unsteady.files_of[8] = [{"id": "b", "filename": "b.zip", "checksum": sent},
                            {"id": "m", "filename": "a.zip", "checksum": "md5:9"}]
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    assert "holds a.zip, which this study did not put there" in capsys.readouterr().err
    assert not any(c[0] == "DELETE" for c in unsteady.calls)


def test_a_draft_found_by_its_title_with_another_licence_zenodo_gave_is_the_study_s(
        tmp_path, unsteady):
    """Hands-on review of the fifteenth commit: a draft whose answer was
    lost, given a licence and kind by Zenodo other than the study's, was
    never described."""
    unsteady.meta = {}
    study, out = _studied(tmp_path, "study")
    real = _Unsteady._list

    def listed(self):
        for draft in type(self).drafts:
            type(self).meta.setdefault(draft["id"], {"title": draft["title"],
                                                     "license": "cc-nc"})
        return real(self)

    _Unsteady._list = listed
    try:
        unsteady.faults = {"POST": ["made"]}
        assert _share_to_zenodo(study, out) == 0
        assert unsteady.meta[7]["title"] == "The 1ABC study"
    finally:
        _Unsteady._list = real


def test_a_draft_deleted_during_its_upload_is_said_gone(tmp_path, unsteady, capsys):
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT file": [404]}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "is gone from zenodo-sandbox (deleted there during this share)" in said
    assert "is kept" not in said


def test_an_empty_answer_to_making_a_draft_says_it_may_be_made(tmp_path, unsteady,
                                                               monkeypatch, capsys):
    study = _study(tmp_path / "study")
    real = zenodo._call

    def empty(method, url, **kwargs):
        if method == "POST" and url.endswith("/api/deposit/depositions"):
            real(method, url, **kwargs)
            return None
        return real(method, url, **kwargs)

    monkeypatch.setattr(zenodo, "_call", empty)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "sharing again looks for it by its title" in capsys.readouterr().err
    monkeypatch.setattr(zenodo, "_call", real)
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert not any(c[0] == "POST" for c in unsteady.calls)


def _empty_answer_to(monkeypatch, method: str, ending: str, *, after: bool = True):
    """Zenodo, asked ``method`` of a URL ending in ``ending``, does it (where
    ``after``) and answers 200 with nothing; the steady ``_call`` returned."""
    real = zenodo._call

    def empty(m, url, **kwargs):
        if m == method and url.endswith(ending):
            if after:
                real(m, url, **kwargs)
            return None
        return real(m, url, **kwargs)

    monkeypatch.setattr(zenodo, "_call", empty)
    return real


def test_an_empty_answer_to_the_description_leaves_it_the_study_s(tmp_path, unsteady,
                                                                  monkeypatch, capsys):
    """Review of the seventeenth commit: one empty answer to the study's
    description left it the person's for good, never written again."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _analysed_again(study, "Words of the second share.")
    real = _empty_answer_to(monkeypatch, "PUT", "/api/deposit/depositions/7")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    monkeypatch.setattr(zenodo, "_call", real)
    report = study / "report" / "report.md"
    report.write_text(report.read_text().replace("Words of the second share.",
                                                 "Words of the third share."))
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "changed there" not in capsys.readouterr().out
    assert "Words of the third share." in unsteady.meta[7]["description"]



def test_an_empty_answer_for_a_new_version_begun_is_a_refusal(tmp_path, unsteady,
                                                             monkeypatch, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    _empty_answer_to(monkeypatch, "GET", "/api/deposit/depositions/8", after=False)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "Zenodo began the new version of record 7 (8)" in capsys.readouterr().err


def test_a_kind_given_a_draft_found_by_its_title_is_the_person_s(tmp_path, unsteady):
    """Hands-on review of the seventeenth commit: a kind the person gave a
    draft whose answer was lost, before the share that found it, was taken
    as Zenodo's and written over."""
    unsteady.meta = {}
    study, out = _studied(tmp_path, "study")
    unsteady.late = len(zenodo._LOOKS_S)
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, out) == 2
    unsteady.faults = {}
    unsteady.meta[7] = {"title": unsteady.drafts[0]["title"], "upload_type": "software"}
    assert _share_to_zenodo(study, out) == 0
    assert unsteady.meta[7]["upload_type"] == "software"


def test_an_empty_answer_for_the_latest_version_is_not_taken_for_none(tmp_path, unsteady,
                                                                     monkeypatch, capsys):
    """Hands-on review of the seventeenth commit: an empty answer for the
    latest version was taken as no newer one, and a new version asked of an
    older version."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.latest = {"id": 9, "metadata": {"version": "10"}, "files": []}
    _empty_answer_to(monkeypatch, "GET", "/api/records/latest", after=False)
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "named no record" in capsys.readouterr().err
    assert not any(c[1].endswith("/actions/newversion") for c in unsteady.calls)


def test_a_licence_alone_on_a_draft_made_and_answered_is_the_person_s(tmp_path, unsteady):
    """Review of the seventeenth commit: a licence the person set alone on a
    draft whose making was answered, to one Zenodo gives others, was written
    over."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.meta[7]["license"] = "cc-zero"
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[7]["license"] == "cc-zero"


def test_a_doi_reserved_for_a_draft_gone_or_unanswered_is_said(tmp_path, unsteady,
                                                              monkeypatch, capsys):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.meta[7]["prereserve_doi"] = {}
    unsteady.faults = {"PUT metadata": [404]}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "is gone from zenodo-sandbox" in capsys.readouterr().err
    unsteady.faults = {}
    _empty_answer_to(monkeypatch, "PUT", "/api/deposit/depositions/7", after=False)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "Zenodo reserved no DOI" in capsys.readouterr().err


def test_the_person_s_file_named_as_the_study_s_archive_is_not_taken_away(tmp_path,
                                                                          unsteady,
                                                                          capsys):
    """Review of the seventeenth commit: a new version's upload failing, the
    person put their own a.zip in place of the record's, and the next share
    deleted it without a word."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    sent = hashlib.md5((tmp_path / "a.zip").read_bytes()).hexdigest()  # noqa: S324
    unsteady.published_ids = {7}
    unsteady.files_of = {n: [{"id": "a", "filename": "a.zip", "checksum": f"md5:{sent}"}]
                         for n in (7, 8)}
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    unsteady.faults = {}
    unsteady.files_of[8] = [{"id": "m", "filename": "a.zip", "checksum": "md5:9"}]
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    assert "holds a.zip, which this study did not put there" in capsys.readouterr().err
    assert not any(c[0] == "DELETE" for c in unsteady.calls)


def test_the_listing_a_share_decides_on_may_take_its_time(monkeypatch):
    asked = []
    monkeypatch.setattr(zenodo, "_call", lambda *a, **k: asked.append(k) or [])
    zenodo._marked("http://z/api/deposit/depositions", "t", "x", None, tries=zenodo.TRIES)
    zenodo._concept_drafts("http://z/api/deposit/depositions", "t",
                           {"conceptrecid": "c7", "id": 7}, None)
    assert [k.get("timeout") for k in asked] == [None, None]


def _md5_of(path: Path) -> str:
    return hashlib.md5(path.read_bytes()).hexdigest()  # noqa: S324


def test_a_new_version_shared_again_holds_its_own_upload_as_the_study_s(tmp_path,
                                                                       unsteady):
    """Review of the eighteenth commit: a new version shared again before it
    was published took the archive the study had just put there, of the
    record's name, for the person's."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    out = tmp_path / "shared.zip"
    assert _share_to_zenodo(study, out) == 0
    unsteady.published_ids = {7}
    unsteady.files_of = {n: [{"id": "s", "filename": "shared.zip",
                              "checksum": _md5_of(out)}] for n in (7, 8)}
    _analysed_again(study, "Words of version two.")
    assert _share_to_zenodo(study, out) == 0
    unsteady.files_of[8] = [{"id": "s2", "filename": "shared.zip", "checksum": _md5_of(out)}]
    report = study / "report" / "report.md"
    report.write_text(report.read_text().replace("Words of version two.", "Words again."))
    assert _share_to_zenodo(study, out) == 0


def test_an_archive_kept_after_an_upload_failed_is_still_the_study_s(tmp_path, unsteady):
    """Review of the eighteenth commit: an upload that failed under the same
    name left the archive the study had put there before taken for the
    person's at every later share."""
    study = _study(tmp_path / "study")
    out = tmp_path / "shared.zip"
    assert _share_to_zenodo(study, out) == 0
    unsteady.files = [{"id": "s", "filename": "shared.zip", "checksum": _md5_of(out)}]
    _analysed_again(study, "Words of the second share.")
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, out) == 2
    unsteady.faults = {}
    assert _share_to_zenodo(study, out) == 0


def test_a_draft_given_zenodo_s_defaults_after_its_answer_is_described(tmp_path, unsteady,
                                                                       monkeypatch):
    """Review of the eighteenth commit: a draft Zenodo gave its documented
    licence only after its answer to making it was taken for the person's,
    and never described."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    real = zenodo._call

    def defaults(method, url, **kwargs):
        answer = real(method, url, **kwargs)
        if method == "POST" and url.endswith("/api/deposit/depositions"):
            unsteady.meta[7].update(license="cc-by", access_right="open")
        return answer

    monkeypatch.setattr(zenodo, "_call", defaults)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[7]["title"] == "The 1ABC study"
    assert unsteady.meta[7]["description"]


def test_an_empty_answer_for_the_latest_version_says_to_share_again(tmp_path, unsteady,
                                                                    monkeypatch, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.submitted = True
    unsteady.latest = {"id": 9, "metadata": {"version": "10"}, "files": []}
    _empty_answer_to(monkeypatch, "GET", "/api/records/latest", after=False)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "share again once Zenodo answers" in capsys.readouterr().err


def test_a_zip_of_the_person_s_named_as_the_new_archive_is_not_uploaded_over(tmp_path,
                                                                            unsteady,
                                                                            capsys):
    """Hands-on review of the eighteenth commit: a zip of the person's on
    the draft, the study's archive renamed to its name, was uploaded over
    and lost."""
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 0
    unsteady.files = [{"id": "c", "filename": "c.zip"},
                      {"id": "m", "filename": "c-final.zip", "checksum": "md5:7777"}]
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "c-final.zip") == 2
    assert "holds c-final.zip, which this study did not put there" in capsys.readouterr().err
    assert not any(c[0] in ("PUT", "DELETE") and "/files/" in c[1] for c in unsteady.calls)


def test_an_empty_answer_for_a_record_is_not_said_to_be_none(tmp_path, unsteady,
                                                              monkeypatch):
    _empty_answer_to(monkeypatch, "GET", "/api/records/7", after=False)
    with pytest.raises(ShareRefused) as refused:
        open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert "has no published record" not in str(refused.value)
    assert "open it again later" in str(refused.value)


def test_a_share_failing_after_its_description_is_written_says_so(tmp_path, unsteady,
                                                                  monkeypatch, capsys):
    """Hands-on review of the eighteenth commit: a disk filling once Zenodo
    had answered the description said it not yet written."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    real_call, real_keep, described = zenodo._call, zenodo._keep, []

    def call(method, url, **kwargs):
        answer = real_call(method, url, **kwargs)
        if method == "PUT" and (kwargs.get("body") or {}).get("metadata", {}) \
                .get("description"):
            described.append(url)
        return answer

    def keep(where, what):
        if described:
            raise ShareRefused("No room.", code="environment.share.too_large")
        return real_keep(where, what)

    monkeypatch.setattr(zenodo, "_call", call)
    monkeypatch.setattr(zenodo, "_keep", keep)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "with the new archive and its description" in said
    assert "not yet written" not in said


def test_a_share_failing_before_the_earlier_archive_goes_warns_of_two(tmp_path, unsteady,
                                                                      monkeypatch, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.files = [{"id": "a", "filename": "a.zip"}]
    real = zenodo._keep
    kept: list = []

    def fills(where, what):
        kept.append(what)
        if what.get("zenodo-sandbox", {}).get("stale") and "uploading" not in \
                what["zenodo-sandbox"]:
            raise ShareRefused("No room.", code="environment.share.too_large")
        return real(where, what)

    monkeypatch.setattr(zenodo, "_keep", fills)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    said = capsys.readouterr().err
    assert "the earlier a.zip beside it, not yet taken away" in said


def test_a_zip_of_the_person_s_on_a_new_version_named_as_an_old_archive_stays(tmp_path,
                                                                             unsteady,
                                                                             capsys):
    """Review of the twentieth commit: a name the study once gave an archive,
    the record no longer holding it, was taken for the study's on a new
    version, and the person's file of that name deleted."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 0
    unsteady.files = [{"id": "b", "filename": "b.zip"}]
    unsteady.faults = {"DELETE": [500] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 2
    unsteady.faults, unsteady.files = {}, None
    # The person deleted b.zip, as advised, and published.
    unsteady.published_ids = {7}
    unsteady.files_of = {n: [{"id": "a", "filename": "a.zip"}] for n in (7, 8)}
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    # On the draft of version 2, the person's own b.zip.
    unsteady.files_of[8] = [{"id": "a", "filename": "a.zip"},
                            {"id": "p", "filename": "b.zip", "checksum": "md5:5"}]
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 2
    assert "holds b.zip, which this study did not put there" in capsys.readouterr().err
    assert not any(c[0] == "DELETE" for c in unsteady.calls)


def test_an_empty_answer_for_a_draft_found_by_its_title_makes_no_second(tmp_path, unsteady,
                                                                       monkeypatch):
    study = _study(tmp_path / "study")
    real = _empty_answer_to(monkeypatch, "POST", "/api/deposit/depositions")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    monkeypatch.setattr(zenodo, "_call", real)
    _empty_answer_to(monkeypatch, "GET", "/api/deposit/depositions/7", after=False)
    unsteady.calls.clear()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert not any(c[0] == "POST" for c in unsteady.calls)
    assert len(unsteady.drafts) == 1


def test_a_record_shared_by_an_earlier_release_is_described_again(tmp_path, unsteady):
    """Review of the twentieth commit: a draft shared by 2.6.0, which kept
    only its id and version, was never described again."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    # As 2.6.0 kept it.
    (study / "shared_to.json").write_text(json.dumps(
        {"zenodo-sandbox": {"id": "7", "version": "1", "doi": "10.5072/zenodo.7"}}))
    _analysed_again(study, "Words of the next share.")
    # Its summary not the study's now, which may be the person's: kept, and
    # said; cleared there, the study's.
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Words of the next share." not in unsteady.meta[7]["description"]
    unsteady.meta[7]["description"] = ""
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Words of the next share." in unsteady.meta[7]["description"]


def test_an_earlier_release_s_archive_is_said_as_perhaps_its(tmp_path, unsteady, capsys):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "old.zip") == 0
    (study / "shared_to.json").write_text(json.dumps(
        {"zenodo-sandbox": {"id": "7", "version": "1", "doi": "10.5072/zenodo.7"}}))
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "new.zip") == 2
    assert "or an earlier release of FastMDXplora put there" in capsys.readouterr().err


def test_an_archive_zenodo_took_away_is_not_said_perhaps_there(tmp_path, unsteady,
                                                              monkeypatch, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.files = [{"id": "a", "filename": "a.zip"}]
    real = zenodo._keep

    def fills(where, what):
        if what.get("zenodo-sandbox", {}).get("stale") == []:
            raise ShareRefused("No room.", code="environment.share.too_large")
        return real(where, what)

    monkeypatch.setattr(zenodo, "_keep", fills)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    said = capsys.readouterr().err
    assert "No room." in said and "may not have taken away" not in said


def test_a_draft_not_read_again_after_its_upload_is_not_described(tmp_path, unsteady,
                                                                  monkeypatch, capsys):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    _empty_answer_to(monkeypatch, "GET", "/api/deposit/depositions/7", after=False)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "reading the draft" in capsys.readouterr().err
    assert not unsteady.meta[7].get("description")


def test_a_draft_gone_since_it_was_listed_is_asked_for_again(tmp_path, unsteady, capsys):
    study, out = _studied(tmp_path, "study")
    unsteady.late = len(zenodo._LOOKS_S)
    unsteady.faults = {"POST": ["made late"]}
    assert _share_to_zenodo(study, out) == 2
    unsteady.faults, unsteady.gone = {}, {7}
    capsys.readouterr()
    assert _share_to_zenodo(study, out) == 0
    assert "is gone from Zenodo since it was listed" in capsys.readouterr().out


def _during_upload(monkeypatch, then) -> None:
    """``then()`` done on Zenodo's site as the study's archive uploads."""
    real = zenodo._call

    def call(method, url, **kwargs):
        answer = real(method, url, **kwargs)
        if kwargs.get("upload") is not None:
            then()
        return answer

    monkeypatch.setattr(zenodo, "_call", call)


def test_a_zip_put_over_the_study_s_during_its_upload_is_said(tmp_path, unsteady,
                                                             monkeypatch, capsys):
    """Hands-on review of the twentieth commit: a zip of the person's put on
    the draft under the archive's name as it uploaded was taken as the
    study's, and the share said to publish it."""
    study = _study(tmp_path / "study")
    _during_upload(monkeypatch, lambda: setattr(unsteady, "files", [
        {"id": "p", "filename": "shared.zip", "checksum": "md5:feed"}]))
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "holds shared.zip with other contents than this study sent" in said
    assert "is kept" not in said


def test_another_zip_put_on_the_draft_during_its_upload_is_said(tmp_path, unsteady,
                                                               monkeypatch, capsys):
    study = _study(tmp_path / "study")
    _during_upload(monkeypatch, lambda: setattr(unsteady, "files", [
        {"id": "s", "filename": "shared.zip"}, {"id": "n", "filename": "next.zip"}]))
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "holds next.zip, put there on Zenodo during this share" in capsys.readouterr().err


def test_every_earlier_archive_left_is_named(tmp_path, unsteady, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.files = [{"id": "a", "filename": "a.zip"}]
    unsteady.faults = {"DELETE": [500] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    unsteady.files = [{"id": "a", "filename": "a.zip"}, {"id": "b", "filename": "b.zip"}]
    unsteady.faults = {"DELETE": [500] * zenodo.TRIES}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "c.zip") == 2
    assert "the earlier a.zip, b.zip, which Zenodo may not" in capsys.readouterr().err


def test_a_lone_zip_not_the_study_s_is_not_said_one_of_two(tmp_path, unsteady, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.files = [{"id": "m", "filename": "a.zip", "checksum": "md5:9"}]
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 2
    said = capsys.readouterr().err
    assert "holds a.zip, which this study did not put there" in said
    assert "two archives" not in said


def _as_2_6_0_kept_it(study: Path) -> None:
    (study / "shared_to.json").write_text(json.dumps(
        {"zenodo-sandbox": {"id": "7", "version": "1", "doi": "10.5072/zenodo.7"}}))


def test_an_earlier_release_s_draft_is_described_after_a_share_that_failed(tmp_path,
                                                                           unsteady):
    """Review of the twenty-second commit: a draft 2.6.0 shared, one share
    failing after its upload, was taken for the person's for good."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _analysed_again(study, "Words of the next share.")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Words of the next share." in unsteady.meta[7]["description"]


def test_words_the_person_added_to_an_earlier_release_s_draft_stay(tmp_path, unsteady):
    """Review of the twenty-second commit: a paragraph the person added to a
    draft 2.6.0 described was written over."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    unsteady.meta[7]["description"] += "<p>Funded by the person's grant.</p>"
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Funded by the person's grant." in unsteady.meta[7]["description"]


def test_an_earlier_release_s_draft_keeps_its_names_and_is_described(tmp_path, unsteady):
    """Hands-on review of the twenty-second commit: a draft 2.6.0 shared with
    other names than this share's was never described again."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    names = [{"name": "Aina, Adekunle"}, {"name": "Other, Author"}]
    unsteady.meta[7]["creators"] = names
    unsteady.meta[7]["description"] = ""
    _analysed_again(study, "Words of the next share.")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Words of the next share." in unsteady.meta[7]["description"]
    assert unsteady.meta[7]["creators"] == names
    _analysed_again_too = study / "report" / "report.md"
    _analysed_again_too.write_text(_analysed_again_too.read_text().replace(
        "Words of the next share.", "Words of the share after."))
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Words of the share after." in unsteady.meta[7]["description"]
    assert unsteady.meta[7]["creators"] == names


def test_a_zip_put_back_under_a_replaced_archive_s_name_is_said(tmp_path, unsteady,
                                                               monkeypatch, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "a.zip") == 0
    unsteady.files = [{"id": "a", "filename": "a.zip"}]
    _during_upload(monkeypatch, lambda: setattr(unsteady, "files", [
        {"id": "b", "filename": "b.zip"},
        {"id": "p", "filename": "a.zip", "checksum": "md5:9"}]))
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    assert "holds a.zip, put there on Zenodo during this share" in capsys.readouterr().err


def test_several_zips_not_the_study_s_are_said_as_several(tmp_path, unsteady, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.files = [{"id": "s", "filename": "shared.zip"}, {"id": "x", "filename": "x.zip"},
                      {"id": "y", "filename": "y.zip"}]
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "Delete them there" in capsys.readouterr().err


def test_a_new_draft_after_an_earlier_release_s_is_described_whole(tmp_path, unsteady):
    """Review of the twenty-third commit: a new draft, made once the draft an
    earlier release shared was deleted, was taken for that release's and
    given no names, kind or keywords."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    unsteady.gone = {7}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[8].get("creators") and unsteady.meta[8].get("upload_type")
    assert not json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"].get(
        "by_earlier")


def test_an_earlier_release_s_undescribed_draft_is_described_whole(tmp_path, unsteady):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.meta[7] = {"title": ""}
    _as_2_6_0_kept_it(study)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[7].get("creators") and unsteady.meta[7].get("keywords")


def test_a_zip_of_the_person_s_on_an_earlier_release_s_draft_is_refused(tmp_path,
                                                                         unsteady, capsys):
    """Review of the twenty-third commit: once a draft an earlier release
    shared was shared here, a zip of the person's of the next archive's
    name was uploaded over."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.files = [{"id": "s", "filename": "shared.zip"},
                      {"id": "p", "filename": "b.zip", "checksum": "md5:5"}]
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "b.zip") == 2
    said = capsys.readouterr().err
    assert "holds b.zip, which this study did not put there" in said
    assert "earlier release" not in said
    assert not any(c[0] in ("PUT", "DELETE") and "/files/" in c[1] for c in unsteady.calls)


def test_words_the_person_writes_once_an_earlier_release_s_draft_is_shared_stay(
        tmp_path, unsteady):
    """Review of the twenty-third commit: a summary the person rewrote on a
    draft an earlier release shared was written over at every share."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    mine = unsteady.meta[7]["description"].replace("<p>", "<p>My own summary. ", 1)
    unsteady.meta[7]["description"] = mine
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[7]["description"] == mine


def test_a_summary_the_person_rewrote_on_an_earlier_release_s_draft_stays(tmp_path,
                                                                          unsteady, capsys):
    """Hands-on review of the twenty-third commit: a summary sentence the
    person rewrote on a draft 2.6.0 shared, in its form, was written over."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    mine = unsteady.meta[7]["description"].replace("<p>", "<p>Corrected: ", 1)
    unsteady.meta[7].update(description=mine, title="My own title")
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert unsteady.meta[7]["description"] == mine
    assert unsteady.meta[7]["title"] == "My own title"
    assert "clear the description there and share again" in capsys.readouterr().out


def test_opening_into_new_folders_on_a_full_disk_is_said(tmp_path, monkeypatch):
    import errno

    from fastmdxplora.sharing import opening

    real = Path.mkdir

    def full(self, *args, **kwargs):
        if "deeper" in str(self):
            raise OSError(errno.ENOSPC, "No space left on device")
        return real(self, *args, **kwargs)

    monkeypatch.setattr(Path, "mkdir", full)
    with pytest.raises(ShareRefused) as refused:
        opening.open_shared(str(tmp_path / "x.zip"), into=tmp_path / "new" / "deeper")
    assert refused.value.refusal.code == "environment.share.too_large"
    assert not (tmp_path / "new").exists()


def test_an_answer_for_the_study_s_draft_naming_none_is_said_once(tmp_path, unsteady,
                                                                  monkeypatch, capsys):
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    real = zenodo._call

    def odd(method, url, **kwargs):
        if method == "GET" and url.endswith("/api/deposit/depositions/7"):
            return {"unexpected": True}
        return real(method, url, **kwargs)

    monkeypatch.setattr(zenodo, "_call", odd)
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    said = capsys.readouterr().err
    assert "is no draft or record" in said and "making the draft" not in said
    assert said.count("as it was") == 0


def test_a_record_answer_that_is_no_record_says_to_open_again(tmp_path, unsteady,
                                                             monkeypatch):
    real = zenodo._call

    def odd(method, url, **kwargs):
        if url.endswith("/api/records/7"):
            return {"unexpected": True}
        return real(method, url, **kwargs)

    monkeypatch.setattr(zenodo, "_call", odd)
    with pytest.raises(ShareRefused) as refused:
        open_shared("10.5072/zenodo.7", workspace=tmp_path / "w")
    assert "open it again later" in str(refused.value)


def test_a_new_version_of_an_earlier_release_s_record_stays_marked_as_such(tmp_path,
                                                                          unsteady,
                                                                          capsys):
    """Review of the twenty-fifth commit: a new version of a record 2.6.0
    shared, its first share stopping before the mark was kept, was from
    then on taken for the person's."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.published_ids = {7}
    _as_2_6_0_kept_it(study)
    unsteady.files_of = {n: [{"id": "s", "filename": "shared.zip"}] for n in (7, 8)}
    assert _share_to_zenodo(study, tmp_path / "new.zip") == 2
    assert json.loads((study / "shared_to.json").read_text())["zenodo-sandbox"].get(
        "by_earlier")
    # 2.6.0's zip deleted there, as said.
    unsteady.files_of[8] = []
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "new.zip") == 0
    assert "changed there" not in capsys.readouterr().out


def test_an_earlier_release_s_draft_whose_description_was_not_taken_is_still_the_study_s(
        tmp_path, unsteady, capsys):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    # As 2.6.0 wrote it, naming its own release.
    unsteady.meta[7]["description"] = re.sub(r"FastMDXplora \S+ as", "FastMDXplora 2.6.0 as",
                                             unsteady.meta[7]["description"])
    unsteady.faults = {"PUT metadata": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "changed there" not in capsys.readouterr().out


def test_an_earlier_release_s_archive_stays_its_own_after_an_upload_failed(tmp_path,
                                                                           unsteady):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    unsteady.files = [{"id": "s", "filename": "shared.zip", "checksum": "md5:2600"}]
    unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    unsteady.faults = {}
    unsteady.files = [{"id": "s", "filename": "shared.zip", "checksum": "md5:2600"}]
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0


def test_no_folder_made_to_open_into_is_left_on_a_full_disk(tmp_path, monkeypatch):
    import errno

    from fastmdxplora.sharing import opening

    real = os.mkdir

    def full(path, *args, **kwargs):
        # As the disk fails: its folders above made, the last one not.
        if str(path).endswith("c") and Path(path).parent.exists():
            raise OSError(errno.ENOSPC, "No space left on device")
        return real(path, *args, **kwargs)

    monkeypatch.setattr(os, "mkdir", full)
    with pytest.raises(ShareRefused):
        opening.open_shared(str(tmp_path / "x.zip"), into=tmp_path / "a" / "b" / "c")
    assert not (tmp_path / "a").exists()


def test_a_gateway_s_404_page_is_not_taken_for_a_draft_gone(tmp_path, unsteady,
                                                           monkeypatch, capsys):
    """Hands-on review of the twenty-fifth commit: a gateway's HTML 404 for
    the study's draft was taken as the draft gone, and a second made."""
    import io
    import urllib.error
    import urllib.request

    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    real = urllib.request.urlopen

    def gateway(request, timeout=None):
        if request.get_method() == "GET" and request.full_url.endswith("/depositions/7"):
            raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {},
                                         io.BytesIO(b"<html>Not Found</html>"))
        return real(request, timeout=timeout)

    monkeypatch.setattr(urllib.request, "urlopen", gateway)
    unsteady.calls.clear()
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 2
    assert "share again once Zenodo answers" in capsys.readouterr().err
    assert not any(c[0] == "POST" for c in unsteady.calls)


def test_an_archive_uploaded_after_many_failed_uploads_is_still_the_study_s(tmp_path,
                                                                            unsteady):
    study = _study(tmp_path / "study")
    out = tmp_path / "shared.zip"
    assert _share_to_zenodo(study, out) == 0
    unsteady.files = [{"id": "s", "filename": "shared.zip", "checksum": _md5_of(out)}]
    report = study / "report" / "report.md"
    for n in range(11):
        report.write_text(report.read_text() + f"\nNote {n}.\n")
        unsteady.faults = {"PUT file": [503] * zenodo.TRIES}
        assert _share_to_zenodo(study, out) == 2
    unsteady.faults = {}
    assert _share_to_zenodo(study, out) == 0


def test_names_changed_on_an_earlier_release_s_draft_leave_its_words_the_study_s(
        tmp_path, unsteady):
    """Hands-on review of the twenty-fifth commit: names the person set on a
    draft 2.6.0 shared, as the share advised, made its description theirs."""
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.meta[7]["creators"] = [{"name": "Aina, Adekunle"}, {"name": "New, Author"}]
    _analysed_again(study, "Words of the next share.")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert "Words of the next share." in unsteady.meta[7]["description"]
    assert len(unsteady.meta[7]["creators"]) == 2


def test_a_licence_removed_on_an_earlier_release_s_draft_stays_removed(tmp_path, unsteady,
                                                                       capsys):
    unsteady.meta = {}
    study = _study(tmp_path / "study")
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    _as_2_6_0_kept_it(study)
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    unsteady.meta[7].update(license="", access_right="closed")
    capsys.readouterr()
    assert _share_to_zenodo(study, tmp_path / "shared.zip") == 0
    assert not unsteady.meta[7].get("license")
    assert "licence and access are kept" in capsys.readouterr().out
