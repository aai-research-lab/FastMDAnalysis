"""Zenodo, for a shared study: a draft made from it, and a record read back.

`--share-to zenodo` makes a draft, never a publication: Zenodo reserves the
draft's DOI first, the archive is packed naming that DOI as its own, the
archive is uploaded, and the draft's description is filled in from the
study's records. The person reads the draft and presses Publish. A study
shared again is a new version of the record it was shared as, kept in
`shared_to.json`.

Opening a study from a DOI or a Zenodo address reads the record
(:func:`record_of`, :func:`archive_of`) and downloads its archive from the
address Zenodo serves files at, checked against the MD5 Zenodo recorded.

Zenodo's REST API, as documented at developers.zenodo.org: depositions
under `/api/deposit/depositions`, a file put into a deposition's bucket,
records under `/api/records/<id>`. `FASTMDXPLORA_ZENODO_URL` points both
sites elsewhere (a test's own server).

Zenodo limits how often it is asked, in requests a minute (429, with when
it may be asked again in `X-RateLimit-Reset`, or how long to wait in
`Retry-After`) and is at times briefly unable to answer (502, 503, 504):
each request is asked again after a wait, up to :data:`TRIES` times. A
request that makes something (a new draft) is asked again only after a 429,
which Zenodo answers before doing anything. A new draft is asked for with a
title of its own (:func:`_marked`) recorded in the study first: where the
answer is lost after Zenodo may have made it (a 5xx, a connection cut), the
draft of that title is looked for and taken as this one, then and at the
next share, so a draft is never made twice and never one the person made.
"""

from __future__ import annotations

import html
import http.client
import json
import os
import re
import secrets
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import Any

from fastmdxplora.sharing import ShareRefused

__all__ = ["SHARED_TO", "SITES", "TRIES", "archive_of", "base_url", "draft", "record_of",
           "wait_before"]

SITES = {"zenodo": "https://zenodo.org", "zenodo-sandbox": "https://sandbox.zenodo.org"}
TOKENS = {"zenodo": "ZENODO_TOKEN", "zenodo-sandbox": "ZENODO_SANDBOX_TOKEN"}
#: Zenodo's DOI prefixes: 10.5281 for records, 10.5072 for the sandbox's.
_PREFIXES = {"10.5281": "zenodo", "10.5072": "zenodo-sandbox"}
#: Where a study keeps the Zenodo record it was shared as.
SHARED_TO = "shared_to.json"
TIMEOUT_S = 120
UPLOAD_TIMEOUT_S = 900
_AGENT = "FastMDXplora (https://github.com/aai-research-lab/FastMDXplora)"
#: How many times a request is asked before Zenodo's answer is taken.
TRIES = 5
#: Answers that say to ask again: too many requests, or briefly unable.
AGAIN = frozenset({429, 502, 503, 504})
#: The longest wait between two tries, whatever `Retry-After` asks.
_LONGEST_WAIT_S = 60.0
_sleep = time.sleep


def wait_before(attempt: int, answer: Any = None) -> float:
    """How long to wait before try ``attempt + 1`` (``attempt`` counting from
    1), by the answer's headers (or its `Retry-After` alone): the seconds
    `Retry-After` asks, else until the time `X-RateLimit-Reset` gives (Unix
    seconds), else 1, 2, 4, 8 s; never over a minute."""
    headers = {"Retry-After": answer} if isinstance(answer, str) or answer is None \
        else answer

    def number(name: str) -> float | None:
        try:
            found = float(str(headers.get(name) or ""))
        except (TypeError, ValueError):
            return None
        return found if found == found and found >= 0 else None  # NaN is not one

    wait = number("Retry-After")
    if wait is None:
        reset = number("X-RateLimit-Reset")
        if reset is not None and reset > 1e9:
            wait = max(reset - time.time(), 0.0) + 1.0
    if wait is None:
        wait = 2.0 ** (attempt - 1)
    return min(wait, _LONGEST_WAIT_S)


class _Progress:
    """A file read as it is sent, saying how far it has gone at each tenth."""

    def __init__(self, path: Path, said: Callable[[str], None] | None) -> None:
        self.handle = path.open("rb")
        self.size = path.stat().st_size
        self.sent = 0
        self.tenth = 0
        self.said = said

    def read(self, amount: int = -1) -> bytes:
        piece = self.handle.read(amount)
        self.sent += len(piece)
        if self.said and self.size >= 10 << 20:
            tenth = self.sent * 10 // self.size
            if tenth > self.tenth and tenth < 10:
                self.tenth = tenth
                self.said(f"Uploaded {self.sent / 1e6:.0f} of {self.size / 1e6:.0f} MB.")
        return piece

    def close(self) -> None:
        self.handle.close()


def base_url(site: str) -> str:
    return (os.environ.get("FASTMDXPLORA_ZENODO_URL") or SITES[site]).rstrip("/")


def _call(method: str, url: str, *, token: str = "", body: Any = None,
          upload: Path | None = None, missing_ok: bool = False,
          said: Callable[[str], None] | None = None, tries: int = TRIES,
          timeout: float | None = None) -> Any:
    """Zenodo's answer as JSON (or None for an empty one), asked again where
    Zenodo says to (:data:`AGAIN`) or could not be reached; refused, naming
    what Zenodo said, where it does not answer as it documents."""
    headers = {"User-Agent": _AGENT, "Accept": "application/json"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    sent: bytes | None = None
    if body is not None:
        sent = json.dumps(body).encode("utf-8")
        headers["Content-Type"] = "application/json"
    if upload is not None:
        headers["Content-Type"] = "application/octet-stream"
        headers["Content-Length"] = str(upload.stat().st_size)
    # Asking again to make something could make it twice: only once Zenodo
    # has said it did nothing (429).
    makes = method == "POST"
    attempt = 0
    while True:
        attempt += 1
        data: Any = _Progress(upload, said) if upload is not None else sent
        request = urllib.request.Request(url, data=data, method=method, headers=headers)
        try:
            # Zenodo answers an upload once it has kept the file, which for
            # gigabytes takes longer than an answer otherwise does.
            with urllib.request.urlopen(request, timeout=timeout or (
                    TIMEOUT_S if upload is None else UPLOAD_TIMEOUT_S)) as answer:
                text = answer.read().decode("utf-8", errors="replace")
            break
        except urllib.error.HTTPError as exc:
            if missing_ok and exc.code in (404, 410):
                return None
            if exc.code in AGAIN and attempt < tries and (not makes or exc.code == 429):
                _again(f"Zenodo answered {method} {url} with HTTP {exc.code}", attempt,
                       exc.headers, said)
                continue
            try:
                answered = exc.read().decode("utf-8", errors="replace")[:400]
            except (OSError, http.client.HTTPException):
                # The answer's own words cut short: said without them.
                answered = ""
            why = {401: "the token was refused", 403: "the token may not do this (its scopes)",
                   404: "no such record or draft",
                   429: f"too many requests, still after {attempt} tries"}.get(
                exc.code, f"HTTP {exc.code}" + (f", still after {attempt} tries"
                                                 if exc.code in AGAIN and attempt > 1 else ""))
            refused = ShareRefused(f"Zenodo answered {method} {url}: {why}. {answered}".strip(),
                                   code="environment.service.unreachable", url=url)
            # Zenodo, or the gateway before it, gave up part way: what was
            # asked may have been done.
            refused.unclear = exc.code >= 500
            raise refused from exc
        except (urllib.error.URLError, OSError, http.client.HTTPException) as exc:
            if attempt < tries and not makes:
                _again(f"The connection to Zenodo at {url} was lost or not made ({exc})",
                       attempt, None, said)
                continue
            refused = ShareRefused(f"The connection to Zenodo at {url} was lost or not made: "
                                   f"{exc}" + (f", still after {attempt} tries" if attempt > 1
                                               else "") + ".",
                                   code="environment.service.unreachable", url=url)
            # Failing as the request was made or sent (urllib says so as a
            # URLError: refused, no such host or network, TLS, a timeout
            # connecting): nothing was done. Lost once it had gone, waiting
            # for or reading the answer: it may have been.
            refused.unclear = not isinstance(exc, urllib.error.URLError)
            raise refused from exc
        finally:
            if isinstance(data, _Progress):
                data.close()
    if not text.strip():
        return None
    try:
        return json.loads(text)
    except ValueError as exc:
        raise ShareRefused(f"Zenodo's answer to {method} {url} was not JSON.",
                           code="environment.service.unreachable", url=url) from exc


# -- reading a record --------------------------------------------------------

def _again(what: str, attempt: int, answer: Any,
           said: Callable[[str], None] | None) -> None:
    wait = wait_before(attempt, answer)
    if said:
        said(f"{what}; asking again in {wait:.0f} s ({attempt + 1} of {TRIES}).")
    _sleep(wait)


def record_of(source: str) -> tuple[str, str]:
    """The site and record a DOI or a Zenodo address names; refused for any
    other (`environment.share.not_a_study`)."""
    text = source.strip()
    text = re.sub(r"^(?:https?://)?(?:dx\.)?doi\.org/", "", text, flags=re.IGNORECASE)
    text = re.sub(r"^(?:https?://)?(?:www\.|sandbox\.)?zenodo\.org/doi/", "", text,
                  flags=re.IGNORECASE)
    text = re.sub(r"^doi:", "", text, flags=re.IGNORECASE)
    found = re.fullmatch(r"(10\.\d{4,9})/zenodo\.(\d+)", text, flags=re.IGNORECASE)
    if found and found.group(1) in _PREFIXES:
        return _PREFIXES[found.group(1)], found.group(2)
    url = urllib.parse.urlparse(source.strip())
    host = (url.hostname or "").lower().removeprefix("www.")
    here = urllib.parse.urlparse(os.environ.get("FASTMDXPLORA_ZENODO_URL") or "")
    sites = {urllib.parse.urlparse(v).hostname: k for k, v in SITES.items()}
    if here.hostname and host == here.hostname:
        sites[host] = "zenodo"
    # A published record's address; a draft's (/uploads/N) opens only for
    # its owner, and is published first.
    numbered = re.search(r"/records?/(\d+)", url.path or "")
    if host in sites and numbered:
        return sites[host], numbered.group(1)
    raise ShareRefused(
        f"{source!r} is not a file here, a Zenodo DOI (10.5281/zenodo.N) or a Zenodo "
        "address (https://zenodo.org/records/N).", code="environment.share.not_a_study")


def archive_of(site: str, record: str,
               said: Callable[[str], None] | None = None) -> dict[str, Any]:
    """The record's archive: its `name`, `size`, `md5`, the `url` it is
    downloaded from, and the record's `doi`, `id` and `version`."""
    base = base_url(site)
    found = _call("GET", f"{base}/api/records/{record}", said=said)
    if not isinstance(found, dict):
        raise ShareRefused(f"Zenodo has no record {record}.",
                           code="environment.share.not_a_study")
    entries = found.get("files") or []
    if isinstance(entries, dict):
        entries = entries.get("entries") or []
        entries = list(entries.values()) if isinstance(entries, dict) else entries
    zips = []
    for entry in entries if isinstance(entries, list) else []:
        name = str(entry.get("key") or entry.get("filename") or "")
        if name.lower().endswith(".zip"):
            checksum = str(entry.get("checksum") or "")
            zips.append({"name": name, "size": int(entry.get("size") or entry.get("filesize")
                                                    or 0),
                         "md5": checksum.split(":", 1)[-1].lower()})
    if len(zips) != 1:
        raise ShareRefused(
            f"Zenodo's record {record} holds {len(zips)} zip files, "
            + ("so there is no shared study to open" if not zips else
               f"({', '.join(z['name'] for z in zips)}); download the study's and open "
               "the file"), code="environment.share.not_a_study")
    archive = zips[0]
    identifier = str(found.get("id") or record)
    # The address Zenodo serves a record's files at; the API's own links
    # have been found to move.
    archive["url"] = (f"{base}/records/{identifier}/files/"
                      f"{urllib.parse.quote(archive['name'])}?download=1")
    metadata = found.get("metadata") if isinstance(found.get("metadata"), dict) else {}
    archive.update(doi=str(found.get("doi") or metadata.get("doi") or ""), id=identifier,
                   version=str(metadata.get("version") or ""), site=site)
    return archive


# -- making a draft ----------------------------------------------------------

#: How long one look among the person's drafts may take.
_LOOK_TIMEOUT_S = 20.0


def _marked(depositions: str, token: str, title: str,
            said: Callable[[str], None] | None) -> list[str] | None:
    """The ids of the person's drafts titled ``title``, the title a new
    draft is asked for with; None where Zenodo did not list them (one try,
    soon given up: the looks are timed by their caller)."""
    try:
        found = _call("GET", f"{depositions}?status=draft&sort=mostrecent&size=100",
                      token=token, missing_ok=True, said=said, tries=1,
                      timeout=_LOOK_TIMEOUT_S)
    except ShareRefused:
        return None
    if not isinstance(found, list):
        return None
    return [str(entry.get("id")) for entry in found if isinstance(entry, dict)
            and isinstance(entry.get("metadata"), dict)
            and entry["metadata"].get("title") == title]


#: The waits before each look among the person's drafts for one whose answer
#: was lost: Zenodo may finish making it after its gateway has given up.
_LOOKS_S = (0.0, 2.0, 4.0, 8.0, 16.0)


def _taken(depositions: str, token: str, title: str, made: list[str],
           said: Callable[[str], None] | None) -> dict[str, Any]:
    """The one draft titled ``title``; refused where there are several."""
    if len(made) > 1:
        raise ShareRefused(
            f"Zenodo holds {len(made)} drafts titled {title!r} ({', '.join(made)}), made by "
            "shares whose answers were lost. Delete all but one among your uploads, then "
            "share again.", code="environment.service.unusable_response")
    if said:
        said(f"Zenodo made draft {made[0]} though its answer was lost; taking that one.")
    return _call("GET", f"{depositions}/{made[0]}", token=token, said=said)


def _new_draft(depositions: str, token: str, title: str,
               said: Callable[[str], None] | None) -> dict[str, Any]:
    """A new draft titled ``title``, made once: where Zenodo's answer to
    making it is lost and it may have been made, the draft of that title is
    looked for and taken; never one of another title."""
    try:
        return _call("POST", depositions, token=token,
                     body={"metadata": {"prereserve_doi": True, "title": title}}, said=said)
    except ShareRefused as exc:
        if not getattr(exc, "unclear", False):
            raise
        for wait in _LOOKS_S:
            if wait:
                if said:
                    said(f"Zenodo's answer was lost; looking again in {wait:.0f} s for the "
                         "draft it may have made.")
                _sleep(wait)
            made = _marked(depositions, token, title, said)
            if made:
                return _taken(depositions, token, title, made, said)
        raise ShareRefused(
            f"{exc} Zenodo may still make the draft (titled {title!r}): sharing again looks "
            "for it first and finishes it, rather than making another.",
            code=exc.refusal.code) from exc


def _creators(authors: list[tuple[str, str]]) -> list[dict[str, str]]:
    out = []
    for name, orcid in authors:
        parts = name.split()
        said = f"{parts[-1]}, {' '.join(parts[:-1])}" if len(parts) > 1 else name
        creator = {"name": said}
        if orcid:
            creator["orcid"] = orcid
        out.append(creator)
    return out


def _metadata(crate: dict[str, Any], *, authors: list[tuple[str, str]], license_id: str,
              system: str, version: int) -> dict[str, Any]:
    root = next(e for e in crate["@graph"] if e.get("@id") == "./")
    software = next((e for e in crate["@graph"] if e.get("@id") == "#fastmdxplora"), {})
    description = (f"<p>{html.escape(str(root.get('description', '')))}</p>"
                   f"<p>A molecular dynamics study packed by FastMDXplora "
                   f"{software.get('version', '')} as an RO-Crate (ro-crate-metadata.json "
                   "lists every file with its SHA-256). Open it with "
                   "<code>fastmdx gui --open &lt;this record's DOI&gt;</code>.</p>")
    related = [{"identifier": "https://github.com/aai-research-lab/FastMDXplora",
                "relation": "isCompiledBy", "scheme": "url", "resource_type": "software"}]
    if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", system or ""):
        related.append({"identifier": f"10.2210/pdb{system.lower()}/pdb",
                        "relation": "isDerivedFrom", "scheme": "doi"})
    keywords = ["molecular dynamics", "FastMDXplora", "RO-Crate"]
    if system:
        keywords.insert(0, system)
    return {"upload_type": "dataset", "title": str(root.get("name") or ""),
            "description": description, "creators": _creators(authors),
            "access_right": "open", "license": license_id.lower(), "keywords": keywords,
            "related_identifiers": related, "version": str(version),
            "prereserve_doi": True}


def _load(path: Path) -> dict[str, Any]:
    try:
        found = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    return found if isinstance(found, dict) else {}


def draft(study: Path, archive: Path, *, site: str, pack: Callable[[str], dict[str, Any]],
          authors: list[tuple[str, str]], license_id: str, system: str,
          said: Callable[[str], None] | None = None) -> dict[str, Any]:
    """A draft on Zenodo of the study packed by ``pack``, which is given the
    DOI the draft reserved (or "") and returns what :func:`share_study`
    does; the draft's `url`, `doi` and `id`. Never published.

    The study is packed once before Zenodo is asked anything, so a study
    that cannot be shared leaves no draft behind; the draft is recorded in
    `shared_to.json` as soon as it exists, so a failure after leaves a
    draft the next share makes again rather than another one; and a draft
    made again loses its old archive only once the new one is written."""
    token = os.environ.get(TOKENS[site], "")
    if not token:
        raise ShareRefused(
            f"--share-to {site} needs a Zenodo token in {TOKENS[site]}, made at "
            f"{SITES[site]}/account/settings/applications/ with the deposit:write scope.",
            code="environment.service.unreachable")
    packed = pack("")
    base = base_url(site)
    depositions = f"{base}/api/deposit/depositions"
    kept = _load(study / SHARED_TO)
    earlier = kept.get(site) if isinstance(kept.get(site), dict) else {}
    # A draft asked for and never known stays the version it was asked as.
    version = int(earlier.get("version") or (0 if earlier.get("id") else 1)) \
        + (1 if earlier.get("id") else 0)
    deposition: dict[str, Any] | None = None
    if earlier.get("id"):
        versioning = False
        try:
            found = _call("GET", f"{depositions}/{earlier['id']}", token=token,
                          missing_ok=True, said=said)
            if not isinstance(found, dict):
                if said:
                    said(f"The draft or record {earlier['id']} this study was shared as is gone "
                         f"from {site}; making a new one.")
                version = 1
            elif not found.get("submitted"):
                # The draft made last time, never published: made again.
                deposition, version = found, int(earlier.get("version") or 1)
            else:
                versioning = True
                newer = _call("POST", f"{depositions}/{found['id']}/actions/newversion",
                              token=token, said=said)
                latest = str(((newer or {}).get("links") or {}).get("latest_draft") or "")
                if not latest:
                    raise ShareRefused("Zenodo made no new version of the record.",
                                       code="environment.service.unreachable")
                deposition = _call("GET", latest, token=token, said=said)
        except ShareRefused as exc:
            began = versioning and getattr(exc, "unclear", False)
            raise ShareRefused(
                f"{exc} "
                + (f"Zenodo may have begun the new version of record {earlier['id']}: look "
                   "among your uploads for its draft, then share again once Zenodo answers."
                   if began else f"The draft or record {earlier['id']} the study was shared as "
                   "is as it was; share again once Zenodo answers."),
                code=exc.refusal.code) from exc
    title = ""
    if deposition is None and earlier.get("asked_as") and not earlier.get("id"):
        # Asked for at a share whose answer was lost: taken if Zenodo made
        # it, asked for again under the same title if not, so that a draft
        # it makes later is seen beside it rather than lost.
        title = str(earlier["asked_as"])
        made = _marked(depositions, token, title, said)
        if made is None:
            raise ShareRefused(
                f"Zenodo did not list your drafts, so whether it made the one asked for at "
                f"the last share (titled {title!r}) is not known; nothing was made. Share "
                "again once Zenodo answers.", code="environment.service.unreachable")
        if made:
            deposition = _taken(depositions, token, title, made, said)
    if deposition is None:
        title = title or f"FastMDXplora: sharing {study.name} ({secrets.token_hex(4)})"
        # Kept before it is asked for, so a draft made whose answer is lost
        # is found again at the next share.
        kept[site] = {"asked_as": title, "version": version}
        (study / SHARED_TO).write_text(json.dumps(kept, indent=2) + "\n", encoding="utf-8")
        deposition = _new_draft(depositions, token, title, said)
    if not isinstance(deposition, dict) or not deposition.get("id"):
        raise ShareRefused("Zenodo made no draft.", code="environment.service.unreachable")
    identifier = str(deposition["id"])
    kept[site] = {"id": identifier, "version": version}
    (study / SHARED_TO).write_text(json.dumps(kept, indent=2) + "\n", encoding="utf-8")
    inherited = deposition.get("metadata") if isinstance(deposition.get("metadata"), dict) \
        else {}
    inherited = {k: v for k, v in inherited.items() if k not in ("doi", "prereserve_doi")}
    doi = str(((deposition.get("metadata") or {}).get("prereserve_doi") or {}).get("doi") or "")
    if not doi:
        reserved = _call("PUT", f"{depositions}/{identifier}", token=token,
                         body={"metadata": {**inherited, "prereserve_doi": True}}, said=said)
        doi = str((((reserved or {}).get("metadata") or {}).get("prereserve_doi") or {})
                  .get("doi") or "")
    if doi:
        packed = pack(f"https://doi.org/{doi}")
    kept[site]["doi"] = doi
    (study / SHARED_TO).write_text(json.dumps(kept, indent=2) + "\n", encoding="utf-8")
    bucket = str((deposition.get("links") or {}).get("bucket") or "")
    if not bucket:
        raise ShareRefused("Zenodo's draft has no bucket to upload into.",
                           code="environment.service.unreachable")
    page = str((deposition.get("links") or {}).get("html")
               or f"{base}/uploads/{identifier}")
    earlier_files = [entry for entry in deposition.get("files") or []
                     if isinstance(entry, dict) and entry.get("id")]
    uploaded, kept_other = False, ""
    try:
        if said:
            said(f"Uploading {archive.name} ({archive.stat().st_size / 1e6:.1f} MB) to "
                 f"{site}...")
        # The new archive first: one of the same name is replaced by it, and
        # the draft keeps the one it had until the new one is there.
        _call("PUT", f"{bucket}/{urllib.parse.quote(archive.name)}", token=token,
              upload=archive, said=said)
        uploaded = True
        for entry in earlier_files:
            name = str(entry.get("filename") or entry.get("key") or "")
            if name != archive.name:
                kept_other = name or "an earlier file"
                # Gone already is what was asked: a try that did it, then
                # lost its answer, is asked again.
                _call("DELETE", f"{depositions}/{identifier}/files/{entry['id']}",
                      token=token, missing_ok=True, said=said)
                kept_other = ""
        ours = _metadata(packed["crate"], authors=authors, license_id=license_id,
                         system=system, version=version)
        _call("PUT", f"{depositions}/{identifier}", token=token,
              body={"metadata": {**inherited, **ours}}, said=said)
    except ShareRefused as exc:
        kept_as = (f"with the new archive and the earlier {kept_other}, which Zenodo did not "
                   "take away (a record of two archives opens as neither: share again, or "
                   "delete it there, before you publish)" if kept_other
                   else "with the new archive, its description not yet written" if uploaded
                   else "with the archive it had" if earlier_files else "with no archive yet")
        raise ShareRefused(
            f"{exc} The draft {page} is kept, {kept_as}: share again once Zenodo answers, "
            "and that draft is finished, not another made.",
            code=exc.refusal.code, url=page) from exc
    if said:
        said(f"Draft {f'of version {version} ' if version > 1 else ''}made: {page}"
             + (f" (DOI {doi}, reserved until you publish)." if doi else ".")
             + " Read it there and press Publish; until then nothing is public.")
    return {"url": page, "doi": doi, "id": identifier, "version": version, **packed}
