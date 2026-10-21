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
next share; a draft made twice that way (one Zenodo listed only late) is
named for the person to delete, and one the person made is never taken.
"""

from __future__ import annotations

import errno
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
            why = {401: "the token was refused (make one with the deposit:write scope, and "
                        "give it in ZENODO_TOKEN, or ZENODO_SANDBOX_TOKEN for the sandbox)",
                   403: "the token may not do this (its scopes: deposit:write, and "
                        "deposit:actions for a new version of a published record)"} \
                if token else {401: "it is not open to be read (restricted, or not "
                                    "published yet: ask its owner for access, or download "
                                    "its zip and open the file)",
                               403: "it is not open to be read (restricted, or not "
                                    "published yet: ask its owner for access, or download "
                                    "its zip and open the file)"}
            why = {**why,
                   404: "no such record or draft",
                   429: f"too many requests, still after {attempt} tries"}.get(
                exc.code, f"HTTP {exc.code}" + (f", still after {attempt} tries"
                                                 if exc.code in AGAIN and attempt > 1 else ""))
            refused = ShareRefused(f"Zenodo answered {method} {url}: {why}. {answered}".strip(),
                                   code="environment.service.unreachable", url=url)
            refused.later = exc.code in AGAIN or exc.code >= 500
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
            refused.later = True
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
    found = _call("GET", f"{base}/api/records/{record}", said=said, missing_ok=True)
    if not isinstance(found, dict):
        raise ShareRefused(f"Zenodo has no published record {record}: look at the DOI or "
                           "address given (a draft's reserved DOI opens once it is published).",
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
            f"Zenodo's record {record} holds {len(zips)} zip files"
            + (", so there is no shared study to open." if not zips else
               f" ({', '.join(z['name'] for z in zips)}): download the study's and open "
               "the file."), code="environment.share.not_a_study")
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
    except ShareRefused as exc:
        if not getattr(exc, "later", False):
            # The token refused, say: no waiting mends it.
            raise
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
    try:
        return _call("GET", f"{depositions}/{made[0]}", token=token, said=said)
    except ShareRefused as exc:
        raise ShareRefused(f"{exc} Draft {made[0]} is as it was; share again once Zenodo "
                           "answers, and it is taken then.", code=exc.refusal.code) from exc


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
            if getattr(exc, "later", False):
                raise ShareRefused(f"{exc} Nothing was made; share again once Zenodo answers.",
                                   code=exc.refusal.code) from exc
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
            "for it first and finishes it; one it makes later still is named for you to "
            "delete.",
            code=exc.refusal.code) from exc


#: What Zenodo itself sets in a draft's description, not what anyone wrote
#: there: told apart from both what the study wrote and what the person did.
_SET_BY_ZENODO = frozenset({"doi", "prereserve_doi", "publication_date", "version"})


def _latest_published(found: dict[str, Any], depositions: str, token: str,
                      said: Callable[[str], None] | None) -> dict[str, Any]:
    """The latest published version of the record ``found`` is a version of,
    as a deposition (as its drafts are read): by its `latest` link, which
    names it as a record; ``found`` itself where it names none."""
    link = str((found.get("links") or {}).get("latest") or "")
    if not link:
        return found
    record = _call("GET", link, token=token, missing_ok=True, said=said)
    newest = str((record or {}).get("id") or "") if isinstance(record, dict) else ""
    if not newest or newest == str(found.get("id")):
        return found
    deposition = _call("GET", f"{depositions}/{newest}", token=token, missing_ok=True,
                       said=said)
    return deposition if isinstance(deposition, dict) else found


def _concept_drafts(depositions: str, token: str, record: dict[str, Any],
                    said: Callable[[str], None] | None) -> list[dict[str, Any]] | None:
    """The person's drafts of a new version of ``record`` (the drafts of its
    concept); None where Zenodo did not list them."""
    concept = str(record.get("conceptrecid") or "")
    if not concept:
        return []
    try:
        found = _call("GET", f"{depositions}?status=draft&sort=mostrecent&size=100",
                      token=token, missing_ok=True, said=said, timeout=_LOOK_TIMEOUT_S)
    except ShareRefused as exc:
        if not getattr(exc, "later", False):
            raise
        return None
    if not isinstance(found, list):
        return None
    return [entry for entry in found if isinstance(entry, dict)
            and str(entry.get("conceptrecid") or "") == concept
            and str(entry.get("id")) != str(record.get("id"))]


def _described(metadata: Any) -> dict[str, Any]:
    return {k: v for k, v in (metadata if isinstance(metadata, dict) else {}).items()
            if k not in _SET_BY_ZENODO}


def _written(metadata: Any) -> str:
    """A draft's description as Zenodo holds it, as one digest: what this
    study wrote, told from what the person changed there since."""
    import hashlib

    return hashlib.sha256(json.dumps(_described(metadata), sort_keys=True).encode()) \
        .hexdigest()


def _sent(metadata: Any) -> str:
    """What the study sends as a draft's description, digested as a draft's
    is (:func:`_written`): a draft written whose answer was lost is told by
    it, and one the person changed since is not."""
    return _written(metadata)


def _bumped(label: str) -> str:
    """The version after ``label``, of its own form: its last number one
    more ("1.0" gives "1.1", "v5" gives "v6", "2.1.5" gives "2.1.6");
    "" where it holds no number."""
    found = list(re.finditer(r"\d+", label))
    if not found:
        return ""
    last = found[-1]
    return f"{label[:last.start()]}{int(last.group()) + 1}{label[last.end():]}"


def _sendings(kept: Any) -> list[str]:
    """The digests sent whose answers were lost, as kept (one, or several)."""
    if isinstance(kept, str):
        return [kept] if kept else []
    return [str(k) for k in kept or [] if k]


def _names(entry: dict[str, Any]) -> set[str]:
    return {str(f.get("filename") or f.get("key") or "") for f in entry.get("files") or []
            if isinstance(f, dict)} - {""}


def _changed_since(record: dict[str, Any], draft: dict[str, Any]) -> str:
    """What a draft of a new version holds that the version it was begun
    from does not: said (its title, its creators, its files...), or ""
    where the draft is as Zenodo began it."""
    was, now = _described(record.get("metadata")), _described(draft.get("metadata"))
    said = [f"its {key}" for key in sorted(set(was) | set(now)) if was.get(key) != now.get(key)]
    if _names(draft) - _names(record):
        said.append("its files")
    return ", ".join(said)


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
              system: str, version: str) -> dict[str, Any]:
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


def _keep(study: Path, kept: dict[str, Any]) -> None:
    """``kept`` written as the study's `shared_to.json`; a disk with no room
    for it refuses, before anything more is asked of Zenodo."""
    temporary = study / f".{SHARED_TO}.{os.getpid()}.tmp"
    try:
        # Whole or not at all: emptied, it would forget the record it names.
        temporary.write_text(json.dumps(kept, indent=2) + "\n", encoding="utf-8")
        os.replace(temporary, study / SHARED_TO)
    except OSError as exc:
        temporary.unlink(missing_ok=True)
        if exc.errno not in (errno.ENOSPC, errno.EDQUOT):
            raise
        raise ShareRefused(
            f"{study / SHARED_TO}, where the study keeps the record it is shared as, could "
            f"not be written: the disk has no room ({exc.strerror}). Make room there and "
            "share again.", code="environment.share.too_large",
            path=str(study / SHARED_TO)) from exc


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
    made again loses its old archive only once the new one is written.

    What the person does on Zenodo's site is theirs: a description they
    changed, files they added, a draft of a new version they began are
    never written over or taken; where the study cannot tell, it leaves
    them and says so."""
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
    found_as = _found(study, kept, earlier, site=site, depositions=depositions, token=token,
                      said=said)
    deposition, version = found_as["deposition"], found_as["version"]
    published, versioned = found_as["published"], found_as["versioned"]
    title = ""
    made_here = False
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
            deposition, made_here = _taken(depositions, token, title, made, said), True
    if deposition is None:
        title = title or f"FastMDXplora: sharing {study.name} ({secrets.token_hex(4)})"
        # Kept before it is asked for, so a draft made whose answer is lost
        # is found again at the next share.
        kept[site] = {"asked_as": title, "version": version}
        _keep(study, kept)
        deposition, made_here = _new_draft(depositions, token, title, said), True
    if not isinstance(deposition, dict) or not deposition.get("id"):
        raise ShareRefused("Zenodo made no draft.", code="environment.service.unreachable")
    identifier = str(deposition["id"])
    same = str(earlier.get("id") or "") == identifier
    # What is known of this draft is kept with it from one share to the next.
    kept[site] = {"id": identifier, "version": version,
                  **{k: earlier[k] for k in ("archive", "written", "doi", "record", "asked_as",
                                             "sending", "record_written", "record_sending",
                                             "stale", "uploading")
                     if same and earlier.get(k)}}
    if versioned:
        kept[site]["record"] = published
        if earlier.get("archive"):
            kept[site]["archive"] = earlier["archive"]
        # A new version begins with the record's description: what the study
        # wrote for the record (not for a draft of a version since deleted).
        recorded = earlier.get("record_written") if found_as["begun_again"] \
            else earlier.get("written")
        if recorded:
            kept[site]["written"] = kept[site]["record_written"] = recorded
        # Sent for the record and its answer lost: what it holds is that.
        sent = earlier.get("record_sending") if found_as["begun_again"] \
            else earlier.get("sending")
        if sent:
            kept[site]["sending"] = kept[site]["record_sending"] = sent
    if made_here and title:
        # The title it was asked for under: until the study describes it, a
        # draft still so titled is the study's to describe, and a draft
        # Zenodo lists under it later is one made twice, to be named.
        kept[site]["asked_as"] = title
    _keep(study, kept)
    held = deposition.get("metadata") if isinstance(deposition.get("metadata"), dict) else {}
    inherited = {k: v for k, v in held.items() if k not in ("doi", "prereserve_doi")}
    page = str((deposition.get("links") or {}).get("html") or f"{base}/uploads/{identifier}")
    earlier_files = [entry for entry in deposition.get("files") or []
                     if isinstance(entry, dict) and entry.get("id")]
    doi = str((held.get("prereserve_doi") or {}).get("doi") or "")
    if not doi:
        reserved = _call("PUT", f"{depositions}/{identifier}", token=token,
                         body={"metadata": {**inherited, "prereserve_doi": True}}, said=said)
        doi = str((((reserved or {}).get("metadata") or {}).get("prereserve_doi") or {})
                  .get("doi") or "")
    if doi:
        try:
            packed = pack(f"https://doi.org/{doi}")
        except ShareRefused as exc:
            holds = "with the archive it had" if earlier_files else "with no archive yet"
            raise ShareRefused(
                f"{exc} The draft {page} is kept, {holds}: share again once it can be "
                "written, and that draft is finished, not another made.",
                code=exc.refusal.code) from exc
    kept[site]["doi"] = doi
    _keep(study, kept)
    bucket = str((deposition.get("links") or {}).get("bucket") or "")
    if not bucket:
        raise ShareRefused("Zenodo's draft has no bucket to upload into.",
                           code="environment.service.unreachable")
    # Only the archive this study put there is taken away; another one
    # keeps the share back, as a record of two archives opens as neither.
    ours = str(earlier.get("archive") or "")
    replacing = {str(n) for n in (earlier.get("stale") or []) if n} \
        if isinstance(earlier.get("stale"), list) else {str(earlier.get("stale") or "")} - {""}
    # Sent at a share whose answer was lost: there, perhaps, and the study's.
    if earlier.get("uploading"):
        replacing.add(str(earlier["uploading"]))
    replacing -= {archive.name}
    stale = [entry for entry in earlier_files
             if str(entry.get("filename") or entry.get("key") or "") == ours
             and ours and ours != archive.name]
    stale += [entry for entry in earlier_files
              if str(entry.get("filename") or entry.get("key") or "") in replacing
              and entry not in stale]
    foreign = sorted(name for name in _names(deposition)
                     if name.lower().endswith(".zip")
                     and name not in {archive.name, ours, *replacing})
    if foreign:
        raise ShareRefused(
            f"The draft {page} holds {', '.join(foreign)}, which this study did not put "
            "there; a record of two archives opens as neither. Delete it there, or keep it "
            "in a record of its own, then share again.",
            code="environment.service.unusable_response", url=page)
    # The description: the study's where it wrote it last or none was
    # written yet; else the person's, kept.
    marker = str(kept[site].get("asked_as") or title or "")
    # What the study wrote last, where it is of this draft (or, for a new
    # version, of the record it is begun from).
    last = str(kept[site].get("written") or "")
    # Written, its answer lost: the draft holds what was sent.
    # Each sent since it was last written, its answer lost: the draft holds
    # one of them, whichever landed.
    sending = _sendings(kept[site].get("sending"))
    untouched = (bool(last) and _written(held) == last) \
        or _sent(held) in sending or (
        not last and not versioned and str(held.get("title") or "") in ("", marker)
        and not held.get("description"))
    uploaded, removing = False, ""
    try:
        if said:
            said(f"Uploading {archive.name} ({archive.stat().st_size / 1e6:.1f} MB) to "
                 f"{site}...")
        # The new archive first: one of the same name is replaced by it, and
        # the draft keeps the one it had until the new one is there.
        kept[site]["uploading"] = archive.name
        _keep(study, kept)
        _call("PUT", f"{bucket}/{urllib.parse.quote(archive.name)}", token=token,
              upload=archive, said=said)
        uploaded = True
        kept[site]["archive"] = archive.name
        kept[site].pop("uploading", None)
        # Each still the study's, until it is taken away.
        kept[site]["stale"] = sorted({str(e.get("filename") or e.get("key")) for e in stale})
        _keep(study, kept)
        for entry in stale:
            removing = str(entry.get("filename") or entry.get("key"))
            # Gone already is what was asked: a try that did it, then lost
            # its answer, is asked again.
            _call("DELETE", f"{depositions}/{identifier}/files/{entry['id']}",
                  token=token, missing_ok=True, said=said)
            kept[site]["stale"] = [n for n in kept[site]["stale"] if n != removing]
            _keep(study, kept)
            removing = ""
        kept[site].pop("stale", None)
        if untouched:
            ours_said = _metadata(packed["crate"], authors=authors, license_id=license_id,
                                  system=system, version=version)
            kept[site]["sending"] = [*_sendings(kept[site].get("sending")),
                                     _sent({**inherited, **ours_said})][-5:]
            _keep(study, kept)
            answer = _call("PUT", f"{depositions}/{identifier}", token=token,
                           body={"metadata": {**inherited, **ours_said}}, said=said)
            kept[site]["written"] = _written((answer or {}).get("metadata")
                                             if isinstance(answer, dict) else {})
            kept[site].pop("sending", None)
        else:
            # The person's words, kept; the version is the study's to say.
            if str(held.get("version") or "") != str(version):
                _call("PUT", f"{depositions}/{identifier}", token=token,
                      body={"metadata": {**inherited, "version": str(version)}}, said=said)
            if said:
                said(f"The draft {page} holds a description changed since this study last "
                     "wrote it: it is kept as it is, and only the archive is new. Write it "
                     "there as you want it before you publish.")
        _keep(study, kept)
    except ShareRefused as exc:
        kept_as = (f"with the new archive and the earlier {removing}, which Zenodo may not "
                   "have taken away (a record of two archives opens as neither: share again, "
                   "or delete it there, before you publish)" if removing
                   else "with the new archive, its description "
                   + ("perhaps " if getattr(exc, "unclear", False) else "") + "not yet written"
                   if uploaded and untouched
                   else "with the new archive" if uploaded
                   else "with the archive it had, and perhaps the new one beside it if Zenodo "
                   "kept it before its answer was lost (a record of two archives opens as "
                   "neither: share again before you publish, and the earlier one is taken "
                   "away)" if earlier_files and getattr(exc, "unclear", False)
                   and archive.name not in _names(deposition)
                   else "with the archive it had, or the new one if Zenodo kept it before its "
                   "answer was lost" if earlier_files and getattr(exc, "unclear", False)
                   else "with the archive it had" if earlier_files else "with no archive yet")
        mend = ("share again once there is room" if exc.refusal.code
                == "environment.share.too_large" else "share again once Zenodo answers")
        raise ShareRefused(
            f"{exc} The draft {page} is kept, {kept_as}: {mend}, and that draft is finished, "
            "not another made.", code=exc.refusal.code, url=page) from exc
    twice = str(kept[site].get("asked_as") or "")
    try:
        others = [n for n in (_marked(depositions, token, twice, said) or [])
                  if n != identifier] if twice else []
    except ShareRefused:
        # The share is done; the stray draft is named at the next one.
        others = []
    if others and said:
        said(f"Zenodo also holds draft{'s' if len(others) > 1 else ''} {', '.join(others)} "
             f"titled {twice!r}, made by a share whose answer was lost and listed only late: "
             f"delete {'them' if len(others) > 1 else 'it'} among your uploads "
             f"({base}/me/uploads).")
    if said:
        said(f"Draft {f'of version {version} ' if str(version) != '1' else ''}made: {page}"
             + (f" (DOI {doi}, reserved until you publish)." if doi else ".")
             + " Read it there and press Publish; until then nothing is public.")
    return {"url": page, "doi": doi, "id": identifier, "version": version, **packed}


def _found(study: Path, kept: dict[str, Any], earlier: dict[str, Any], *, site: str,
           depositions: str, token: str,
           said: Callable[[str], None] | None) -> dict[str, Any]:
    """The draft the study was shared as, made again (its `deposition`), or
    a draft of a new version of the record it was published as, begun now;
    None where a new draft is to be made; with the `version` it is, the
    `published` record a new version is of, and whether it is one
    (`versioned`)."""
    out: dict[str, Any] = {"deposition": None, "version": str(earlier.get("version") or "1"),
                           "published": str(earlier.get("record") or ""), "versioned": False,
                           "begun_again": False}
    if not earlier.get("id"):
        return out
    published = out["published"]
    again_from = False
    try:
        found = _call("GET", f"{depositions}/{earlier['id']}", token=token, missing_ok=True,
                      said=said)
        if not isinstance(found, dict) and published and published != str(earlier["id"]):
            again_from = True
            # The draft of a new version deleted: begun again from its record.
            found = _call("GET", f"{depositions}/{published}", token=token, missing_ok=True,
                          said=said)
            if isinstance(found, dict) and found.get("submitted") and said:
                said(f"The draft {earlier['id']} of a new version is gone from {site}; "
                     f"beginning a new version of record {published} again.")
        if not isinstance(found, dict):
            if said:
                said(f"The draft or record {earlier['id']} this study was shared as is gone "
                     f"from {site}; making a new one.")
            out["version"] = "1"
            return out
        if not found.get("submitted"):
            # The draft made last time, never published: made again.
            out["deposition"] = found
            return out
    except ShareRefused as exc:
        again = "; share again once Zenodo answers." if getattr(exc, "later", False) else "."
        raise ShareRefused(f"{exc} The draft or record {earlier['id']} the study was shared "
                           f"as is as it was{again}", code=exc.refusal.code) from exc
    # Published: a new version, begun from the latest one published (which
    # may be one the person published by hand since), and counted on from it.
    try:
        newest = _latest_published(found, depositions, token, said)
        published = str(newest.get("id") or found["id"])
        # The next of the latest version's own label; where it holds no
        # number, the study's own count: the version its deleted draft was
        # to be, or the one after the one it published.
        mine = str(earlier.get("version") or "1")
        version = _bumped(str((newest.get("metadata") or {}).get("version") or "")) or (
            mine if again_from else _bumped(mine) or "2")
        begun = _concept_drafts(depositions, token, newest, said)
        if begun is None:
            raise ShareRefused(
                f"Zenodo did not list your drafts, so whether a new version of record "
                f"{published} is begun already is not known; nothing was changed. Share again "
                "once Zenodo answers.", code="environment.service.unreachable")
        if begun:
            # Zenodo keeps one draft of a new version, and gives it to whoever
            # asks: one there already is the person's, or one a share began
            # and lost the answer to; either way left alone.
            draft_id = str(begun[0].get("id"))
            raise ShareRefused(
                f"Zenodo holds a draft of a new version of record {published} ({draft_id}), "
                "begun on its site or by a share whose answer was lost: publish it, or "
                "delete it there, then share again.",
                code="environment.service.unusable_response",
                url=str((begun[0].get("links") or {}).get("html") or ""))
        newer = _call("POST", f"{depositions}/{published}/actions/newversion", token=token,
                      said=said)
        latest = str(((newer or {}).get("links") or {}).get("latest_draft") or "")
        if not latest:
            raise ShareRefused("Zenodo made no new version of the record.",
                               code="environment.service.unreachable")
        deposition = _call("GET", latest, token=token, said=said)
    except ShareRefused as exc:
        if exc.refusal.code == "environment.service.unusable_response" \
                or "did not list your drafts" in str(exc):
            raise
        began = getattr(exc, "unclear", False)
        again = "; share again once Zenodo answers." if getattr(exc, "later", False) else "."
        raise ShareRefused(
            f"{exc} "
            + (f"Zenodo may have begun the new version of record {published}: the next share "
               "finds its draft and asks you to publish or delete it, as it cannot tell it "
               "from one you began." if began else f"Record {published} is as it was{again}"),
            code=exc.refusal.code) from exc
    changed = _changed_since(newest, deposition)
    if changed:
        # The one draft Zenodo keeps, begun by the person and not listed:
        # given back to this share, and left alone.
        raise ShareRefused(
            f"Zenodo holds a draft of a new version of record {published} "
            f"({latest.rsplit('/', 1)[-1]}) that has been changed ({changed}): publish or "
            "delete it there, then share again.",
            code="environment.service.unusable_response", url=latest)
    out.update(deposition=deposition, version=version, published=published, versioned=True,
               begun_again=again_from)
    return out
