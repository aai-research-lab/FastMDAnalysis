"""Fetch real structures from RCSB through the retry, with failures put in.

`tests/test_a_fetch_is_tried_again.py` answers every request from a fake
RCSB. This asks the real one, from the machine a study runs on, so what is
measured is the retry against the service whose dropped connection it was
written for (6B73 on the nightly validation of 2.6.0): that a request which
drops once is made again and the file arrives, that one which drops each
time is given up after three tries and about ten seconds, that a machine
with no route out is told so at once, and that an entry with no PDB file
still comes as mmCIF.

A failure is put in by answering the first requests with it and passing
the rest to RCSB. A clean fetch is made too, several times, to see whether
RCSB drops one on its own.

Usage::

    python scripts/check_rcsb_fetch.py --output runs/rcsb-fetch
    python scripts/check_rcsb_fetch.py --output runs/rcsb-fetch --repeats 20

Writes ``result.json`` and ``result.md`` to the output, and exits 1 where a
check fails.
"""
from __future__ import annotations

import argparse
import json
import logging
import platform
import sys
import time
import urllib.error
import urllib.request
from contextlib import contextmanager
from pathlib import Path

sys.path.insert(0, "src")

from fastmdxplora.refusals import StudyError  # noqa: E402
from fastmdxplora.setup import pipeline  # noqa: E402

PDB_ID = "6B73"
#: An entry deposited as mmCIF only, so the PDB request is answered 404.
#: Such entries are those too large for the PDB format, so the mmCIF that
#: arrives is refused at conversion, as it should be; what is checked is
#: that it arrived.
CIF_ONLY_ID = "5T2A"


def _reset() -> OSError:
    return urllib.error.URLError(ConnectionResetError(104, "Connection reset by peer"))


def _busy(url: str) -> urllib.error.HTTPError:
    return urllib.error.HTTPError(url, 503, "Service Unavailable", {}, None)


@contextmanager
def _rcsb(failures):
    """RCSB, with the first requests answered by ``failures`` (an entry of
    None passes that request on) and every request after by RCSB itself.
    ``failures`` may be a callable, which answers every request.

    Each request is caught where the pipeline makes it: `pipeline._download`
    where there is one (it fetches with a timeout), `urlretrieve` before it."""
    owner, name = ((pipeline, "_download") if hasattr(pipeline, "_download")
                   else (urllib.request, "urlretrieve"))
    real = getattr(owner, name)
    calls: list[str] = []
    waits: list[float] = []
    queue = None if callable(failures) else list(failures)

    def urlretrieve(url, dest, *args, **kwargs):
        calls.append(url)
        failure = failures(url) if queue is None else (queue.pop(0) if queue else None)
        if failure is not None:
            raise failure
        return real(url, dest, *args, **kwargs)

    sleep = pipeline._sleep

    def wait(seconds):
        waits.append(seconds)
        sleep(seconds)

    setattr(owner, name, urlretrieve)
    pipeline._sleep = wait
    try:
        yield calls, waits
    finally:
        setattr(owner, name, real)
        pipeline._sleep = sleep


def _fetch(pdb_id: str, dest: Path):
    """Fetch, and say what came: the file's size, or the refusal."""
    started = time.monotonic()
    try:
        pipeline._fetch_pdb_from_rcsb(pdb_id, dest)
        refusal = None
    except StudyError as exc:
        refusal = exc
    return {
        "seconds": round(time.monotonic() - started, 2),
        "bytes": dest.stat().st_size if dest.exists() else 0,
        "atoms": _atoms(dest),
        "refusal": str(refusal) if refusal else None,
        "code": getattr(refusal, "code", None),
    }


def _atoms(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open(errors="replace") as f:
        return sum(1 for line in f if line.startswith(("ATOM", "HETATM")))


def check_a_reset_once(out: Path) -> dict:
    with _rcsb([_reset()]) as (calls, waits):
        got = _fetch(PDB_ID, out / "reset-once" / "input.pdb")
    return {
        "name": "A connection reset once is asked again and the file comes",
        "passed": got["refusal"] is None and got["atoms"] > 0
        and len(calls) == 2 and waits == [2.0],
        "requests": len(calls), "waits": waits, **got,
    }


def check_busy_once(out: Path) -> dict:
    url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
    with _rcsb([_busy(url)]) as (calls, waits):
        got = _fetch(PDB_ID, out / "busy-once" / "input.pdb")
    return {
        "name": "A 503 once is asked again and the file comes",
        "passed": got["refusal"] is None and got["atoms"] > 0
        and len(calls) == 2 and waits == [2.0],
        "requests": len(calls), "waits": waits, **got,
    }


def check_dropped_each_time(out: Path) -> dict:
    with _rcsb(lambda url: _reset()) as (calls, waits):
        got = _fetch(PDB_ID, out / "dropped" / "input.pdb")
    return {
        "name": "A connection dropped each time is given up after three tries",
        "passed": got["refusal"] is not None and len(calls) == 3
        and waits == [2.0, 8.0] and "in 3 tries" in got["refusal"]
        and got["code"] == "environment.service.unreachable"
        and got["bytes"] == 0 and got["seconds"] < 20,
        "requests": len(calls), "waits": waits, **got,
    }


def check_no_route(out: Path) -> dict:
    # A proxy nothing listens on: the connection is refused, as a machine
    # with no route out refuses it, and is not tried again. urllib keeps the
    # opener it built first, proxies and all, so one is put in its place.
    urllib.request.install_opener(urllib.request.build_opener(
        urllib.request.ProxyHandler({"https": "http://127.0.0.1:9"})))
    try:
        with _rcsb([]) as (calls, waits):
            got = _fetch(PDB_ID, out / "no-route" / "input.pdb")
    finally:
        urllib.request.install_opener(None)
    return {
        "name": "No route out is not asked again and is said plainly",
        "passed": got["refusal"] is not None and len(calls) == 1 and waits == []
        and "no route to the internet" in got["refusal"] and got["seconds"] < 5,
        "requests": len(calls), "waits": waits, **got,
    }


def _cif_arrived(got: dict, dest: Path) -> bool:
    """The mmCIF came whole: written as PDB, or refused as too large for
    PDB only after it was read, which needs the whole file."""
    cif = dest.with_suffix(".cif")
    got["cif_bytes"] = cif.stat().st_size if cif.exists() else 0
    if got["refusal"] is None:
        return got["atoms"] > 0
    return got["code"] == "setup.structure.too_large" and got["cif_bytes"] > 0


def check_cif_only(out: Path) -> dict:
    dest = out / "cif-only" / "input.pdb"
    with _rcsb([]) as (calls, waits):
        got = _fetch(CIF_ONLY_ID, dest)
    return {
        "name": "An entry with no PDB file comes as mmCIF, asked once",
        "passed": _cif_arrived(got, dest)
        and [c.rsplit(".", 1)[-1] for c in calls] == ["pdb", "cif"] and waits == [],
        "requests": len(calls), "waits": waits, **got,
    }


def check_cif_dropped_once(out: Path) -> dict:
    dest = out / "cif-dropped-once" / "input.pdb"
    with _rcsb([None, _reset()]) as (calls, waits):
        got = _fetch(CIF_ONLY_ID, dest)
    return {
        "name": "An mmCIF request that drops once is asked again",
        "passed": _cif_arrived(got, dest)
        and [c.rsplit(".", 1)[-1] for c in calls] == ["pdb", "cif", "cif"]
        and waits == [2.0],
        "requests": len(calls), "waits": waits, **got,
    }


def check_clean_repeats(out: Path, repeats: int) -> dict:
    fetches = []
    with _rcsb([]) as (calls, waits):
        for i in range(repeats):
            fetches.append(_fetch(PDB_ID, out / "clean" / f"{i:02d}" / "input.pdb"))
    retried = len(calls) - repeats
    return {
        "name": f"{PDB_ID} fetched cleanly {repeats} times",
        "passed": all(f["refusal"] is None and f["atoms"] > 0 for f in fetches),
        "requests": len(calls), "retries_seen": retried, "waits": waits,
        "seconds": [f["seconds"] for f in fetches],
        "refusals": [f["refusal"] for f in fetches if f["refusal"]],
    }


def _run(check, *args) -> dict:
    """A check, or what it raised: one that raises is a failure to read,
    not a reason to leave the others unmade."""
    try:
        return check(*args)
    except Exception as exc:  # noqa: BLE001
        logging.exception("%s raised", check.__name__)
        return {"name": check.__name__, "passed": False, "requests": 0,
                "waits": [], "raised": f"{type(exc).__name__}: {exc}"}


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--repeats", type=int, default=10,
                   help="Clean fetches of 6B73 to make (default 10).")
    args = p.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=False)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(name)s %(message)s",
        handlers=[logging.FileHandler(args.output / "fetch.log"),
                  logging.StreamHandler()])

    checks = []
    for check in (check_a_reset_once, check_busy_once, check_dropped_each_time,
                  check_no_route, check_cif_only, check_cif_dropped_once):
        checks.append(_run(check, args.output))
    checks.append(_run(check_clean_repeats, args.output, args.repeats))
    import subprocess
    record = {
        "host": platform.node(),
        "commit": subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                                 text=True).stdout.strip() or None,
        "python": platform.python_version(),
        "when": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "passed": all(c["passed"] for c in checks),
        "checks": checks,
    }
    (args.output / "result.json").write_text(json.dumps(record, indent=2) + "\n")
    lines = [f"# RCSB fetch on {record['host']}", "",
             f"Commit {record['commit']}, Python {record['python']}, "
             f"{record['when']}", "",
             "| Check | Result | Requests | Waits (s) |", "|---|---|---|---|"]
    for c in checks:
        lines.append(f"| {c['name']} | {'pass' if c['passed'] else '**FAIL**'} | "
                     f"{c['requests']} | {', '.join(f'{w:g}' for w in c['waits']) or '-'} |")
    (args.output / "result.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if record["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
