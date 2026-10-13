"""How well a paper's MD studies are read, against a truth made by hand.

The set is twenty papers (``PAPERS``), chosen by two researchers working
apart and agreed by both: open access where it can be, methods stated in
full or with gaps on purpose, force fields this software runs and ones it
does not, results with errors, inputs deposited where they are. Each was
read by two readers apart, into the fields of
:mod:`fastmdxplora.paper.fields`; :func:`reconcile` keeps what the two
agree on as the truth and marks the rest contested, and only the agreed
part is scored. The truth is in ``preregistration/paper-reading/``, with
the registration that fixes, before any AI model's reading of these papers
was seen, what is counted and what is claimed.

What is counted, for each paper the feature reads (:func:`score`):

- **studies**: the truth's studies the reading found, and the reading's
  studies the truth has;
- **each agreed field** of a study found: ``correct`` (the value used is the
  truth's), ``wrong`` (another value used), ``missed`` (the truth states it,
  nothing was used), ``invented`` (the truth says the paper does not state
  it, and a value was used), ``agreed_absent`` (neither);
- **each agreed result**: found with the truth's value, or not.

Run, with an AI model chosen (`fastmdx agent model`) and the internet::

    python -m fastmdxplora.validation.paper_reading --out paper_reading.json

**A paper not read, or not found** (the registration's rule, set down
before check 5's second sitting). A refusal by the AI model's provider
(its account's limit, billing, a refused key, overload or a server error
after the retries) means the paper was not read (``not_read``): it is read
in a later sitting, given that sitting's file with ``--earlier``, and the
claims are not judged until every paper is read. So is a paper not
reached, which the rule does not name since it is no reading either:
another service's refusal, a connection cut, a paper not open where it was
fetched, a file not there or not a paper. FastMDXplora's own failure to read an answer it was given counts the
paper's studies as not found (``not_found``) and the paper is not read
again. A refusal that will refuse every paper after it (the account's
limit, no AI model chosen, no key) stops the sitting there. The file is
written again after each paper, so a sitting stopped keeps what it read.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from datetime import datetime, timezone
from http.client import HTTPException
from pathlib import Path
from typing import Any

from fastmdxplora.paper.fields import BY_NAME, FIELDS

__all__ = ["PAPERS", "TRUTH", "DESCRIPTIONS", "THRESHOLDS", "reconcile", "same_value",
           "match_studies", "score", "claims_met", "status_of", "main"]

#: Where the truth is kept, beside the registration.
TRUTH = Path(__file__).resolve().parents[3] / "preregistration" / "paper-reading" / "truth"

#: The set: key, DOI, how it is reached, its licence, and what it tests.
PAPERS: tuple[dict[str, str], ...] = (
    {"key": "seekr2_2022", "doi": "10.1021/acs.jcim.2c00501", "route": "PMC9277580",
     "licence": "CC BY 4.0", "tests": "trypsin-benzamidine in OpenMM; water model only in the deposit; milestoning refused"},
    {"key": "ansari_2022", "doi": "10.1038/s41467-022-33104-3", "route": "PMC9481606",
     "licence": "CC BY 4.0", "tests": "trypsin-benzamidine by OPES; settings in the SI; temperature in a table"},
    {"key": "plattner_2015", "doi": "10.1038/ncomms8653", "route": "PMC4506540",
     "licence": "CC BY 4.0", "tests": "settings by reference to another paper; runs reused from it"},
    {"key": "buch_2011", "doi": "10.1073/pnas.1103547108", "route": "PMC3121846",
     "licence": "free to read, not open access", "tests": "the not-open path: refused by DOI, read from its PDF"},
    {"key": "lindorff_larsen_2012", "doi": "10.1371/journal.pone.0032131", "route": "PMC3285199",
     "licence": "CC BY 4.0", "tests": "eight force fields per system, refused one by one; simulated tempering"},
    {"key": "pereira_2021", "doi": "10.1371/journal.pone.0247841", "route": "PMC7906464",
     "licence": "CC BY 4.0", "tests": "mutants; methods with gaps; units missing from results"},
    {"key": "ilter_2025", "doi": "10.1021/acs.jcim.4c01964", "route": "PMC11776055",
     "licence": "CC BY 4.0", "tests": "OpenMM and CHARMM36m in full; a predicted starting model; QM/MM"},
    {"key": "fischer_2024", "doi": "10.1021/acs.jctc.3c01106", "route": "PMC10938642",
     "licence": "CC BY 4.0", "tests": "peptide folding; ff19SB with OPC; Berendsen barostat; no ions"},
    {"key": "rieloff_2021", "doi": "10.3390/ijms221810174", "route": "PMC8470740",
     "licence": "CC BY 4.0", "tests": "phosphorylated IDPs; TIP4P-D refused; results with errors"},
    {"key": "pantsar_2018", "doi": "10.1371/journal.pcbi.1006458", "route": "PMC6147662",
     "licence": "CC BY 4.0", "tests": "OPLS-AA refused; replica counts in the SI; deposited trajectories"},
    {"key": "love_2023", "doi": "10.1021/acs.jctc.3c00233", "route": "PMC10339676",
     "licence": "CC BY-NC-ND 4.0", "tests": "48 DNA studies; OL21 and Tumuc1; HMR without its mass"},
    {"key": "montgomery_2024", "doi": "10.1021/acs.jpcb.4c01634", "route": "PMC11301690",
     "licence": "CC BY 4.0", "tests": "umbrella sampling in a bilayer in OpenMM; Drude refused"},
    {"key": "siwy_2017", "doi": "10.1371/journal.pcbi.1005314", "route": "PMC5279813",
     "licence": "CC BY 4.0", "tests": "replica exchange refused; CHARMM22* and OPLS-AA"},
    {"key": "larsen_2021", "doi": "10.1371/journal.pcbi.1008807", "route": "PMC8491906",
     "licence": "CC BY 4.0", "tests": "Martini and coarse-grained FEP refused beside an all-atom study"},
    {"key": "lee_2016", "doi": "10.1021/acs.jctc.5b00935", "route": "PMC4712441",
     "licence": "ACS AuthorChoice", "tests": "bilayers in five programs, OpenMM's among them"},
    {"key": "atlas_2024", "doi": "10.1093/nar/gkad1084", "route": "PMC10767941",
     "licence": "CC BY 4.0", "tests": "one protocol for many systems; GROMACS settings to OpenMM"},
    {"key": "lu_huang_2022", "doi": "10.1038/s42003-022-03818-7", "route": "PMC9388646",
     "licence": "CC BY 4.0", "tests": "OpenMM runs beside bias-exchange metadynamics, refused"},
    {"key": "aldeghi_2016", "doi": "10.1039/c5sc02678d", "route": "PMC4700411",
     "licence": "CC BY 3.0", "tests": "alchemical free energies refused; docked poses"},
    {"key": "moore_2023", "doi": "10.1002/prot.26432", "route": "PMC10092333",
     "licence": "CC BY 4.0", "tests": "six force-field and water pairs; GROMOS and OPLS refused"},
    {"key": "mdcath_2024", "doi": "10.1038/s41597-024-04140-z", "route": "PMC11604666",
     "licence": "CC BY 4.0", "tests": "a dataset at five temperatures; CHARMM22* refused"},
)

#: Relative agreement two numbers need to be one value: a paper's numbers
#: converted to one unit by two readers.
_TOLERANCE = 1e-3


def _reduced(text: Any) -> str:
    from fastmdxplora.paper.quotes import squash

    return squash(str(text))[0]


def same_value(name: str, a: Any, b: Any) -> bool:
    """Whether two readings of one field are one value: numbers within a
    thousandth, names naming the same things (:func:`fastmdxplora.paper.
    mapping.recognized`), descriptions sharing their words."""
    from fastmdxplora.paper.mapping import recognized

    kind = BY_NAME[name].kind
    if kind in ("temperature", "pressure", "time", "timestep", "length", "concentration",
                "mass", "count"):
        try:
            x, y = float(a), float(b)
        except (TypeError, ValueError):
            return False
        return math.isclose(x, y, rel_tol=_TOLERANCE, abs_tol=1e-9)
    if kind == "flag":
        return bool(a) == bool(b)
    if name == "mutations":
        from fastmdxplora.paper.extract import mutation_forms

        def forms(value: Any) -> set[str]:
            items = value if isinstance(value, list) else re.split(r"[,;\s]+", str(value))
            return {sorted(mutation_forms(str(item)))[0] for item in items if str(item).strip()}
        return forms(a) == forms(b)
    if name == "pdb_id":
        def codes(value: Any) -> set[str]:
            return {code.upper() for code in re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", str(value))}
        return codes(a) == codes(b) and bool(codes(a))
    said_a, said_b = recognized(name, str(a)), recognized(name, str(b))
    if said_a is not None and (said_a or said_b):
        return said_a == said_b
    words_a = set(re.findall(r"[a-z0-9]{3,}", str(a).lower()))
    words_b = set(re.findall(r"[a-z0-9]{3,}", str(b).lower()))
    if not words_a or not words_b:
        return _reduced(a) == _reduced(b)
    return len(words_a & words_b) / min(len(words_a), len(words_b)) >= 0.5


def _stated(study: dict[str, Any], name: str) -> dict[str, Any] | None:
    record = (study.get("fields") or {}).get(name)
    return record if isinstance(record, dict) and record.get("status") == "stated" else None


#: The fields that tell studies apart, for matching one reader's to the other's.
_IDENTITY = ("pdb_id", "mutations", "protein_forcefield", "nucleic_forcefield",
             "water_model", "temperature", "method", "membrane", "lipid_forcefield",
             "ligands", "engine", "timestep")


def _likeness(a: dict[str, Any], b: dict[str, Any], penalize: bool = True) -> float:
    score = 0.0
    for name in _IDENTITY:
        ra, rb = _stated(a, name), _stated(b, name)
        if ra and rb:
            if same_value(name, ra.get("value"), rb.get("value")):
                score += 2.0
            elif penalize:
                score -= 2.0
    words_a = set(re.findall(r"[a-z0-9]{2,}", str(a.get("label", "")).lower()))
    words_b = set(re.findall(r"[a-z0-9]{2,}", str(b.get("label", "")).lower()))
    if words_a and words_b:
        score += 4.0 * len(words_a & words_b) / max(len(words_a | words_b), 1)
    return score


def match_studies(first: list[dict[str, Any]], second: list[dict[str, Any]],
                  least: float = 1.0, penalize: bool = True) -> list[tuple[int, int]]:
    """Pairs of indices, one study of each list, the likeliest first, each
    study in one pair at most, none liked less than ``least``. With
    ``penalize``, a value two studies give differently counts against them
    (two readers telling their studies apart); without, only what they
    share counts (a reading paired with the truth, whose wrong values must
    be counted, not make its study unpaired)."""
    scored = sorted(((_likeness(a, b, penalize), i, j) for i, a in enumerate(first)
                     for j, b in enumerate(second)), reverse=True)
    used_i: set[int] = set()
    used_j: set[int] = set()
    pairs = []
    for value, i, j in scored:
        if value < least or i in used_i or j in used_j:
            continue
        used_i.add(i)
        used_j.add(j)
        pairs.append((i, j))
    return sorted(pairs)


def _claims_agree(a: dict[str, Any], b: dict[str, Any]) -> bool:
    try:
        return (a.get("quantity") == b.get("quantity")
                and math.isclose(float(a["value"]), float(b["value"]), rel_tol=_TOLERANCE))
    except (KeyError, TypeError, ValueError):
        return False


def reconcile(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    """Two readers' readings of one paper as one truth: the studies both
    found, each field both state with one value (or both leave unstated)
    agreed, the rest contested; the results both list, agreed. A field not
    in ``fields`` and not ``contested`` is agreed not stated. Quotes are
    left out: the truth keeps values, never the paper's text."""
    studies_a, studies_b = a.get("studies") or [], b.get("studies") or []
    pairs = match_studies(studies_a, studies_b)
    out = []
    agreed_fields = contested_fields = 0
    for i, j in pairs:
        sa, sb = studies_a[i], studies_b[j]
        fields: dict[str, Any] = {}
        contested = []
        for field in FIELDS:
            ra, rb = _stated(sa, field.name), _stated(sb, field.name)
            status_a = ((sa.get("fields") or {}).get(field.name) or {}).get("status", "not_stated")
            status_b = ((sb.get("fields") or {}).get(field.name) or {}).get("status", "not_stated")
            if ra and rb and same_value(field.name, ra.get("value"), rb.get("value")):
                fields[field.name] = {"status": "stated", "value": ra.get("value"),
                                      "unit": ra.get("unit"), "where": ra.get("where")}
                agreed_fields += 1
            elif not ra and not rb:
                # Both say the paper's text does not give it (by reference,
                # in the supporting information, or not at all); kept only
                # where they say it the same way, absence meaning not stated.
                if status_a == status_b and status_a != "not_stated":
                    fields[field.name] = {"status": status_a}
                agreed_fields += 1
            else:
                contested.append(field.name)
                contested_fields += 1
        claims = []
        for claim in sa.get("claims") or []:
            if any(_claims_agree(claim, other) for other in sb.get("claims") or []):
                claims.append({key: claim.get(key) for key in (
                    "quantity", "what", "value", "error", "error_kind", "n", "unit", "where")})
        out.append({"id": f"{sa.get('id')}/{sb.get('id')}",
                    "label": sa.get("label") or sb.get("label"),
                    "fields": fields, "contested": contested, "claims": claims})
    matched_a = {i for i, _ in pairs}
    matched_b = {j for _, j in pairs}
    return {
        "key": a.get("key") or b.get("key"),
        "doi": a.get("doi") or b.get("doi"),
        "pmcid": a.get("pmcid") or b.get("pmcid"),
        "studies": out,
        "only_a": [s.get("label") for i, s in enumerate(studies_a) if i not in matched_a],
        "only_b": [s.get("label") for j, s in enumerate(studies_b) if j not in matched_b],
        "agreement": {"studies_a": len(studies_a), "studies_b": len(studies_b),
                      "studies_matched": len(pairs), "fields_agreed": agreed_fields,
                      "fields_contested": contested_fields},
    }


#: Fields that are descriptions rather than choices: kept in the truth, not
#: scored, since two wordings of one description are not told apart well.
DESCRIPTIONS = ("system", "structure_source", "ligands", "protonation", "minimization",
                "equilibration_restraints", "method_details", "chains")


def _used(record: dict[str, Any] | None) -> bool:
    return isinstance(record, dict) and record.get("status") == "stated"


def score(truth: dict[str, Any], reading: dict[str, Any]) -> dict[str, Any]:
    """A reading of one paper counted against its truth (see the module).

    A study the reading has and the truth does not is counted too: each
    value it uses is ``invented``. A field the truth says the paper gives
    only by reference or in its supporting information is not counted
    (``elsewhere``): a reading given the supporting information may find
    it there, and the truth does not know its value."""
    counts = {"correct": 0, "wrong": 0, "missed": 0, "invented": 0, "agreed_absent": 0,
              "elsewhere": 0}
    by_field: dict[str, dict[str, int]] = {}
    truths, everything = truth.get("studies") or [], reading.get("studies") or []
    # A study of the reading named more like one only one reader found than
    # like any of the truth's is neither the truth's nor made up: set aside.
    one_reader = [str(label) for label in (truth.get("only_a") or []) + (truth.get("only_b") or [])]
    found, unsure = [], 0
    for study in everything:
        label = str(study.get("label") or "")
        nearest_one = max((_label_share(label, other) for other in one_reader), default=0.0)
        nearest_truth = max((_label_share(label, str(t.get("label") or "")) for t in truths),
                            default=0.0)
        if nearest_one >= 0.5 and nearest_one > nearest_truth:
            unsure += 1
        else:
            found.append(study)
    pairs = match_studies(truths, found, penalize=False)
    claims_agreed = claims_found = 0

    def count(name: str, outcome: str) -> None:
        counts[outcome] += 1
        by_field.setdefault(name, dict.fromkeys(counts, 0))[outcome] += 1

    for i, j in pairs:
        expected, got = truths[i], found[j]
        contested = set(expected.get("contested") or [])
        for field in FIELDS:
            name = field.name
            if name in DESCRIPTIONS or name in contested:
                continue
            record = (expected.get("fields") or {}).get(name) or {"status": "not_stated"}
            mine = (got.get("fields") or {}).get(name)
            if record.get("status") in ("by_reference", "in_si"):
                outcome = "elsewhere"
            elif record.get("status") == "stated":
                if not _used(mine):
                    outcome = "missed"
                elif same_value(name, record.get("value"), mine.get("value")):
                    outcome = "correct"
                else:
                    outcome = "wrong"
            else:
                outcome = "invented" if _used(mine) else "agreed_absent"
            count(name, outcome)
        for claim in expected.get("claims") or []:
            claims_agreed += 1
            if any(_claims_agree(claim, other) for other in got.get("claims") or []
                   if other.get("status") == "stated"):
                claims_found += 1
    paired = {j for _, j in pairs}
    for j, got in enumerate(found):
        if j in paired:
            continue
        label = str(got.get("label") or "")
        if any(_label_share(label, other) >= 0.5 for other in one_reader):
            unsure += 1  # as like a study one reader found as any: not counted
            continue
        for field in FIELDS:
            if field.name not in DESCRIPTIONS and _used((got.get("fields") or {}).get(field.name)):
                count(field.name, "invented")
    return {"key": truth.get("key"), "studies_truth": len(truths),
            "studies_read": len(everything),
            "studies_matched": len(pairs), "studies_one_reader_found": unsure,
            "fields": counts, "by_field": by_field,
            "claims_agreed": claims_agreed, "claims_found": claims_found}


def _label_share(a: str, b: str) -> float:
    """The words two studies' labels share, over the words of both."""
    words_a = set(re.findall(r"[a-z0-9]{2,}", a.lower()))
    words_b = set(re.findall(r"[a-z0-9]{2,}", b.lower()))
    return len(words_a & words_b) / len(words_a | words_b) if words_a and words_b else 0.0


#: What is claimed (``preregistration/paper-reading.md``), each a test of
#: the pooled counts.
THRESHOLDS = {"values_used_right": 0.98, "stated_values_found": 0.75, "studies_found": 0.80}


#: A refusal's code that is no reading: a service's status (the provider's
#: above all), or a paper not reached (not open where it was fetched, a
#: file not there or not a paper). The paper was not read.
_NOT_READ = ("environment.service.", "environment.model.", "environment.credentials.",
             "environment.paper.not_open", "environment.paper.unreadable", "environment.path.")

#: A refusal's words that are the provider's account, not this paper: every
#: paper after it would be refused the same way. A limit per minute is not
#: one: it passes.
_ACCOUNT = re.compile(r"\((?:401|402|403)\)|usage limits?|spend(?:ing)? limits?|credit balance|"
                      r"insufficient.quota|billing|insufficient.funds", re.I)


def status_of(result: dict[str, Any]) -> str:
    """What one paper's result is: ``read``, ``not_found`` (FastMDXplora
    could not read it: its studies count as not found) or ``not_read`` (a
    service's refusal: read it again in a later sitting). A result from
    before results said so is told by its code and words: FastMDXplora's
    own failure then carried a service's code, and says "FastMDXplora
    could not read the paper"."""
    if "fields" in result:
        return "read"
    if result.get("status") in ("not_found", "not_read"):
        return str(result["status"])
    said = str(result.get("refused") or "")
    if said.startswith(("FastMDXplora could not read the paper", "FastMDXplora found no MD")):
        return "not_found"
    return "not_read" if str(result.get("code") or "").startswith(_NOT_READ) else "not_found"


def _stops_the_sitting(code: str, said: str) -> bool:
    """Whether a refusal will refuse every paper after it: no AI model
    chosen, no key, or the provider's account (its limit, billing, a
    refused key)."""
    if code in ("environment.model.unset", "environment.credentials.absent"):
        return True
    return code.startswith(_NOT_READ) and bool(_ACCOUNT.search(said))


def claims_met(results: list[dict[str, Any]],
               registered: list[str] | tuple[str, ...] | None = None) -> dict[str, Any]:
    """The registered claims, each with its number and whether it is met,
    from the papers' counts pooled. Not judged (``value`` and ``met``
    None) while a paper of ``registered`` (the registered set, ``PAPERS``,
    by default) is not read, or not in ``results``: the papers are named in
    ``not_read``."""
    total = {key: sum(r.get("fields", {}).get(key, 0) for r in results)
             for key in ("correct", "wrong", "missed", "invented", "agreed_absent", "elsewhere")}
    statuses = {str(r.get("key")): status_of(r) for r in results}
    keys = [paper["key"] for paper in PAPERS] if registered is None else list(registered)
    not_read = [key for key in keys if statuses.get(key, "not_read") == "not_read"]
    # A paper not found counts its studies as not found, but the one that
    # must be refused by its DOI and was not given as a file.
    truth_studies = sum(r.get("studies_truth", 0) for r in results
                        if statuses[str(r.get("key"))] != "not_read"
                        and ("fields" in r or r.get("key") != "buch_2011"))
    read_studies = sum(r.get("studies_read", 0) for r in results if "fields" in r)
    matched = sum(r.get("studies_matched", 0) for r in results if "fields" in r)

    def ratio(top: int, bottom: int) -> float | None:
        return round(top / bottom, 4) if bottom and not not_read else None

    def judged(value: float | None, threshold: float) -> dict[str, Any]:
        return {"value": value, "met": None if value is None else value >= threshold}

    used_right = ratio(total["correct"], total["correct"] + total["wrong"] + total["invented"])
    found = ratio(total["correct"], total["correct"] + total["wrong"] + total["missed"])
    studies = ratio(matched, truth_studies)
    buch = next((r for r in results if r.get("key") == "buch_2011"), None)
    return {
        "values_used_right": judged(used_right, THRESHOLDS["values_used_right"]),
        "stated_values_found": judged(found, THRESHOLDS["stated_values_found"]),
        "studies_found": judged(studies, THRESHOLDS["studies_found"]),
        "not_open_refused": {"value": (buch or {}).get("refused_by_doi"),
                             "met": None if not_read else bool(buch) and buch.get(
                                 "refused_by_doi") == "environment.paper.not_open"
                             and "fields" in buch},
        "studies_read_that_are_the_truth_s": ratio(matched, read_studies),
        "not_found": [key for key, status in statuses.items() if status == "not_found"],
        "not_read": not_read,
        "totals": total,
    }


def _reconcile_folders(first: Path, second: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for path in sorted(first.glob("*.json")):
        other = second / path.name
        if not other.is_file():
            continue
        truth = reconcile(json.loads(path.read_text(encoding="utf-8")),
                          json.loads(other.read_text(encoding="utf-8")))
        (out / path.name).write_text(json.dumps(truth, indent=None, ensure_ascii=False,
                                                separators=(",", ":")) + "\n",
                                     encoding="utf-8")
        agreement = truth["agreement"]
        print(f"{truth['key']}: {agreement['studies_matched']} studies matched of "
              f"{agreement['studies_a']} and {agreement['studies_b']}; "
              f"{agreement['fields_agreed']} fields agreed, "
              f"{agreement['fields_contested']} contested")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m fastmdxplora.validation.paper_reading",
        description="Read the registered papers with the AI model chosen, and count the "
                    "readings against the truth.")
    parser.add_argument("--papers", default="all",
                        help="Which papers, by key, comma-separated; all by default.")
    parser.add_argument("--out", default="paper_reading.json", help="Where the counts go.")
    parser.add_argument("--file", action="append", default=[], metavar="KEY=PATH",
                        help="Read this paper from a file here rather than fetching it, "
                             "as the one not open for programs must be; repeatable.")
    parser.add_argument("--earlier", action="append", default=[], metavar="FILE",
                        help="An earlier sitting's file: its papers read or not found are "
                             "kept and not read again, its papers not read are read now; "
                             "repeatable.")
    parser.add_argument("--reconcile", nargs=3, metavar=("FIRST", "SECOND", "OUT"),
                        help="Make the truth from two readers' folders, and stop.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.reconcile:
        first, second, out = (Path(p) for p in args.reconcile)
        _reconcile_folders(first, second, out)
        return 0
    from fastmdxplora.paper.studies import studies_in
    from fastmdxplora.refusals import CodedError

    keys = {paper["key"] for paper in PAPERS}
    wanted = None if args.papers == "all" else set(args.papers.split(","))
    files: dict[str, str] = {}
    for given in args.file:
        key, equals, path = given.partition("=")
        if not equals or key not in keys or not path:
            print(f"--file {given}: give KEY=PATH, with KEY one of {', '.join(sorted(keys))}.",
                  file=sys.stderr)
            return 2
        files[key] = path
    if wanted and wanted - keys:
        print(f"--papers: no paper {', '.join(sorted(wanted - keys))}.", file=sys.stderr)
        return 2
    kept, sitting, was_not_read = _earlier(args.earlier)
    if kept is None:
        return 2

    def to_read(key: str) -> bool:
        return (not wanted or key in wanted) and (
            key not in kept or status_of(kept[key]) == "not_read")

    if to_read("buch_2011") and "buch_2011" not in files:
        print("buch_2011 is read from its PDF, as registered: give it with "
              "--file buch_2011=PATH, or leave it out with --papers.", file=sys.stderr)
        return 2
    out = Path(args.out)
    # Papers read or not found before, and those not read before that this
    # sitting does not try, kept as they were.
    results = [result for key, result in kept.items() if not to_read(key)]
    stopped: dict[str, Any] | None = None
    for paper in PAPERS:
        if wanted and paper["key"] not in wanted:
            continue
        earlier = kept.get(paper["key"])
        if earlier is not None and status_of(earlier) != "not_read":
            continue  # read, or not found, in an earlier sitting: not read again
        truth = json.loads((TRUTH / f"{paper['key']}.json").read_text(encoding="utf-8"))
        if stopped is not None:
            results.append({"key": paper["key"], "status": "not_read", "sitting": sitting,
                            "refused": f"Not tried: the sitting stopped at {stopped['key']}.",
                            "code": stopped["code"],
                            "studies_truth": len(truth.get("studies") or [])})
            continue
        refused_by_doi = None
        if paper["licence"].startswith("free to read"):
            # Not licensed for programs: fetching it by its DOI must be refused.
            from fastmdxplora.paper.fetch import open_paper
            try:
                open_paper(paper["doi"])
                refused_by_doi = False
            except CodedError as exc:
                refused_by_doi = exc.code
            print(f"{paper['key']}: by its DOI, {refused_by_doi or 'fetched'}")
        try:
            _paper, reading = studies_in(files.get(paper["key"], paper["route"]), use_kept=False)
        except (CodedError, OSError, HTTPException) as exc:
            # A connection cut part way is the network's: not read.
            code = getattr(exc, "code", None) if isinstance(exc, CodedError) else None
            code = code or "environment.service.unreachable"
            result = {"key": paper["key"], "refused": str(exc) or type(exc).__name__,
                      "code": code, "refused_by_doi": refused_by_doi, "sitting": sitting,
                      "studies_truth": len(truth.get("studies") or [])}
            result["status"] = status_of(result)
            results.append(result)
            print(f"{paper['key']}: {result['status'].replace('_', ' ')} ({code}): "
                  f"{result['refused']}")
            if result["status"] == "not_read" and _stops_the_sitting(code, str(exc)):
                stopped = result
                print("Every paper after it would be refused the same way: the sitting "
                      "stops here.")
            _write(out, results, sitting, was_not_read)
            continue
        counted = score(truth, reading)
        counted["model"] = reading.get("model")
        counted["sitting"] = sitting
        if refused_by_doi is not None:
            counted["refused_by_doi"] = refused_by_doi
        results.append(counted)
        fields = counted["fields"]
        print(f"{paper['key']}: {counted['studies_matched']}/{counted['studies_truth']} studies; "
              f"{fields['correct']} correct, {fields['wrong']} wrong, {fields['missed']} missed, "
              f"{fields['invented']} invented")
        _write(out, results, sitting, was_not_read)
    claims = _write(out, results, sitting, was_not_read)
    if claims["not_read"]:
        print(f"Not judged: {len(claims['not_read'])} paper(s) not read: "
              f"{', '.join(claims['not_read'])}. Read them in a later sitting:")
        print(f"  python -m fastmdxplora.validation.paper_reading --earlier {out} "
              f"--out <the next sitting's file>" + "".join(
                  f" --file {key}={path}" for key, path in files.items()))
    else:
        for name in ("values_used_right", "stated_values_found", "studies_found",
                     "not_open_refused"):
            print(f"{name}: {claims[name]['value']} "
                  f"({'met' if claims[name]['met'] else 'not met'})")
    print(f"Totals: {claims['totals']}. Written to {out}.")
    return 0


def _earlier(paths: list[str]) -> tuple[dict[str, dict[str, Any]] | None, int, list[str]]:
    """The papers' results from earlier sittings' files, given in the order
    they ran: of each paper, its first result read or not found, else its
    last not read; this sitting's number; and the papers not read before."""
    kept: dict[str, dict[str, Any]] = {}
    sitting = 1
    for path in paths:
        try:
            earlier = json.loads(Path(path).read_text(encoding="utf-8"))
            results = earlier["results"]
            said = int(earlier.get("sittings") or 1)
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
            print(f"--earlier {path}: not a sitting's file ({exc}).", file=sys.stderr)
            return None, 0, []
        sitting = max(sitting + 1, said + 1)
        for result in results:
            if not isinstance(result, dict) or not result.get("key"):
                continue
            key = str(result["key"])
            # A paper read, or not found, is not read again: its first such
            # result stands, whatever a later file holds.
            if key not in kept or status_of(kept[key]) == "not_read":
                kept[key] = result
    return kept, sitting, [paper["key"] for paper in PAPERS if paper["key"] in kept
                           and status_of(kept[paper["key"]]) == "not_read"]


def _write(out: Path, results: list[dict[str, Any]], sitting: int,
           was_not_read: list[str]) -> dict[str, Any]:
    """Write the sitting's file whole (through a file beside it, so a
    sitting stopped part way leaves the last one written), and give its
    claims."""
    claims = claims_met(results)
    read_again = [key for key in was_not_read
                  if any(r.get("key") == key and r.get("sitting") == sitting
                         and status_of(r) != "not_read" for r in results)]
    body = json.dumps({
        "when": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sittings": sitting, "read_again": read_again,
        "results": results, "totals": claims["totals"], "claims": claims}, indent=1) + "\n"
    part = out.with_name(out.name + ".part")
    part.write_text(body, encoding="utf-8")
    part.replace(out)
    return claims


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
