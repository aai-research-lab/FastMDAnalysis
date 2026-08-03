#!/usr/bin/env python
"""Compare setup success under the previous and proposed default force fields.

FastMDXplora 2.2.0 proposes changing the default force field from CHARMM36 to
amber-openff, because the latter can parameterize a bound ligand and the
former cannot. One structure, 4W52, was verified to work. Four tests then
failed, one of them a PDBFixer-prepared peptide, which showed the change is
not free: AMBER's residue templates are stricter than CHARMM's.

This script measures the cost instead of guessing at it. It runs setup only,
which is seconds per structure rather than hours, under both force fields, and
reports every structure whose outcome differs.

    python scripts/compare_forcefield_defaults.py --output-root /tmp/ffcmp

Reading the result:

* Failures confined to degenerate inputs (single residues, fragments) are
  acceptable; fix the fixtures and ship.
* Ordinary proteins failing under AMBER that worked under CHARMM36 mean the
  default should not change yet, however convenient it would be.

A refusal from the heterogen classifier is not a force-field failure and is
reported separately: those are the same under either force field.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

#: A spread chosen to exercise the cases where the two force fields differ:
#: termini, non-standard residues, metals, nucleic acids, and glycans. Adjust
#: freely; the comparison is what matters, not this particular list.
DEFAULT_STRUCTURES = [
    # small, well-behaved proteins
    "1L2Y",   # Trp-cage, 20 residues
    "1UBQ",   # ubiquitin
    "1CRN",   # crambin, three disulfides
    "1AKI",   # hen lysozyme
    "1A2P",   # barnase, three chains
    # the validation pair
    "4W51",   # T4 lysozyme L99A, apo
    "4W52",   # T4 lysozyme L99A with benzene
    # protein-ligand
    "1STP",   # streptavidin with biotin
    "1HVR",   # HIV-1 protease with an inhibitor
    "6LU7",   # SARS-CoV-2 main protease with an inhibitor
    # metals and cofactors
    "1ZNI",   # insulin, zinc
    "4HHB",   # haemoglobin, haem: expected to be refused by the classifier
    # nucleic acid
    "1BNA",   # B-DNA dodecamer
    # glycosylated
    "1IGT",   # antibody, N-linked glycans: expected to be refused
]


@dataclass
class Outcome:
    """What happened to one structure under one force field."""

    ok: bool
    seconds: float
    detail: str = ""

    @property
    def refused(self) -> bool:
        """A deliberate refusal, not a force-field failure."""
        markers = (
            "does not determine what should be simulated",
            "needs dedicated parameters",
            "a sugar can be",
            "covalently bonded",
        )
        return any(m in self.detail for m in markers)


@dataclass
class Row:
    structure: str
    results: dict[str, Outcome] = field(default_factory=dict)

    @property
    def verdict(self) -> str:
        old, new = self.results["charmm36"], self.results["amber-openff"]
        if old.refused or new.refused:
            return "refused"
        if old.ok and new.ok:
            return "both ok"
        if old.ok and not new.ok:
            return "REGRESSION"
        if not old.ok and new.ok:
            return "improved"
        return "both fail"


def run_setup(structure: str, forcefield: str, root: Path) -> Outcome:
    """Run the setup phase alone and report whether it produced a system."""
    out = root / forcefield / structure
    started = time.time()
    process = subprocess.run(
        [
            "fastmdx", "setup",
            "--system", structure,
            "--output", str(out),
            "--forcefield", forcefield,
        ],
        capture_output=True,
        text=True,
        timeout=900,
    )
    elapsed = time.time() - started
    produced = (out / "setup" / "system.xml").is_file()
    detail = ""
    if not produced:
        tail = (process.stderr or process.stdout).strip().splitlines()
        detail = " ".join(tail[-4:])[:400] if tail else "no output"
    return Outcome(ok=produced, seconds=elapsed, detail=detail)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structures", nargs="*", default=None,
                        help="PDB identifiers; defaults to a representative spread")
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--json", type=Path, help="write the full result table here")
    args = parser.parse_args()

    structures = args.structures or DEFAULT_STRUCTURES
    args.output_root.mkdir(parents=True, exist_ok=True)

    rows: list[Row] = []
    print(f"Comparing setup under charmm36 and amber-openff, "
          f"{len(structures)} structures\n")
    for structure in structures:
        row = Row(structure=structure)
        for forcefield in ("charmm36", "amber-openff"):
            row.results[forcefield] = run_setup(structure, forcefield, args.output_root)
        rows.append(row)
        old, new = row.results["charmm36"], row.results["amber-openff"]
        print(f"  {structure:6} charmm36={'ok ' if old.ok else 'FAIL':4} "
              f"amber-openff={'ok ' if new.ok else 'FAIL':4}  {row.verdict}")
        if row.verdict == "REGRESSION":
            print(f"         {new.detail[:150]}")

    regressions = [r for r in rows if r.verdict == "REGRESSION"]
    refusals = [r for r in rows if r.verdict == "refused"]
    both_ok = [r for r in rows if r.verdict == "both ok"]

    print("\n" + "=" * 64)
    print(f"  work under both force fields : {len(both_ok)}")
    print(f"  refused by the classifier    : {len(refusals)}"
          f"  ({', '.join(r.structure for r in refusals) or 'none'})")
    print(f"  REGRESSIONS under amber      : {len(regressions)}"
          f"  ({', '.join(r.structure for r in regressions) or 'none'})")
    print("=" * 64)

    if regressions:
        print("\nStructures that worked under CHARMM36 and do not under AMBER.")
        print("If these are ordinary proteins, the default should not change yet.")
    else:
        print("\nNo regressions: changing the default costs nothing on this set.")

    if args.json:
        args.json.write_text(json.dumps(
            [{"structure": r.structure, "verdict": r.verdict,
              **{ff: {"ok": o.ok, "seconds": round(o.seconds, 1), "detail": o.detail}
                 for ff, o in r.results.items()}}
             for r in rows], indent=2), encoding="utf-8")
        print(f"\nFull table written to {args.json}")

    return 1 if regressions else 0


if __name__ == "__main__":
    sys.exit(main())
