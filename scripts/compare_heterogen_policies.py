#!/usr/bin/env python
"""Compare the `drop` and `auto` heterogen policies, to decide the default.

`drop` removes every non-standard residue and prepares the protein alone.
`auto` decides per component: solvent and crystallization additives are
discarded, coordinated metals kept, a bound ligand retrieved and
parameterized, and anything the structure does not determine is refused.

`auto` is better science when it works. The question this answers is what it
costs: whether structures that prepare cleanly under `drop` stop under `auto`,
and whether the ones that stop are stopping for good reasons.

    python scripts/compare_heterogen_policies.py --output-root /tmp/hetcmp

Both policies run on the same force field, so nothing here is about the force
field. Setup only, seconds per structure.

Reading the result:

* `auto` refusing a haem protein is correct: a haem needs haem parameters.
* `auto` refusing an ordinary protein-ligand complex is not, and means the
  default should stay `drop`.
* `auto` succeeding where `drop` succeeded, plus preparing the ligand, is the
  case for changing the default.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

#: Chosen to span what the classifier must handle: plain proteins, bound
#: ligands, coordinated metals, cofactors, glycans, nucleic acid, and
#: crystallization additives with nothing else present.
DEFAULT_STRUCTURES = [
    "1L2Y",   # Trp-cage: nothing but protein
    "1UBQ",   # ubiquitin
    "1CRN",   # crambin, disulfides
    "1AKI",   # lysozyme
    "4W51",   # T4 lysozyme L99A, apo
    "4W52",   # T4 lysozyme L99A + benzene: the case auto exists for
    "1STP",   # streptavidin + biotin, a tight pocket
    "1HVR",   # HIV-1 protease + inhibitor, modified cysteine
    "6LU7",   # SARS-CoV-2 Mpro + covalent inhibitor
    "1A2P",   # barnase, zinc
    "1ZNI",   # insulin, six zincs
    "4HHB",   # haemoglobin, haem: should refuse
    "1IGT",   # antibody, N-glycans: should refuse
    "1BNA",   # B-DNA
]

REFUSAL = re.compile(
    r"does not determine what should be simulated|needs dedicated parameters|"
    r"a sugar can be|covalently bonded|cannot parameterize small molecules|"
    r"can only be retrieved"
)


@dataclass
class Outcome:
    ok: bool
    seconds: float
    detail: str = ""
    ligands: tuple[str, ...] = ()

    @property
    def refused(self) -> bool:
        return bool(REFUSAL.search(self.detail))


@dataclass
class Row:
    structure: str
    results: dict[str, Outcome] = field(default_factory=dict)

    @property
    def verdict(self) -> str:
        drop, auto = self.results["drop"], self.results["auto"]
        if drop.ok and auto.ok:
            return "richer" if auto.ligands else "same"
        if drop.ok and auto.refused:
            return "auto refuses"
        if drop.ok and not auto.ok:
            return "AUTO BREAKS"
        if not drop.ok and auto.ok:
            return "auto better"
        return "both fail"


def run_setup(structure: str, policy: str, root: Path, forcefield: str) -> Outcome:
    out = root / policy / structure
    started = time.time()
    process = subprocess.run(
        ["fastmdx", "setup", "--system", structure, "--output", str(out),
         "--heterogens", policy, "--forcefield", forcefield],
        capture_output=True, text=True, timeout=900,
    )
    elapsed = time.time() - started
    combined = (process.stdout or "") + (process.stderr or "")
    produced = (out / "setup" / "system.xml").is_file()

    ligands: tuple[str, ...] = ()
    ligand_dir = out / "setup" / "ligands"
    if ligand_dir.is_dir():
        ligands = tuple(sorted(p.stem for p in ligand_dir.glob("*.sdf")))

    detail = ""
    if not produced:
        lines = [l for l in combined.splitlines() if "ERROR" in l or " - " in l]
        detail = " ".join(lines[-3:])[:500] if lines else combined[-300:]
    return Outcome(ok=produced, seconds=elapsed, detail=detail, ligands=ligands)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("structures", nargs="*", default=None)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--forcefield", default="amber-openff",
                        help="the same one for both policies")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    structures = args.structures or DEFAULT_STRUCTURES
    args.output_root.mkdir(parents=True, exist_ok=True)

    rows: list[Row] = []
    print(f"Comparing heterogens=drop against heterogens=auto on "
          f"{args.forcefield}, {len(structures)} structures\n")
    for structure in structures:
        row = Row(structure=structure)
        for policy in ("drop", "auto"):
            row.results[policy] = run_setup(structure, policy,
                                            args.output_root, args.forcefield)
        rows.append(row)
        drop, auto = row.results["drop"], row.results["auto"]
        ligands = f"  ligands: {', '.join(auto.ligands)}" if auto.ligands else ""
        print(f"  {structure:6} drop={'ok ' if drop.ok else 'FAIL':4} "
              f"auto={'ok ' if auto.ok else 'FAIL':4}  {row.verdict}{ligands}")
        if row.verdict in {"AUTO BREAKS", "auto refuses"}:
            print(f"         {auto.detail[:150]}")

    breaks = [r for r in rows if r.verdict == "AUTO BREAKS"]
    refuses = [r for r in rows if r.verdict == "auto refuses"]
    richer = [r for r in rows if r.verdict == "richer"]
    same = [r for r in rows if r.verdict == "same"]

    print("\n" + "=" * 68)
    print(f"  auto prepares a ligand drop would have discarded : {len(richer)}"
          f"  ({', '.join(r.structure for r in richer) or 'none'})")
    print(f"  identical outcome                                : {len(same)}")
    print(f"  auto refuses, with a reason                      : {len(refuses)}"
          f"  ({', '.join(r.structure for r in refuses) or 'none'})")
    print(f"  auto BREAKS what drop prepared                   : {len(breaks)}"
          f"  ({', '.join(r.structure for r in breaks) or 'none'})")
    print("=" * 68)

    if breaks:
        print("\nauto fails where drop succeeded, and not by refusing: these are")
        print("bugs, not judgements. The default should stay drop.")
    elif refuses:
        print("\nauto refuses some structures that drop prepared. That is the")
        print("intended behaviour where the refusal is correct: check each")
        print("reason above. A refusal on an ordinary complex argues against")
        print("changing the default; one on a haem or a glycan does not.")
    else:
        print("\nauto costs nothing on this set and prepares ligands drop would")
        print("have silently discarded. The default can change.")

    if args.json:
        args.json.write_text(json.dumps(
            [{"structure": r.structure, "verdict": r.verdict,
              **{p: {"ok": o.ok, "seconds": round(o.seconds, 1),
                     "ligands": list(o.ligands), "detail": o.detail}
                 for p, o in r.results.items()}}
             for r in rows], indent=2), encoding="utf-8")
        print(f"\nFull table written to {args.json}")

    return 1 if breaks else 0


if __name__ == "__main__":
    sys.exit(main())
