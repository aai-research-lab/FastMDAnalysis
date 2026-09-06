#!/usr/bin/env python
"""Compute V4's pre-registered statistic from the run's own output.

Reads `bfactor_comparison.dat` -- residue, simulated RMSF (nm), implied
RMSF (nm) -- and applies the residue set and statistic fixed in
`V4-experimental-anchor.md` before the result existed. Nothing here chooses
anything: the exclusions and the threshold are arguments with the
pre-registered values as defaults, so running it with no flags reproduces
the registered analysis and any deviation from it has to be typed out.

    python check_v4.py runs/v4-ubiquitin-anchors/analysis/bfactor_comparison/bfactor_comparison.dat
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

FIRST, LAST = 1, 72          # the folded beta-grasp domain
THRESHOLD = 0.60             # pre-registered, on Spearman rho


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dat", type=Path)
    parser.add_argument("--first", type=int, default=FIRST)
    parser.add_argument("--last", type=int, default=LAST)
    parser.add_argument("--threshold", type=float, default=THRESHOLD)
    args = parser.parse_args()

    table = np.loadtxt(args.dat)
    if table.ndim != 2 or table.shape[1] < 3:
        print(f"{args.dat} does not have the three columns this expects "
              f"(residue, simulated, implied); got shape {table.shape}",
              file=sys.stderr)
        return 2

    residue, simulated, implied = table[:, 0], table[:, 1], table[:, 2]
    keep = (residue >= args.first) & (residue <= args.last)
    dropped = int((~keep).sum())
    if keep.sum() < 3:
        print("fewer than three residues survive the range", file=sys.stderr)
        return 2

    r = float(np.corrcoef(simulated[keep], implied[keep])[0, 1])
    rho = float(spearmanr(simulated[keep], implied[keep]).statistic)
    ratio = float(np.mean(simulated[keep]) / np.mean(implied[keep]))

    registered = (args.first, args.last, args.threshold) == (FIRST, LAST, THRESHOLD)
    print(f"file            {args.dat}")
    print(f"residues        {int(keep.sum())} in {args.first}-{args.last} "
          f"({dropped} outside the range, excluded)")
    print(f"Spearman rho    {rho:+.3f}   "
          f"{'MEETS' if rho >= args.threshold else 'BELOW'} the "
          f"{args.threshold:.2f} threshold")
    print(f"Pearson r       {r:+.3f}   (reported, no threshold)")
    print(f"mean ratio      {ratio:.3f}   simulated / implied -- "
          f"{'above 1 as predicted' if ratio > 1 else 'BELOW 1, the informative failure'}")
    if not registered:
        print("\n!! these are not the pre-registered settings "
              f"({FIRST}-{LAST}, threshold {THRESHOLD:.2f}). Say so wherever "
              "the number is quoted.")

    # Where the disagreement sits, so a shortfall can be described rather
    # than only reported.
    if rho < args.threshold:
        gap = simulated[keep] - implied[keep]
        order = np.argsort(-np.abs(gap))[:8]
        print("\nlargest deviations (residue, simulated, implied, nm):")
        for i in order:
            print(f"  {int(residue[keep][i]):4d}  {simulated[keep][i]:.4f}  "
                  f"{implied[keep][i]:.4f}  {gap[i]:+.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
