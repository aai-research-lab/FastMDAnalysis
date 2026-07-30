"""Structure metadata helpers for the live dashboard.

The dashboard needs lightweight, failure-isolated metadata for whichever PDB
is being displayed.  The helpers in this module deliberately avoid optional
scientific dependencies: they parse fixed-width PDB records, cache the result
by file modification time, and always return JSON-serialisable dictionaries.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Any

from fastmdxplora.live.ligand_detection import (
    AMINO_ACID_RESNAMES,
    ION_RESNAMES,
    WATER_RESNAMES,
    detect_ligands,
)

# Keep the lightweight metadata scanner bounded. The viewer now prefers the
# prepared solute PDB, so normal dashboard runs avoid parsing a full solvated
# topology while retaining the original safety limit for unusually large files.
_MAX_PDB_BYTES_FOR_INLINE_SCAN = 8 * 1024 * 1024
_COORD_SANITY_LIMIT = 1_000_000.0  # angstrom


def count_structure(path: str | Path) -> dict[str, Any]:
    """Return chain/residue/atom/ligand counters for a PDB file.

    Results are cached by absolute path, file size, and nanosecond mtime.  A
    fresh dictionary is returned to callers so API code may safely enrich the
    payload without mutating the cached value.
    """

    p = Path(path)
    if not p.is_file():
        return {"valid": False, "reason": "missing", "path": p.as_posix()}
    try:
        stat = p.stat()
    except OSError as exc:
        return {"valid": False, "reason": f"stat-error: {exc}", "path": p.as_posix()}
    if stat.st_size > _MAX_PDB_BYTES_FOR_INLINE_SCAN:
        return {
            "valid": False,
            "reason": "too-large",
            "path": p.as_posix(),
            "size": stat.st_size,
            "max_inline_bytes": _MAX_PDB_BYTES_FOR_INLINE_SCAN,
        }
    result = _count_structure_cached(
        str(p.resolve()),
        int(stat.st_mtime_ns),
        int(stat.st_size),
    )
    return dict(result)


@lru_cache(maxsize=32)
def _count_structure_cached(
    path_string: str,
    _mtime_ns: int,
    file_size: int,
) -> dict[str, Any]:
    p = Path(path_string)
    protein_chains: set[str] = set()
    all_chains: set[str] = set()
    protein_residue_keys: set[tuple[str, str, str]] = set()
    non_protein_residue_keys: set[tuple[str, str, str]] = set()
    water_residue_keys: set[tuple[str, str, str]] = set()
    ion_residue_keys: set[tuple[str, str, str]] = set()

    atoms = 0
    protein_atoms = 0
    hetatm_atoms = 0
    water_atoms = 0
    ion_atoms = 0
    min_coords = [float("inf"), float("inf"), float("inf")]
    max_coords = [float("-inf"), float("-inf"), float("-inf")]

    try:
        with p.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.startswith(("ATOM", "HETATM")):
                    continue

                atoms += 1
                is_hetatm = line.startswith("HETATM")
                if is_hetatm:
                    hetatm_atoms += 1

                chain_id = line[21:22].strip() or "A"
                resname = line[17:20].strip().upper()
                resi = line[22:27].strip()
                key = (chain_id, resname, resi)
                all_chains.add(chain_id)

                # Parse coordinates before classifying so the overall extent
                # remains meaningful even for solvated structures.
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if all(
                        -_COORD_SANITY_LIMIT <= value <= _COORD_SANITY_LIMIT
                        for value in (x, y, z)
                    ):
                        for idx, value in enumerate((x, y, z)):
                            min_coords[idx] = min(min_coords[idx], value)
                            max_coords[idx] = max(max_coords[idx], value)
                except (ValueError, IndexError):
                    pass

                if resname in WATER_RESNAMES:
                    water_atoms += 1
                    water_residue_keys.add(key)
                    continue
                if resname in ION_RESNAMES:
                    ion_atoms += 1
                    ion_residue_keys.add(key)
                    continue
                if resname in AMINO_ACID_RESNAMES:
                    protein_atoms += 1
                    protein_chains.add(chain_id)
                    protein_residue_keys.add(key)
                    continue

                # Any remaining residue is a ligand/cofactor candidate.  The
                # detector applies its own cofactor exclusions.
                non_protein_residue_keys.add(key)
    except OSError as exc:
        return {"valid": False, "reason": f"read-error: {exc}", "path": p.as_posix()}

    ligand_info = detect_ligands(non_protein_residue_keys)
    if min_coords[0] == float("inf"):
        bounding_extent = [0.0, 0.0, 0.0]
        centroid: list[float | None] = [None, None, None]
    else:
        bounding_extent = [max_coords[i] - min_coords[i] for i in range(3)]
        centroid = [(min_coords[i] + max_coords[i]) / 2 for i in range(3)]

    return {
        "valid": True,
        "path": p.as_posix(),
        "size": file_size,
        "atoms": atoms,
        "protein_atoms": protein_atoms,
        "hetatm_atoms": hetatm_atoms,
        "water_atoms": water_atoms,
        "ion_atoms": ion_atoms,
        "chains": sorted(protein_chains or all_chains),
        "all_chains": sorted(all_chains),
        "n_chains": len(protein_chains or all_chains),
        "protein_residues": len(protein_residue_keys),
        "water_residues": len(water_residue_keys),
        "ions": len(ion_residue_keys),
        "ligand_residues_total": len(ligand_info["instances"]),
        "ligand_resnames": ligand_info["resnames"],
        "ligand_instances": ligand_info["instances"],
        "extents_angstrom": [round(value, 2) for value in bounding_extent],
        "centroid_angstrom": [
            round(value, 2) if value is not None else None for value in centroid
        ],
    }


def ligand_atom_counts(path: str | Path) -> dict[str, int]:
    """Return ligand residue-name to atom-count mappings.

    Both ``ATOM`` and ``HETATM`` records are considered because some input
    generators encode small molecules using ``ATOM`` records.
    """

    p = Path(path)
    if not p.is_file():
        return {}
    counts: Counter[str] = Counter()
    try:
        with p.open("r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                resname = line[17:20].strip().upper()
                if (
                    resname in WATER_RESNAMES
                    or resname in ION_RESNAMES
                    or resname in AMINO_ACID_RESNAMES
                ):
                    continue
                counts[resname] += 1
    except OSError:
        return {}
    return dict(counts)
