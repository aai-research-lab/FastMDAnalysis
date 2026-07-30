"""Protein preview generation for the local live dashboard."""

from __future__ import annotations

import json
import shutil
import subprocess
from colorsys import hsv_to_rgb
from pathlib import Path
from typing import Any


PREVIEW_CANDIDATES = (
    "report/dashboard_assets/protein_preview.png",
    "simulation/protein_preview.png",
)

# Prefer the prepared solute for browser rendering.  The simulation topology
# normally contains the entire solvent box and can easily be tens of megabytes,
# which makes a browser viewer slow and also obscures the protein/ligand.  Live
# frames and trajectory playback still use the full simulation topology where
# atom-for-atom correspondence is required.
STRUCTURE_CANDIDATES = (
    "setup/prepared.pdb",
    "simulation/final.pdb",
    "simulation/topology.pdb",
    "setup/topology.pdb",
    "setup/solvated.pdb",
)

AMINO_ACID_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}


def protein_preview_payload(root: str | Path, *, regenerate: bool = False) -> dict[str, Any]:
    """Return preview status, generating an image from available structure if needed."""
    project_root = Path(root)
    structure = find_structure(project_root)
    existing = _find_existing_preview(project_root)
    if structure is None and existing is not None and not regenerate:
        return _available(project_root, existing, static_mode="existing", message="Existing protein preview found.")
    if structure is None:
        return {
            "available": False,
            "static_available": False,
            "viewer_available": False,
            "message": "No topology/PDB found yet.",
        }

    if existing is not None and not regenerate:
        return _available(project_root, existing, static_mode="pymol", message="PyMOL render found.")

    output_path = _preview_output_path(project_root)
    pymol = _find_pymol_executable(project_root)
    if pymol:
        try:
            _render_with_pymol(pymol, structure, output_path)
            return _available(project_root, output_path, static_mode="pymol", message="PyMOL render generated.")
        except Exception as exc:  # noqa: BLE001
            return _viewer_only(project_root, f"PyMOL preview unavailable. Showing schematic fallback. {exc}")
    return {
        **_viewer_only(
            project_root,
            "PyMOL preview unavailable. Showing schematic fallback. Install pymol-open-source "
            "or run from the fastmdx-local micromamba environment.",
        ),
    }


def _find_pymol_executable(root: Path) -> str | None:
    direct = shutil.which("pymol")
    if direct:
        return direct
    candidates: list[Path] = []
    for base in [Path.cwd(), *root.resolve().parents]:
        candidates.append(base / ".micromamba-root" / "envs" / "fastmdx-local" / "bin" / "pymol")
    for path in candidates:
        if path.is_file():
            return path.as_posix()
    return None


def _find_existing_preview(root: Path) -> Path | None:
    for rel in PREVIEW_CANDIDATES:
        path = root / rel
        if path.is_file():
            return path
    return None


def find_structure(root: str | Path) -> Path | None:
    """Return the best available PDB/topology file for live preview rendering."""
    return _find_structure(Path(root))


def _find_structure(root: Path) -> Path | None:
    for rel in STRUCTURE_CANDIDATES:
        path = root / rel
        if path.is_file():
            return path

    manifest = _load_json(root / "manifest.json")
    system = manifest.get("system")
    if system:
        candidates = [Path(str(system)), root / str(system)]
        for path in candidates:
            if path.is_file() and path.suffix.lower() == ".pdb":
                return path
    return None


def _load_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _preview_output_path(root: Path) -> Path:
    report_assets = root / "report" / "dashboard_assets"
    if report_assets.exists():
        return report_assets / "protein_preview.png"
    return root / "simulation" / "protein_preview.png"


def _available(root: Path, path: Path, *, static_mode: str, message: str) -> dict[str, Any]:
    rel = path.relative_to(root).as_posix()
    image_url = f"/artifacts/{rel}?v={int(path.stat().st_mtime)}"
    structure_info = _structure_info(root)
    return {
        "available": True,
        "static_available": True,
        "static_mode": static_mode,
        "static_image_url": image_url,
        "path": rel,
        "href": image_url,
        "image_url": image_url,
        **structure_info,
        "message": message,
    }


def _viewer_only(root: Path, message: str) -> dict[str, Any]:
    return {
        "available": True,
        "static_available": False,
        "static_mode": None,
        "static_image_url": None,
        "path": None,
        "href": None,
        "image_url": None,
        **_structure_info(root),
        "message": message,
    }


def _structure_info(root: Path) -> dict[str, Any]:
    structure = find_structure(root)
    structure_rel = None
    structure_url = None
    if structure is not None:
        try:
            structure_rel = structure.relative_to(root).as_posix()
            structure_url = f"/structure/topology.pdb?v={int(structure.stat().st_mtime)}"
        except ValueError:
            structure_rel = structure.as_posix()
    return {
        "structure_path": structure_rel,
        "structure_url": structure_url,
        "structure_available": structure_url is not None,
        "viewer_available": structure_url is not None and viewer_asset_available(),
        "viewer_mode": "3dmol" if structure_url is not None and viewer_asset_available() else None,
        "fallback_available": structure_url is not None,
        "fallback_mode": "schematic" if structure_url is not None else None,
    }


def viewer_asset_available() -> bool:
    """Return True when the bundled 3Dmol browser asset is available."""
    return (Path(__file__).with_name("static") / "3Dmol-min.js").is_file()


def _render_with_pymol(pymol: str, structure: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    script = output_path.with_suffix(".pml")
    residue_colors = _residue_color_commands(structure)
    script.write_text(
        "\n".join(
            [
                "reinitialize",
                "viewport 1800, 1400",
                f"load {structure.as_posix()}, prot",
                "hide everything",
                "show cartoon, prot",
                "spectrum count, rainbow, prot, byres=1",
                *residue_colors,
                "bg_color white",
                "set ray_opaque_background, off",
                "set antialias, 2",
                "set ray_trace_mode, 1",
                "set cartoon_fancy_helices, 1",
                "set cartoon_smooth_loops, 1",
                "set cartoon_tube_radius, 0.45",
                "set cartoon_sampling, 14",
                "orient prot",
                "center prot",
                "zoom prot, 1.8",
                "ray 1800, 1200",
                f"png {output_path.as_posix()}, width=1800, height=1200, dpi=300, ray=1",
                "quit",
            ]
        ),
        encoding="utf-8",
    )
    try:
        subprocess.run(
            [pymol, "-cq", str(script)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=45,
        )
    finally:
        try:
            script.unlink()
        except OSError:
            pass
    if not output_path.is_file():
        raise RuntimeError("PyMOL did not write a preview image")


def _residue_color_commands(structure: Path) -> list[str]:
    residues: list[tuple[str, str]] = []
    seen: set[tuple[str, str]] = set()
    try:
        lines = structure.read_text(encoding="utf-8", errors="ignore").splitlines()
    except OSError:
        return []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip().upper()
        if resname not in AMINO_ACID_RESNAMES:
            continue
        chain = line[21].strip()
        resi = line[22:27].strip()
        if not resi:
            continue
        key = (chain, resi)
        if key not in seen:
            seen.add(key)
            residues.append(key)
    if len(residues) < 2:
        return []
    commands: list[str] = []
    denom = max(1, len(residues) - 1)
    for index, (chain, resi) in enumerate(residues):
        hue = 0.72 * (index / denom)
        red, green, blue = hsv_to_rgb(hue, 0.78, 1.0)
        color_name = f"fastmdx_res_{index + 1}"
        selection = f"prot and resi {resi}"
        if chain:
            selection += f" and chain {chain}"
        commands.append(f"set_color {color_name}, [{red:.3f}, {green:.3f}, {blue:.3f}]")
        commands.append(f"color {color_name}, {selection}")
    return commands
