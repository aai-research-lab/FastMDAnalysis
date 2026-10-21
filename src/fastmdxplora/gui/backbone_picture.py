"""A picture of a study's protein for its card on All studies, where the
study has plotted no figure yet.

The backbone as its Cα trace (a nucleic acid's by its phosphorus atoms),
read from the structure the Viewer would render first, seen down the
axis along which the molecule is thinnest so its widest face is shown,
coloured from the N terminus to the C terminus (viridis, its lightest
fifth left out to read on white) and shaded by depth, the nearer part
darker and wider. The trace is broken at a chain's end and wherever two
consecutive atoms are farther apart than a bond allows (a missing loop).
Written as SVG from the PDB's own lines, in milliseconds and with no
browser, and kept in memory by the file's size and time, never in the
study's folder.

Every file is read through `record_access`, so where the read is held to
the study (All studies' cards, `/api/study-thumbnail`) a structure linked
out of it, or one its records name outside it, is not read.
"""

from __future__ import annotations

import threading
from collections import OrderedDict
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np

from fastmdxplora import record_access

__all__ = ["backbone_points", "backbone_svg", "molecule_structure", "start_outside",
           "structure_for_picture"]

WIDTH, HEIGHT, MARGIN = 480, 200, 14
#: Consecutive Cα atoms are 3.8 Å apart; phosphorus atoms about 6 to 7.
_GAP = {"CA": 4.3, "P": 8.0}
#: Above this many atoms the trace takes every k-th, so the picture stays
#: small to send.
MOST_POINTS = 1500
#: Files larger than this are not read for a card.
MOST_BYTES = 200 * 1024 * 1024
#: Viridis.
_STOPS = ((0.0, (0x44, 0x01, 0x54)), (0.2, (0x41, 0x44, 0x87)), (0.4, (0x2A, 0x78, 0x8E)),
          (0.6, (0x22, 0xA8, 0x84)), (0.8, (0x7A, 0xD1, 0x51)), (1.0, (0xFD, 0xE7, 0x25)))
#: Each scheme's part of viridis, its background and the colour of N and C:
#: on paper its lightest fifth is left out, to read on white; on the dark
#: scheme its darkest, to read on the page, which shows through.
SCHEMES = {"light": ((0.0, 0.8), "#ffffff", "#444444"),
           "dark": ((0.2, 1.0), "none", "#c8c8c8")}
#: Pictures kept in memory, the oldest let go first.
MOST_KEPT = 256
_KEPT: OrderedDict[tuple[str, int, int, str], str | None] = OrderedDict()
_LOCK = threading.Lock()


def structure_for_picture(folder: str | Path) -> Path | None:
    """The structure a card's picture is made from: the study's own, as the
    Viewer would find it (`protein_preview.find_structure`), else that of
    the first of its runs with one."""
    base = Path(folder)
    found = _structure_in(base)
    if found is not None and found.suffix.lower() == ".pdb":
        return found
    runs = base / "runs"
    for run in sorted(record_access.iterdir(runs)) if record_access.is_dir(runs) else []:
        found = _structure_in(run) if record_access.is_dir(run) else None
        if found is not None and found.suffix.lower() == ".pdb":
            return found
    return None


def _structure_in(root: Path) -> Path | None:
    """`protein_preview.find_structure`, its files read through
    `record_access`: the same candidates in the same order, then the
    structure the study's Manifest names, where it can be read (one outside
    a held read is not looked at)."""
    import json

    from fastmdxplora.gui.protein_preview import STRUCTURE_CANDIDATES

    for rel in STRUCTURE_CANDIDATES:
        if record_access.is_file(root / rel):
            return root / rel
    try:
        manifest = json.loads(record_access.read_text(root / "manifest.json", encoding="utf-8"))
    except (OSError, ValueError):
        return None
    system = manifest.get("system") if isinstance(manifest, dict) else None
    if system:
        for path in (Path(str(system)), root / str(system)):
            if (path.suffix.lower() == ".pdb" and record_access.can_read(path)
                    and record_access.is_file(path)):
                return path
    return None


def _pieces(path: Path) -> tuple[list[np.ndarray], str]:
    """The backbone's atoms as pieces of consecutive residues, in Å."""
    atoms: dict[str, list[tuple[str, str, tuple[float, float, float]]]] = {"CA": [], "P": []}
    seen: set[tuple[str, str, str]] = set()
    with record_access.open_text(path, encoding="utf-8", errors="replace") as lines:
        for line in lines:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                continue
            name = line[12:16]
            # " CA " is an alpha carbon; a calcium ion is "CA  ".
            kind = "CA" if name == " CA " else "P" if name.strip() == "P" else None
            if kind is None or line[16] not in " A":
                continue
            residue = (kind, line[21], line[22:27])
            if residue in seen:
                continue
            seen.add(residue)
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            atoms[kind].append((line[21], line[22:27], xyz))
    kind = "CA" if len(atoms["CA"]) >= 3 else "P"
    chosen = atoms[kind]
    pieces: list[list[tuple[float, float, float]]] = []
    for k, (chain, _, xyz) in enumerate(chosen):
        joined = k > 0 and chosen[k - 1][0] == chain and float(np.linalg.norm(
            np.subtract(xyz, chosen[k - 1][2]))) <= _GAP[kind]
        if joined:
            pieces[-1].append(xyz)
        else:
            pieces.append([xyz])
    return [np.array(piece) for piece in pieces], kind


def _colour(t: float, scheme: str = "light") -> str:
    low_end, high_end = SCHEMES[scheme][0]
    t = low_end + min(max(t, 0.0), 1.0) * (high_end - low_end)
    (a, low), (b, high) = next(pair for pair in zip(_STOPS, _STOPS[1:]) if t <= pair[1][0])
    f = (t - a) / (b - a)
    return "#%02x%02x%02x" % tuple(round(lo + f * (hi - lo)) for lo, hi in zip(low, high))


def _rendered(path: Path, scheme: str = "light") -> str | None:
    pieces, kind = _pieces(path)
    count = sum(len(piece) for piece in pieces)
    if count < 3:
        return None
    step = max(1, -(-count // MOST_POINTS))
    pieces = [piece[::step] if len(piece) > step else piece for piece in pieces]
    every = np.concatenate(pieces)
    centre = every.mean(axis=0)
    _, _, axes = np.linalg.svd(every - centre, full_matrices=False)
    across, up = axes[0], axes[1]
    # The N terminus on the left, and the depth axis right-handed, so the
    # nearer part is the one facing the reader.
    if (every[0] - centre) @ across > (every[-1] - centre) @ across:
        across = -across
    toward = np.cross(across, up)
    flat = [np.stack([(piece - centre) @ across, -((piece - centre) @ up),
                      (piece - centre) @ toward], axis=1) for piece in pieces]
    allflat = np.concatenate(flat)
    low, high = allflat.min(axis=0), allflat.max(axis=0)
    span = np.maximum(high - low, 1e-6)
    scale = min((WIDTH - 2 * MARGIN) / span[0], (HEIGHT - 2 * MARGIN) / span[1])
    offset = np.array([(WIDTH - scale * span[0]) / 2, (HEIGHT - scale * span[1]) / 2])
    width = 3.0 if count < 150 else 2.2 if count < 600 else 1.4
    segments = []
    done = 0
    total = max(len(allflat) - 1, 1)
    for piece in flat:
        xy = (piece[:, :2] - low[:2]) * scale + offset
        near = (piece[:, 2] - low[2]) / span[2]
        for k in range(len(piece) - 1):
            depth = float((near[k] + near[k + 1]) / 2)
            segments.append((depth, (*xy[k], *xy[k + 1]), (done + k + 0.5) / total))
        done += len(piece)
    segments.sort(key=lambda kept: kept[0])
    title = (f"The backbone of {count:,} residues ({'Cα' if kind == 'CA' else 'phosphorus'} "
             "trace), coloured from the N terminus (purple) to the C terminus "
             f"({'green' if scheme == 'light' else 'yellow'}), the nearer part "
             f"{'darker' if scheme == 'light' else 'brighter'}")
    _, background, ink = SCHEMES[scheme]
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}" '
           f'width="{WIDTH}" height="{HEIGHT}" role="img"><title>{escape(title)}</title>',
           f'<rect width="{WIDTH}" height="{HEIGHT}" fill="{background}"/>',
           '<g fill="none" stroke-linecap="round">']
    for depth, (x1, y1, x2, y2), t in segments:
        out.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                   f'stroke="{_colour(t, scheme)}" stroke-opacity="{0.35 + 0.65 * depth:.2f}" '
                   f'stroke-width="{width * (0.7 + 0.5 * depth):.2f}"/>')
    out.append("</g>")
    first = (flat[0][0, :2] - low[:2]) * scale + offset
    last = (flat[-1][-1, :2] - low[:2]) * scale + offset
    for (x, y), label in ((first, "N"), (last, "C")):
        out.append(f'<text x="{min(max(x, 8), WIDTH - 8):.1f}" '
                   f'y="{min(max(y - 6, 11), HEIGHT - 3):.1f}" font-family="sans-serif" '
                   f'font-size="11" text-anchor="middle" fill="{ink}">{label}</text>')
    out.append("</svg>")
    return "".join(out)


def backbone_svg(folder: str | Path, scheme: str = "light") -> str | None:
    """The card's picture of the study's backbone as SVG, in the scheme
    given, or None where the study has no structure with a backbone to show."""
    scheme = scheme if scheme in SCHEMES else "light"
    path = structure_for_picture(folder)
    if path is None:
        return None
    try:
        stat = record_access.stat(path)
    except OSError:
        return None
    if stat.st_size > MOST_BYTES:
        return None
    key = (str(path.resolve()), stat.st_size, stat.st_mtime_ns, scheme)
    with _LOCK:
        if key in _KEPT:
            _KEPT.move_to_end(key)
            return _KEPT[key]
    try:
        picture = _rendered(path, scheme)
    except (OSError, ValueError, np.linalg.LinAlgError):
        picture = None
    with _LOCK:
        _KEPT[key] = picture
        while len(_KEPT) > MOST_KEPT:
            _KEPT.popitem(last=False)
    return picture


# ---------------------------------------------------------------------------
# The molecule on All studies' rows and cards
# ---------------------------------------------------------------------------

#: The points a row's or a card's picture is sent: enough for a protein's
#: shape at 260 px, few enough that two hundred studies list quickly.
LIST_POINTS = 360
#: A structure larger than this is not read for the list (a solvated
#: system written whole): the saved frames leave the water out.
LIST_BYTES = 40 * 1024 * 1024
_POINTS_KEPT: OrderedDict[tuple[str, int, int], dict | None] = OrderedDict()


def molecule_structure(folder: str | Path) -> tuple[Path, str] | None:
    """The structure a study's molecule is drawn from on All studies, and
    what it is: the latest frame its run wrote (`simulation/live_frame.pdb`,
    the run's last where it has ended), else the structure the Viewer draws
    first, else the file it is to start from where nothing ran yet. For a
    study of several runs, that of its first run with one."""
    base = Path(folder)
    runs = base / "runs"
    for where in [base] + (sorted(p for p in record_access.iterdir(runs)
                                  if record_access.is_dir(p))
                           if record_access.is_dir(runs) else []):
        latest = where / "simulation" / "live_frame.pdb"
        if record_access.is_file(latest):
            return latest, "latest frame"
    found = structure_for_picture(base)
    if found is not None:
        return found, "structure"
    path = _start_named(base)
    if path is not None and record_access.can_read(path) and record_access.is_file(path):
        return path, "start"
    return None


def _start_named(base: Path) -> Path | None:
    """The .pdb the study's Config names to start from, where it names one."""
    from fastmdxplora.gui.workspace import _read_yaml

    config = _read_yaml(base / "resolved_config.yml") or _read_yaml(base / "exploration.yml")
    systems = config.get("systems") if isinstance(config, dict) else None
    first = (systems[0] if isinstance(systems, list) and systems
             and isinstance(systems[0], dict) else {})
    named = str(first.get("system") or "")
    if not named.lower().endswith(".pdb"):
        return None
    path = Path(named).expanduser()
    return path if path.is_absolute() else base / path


def start_outside(folder: str | Path) -> str:
    """The name of the .pdb a study is to start from where it is outside
    the study, which a read held to the study does not open (nor where it
    leads); '' otherwise. Told without opening it."""
    import os

    path = _start_named(Path(folder))
    if path is None or record_access.can_read(path) or not os.path.lexists(path):
        return ""
    return path.name


def backbone_points(folder: str | Path) -> dict | None:
    """The study's backbone as the page draws it at any size: each piece's
    points in the plane of the molecule's two widest axes (``x`` from 0 to
    ``aspect``, ``y`` from 0 to 1, N terminus on the left), each with its
    depth from 0 (far) to 1 (near) and its place along the chain from 0 (N)
    to 1 (C). ``from`` says which structure (`molecule_structure`). Kept in
    memory by the file's size and time; None where there is none to draw."""
    found = molecule_structure(folder)
    if found is None:
        return None
    path, source = found
    try:
        stat = record_access.stat(path)
    except OSError:
        return None
    if stat.st_size > LIST_BYTES:
        return None
    key = (str(path.resolve()), stat.st_size, stat.st_mtime_ns)
    with _LOCK:
        if key in _POINTS_KEPT:
            _POINTS_KEPT.move_to_end(key)
            kept = _POINTS_KEPT[key]
            return None if kept is None else {**kept, "from": source}
    try:
        drawn = _points(path)
    except (OSError, ValueError, np.linalg.LinAlgError):
        drawn = None
    with _LOCK:
        _POINTS_KEPT[key] = drawn
        while len(_POINTS_KEPT) > MOST_KEPT:
            _POINTS_KEPT.popitem(last=False)
    return None if drawn is None else {**drawn, "from": source}


def _points(path: Path) -> dict | None:
    pieces, kind = _pieces(path)
    count = sum(len(piece) for piece in pieces)
    if count < 3:
        return None
    step = max(1, -(-count // LIST_POINTS))
    pieces = [piece[::step] if len(piece) > step else piece for piece in pieces]
    every = np.concatenate(pieces)
    centre = every.mean(axis=0)
    _, _, axes = np.linalg.svd(every - centre, full_matrices=False)
    across, up = axes[0], axes[1]
    if (every[0] - centre) @ across > (every[-1] - centre) @ across:
        across = -across
    toward = np.cross(across, up)
    flat = [np.stack([(piece - centre) @ across, -((piece - centre) @ up),
                      (piece - centre) @ toward], axis=1) for piece in pieces]
    allflat = np.concatenate(flat)
    low, high = allflat.min(axis=0), allflat.max(axis=0)
    span = np.maximum(high - low, 1e-6)
    total = max(len(allflat) - 1, 1)
    out, done = [], 0
    for piece in flat:
        scaled = (piece - low) / span[1]
        depth = (piece[:, 2] - low[2]) / span[2]
        out.append([[round(float(x), 3), round(float(y), 3), round(float(d), 2),
                     round((done + k) / total, 3)]
                    for k, ((x, y, _), d) in enumerate(zip(scaled, depth))])
        done += len(piece)
    return {"pieces": out, "aspect": round(float(span[0] / span[1]), 3), "residues": count,
            "atom": "CA" if kind == "CA" else "P"}
