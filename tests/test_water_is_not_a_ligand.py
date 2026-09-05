"""Retained crystallographic water has no chemistry to fetch either.

The sibling of `test_an_ion_is_not_a_ligand`, and found the same way --
by running. `keep_water: true` marks water SIMULATE, which routes it into
the ligand path, and the exemption written directly above for ions was
never extended to it.

Both outcomes were wrong and the quieter one was worse. From a local file
the run refused, advising the reader to `curl` chemistry for water:

    heterogens: auto identified components to simulate (HOH x887), but
    their chemistry can only be retrieved for a structure given by PDB
    identifier ... curl -O .../XXX_ideal.sdf

From a PDB identifier there is no refusal at all. `HOH` is a real CCD
component, so the fetch succeeds, and the retained waters are
parameterised by openff as a small molecule -- a second water model, in
the same box as the solvent poured around them, with nothing said. That
is the failure mode this repository's own notes call the worst kind: a
wrong answer that completes.

The water model is part of the force field. `amber14/tip3p.xml` carries
HOH directly, so `HOH_ideal.sdf` describes a molecule that is already
parameterised.
"""

from __future__ import annotations

import io
import logging
from pathlib import Path

import pytest

from fastmdxplora.setup.heterogens import WATER_NAMES, resolve
from fastmdxplora.setup.pipeline import _auto_ligands


def _structure(tmp_path: Path, *, waters: int = 3, ligand: bool = False,
               water_name: str = "HOH") -> Path:
    """A short peptide, some crystallographic waters, optionally a ligand."""
    lines, n = [], 0
    for i, x in enumerate((10.0, 13.8), 1):
        n += 1
        lines.append(f"ATOM  {n:5d}  CA  ALA A{i:4d}    "
                     f"{x:8.3f}{10.0:8.3f}{10.0:8.3f}  1.00  0.00           C")
    for w in range(waters):
        n += 1
        lines.append(f"HETATM{n:5d}  O   {water_name:>3s} A{600 + w:4d}    "
                     f"{20.0 + 3.5 * w:8.3f}{10.0:8.3f}{10.0:8.3f}"
                     f"  1.00  0.00           O")
    if ligand:
        for i, dx in enumerate((0.0, 1.4), 1):
            n += 1
            lines.append(f"HETATM{n:5d}  C{i}  BNZ A 500    "
                         f"{40.0 + dx:8.3f}{10.0:8.3f}{10.0:8.3f}"
                         f"  1.00  0.00           C")
    path = tmp_path / "system.pdb"
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")
    return path


KEEP = {"heterogens": "auto", "keep_water": True}


class TestWaterIsKeptWithoutBeingFetched:

    def test_keeping_water_from_a_local_file_is_not_refused(
            self, tmp_path: Path) -> None:
        """The whole failure. A pure water box could not be prepared at all,
        and neither could any structure whose waters were asked for."""
        assert _auto_ligands(KEEP, _structure(tmp_path), tmp_path, None) == []

    def test_it_says_what_it_kept_and_why(self, tmp_path: Path) -> None:
        """Captured from the module's own logger rather than through
        `caplog`. `caplog` attaches to the root and depends on propagation
        surviving whatever the rest of the suite did to logging: this test
        passed alone and on a clean 3.9, and came back with an empty
        `caplog.text` inside the full run on CI. A handler on the named
        logger does not care what else has happened."""
        logger = logging.getLogger("fastmdx.setup")
        stream = io.StringIO()
        handler = logging.StreamHandler(stream)
        handler.setLevel(logging.INFO)
        previous = logger.level
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        try:
            _auto_ligands(KEEP, _structure(tmp_path), tmp_path, None)
        finally:
            logger.removeHandler(handler)
            logger.setLevel(previous)
        said = stream.getvalue()
        assert "HOH" in said
        assert "water model is part of the force field" in said

    def test_no_chemistry_is_fetched_even_with_an_entry_to_look_it_up(
            self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """The silent half. With an identifier the fetch succeeds, and the
        retained waters come back parameterised as a small molecule -- a
        different water model from the solvent, in the same box."""
        import fastmdxplora.setup.ccd as ccd

        def refuse(*args, **kwargs):  # pragma: no cover - must not be reached
            raise AssertionError(
                "chemistry was fetched for water; the force field already "
                "carries it")

        monkeypatch.setattr(ccd, "fetch_chemistry", refuse)
        assert _auto_ligands(KEEP, _structure(tmp_path), tmp_path, "1ABC") == []

    @pytest.mark.parametrize("name", sorted(WATER_NAMES))
    def test_every_name_the_classifier_calls_water(
            self, tmp_path: Path, name: str) -> None:
        """`WATER_NAMES` holds eight spellings. Exempting only HOH would
        leave a GROMACS-written SOL or an AMBER WAT taking the ligand path,
        which is the shape of AUD9 -- one naming convention covered and the
        rest silently not."""
        structure = _structure(tmp_path, water_name=name)
        kept = [d for d in resolve(structure, keep_water=True)
                if d.resname == name]
        if not kept:
            pytest.skip(f"{name} is not read as a residue by the parser here")
        assert _auto_ligands(KEEP, structure, tmp_path, None) == []


class TestTheRestOfThePathIsUnchanged:

    def test_an_organic_ligand_still_refuses_without_an_entry(
            self, tmp_path: Path) -> None:
        with pytest.raises(ValueError, match="BNZ"):
            _auto_ligands(KEEP, _structure(tmp_path, ligand=True),
                          tmp_path, None)

    def test_water_is_not_named_in_that_refusal(self, tmp_path: Path) -> None:
        """Naming it sends somebody after `HOH_ideal.sdf`, which is the
        fault this file is about."""
        with pytest.raises(ValueError) as caught:
            _auto_ligands(KEEP, _structure(tmp_path, ligand=True),
                          tmp_path, None)
        first = str(caught.value).split("\n")[0]
        assert "BNZ" in first and "HOH" not in first

    def test_discarding_water_is_still_the_default(
            self, tmp_path: Path) -> None:
        """Without `keep_water` the classifier drops it as solvent, so it
        never reaches this function at all."""
        assert _auto_ligands({"heterogens": "auto"}, _structure(tmp_path),
                             tmp_path, None) == []
