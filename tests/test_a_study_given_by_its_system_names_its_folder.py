"""A study given by its system alone, with no output folder, is named by it.

`FastMDXplora(system=...)` named the folder from `self.config`, an attribute
it never sets, so the call `docs/api.md` gives stopped with an
AttributeError whenever `output_dir` was left to its default.
"""

from __future__ import annotations

import re

from fastmdxplora import FastMDXplora


def test_no_output_dir_names_the_folder_by_the_system(tmp_path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    structure = tmp_path / "trp cage.pdb"
    structure.write_text(
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "END\n", encoding="utf-8")
    study = FastMDXplora(system=str(structure))
    assert re.fullmatch(r"fastmdxplora_trp-cage_study_\d{14}", str(study.output_dir))
    assert (tmp_path / study.output_dir).is_dir()
