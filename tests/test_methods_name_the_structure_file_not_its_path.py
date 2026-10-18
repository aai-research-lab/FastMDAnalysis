"""The report's Methods name the structure file, not the folder it was in.

"Starting coordinates were taken from `/tmp/.../1VII_villin.pdb`" published
a path on one computer; the file is identified by its name, its digest and
the entry its header names.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fastmdxplora.report.methods import methods_paragraphs


def _sentence(system: str, **structure) -> str:
    text = methods_paragraphs(Path("."), {"input": {
        "system": system, "form": "pdb_file",
        "structure": {"form": "pdb_file", **structure} if structure else None}}, {})
    start = text.find("Starting coordinates")
    return text[start:text.find("\n", start)]


@pytest.mark.parametrize("system", [
    "/tmp/pytest-of-someone/pytest-3/studies0/1VII_villin.pdb",
    "/Users/kunle/aai-research-lab/structures/1VII_villin.pdb",
    r"C:\Users\kunle\structures\1VII_villin.pdb",
    "structures/1VII_villin.pdb",
])
def test_a_path_is_said_as_the_file_s_name(system):
    said = _sentence(system)
    assert said.startswith("Starting coordinates were taken from `1VII_villin.pdb`")
    assert "/" not in said and "\\" not in said


def test_the_header_s_entry_is_still_said_beside_the_name():
    text = methods_paragraphs(Path("."), {"input": {
        "system": "/home/kunle/work/prepped.pdb", "form": "pdb_file",
        "structure": {"form": "pdb_file", "entry": "1VII", "sha256": "e" * 64}}}, {})
    assert "`prepped.pdb`, whose header identifies it as Protein Data Bank entry 1VII" in text
    assert "/home/kunle" not in text
