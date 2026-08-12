"""A run recorded the string somebody typed, and called it the input.

`system: 4hhb_cleaned.pdb` is an answer for about a week. The file gets moved,
the directory tidied, a second copy made under a name that sorts better --
and the run no longer says which structure it used. The report reads the same
field, so a methods section ends up stating that coordinates came from a
filename on somebody's laptop, which is not something a reader can check.

Six months on, the question asked of a finished study is *which entry was
this, and when*. Nothing in the manifest answered either.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from fastmdxplora.provenance import (
    described_structure,
    structure_provenance,
)
from fastmdxplora.report.methods import methods_paragraphs
from fastmdxplora.setup.pipeline import run as setup_run

MINI = (
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
    "END\n"
)
HEADER = "HEADER    OXYGEN TRANSPORT                        07-MAR-84   4HHB\n"


@pytest.fixture
def deposited(tmp_path: Path) -> Path:
    """A structure that has been through a preparation step and kept its
    header -- which is the ordinary case, and the one the path fails."""
    p = tmp_path / "4hhb_cleaned.pdb"
    p.write_text(HEADER + MINI)
    return p


@pytest.fixture
def built_by_hand(tmp_path: Path) -> Path:
    p = tmp_path / "built.pdb"
    p.write_text(MINI)
    return p


class TestTheDigestIdentifiesTheStructure:
    """The path is the field that stops being true. The bytes do not."""

    def test_the_same_file_gives_the_same_digest(self, deposited: Path) -> None:
        first = structure_provenance(str(deposited), "pdb_file", deposited)
        second = structure_provenance(str(deposited), "pdb_file", deposited)
        assert first["sha256"] == second["sha256"]

    def test_a_file_edited_since_no_longer_matches(self, deposited: Path) -> None:
        """What makes the digest worth recording: a run and a file either
        agree or they do not, and this is how somebody finds out."""
        before = structure_provenance(str(deposited), "pdb_file", deposited)
        deposited.write_text(HEADER + MINI.replace("0.000   0.000", "0.001   0.000"))
        after = structure_provenance(str(deposited), "pdb_file", deposited)
        assert before["sha256"] != after["sha256"]

    def test_the_record_is_a_snapshot_not_a_pointer(self, deposited: Path) -> None:
        """Taken when the file arrived, so what the run then did to it -- the
        chain selection that rewrites input.pdb -- does not rewrite history."""
        recorded = structure_provenance(str(deposited), "pdb_file", deposited)
        digest = recorded["sha256"]
        deposited.write_text("ATOM      1  N   GLY B   1       0.0   0.0   0.0\nEND\n")
        assert recorded["sha256"] == digest


class TestTheEntryTheFileNamesItself:
    """The answer to "which entry" for exactly the case where the path has
    stopped giving one."""

    def test_a_prepared_file_still_names_its_entry(self, deposited: Path) -> None:
        got = structure_provenance(str(deposited), "pdb_file", deposited)
        assert got["entry"] == "4HHB"

    def test_the_deposition_date_comes_with_it(self, deposited: Path) -> None:
        got = structure_provenance(str(deposited), "pdb_file", deposited)
        assert got["deposited"] == "07-MAR-84"

    def test_a_file_with_no_header_claims_no_entry(self, built_by_hand: Path) -> None:
        """A structure built by hand is not a deposited one, and guessing an
        entry from a four-character filename is how a methods section comes
        to state a false provenance."""
        got = structure_provenance(str(built_by_hand), "pdb_file", built_by_hand)
        assert "entry" not in got and "deposited" not in got

    def test_an_mmcif_entry_id_is_read_too(self, tmp_path: Path) -> None:
        p = tmp_path / "s.cif"
        p.write_text("data_4HHB\n#\n_entry.id   4HHB\n#\n")
        assert structure_provenance(str(p), "pdb_file", p)["entry"] == "4HHB"

    def test_the_search_is_bounded(self, tmp_path: Path) -> None:
        """Only the header is read, so the cost does not follow the size of
        the structure. A HEADER line a megabyte in is not a header."""
        p = tmp_path / "big.pdb"
        p.write_text(MINI * 40_000 + HEADER)
        assert "entry" not in structure_provenance(str(p), "pdb_file", p)


class TestWhatWasGivenAndWhatItResolvedTo:
    def test_a_local_file_records_the_absolute_path(self, deposited: Path) -> None:
        """The given path is resolved against a working directory nobody
        records and which is not where the run is read from later."""
        got = structure_provenance("./" + deposited.name, "pdb_file", deposited)
        assert got["given"] == "./" + deposited.name
        assert Path(got["path"]).is_absolute()

    def test_a_fetch_records_the_url_it_used(self, deposited: Path) -> None:
        got = structure_provenance("1ubq", "pdb_id", deposited)
        assert got["url"] == "https://files.rcsb.org/download/1UBQ.pdb"

    def test_when_is_recorded(self, deposited: Path) -> None:
        """Deposited entries are revised. A run made before a revision used
        different coordinates from one made after it, and the identifier
        alone does not say which side of it a run falls on."""
        got = structure_provenance("1ubq", "pdb_id", deposited)
        assert got["retrieved_at"].startswith("20")

    def test_no_file_is_recorded_as_no_record(self, tmp_path: Path) -> None:
        """A record of nothing is worse than no record: it reads as a
        structure that was there."""
        assert structure_provenance("x", "pdb_file", tmp_path / "absent.pdb") is None

    def test_described_reads_as_a_sentence(self, deposited: Path) -> None:
        said = described_structure(structure_provenance(str(deposited), "pdb_file", deposited))
        assert "4HHB" in said and "sha256" in said

    def test_described_says_nothing_about_nothing(self) -> None:
        assert described_structure(None) is None


class TestTheManifestCarriesIt:
    @staticmethod
    def _run(tmp_path: Path, system: str) -> dict:
        orch = MagicMock()
        orch.system = system
        orch.output_dir = tmp_path / "proj"
        orch.output_dir.mkdir()
        orch._presenter = None
        out = orch.output_dir / "setup"
        out.mkdir()
        try:
            setup_run(orchestrator=orch, output_dir=out)
        except Exception:  # noqa: BLE001
            # Preparing a two-atom alanine is allowed to fail, and does on
            # any machine with OpenMM installed: PDBFixer invents the rest
            # of the residue and no template matches the result. That is
            # fine, because it is not what is under test. The record is
            # taken when the structure arrives, and the pipeline writes the
            # manifest before raising -- so the property holds whatever
            # became of preparation, which is the property.
            pass
        return json.loads((out / "setup_parameters.json").read_text(encoding="utf-8"))

    def test_a_run_records_which_structure_it_used(
            self, tmp_path: Path, deposited: Path) -> None:
        manifest = self._run(tmp_path, str(deposited))
        structure = manifest["input"]["structure"]
        assert structure["entry"] == "4HHB"
        assert len(structure["sha256"]) == 64

    def test_the_fields_the_report_reads_are_untouched(
            self, tmp_path: Path, deposited: Path) -> None:
        """The report and the gui both read `system` and `form`; this record
        is added beside them, not in place of them."""
        manifest = self._run(tmp_path, str(deposited))
        assert manifest["input"]["system"] == str(deposited)
        assert manifest["input"]["form"] == "pdb_file"

    def test_a_run_that_fails_in_preparation_still_records_it(
            self, tmp_path: Path, deposited: Path, monkeypatch) -> None:
        """The failures where "which structure was this" is the first
        question are the failures. Stage 1 wrote the manifest before
        raising; stage 3 caught only a missing dependency, so a resolved
        structure that failed to prepare left no manifest at all.

        The failure is forced rather than provoked, so this holds on a
        machine without OpenMM too -- the two differed on it for a week."""
        import fastmdxplora.setup.prepare as prepare_mod

        def refuse(*args, **kwargs):
            raise ValueError("No template found for residue 0 (ALA).")

        monkeypatch.setattr(prepare_mod, "prepare_system", refuse)

        orch = MagicMock()
        orch.system = str(deposited)
        orch.output_dir = tmp_path / "proj"
        orch.output_dir.mkdir()
        orch._presenter = None
        out = orch.output_dir / "setup"
        out.mkdir()
        with pytest.raises(ValueError, match="No template"):
            # fixed_pdb skips PDBFixer, so stage 3 is reached on a machine
            # without it too.
            setup_run(orchestrator=orch, output_dir=out,
                      fixed_pdb=str(deposited))

        manifest = json.loads(
            (out / "setup_parameters.json").read_text(encoding="utf-8"))
        assert manifest["input"]["structure"]["entry"] == "4HHB"
        assert any("preparation failed" in note for note in manifest["notes"])

    def test_a_sequence_records_no_structure(self, tmp_path: Path) -> None:
        """Nothing was read, so there is nothing to describe -- and saying so
        is different from the field being absent because nobody set it."""
        manifest = self._run(tmp_path, "MKTAYIAKQRQISFVKSHFSRQ")
        assert manifest["input"]["form"] == "sequence"
        assert manifest["input"]["structure"] is None


class TestTheMethodsSectionSaysSomethingCheckable:
    """The reason the record is worth keeping: this sentence is printed in a
    paper, and it used to name a file on somebody's laptop."""

    @staticmethod
    def _sentence(setup: dict) -> str:
        text = methods_paragraphs(Path("."), setup, {})
        return text[text.find("Starting coordinates"):][:320]

    def test_a_local_file_names_the_entry_its_header_claims(self) -> None:
        said = self._sentence({"input": {
            "system": "4hhb_cleaned.pdb", "form": "pdb_file",
            "structure": {"form": "pdb_file", "entry": "4HHB",
                          "deposited": "07-MAR-84", "sha256": "d" * 64,
                          "retrieved_at": "2026-08-11T18:03:22+00:00"}}})
        assert "4hhb_cleaned.pdb" in said
        assert "entry 4HHB" in said and "07-MAR-84" in said

    def test_it_is_attributed_rather_than_asserted(self) -> None:
        """A prepared structure keeps the header of the entry it began as.
        That is what makes this useful and also what stops it being a claim
        that the file is that entry."""
        said = self._sentence({"input": {
            "system": "prepped.pdb", "form": "pdb_file",
            "structure": {"form": "pdb_file", "entry": "4HHB",
                          "sha256": "d" * 64, "retrieved_at": "2026-08-11T00:00:00+00:00"}}})
        assert "header identifies it as" in said

    def test_a_fetched_entry_says_when(self) -> None:
        said = self._sentence({"input": {
            "system": "1UBQ", "form": "pdb_id",
            "structure": {"form": "pdb_id", "sha256": "9" * 64,
                          "retrieved_at": "2026-08-11T18:03:22+00:00"}}})
        assert "retrieved 2026-08-11" in said

    def test_the_digest_reaches_the_reader(self) -> None:
        said = self._sentence({"input": {
            "system": "built.pdb", "form": "pdb_file",
            "structure": {"form": "pdb_file", "sha256": "abc123def456" + "0" * 52,
                          "retrieved_at": "2026-08-11T18:03:22+00:00"}}})
        assert "abc123def456" in said

    def test_a_manifest_from_before_this_still_reads(self) -> None:
        """Studies already on disk have no such record, and a report built
        from one must not be worse than it was."""
        said = self._sentence({"input": {"system": "1UBQ", "form": "pdb_id"}})
        assert said.startswith(
            "Starting coordinates were taken from the Protein Data Bank entry 1UBQ.")

    def test_no_structure_block_does_not_invent_one(self) -> None:
        said = self._sentence({"input": {"system": "x.pdb", "form": "pdb_file",
                                         "structure": None}})
        assert "SHA-256" not in said
