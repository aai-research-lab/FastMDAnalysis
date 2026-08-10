"""A path that names a structure and finds none says which.

The classifier tested existence and file type in one condition, so a `.pdb`
that was not there fell through to a message saying a `.pdb` was expected --
advising the author to supply what they had just supplied. A path written one
directory off therefore read as a complaint about the file format, and cost a
run start before anyone looked at the path itself.

The path is read from where `fastmdx` was run, not from the config that names
it, which is the thing worth saying and the thing the old message could not.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import fastmdxplora.setup.pipeline as pipeline

_classify = next(v for k, v in vars(pipeline).items()
                 if callable(v) and "classify" in k.lower())


class TestAMissingStructureIsItsOwnFailure:
    @pytest.mark.parametrize("name", ["absent.pdb", "gone.cif", "nope.pdbx"])
    def test_it_raises_file_not_found(self, name: str) -> None:
        with pytest.raises(FileNotFoundError):
            _classify(name)

    def test_it_does_not_ask_for_what_was_given(self) -> None:
        """The old message told a `.pdb` author to supply a `.pdb`."""
        with pytest.raises(FileNotFoundError) as caught:
            _classify("absent.pdb")
        assert "Expected the path to a" not in str(caught.value)

    def test_it_names_the_directory_paths_are_read_from(self) -> None:
        with pytest.raises(FileNotFoundError) as caught:
            _classify("absent.pdb")
        message = str(caught.value)
        assert str(Path.cwd()) in message
        assert "not from the config file" in message


class TestTheOtherAnswersAreUnchanged:
    def test_a_structure_that_exists_is_a_file(self, tmp_path: Path) -> None:
        here = tmp_path / "system.pdb"
        here.write_text("END\n", encoding="utf-8")
        assert _classify(str(here)) == "pdb_file"

    def test_four_characters_are_an_identifier(self) -> None:
        assert _classify("1L2Y") == "pdb_id"

    def test_letters_are_a_sequence(self) -> None:
        """Recognised and refused further on, because building a structure
        from one needs a predictor this software does not carry."""
        assert _classify("ACDEFGHIK") == "sequence"

    def test_something_unrecognisable_still_says_what_is_accepted(self) -> None:
        """That message is right when the input names no format at all."""
        with pytest.raises(ValueError, match="Could not classify"):
            _classify("nonsense!")

    def test_nothing_at_all_is_refused(self) -> None:
        with pytest.raises(ValueError, match="requires a system input"):
            _classify(None)
