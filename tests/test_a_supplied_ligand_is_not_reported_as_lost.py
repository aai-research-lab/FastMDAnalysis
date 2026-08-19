"""A ligand you supply yourself is re-added, and the log must say so.

`pdbfix` splits heterogens into kept and discarded, and warns only about the
discarded: warning about one that was re-added "would describe a loss that did
not happen". The split reads `_reinstated_heterogens`, which only
`_auto_ligands` populated -- and `_auto_ligands` is skipped when the caller
passes `ligand:` themselves. So the offline path, where an SDF is handed over
by hand because there is no network, got:

    WARNING  Removed heterogens: BEN. Pass --setup-keep-heterogens ...
    INFO     Loaded ligand BEN from BEN_ideal.sdf (net charge=0 ...)

four seconds apart, about the same component.
"""
import pytest

from fastmdxplora.setup.pipeline import _explicit_ligand_resnames


def test_a_named_ligand_is_recorded_as_reinstated():
    params = {"ligand": "BEN_ideal.sdf", "ligand_name": "BEN"}
    assert _explicit_ligand_resnames(params) == ("BEN",)


def test_several_ligands_are_all_recorded():
    params = {"ligand": ["a.sdf", "b.sdf"], "ligand_name": ["BEN", "BNZ"]}
    assert _explicit_ligand_resnames(params) == ("BEN", "BNZ")


def test_the_component_code_is_normalised():
    """pdbfix compares against upper-cased names."""
    assert _explicit_ligand_resnames(
        {"ligand": "x.sdf", "ligand_name": "ben"}) == ("BEN",)


@pytest.mark.parametrize("named", [None, "", []])
def test_an_unnamed_ligand_vouches_for_nothing(named):
    """Silence beats a guess: an unnamed file could be anything, and
    claiming a component was kept would suppress a warning about one that
    really was dropped."""
    assert _explicit_ligand_resnames(
        {"ligand": "x.sdf", "ligand_name": named}) == ()


def test_the_recorded_names_satisfy_the_pdbfix_split():
    """The end the fix exists for: what is recorded must be what pdbfix
    checks against, which is an upper-cased set membership test."""
    recorded = _explicit_ligand_resnames(
        {"ligand": "BEN_ideal.sdf", "ligand_name": "BEN"})
    accounted = set(recorded)
    removed = {"BEN": 1, "EPE": 2}
    discarded = {n: c for n, c in removed.items() if n.upper() not in accounted}
    assert "BEN" not in discarded
    assert discarded == {"EPE": 2}
