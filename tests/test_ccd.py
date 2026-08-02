"""Tests for obtaining the chemistry a crystal structure leaves out.

None of these touch the network. The fetcher is injected, so the tests
describe the module's contract rather than RCSB's availability.
"""

from __future__ import annotations

import pytest

from fastmdxplora.setup import ccd

# Reading chemical definitions and detecting ionizable groups needs RDKit,
# which ships in the optional ligand extra. CI installs the test extra only,
# so these skip there rather than failing.
pytest.importorskip("rdkit", reason="requires the [ligand] extra")


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path, monkeypatch):
    """Never read or write the developer's real cache during tests."""
    monkeypatch.setenv("FASTMDXPLORA_CACHE_DIR", str(tmp_path / "cache"))


def _sdf(smiles: str, name: str) -> str:
    """Build an SDF the way ModelServer would return one."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    molecule = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(molecule, randomSeed=1)
    block = Chem.MolToMolBlock(molecule)
    return f"{name}\n" + block.split("\n", 1)[1] + "$$$$\n"


def _serving(text: str):
    """A fetcher that returns ``text`` and counts how often it was called."""
    calls: list[str] = []

    def fetcher(url: str, timeout: int = 30) -> str:
        calls.append(url)
        return text

    fetcher.calls = calls  # type: ignore[attr-defined]
    return fetcher


class TestRetrieval:
    def test_a_component_is_fetched_at_its_pose(self) -> None:
        fetcher = _serving(_sdf("c1ccccc1", "BNZ"))
        chemistry = ccd.fetch_chemistry(
            "4W52", "BNZ", chain="A", resseq=201, fetcher=fetcher
        )

        assert chemistry.resname == "BNZ"
        assert chemistry.formal_charge == 0
        assert chemistry.from_cache is False
        # The request must identify the copy, not merely the component, or a
        # structure with several copies would get the first one every time.
        url = fetcher.calls[0]
        assert "auth_asym_id=A" in url
        assert "auth_seq_id=201" in url
        assert "encoding=sdf" in url

    def test_a_second_run_needs_no_network(self) -> None:
        """The cache is the offline story; the dictionary is too large to ship."""
        fetcher = _serving(_sdf("c1ccccc1", "BNZ"))
        ccd.fetch_chemistry("4W52", "BNZ", chain="A", resseq=201, fetcher=fetcher)
        again = ccd.fetch_chemistry("4W52", "BNZ", chain="A", resseq=201, fetcher=fetcher)

        assert again.from_cache is True
        assert len(fetcher.calls) == 1

    def test_copies_are_cached_separately(self) -> None:
        """Two copies sit at different poses and are different files."""
        fetcher = _serving(_sdf("c1ccccc1", "BNZ"))
        ccd.fetch_chemistry("4W52", "BNZ", chain="A", resseq=201, fetcher=fetcher)
        ccd.fetch_chemistry("4W52", "BNZ", chain="B", resseq=201, fetcher=fetcher)
        assert len(fetcher.calls) == 2


class TestRefusals:
    def test_being_offline_with_nothing_cached_stops(self) -> None:
        def offline(url: str, timeout: int = 30) -> str:
            raise ccd.ChemistryUnavailableError("network unreachable")

        with pytest.raises(ccd.ChemistryUnavailableError, match="not cached"):
            ccd.fetch_chemistry("9def", "XYZ", chain="A", resseq=9, fetcher=offline)

    def test_an_empty_response_stops(self) -> None:
        with pytest.raises(ccd.ChemistryUnavailableError, match="no usable chemistry"):
            ccd.fetch_chemistry(
                "9def", "XYZ", chain="A", resseq=9,
                fetcher=lambda url, timeout=30: "",
            )

    def test_a_partially_resolved_component_stops(self) -> None:
        """Otherwise the missing atoms would be placed at invented positions."""
        text = _sdf("c1ccccc1", "BNZ")
        with pytest.raises(ccd.ChemistryUnavailableError, match="partially resolved"):
            ccd.fetch_chemistry(
                "9abc", "BNZ", chain="A", resseq=1,
                expected_heavy_atoms=4,
                fetcher=lambda url, timeout=30: text,
            )


class TestProtonation:
    """A ligand's pKa in a pocket is a property of the complex, not the ligand."""

    def _chemistry(self, smiles: str, code: str):
        return ccd.fetch_chemistry(
            "9xyz", code, chain="A", resseq=1,
            fetcher=lambda url, timeout=30: _sdf(smiles, code),
        )

    def test_a_ligand_with_no_ionizable_group_proceeds(self) -> None:
        chemistry = self._chemistry("c1ccccc1", "BNZ")
        assert chemistry.titratable_groups == ()
        ccd.require_determined_protonation(chemistry, 7.0)  # must not raise

    @pytest.mark.parametrize(
        ("smiles", "code", "group"),
        [
            ("CC(=O)Oc1ccccc1C(=O)O", "AIN", "carboxylic acid"),
            ("NCCO", "ETA", "amine"),
            ("c1cnc[nH]1", "IMD", "imidazole"),
            ("OP(=O)(O)OCC1OC(O)C(O)C1O", "R5P", "phosphate"),
        ],
    )
    def test_a_titratable_ligand_stops(self, smiles, code, group) -> None:
        chemistry = self._chemistry(smiles, code)
        assert chemistry.titratable_groups, f"{code} should be flagged"

        with pytest.raises(ccd.ProtonationUndeterminedError) as excinfo:
            ccd.require_determined_protonation(chemistry, 7.4)
        assert group in str(excinfo.value)

    def test_the_refusal_explains_why_the_ligand_alone_cannot_answer(self) -> None:
        """The message has to teach, or the user will just force it through."""
        chemistry = self._chemistry("CC(=O)Oc1ccccc1C(=O)O", "AIN")
        with pytest.raises(ccd.ProtonationUndeterminedError) as excinfo:
            ccd.require_determined_protonation(chemistry, 7.4)
        message = str(excinfo.value)
        assert "binding site" in message
        assert "pKa" in message
        assert "SDF or MOL2" in message


class TestIdentity:
    """The query names a position; the response must still be checked."""

    def test_a_different_component_is_rejected(self) -> None:
        """A wrong residue number returns a real molecule, just the wrong one.

        Without this check the atom-count mismatch that follows reads as
        partial occupancy, which sends the user looking in the wrong place.
        """
        hepes = _sdf("OCCN1CCN(CCS(=O)(=O)O)CC1", "EPE")
        with pytest.raises(ccd.ChemistryUnavailableError, match="different component"):
            ccd.fetch_chemistry(
                "4W52", "BNZ", chain="A", resseq=201,
                fetcher=lambda url, timeout=30: hepes,
            )

    def test_the_request_names_the_component(self) -> None:
        fetcher = _serving(_sdf("c1ccccc1", "BNZ"))
        ccd.fetch_chemistry("4W52", "BNZ", chain="A", resseq=201, fetcher=fetcher)
        assert "label_comp_id=BNZ" in fetcher.calls[0]


class TestHydrogens:
    """ModelServer returns what the experiment resolved: heavy atoms only."""

    RAW_BENZENE = """BNZ
  ModelServer 0.9.13

  6  6  0  0  0  0  0  0  0  0  0
  -32.9690    6.1960    2.8770 C   0  0  0  0  0  0  0  0  0  0  0  0
  -32.9450    7.0460    3.9730 C   0  0  0  0  0  0  0  0  0  0  0  0
  -33.7190    6.7980    5.1130 C   0  0  0  0  0  0  0  0  0  0  0  0
  -34.5400    5.6800    5.1430 C   0  0  0  0  0  0  0  0  0  0  0  0
  -34.5640    4.8300    4.0470 C   0  0  0  0  0  0  0  0  0  0  0  0
  -33.7900    5.0780    2.9070 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END
$$$$
"""

    def test_hydrogens_are_added_to_what_is_retrieved(self) -> None:
        """Otherwise benzene is parameterized as a six-carbon radical.

        This is the exact payload RCSB returned for 4W52 BNZ: six atoms, six
        bonds, no hydrogens. A force field handed that describes a molecule
        that does not exist.
        """
        chemistry = ccd.fetch_chemistry(
            "4W52", "BNZ", chain="A", resseq=200,
            expected_heavy_atoms=6,
            fetcher=lambda url, timeout=30: self.RAW_BENZENE,
        )
        assert chemistry.n_atoms == 12

    def test_the_cached_file_is_chemically_correct(self) -> None:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors

        chemistry = ccd.fetch_chemistry(
            "4W52", "BNZ", chain="A", resseq=200,
            fetcher=lambda url, timeout=30: self.RAW_BENZENE,
        )
        molecule = Chem.MolFromMolFile(str(chemistry.path), removeHs=False)

        assert molecule is not None
        assert rdMolDescriptors.CalcMolFormula(molecule) == "C6H6"
        assert Chem.MolToSmiles(Chem.RemoveHs(molecule)) == "c1ccccc1"
        # A radical here would mean the hydrogens never arrived.
        assert sum(a.GetNumRadicalElectrons() for a in molecule.GetAtoms()) == 0
        # The crystallographic pose must survive the round trip.
        assert molecule.GetNumConformers() == 1

    def test_the_heavy_atom_check_still_sees_the_structure(self) -> None:
        """Validation runs before hydrogens, against what was resolved."""
        with pytest.raises(ccd.ChemistryUnavailableError, match="partially resolved"):
            ccd.fetch_chemistry(
                "4W52", "BNZ", chain="A", resseq=200,
                expected_heavy_atoms=4,
                fetcher=lambda url, timeout=30: self.RAW_BENZENE,
            )
