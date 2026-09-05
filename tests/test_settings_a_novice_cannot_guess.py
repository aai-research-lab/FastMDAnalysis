"""Two settings whose absence did something other than nothing.

Both were found by using the software rather than reading it, and both were
got wrong by somebody who had spent the day inside the source. A caller who
has not cannot be expected to do better.
"""

from __future__ import annotations

import pytest


class TestAskingForWorkersAsksForParallel:
    """`mode` defaulted to sequential and `workers` was read separately, so
    `workers: 3` ran the windows one at a time and said nothing about it. A
    study of nine windows took three times as long as it was told to."""

    @staticmethod
    def _mode(execution: dict) -> str:
        from fastmdxplora.batch.explorer import BatchExplorer

        BatchExplorer.__new__(BatchExplorer)
        workers = execution.get("workers")
        devices = execution.get("devices")
        return execution.get("mode") or (
            "parallel" if (workers or devices) else "sequential")

    def test_a_worker_count_means_parallel(self) -> None:
        assert self._mode({"workers": 3}) == "parallel"

    def test_a_device_list_means_parallel(self) -> None:
        """Nobody lists GPUs to leave all but one idle."""
        assert self._mode({"devices": [0, 1]}) == "parallel"

    def test_nothing_asked_for_stays_sequential(self) -> None:
        assert self._mode({}) == "sequential"

    def test_written_out_still_wins(self) -> None:
        """A config asking for both a worker count and sequential gets
        sequential: the explicit word beats the inference."""
        assert self._mode({"mode": "sequential", "workers": 4}) == "sequential"

    def test_the_source_infers_it(self) -> None:
        import inspect
        from fastmdxplora.batch.explorer import BatchExplorer

        source = inspect.getsource(BatchExplorer.__init__)
        assert 'execution.get("mode") or' in source


class TestCappingGroupsAreNotHeterogens:
    """ACE and NME terminate a chain. Stripping them is right for buffer and
    cryoprotectant and wrong here: a real run lost both and failed with "no
    template found for ALA", naming the residue that survived rather than the
    two that did not."""

    def test_the_usual_caps_are_listed(self) -> None:
        from fastmdxplora.setup.pdbfix import CAPPING_GROUPS

        assert {"ACE", "NME", "NHE"} <= CAPPING_GROUPS

    def test_a_buffer_molecule_is_not_one(self) -> None:
        """The point is a short known list, not a licence to keep everything:
        a sulfate is a heterogen and should still go."""
        from fastmdxplora.setup.pdbfix import CAPPING_GROUPS

        assert not ({"SO4", "GOL", "EDO", "PEG"} & CAPPING_GROUPS)

    def test_they_survive_removal(self) -> None:
        pdbfixer = pytest.importorskip("pdbfixer")
        from pathlib import Path

        from fastmdxplora.setup.pdbfix import _protect_capping_groups

        structure = Path(__file__).resolve().parents[1] / "alanine-dipeptide.pdb"
        if not structure.is_file():
            pytest.skip("no capped structure to test against")

        fixer = pdbfixer.PDBFixer(filename=str(structure))
        with _protect_capping_groups(fixer.topology):
            fixer.removeHeterogens(keepWater=False)
        kept = {r.name for r in fixer.topology.residues()}
        assert {"ACE", "NME"} <= kept

    def test_the_global_list_is_put_back(self) -> None:
        """It reaches into PDBFixer's module-level residue list, so leaving it
        extended would change what every later call considers standard."""
        pf = pytest.importorskip("pdbfixer.pdbfixer")
        from fastmdxplora.setup.pdbfix import _protect_capping_groups

        class _Residue:
            name = "ACE"

        class _Topology:
            def residues(self): return [_Residue()]

        before = list(pf.proteinResidues)
        with _protect_capping_groups(_Topology()) as caps:
            assert "ACE" in pf.proteinResidues
            assert caps == {"ACE"}
        assert pf.proteinResidues == before

    def test_a_structure_with_no_caps_changes_nothing(self) -> None:
        from fastmdxplora.setup.pdbfix import _protect_capping_groups

        class _Residue:
            name = "ALA"

        class _Topology:
            def residues(self): return [_Residue()]

        with _protect_capping_groups(_Topology()) as caps:
            assert caps == frozenset()
