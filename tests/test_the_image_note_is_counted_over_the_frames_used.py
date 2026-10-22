"""The image check's note on an analysis is counted over the frames that
analysis uses: those from its ``start`` (RMSF, ``ligand_rmsf``,
``bfactor_comparison``), the energy rows each thermodynamic mean averages,
and, for an analysis without a main mean, the frames its means use.

An RMSF over frames 149 to 299, none of which came within the cutoff, said
"100 of 300 frames used (33%)" and that its shape may be biased.
"""

from __future__ import annotations

import json

import mdtraj as md
import numpy as np
import pytest

import fastmdxplora.analysis  # noqa: F401 -- registers the analyses


def _beads(n: int = 300, reaching: int = 100, box: float = 3.0) -> md.Trajectory:
    """Ten alpha carbons bonded in a chain, 10 ps a frame: stretched along x
    (its ends 0.12 nm from each other's image) in the first ``reaching``
    frames, a compact coil after."""
    rng = np.random.default_rng(3)
    top = md.Topology()
    chain = top.add_chain()
    previous = None
    for i in range(10):
        atom = top.add_atom("CA", md.element.carbon,
                            top.add_residue("ALA", chain, resSeq=i + 1))
        if previous is not None:
            top.add_bond(previous, atom)
        previous = atom
    t = np.arange(10) * 1.7
    coil = np.c_[0.25 * np.cos(t), 0.25 * np.sin(t), 0.15 * np.arange(10)]
    straight = np.c_[np.arange(10) * 0.32, np.zeros(10), np.zeros(10)]
    xyz = np.asarray([(straight if f < reaching else coil)
                      + rng.normal(scale=0.01, size=(10, 3)) + 1.0 for f in range(n)])
    traj = md.Trajectory(xyz.astype(np.float32), top, time=np.arange(n) * 10.0)
    traj.unitcell_vectors = np.repeat(np.diag([box] * 3)[None], n, axis=0).astype(np.float32)
    return traj


def _orchestrated(folder, traj, include, options=None):
    from fastmdxplora.analysis import AnalysisOrchestrator

    traj[0].save_pdb(str(folder / "top.pdb"))
    traj.save_dcd(str(folder / "traj.dcd"))
    results = AnalysisOrchestrator(str(folder / "traj.dcd"), str(folder / "top.pdb"),
                                   output_dir=folder / "analysis",
                                   saving_interval_ps=10.0).run(include=include,
                                                                options=options or {})
    assert all(r.status == "ok" for r in results.values()), {
        n: r.message for n, r in results.items()}
    return folder / "analysis"


def _findings(folder, name):
    return json.loads((folder / name / "options.json").read_text())["findings"]


class TestFromItsStart:
    def test_rmsf_from_past_the_frames_that_met_it_carries_no_note(self, tmp_path):
        folder = _orchestrated(tmp_path, _beads(), ["rmsf", "rg"],
                               {"rmsf": {"start": 1.5}})
        assert _findings(folder, "image_check")["solute_image"]["frames_under_cutoff"] == 100
        found = _findings(folder, "rmsf")
        assert found["frames"]["first_frame"] > 100
        assert "solute_image" not in found

    def test_rmsf_from_among_them_counts_from_its_start(self, tmp_path):
        folder = _orchestrated(tmp_path, _beads(), ["rmsf"], {"rmsf": {"start": 0.5}})
        found = _findings(folder, "rmsf")
        first = found["frames"]["first_frame"]
        assert 0 < first < 100
        assert found["solute_image"]["frames_used"] == 300 - first
        assert found["solute_image"]["frames_under_cutoff"] == 100 - first

    @pytest.mark.parametrize("name", ["rmsf", "ligand_rmsf", "bfactor_comparison"])
    def test_each_reads_the_frames_its_record_says(self, name):
        """The three set their first frame in `findings.frames` (from
        `starting_frame.first_frame`), and the note reads it there."""
        from fastmdxplora.analysis.image_check import solute_image_record, solute_atoms
        from fastmdxplora.analysis.orchestrator import get_analysis_class

        traj = _beads()
        series, finding = solute_image_record(traj, solute_atoms(traj.topology, None)[0])
        analysis = get_analysis_class(name)(
            **({"ligand_resname": "LIG"} if name == "ligand_rmsf" else {}))
        analysis.solute_image, analysis.solute_image_series = finding, series
        analysis.findings["frames"] = {"first_frame": 149, "n_frames_given": 300}
        assert analysis._image_note_for(traj, None) is None
        analysis.findings["frames"] = {"first_frame": 60, "n_frames_given": 300}
        note = analysis._image_note_for(traj, None)
        assert note["frames_used"] == 240 and note["frames_under_cutoff"] == 40
