"""The parts of the cross-tool comparison that run without the reference.

Its ProLIF-driven commands cannot execute here: the reference tools are
deliberately not dependencies, and there is no finished run in a test
environment. What can be exercised is everything the ten rounds of
hardening went into -- finding artifacts whose recorded paths belong to
another machine, reading a table whose format nobody wrote down, and
turning per-pair rows into a residue-level bracket. Each of those was a
cluster-side failure, and each is a pure function.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from fastmdxplora.validation.cross_tool import (
    _numeric_column,
    _occupancy_from_csv,
    find_artifact,
    load_manifest,
    our_kind_family,
    residue_label,
    trajectory_and_topology,
)


class TestJoiningTwoToolsTables:
    """Neither tool writes a residue the way the other does."""

    @pytest.mark.parametrize("raw,expected", [
        ("ASP189", "ASP189"), ("ASP 189 A", "ASP189"), ("ASP189.A", "ASP189"),
        ("ASP189-OD1", "ASP189"), ("asp189", "ASP189"), ("TRP215.A", "TRP215"),
    ])
    def test_every_decoration_reduces_to_one_label(self, raw, expected):
        assert residue_label(raw) == expected

    def test_something_unrecognisable_is_passed_through(self):
        """Rather than dropped, which would silently shorten a table."""
        assert residue_label("LIG") == "LIG"

    @pytest.mark.parametrize("kind,family", [
        ("hbond_donor", "hbond"), ("HydrogenBond", "hbond"),
        ("salt_bridge", "salt_bridge"), ("Ionic", "salt_bridge"),
        ("pi_stacking", "pi_stacking"), ("Hydrophobic", "hydrophobic"),
        ("XBond_halogen", "halogen"), ("MetalDonor", "metal"),
    ])
    def test_kinds_map_onto_families(self, kind, family):
        assert our_kind_family(kind) == family

    def test_an_unfamiliar_kind_is_excluded_rather_than_guessed(self):
        """A water bridge is not any of the families being compared.

        Mapping it onto the nearest one would put a row in the table that
        the other tool never produced.
        """
        assert our_kind_family("water_bridge") is None
        assert our_kind_family("something novel") is None


class TestFindingArtifactsAfterTheFerry:
    """A manifest records the cluster's absolute paths.

    They exist nowhere on the laptop the run was copied to, which is how
    the first real comparison failed. An absolute path is re-rooted under
    the local directory by everything after the run's own name in it.
    """

    def _run(self, tmp_path: Path) -> Path:
        run = tmp_path / "benchmark-3ptb-bound"
        (run / "simulation").mkdir(parents=True)
        (run / "simulation" / "production.dcd").write_bytes(b"x")
        (run / "setup").mkdir()
        (run / "setup" / "topology.pdb").write_text("ATOM\n", encoding="utf-8")
        return run

    def test_a_cluster_path_is_re_rooted(self, tmp_path):
        run = self._run(tmp_path)
        manifest = {"artifacts": {"traj": (
            "/scratch/aaina/projects/x/benchmark-3ptb-bound/simulation/"
            "production.dcd")}}
        found = find_artifact(run, manifest, "production", ".dcd")

        assert found is not None
        assert found == run / "simulation" / "production.dcd"

    def test_a_relative_path_is_tried_as_given(self, tmp_path):
        run = self._run(tmp_path)
        manifest = {"traj": "simulation/production.dcd"}
        assert find_artifact(run, manifest, ".dcd") == (
            run / "simulation" / "production.dcd")

    def test_a_recorded_path_that_does_not_exist_is_not_returned(self, tmp_path):
        """So the caller falls through to the conventional layout.

        Returning a path that is merely recorded would fail later, further
        from the cause.
        """
        run = self._run(tmp_path)
        manifest = {"traj": "/elsewhere/gone/production.dcd"}
        assert find_artifact(run, manifest, "gone") is None

    def test_all_needles_must_match(self, tmp_path):
        run = self._run(tmp_path)
        manifest = {"a": "simulation/production.dcd"}
        assert find_artifact(run, manifest, "production", ".pdb") is None

    def test_the_conventional_layout_is_the_fallback(self, tmp_path, capsys):
        run = self._run(tmp_path)
        (run / "manifest.json").write_text(json.dumps({}), encoding="utf-8")
        traj, top = trajectory_and_topology(run, load_manifest(run))

        assert traj.name == "production.dcd"
        assert top.name == "topology.pdb"
        assert "fallback" in capsys.readouterr().out

    def test_a_directory_without_a_manifest_says_so(self, tmp_path):
        with pytest.raises(SystemExit, match="run directory"):
            load_manifest(tmp_path)


class TestReadingATableNobodyWroteDown:
    """The interaction table is `.dat`, comment-led, and pre-aggregated.

    Each of those was discovered by a failure: the extension is not `.csv`,
    the header may be behind a `#`, the numbers may be in scientific
    notation, and the rows are per atom pair rather than per frame.
    Reading it as per-frame records printed 100.0 for every pair that
    existed at all, which was the whole of the first deltas table.
    """

    def _write(self, tmp_path: Path, text: str) -> Path:
        path = tmp_path / "pl_interactions.dat"
        path.write_text(text, encoding="utf-8")
        return path

    def test_an_aggregated_table_becomes_a_residue_bracket(self, tmp_path):
        path = self._write(tmp_path, (
            "kind,ligand_atom,protein_atom,frames_present,frames_total,"
            "occupancy,episodes,standard_error,well_sampled,residue\n"
            "hydrophobic,3224,2724,804,2000,0.402,281,0.029,True,VAL213\n"
            "hydrophobic,3225,2724,601,2000,0.3005,279,0.027,True,VAL213\n"
            "hydrophobic,3223,2728,8,2000,0.004,5,0.028,True,VAL213\n"))
        occupancy = _occupancy_from_csv(path)

        value = occupancy[("VAL213", "hydrophobic")]
        low = value[0] if isinstance(value, tuple) else value
        assert low == pytest.approx(40.2, abs=0.1), (
            "the floor is the largest single pair, not the sum")

    def test_pairs_on_one_protein_atom_do_not_add(self, tmp_path):
        """Ring carbons touching one atom fire in the same frames.

        Summing them convicted a passing negative control at 24.8% where
        the union was about 3%.
        """
        path = self._write(tmp_path, (
            "kind,ligand_atom,protein_atom,frames_present,frames_total,"
            "occupancy,episodes,standard_error,well_sampled,residue\n"
            + "".join(
                f"hydrophobic,{3220 + i},2724,60,2000,0.03,50,0.02,True,LEU118\n"
                for i in range(6))))
        value = _occupancy_from_csv(path)[("LEU118", "hydrophobic")]
        high = value[1] if isinstance(value, tuple) else value

        assert high == pytest.approx(3.0, abs=0.2), (
            "six pairs on one atom are one contact, not six")


class TestPickingTheColumnThatHoldsTheNumbers:
    def test_a_headerless_column_is_read(self, tmp_path):
        path = tmp_path / "rmsd.dat"
        path.write_text("1.0e-01\n2.0e-01\n3.0e-01\n", encoding="utf-8")
        values = _numeric_column(path, "rmsd")
        assert len(values) == 3
        assert values[0] == pytest.approx(0.1)

    def test_a_commented_header_is_not_taken_for_data(self, tmp_path):
        """`#` leads the header, and its tokens are not floats.

        Detecting the header by trying to parse the first line is what
        makes scientific notation dangerous: 'e' looks like a name.
        """
        path = tmp_path / "rg.dat"
        path.write_text("# frame rg\n0 1.67\n1 1.68\n", encoding="utf-8")
        values = _numeric_column(path, "rg")
        assert len(values) == 2


class TestReadingPerFrameRecords:
    """The other shape the interaction record can arrive in.

    A JSON list of per-frame contacts is the shape a future per-frame
    export gives, and it is the one where a residue's occupancy is exact
    rather than bracketed, because the frames are there to be unioned.
    """

    def test_frames_are_unioned_per_residue(self):
        from fastmdxplora.validation.cross_tool import _occupancy_from_records

        records = (
            [{"frame": f, "kind": "hydrophobic", "residue": "LEU118"}
             for f in range(0, 40)]
            + [{"frame": f, "kind": "hydrophobic", "residue": "LEU118"}
               for f in range(20, 60)]
            + [{"frame": f, "kind": "hbond", "residue": "ASP189"}
               for f in range(0, 100)]
        )
        occupancy = _occupancy_from_records(records)

        # 100 distinct frames appear; LEU118 is present in 60 of them.
        assert occupancy[("LEU118", "hydrophobic")] == pytest.approx(60.0)
        assert occupancy[("ASP189", "hbond")] == pytest.approx(100.0)

    def test_records_under_a_key_are_found(self):
        from fastmdxplora.validation.cross_tool import _occupancy_from_records

        payload = {"contacts": [
            {"frame": 0, "kind": "salt_bridge", "residue": "ASP189"},
            {"frame": 1, "kind": "salt_bridge", "residue": "ASP189"},
        ]}
        assert _occupancy_from_records(payload)[
            ("ASP189", "salt_bridge")] == pytest.approx(100.0)

    def test_a_shape_it_does_not_recognise_yields_nothing(self):
        """So the caller tries the next candidate file.

        Returning a partial reading of an unfamiliar structure would put
        numbers in the table that no file actually contains.
        """
        from fastmdxplora.validation.cross_tool import _occupancy_from_records

        assert _occupancy_from_records({"summary": {"mean": 1.0}}) == {}
        assert _occupancy_from_records([{"kind": "hbond"}]) == {}
        assert _occupancy_from_records([]) == {}


class TestFindingTheLigand:
    def test_it_reads_the_resolved_configuration(self, tmp_path):
        import yaml

        from fastmdxplora.validation.cross_tool import ligand_resname

        (tmp_path / "resolved_config.yml").write_text(
            yaml.safe_dump({"setup": {"ligand_name": "BEN"}}),
            encoding="utf-8")
        assert ligand_resname(tmp_path) == "BEN"

    def test_without_one_it_says_what_to_pass(self, tmp_path):
        from fastmdxplora.validation.cross_tool import ligand_resname

        with pytest.raises(SystemExit, match="--ligand"):
            ligand_resname(tmp_path)
