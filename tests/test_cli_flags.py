"""Tests for the per-phase CLI flag plumbing.

These verify that:
  1. Each per-phase subcommand exposes its phase's flags.
  2. Flag values flow into the phase's run() kwargs (not silently dropped).
  3. The `explore` subcommand exposes prefixed equivalents and routes them
     to the right phase via the orchestrator's options dict.
  4. Unset flags do NOT appear as None in kwargs (which would override
     the phase's own defaults).
"""

from __future__ import annotations

import argparse

import pytest

from fastmdxplora.cli.main import (
    _PHASE_SPEC,
    _build_parser,
    _harvest_phase_options,
    _normalize_analysis_options,
)


# ===========================================================================
# Parser surface — every phase flag is actually registered
# ===========================================================================
class TestPerPhaseFlagsExist:
    """Each per-phase subcommand should expose every flag from its options list."""

    @pytest.mark.parametrize("phase,opts_attr,prefix", [
        ("setup",    "_SETUP_OPTIONS",      ""),
        ("simulate", "_SIMULATION_OPTIONS", ""),
        ("analyze",  "_ANALYSIS_OPTIONS",   ""),
        ("report",   "_REPORT_OPTIONS",     ""),
    ])
    def test_phase_has_all_its_flags(self, phase, opts_attr, prefix):
        parser = _build_parser()
        # Argparse's introspection: get the subparser action then look up
        # the subcommand's parser.
        subparsers_action = next(
            a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
        )
        sub_parser = subparsers_action.choices[phase]
        registered = {a.option_strings[0] for a in sub_parser._actions
                      if a.option_strings}

        opts_list = _PHASE_SPEC[phase][0]
        for cli_suffix, _, _ in opts_list:
            expected = f"--{cli_suffix}"
            assert expected in registered, (
                f"{phase} subcommand missing flag {expected!r}; "
                f"registered: {sorted(registered)}"
            )


class TestExploreHasPrefixedFlags:
    """The `explore` subcommand should expose every phase's flag under a prefix."""

    def test_explore_has_setup_prefixed_flags(self):
        parser = _build_parser()
        subparsers_action = next(
            a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
        )
        explore = subparsers_action.choices["explore"]
        registered = {a.option_strings[0] for a in explore._actions
                      if a.option_strings}

        # Spot-check: setup options should have --setup- prefix
        for sample_setup_flag in ("--setup-ph", "--setup-keep-water"):
            assert sample_setup_flag in registered

    def test_explore_has_simulation_prefixed_flags(self):
        parser = _build_parser()
        subparsers_action = next(
            a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
        )
        explore = subparsers_action.choices["explore"]
        registered = {a.option_strings[0] for a in explore._actions
                      if a.option_strings}

        for flag in ("--simulate-duration-ns", "--simulate-platform",
                     "--simulate-nvt-duration-ns", "--simulate-npt-duration-ns"):
            assert flag in registered

    def test_explore_has_analysis_prefixed_flags(self):
        parser = _build_parser()
        subparsers_action = next(
            a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
        )
        explore = subparsers_action.choices["explore"]
        registered = {a.option_strings[0] for a in explore._actions
                      if a.option_strings}

        for flag in ("--analyze-analyses", "--analyze-selection", "--analyze-stride"):
            assert flag in registered

    def test_explore_has_report_prefixed_flags(self):
        parser = _build_parser()
        subparsers_action = next(
            a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
        )
        explore = subparsers_action.choices["explore"]
        registered = {a.option_strings[0] for a in explore._actions
                      if a.option_strings}

        for flag in ("--report-title", "--report-no-slides"):
            assert flag in registered


# ===========================================================================
# Harvest — values flow into kwargs and Nones are dropped
# ===========================================================================
class TestHarvestPhaseOptions:
    def test_set_values_flow_through(self):
        parser = _build_parser()
        args = parser.parse_args([
            "simulate", "-system", "/tmp/x.pdb",
            "--duration-ns", "50.0",
            "--platform", "CUDA",
            "--temperature-K", "310.0",
        ])
        from fastmdxplora.cli.main import _SIMULATION_OPTIONS
        kwargs = _harvest_phase_options(args, _SIMULATION_OPTIONS)
        assert kwargs["duration_ns"] == 50.0
        assert kwargs["platform"] == "CUDA"
        assert kwargs["temperature_K"] == 310.0

    def test_unset_values_dropped_not_none(self):
        """Argparse defaults to None for unset args. Those must NOT reach run().

        If they did, the phase's own DEFAULTS table would be overridden with
        None and every phase would crash on type coercion.
        """
        parser = _build_parser()
        args = parser.parse_args(["simulate", "-system", "/tmp/x.pdb",
                                  "--duration-ns", "10.0"])
        from fastmdxplora.cli.main import _SIMULATION_OPTIONS
        kwargs = _harvest_phase_options(args, _SIMULATION_OPTIONS)
        # Only duration_ns was supplied; everything else must be absent.
        assert kwargs == {"duration_ns": 10.0}
        assert "platform" not in kwargs
        assert "temperature_K" not in kwargs

    def test_explore_routes_to_correct_phase(self):
        """--setup-ph lands in options['setup']; --simulate-duration-ns in options['simulation']."""
        parser = _build_parser()
        args = parser.parse_args([
            "explore", "-system", "/tmp/x.pdb",
            "--setup-ph", "6.5",
            "--simulate-duration-ns", "100.0",
            "--analyze-selection", "name CA",
            "--report-title", "Test Run",
        ])

        # Mimic what _cmd_explore does
        options: dict = {}
        from fastmdxplora.cli.main import _PHASE_SPEC, _PHASE_TO_ORCH
        for phase, (opts, _) in _PHASE_SPEC.items():
            h = _harvest_phase_options(args, opts, dest_prefix=phase)
            if h:
                options[_PHASE_TO_ORCH[phase]] = h

        assert options["setup"]["ph"] == 6.5
        assert options["simulation"]["duration_ns"] == 100.0
        assert options["analysis"]["selection"] == "name CA"
        assert options["report"]["title"] == "Test Run"

    def test_explore_with_nothing_supplied_gives_empty_options(self):
        parser = _build_parser()
        args = parser.parse_args(["explore", "-system", "/tmp/x.pdb"])

        options: dict = {}
        from fastmdxplora.cli.main import _PHASE_SPEC, _PHASE_TO_ORCH
        for phase, (opts, _) in _PHASE_SPEC.items():
            h = _harvest_phase_options(args, opts, dest_prefix=phase)
            if h:
                options[_PHASE_TO_ORCH[phase]] = h
        assert options == {}


# ===========================================================================
# Boolean flags — store_false flags should land in kwargs when set
# ===========================================================================
class TestBooleanFlags:
    def test_no_minimize_lands_as_false(self):
        """--no-minimize sets minimize=False (action=store_false)."""
        parser = _build_parser()
        args = parser.parse_args(["simulate", "-system", "/tmp/x.pdb",
                                  "--no-minimize"])
        from fastmdxplora.cli.main import _SIMULATION_OPTIONS
        kwargs = _harvest_phase_options(args, _SIMULATION_OPTIONS)
        assert kwargs["minimize"] is False

    def test_no_minimize_absent_means_kwarg_absent(self):
        """When --no-minimize is NOT passed, minimize should NOT be in kwargs."""
        parser = _build_parser()
        args = parser.parse_args(["simulate", "-system", "/tmp/x.pdb"])
        from fastmdxplora.cli.main import _SIMULATION_OPTIONS
        kwargs = _harvest_phase_options(args, _SIMULATION_OPTIONS)
        # Not passed → not in kwargs → phase falls through to its default True
        assert "minimize" not in kwargs

    def test_report_no_slides(self):
        parser = _build_parser()
        args = parser.parse_args(["report", "-system", "/tmp/x.pdb",
                                  "--no-slides"])
        from fastmdxplora.cli.main import _REPORT_OPTIONS
        kwargs = _harvest_phase_options(args, _REPORT_OPTIONS)
        assert kwargs["slides"] is False


# ===========================================================================
# Multi-valued flags (nargs='+')
# ===========================================================================
class TestMultiValuedFlags:
    def test_analyses_subset(self):
        parser = _build_parser()
        args = parser.parse_args([
            "analyze", "-system", "/tmp/x.pdb",
            "--analyses", "rmsd", "rmsf", "rg",
        ])
        from fastmdxplora.cli.main import _ANALYSIS_OPTIONS
        kwargs = _harvest_phase_options(args, _ANALYSIS_OPTIONS)
        assert kwargs["include"] == ["rmsd", "rmsf", "rg"]

    def test_force_field_multi(self):
        parser = _build_parser()
        args = parser.parse_args([
            "setup", "-system", "/tmp/x.pdb",
            "--force-field", "amber14-all.xml", "amber14/tip3pfb.xml",
        ])
        from fastmdxplora.cli.main import _SETUP_OPTIONS
        kwargs = _harvest_phase_options(args, _SETUP_OPTIONS)
        assert kwargs["force_field"] == ["amber14-all.xml", "amber14/tip3pfb.xml"]


class TestAnalysisCompatibility:
    def test_fastmdanalysis_profile_sets_paper_defaults(self):
        parser = _build_parser()
        args = parser.parse_args([
            "analyze", "-system", "/tmp/bpti.pdb",
            "--compat", "fastmdanalysis",
        ])
        from fastmdxplora.cli.main import _ANALYSIS_OPTIONS
        kwargs = _normalize_analysis_options(
            _harvest_phase_options(args, _ANALYSIS_OPTIONS)
        )
        assert kwargs["scope"] == "protein"
        assert kwargs["selection"] == "protein"
        assert kwargs["stride"] == 2
        assert kwargs["include"] == [
            "rmsd", "rmsf", "rg", "hbonds", "sasa", "ss", "dimred", "cluster"
        ]
        assert kwargs["options"]["rmsd"] == {"ref": 0, "align": True}
        assert kwargs["options"]["rmsf"] == {"ref": 0, "per_residue": False}
        assert kwargs["options"]["sasa"] == {
            "mode": "total", "probe_radius": 0.14, "n_sphere_points": 960
        }
        assert kwargs["options"]["hbonds"] == {
            "method": "baker_hubbard",
            "freq": 0.1,
            "candidate_freq": 0.0,
            "count_multiplier": 2,
        }
        assert kwargs["options"]["dimred"] == {"methods": ["pca"], "n_components": 2}
        assert kwargs["options"]["cluster"] == {
            "methods": ["hierarchical"], "n_clusters": 6
        }

    def test_explicit_analysis_values_override_profile(self):
        parser = _build_parser()
        args = parser.parse_args([
            "analyze", "-system", "/tmp/bpti.pdb",
            "--compat", "fastmdanalysis",
            "--stride", "4",
            "--cluster-methods", "kmeans",
            "--cluster-n-clusters", "3",
        ])
        from fastmdxplora.cli.main import _ANALYSIS_OPTIONS
        kwargs = _normalize_analysis_options(
            _harvest_phase_options(args, _ANALYSIS_OPTIONS)
        )
        assert kwargs["stride"] == 4
        assert kwargs["options"]["cluster"] == {
            "methods": ["kmeans"], "n_clusters": 3
        }


# ===========================================================================
# Choices enforcement (argparse-level)
# ===========================================================================
class TestChoicesValidation:
    def test_platform_choice_rejects_invalid(self):
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["simulate", "-system", "/tmp/x.pdb",
                               "--platform", "NotAPlatform"])

    def test_box_shape_choice(self):
        parser = _build_parser()
        # Valid value passes
        args = parser.parse_args(["setup", "-system", "/tmp/x.pdb",
                                  "--box-shape", "dodecahedron"])
        assert args.box_shape == "dodecahedron"
        # Invalid value rejected
        with pytest.raises(SystemExit):
            parser.parse_args(["setup", "-system", "/tmp/x.pdb",
                               "--box-shape", "tetrahedron"])


# ===========================================================================
# Dual-dash input flags — both -system and --system must work
# ===========================================================================
class TestSystemInputFlag:
    """The system input flag accepts three forms: -s, -system, --system.

    -s is the GNU short option; -system is the GROMACS/AMBER/NAMD-style
    single-dash long flag; --system is the GNU double-dash long flag.
    There is no separate --pdb-id flag — PDB IDs go through --system and
    are auto-detected by the setup classifier.
    """

    @pytest.mark.parametrize("flag", ["-s", "-system", "--system"])
    def test_system_all_three_forms(self, flag):
        parser = _build_parser()
        args = parser.parse_args(["setup", flag, "protein.pdb"])
        assert args.system == "protein.pdb"

    def test_pdb_id_passes_through_system(self):
        """A 4-char PDB ID is a valid --system value (no separate flag)."""
        parser = _build_parser()
        args = parser.parse_args(["setup", "--system", "1L2Y"])
        assert args.system == "1L2Y"

    def test_no_pdb_id_flag_exists(self):
        """The --pdb-id flag was removed; --system handles the ID case."""
        parser = _build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["setup", "--pdb-id", "1L2Y"])

    def test_single_dash_system_not_split_into_chars(self):
        """Guard against argparse misreading -system as stacked short flags."""
        parser = _build_parser()
        args = parser.parse_args(["analyze", "-system", "x.pdb"])
        assert args.system == "x.pdb"

    def test_all_forms_on_explore(self):
        """All three forms work on the explore subcommand too."""
        parser = _build_parser()
        a1 = parser.parse_args(["explore", "-s", "x.pdb"])
        a2 = parser.parse_args(["explore", "-system", "x.pdb"])
        a3 = parser.parse_args(["explore", "--system", "x.pdb"])
        assert a1.system == a2.system == a3.system == "x.pdb"
