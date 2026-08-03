"""An option a user can set must be settable from every interface.

FastMDXplora offers four routes to the same workflow: the command line, a
configuration file, the Python API, and the browser. An option reachable from
only some of them is a feature that exists for some users and not others, and
`setup.heterogens` shipped that way: the automatic protein-ligand preparation
it controls could not be reached from the GUI at all.

These tests check the routes agree, so the gap fails here rather than
appearing as a missing control someone eventually notices.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from fastmdxplora.config.schema import all_schemas

GUI = Path(__file__).resolve().parents[1] / "src" / "fastmdxplora" / "gui"

#: Options deliberately absent from the browser, with the reason.
GUI_EXEMPT = {
    # Paths and file handles: the browser cannot meaningfully offer a
    # filesystem path on a machine it may not be running on.
    "fixed_pdb", "ligand", "trajectory", "topology",
    # Expert force-field escape hatches, set through a config file.
    "force_field", "ligand_forcefield", "ligand_net_charge", "ligand_name",
    # Batch and reporting concerns outside the builder's scope.
    "region_highlights", "comparison", "options", "include", "exclude",
}



def _parse(argv):
    """Parse a command line the way the CLI does, without running anything."""
    from fastmdxplora.cli.main import _build_parser

    return _build_parser().parse_args(argv)


class TestHeterogensReachesEveryInterface:
    """The option behind 2.2.0's headline feature."""

    def test_the_schema_defines_it(self) -> None:
        setup = all_schemas()["setup"]
        assert "heterogens" in setup.field_names()

    def test_the_cli_exposes_it(self) -> None:
        """Under both spellings: --heterogens on setup, --setup-heterogens on explore."""
        from fastmdxplora.cli.main import _SETUP_OPTIONS

        names = {cli_suffix for cli_suffix, _kwarg, _kw in _SETUP_OPTIONS}
        assert "heterogens" in names

        parsed = _parse(["explore", "--system", "4W52", "--setup-heterogens", "auto"])
        assert getattr(parsed, "setup__heterogens", None) == "auto"

        parsed = _parse(["setup", "--system", "4W52", "--heterogens", "auto"])
        assert getattr(parsed, "heterogens", None) == "auto"

    def test_a_config_file_carries_it(self, tmp_path) -> None:
        from fastmdxplora.config.loader import load_config_file

        config = tmp_path / "study.yml"
        config.write_text(
            "system: 4W52\nsetup:\n  heterogens: auto\n", encoding="utf-8"
        )
        loaded = load_config_file(str(config))
        assert loaded["setup"]["heterogens"] == "auto"

    def test_the_browser_offers_all_three_choices(self) -> None:
        markup = (GUI / "templates" / "dashboard.html").read_text(encoding="utf-8")
        assert 'id="builder-heterogens"' in markup
        for choice in ("drop", "auto", "keep"):
            assert f'value="{choice}"' in markup, f"the {choice} option is missing"

    def test_the_browser_sends_it_and_the_server_validates_it(self) -> None:
        script = (GUI / "static" / "simulation-builder.js").read_text(encoding="utf-8")
        assert "heterogens:" in script, "the builder must send the choice"

        backend = (GUI / "exploration.py").read_text(encoding="utf-8")
        assert '"drop", "keep", "auto"' in backend, "the server must validate it"

    def test_the_gui_puts_it_in_the_command_and_the_saved_config(self) -> None:
        backend = (GUI / "exploration.py").read_text(encoding="utf-8")
        assert "--setup-heterogens" in backend, "the launched command must carry it"
        assert 'setup_block["heterogens"]' in backend, (
            "a saved config must carry it, or a study designed in the browser "
            "cannot be reproduced on a cluster"
        )


class TestDefaultsAgree:
    """One default, however the option is reached."""

    @pytest.mark.parametrize("option", ["heterogens", "ph", "forcefield"])
    def test_the_phase_default_matches_the_schema(self, option) -> None:
        from fastmdxplora.setup.pipeline import DEFAULTS

        schema = {f.name: f.default for f in all_schemas()["setup"].fields}
        assert DEFAULTS[option] == schema[option], (
            f"{option} defaults to {DEFAULTS[option]!r} in the phase but "
            f"{schema[option]!r} in the schema"
        )

    def test_centre_of_mass_motion_is_removed_by_default(self) -> None:
        """OpenMM removes it by default; silently differing would surprise."""
        from fastmdxplora.setup.pipeline import DEFAULTS

        assert DEFAULTS["remove_cm_motion"] is True
        schema = {f.name: f.default for f in all_schemas()["setup"].fields}
        assert schema["remove_cm_motion"] is True


class TestNoOptionIsQuietlyUnreachable:
    """A broad sweep, so the next gap is caught without being looked for."""

    def test_every_builder_setup_option_is_a_real_schema_field(self) -> None:
        """The browser must not offer settings the schema does not know."""
        import re

        script = (GUI / "static" / "simulation-builder.js").read_text(encoding="utf-8")
        # The payload's setup section, which is the object literal following
        # the key. Slice past the key itself so it is not read as a field.
        start = script.index("setup: {") + len("setup: {")
        block = script[start: script.index("},", start)]
        offered = set(re.findall(r"^\s*([a-z_]+):", block, re.M))

        known = set(all_schemas()["setup"].field_names())
        unknown = offered - known
        assert not unknown, (
            f"the builder sends settings the schema does not define: {unknown}"
        )


class TestAutoResolvesConsistently:
    """`auto` must mean one stack, not a choice made per structure."""

    def test_auto_is_ligand_capable(self) -> None:
        """Otherwise the default combination refuses on any bound ligand."""
        from fastmdxplora.setup.forcefields import resolve_forcefield

        assert resolve_forcefield("auto").supports_ligand is True

    def test_auto_does_not_depend_on_the_structure(self) -> None:
        """An apo run and its holo partner must share a protein force field.

        Choosing per-system would give 4W51 CHARMM36 and 4W52 AMBER, and the
        comparison between them, which is usually the whole point, would be
        meaningless.
        """
        from fastmdxplora.setup.forcefields import AUTO_FORCEFIELD, resolve_forcefield

        assert resolve_forcefield("auto").name == AUTO_FORCEFIELD
        assert resolve_forcefield(None).name == AUTO_FORCEFIELD

    def test_every_interface_accepts_auto(self) -> None:
        from fastmdxplora.gui.exploration import _FORCEFIELDS

        assert "auto" in _FORCEFIELDS

        parsed = _parse(["explore", "--system", "X", "--setup-forcefield", "auto"])
        assert parsed.setup__forcefield == "auto"

    def test_naming_a_force_field_still_overrides(self) -> None:
        from fastmdxplora.setup.forcefields import resolve_forcefield

        assert resolve_forcefield("charmm36").name == "charmm36"

    def test_the_defaults_work_together(self) -> None:
        """heterogens=auto needs a force field that can take a ligand."""
        from fastmdxplora.setup.forcefields import resolve_forcefield
        from fastmdxplora.setup.pipeline import DEFAULTS

        # heterogens stays opt-in for now: the classifier keeps coordinated
        # metals that the clash check then rejects, so making it the default
        # would break structures that prepare cleanly today. The force field
        # is nonetheless ligand-capable, so turning the policy on needs no
        # second flag.
        assert DEFAULTS["heterogens"] == "drop"
        assert resolve_forcefield(DEFAULTS["forcefield"]).supports_ligand, (
            "--setup-heterogens auto must work without also naming a force "
            "field, or the feature needs two flags instead of one"
        )
