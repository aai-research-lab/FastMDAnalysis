"""The plumed setting is asked for as a script, and the config carries it.

The dashboard drew `plumed` with the same YAML textarea the other
enhanced-sampling blocks use, which meant writing `script: |` with
block-scalar indentation around a working .dat file -- YAML inside a box,
where one stray tab changes a PLUMED input. And the setting that exists to
name a file on the machine that runs the study was the one path in the form
without the Browse control every other path field has.

Now the schema payload marks the field as a script; the page draws an
on-switch, a path wired to the shared picker, and a place to write the
script itself; and whichever way it arrives, the config carries the one
`script` string the engine already knows how to read (a single line naming
a file that exists is read from disk; anything else is the script).
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from fastmdxplora.config.languages import (
    UntranslatableSetting,
    cli_command,
    python_script,
)
from fastmdxplora.config.schema import PHASE_SCHEMAS
from fastmdxplora.gui.browse import KINDS
from fastmdxplora.gui.config_builder import _coerce, build_config, config_yaml
from fastmdxplora.gui.schema_payload import field_payload

JS = Path("src/fastmdxplora/gui/static/run-builder.js").read_text(
    encoding="utf-8")

PLUMED = {
    "enabled": True,
    "script": "d: DISTANCE ATOMS=1,2\nPRINT ARG=d FILE=colvar\n",
}


def _plumed_field():
    return next(
        f for f in PHASE_SCHEMAS["simulation"].fields if f.name == "plumed")


def _state():
    return {
        "output": "out",
        "include": ["simulation"],
        "system": "protein.pdb",
        "simulation": {"plumed": dict(PLUMED)},
    }


class TestTheFieldIsAskedForAsAScript:
    def test_the_payload_calls_it_a_script(self):
        # The payload is the contract between schema and page: the page keys
        # its drawing on `control`, so this is the fact that turns the YAML
        # textarea into the composite.
        assert field_payload(_plumed_field())["control"] == "script"

    def test_the_picker_knows_what_a_plumed_file_is(self):
        # plumed.dat is the convention; .plumed is what this package itself
        # writes next to a run, so one study's script can be picked up by
        # the next from the same dialog.
        assert ".dat" in KINDS["plumed"]
        assert ".plumed" in KINDS["plumed"]

    def test_the_page_draws_the_composite(self):
        assert 'field.control === "script"' in JS
        # The path input is marked for the shared picker rather than wired
        # to a second browsing mechanism; attachAll() finds it by this mark
        # after every render. The kind and the id are derived from the
        # field, so no setting is named in the page -- the rule the rest of
        # the form lives by.
        assert "path.dataset.picks = field.name" in JS
        assert 'path.id = "builder-" + field.name + "-path"' in JS

    def test_the_picker_kind_is_the_field_name(self):
        # The derivation above only works because the KINDS table is keyed
        # to agree; this is the coupling, pinned.
        assert _plumed_field().name in KINDS

    def test_the_precedence_rule_is_stated_where_it_applies(self):
        # Two slots feed one `script` key. The rule -- written text wins --
        # is told to the person at the moment both are filled, not left to
        # be discovered from the config file afterwards.
        assert "clear it to use the file path" in JS

    def test_an_empty_control_writes_nothing(self):
        # An on-switch over two empty slots is not a study anybody
        # described, so the control refuses to put `enabled: true` with no
        # script into a config.
        assert "if (!text.trim() && !file) return null;" in JS


class TestTheHelpNamesBothRoutes:
    def test_the_help_speaks_to_the_person_choosing(self):
        help_text = _plumed_field().help
        assert "travels inside the config" in help_text
        assert ".dat" in help_text
        assert "openmm-plumed" in help_text


class TestOneScriptStringLeavesThePage:
    def test_a_dict_from_the_page_passes_through_untouched(self):
        # The composite sends a real object, unlike the mapping textareas
        # whose YAML text is parsed server-side. _coerce must recognise it
        # as already being what the field wants.
        assert _coerce(dict(PLUMED), _plumed_field()) == PLUMED

    def test_the_config_carries_the_script(self):
        config = build_config(_state())
        assert config["simulation"]["plumed"] == PLUMED

    def test_a_full_config_carries_it_too(self):
        config = build_config(_state(), full=True)
        assert config["simulation"]["plumed"] == PLUMED


class TestTheOtherLanguagesAnswerForIt:
    def test_the_command_line_refuses_it_by_name(self):
        # There is no flag that could carry a script, and the refusal names
        # the setting rather than dropping it -- which is what routes the
        # YAML header onto its keep-this-as-a-config-file parenthetical.
        with pytest.raises(UntranslatableSetting) as refusal:
            cli_command(build_config(_state()))
        assert "simulate.plumed" in str(refusal.value)

    def test_the_yaml_says_so_and_still_carries_the_script(self):
        result = config_yaml(_state())
        assert result["ok"], result["error"]
        assert result["command"] is None
        assert "keep this study as a config file" in result["yaml"]
        parsed = yaml.safe_load(result["yaml"])
        assert parsed["simulation"]["plumed"]["enabled"] is True
        assert "DISTANCE" in parsed["simulation"]["plumed"]["script"]

    def test_the_python_script_still_translates(self):
        result = config_yaml(_state())
        script = result["script"]
        compile(script, "<fastmdxplora_study.py>", "exec")
        assert "plumed" in script
        assert "DISTANCE" in script
