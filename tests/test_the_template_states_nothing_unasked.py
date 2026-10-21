"""The template `fastmdx config` writes states nothing nobody asked for.

Its top-level provenance, budget and decisions were written live with their
examples. A study made from it by giving only a system was recorded as
drafted by an AI model (`agent: assisted`, a named model), ran in stages
under a 40 hour ceiling, and its report said "`setup.ph` = `7.4`: the assay
buffer is pH 7.0 (decided by person)". Its help lines were written after the
same `# ` as its settings, so one taken off to use a block made "Integrator:
langevin_middle, ..." a key of its own, and two settings that cannot stand
together (`include` and `exclude`, `pressure_bar` and `pressure_atm`) were
offered side by side with nothing to say so.
"""

from __future__ import annotations

import yaml

from fastmdxplora.cli.main import main as cli_main
from fastmdxplora.config import generate_template, validate_config
from fastmdxplora.config.schema import PHASE_SCHEMAS


def _given_a_system(text: str) -> dict:
    data = yaml.safe_load(text)
    data["systems"] = [{"id": "trp", "system": "1L2Y"}]
    validate_config(data, require_systems=True)
    return data


def test_a_template_given_only_a_system_states_no_provenance_budget_or_decisions():
    data = _given_a_system(generate_template())
    for key in ("agent", "agent_model", "budget_hours", "decisions",
                "include_phase", "exclude_phase"):
        assert key not in data, key
    # Shown, for a person to uncomment.
    text = generate_template()
    assert "# agent_model: " in text and "# budget_hours: 40" in text
    assert "# decisions: " in text


def test_the_written_file_says_the_same(tmp_path):
    out = tmp_path / "study.yml"
    assert cli_main(["config", "-f", str(out)]) == 0
    data = yaml.safe_load(out.read_text(encoding="utf-8"))
    assert not {"agent", "agent_model", "budget_hours", "decisions"} & set(data)


def _uncommented_once(text: str) -> str:
    """The template with its phase headers and every line inside a phase
    block uncommented once, as a person using the blocks would: one `# `
    taken off each."""
    headers = {f"# {phase}:" for phase in PHASE_SCHEMAS}
    lines = []
    for line in text.splitlines():
        if line in headers:
            lines.append(line[2:])
        elif line.startswith("  # "):
            lines.append("  " + line[4:])
        else:
            lines.append(line)
    return "\n".join(lines)


def test_help_stays_a_comment_when_a_block_is_uncommented():
    text = generate_template()
    assert "  ## Integrator: langevin_middle" in text
    opened = yaml.safe_load(_uncommented_once(text))
    assert "Integrator" not in opened["simulation"]
    # Every key the uncommented blocks hold is a setting of its block.
    for phase, schema in PHASE_SCHEMAS.items():
        names = {f.name for f in schema.fields}
        assert set(opened[phase]) <= names, (phase, set(opened[phase]) - names)


def test_each_setting_of_a_pair_names_the_one_it_replaces():
    text = generate_template()
    for line in ("# exclude_phase: [report]",
                 "  # exclude: [dimred]",
                 "  # pressure_atm: 1.0",
                 "  # production_steps: 1000000"):
        found = [row for row in text.splitlines() if row.startswith(line)]
        assert found, line
        assert "instead of " in found[0] and "not beside it" in found[0], found[0]


def test_an_example_value_says_it_is_one():
    text = generate_template()
    row = next(r for r in text.splitlines() if r.startswith("  # nonbonded_cutoff_nm:"))
    assert "an example; unset by default" in row


def test_the_minimal_template_offers_no_sequence():
    assert "or sequence" not in generate_template(minimal=True)
