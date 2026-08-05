"""Why each step happens, said while it happens.

Molecular dynamics has a lot of steps that are obvious once you know them and
opaque before that. A pipeline that does all of it silently is faster to use
and teaches nothing: somebody running their first simulation ends up with a
trajectory they cannot defend.
"""

from __future__ import annotations

import io

import pytest


class TestTheExplanationsThemselves:
    def test_every_one_says_why_rather_than_what(self) -> None:
        """The step already says what happened. An explanation that repeats it
        has added nothing."""
        from fastmdxplora.explain import EXPLANATIONS

        assert EXPLANATIONS
        for key, entry in EXPLANATIONS.items():
            assert len(entry.why) > 120, (
                f"{key}: too short to explain anything"
            )
            assert len(entry.why) < 700, (
                f"{key}: this is documentation, not a step note"
            )

    def test_a_reference_is_a_real_citation_or_absent(self) -> None:
        """Inventing one to look thorough would be worse than having none."""
        from fastmdxplora.explain import EXPLANATIONS

        for key, entry in EXPLANATIONS.items():
            if entry.reference is None:
                continue
            assert any(char.isdigit() for char in entry.reference), (
                f"{key}: a reference should carry a year"
            )
            assert "," in entry.reference, (
                f"{key}: a reference should name authors and a venue"
            )

    def test_an_unknown_key_returns_nothing_rather_than_raising(self) -> None:
        """A missing explanation should not stop a run. The guard below fails
        on one, so it does not go unnoticed either."""
        from fastmdxplora.explain import explain

        assert explain("no_such_step") is None

    def test_it_wraps_for_a_terminal(self) -> None:
        from fastmdxplora.explain import EXPLANATIONS

        entry = next(iter(EXPLANATIONS.values()))
        for line in entry.as_text().splitlines():
            assert len(line) <= 100, "an explanation should not run off screen"


class TestKeysRatherThanMatching:
    """A call site names the explanation it wants, so one cannot drift onto
    the wrong step and a step cannot quietly lose its explanation."""

    def test_every_key_a_call_site_asks_for_exists(self) -> None:
        import re
        from pathlib import Path

        from fastmdxplora.explain import EXPLANATIONS

        package = Path(__file__).resolve().parents[1] / "src" / "fastmdxplora"
        asked = set()
        for path in package.rglob("*.py"):
            for match in re.finditer(r'explain=["\']([a-z_]+)["\']',
                                     path.read_text(encoding="utf-8")):
                asked.add(match.group(1))

        missing = sorted(asked - set(EXPLANATIONS))
        assert not missing, f"these steps ask for explanations that do not exist: {missing}"

    def test_the_presenter_takes_a_key_not_a_string(self) -> None:
        """Matching on the step's message would have made a mismatch
        silent."""
        import inspect

        from fastmdxplora.utils.presenter import SessionPresenter

        signature = inspect.signature(SessionPresenter.step)
        assert "explain" in signature.parameters


class TestItPrintsAndCanBeSilenced:
    @staticmethod
    def _run(explain_enabled):
        from fastmdxplora.utils.presenter import SessionPresenter

        out = io.StringIO()
        presenter = SessionPresenter(stream=out, explain=explain_enabled)
        presenter.step("Fixed PDB with PDBFixer (pH=7.4)", explain="protonation")
        return out.getvalue()

    def test_on_by_default(self) -> None:
        import inspect

        from fastmdxplora.utils.presenter import SessionPresenter

        assert inspect.signature(
            SessionPresenter.__init__).parameters["explain"].default is True

    def test_it_explains_when_on(self) -> None:
        printed = self._run(True)
        assert "X-rays do not see hydrogens" in printed
        assert "PROPKA3" in printed

    def test_it_is_silent_when_off(self) -> None:
        printed = self._run(False)
        assert "Fixed PDB" in printed, "the step itself still appears"
        assert "X-rays" not in printed

    def test_a_step_with_no_key_explains_nothing(self) -> None:
        from fastmdxplora.utils.presenter import SessionPresenter

        out = io.StringIO()
        SessionPresenter(stream=out, explain=True).step("Wrote report.md")
        assert out.getvalue().strip() == "✓ Wrote report.md".replace("✓ ", "✓ ")


class TestItReachesTheInterfaces:
    def test_the_setting_is_declared(self) -> None:
        from fastmdxplora.config.schema import TOP_LEVEL_KEYS

        assert "explain" in TOP_LEVEL_KEYS

    def test_it_defaults_to_on(self) -> None:
        from fastmdxplora.config.schema import TOP_LEVEL

        field = TOP_LEVEL.get("explain")
        assert field.default is True

    def test_the_command_line_can_turn_it_off(self) -> None:
        from fastmdxplora.cli.main import _build_parser

        args = _build_parser().parse_args(
            ["explore", "--system", "1UBQ", "--no-explain"])
        assert args.explain is False

    def test_and_leaves_it_on_otherwise(self) -> None:
        from fastmdxplora.cli.main import _build_parser

        args = _build_parser().parse_args(["explore", "--system", "1UBQ"])
        assert args.explain is True


class TestWhatTheStepsThemselvesSay:
    """A step's message is read far more often than any explanation under it,
    and two of them were wrong."""

    def test_the_solvation_message_has_no_nested_parentheses(self) -> None:
        """It read "Solvating (box=cube) (padding=1.00 nm, ...)" -- a
        parenthesis inside a parenthesis, introduced when the message learned
        to choose between solvating and embedding."""
        import inspect

        from fastmdxplora.setup import prepare

        source = inspect.getsource(prepare.prepare_system)
        assert "Solvating (box=" not in source
        assert "padding=%.2f nm, ions" in source

    def test_the_force_field_step_says_what_it_resolved_to(self) -> None:
        """It said "(auto)", which tells a reader nothing -- while the banner
        two inches above said amber-openff (auto)."""
        import inspect

        from fastmdxplora.setup import pipeline

        source = inspect.getsource(pipeline)
        assert "AUTO_FORCEFIELD" in source, (
            "the step should name the force field auto resolved to"
        )

    def test_no_step_message_nests_parentheses(self) -> None:
        """This mistake has been made twice: once when the solvation message
        learned to choose between solvating and embedding, and again in the
        fix for it, when a label already carrying "(auto)" was passed to a
        step that wraps its label in parentheses.

        Checking the shape rather than the instance, because the instance is
        not what recurs.
        """
        import re
        from pathlib import Path

        package = Path(__file__).resolve().parents[1] / "src" / "fastmdxplora"
        offenders = []
        for path in package.rglob("*.py"):
            for line in path.read_text(encoding="utf-8").splitlines():
                if ".step(" not in line and "presenter.info(" not in line:
                    continue
                # A message that both wraps a value in parentheses and builds
                # that value with parentheses of its own.
                if re.search(r'\(\{[a-z_]+\}\)', line) and "(" in line:
                    text = re.search(r'f?"([^"]+)"', line)
                    if text and text.group(1).count("(") > 1:
                        offenders.append(f"{path.name}: {text.group(1)[:60]}")
        assert not offenders, f"these step messages nest parentheses: {offenders}"
