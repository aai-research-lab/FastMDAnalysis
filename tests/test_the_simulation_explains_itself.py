"""The phase where the choices are made says why they were made.

`explain.py` opens by saying a pipeline that does all this silently "teaches
nothing". That held for setup, which explains its protonation and its
solvation, and did not hold here: the entries for minimisation, NVT, NPT,
production and the ensemble existed and were never printed, because the
runner announces through a callback that carried a message and nothing else.
A real run showed two explanations during setup and then went quiet for the
twenty minutes where the ensemble is actually chosen.
"""

from __future__ import annotations

import io

import pytest

from fastmdxplora.explain import EXPLANATIONS
from fastmdxplora.simulation.runner import _explanation_for
from fastmdxplora.utils.presenter import SessionPresenter


def _presenter(explain: bool = True) -> SessionPresenter:
    return SessionPresenter(stream=io.StringIO(), explain=explain)


class TestAStageKnowsItsExplanation:
    @pytest.mark.parametrize("label,key", [
        ("NVT equilibration", "nvt"),
        ("NPT equilibration", "npt"),
        ("Production", "production"),
    ])
    def test_the_label_finds_it(self, label: str, key: str) -> None:
        assert _explanation_for(label) == key

    def test_a_stage_with_none_gets_none(self) -> None:
        """Not every message is a stage, and inventing an explanation for
        one would put it beside the wrong thing."""
        assert _explanation_for("Loading System, State, and topology") is None

    def test_every_key_it_names_exists(self) -> None:
        """A typo here would be a stage that silently stopped explaining
        itself, which is the failure this whole test file is about."""
        for label in ("NVT equilibration", "NPT equilibration", "Production"):
            assert _explanation_for(label) in EXPLANATIONS


class TestThePresenterCanPrintOneOnItsOwn:
    """Stages announce through `info`, which prints no status icon -- "NVT
    equilibration: 250,000 steps" is not a thing that succeeded or failed --
    so the explanation had to be reachable from there and not only from
    `step`."""

    def test_info_carries_it(self) -> None:
        pres = _presenter()
        pres.info("NVT equilibration: 50,000 steps", explain="nvt")
        assert "held fixed" in pres.stream.getvalue()

    def test_it_can_be_printed_alone(self) -> None:
        pres = _presenter()
        pres.explanation("ensemble_choice")
        assert "NVT production is also legitimate" in pres.stream.getvalue()

    def test_no_key_prints_nothing(self) -> None:
        pres = _presenter()
        pres.info("Loading System, State, and topology")
        assert "\n\n" not in pres.stream.getvalue()

    def test_an_unknown_key_is_not_an_error(self) -> None:
        pres = _presenter()
        pres.explanation("no_such_entry")
        assert pres.stream.getvalue() == ""

    def test_explanations_off_prints_the_message_only(self) -> None:
        pres = _presenter(explain=False)
        pres.info("NVT equilibration: 50,000 steps", explain="nvt")
        out = pres.stream.getvalue()
        assert "NVT equilibration" in out and "held fixed" not in out


class TestTheRunnerReachesThePresenter:
    """The two ends were both present and not joined: the runner had no way
    to name an explanation, and `info` had no way to print one."""

    def test_the_stage_runners_call_it(self) -> None:
        from pathlib import Path
        import fastmdxplora.simulation.runner as runner

        source = Path(runner.__file__).read_text(encoding="utf-8")
        assert source.count("on_explain(_explanation_for(label))") == 2, (
            "both stage runners announce, and both should explain")
        assert 'on_explain("minimize")' in source
        assert 'on_explain("ensemble_choice")' in source

    def test_the_pipeline_passes_a_callback(self) -> None:
        from pathlib import Path
        import fastmdxplora.simulation.pipeline as pipeline

        source = Path(pipeline.__file__).read_text(encoding="utf-8")
        assert "on_explain=_explain" in source
        assert "presenter.explanation(key)" in source

    def test_it_is_a_separate_callback_from_progress(self) -> None:
        """`on_progress` is defined as one argument everywhere it is passed,
        including in tests. Widening it would mean catching TypeError, which
        also swallows a TypeError raised inside the callback itself."""
        import inspect
        from fastmdxplora.simulation.runner import _run_md_stage

        params = inspect.signature(_run_md_stage).parameters
        assert "on_explain" in params
        assert params["on_explain"].default is None


class TestTheWholeChainRenders:
    def test_a_stage_and_its_reason_arrive_together(self) -> None:
        pres = _presenter()
        label = "NPT equilibration"
        pres.info(f"{label}: 500,000 steps")
        pres.explanation(_explanation_for(label))
        out = pres.stream.getvalue()
        assert out.index("500,000 steps") < out.index("P replaces V")

    def test_every_simulation_stage_has_something_to_say(self) -> None:
        for key in ("minimize", "nvt", "npt", "production",
                    "ensemble_choice"):
            assert key in EXPLANATIONS, f"{key} is announced and undefined"
