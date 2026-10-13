"""The turns before the last twelve are kept as what was decided and what is
still open, and the Agent's plan is kept across turns as a checklist.

The turns before the last twelve were kept by their first sentence, so
"Thanks, that helps. Use charmm36." lost the force field fourteen turns
later. Now the person's words are kept whole, a config the Agent wrote as
the settings it set or changed, a question with whether the next message
answered it, each with the turn it came from (`agent/context.py`). The plan
the Agent keeps (`keep_plan`) is shown under its reply as a checklist and
sent back with the next message, so the twentieth turn knows which step is
next; it does nothing on its own.
"""

from __future__ import annotations

import json
import types
from pathlib import Path

from fastmdxplora.agent.context import MOST_PLAN_STEPS, earlier_kept, plan_of, plan_said
from fastmdxplora.agent.turns import ToolCall, Turn, Usage

STATIC = Path(__file__).resolve().parents[1] / "src" / "fastmdxplora" / "gui" / "static"

ASIDES = [{"role": "user" if n % 2 == 0 else "agent", "text": f"Aside {n}. More of it."}
          for n in range(12)]


def _settled() -> list[dict[str, str]]:
    return [
        {"role": "user", "text": "Simulate 1UBQ."},
        {"role": "agent", "text": "Which force field?\n\nCandidates: amber14; charmm36."},
        {"role": "user", "text": "Thanks, that helps. Use charmm36 for the force field."},
        {"role": "agent", "text": "Noted.\n\nWrote a config:\nsystems:\n- system: 1UBQ\n"
                                  "setup:\n  forcefield: charmm36\nsimulation:\n"
                                  "  duration_ns: 10\nagent: assisted\n"},
        {"role": "user", "text": "Make it 310 K."},
        {"role": "agent", "text": "Wrote a config:\nsystems:\n- system: 1UBQ\nsetup:\n"
                                  "  forcefield: charmm36\nsimulation:\n  duration_ns: 10\n"
                                  "  temperature_K: 310\n"},
    ]


def test_a_choice_in_a_second_sentence_is_kept_with_its_turn():
    kept = earlier_kept(_settled() + ASIDES)
    assert kept.startswith("## Earlier in this conversation\n")
    assert "- Turn 3, the person: Thanks, that helps. Use charmm36 for the force field." in kept
    assert "- Turn 2, the Agent asked: Which force field? Candidates: amber14; charmm36." in kept
    # A config as the settings it set, then as what changed.
    assert ("- Turn 4, the Agent wrote config 1: systems 1UBQ; setup.forcefield "
            "charmm36; simulation.duration_ns 10.") in kept
    assert ("- Turn 6, the Agent wrote config 2, changing simulation.temperature_K "
            "unset to 310.") in kept
    # The twelve kept whole are not in it.
    assert "Aside" not in kept
    assert earlier_kept(ASIDES) == ""


def test_a_question_moved_past_is_still_open():
    thread = [
        {"role": "user", "text": "Simulate lysozyme."},
        {"role": "agent", "text": "Which one?\n\nCandidates: 1AKI; 2LZM."},
        {"role": "user", "text": "Before that, does the Agent keep my conversations?"},
    ] + ASIDES
    kept = earlier_kept(thread)
    assert "### Still open" in kept
    assert ("- Turn 2, the Agent asked: Which one? Candidates: 1AKI; 2LZM. The person's "
            "next message asked something else, so it may still be open.") in kept
    # Answered by naming a candidate: not open, and said as asked.
    answered = earlier_kept(thread[:2] + [{"role": "user", "text": "2LZM, the T4 one?"}]
                            + ASIDES)
    assert "### Still open" not in answered
    assert "- Turn 2, the Agent asked: Which one?" in answered
    # Answered later, after a question of the person's own: not open.
    later = earlier_kept(thread[:3] + [
        {"role": "agent", "text": "Yes, each is kept with its study."},
        {"role": "user", "text": "Then 2LZM."}] + ASIDES)
    assert "### Still open" not in later
    # Nor where a config was written after it, which settled what it asked.
    written = earlier_kept(thread[:3] + [
        {"role": "agent", "text": "Wrote a config:\nsystems:\n- system: 2LZM\n"}] + ASIDES)
    assert "the Agent asked: Which one?" in written and "may still be open" not in written


def test_a_config_not_run_is_open_and_a_run_closes_it():
    thread = _settled() + ASIDES
    assert "- Config 2 (turn 6) is written and not run." in earlier_kept(thread)
    ran = _settled() + [{"role": "agent", "text": "DO: run"}] + ASIDES
    kept = earlier_kept(ran)
    assert "is written and not run" not in kept
    assert "- Turn 7, the action named: run." in kept


def test_a_run_started_from_the_page_closes_it_and_a_run_declined_does_not():
    started = _settled() + [{"role": "agent", "text": "Started version 2."}] + ASIDES
    kept = earlier_kept(started)
    assert "is written and not run" not in kept
    assert "- Turn 7, the page: Started version 2." in kept
    declined = _settled() + [{"role": "agent", "text": "DO: run"},
                             {"role": "agent", "text": "Not run."}] + ASIDES
    assert "- Config 2 (turn 6) is written and not run." in earlier_kept(declined)
    # A fix declined after a run is not the run declined.
    fix_declined = _settled() + [
        {"role": "agent", "text": "DO: run"}, {"role": "user", "text": "it failed, why"},
        {"role": "agent", "text": "Run the fix: raise the box padding?"},
        {"role": "user", "text": "no"}, {"role": "agent", "text": "Not run."}] + ASIDES
    assert "is written and not run" not in earlier_kept(fix_declined)
    # Declined with the person's answer between: the run named is undone.
    answered_no = _settled() + [{"role": "agent", "text": "DO: run"},
                                {"role": "user", "text": "no"},
                                {"role": "agent", "text": "Not run."}] + ASIDES
    assert "- Config 2 (turn 6) is written and not run." in earlier_kept(answered_no)
    # "run the fix" is another action, not this config run.
    fixed = _settled() + [{"role": "agent", "text": "DO: run the fix"}] + ASIDES
    assert "- Config 2 (turn 6) is written and not run." in earlier_kept(fixed)


def test_a_config_with_aliases_or_too_long_is_not_read_back():
    import time

    bomb = "Wrote a config:\n" + "".join(
        f"l{i}: &l{i} [" + ", ".join([f"*l{i - 1}" if i else "x"] * 8) + "]\n"
        for i in range(10))
    started = time.monotonic()
    kept = earlier_kept([{"role": "agent", "text": bomb}] + ASIDES)
    assert time.monotonic() - started < 2
    assert "wrote config" not in kept
    deep = "Wrote a config:\na: " + "[" * 5000
    assert "wrote config" not in earlier_kept([{"role": "agent", "text": deep}] + ASIDES)
    long = "Wrote a config:\n" + "".join(f"k{i}: {i}\n" for i in range(5000))
    assert "wrote config" not in earlier_kept([{"role": "agent", "text": long}] + ASIDES)


def test_a_long_message_is_cut_and_says_so():
    long = "word " * 400
    kept = earlier_kept([{"role": "user", "text": long}] + ASIDES)
    line = next(line for line in kept.splitlines() if line.startswith("- Turn 1"))
    assert line.endswith("... (cut)") and len(line) < 650


def test_both_ways_of_asking_carry_the_kept_turns_and_the_old_way_stays():
    from fastmdxplora.agent.conversation import _history
    from fastmdxplora.agent.propose import earlier_turns, prompt_for

    thread = _settled() + ASIDES
    prompt = prompt_for("write it", history=thread)
    assert "Use charmm36 for the force field." in prompt
    assert "Use charmm36" in _history(thread)[0]["text"]
    # Without the context: each turn by its first sentence, as before.
    old = earlier_turns(thread, context=False)
    assert "- Person: Thanks, that helps." in old and "charmm36 for" not in old
    assert "Use charmm36" not in _history(thread, context=False)[0]["text"]


def test_a_plan_is_kept_bounded_and_said_as_a_checklist():
    steps = plan_of([{"step": " Write  the study ", "done": True}, "Run it", {"step": ""},
                     {"step": "x" * 500}, 7])
    assert steps[0] == {"step": "Write the study", "done": True}
    assert steps[1] == {"step": "Run it", "done": False}
    assert len(steps) == 3 and len(steps[2]["step"]) == 200
    assert len(plan_of(["s"] * 40)) == MOST_PLAN_STEPS
    assert plan_of({"steps": ["one"]}) == ({"step": "one", "done": False},)
    assert plan_of("not a plan") is None and plan_of([]) is None
    said = plan_said(steps)
    assert said.startswith("## Your plan so far\n")
    assert "1. [x] Write the study" in said and "2. [ ] Run it" in said
    assert plan_said(None) == ""


def _turns(*turns):
    left = list(turns)
    seen: list = []

    def complete(prompt: str, **_: object) -> str:
        raise AssertionError("asked in text")

    def turn(system, messages, tools, **_: object):
        seen.append((system, list(messages), [t.name for t in tools]))
        return left.pop(0)

    complete.turn = turn
    complete.seen = seen
    return complete


def test_a_plan_kept_beside_an_answer_comes_back_with_it():
    from fastmdxplora.agent.propose import propose_config

    plan = ToolCall("p1", "keep_plan", {"steps": [{"step": "Write the study", "done": True},
                                                  {"step": "Run it"}]})
    complete = _turns(Turn(text="I wrote it. Say run when you have read it.",
                           calls=(plan,), usage=Usage()))
    proposal = propose_config("set it up, run it, then analyse it", complete,
                              plan=[{"step": "Write the study"}])
    assert proposal.answer.startswith("I wrote it.")
    assert proposal.steps == ({"step": "Write the study", "done": True},
                              {"step": "Run it", "done": False})
    system, messages, tools = complete.seen[0]
    assert "keep_plan" in tools and "`keep_plan`" in system
    # The plan given back, in the message: what the AI model keeps to.
    assert "## Your plan so far" in messages[-1]["text"]
    assert "1. [ ] Write the study" in messages[-1]["text"]


def test_a_plan_kept_alone_is_answered_and_the_reply_comes_next():
    from fastmdxplora.agent.propose import propose_config

    plan = ToolCall("p1", "keep_plan", {"steps": ["Find the structure", "Write the study"]})
    complete = _turns(Turn(text="", calls=(plan,), usage=Usage()),
                      Turn(text="Which lysozyme do you mean?", calls=(), usage=Usage()))
    proposal = propose_config("simulate lysozyme, run it and analyse it", complete)
    assert proposal.answer == "Which lysozyme do you mean?"
    assert [s["step"] for s in proposal.steps] == ["Find the structure", "Write the study"]
    results = complete.seen[1][1][-1]
    assert results["role"] == "results"
    assert results["results"][0]["content"] == "Kept: 2 steps."
    assert not proposal.attempts


def test_a_plan_not_given_as_steps_is_not_kept_and_says_why():
    from fastmdxplora.agent.propose import propose_config

    plan = ToolCall("p1", "keep_plan", {"steps": "first this, then that"})
    complete = _turns(Turn(text="", calls=(plan,), usage=Usage()),
                      Turn(text="Done.", calls=(), usage=Usage()))
    proposal = propose_config("plan it", complete)
    assert proposal.steps is None
    said = complete.seen[1][1][-1]["results"][0]
    assert said["is_error"] and "Not kept" in said["content"]


def test_the_text_protocol_is_not_told_of_a_plan():
    from fastmdxplora.agent.propose import _instructions

    assert "keep_plan" not in _instructions("text")
    assert "keep_plan" in _instructions("tools")


def test_the_page_sends_the_plan_and_gets_the_new_one(monkeypatch, tmp_path):
    import fastmdxplora.agent as agent_module
    from fastmdxplora.agent.propose import Proposal
    from fastmdxplora.gui import agent_panel

    given = {}

    def propose(request, complete, **kw):
        given.update(kw)
        return Proposal(config=None, attempts=(), answer="Next, run it.",
                        steps=({"step": "Write it", "done": True},
                               {"step": "Run it", "done": False}))

    monkeypatch.setattr(agent_module, "completion_for", lambda *a, **k: object())
    monkeypatch.setattr(agent_module, "propose_config", propose)
    monkeypatch.setattr(agent_panel, "GROW_IN_BACKGROUND", False)
    runtime = types.SimpleNamespace(exploration_root=str(tmp_path), active_root=None)
    answer = agent_panel.propose_endpoint(
        {"request": "go on", "steps": [{"step": "Write it", "done": True}, "Run it", 3]},
        runtime, memory_store=None)
    assert given["plan"] == ({"step": "Write it", "done": True},
                             {"step": "Run it", "done": False})
    assert answer["steps"] == [{"step": "Write it", "done": True},
                               {"step": "Run it", "done": False}]
    # A config's own `plan` (the config said in words) is another thing.
    assert "plan" not in answer
    json.dumps(answer)


def test_the_panel_shows_the_plan_as_a_checklist_and_sends_it_back():
    js = (STATIC / "agent-panel.js").read_text(encoding="utf-8")
    assert "function stepsCard(box, steps)" in js
    assert "steps: currentSteps" in js
    # Kept with the reply's entry, restored from it, and cleared with the thread.
    assert "if (steps) entry.steps = steps;" in js
    assert "var kept = stepsOf(e.steps);" in js
    assert js.count("currentSteps = null;") >= 2
    # Ticks are shown, never pressed: the plan does nothing.
    assert "tick.disabled = true;" in js
