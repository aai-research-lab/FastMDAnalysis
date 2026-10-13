# The Agent's evaluation, version 2: pre-registered

**Written 2026-10-13, before any reply from an AI model to the cases it
adds was seen.** The harness is `fastmdxplora.validation.agent_eval`, set
`agent-eval`, version 2. Version 2 is version 1 (`agent-eval.md`) with one
kind added, `context`: every case of version 1, its messages, what it must
come to, how it is judged and what is claimed of it, stays as registered
there, and a recording of version 1 replays under version 2. This file
fixes what the new kind asks, how it is judged and what is claimed of it. A
change to any of them after the first reply is seen is a new version.
Results are added under **Results** at the end, as found, and nothing above
that heading changes with them.

## Why

The Agent now knows more of where the person is (the round `agent-context`):
the page they are on and what is open in it, what they name with @ in a
message, and the turns before the last twelve, kept as what was decided
and what is still open rather than by each turn's first sentence. Each of
these is meant to make a reply right that was wrong without it. Version 1
cannot show it: none of its cases has a page, an @, or more than twelve
turns. The `context` kind asks for exactly what each part adds, and is
asked twice, with the Agent's context and without it, so what it adds is
counted against the same AI model on the same day.

## How a case is asked

As in version 1, with these additions for a `context` case:

- **Earlier turns.** The turns the case lists come before its message, in
  the conversation as the page keeps it (a question with its candidates on
  a line `Candidates: ...`). Each far case has sixteen: four that settle or
  ask something, then twelve about the software (`_ASIDES`), so what was
  settled is no longer among the last twelve kept whole.
- **The page.** Where the case lists the page's hints, they are given to
  the look at what the page shows (`current_view`), `{study}` standing for
  the finished study given, and the open study is that study.
- **The folder of studies.** With a study given, the folder it is in is
  the folder of studies, as the GUI started there would have it.

**With the Agent's context** (the default), each message is asked as the
page asks it on the branch that adds it: the earlier turns kept as decided
and open, what is named with @ resolved to what it names, and the page's
hints whole.

**Without it** (`--without-context`), the Agent as it is on that branch
with these three parts left out: the earlier turns each by its first
sentence, as before; nothing named with @ resolved (the @ is in the message
as typed); and the page's hints only its name and the study open, which is
what the page sent before from every page while nothing had been shown in
the Viewer (once something had, its frame and view went too). The rest of
the branch (the plan, the look at a study's file) is in both arms.

Both are asked by the same AI model, the one the Agent is set to, three
times each, on the same day, the one after the other.

A recording replays as version 1's do. What is named with @ is resolved
from the folders on the machine replaying it, so where the study is not
there the requests of the two @ cases are counted as asked differently;
their verdicts replay all the same.

## The cases

Eight, in three groups. The four that read a finished study are asked only
where one is given (`--study`), as version 1's `records` cases are.

### What was settled many turns before (judged as `write`)

- `ctx_far_forcefield`: earlier "I want to simulate ubiquitin, 1UBQ, in
  water.", the Agent asking "Which force field should it use?" with the
  candidates amber14 and charmm36, then "Thanks, that helps. Use charmm36
  for the force field.", and "Understood. I will use it when I write the
  study."; the message "Now write the study for it, 10 ns." Must hold
  `setup.forcefield` charmm36 and `simulation.duration_ns` 10, of 1UBQ.
- `ctx_far_temperature`: earlier "I will be running trp-cage studies this
  week.", "Good. Tell me what each should hold and I will write it.",
  "Fine. Every run should be at 310 K, not the default.", "Understood. I
  will keep to that."; the message "Write the study for 1L2Y, 5 ns." Must
  hold `simulation.temperature_K` 310 and `simulation.duration_ns` 5, of
  1L2Y.
- `ctx_far_ph`: earlier "I am reproducing a paper's ubiquitin runs.",
  "Which conditions does the paper use?", "It is a standard setup. The
  protonation is set at pH 6.0, as the paper did.", "Understood. I will set
  that when I write it."; the message "Write it for 1UBQ, 10 ns." Must hold
  `setup.ph` 6.0 and `simulation.duration_ns` 10, of 1UBQ.

In each, what decides the study is in the second sentence of the person's
turn, and the Agent's turns after it do not repeat it, so the first
sentence of each turn does not hold it.

### A question still open (judged as `ask`)

- `ctx_open_question`: earlier "Simulate lysozyme in water.", the Agent
  saying "There are two common lysozyme structures in the PDB. Hen
  egg-white (1AKI) or T4 (2LZM): which one should it be?" with the
  candidates 1AKI and 2LZM, then
  "Before that, does the Agent keep my conversations?" and "Yes. Each is
  kept with the study it is about, or in your settings folder where it is
  about none."; the message "Go ahead and write it, 10 ns." Passes as a
  question whose words or choices name 1AKI, 2LZM, hen or T4. A config of
  either is the AI model choosing for the person, and fails. Without the
  context, the Agent's turn is kept as its first sentence, which says there
  are two structures and names neither, nor that a question is open.

### What the page shows, and what is named with @ (judged as `records`)

- `ctx_file_open`: on the Files page, the study open, its
  `analysis/rg/options.json` open in the file reader; the message "What
  mean does the file I have open record?" The mean of `rg`.
- `ctx_report_section`: on the Report page, the study open, its section
  "RMSD" being read; the message "What value does this section give?" The
  mean of `rmsd`.
- `ctx_mention_study`: on the Agent page, no study open; the message "What
  mean radius of gyration did @{study_name} record?", `{study_name}` the
  name of the given study's folder. The mean of `rg`.
- `ctx_mention_analysis`: on the Overview, the study open; the message
  "What did @rmsd come to?" The mean of `rmsd`.

Each asks for the `ANSWER:` line as version 1's `records` cases do, and is
judged as they are.

## How a case is judged

On what the last message came to, as version 1 judges the kind named:

- **write:** as version 1's `write`: an accepted config whose runs hold
  each setting listed and name exactly the structures listed.
- **records:** as version 1's `records`: the reply's `ANSWER:` line, its
  first number within one standard error of the mean the study records (or
  within half a unit of the reply's last place, as there).
- **ask:** a question whose words or choices name one of those listed. An
  accepted config fails, as does an answer, an action or nothing accepted.

The rest is as version 1: a trial the provider fails is asked again twice
at most and is not judged if it never completes; any other error fails it.

## What is claimed

Counted per trial, each arm 24 trials (the eight cases three times), the
twelve that read a study only where one is given. Where trials a claim
names were not judged, the claim is not decided unless those judged decide
it.

1. **With the Agent's context, at least 18 of the 24 `context` trials
   pass.**
2. **The context adds at least 6:** with it, at least 6 more of the 24
   trials pass than without it, asked by the same AI model on the same
   day.
3. **No choice made for the person:** with the Agent's context, none of the
   3 `ctx_open_question` trials comes to a config.

The thresholds were set by judgment, before any reply was seen: three
trials in four for the kind, and a gain of one case in four for what the
context adds. The two cases named with @ can pass without the context too
(a study's folder can be listed, and the page always named the study open),
so claim 2 rests mostly on the other six.

**Checked, not claimed.** Each case's counts with and without, side by
side; which looks were taken; whether `current_view` was the look that read
the file or the section.

**What is not measured here.** The plan kept as a checklist (it does
nothing, and is held by the suite's tests); what the Analysis page, the
Config Builder and the Viewer's selection add, whose hints come in a later
round (`agent-context-2`) and are for a later version; whether a resolved
@ helps where the AI model would have found the record anyway.

## Running it

On a machine with an AI model chosen, from the branch that adds the
context, with a finished study whose RMSD and Rg means are determined:

```bash
python -m fastmdxplora.validation.agent_eval --kind context --repeats 3 --study <the study> --out agent_eval_context.json
python -m fastmdxplora.validation.agent_eval --kind context --repeats 3 --study <the study> --without-context --out agent_eval_context_without.json
```

Bring both files back: they are reported under **Results** and kept as
recordings the suite replays.

## Results

None yet.
