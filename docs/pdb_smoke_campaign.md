# PDB smoke campaigns

*For maintainers.*

A smoke campaign runs FastMDXplora over many PDB structures at once, short, to
find what breaks. Real structures are more varied than any test suite: unusual
residues, missing density, exotic cofactors, chain breaks in awkward places.
Running fifty of them for ten picoseconds each finds failures that fifty
unit tests will not.

This is for hardening a release, not for ordinary CI.

```bash
python scripts/run_pdb_smoke_campaign.py \
  --output-root runs/campaign \
  1UBQ 1L2Y 4LYT 181L 1AFO
```

---

## Choosing what to run

Positional arguments are PDB identifiers or paths to local `.pdb` / `.cif`
files. For anything more than a handful, keep them in a file:

```bash
python scripts/run_pdb_smoke_campaign.py \
  --output-root runs/campaign \
  --input-list examples/pdb_list.txt
```

One entry per line; `#` starts a comment.

A good list mixes the easy and the awkward: a small soluble protein, one with
a ligand, one with a metal, one with a nucleic acid, one with missing loops,
one that is large. The point is coverage of what breaks, not of what works.

---

## Keeping it short and bounded

The defaults are already short. The flags that matter are the ones that stop a
single structure eating the campaign:

| | |
|---|---|
| `--preset gentle` | conservative simulation settings |
| `--nvt-steps`, `--npt-steps`, `--production-steps` | shorter still |
| `--max-input-mb` | skip a structure whose file is larger |
| `--max-setup-atoms` | skip one that solvates to more atoms than this |
| `--platform` | force CPU or CUDA |
| `--no-report` | skip the reporting phase |

The two limits are what keep a campaign finishing overnight. A ribosome will
otherwise consume the time budget for everything after it.

---

## Letting it finish

```bash
python scripts/run_pdb_smoke_campaign.py \
  --output-root runs/campaign \
  --input-list examples/pdb_list.txt \
  --continue-on-error
```

`--continue-on-error` is the point of a campaign: one structure failing is a
result, not a reason to stop. `--stop-on-error` does the opposite when you are
chasing one specific failure.

---

## Reading the results

Two summaries are written to the output root:

- **`campaign_summary.csv`** — one row per structure, for sorting and
  filtering
- **`campaign_summary.json`** — the same with the detail

Each structure gets one of these:

| | |
|---|---|
| `ok` | ran through |
| `expected_limitation` | refused something it should refuse — a structure the software correctly declines |
| `validation_failed` | produced output that did not pass its own checks |
| `failed` | a phase raised |
| `skipped` | over a size limit |
| `error` | the campaign itself had a problem |

Each run's own directory is beside the summaries, so a failure can be opened
and read like any other run.

---

## What counts as a bug

**`expected_limitation` is not a bug.** A structure with an unparameterisable
cofactor, a charge that cannot be settled, or geometry the force field has no
templates for *should* be refused. A campaign with several of these is
working.

**`validation_failed` usually is.** It means a phase produced output and the
output was wrong, which is worse than refusing.

**`failed` needs reading.** A phase raising with a clear message naming what it
could not do is close to `expected_limitation`. A phase raising with a
traceback from inside a library is a bug, and the message is where to start.

The distinction is the whole reason the statuses are separate: software that
refuses cleanly and software that breaks look the same in a pass/fail count.
