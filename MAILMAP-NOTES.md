# Contributor identities: what `.mailmap` fixes, and what it cannot

Nine author identities appear on `main` for four people. Three of them are
the same person on different machines, and one is a misconfigured
environment. `.mailmap` collapses them for every tool that reads it:

```
$ git shortlog -sne main          # before                    # after
  525  Adekunle Aina <kunleaina@gmail.com>            595  Adekunle Aina <kunleaina@gmail.com>
   53  Tomato_Cultivator <...Paradoxicaly@...>         76  Derrick K <...Paradoxicaly@...>
   25  aina <kunleaina@gmail.com>                      20  princeote <princeote12@gmail.com>
   23  Derrick K <...Paradoxicaly@...>                 10  Victory0102 <victoryunegbu2010@gmail.com>
   22  Adekunle Aina <aaina@csudh.edu>
   20  princeote <princeote12@gmail.com>
   18  Adekunle Aina <aaina@Adekunles-MBP.lan>
   10  Victory0102 <victoryunegbu2010@gmail.com>
    5  Review <review@localhost>
```

## What it does not fix

**GitHub's contributors graph.** That graph is built from the default
branch by matching each commit's author email against an email registered
on a GitHub account. GitHub's documentation for it does not mention
`.mailmap`, and the community question asking exactly this has no
authoritative answer. Do not assume the sidebar will change.

That splits the three non-canonical addresses into two kinds:

| address | commits | deliverable? | fix |
|---|---|---|---|
| `aaina@csudh.edu` | 22 | yes | add it to the GitHub account — links retroactively, no history change |
| `aaina@Adekunles-MBP.lan` | 18 | no — git's `user@hostname` fallback | history rewrite only |
| `review@localhost` | 5 authored, 9 committed | no | history rewrite only |

## The rewrite, and why it is not done yet

Rewriting is straightforward — `git filter-repo` takes this same file:

```bash
pip install git-filter-repo
git clone git@github.com:aai-research-lab/FastMDXplora.git fresh && cd fresh
git filter-repo --mailmap ../FastMDXplora/.mailmap
git remote add origin git@github.com:aai-research-lab/FastMDXplora.git
git push --force origin main
```

**The cost is the provenance ledger.** Every commit from 2026-08-19 onward
gets a new SHA, and the validation programme records which commit produced
each result:

```
v6-guardrails      998f682   2.5.6.dev31+g998f682c3
v3-tip3p-density   998f682   (finished under af6c859)
```

Those SHAs appear in `.queue-state`, in `run.out`, in every run's
`manifest.json`, in the `setuptools_scm` version string baked into each
result, and in the manuscript's provenance column. After a rewrite they
name commits that do not exist, and a reader who tries to check one gets
nothing. That is precisely the failure the ledger was built to prevent.

**So the order is: finish the programme, then rewrite, then remap, then
deposit.** `git filter-repo` writes `.git/filter-repo/commit-map` — two
columns, old SHA and new — which makes the remap mechanical and auditable
rather than a hand edit:

```bash
# after the rewrite, from the fresh clone
awk 'NR>1 && $1 != $2 {print substr($1,1,7), substr($2,1,7)}' \
    .git/filter-repo/commit-map > /tmp/sha-map
# then rewrite the short SHAs in the ledger and logs
while read -r old new; do
  sed -i "s/\\b$old/$new/g" ~/aai-research-lab/.queue-state \
                            ~/aai-research-lab/run.out
done < /tmp/sha-map
```

Keep `commit-map` itself in the Zenodo deposit. It is the evidence that the
pre-rewrite SHAs quoted anywhere else correspond to the post-rewrite
history, and without it the rewrite destroys the audit trail rather than
tidying it.

## One decision inside the mapping

`review@localhost` committed four commits that **Derrick authored** (GUI
overlay and test work, 2026-08-19). Mapping that identity leaves Derrick as
their author and makes Adekunle Aina their committer — which is what
running the commit means. Author credit is unchanged.
