# Scripts

Helper scripts for development, benchmarking, and release engineering.

- `make_mark.py` makes FastMDXplora's mark, a three-atom molecule going
  right with each atom's trail behind it, from its atoms, bonds and trails.
  It writes the tab's icon (`src/fastmdxplora/gui/static/fastmdx-mark.svg`),
  prints its shapes for `ICONS["mark"]` in
  `src/fastmdxplora/gui/sidebar_icons.py`, checks that both match it with
  `--check`, and writes a page of the mark at each size in both schemes
  with `--preview FILE`.
- `make_demo.py` packages a finished study made from
  `src/fastmdxplora/demo/3ptb.yml` (run on a GPU) as the demo study in
  `src/fastmdxplora/demo/3ptb/`: what the GUI shows, without what only a
  run needs, and `demo.json` recording the release and the folder it was
  made in. `fastmdx gui --demo DIR` copies it into DIR and opens it.
- `run_validation_campaign.sh` runs, on one machine and one stage after
  another, what a release is held to: `fastmdx info`, the fetch's tests
  and `check_rcsb_fetch.py`, the tests that need a real OpenMM, the
  network-marked tests and corpus, a CUDA smoke campaign over
  `examples/pdb_list.txt`, and the suite. Each stage waits for a free GPU
  and writes its log, its JUnit XML and a row of `registry.csv` to
  `--output`, which must not exist yet.
- `check_rcsb_fetch.py` asks the real RCSB, with a dropped connection, a
  503, a drop on every try and no route out put in front of it, and
  fetches an mmCIF-only entry and 6B73 cleanly, so the fetch's retry is
  seen working from the machine a study runs on.
