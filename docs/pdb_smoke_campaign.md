# PDB Smoke Campaigns

This is a staged hardening workflow for finding FastMDXplora setup,
simulation, analysis, and report failures across diverse PDB structures. It is
not intended to test every PDB entry in normal CI.

## Staged Strategy

1. Run unit tests and campaign-runner dry tests.
2. Run a fake/missing-file campaign to verify summaries and error handling.
3. Run one or two real tiny proteins, such as `1L2Y` and `1CRN`.
4. Run the starter list in `examples/pdb_list.txt` and classify failures.
5. Expand to 10-20 hand-picked PDBs per category after the starter set is clean.
6. Move to HPC batches only after the small campaign has clear failure buckets.
7. Ask before scaling to hundreds or thousands of structures.

## Curated Starter Set

| PDB | Why included | Expected difficulty | Expected behavior | Failure counts as |
| --- | --- | --- | --- | --- |
| `1L2Y` | Tiny Trp-cage protein | Low | Setup and gentle simulation should pass with OpenMM/PDBFixer | Likely code bug if it fails with an internal Python error |
| `1CRN` | Small stable crambin | Low | Basic protein smoke test | Likely code bug if standard artifacts are missing |
| `1UBQ` | Compact mixed alpha/beta protein | Low | Normal protein baseline | Likely code bug if setup/simulation crashes internally |
| `2MHR` | Alpha-helical protein | Low/medium | Helical protein baseline | Likely code bug unless input repair fails clearly |
| `1TEN` | Beta-rich fibronectin domain | Low/medium | Beta-sheet topology baseline | Likely code bug unless preparation reports a clear input issue |
| `1AKE` | Larger mixed alpha/beta enzyme | Medium | Should run, but slower than tiny proteins | Expected limitation if too large for smoke settings |
| `1A8O` | Selenomethionine/nonstandard residue | Medium | May require residue replacement or template support | Expected limitation if template error is clear |
| `1HHP` | Protein with ligand/heterogens | Medium/high | May need explicit ligand parameterization | Expected limitation if ligand chemistry is unsupported |
| `4HHB` | Multi-chain hemoglobin with heme/iron | High | Cofactor/metal support may fail | Expected limitation if heme/metal templates are unsupported |
| `1BNA` | DNA duplex | High | Nucleic acid support may be outside current workflow | Expected limitation unless nucleic-acid support is claimed |
| `1FAT` | Multi-chain lectin with heterogens/waters | Medium | Tests chain and heterogen handling | Needs review; code bug if artifacts contradict success |
| `1A3N` | Alternate locations/multi-chain hemoglobin family | Medium/high | Tests altloc behavior | Expected limitation if altloc/template issue is clear |
| `2PTC` | Protein-protein complex | Medium | Multi-chain setup and simulation smoke | Likely code bug if normal protein complex fails internally |
| `1CAG` | Large repeat-rich fragment | High | Should be handled carefully or skipped | Expected limitation if too large for smoke settings |
| `6VXX` | Large membrane/spike glycoprotein | Very high | Not appropriate for tiny smoke settings | Expected limitation/skip candidate |

## Local Commands

Use a small dry run first:

```bash
python scripts/run_pdb_smoke_campaign.py \
  --output-root runs/pdb_smoke_missing \
  missing_file.pdb \
  --continue-on-error
```

Run two real smoke cases:

```bash
python scripts/run_pdb_smoke_campaign.py \
  --output-root runs/pdb_smoke_tiny \
  --preset gentle \
  1L2Y 1CRN \
  --continue-on-error
```

Run the curated starter list:

```bash
python scripts/run_pdb_smoke_campaign.py \
  --input-list examples/pdb_list.txt \
  --output-root runs/pdb_smoke_starter \
  --preset gentle \
  --continue-on-error
```

If DNS is broken, download PDB files elsewhere and use local paths:

```bash
mkdir -p local_pdbs
curl -L -o local_pdbs/1L2Y.pdb https://files.rcsb.org/download/1L2Y.pdb
curl -L -o local_pdbs/1CRN.pdb https://files.rcsb.org/download/1CRN.pdb
printf '%s\n' local_pdbs/1L2Y.pdb local_pdbs/1CRN.pdb > local_pdb_files.txt
```

```bash
python scripts/run_pdb_smoke_campaign.py \
  --input-list local_pdb_files.txt \
  --output-root runs/pdb_smoke_local \
  --preset gentle \
  --continue-on-error
```

The campaign writes `campaign_summary.csv` / `campaign_summary.json` under
the output root, plus compatibility copies named `summary.csv` /
`summary.json`. Each protein gets its own output directory.

## Windows PowerShell Commands

If the installed `fastmdx` command is not on PATH, use
`python -m fastmdxplora.cli.main info` to verify the package and run the
campaign script with `python`.

Download two local PDB inputs:

```powershell
mkdir local_pdbs -Force
curl.exe -L -o local_pdbs\1L2Y.pdb https://files.rcsb.org/download/1L2Y.pdb
curl.exe -L -o local_pdbs\1CRN.pdb https://files.rcsb.org/download/1CRN.pdb
python -m fastmdxplora.cli.main info
```

Run a small professor/demo output campaign with the gentle preset:

```powershell
python scripts\run_pdb_smoke_campaign.py `
  --output-root local_runs\professor_outputs `
  --preset gentle `
  --continue-on-error `
  local_pdbs\1L2Y.pdb local_pdbs\1CRN.pdb
```

Verify the expected outputs:

```powershell
Get-ChildItem local_runs\professor_outputs -Recurse |
  Where-Object { $_.Name -match 'report\.md|slides\.pptx|project_bundle\.zip|production\.dcd|state_minimized\.xml|summary\.csv|summary\.json|manifest\.json|\.png|\.dat' } |
  Select-Object FullName
Get-Content local_runs\professor_outputs\summary.csv
```

Package outputs for review:

```powershell
New-Item -ItemType Directory -Force private_reports
Compress-Archive -Path local_runs\professor_outputs, README.md, docs\pdb_smoke_campaign.md, examples\pdb_list.txt -DestinationPath private_reports\FastMDXplora_professor_outputs.zip -Force
```

`local_runs/`, `local_pdbs/`, and `private_reports/` are intended for local
artifacts and are ignored by git.

## GPU and offline commands

Run the campaign on a verified GPU environment only after a short real
OpenMM smoke test. Do not assume a site module name, scheduler, checkout path,
or internet access. See [Production and GPU runs](production.md) for generic
remote-GPU and offline-host procedures.

For DNS-restricted or offline machines, stage local files and point the list
at paths:

```text
/scratch/$USER/pdbs/1L2Y.pdb
/scratch/$USER/pdbs/1CRN.pdb
/scratch/$USER/pdbs/1UBQ.pdb
```

Then run from the verified environment:

```bash
micromamba run -n <verified-env> python scripts/run_pdb_smoke_campaign.py \
  --input-list local_pdb_files.txt \
  --output-root ~/runs/pdb_smoke_campaign_local \
  --preset gentle \
  --continue-on-error
```

## Failure Classification

The campaign summarizes failures as:

- `DNS/download failure`
- `unsupported residue/template failure`
- `ligand unsupported`
- `metal unsupported`
- `missing atoms/residues issue`
- `bad geometry/clash`
- `solvation/box issue`
- `OpenMM NaN`
- `analysis failure`
- `report generation failure`
- `code exception/bug`
- `missing dependency`
- `missing input file`
- `unknown`

## What Counts As A Bug

Mark a failure as a likely code bug when:

- A normal protein with no unusual chemistry crashes with an internal Python error.
- A phase reports success but expected artifacts are missing.
- `manifest.json` says success but files are absent.
- Particle counts mismatch between setup artifacts.
- Validation catches NaN/Inf after a phase reported success.
- Analysis or report generation crashes on valid simulation outputs.

Mark it as an expected limitation or input issue when:

- An unsupported ligand, metal, or nonstandard residue causes a clear template error.
- DNS prevents downloading.
- Missing experimental data cannot be repaired.
- The protein is too large for smoke-test settings.
- Ligand chemistry is not parameterized.

## Optional Integration Test

Normal tests do not run real MD. To opt into tiny real OpenMM/PDBFixer smoke
tests:

```bash
FASTMDX_RUN_OPENMM_TESTS=1 pytest tests/test_pdb_smoke_campaign.py
```
