# Changelog
All notable changes to this project will be documented in this file.

This format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]
### Added
- **Fraction of Native Contacts (Q-Value) Analysis**: New `qvalue` module implementing the Best-Hummer-Eaton Q metric for measuring protein folding state and native structure retention.
  - User-configurable parameters: `reference_frame`, `beta_const`, `lambda_const`, `native_cutoff`
  - Automatic native contact identification and metadata reporting
  - Integration with multi-analysis orchestrator and slide generation
  - Full CLI support with `fastmda qvalue` command
- (Add upcoming features here.)
### Changed
- (Add behavior/CLI/API changes here.)
### Fixed
- (Add bug fixes here.)

## [1.1.0] - 2026-02-15
### Added
- **Dihedral analysis** module with Ramachandran plots, including error bars and per-residue frame plots; expanded CLI support and documentation. (#16, #18, #19, #20)
- **Fraction of Native Contacts (Q-value)** analysis module; integrated into analysis discovery and workflow execution.
- **Adaptive plotting helpers** for slide-ready figures and expanded plotting regression coverage (including slide style helpers).
- **Permissive options passthrough** and test coverage for the options forwarder.
- **CODE_OF_CONDUCT.md**.

### Changed
- **MDS dimensionality reduction** updated for robust scikit-learn compatibility via API detection and exhaustive fallback parameter handling; avoids deprecated parameters across versions.
- Improved clustering and dendrogram plot styling (layout/whitespace reduction; x-axis and label cleanup).
- Enhanced analysis outputs with `compute_stat` summaries and refined plot bounds/defaults. (#21)
- Documentation updates for plotting and output formats; README improvements including conda-forge installation notes.
- Metadata updates for Zenodo DOI (pyproject.toml, CITATION, README).

### Fixed
- Fixed dihedral residue selection (`--residues`) to properly restrict computation; corrected option aliases and added coverage. (#16, #18)
- Fixed CLI `-o/--output` flag issue. (#7)
- Improved warning handling: filtered pyparsing deprecation warnings from MDTraj; fixed pytest warning filter compatibility across pyparsing versions.
- Updated plotting tests to match the current plotting API; removed outdated clustering utilities test file.
- CI/test maintenance for scikit-learn compatibility and warning suppression.

### Chores
- Version bump to 1.1.0.
- Ignored local validation outputs and example data in version control.

## [1.0.0] - 2025-11-02
### Added
- **First stable release of FastMDAnalysis**
- Version header in log files.

## [0.0.3] - 2025-10-27
### Added
- **First public release of FastMDAnalysis** (non-testing).
- Command-line workflow: `fastmda analyze --traj <traj> --top <top> --include rmsd rg cluster`
- Core analyses: RMSD, Rg, RMSF, hydrogen bonds, secondary structure, SASA, clustering, dimensionality reduction.
- Publication-quality figure generation.
- Automatic slide generation.
- Initial documentation/README.
- Packaging and metadata for citation (`CITATION.cff`) and licensing.

---

[Unreleased]: https://github.com/aai-research-lab/FastMDAnalysis/compare/v1.1.0...HEAD
[1.1.0]: https://github.com/aai-research-lab/FastMDAnalysis/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/aai-research-lab/FastMDAnalysis/compare/v0.0.3...v1.0.0
[0.0.3]: https://github.com/aai-research-lab/FastMDAnalysis/releases/tag/v0.0.3
