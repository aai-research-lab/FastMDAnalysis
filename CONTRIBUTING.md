# Contributing to FastMDXplora

Thank you for your interest in contributing to FastMDXplora. We welcome
contributions of all kinds — bug reports, feature requests, documentation
improvements, and code.

## Getting started

```bash
# Clone the repository
git clone https://github.com/aai-research-lab/FastMDXplora.git
cd FastMDXplora

# Create and activate a virtual environment
python -m venv .venv
source .venv/bin/activate   # On Windows: .venv\Scripts\activate

# Install in editable mode with development dependencies
pip install -e ".[dev]"

# Verify the install
fastmdx --version
fastmdx info
pytest
```

## Development workflow

1. Open an issue describing the change you'd like to make (skip this for
   trivial fixes such as typos).
2. Fork the repository and create a topic branch from `main`.
3. Make your changes. Write or update tests under `tests/`.
4. Run the test suite locally (`pytest`).
5. Run the linter (`ruff check src tests`) and ensure it passes.
6. Open a pull request against `main`.

## Coding conventions

- **Python ≥ 3.9.** Use modern type hints and the standard library where possible.
- **`src/` layout.** All package code lives under `src/fastmdxplora/`.
- **Docstrings.** Public functions and classes get NumPy-style docstrings.
- **Tests required for new functionality.** Smoke tests at minimum; full
  numerical/equivalence tests for any analytical code migrated from
  FastMDXplora version 1.
- **Lazy imports for heavy optional dependencies** (OpenMM, PDBFixer,
  python-pptx, etc.).
- **Consistent output structure.** Every phase writes its outputs to
  `output_dir/<phase>/`, including a `*_parameters.json` manifest.

## Reporting issues

Please include:

- The output of `fastmdx info`
- The command line you ran (or the Python code snippet)
- The full error message and traceback
- A minimal reproducing example if possible

## Code of conduct

By participating in this project, you agree to abide by our
[Code of Conduct](CODE_OF_CONDUCT.md).

## License

By contributing, you agree that your contributions will be licensed under
the project's [MIT License](LICENSE).
