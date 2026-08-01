# API reference

FastMDXplora's public Python API is centered on two classes, both importable
from the top-level package.

```python
import fastmdxplora as fastmdx

study = fastmdx.FastMDXplora(config="study.yml")
results = study.explore()
```

The Python API runs the same setup, simulation, analysis, and report engine as
the CLI. It does not need a separate Python dashboard implementation: enable
simulation telemetry in `options`, then serve the resulting output with
`fastmdx gui --output <run>` from another terminal. See the
[Beginner's guide](getting_started.md) for a complete example.

## FastMDXplora

The project-level orchestrator: the entry point for running a complete study
(any subset of the four phases) from Python.

```{eval-rst}
.. autoclass:: fastmdxplora.FastMDXplora
   :members:
   :undoc-members:
   :show-inheritance:
```

## AnalysisOrchestrator

Coordinates the analysis phase: trajectory loading, the analysis plan
(including auto-detected protein-ligand analyses), and writing results.

```{eval-rst}
.. autoclass:: fastmdxplora.AnalysisOrchestrator
   :members:
   :undoc-members:
   :show-inheritance:
```

## Package metadata

```{eval-rst}
.. autodata:: fastmdxplora.__version__
.. autodata:: fastmdxplora.__citation__
.. autodata:: fastmdxplora.__doi__
```
